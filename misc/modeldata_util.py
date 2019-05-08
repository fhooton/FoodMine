import pandas as pd
import numpy as np
import pickle
import spacy
from spacy.matcher import Matcher
import en_core_web_sm
import re
from time import time, sleep
import urllib.request as request
from lxml import etree

# For Microsoft Acedemic Graph
import http.client, urllib.request, urllib.parse, urllib.error, base64
import json

# For PDF Extraction
import chemdataextractor as cd
import tabula


# Takes in data from paper classifier and organizes/reforms sed data into the appropriate form for the model class
class ModelData():

	def __init__(self, structure):
		self.structure = structure

		self.nlp = en_core_web_sm.load()

		dicts = pickle.load(open('dicts.pkl', 'rb'))

		food_dict = list(set( dicts['food'] ))
		chem_dict = list(set( dicts['chem'] ))
		sci_dict = list(set( dicts['sci_name'] ))

		gen_dict = ['food', 'meat', 'vegetable', 'database']
		#food_dict = ['garlic', 'tomato']
		#chem_dict = ['pinoresinol', 'lariciresinol', 'isoflavones']

		search_dict = ['garlic', 'allium sativum', 'a. sativum']
		#search_dict = ['cocao', 'theobroma cacao', 't. cacao']
		self.search_matcher = self.__matcher_from_list__(search_dict)

		# Creates specialized matchers for food
		self.gen_matcher = self.__matcher_from_list__(gen_dict)
		self.food_matcher = self.__matcher_from_list__(food_dict)
		self.chem_matcher = self.__matcher_from_list__(chem_dict)
		self.sci_matcher = self.__matcher_from_list__(sci_dict)

		self.eval_sci_matches = True

	# Organizes data (in form pandas dataframe) into appropriate form for model
	def build_data(self, input_data, is_traindata = False):

		self.input_data = input_data
		self.is_traindata = is_traindata

		# Restructers data according to the specified data format
		if self.structure == 'search_selection':
			self.__search_selection__()

		elif self.structure == 'pdf':
			self.__pdf_struct__()

		else:
			print("Please enter acceptable structure")
			return


	# The initial structure of data for the model
	def __search_selection__(self):
		
		self.__init_preloaded_pubmed_vars__()

		self.data = pd.DataFrame()
		self.meta_data = pd.DataFrame()
		
		start = time()
		
		for idx, row in self.input_data.iterrows():

			if not idx % 200:
				print('# rows:', idx, round(time() - start) / 60, "min")

			# Detects the presence of measurement methodology word in abstract
			measurement_methods = ['spectrometry', 'chromatography', 'spectrophotometry']
			measurement_detection = self.__detect_measurement_methodology__(row['abstract'], row['mesh_terms'], measurement_methods)

			if row['abstract'] == None:
				row['abstract'] = ''

			# Mesh_terms might come in as a string from a csv file. Need to convert to list if so
			if not isinstance(row['mesh_terms'],list):
				row['mesh_terms'] = self.__mesh_string_to_list__(row['mesh_terms'])
			else:
				row['mesh_terms'] = row['mesh_terms']

			stc_abstract, _ = self.__matches__(self.search_matcher, row['abstract'])
			stc_mesh_terms, _ = self.__matches__(self.search_matcher, row['mesh_terms'])
			search_term_count = stc_abstract + stc_mesh_terms

			# Detects the presence of measurement methodology word in abstract
			#measurement_methods = ['spectrometry', 'chromatography', 'spectrophotometry']
			#measurement_detection = self.__detect_measurement_methodology__(row['abstract'], row['mesh_terms'], measurement_methods)
			
			pubmed_pack = self.__get_pubmed_features__(row)

			data_row = {
				'PMID' : row['PMID'],
				'search_term_count' : search_term_count
			}

			data_row.update(pubmed_pack)

			# Use variable classes if data is training data, otherwise just mark as useful for search
			if self.is_traindata:
				data_row['class'] = row['is_useful']
			else:
				data_row['class'] = 1

			# Add the existance of measurement methods as different features
			for method in measurement_methods:
				data_row[method] = measurement_detection[method]

			self.data = self.data.append(data_row, ignore_index = True)

			meta_row = {
				'PMID' : row['PMID']
				#'gen_terms' : gen_terms,
				#'food_terms' : food_terms,
				#'sci_terms' : sci_terms,
				#'chem_terms' : chem_terms
			}

			self.meta_data = self.meta_data.append(meta_row, ignore_index = True)

		self.data['PMID'] = self.data['PMID'].astype('int32')
		self.meta_data['PMID'] = self.meta_data['PMID'].astype('int32')

		self.data = self.data.fillna(0).set_index('PMID', drop = True)

		for col in self.data.columns:
			self.data[col] = pd.to_numeric(self.data[col])
		#self.meta_data = self.meta_data.set_index('PMID', drop = True)

		self.__normalize_feature_vectors__()

		#dictionary_condition = '((row["gen_term_count"] > 0) | (row["food_term_count"] > 0) | (row["chem_term_count"] > 0))'
		dictionary_condition = '((int(row["gen_term_count"]) > 0) + (int(row["food_term_count"]) > 0) + (int(row["chem_term_count"]) > 0) + (int(row["sci_term_count"]) > 0) > 1)'
		
		# Dynamically builds the measurement conditions based on different measurement methods
		measurement_condition = ''
		for method in measurement_methods:
			
			if measurement_condition == '':
				measurement_condition += '((int(row["' + method + '"]) == 1)'
			else:
				measurement_condition += ' | (int(row["' + method + '"]) == 1)'

		measurement_condition += ')'

		# Creates the logical pattern to be used as a filter function in the Model class
		self.extract_pattern = {
			#'(row["has_spectrometry"] == 1) | (row["has_chromatography"] == 1)' : 'True'
			dictionary_condition + ' & ' +  measurement_condition : 'True'
		}

	# Retrieves features from pubmed document information
	def __get_pubmed_features__(self, row):
		
		qualifier_perc = self.__qualifier_percent__(row)

		num_mesh_terms = len(row['mesh_terms'])
		abstract_length = len(row['abstract'].split())

		num_mass_in_abstract = len(re.findall(self.mass_unit_regex, row['abstract']))

		mesh_stats = self.__get_mesh_stats__(row['mesh_terms'])

		num_citations = self.__get_num_citations__(row['paper'])

		# Returns the term count and terms based on specific matcher
		gtc_abstract, gt_abstract = self.__matches__(self.gen_matcher, row['abstract'])
		gtc_mesh_terms, gt_mesh_terms = self.__matches__(self.gen_matcher, row['mesh_terms'])

		ftc_abstract, ft_abstract = self.__matches__(self.food_matcher, row['abstract'])
		ftc_mesh_terms, ft_mesh_terms = self.__matches__(self.food_matcher, row['mesh_terms'])

		ctc_abstract, ct_abstract = self.__matches__(self.chem_matcher, row['abstract'])
		#ctc_abstract, ct_abstract = self.__PubTator__(row['PMID'])
		ctc_mesh_terms, ct_mesh_terms = self.__matches__(self.chem_matcher, row['mesh_terms'])

		# int's of counts for dictionary terms in mesh terms and abstracts
		gen_term_count = gtc_abstract + gtc_mesh_terms
		food_term_count = ftc_abstract + ftc_mesh_terms
		chem_term_count = ctc_abstract + ctc_mesh_terms

		gen_terms = list(set( gt_abstract + gt_mesh_terms ))
		food_terms = list(set( ft_abstract + ft_mesh_terms ))
		chem_terms = list(set( ct_abstract + ct_mesh_terms ))

		if self.eval_sci_matches == True:
			stc_abstract, st_abstract = self.__matches__(self.sci_matcher, row['abstract'])
			stc_mesh_terms, st_mesh_terms = self.__matches__(self.sci_matcher, row['mesh_terms'])

			sci_term_count = stc_abstract + stc_mesh_terms
			sci_terms = list(set( st_abstract + st_mesh_terms ))

		pubmed_pack = {
			'gen_term_count' : gen_term_count,
			'food_term_count' : food_term_count,
			'sci_term_count' : sci_term_count,
			'chem_term_count' : chem_term_count,
			'qualifier_perc' : qualifier_perc,
			#'num_mesh_terms' : num_mesh_terms,
			#'abstract_length' : abstract_length,
			'num_mass_in_abstract' : num_mass_in_abstract,
			'num_citations' : num_citations
		}

		if len(mesh_stats) > 0:
			pubmed_pack.update(mesh_stats)

		return pubmed_pack

	# Calculate mesh term statistics from mesh description file
	def __get_mesh_stats__ (self, mesh_terms):
		
		if len(mesh_terms) == 0:
			return {}
		
		stat_df = pd.DataFrame({'mesh_term' : mesh_terms})

		stat_df = stat_df.merge(self.mesh_tree, how='left', on='mesh_term')

		sub_branch_levels = stat_df[stat_df.sub_branch_level.notnull()].sub_branch_level.tolist()

		analytical_sub_branch_count = self.sub_branch_counts(stat_df)

		branches = self.branch_counts(stat_df.main_branch.tolist())

		mean_branch = np.mean(sub_branch_levels)

		stat_pack = {
			'mean_mesh_branch_level' : mean_branch,
			'analytical_sub_branch_count' : analytical_sub_branch_count
		}

		stat_pack.update(branches)

		return stat_pack
	
	# Counts the top level branches in a set of mesh terms
	def branch_counts(self, branches):
		#b = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'V', 'Z']
		b = ['D', 'E', 'J']
		return {k : branches.count(k) for k in b}

	# Counts the number of occurances for specified sublevel
	def sub_branch_counts(self, stat_df, level = 1):
		# Chemistry Techniques, Analytical in mesh tree
		branch_to_compare = 'E05.196'
		count = 0

		for idx,row in stat_df.fillna('nan').iterrows():
			if len(row['mesh_code'].split('.')) > level:
				branch_at_level = '.'.join(row['mesh_code'].split('.')[0:level+1])

				if branch_at_level == branch_to_compare:
					count += 1

		return count

	# Uses Microsoft Acedemic Graph to get number of citations
	def __get_num_citations__(self, paper):
		# Get key from https://labs.cognitive.microsoft.com/en-us/project-academic-knowledge
		self.msft_apikey = '4c1e37e05c944c598b98c530eda9aaeb'

		# Code snippits at: https://dev.labs.cognitive.microsoft.com/docs/services/56332331778daf02acc0a50b/operations/56332331778daf06340c9666

		# Mode can be 'calchistogram', 'evaluate', 'interpret', or 'similarity'
		mode = 'interpret'

		query = paper

		headers = {
			# Request headers
			'Ocp-Apim-Subscription-Key': self.msft_apikey,
		}

		params = urllib.parse.urlencode({
			# Request parameters
			'query': query,
			#'complete': '0',
			'count': '10',
			#'offset': '{number}',
			#'timeout': '{number}',
			'model': 'latest',
			'attributes': 'Id',
		})

		try:
			conn = http.client.HTTPSConnection('api.labs.cognitive.microsoft.com')
			conn.request("GET", "/academic/v1.0/" + mode + "?%s" % params, "{body}", headers)
			response = conn.getresponse()
			data = response.read()
			#print(data)
			conn.close()
		except Exception as e:
			print("[Errno {0}] {1}".format(e.errno, e.strerror))

		loaded_interpret = json.loads(data)

		# If there are issiues with retrieving info, like no interpretations returned, return 0
		try:
			paper_query = loaded_interpret['interpretations'][0]['rules'][0]['output']['value']
			num_citations = self.__get_msft_attribute__(paper_query, 'CC')
		except:
			return 0

		return num_citations


	# Gets a attribute from paper of microsoft acedemic graph
	def __get_msft_attribute__(self, query, attr):
		# Mode can be 'calchistogram', 'evaluate', 'interpret', or 'similarity'
		mode = 'evaluate'

		headers = {
			# Request headers
			'Ocp-Apim-Subscription-Key': self.msft_apikey,
		}

		params = urllib.parse.urlencode({
			# Request parameters
			'expr': query,
			#'complete': '0',
			'count': '10',
			#'offset': '{number}',
			#'timeout': '{number}',
			'model': 'latest',
			'attributes': 'CC',
		})

		try:
			conn = http.client.HTTPSConnection('api.labs.cognitive.microsoft.com')
			conn.request("GET", "/academic/v1.0/" + mode + "?%s" % params, "{body}", headers)
			response = conn.getresponse()
			data = response.read()
			#print(data)
			conn.close()
		except Exception as e:
			print("[Errno {0}] {1}".format(e.errno, e.strerror))

		loaded_eval = json.loads(data)

		# See for all attributes: https://docs.microsoft.com/en-us/azure/cognitive-services/academic-knowledge/paperentityattributes
		num_citations = loaded_eval['entities'][0][attr]

		return num_citations


	# Determines what measurement methodolgies are present in paper details
	def __detect_measurement_methodology__(self, abstract, mesh_terms, measurement_methods):
		
		detection = {}

		if abstract == None:
			return {m : 0 for m in measurement_methods}
		elif mesh_terms == None:
			return {m : 0 for m in measurement_methods}

		for method in measurement_methods:
			
			# Returns 1 if word in statement, 0 otherwise
			detection[method] = self.__detect_word__(abstract, method)

			# also checks the mesh terms if the measurement method doesn't occur in the abstract
			if detection[method] == 0:
				detection[method] = self.__detect_word__(mesh_terms, method)

		return detection


	# Determines if a word is present in paper abstract or mesh terms
	def __detect_word__(self, text, word):

		if isinstance(text, list):
			text = " ".join(text)

		matcher = Matcher(self.nlp.vocab)

		pattern = [{'LOWER' : word}]
		
		# Adds the pattern to look for to the Matcher object
		matcher.add(word, None, pattern)

		# Specifically determines whether or not there is a match
		match = matcher(self.nlp(text))

		if len(match) > 0:
			return 1

		else:
			return 0


	# Counts the number of times a match appears in a text and return list of matches
	def __matches__(self, matcher, text):
		if isinstance(text, list):
			text = " ".join(text)

		text = self.nlp(text)

		matches = matcher(text)

		match_count = len(matches)

		match_terms = [self.nlp.vocab.strings[match_id] for match_id, _, _ in matches]

		return match_count, match_terms


	# Creates a matcher object that identifies lists of words for incoming text
	def __matcher_from_list__(self, dictionary):
		
		matcher = Matcher(self.nlp.vocab)

		for term in dictionary:

			#pattern = [{'LOWER' : term, 'OP' : '?'}, {'LOWER' : term + 's', 'OP' : '?'}]
			pattern = [{'LOWER' : term}]

			matcher.add(term, None, pattern)

		return matcher


	# Uses pubtator tool to return chemicals and chemical counts based on abstract
	# https://www.ncbi.nlm.nih.gov/research/bionlp/APIs/usage/
	def __PubTator__(self, PMID):
		concept = 'Chemical'

		url_end = concept + '/' + str(PMID) + '/BioC'

		url = 'https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/' + url_end

		with request.urlopen(url) as response:
			xml = response.read()

		root = etree.fromstring(xml)

		elements = root.findall('.//annotation')

		chem_terms = [i.text for i in elements]

		return len(chem_terms), chem_terms

	# Calculate the percentage of pre-determined qualifiers present in a paper's mesh terms
	def __qualifier_percent__(self, row):

		# Tope qualifiers that should be the most useful
		quals = ['analysis', 'chemistry', 'methods', 'metabolism']

		qual_count = sum( [row['qual_terms'].count(q) for q in quals] )

		if len(row['mesh_terms']) == 0:
			return 0
		else:
			return qual_count / len(row['mesh_terms'])

	# Creates different feature vectors for files with pdfs
	def __pdf_struct__(self):
		self.__get_pdf__filenames__()

		self.data = pd.DataFrame()
		self.meta_data = pd.DataFrame()

		self.__init_preloaded_pdf_vars__()
		self.__init_preloaded_pubmed_vars__()
		
		start = time()
		
		for idx, row in self.input_data.iterrows():

			if not idx % 200:
				print('# rows:', idx, round(time() - start) / 60, "min")
			
			# Detects the presence of measurement methodology word in abstract
			measurement_methods = ['spectrometry', 'chromatography', 'spectrophotometry']
			measurement_detection = self.__detect_measurement_methodology__(row['abstract'], row['mesh_terms'], measurement_methods)

			if row['abstract'] == None:
				row['abstract'] = ''

			# Mesh_terms might come in as a string from a csv file. Need to convert to list if so
			if not isinstance(row['mesh_terms'],list):
				row['mesh_terms'] = self.__mesh_string_to_list__(row['mesh_terms'])
			else:
				row['mesh_terms'] = row['mesh_terms']

			stc_abstract, _ = self.__matches__(self.search_matcher, row['abstract'])
			stc_mesh_terms, _ = self.__matches__(self.search_matcher, row['mesh_terms'])

			search_term_count = stc_abstract + stc_mesh_terms
			
			data_row = {
				'PMID' : row['PMID'],
				'search_term_count' : search_term_count
			}

			pdf_pack = self.__get_pdf_features__(row['filename'])
			data_row.update(pdf_pack)

			pubmed_pack = self.__get_pubmed_features__(row)
			data_row.update(pubmed_pack)

			# Use variable classes if data is training data, otherwise just mark as useful for search
			if self.is_traindata:
				data_row['class'] = row['is_useful']
			else:
				data_row['class'] = 1

			
			# Add the existance of measurement methods as different features
			for method in measurement_methods:
				data_row[method] = measurement_detection[method]
			

			self.data = self.data.append(data_row, ignore_index = True)
			

			meta_row = {
				'PMID' : row['PMID'],
			}

			self.meta_data = self.meta_data.append(meta_row, ignore_index = True)

		self.data['PMID'] = self.data['PMID'].astype('int32')
		self.meta_data['PMID'] = self.meta_data['PMID'].astype('int32')

		# Need to fill NA values with something for scikit to works
		self.data = self.data.fillna(0).set_index('PMID', drop = True)

		# Some features that are from other files are of type string, so need to set to numeric
		for col in self.data.columns:
			self.data[col] = pd.to_numeric(self.data[col]) 

		#self.__normalize_feature_vectors__()
		#self.meta_data = self.meta_data.set_index('PMID', drop = True)
		

	# Retrieves featrues based on reading the pdf
	def __get_pdf_features__(self, file):
		
		# Some files might not have header information
		try:
			headers = self.table_headers[file]
		except:
			headers = []
			print("Error retreiving header info. Perhaps no headers in pre-extracted header dictionary")

		# Some files might not have table information
		try:
			tables = self.table_data.set_index('file', drop = True).loc[file, 'tables']
		except:
			tables = []
			print("Error retreiving table info. Perhaps no tables in pre-extracted table dataframe")

		try:
			figure_headers = self.figure_headers[file]
		except:
			figure_headers = []
			print("Error retreiving figure header info. Perhaps no figure headers in pre-extracted header dictionary")


		# Calculates the total number of times a search term appeared in headers in a paper
		num_search_terms_in_headers = sum([
			len(self.__matches__(self.search_matcher, header)) for header in headers
		])

		num_mass_in_headers = sum([
			len(re.findall(self.mass_unit_regex, header)) for header in headers
		])

		# Retrieves dict associated with feature vectors from tabula extraction
		table_featuers = self.__sum_precomputed_table_features__(file)

		pdf_pack = {
			'num_tables' : round((len(headers) + len(tables)) / 2),
			'num_figures' : len(figure_headers),
			'num_search_terms_in_headers' : num_search_terms_in_headers,
			'num_mass_in_headers' : num_mass_in_headers
		}

		# Not all files have tables associated, and therefore they don't have table vectors either
		if len(table_featuers) is not 0:
			# List out the keys for the table features to be included in pdf_pack
			table_featuers_to_include = {
				'perc_numbers_in_tables' : table_featuers['perc_cells_numbers'],
				'perc_meanvar_in_tables' : table_featuers['perc_mean_variance_cells'],
				'search_count_in_tables' : table_featuers['search_count'],
				'perc_verbs_in_tables' : table_featuers['perc_verbs']
			}

			pdf_pack.update(table_featuers_to_include)

		return pdf_pack

	# Initializes data dictionarys, dataframes, and matchers for pdf files
	def __init_preloaded_pdf_vars__(self):
		# Read in files where these table headers and tables are already extracted to save time
		self.table_headers = pickle.load( open( "C:/Users/forresthooton/Documents/NetSci Research/Foodome Project/food_chemical_db/misc_save/table_headers.pkl", "rb" ) )
		self.figure_headers = pickle.load( open( "C:/Users/forresthooton/Documents/NetSci Research/Foodome Project/food_chemical_db/misc_save/figure_headers.pkl", "rb" ) )
		self.table_data = pickle.load( open( 'C:/Users/forresthooton/Documents/NetSci Research/Foodome Project/food_chemical_db/misc_save/table_data.pkl', "rb" ) )
		self.table_vectors = pickle.load( open( 'C:/Users/forresthooton/Documents/NetSci Research/Foodome Project/food_chemical_db/misc_save/table_vectors.pkl', "rb" ) )

		self.mass_unit_regex = r'[kmuµ(100)]?g[ ( per )/\\][(100)(100 )]?[kmuµ]?[gLl]'

	def __init_preloaded_pubmed_vars__(self):
		self.mesh_tree = pickle.load( open( "C:/Users/forresthooton/Documents/NetSci Research/Foodome Project/food_chemical_db/mesh_tree_counts.pkl", "rb" ) )
		self.mass_unit_regex = r'[kmuµ(100)]?g[ ( per )/\\][(100)(100 )]?[kmuµ]?[gLl]'

	# Combines values of different table from precomputed self.table_features variable into single dict
	def __sum_precomputed_table_features__(self, file):
		try:
			keys = list( self.table_vectors.loc[file, 'feature_vectors'][0].keys())
		except:
			return {}
		
		vectors_for_file = np.array([np.array(list(dic.values())) for dic in self.table_vectors.loc[file,'feature_vectors']])
		
		column_sumes = np.sum(vectors_for_file, axis = 0)

		tf_dict = { keys[i] : column_sumes[i] for i in range(len(column_sumes)) }

		return tf_dict
		

	# Appends filenames to input data df
	def __get_pdf__filenames__(self):
		# Filenames are stored in specifc file
		filename_df = pd.read_csv('C:/Users/forresthooton/Dropbox/PDFs/garlic_search.csv', encoding = 'latin1')[['PMID', 'filename']]
		# Only add filename column
		self.input_data = self.input_data.merge(filename_df, on = 'PMID')


	def __normalize_feature_vectors__(self):
		
		# Compute Normalization for each column
		for col in self.data.columns:
			# Avoid nomralizing class feature
			if col == 'class':
				continue

			mean = self.data[col].mean()
			stdev = self.data[col].std()

			if stdev == 0:
				norm = lambda x: (x - mean)
			else:
				norm = lambda x: (x - mean) / stdev

			self.data[col] = self.data[col].apply(norm)

	# Get random subset with replacement
	def sample(self, n = 10, replace = True):
		self.data = self.data.sample(n=n, replace=replace)

		#self.meta_data = self.meta_data.iloc[self.data.index]

	def __mesh_string_to_list__(self,mesh_string):
		mesh_terms = mesh_string.lstrip('[').rstrip(']').split('\'')
		mesh_terms = [t.strip() for t in mesh_terms if t.strip() != ',' and t != '']
		mesh_terms = [
			t if t != ', None' or t != 'None,' else None for t in mesh_terms
		]

		return mesh_terms


	#def copy(self):
	#	return self


