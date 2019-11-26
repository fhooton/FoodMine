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


class Filter():


	def __init__(self):

		self.nlp = en_core_web_sm.load()

		dicts = pickle.load(open('data/dicts.pkl', 'rb'))

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


	# Appends paper based on whether or not it fits some criteria defined in the ModelData class
	def filter(self, data):

		self.input_data = data

		self.__build_features__()

		selected_articles = pd.DataFrame()

		for PMID, row in self.data.iterrows():

			if self.__eval_conditionals__(row):

				selected_articles = selected_articles.append({'PMID' : PMID}, ignore_index = True)

		selected_articles['PMID'] = selected_articles['PMID'].astype('int32')

		return selected_articles


	# The initial structure of data for the model
	def __build_features__(self):
		
		self.data = pd.DataFrame()
		
		start = time()
		
		for idx, row in self.input_data.iterrows():

			if not idx % 200:
				print('# rows:', idx, round(time() - start) / 60, "min")

			# Detects the presence of measurement methodology word in abstract
			measurement_methods = ['spectrometry', 'chromatography', 'spectrophotometry']
			measurement_detection = self.__detect_measurement_methodology__(row['abstract'], row['mesh_terms'], measurement_methods)

			if row['abstract'] == None:
				row['abstract'] = ''

			pubmed_features = self.__get_pubmed_features__(row)

			data_row = {
				'PMID' : row['PMID']
			}

			data_row.update(pubmed_features)

			# Add the existance of measurement methods as different features
			for method in measurement_methods:
				data_row[method] = measurement_detection[method]

			self.data = self.data.append(data_row, ignore_index = True)

		self.data['PMID'] = self.data['PMID'].astype('int32')

		self.data = self.data.fillna(0).set_index('PMID', drop = True)

		for col in self.data.columns:
			self.data[col] = pd.to_numeric(self.data[col])

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

		abstract = ' '.join(row['abstract'])
		mesh_terms = ' '.join(row['mesh_terms'])

		text = ' '.join([abstract, mesh_terms])
		
		gen_term_count = self.__matches__(self.gen_matcher, text)
		food_term_count = self.__matches__(self.food_matcher, text)
		chem_term_count = self.__matches__(self.chem_matcher, text)

		if self.eval_sci_matches == True:
			sci_term_count = self.__matches__(self.sci_matcher, text)

		pubmed_features = {
			'gen_term_count' : gen_term_count,
			'food_term_count' : food_term_count,
			'sci_term_count' : sci_term_count,
			'chem_term_count' : chem_term_count
		}

		return pubmed_features


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

		text = self.nlp(text)

		matches = matcher(text)

		match_count = len(matches)

		return match_count


	# Creates a matcher object that identifies lists of words for incoming text
	def __matcher_from_list__(self, dictionary):
		
		matcher = Matcher(self.nlp.vocab)

		for term in dictionary:

			#pattern = [{'LOWER' : term, 'OP' : '?'}, {'LOWER' : term + 's', 'OP' : '?'}]
			pattern = [{'LOWER' : term}]

			matcher.add(term, None, pattern)

		return matcher


	def __eval_conditionals__(self, row):

		# extract_pattern is a dictionary in ModelData where the keys are conditions
		# and values actions if the condition is met
		for condition, action in self.extract_pattern.items():

			if eval(condition):
				return eval(action)





