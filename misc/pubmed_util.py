import urllib.request as request
from lxml import etree
import pandas as pd
import numpy as np
import math
from time import time

# Imports from directory
import model_util


# Recieves search terms to classify usefulness of pubmed documents from search results
def filter_results(search_terms):

	doc_ids = search_pubmed(search_terms)

	print('ids', len(doc_ids))

	start = time()

	print('here 1')
	# pandas dataframe of document info
	doc_info = retrieve_doc_info(doc_ids)
	print(search_terms, "Document Info", "(" + str(len(doc_info)) + " entries)", "Retrieved in", (time() - start) / 60, "min")

	start = time()
	model_data = ModelData('search_selection')
	model_data.build_data(doc_info)
	print("Data Converted in", (time() - start) / 60, "min")

	output = model(model_data)

	output_info = retrieve_doc_info(output['PMID'].tolist())

	#output_info.to_pickle('temp_sub_search.pkl')

	#print(output_info['PMID'])

	return output_info

# Enters search terms into pubmed database to return document ID's
def search_pubmed(search_terms):

	url = construct_url(search_terms, 'search')

	with request.urlopen(url) as response:
		xml = response.read()

	root = etree.fromstring(xml)

	# Recursively gets all objects where the tag is Id
	elements = root.findall('.//Id')

	# Converts all lxml objects to their text values
	ids = [i.text for i in elements]

	return ids


# Retrieves document (paper) info using pubmed paper ids
def retrieve_doc_info(ids):
	# Can't query too much in a single query, so divides larger id lists into seperate queries
	num_loops = int(math.ceil(len(ids) / 100))

	print('here 2')
	# Have to split requests larger than 100 documents to keep it within url size
	ids = divide_list(ids, num_loops)

	documents = []

	print('here 3')
	# Retrieves xml data from pubmed
	for i in ids:
		url = construct_url(i, 'document')

		with request.urlopen(url) as response:
			xml = response.read()

		root = etree.fromstring(xml)

		documents = documents + root.findall('PubmedArticle')

	info = pd.DataFrame()

	print('here 4')
	for document in documents:

		doc_id = int(document.find('.//PMID').text)
		
		paper = document.find('.//ArticleTitle').text

		journal = document.find('.//Title').text

		year = document.find('.//Year').text

		if document.find('.//AbstractText') is not None:
			abstract = document.find('.//AbstractText').text
		else:
			abstract = None

		mesh_terms = []
		mesh_UIds = []
		qual_terms = []
		qual_UIds = []

		for mesh_section in document.findall('.//MeshHeading'):
			mesh_terms.append(mesh_section.find('.//DescriptorName').text)
			mesh_UIds.append(mesh_section.find('.//DescriptorName').attrib['UI'])

			if mesh_section.find('.//QualifierName') is not None:
				qual_terms.append(mesh_section.find('.//QualifierName').text)
				qual_UIds.append(mesh_section.find('.//QualifierName').attrib['UI'])
			else:
				qual_terms.append(None)
				qual_UIds.append(None)

		new_row = {
			'PMID' : doc_id,
			'paper' : paper,
			'journal' : journal,
			'year' : year,
			'abstract' : abstract,
			'mesh_terms' : mesh_terms,
			'mesh_UIds' : mesh_UIds,
			'qual_terms' : qual_terms,
			'qual_UIds' : qual_UIds,
			'webpage' : 'https://www.ncbi.nlm.nih.gov/pubmed/' + str(doc_id)
		}
		
		info = info.append(new_row, ignore_index = True)

	print('here 5')
	info['PMID'] = info['PMID'].astype('int32')

	return info.reset_index(drop = True)


def pubchem_synonym_info(chem_name):
	url = construct_url(chem_name, 'synonym')
	
	with request.urlopen(url) as response:
		xml = response.read()

	root = root = etree.fromstring(xml)
	compound_id = float(root.findall(".//{http://pubchem.ncbi.nlm.nih.gov/pug_rest}CID")[0].xpath('.//text()')[0])
	compound_name = str(root.findall(".//{http://pubchem.ncbi.nlm.nih.gov/pug_rest}Synonym")[0].xpath('.//text()')[0])

	return compound_id, compound_name


def pubchem_SMILE(chem_id):
	url = construct_url(chem_id, 'SMILE')

	with request.urlopen(url) as response:
		xml = response.read()

	root = root = etree.fromstring(xml)
	SMILE = root.findall(".//{http://pubchem.ncbi.nlm.nih.gov/pug_rest}CanonicalSMILES")[0].xpath('.//text()')[0]

	return SMILE


# Constructs appropriate url for pubmed api from search terms
def construct_url(url_input, query_type, num_results = 1000000):
	
	# Constructs url for search query
	if query_type == 'search':
		base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term='

		adjusted_terms = [s.replace(" ", "%20") for s in url_input]

		term_url = '%20AND%20'.join(adjusted_terms)

		# Maybe include later
		results_num_url = '&retmax=' + str(num_results)

		return base_url + term_url + results_num_url

	elif query_type == 'document':
		base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id='

		doc_urls = ""
		for i in url_input:
			if isinstance(i, str): 
				doc_urls = doc_urls + "," + i
			else:
				doc_urls = doc_urls + "," + str(i)

		url = base_url + doc_urls.lstrip(",") + '&retmode=xml'

		return url

	elif query_type == 'synonym':
		url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + url_input + "/synonyms/XML"
		return url

	elif query_type == 'SMILE':
		url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + url_input + "/property/CanonicalSMILES/XML"
		return url


# Divides doc ids for larger paper queries in retrieve_doc_info()
def divide_list(ids, num_divisions):

	split_ids = np.array_split(np.asarray(ids), num_divisions)
	split_ids = [np.ndarray.tolist(split_ids[i]) for i in range(len(split_ids))]

	return split_ids