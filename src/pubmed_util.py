# Author: Forrest Hooton


import requests
from lxml import etree
import pandas as pd
import numpy as np
import math
import time

# Imports from directory
from .filter import Filter


def filter_results(search_terms):
	"""
		Receives search terms to classify usefulness of pubmed documents from search results.

		Parameters
		-----------------
		search_terms : list
			List of terms to include in PubMed query


		Returns
		-----------------
		output_info: pd.DataFrame
			PubMed entries that met the specified criteria with metadata attached
	"""

	# Retrieves doc_ids of search terms
	doc_ids = search_pubmed(search_terms)

	print('ids', len(doc_ids))

	start = time.time()

	# pandas dataframe of document info
	doc_info = retrieve_doc_info(doc_ids)
	print(search_terms, "Document Info", "(" + str(len(doc_info)) + " entries)", "Retrieved in", (time.time() - start) / 60, "min")

	# Filters out results deemed irrelevant
	start = time.time()
	filt = Filter()
	output = filt.filter(doc_info)
	print("Data Converted in", (time.time() - start) / 60, "min")

	# Reincorporate metadata
	output_info = retrieve_doc_info(output['PMID'].tolist())

	return output_info


def search_pubmed(search_terms):
	"""
		Receives search terms to classify usefulness of pubmed documents from search results.

		Parameters
		-----------------
		search_terms : list
			List of terms to include in PubMed query


		Returns
		-----------------
		output_info: pd.DataFrame
			PubMed entries that met the specified criteria with metadata attached
	"""

	# Gets url to retrieve information from search
	url = construct_url(search_terms, 'search')

	xml = __safe_urlopen__(url)

	root = etree.fromstring(xml)

	# Recursively gets all objects where the tag is Id
	elements = root.findall('.//Id')

	# Converts all lxml objects to their text values
	ids = [i.text for i in elements]

	return ids


def retrieve_doc_info(ids):
	"""
		Retrieves document (paper) info using pubmed paper ids.

		Parameters
		-----------------
		ids : list
			List of PubMed ids to retrieve PubMed entry information


		Returns
		-----------------
		info: pd.DataFrame
			PubMed ids with metadata attached
	"""

	# Can't query too much in a single query, so divides larger id lists into separate queries
	num_loops = int(math.ceil(len(ids) / 100))

	# Have to split requests larger than 100 documents to keep it within url size
	ids = divide_list(ids, num_loops)

	documents = []

	# Retrieves xml data from pubmed
	for i in ids:
		url = construct_url(i, 'document')

		xml = __safe_urlopen__(url)

		root = etree.fromstring(xml)

		documents = documents + root.findall('PubmedArticle')

	info = pd.DataFrame()

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

	info['PMID'] = info['PMID'].astype('int32')

	return info.reset_index(drop = True)


def pubchem_synonym_info(chem_name):
	"""
		Retrieves compound id and first compound synonym name from PubChem based on a queried chemical.

		Parameters
		-----------------
		chem_name : string
			Name of chemical


		Returns
		-----------------
		compound_id : int
			PubChem ID of queried chemical

		compound_name : string
			First name of compound listed in PubChem synonyms list
	"""

	# Creates synonym query url
	url = construct_url(chem_name, 'synonym')
	
	xml = __safe_urlopen__(url)

	root = etree.fromstring(xml)

	# Extracts the compound ID and first synonym name
	compound_id = float(root.findall(".//{http://pubchem.ncbi.nlm.nih.gov/pug_rest}CID")[0].xpath('.//text()')[0])
	compound_name = str(root.findall(".//{http://pubchem.ncbi.nlm.nih.gov/pug_rest}Synonym")[0].xpath('.//text()')[0])

	return compound_id, compound_name


def pubchem_SMILE(chem_id):
	"""
		Retrieves compound SMILE using pubchem ID.

		Parameters
		-----------------
		chem_id : int
			Compound pubchem ID


		Returns
		-----------------
		SMILE : string
			SMILE corresponding to input compound ID
	"""
	url = construct_url(chem_id, 'SMILE')

	xml = __safe_urlopen__(url)

	root = root = etree.fromstring(xml)

	# Extracts the compound SMILE
	SMILE = root.findall(".//{http://pubchem.ncbi.nlm.nih.gov/pug_rest}CanonicalSMILES")[0].xpath('.//text()')[0]

	return SMILE


# Constructs appropriate url for pubmed api from search terms
def construct_url(url_input, query_type, num_results = 1000000):
	"""
		Constructs the url for various pubchem queries

		Parameters
		-----------------
		url_input : string or list depending on query_type
			The input for a pubchem query such as pubchem id, pubmed id, or list of search terms

		query_type : string
			Specifies which search method should be used... can be "search", "document", "synonym", or "SMILE"

		num_results : int
			Top number of results to keep from a search query


		Returns
		-----------------
		url : string
			Url to pass to request
	"""
	
	# Constructs url for search query from list of search terms
	if query_type == 'search':
		base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term='

		# Replace spaces with something that works with as a url
		adjusted_terms = [s.replace(" ", "%20") for s in url_input]

		# Join separate search queries
		term_url = '%20AND%20'.join(adjusted_terms)

		# Cap the number of results
		results_num_url = '&retmax=' + str(num_results)

		return base_url + term_url + results_num_url

	# Constructs url for document query from list of pubmed document ids
	elif query_type == 'document':
		base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id='

		# Handles either string or integer pubmed ids
		doc_urls = ""
		for i in url_input:
			if isinstance(i, str): 
				doc_urls = doc_urls + "," + i
			else:
				doc_urls = doc_urls + "," + str(i)

		url = base_url + doc_urls.lstrip(",") + '&retmode=xml'

		return url
	
	# Constructs url for synonym search from a chemical string
	elif query_type == 'synonym':
		url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + url_input + "/synonyms/XML"
		return url

	# Constructs url for SMILE retrival from a pubchem id string
	elif query_type == 'SMILE':
		url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + url_input + "/property/CanonicalSMILES/XML"
		return url

	else:
		print('Please enter a valid query type ("search", "document", "synonym", or "SMILE")')


# Divides doc ids for larger paper queries in retrieve_doc_info()
def divide_list(ids, num_divisions):
	"""
		Splits a single list of ids into num_divisions numbers of separate lists

		Parameters
		-----------------
		ids : list
			list of pubmed ids

		num_divisions : int
			Number of divisions in which to partition lists


		Returns
		-----------------
		split_ids : 2d list
			list of lists of pubchem ids
	"""

	split_ids = np.array_split(np.asarray(ids), num_divisions)
	split_ids = [np.ndarray.tolist(split_ids[i]) for i in range(len(split_ids))]

	return split_ids


def __safe_urlopen__(url):
    """
        Retrieves information from url query without throwing errors, such as for values that do not
        exist or excessive querying. Designed for Pubchem and Pubmed apis

        Input
        ----------------------------------------------------------------
        url : str
            url to query

        Returns
        ----------------------------------------------------------------
        response.content : str (maybe bytes) or None
            Returns the response of a url query if it exists, else None
    """
    try:
        response = requests.get(url)
    except TimeoutError:
        time.sleep(.5)
        return __safe_urlopen__(url)

    if response.status_code == 200: # Successful
        return response.content

    elif response.status_code == 429: # Too many requests
        # print('Retrying...')
        time.sleep(.5)
        return __safe_urlopen__(url)

    elif response.status_code == 502: # Bad Gateway
        # print('Retrying...')
        time.sleep(1)
        return __safe_urlopen__(url)

    elif response.status_code == 404: # PUGREST.NotFound (aka doesn't exist)
        return None