# Author: Forrest Hooton


import pandas as pd
import numpy as np
import pickle
import spacy
from spacy.matcher import Matcher
import en_core_web_sm
import re
from time import time, sleep
from lxml import etree
from tqdm import tqdm


class Filter():

	"""
		A class that calculates features for pubmed entry filtration
	"""


	def __init__(self):

		# Loads spacy language model
		self.nlp = en_core_web_sm.load()

		# Loads vocabulary from predefined dictionaries of words
		dicts = pickle.load(open('data/dicts.pkl', 'rb'))

		food_dict = list(set( dicts['food'] ))
		chem_dict = list(set( dicts['chem'] ))
		sci_dict = list(set( dicts['sci_name'] ))

		# Sets a dictionary to capture broad, possibly relevant concepts
		gen_dict = ['food', 'meat', 'vegetable', 'database']

		# Creates specialized matchers for different dictionaries
		self.gen_matcher = self.__matcher_from_list__(gen_dict)
		self.food_matcher = self.__matcher_from_list__(food_dict)
		self.chem_matcher = self.__matcher_from_list__(chem_dict)
		self.sci_matcher = self.__matcher_from_list__(sci_dict)

		self.eval_sci_matches = True

		# Specifies measurement methods to search for and creates matchers for those methods
		measurement_methods = ['spectrometry', 'chromatography', 'spectrophotometry']
		self.measurement_matchers = {m : self.__matcher_from_list__([m]) for m in measurement_methods}


	def filter(self, data):
		"""
			Takes in search from PubMed with the specified features from pubmed_util and filters results based one
			pre-specified criteria.

			Parameters
			-----------------
			data : pd.DataFrame
				search results data from PubMed


			Returns
			-----------------
			selected_articles: pd.DataFrame
				PubMed entries that met the specified criteria and were not filtered out.
		"""

		self.input_data = data

		# Creates word counts for filtration
		self.build_features()

		selected_articles = pd.DataFrame(columns=['PMID'])

		# Iterates over each PubMed entry and filters out results that don't meed criteria
		for PMID, row in self.data.iterrows():

			if self.__eval_conditionals__(row):

				selected_articles = selected_articles.append({'PMID' : PMID}, ignore_index = True)

		selected_articles['PMID'] = selected_articles['PMID'].astype('int32')

		return selected_articles


	def build_features(self, input_data=None, is_traindata=False):
		"""
			Creates features. In this instance word count frequencies.

			Parameters
			-----------------
			None


			Returns
			-----------------
			None
		"""
		
		if input_data is not None:
			self.input_data = input_data

		self.data = pd.DataFrame()
		
		start = time()
		print('Creating features...')

		# Builds features for each row of data
		for idx, row in tqdm(self.input_data.iterrows()):

			# if not idx % 200:
			# 	print('# rows:', idx, round(time() - start) / 60, "min")

			if row['abstract'] == None:
				abstract = ''
			else:
				abstract = row['abstract']
			
			mesh_terms = ' '.join(row['mesh_terms'])

			# Combine all info into a single string
			text = ' '.join([abstract, mesh_terms])


			# Returns 1 per method if the measurement method is in the abstract or mesh terms, and 0 otherwise
			measurement_detection = {
				method : self.__detect_word_presence__(text, matcher)
				for method, matcher in self.measurement_matchers.items()
			}

			# Retrieves features from pubmed
			pubmed_features = self.__get_pubmed_features__(text)

			data_row = {
				'PMID' : row['PMID']
			}

			data_row.update(pubmed_features)

			# Use variable classes if data is training data, otherwise just mark as useful for search
			if is_traindata:
				data_row['class'] = row['is_useful']
			else:
				data_row['class'] = 1

			# Add the existence of measurement methods as different features
			for method in self.measurement_matchers.keys():
				data_row[method] = measurement_detection[method]

			# Updates data with new information
			self.data = self.data.append(data_row, ignore_index = True)

		self.data['PMID'] = self.data['PMID'].astype('int32')

		self.data = self.data.fillna(0).set_index('PMID', drop = True)

		for col in self.data.columns:
			self.data[col] = pd.to_numeric(self.data[col])

		# Specifies criteria for filtration
		dictionary_condition = '((int(row["gen_term_count"]) > 0) + (int(row["food_term_count"]) > 0) + (int(row["chem_term_count"]) > 0) + (int(row["sci_term_count"]) > 0) > 1)'
		
		# Dynamically builds the measurement conditions based on different measurement methods
		measurement_condition = ''
		for method in self.measurement_matchers.keys():
			
			if measurement_condition == '':
				measurement_condition += '((int(row["' + method + '"]) == 1)'
			else:
				measurement_condition += ' | (int(row["' + method + '"]) == 1)'

		measurement_condition += ')'

		# Creates the logical pattern to be used as a filter function in the Model class
		self.extract_pattern = {
			dictionary_condition + ' & ' +  measurement_condition : 'True'
		}


	def __get_pubmed_features__(self, text):
		"""
			Calculates count frequency features from PubMed entry abstract and mesh terms.

			Parameters
			-----------------
			row : pd.Series
				Single row containing information on PubMed entry


			Returns
			-----------------
			pubmed_features : dict
				Dictionary containing the count frequencies for general, food, scientific name, and chemical terms.
		"""
		
		# Count the terms that occur in each dictionary
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


	# Determines if a word is present in paper abstract or mesh terms
	def __detect_word_presence__(self, text, matcher):
		"""
			Calculates count frequency features from PubMed entry abstract and mesh terms.

			Parameters
			-----------------
			text : string
				Text to search for specified matcher pattern
			
			matcher : spacy Matcher object
				Specific matcher with defined rules to use with the text input

			Returns
			-----------------
			presence of pattern : int
				Returns 1 if there is any presence of the matcher pattern, and 0 otherwise.
		"""

		# Retrieves the matches from text
		match = matcher(self.nlp(text))

		# Returns 1 if there are any matches, else 0
		if len(match) > 0:
			return 1

		else:
			return 0


	def __matches__(self, matcher, text):
		"""
			Calculates count frequency features from PubMed entry abstract and mesh terms.

			Parameters
			-----------------
			matcher : spacy Matcher object
				Specific matcher with defined rules to use with the text input

			text : string
				Text to search for specified matcher pattern

			Returns
			-----------------
			match_count : int
				Frequency of matches between the text input and matcher pattern
		"""

		text = self.nlp(text)

		matches = matcher(text)

		match_count = len(matches)

		return match_count


	def __matcher_from_list__(self, dictionary):
		"""
			Creates the spacy matcher object using an input dictionary of terms

			Parameters
			-----------------
			dictionary : list
				List of terms to be included as matching patterns


			Returns
			-----------------
			matcher : spacy Matcher object
				Frequency of matches between the text input and matcher pattern
		"""

		matcher = Matcher(self.nlp.vocab)

		for term in dictionary:

			#pattern = [{'LOWER' : term, 'OP' : '?'}, {'LOWER' : term + 's', 'OP' : '?'}]
			pattern = [{'LOWER' : term}]

			matcher.add(term, None, pattern)

		return matcher


	def __eval_conditionals__(self, row):
		"""
			Creates the spacy matcher object using an input dictionary of terms

			Parameters
			-----------------
			row : pd.Series
				Feature row with which to evaluate conditionals

			Returns
			-----------------
			eval(action) : 
				Evaluation of action when an extraction pattern condition is met
		"""
		
		for condition, action in self.extract_pattern.items():

			if eval(condition):
				return eval(action)





