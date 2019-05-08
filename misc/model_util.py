import pandas as pd
import numpy as np
from imblearn.over_sampling import SMOTE
import pickle

# sklearn models
from sklearn import svm
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier

# sklearn metrics
from sklearn.model_selection import cross_validate, train_test_split
from sklearn.metrics import classification_report,confusion_matrix

# sklearn feature selection
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2

# Trains and loads models from ModelData inputs
class Model():

	def __init__(self, model_type):

		self.model_type = model_type
		self.is_train = False

		#self.model = 'undeclared'

	# Wrapper function for all models that reroutes data based on model designation in params
	def __call__(self, ModelData, loadmodel = False):

		# Extracts based on pre specified patterns in ModelData
		if self.model_type == 'simple_extract':
			
			self.data = ModelData
			selected_articles = self.__simple_extract__()

			return selected_articles

		self.data_structure = ModelData.structure

		self.Y_data = ModelData.data['class']
		self.X_data = ModelData.data.drop('class', axis = 1)

		#self.__auto_feature_selection__()
		
		filename = self.__make_filepath__()

		# Initialize the model differently based on whether or not it should be trained or loaded
		if self.is_train:
			self.__get_model__()

			smote_x, smote_y = SMOTE().fit_sample(self.X_data, self.Y_data)
			self.X_data = pd.DataFrame(smote_x, columns = self.X_data.columns)
			#print(self.X_data.shape)
			self.Y_data = pd.Series(smote_y)

			self.__report__()

			self.model.fit(self.X_data, self.Y_data)

			pickle.dump(self.model, open(filename, 'wb'))

		else:

			if loadmodel is True:
				self.load_model()

			elif self.model is 'undeclared':
				print('Please run a command to initialize a model')


	

	# Returns initialized model based on model type
	def __get_model__(self):

		if self.model_type == 'logistic_regression':
			self.model = LogisticRegression()

		elif self.model_type == 'svm':
			self.model = svm.LinearSVC()

		elif self.model_type == 'decision_tree':
			self.model = DecisionTreeClassifier()

		elif self.model_type == 'random_forest':
			self.model = RandomForestClassifier(n_estimators = 100, criterion = 'entropy')

		elif self.model_type == 'xgboost':
			self.model = XGBClassifier()

		elif self.model_type == 'nn':
			hidden_layer_sizes = (len(self.X_data.values.tolist()[0]) * 2, len(self.X_data.values.tolist()[0]) * 2)
			self.model = MLPClassifier(solver='lbfgs', activation = 'logistic', alpha=1e-5, hidden_layer_sizes=hidden_layer_sizes)

		elif self.model_type == 'k_neighbors':
			self.model = KNeighborsClassifier(n_neighbors=5)

		else:
			print('Please declare the model with a valid model type:')
			print('logistic_regression\tsvm\trandom_forest\tnn\tk_neighbors')


	def load_model(self, path):
		try:
			self.model = pickle.load(open(path, 'rb'))
		except:
			print('Error loading model. File may not exist.')

	# Automactically filter features
	def __auto_feature_selection__(self):
		# Need to select the columns for the training data, but use the selected columns for other data
		if self.is_train:
			x_new = SelectKBest(chi2, k=15).fit_transform(self.X_data, self.Y_data)

			x_new_df = pd.DataFrame(x_new)

			col_subset = []

			for col in x_new_df.columns:
				for col2 in self.X_data.columns:

					if x_new_df[col].tolist() == self.X_data[col2].tolist():
						col_subset.append(col2)

			self.X_data = self.X_data[col_subset]
			self.filtered_cols = col_subset

		else:
			self.X_data = self.X_data[self.filtered_cols]

		print("Filtered Features\n", self.filtered_cols)



	# Function takes in training data to train model type
	# Meant to be used in conjuction with load_training_data()
	def train(self, temp):
		self.is_train = True

		#self.__call__(self.training_data)
		self.__call__(temp)

		self.is_train = False


	# Appends paper based on whether or not it fits some criteria defined in the ModelData class
	def __simple_extract__(self):

		selected_articles = pd.DataFrame()

		for PMID, row in self.data.data.iterrows():

			if self.__eval_conditionals__(row):

				selected_articles = selected_articles.append({'PMID' : PMID}, ignore_index = True)

		selected_articles['PMID'] = selected_articles['PMID'].astype('int32')

		# Add metadata to outgoing selected articles
		selected_articles = pd.merge(selected_articles, self.data.meta_data, how = 'left', on = 'PMID')

		return selected_articles


	# Evaluates criteria in ModelData class
	def __eval_conditionals__(self, row):

		# extract_pattern is a dictionary in ModelData where the keys are conditions
		# and values actions if the condition is met
		for condition, action in self.data.extract_pattern.items():

			if eval(condition):
				return eval(action)


	# Creates the filepath to save and load models
	def __make_filepath__(self):
		return 'saved_models/' + self.model_type + '-' + self.data_structure + '.pkl'


	# Helper function for all models to keep track of reporting for best practices and research purposes
	def __report__(self):
		print("Sample Size:", len(self.Y_data), '\n')

		scores = cross_validate(self.model, self.X_data, self.Y_data, cv=5)

		test_scores = scores['test_score']
		times = scores['fit_time']

		print("avg score:", np.mean(test_scores))
		print("avg time:", np.mean(times))

		x_train, x_test, y_train, y_test = self.balanced_train_test_split(self.X_data.values.tolist(), self.Y_data.tolist())

		"""
		print('\nBalanced data partition')

		for clas in list(set( self.Y_data.tolist() )):
			print(clas)
			print('Training:', y_train.count(clas) / len(y_train))
			print('Testing:', y_test.count(clas) / len(y_test))

		print()
		"""

		model = self.model

		model.fit(x_train, y_train)

		predictions = model.predict(x_test)

		print('Confusion Matrix')
		print(confusion_matrix(y_test,predictions), '\n')

		print('Classification Report')
		print(classification_report(y_test,predictions))

	def balanced_train_test_split(self, x, y, test_size = .4):
		x_train, x_test, y_train, y_test = [], [], [], []

		for clas in list(set( y )):
			x_temp = [x[i] for i, e in enumerate(y) if e == clas]
			y_temp = [clas] * len(x_temp)

			t_x_train, t_x_test, t_y_train, t_y_test = train_test_split(x_temp, y_temp, test_size = test_size)

			x_train += t_x_train
			x_test += t_x_test
			y_train += t_y_train
			y_test += t_y_test

		return np.array(x_train), np.array(x_test), np.array(y_train), np.array(y_test)


	def load_training_data(self, model_data):
		#self.training_data = pd.read_csv("C:/Users/forresthooton/Dropbox/PDFs/dbf/garlic_scoring.csv", encoding='latin1')

		training_df = pd.read_csv("C:/Users/forresthooton/Dropbox/PDFs/dbf/garlic_scoring.csv", encoding='latin1')
		
		# Only keep data that has already been labeled
		training_df = training_df[training_df['is_useful'].notnull()]
		#training_df = training_df[training_df['is_useful'] != '?']

		# Option to change extra class of useful but not quantified to 1
		training_df['is_useful'] = training_df['is_useful'].replace(2, 1, regex=True)

		"""
		phenol_data = pd.read_csv('phenol_explorer_training_data.csv')
		phenol_data['abstract'] = phenol_data['abstract'].fillna('')
		training_df = pd.concat([training_df, phenol_data], axis=0, ignore_index = True)
		#return training_df
		"""

		"""
		# Handle code differently if phenol explorer training data is included
		if len(training_df[training_df['is_useful'] != '0']) > len(training_df[training_df['is_useful'] == '0']):
			t1 = training_df[training_df['is_useful'] == '0']
			t2 = training_df[
				(training_df['is_useful'] != '0') & (training_df['phenol_td'] == True)
			].copy().sample(n=len(t1))

			# Create seperate holdout set for data from search to make sure learning transfers correctly
			self.pos_holdout = training_df[(training_df['is_useful'] != '0') & (training_df['phenol_td'].isnull())]
		
		else:
			t1 = training_df[training_df['is_useful'] != '0']
			t2 = training_df[training_df['is_useful'] == '0'].copy().sample(n=round(len(t1)*1.5))

		training_df = pd.concat([t1, t2], axis=0, ignore_index = True)
		#"""

		self.training_data = model_data
		self.training_data.build_data(training_df, is_traindata = True)
