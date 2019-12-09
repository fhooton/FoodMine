# FoodMine

FoodMine is a food database project from the Barabasi Lab of the Network Science Institute (Northeastern University). The objective of this project is to mine research litereature to uncover chemical contents in food. The code in this repository was used to generate the food pilots, plots, and other information in **FoodMine: Exploring Food Contents in Scientific Literature** (Hooton, Menichetti, Barabási). Additionally, we have included the data from out initial experiments.

# Files

#### Paper_Screening.ipynb
Notebook to search through the PubMed database and filter out search results.

#### Data_Statistics.ipynb
Notebook to compare clean raw data, compare databases, and vizualize comparisons.

#### Molecule_Embedding.ipynb
Notebook to retrive molecule smiles, embed compounds, and vizualize embeddings.

#### Paper_Citations.ipynb
Notebook to analyze the citation overlap between CTD and the papers gathered in FoodMine.

#### misc
Miscellaneous functions and classes to serve as helpers in notebooks. 

* **pubmed_util.py**
Holds functions to interact with PubMed API for the purposes of our research.

* **filter.py**
Contains class Filter to filter out PubMed search results.

* **collected_data_handling.py**
Helper functions to clean and aggregate information into FoodMine database.

* **tools/chemidr**
Suite of tools developed to interact with the PubChem API.

#### data
Folder that holds raw data from paper data collection. **Note that several files used in our analysis are too large to be uploaded to github**

