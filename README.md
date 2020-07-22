# FoodMine

FoodMine is a food database project from the Barabasi Lab of the Network Science Institute (Northeastern University). The objective of this project is to mine research litereature to uncover chemical contents in food. The code in this repository was used to generate the food pilots, plots, and other information in [FoodMine: Exploring Food Contents in Scientific Literature](https://www.biorxiv.org/content/10.1101/2019.12.17.880062v1) (Hooton, Menichetti, Barabási). Additionally, we have included the data from out initial experiments.

Run the [python notebooks](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/) in the following order to replicate the FoodMine process.
 - Paper_Screening.ipynb
 - Data_statistics.ipynb (this notebook is key to running other files as it downloads relevant data hosted online)
 - Molecule_Embedding.ipynb
 - Paper_Citations.ipynb (this notebook is not runnable without a Microsoft Acedemic Graph API key)
 - Phenol_Explorer_Comparison.ipynb

# Files

#### Paper_Screening.ipynb
Notebook to search through the PubMed database and filter out search results.

#### Data_Statistics.ipynb
Notebook to compare clean raw data, compare databases, and visualize comparisons. This executes the majority of the analysis discussed in FoodMine. It begins by assigning keys to compounds, then calculates comparisons between databases, and finally generates paper visualizations.

#### Molecule_Embedding.ipynb
Notebook to retrive molecule smiles, embed compounds, and vizualize embeddings.

#### Paper_Citations.ipynb
Notebook to analyze the citation overlap between CTD and the papers gathered in FoodMine.

#### Phenol_Explorer_Comparison.ipynb
Notebook comparing information available on [Phenol Explorer](http://phenol-explorer.eu/) to data collected in FoodMine.

#### src
Miscellaneous functions and classes to serve as helpers in notebooks. 

* **pubmed_util.py**
Holds functions to interact with PubMed API for the purposes of our research.

* **filter.py**
Contains class Filter to filter out PubMed search results.

* **collected_data_handling.py**
Helper functions to clean and aggregate information into FoodMine database.

* **tools/chemidr**
Suite of tools developed to interact with the PubChem API.

* **file_downloader.py**
Scripts to download data folders hosted on Google Drive (large file sites prevents storing the files on Git).

#### data
Folder that holds raw data from paper data collection. **Note that several files used in our analysis are too large to be uploaded to github**

Select data files:

* (garlic/cocoa)_data.csv
Files that contains extracted information such as chemcials, PMID's, origional units, and sample metadata

* FoodMine_Output/fm_(garlic/cocoa).pkl
Files that contains the proccesed data that constitutes the FoodMine final data stores

* FoodMine_Output/compound_names_(garlic/cocoa).pkl
Files that contains processed compound names from literature and it's designated id in the FoodMine process