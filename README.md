![Logo](https://github.com/fhooton/FoodMine/blob/master/images/NetSci_Logo.png)

# FoodMine

FoodMine is a food database project from the [Barabasi Lab](https://www.barabasilab.com/) of the [Network Science Institute](https://www.networkscienceinstitute.org/). The objective of this project is to mine research litereature to uncover chemical contents in food. The code in this repository was used to generate the food pilots, plots, and other information in [FoodMine: Exploring Food Contents in Scientific Literature](https://www.biorxiv.org/content/10.1101/2019.12.17.880062v1) (Hooton, [Menichetti](https://www.barabasilab.com/people/giulia-menichetti), [Barabási](barabasilab.com/people/laszlo-barabasi)). Additionally, we have included the data from out initial experiments.

> **FoodMine Abstract**
> 
> Thanks to the many chemical and nutritional components it carries, diet critically affects human health. However, the currently available comprehensive databases on food composition cover only a tiny fraction of the total number of chemicals present in our food, focusing on the nutritional components essential for our health. Indeed, thousands of other molecules, many of which have well documented health implications, remain untracked. To explore the body of knowledge available on food composition, we built FoodMine, an algorithm that uses natural language processing to identify papers from PubMed that potentially report on the chemical composition of garlic and cocoa. After extracting from each paper information on the reported quantities of chemicals, we find that the scientific literature carries extensive information on the detailed chemical components of food that is currently not integrated in databases. Finally, we use unsupervised machine learning to create chemical embeddings, finding that the chemicals identified by FoodMine tend to have direct health relevance, reflecting the scientific community’s focus on health-related chemicals in our food.

> ![F3](https://github.com/fhooton/FoodMine/blob/master/images/database-comp.png)
> **Figure 2:** Number of Unique Compounds Recovered by FoodMine, USDA, and FooDB. The plots show the number of unique compounds reported by USDA, FooDB, and FoodMine. The columns display 1) the total number of unique quantified compounds in each database, 2) the total number of unique unquantified compounds in each database, and 3) the number of quantified compounds retrieved by FoodMine and never reported before in USDA or FooDB.

> ![F5](https://github.com/fhooton/FoodMine/blob/master/images/cocoa_TSNE_replot.png)
> **Figure 5:** TSNE Dimensionality Reduction of Chemical Embeddings with Health Associations. TSNE plots of Mol2Vec chemical embeddings for garlic (A, B, and C) and cocoa (D, E, and F). The colors of each data point encode the number of health implications associated with the compounds based on the CTD database. Dark gray represents chemicals with 0 health associations. We show chemicals catalogued by each studied database for FoodMine (A & D), USDA (B & E), and FooDB (C & F). Markers are filled if the database contains the chemical, and empty if it does not.

# Setup

Getting started using a conda env:

```shell
conda create --name foodmine python=3.6
conda activate foodmine
```

Installing packages:

```shell
pip install -r requirements.txt
./install_additional_packages.sh
```

Run the [python notebooks](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/) in the following order to replicate the FoodMine process.
 - Paper_Screening.ipynb
 - Data_statistics.ipynb (this notebook is key to running other downstream files as it downloads relevant data hosted online)
 - Molecule_Embedding.ipynb
 - Paper_Citations.ipynb (this notebook is not runnable without a Microsoft Acedemic Graph API key)
 - Phenol_Explorer_Comparison.ipynb

# Files

#### Paper_Screening.ipynb
Notebook to search through the PubMed database and filter out search results.

#### Data_Statistics.ipynb
Notebook to compare clean raw data, compare databases, and visualize comparisons. This executes the majority of the analysis discussed in FoodMine. It begins by assigning keys to compounds, then calculates comparisons between databases, and finally generates paper visualizations.

#### Molecule_Embedding.ipynb
Notebook to retrieve molecule smiles, embed compounds, and visualize embeddings.

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
Folder that holds raw data from paper data collection.

Select data files:

* data/(garlic/cocoa)_data.csv
Files that contains extracted information such as chemcials, PMID's, origional units, and sample metadata

* FoodMine_Output/fm_(garlic/cocoa).pkl
Files that contains the processed data that constitutes the FoodMine final data stores

* FoodMine_Output/compound_names_(garlic/cocoa).pkl
Files that contains processed compound names from literature and it's designated id in the FoodMine process

**Note: some of the matching procedures leveraged undisclosed datasets under development for future projects. For more info contact Giulia Menichetti at g.menichetti@northeastern.edu.**