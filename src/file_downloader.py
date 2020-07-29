import os
from shutil import copyfile
import gdown
from zipfile import ZipFile

cwd = os.getcwd()

def unzip_file(dirzip):
    with ZipFile(dirzip, 'r') as zip: 
        zip.extractall() 

def download_core_data():
    url = 'https://drive.google.com/uc?id=1urQmygS0vrbWN_4nU6R44q85XkPomYLD'
    output = os.path.join(cwd, "data.zip")
    gdown.download(url, output, quiet=False)

    unzip_file(output)

def download_intermediate_data():
    url = 'https://drive.google.com/uc?id=1Enc3FOXDb8R2gGGVnn73FOQEAYqePn9M'
    output = os.path.join(cwd, "misc_save.zip")
    gdown.download(url, output, quiet=False)

    unzip_file(output)

def download_chemidr_intermediate_data():
    url = 'https://drive.google.com/uc?id=1vl5TSitI30bT5V_wkFfil8mqvQzvFhNn'
    output = os.path.join(cwd, "src/tools/intermediate_data.zip")
    
    gdown.download(url, output, quiet=False)

    with ZipFile(output, 'r') as zip: 
        zip.extractall(path = os.path.join(cwd, "src/tools"))

def make_copypath_and_copy(src, dst):
    dst = f'{dst}/{src}'
    src = f'data/{src}'

    print(src,dst)
    copyfile(src, dst)


def copy_chemidr_files():
    os.mkdir('src/tools/data')
    files = ['compounds.csv', 'compound_synonymssql.csv', 'contentssql.csv', 'usda_raw_garlic.csv']

    for f in files:
        make_copypath_and_copy(f, 'src/tools/data')


def download_all_data():
    download_core_data()
    download_intermediate_data()
    download_chemidr_intermediate_data()
    copy_chemidr_files()