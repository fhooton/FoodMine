# import requests

# def download_file_from_google_drive(id, destination):
#     URL = "https://drive.google.com/uc?export=download"

#     session = requests.Session()

#     response = session.get(URL, params = { 'id' : id }, stream = True)
#     token = get_confirm_token(response)

#     if token:
#         params = { 'id' : id, 'confirm' : token }
#         response = session.get(URL, params = params, stream = True)

#     save_response_content(response, destination)    

# def get_confirm_token(response):
#     for key, value in response.cookies.items():
#         if key.startswith('download_warning'):
#             return value

#     return None

# def save_response_content(response, destination):
#     CHUNK_SIZE = 32768

#     with open(destination, "wb") as f:
#         for chunk in response.iter_content(CHUNK_SIZE):
#             if chunk: # filter out keep-alive new chunks
#                 f.write(chunk)

# def download_data():
#     file_id = '1qUSKg-KepGD50L9HLtdMQM6PIP9N2qwQ'
#     destination = r"C:\Users\forresthooton\Documents\NetSci Research\Foodome Project\FoodMine\src\misc_save"
#     download_file_from_google_drive(file_id, destination)

# if __name__ == "__main__":
#     file_id = '1qUSKg-KepGD50L9HLtdMQM6PIP9N2qwQ'
#     destination = 'misc_save'
#     download_file_from_google_drive(file_id, destination)

# https://drive.google.com/drive/folders/1qUSKg-KepGD50L9HLtdMQM6PIP9N2qwQ?usp=sharing
# https://drive.google.com/file/d/1Enc3FOXDb8R2gGGVnn73FOQEAYqePn9M/view?usp=sharing

import gdown
# importing required modules 
from zipfile import ZipFile

def unzip_file(dirzip):
    with ZipFile(dirzip, 'r') as zip: 
        zip.extractall() 

def download_core_data():
    url = 'https://drive.google.com/uc?id=1urQmygS0vrbWN_4nU6R44q85XkPomYLD'
    output = r"C:\Users\forresthooton\Documents\NetSci Research\Foodome Project\FoodMine\src\misc_save.zip"
    gdown.download(url, output, quiet=False)

    unzip_file(output)

def download_intermediate_data():
    url = 'https://drive.google.com/uc?id=1Enc3FOXDb8R2gGGVnn73FOQEAYqePn9M'
    output = r"C:\Users\forresthooton\Documents\NetSci Research\Foodome Project\FoodMine\src\misc_save.zip"
    gdown.download(url, output, quiet=False)

    unzip_file(output)

def download_all_data():
    download_core_data()
    download_intermediate_data()