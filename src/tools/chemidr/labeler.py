"""
    Author: Forrest Hooton
    Purpose: Label given chemical strings with corresponding id's in popular databases (pubchem and foodb)

    Primary Functions:
        get_compound_pubchem_info: returns pubchem id and name given a string
        append_foodb_id: adds column with foodb id's given a dataframe with chemical strings
        append_pubchem_id: adds column with pubchem ids's given a dataframe with chemical strings
        id_searcher: adds columns with foodb and pubchem id's, along with a composite column

    Example Usage:
        test = pd.DataFrame({'chem' : ['allicin', 'diallyl disulphide']})
        id_searcher(test, 'chem')
"""

import pandas as pd
import numpy as np
import time
import requests
from lxml import etree
import json
import math
import pickle
import re
from fuzzywuzzy import fuzz

from .id_map import cids2inchis

import os
cwd = os.path.dirname(__file__) # get current location of script
package_path = 'tools'.join(cwd.split('tools')[:-1]) + 'tools' # get path to head of Chemidr package


# Filepath wrapper to make data load-able from Chemidr
def __make_fp__(fp):
    return f'{package_path}/{fp}'
    

# Convert text to greet letter using *~greek letter name~*
def __greek_letter_converter__(chem, convert_letter = True):
    if convert_letter:
        chem = chem.replace('*alpha*', 'α')
        chem = chem.replace('*beta*', 'β')
        chem = chem.replace('*gamma*', 'γ')
        chem = chem.replace('*rho*', 'ρ')
        chem = chem.replace('*delta*', 'δ')
    else:
        chem = chem.replace('*alpha*', 'alpha')
        chem = chem.replace('*beta*', 'beta')
        chem = chem.replace('*gamma*', 'gamma')
        chem = chem.replace('*rho*', 'rho')
        chem = chem.replace('*delta*', 'delta')
    return chem


def __clean_term__(term, convert_letter = True, w_space = True, is_url=True):
    """
        Prepares an input term to be queried in a url

        Input
        ----------------------------------------------------------------
        term : str
            term to clean
        convert_letter : bool (default True)
            Whether or not to convert letter to greek representation
        w_space : bool (default True)
            keep space in term, removes spaces if false
        is_url : bool (default True)
            Replaces spaces with %20 if true

        Returns
        ----------------------------------------------------------------
        term : str
            same term as input after cleaning process

    """
    term = term.lower().strip()
    
    term = __greek_letter_converter__(term, convert_letter=convert_letter)
    
    # Keeps spaces in string
    if w_space:
        if is_url:
            term = term.replace(' ', '%20') # To replace ' ' in request
        else:
            pass
    else:
        term = term.replace(' ', '')
    return term


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
    response = requests.get(url)

    if response.status_code == 200: # Successful
        return response.content

    elif response.status_code == 429: # Too many requests
        time.sleep(.5)
        return __safe_urlopen__(url)

    elif response.status_code == 503: # PUGREST.ServerBusy
        time.sleep(1)
        return __safe_urlopen__(url)

    elif response.status_code == 404: # PUGREST.NotFound (aka doesn't exist)
        return None


def __exact_retrieval__(req):
    """
        retrieves pubchem synonym information using exact string matches and extracts the 
        compound id and primary compound name

        Input
        ----------------------------------------------------------------
        req : str
            chemical string to query in pubmed

        Returns
        ----------------------------------------------------------------
        compound_id : float
            pubchem compound id of req chemical
        compound_name : str
            compound name of primary compound for synonym entry
    """
    # Ex. https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/allicin/synonyms/XML
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{req}/synonyms/XML"

    xml = __safe_urlopen__(url)

    if xml is None:
        return np.nan, np.nan
    
    # Extracts information from xml string
    root = etree.fromstring(xml)
    compound_id = float(root.findall(".//{http://pubchem.ncbi.nlm.nih.gov/pug_rest}CID")[0].xpath('.//text()')[0])
    compound_name = str(root.findall(".//{http://pubchem.ncbi.nlm.nih.gov/pug_rest}Synonym")[0].xpath('.//text()')[0])
    
    return compound_id, compound_name


def __is_nutrient__(term):
    general_terms = ['protein', 'fat', 'sugar', 'carbohydrate', 'fiber', 'ash']
    
    if sum([1 for t in general_terms if t in term]) > 0:
        return True
    else:
        return False


def __clean_compound_name__(s):
    s = s.lower()

    # Remove +/- from string
    s = re.sub('[\+\-]', '', s)

    # Remove trailing s in string
    s = re.sub('s$', '', s)

    # Remove cis/trans from string
    s = re.sub('[\s\-\(]cis[\s\-\)]|[\s\-\(]trans[\s\-\)]', '', s)

    # Remove any single alpha characters surrounded by whitespace or '-'
    s = re.sub('[\s\-\(][a-z][\s\-\)]|^[a-z][\s\-]', '-', s)

    # Replace hyphen with space
    s = re.sub('-', ' ', s)

    # Replace stretch of whitespace with space
    s = re.sub('\s+', ' ', s)

    s = s.strip()

    return s


def __complex_string_equivalence__(s1, s2):
    # Remove parenthesis and anything in between
    s1 = __clean_compound_name__(s1)
    s2 = __clean_compound_name__(s2)

    # Calculates fuzz ratio (Levenshtein Distance)
    if fuzz.ratio(''.join(s1.split()), ''.join(s2.split())) >= 98:
        return True
    else:
        return False


# Examples https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term=allicin&retmode=json
# https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/synonyms/json
def __compound_search__(req):
    """
        Searches PubChem database with search term (req) to retrieve information for in-exact matches

        Input
        ----------------------------------------------------------------
        req : str
            chemical string to query in pubmed

        Returns
        ----------------------------------------------------------------
        cid : float
            pubchem compound id of req chemical
        name : str
            compound name of primary compound for synonym entry
    """
    req = req.lower()

    # Url to search for compound
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term={req}&retmode=json"

    r = __safe_urlopen__(url)

    if r is None:
        return np.nan, np.nan

    req = re.sub('%20', ' ', req)

    # Return np.nan for no search results
    if int(json.loads(r)['esearchresult']['count']) == 0:
        return np.nan, np.nan

    cid = float(json.loads(r)['esearchresult']['idlist'][0])

    # Url to get the name associated with a cid
    name_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/json"

    # a cid did not have synonyms for some reason
    jsn = __safe_urlopen__(name_url)

    if jsn is None:
        return np.nan, np.nan

    name = json.loads(jsn)['InformationList']['Information'][0]['Synonym'][0]

    # Return the cid and name if the search retrieves 1 result and passes string equivalence
    if int(json.loads(r)['esearchresult']['count']) == 1:
        if __complex_string_equivalence__(req, name) == True:
            return cid, name
        else:
            return np.nan, np.nan

    # Returns cid and name if first result passes the complex equivalence test with the req
    elif int(json.loads(r)['esearchresult']['count']) > 1:
        if __complex_string_equivalence__(req, name) == True:
            return cid, name
        else:
            return np.nan, np.nan


# try-except wrapper to retrieve pubchem info
def get_compound_pubchem_info(chem):
    if __is_nutrient__(chem):
        return np.nan, np.nan
        
    req = __clean_term__(chem, convert_letter=False)
    compound_id, compound_name = __exact_retrieval__(req)
    
    if not math.isnan(compound_id):
        return compound_id, compound_name

    req = __clean_term__(chem, convert_letter=False, w_space=False)
    compound_id, compound_name = __exact_retrieval__(req)

    if not math.isnan(compound_id):
        return compound_id, compound_name

    req = __clean_term__(chem, convert_letter=False)
    compound_id, compound_name = __compound_search__(req)
    
    if not math.isnan(compound_id):
        return compound_id, compound_name

    req = __clean_term__(__clean_compound_name__(chem), convert_letter=False)
    compound_id, compound_name = __compound_search__(req)

    return compound_id, compound_name


# Evaluates if a cell already contains a key. If not, and if there is a match, insert a key
def __check_key__(row, id_col, str_col, input_dict):
    if np.isnan(row[id_col]):
        if row[str_col] in input_dict:
            return input_dict[row[str_col]]
        else:
            return np.nan
    else:
        return row[id_col]
    

def append_foodb_id(df, chem_key, load_ids=True):
    """
        label chemicals with foodb id where there is a match to chemicals in foodb files under the data directory

        Input
        ----------------------------------------------------------------
        df : pd.DataFrame
            dataframe containing chemical names
        chem_key : str
            column name with the chemical strings
        load_ids : bool (default True)
            Whether or not to use pre-altered pubchem data. Want to use true after the program has run once to reduce runtime

        Returns
        ----------------------------------------------------------------
        df : pd.DataFrame
            input dataframe with a foodb_id column
    """
    df[chem_key] = df[chem_key].str.strip().str.lower()
    
    # Read in foodb primary compound names
    fdb_compounds = pd.read_csv(__make_fp__('data/compounds.csv'), encoding='latin1')[['id', 'name']]
    fdb_compounds.name = fdb_compounds.name.str.strip().str.lower()
    fdb_compounds.rename(columns={'id' : 'foodb_id'}, inplace=True)
    
    # Merge on matching names
    df = df.merge(fdb_compounds, how = 'left', left_on = chem_key, right_on = 'name')
    
    if load_ids:
        with open(__make_fp__('intermediate_save/fdb_source_strings.pkl'), 'rb') as f:
            input_dict = pickle.load(f)
    else:
        # Dataframe with contents of foodb
        foodb = pd.read_csv(__make_fp__('data/contentssql.csv'))

        # Gets the subset of the database pertaining to garlic
        foodb = foodb[['source_id', 'orig_source_name']].drop_duplicates()

        # Transforms all the chemical names to lowercase for syncing
        foodb.orig_source_name = foodb.orig_source_name.str.strip().str.lower()

        # Creates a list of the unique chemicals in garlic from foodb
        input_dict = {row['orig_source_name'] : row['source_id'] for _, row in foodb.iterrows()}
        
        with open(__make_fp__('intermediate_save/fdb_source_strings.pkl'), 'wb') as f:
            pickle.dump(input_dict, f)
    
    df.foodb_id = df.apply(__check_key__, id_col='foodb_id', str_col=chem_key, input_dict=input_dict, axis=1)
    
    if load_ids:
        with open(__make_fp__('intermediate_save/usda_ids.pkl'), 'rb') as f:
            input_dict = pickle.load(f)
    else:
        usda = pd.read_csv(__make_fp__('data/usda_raw_garlic.csv'), encoding = 'latin1')
        usda.nut_desc = usda.nut_desc.str.strip().str.lower()
        input_dict = {row['nut_desc'] : row['chem_id'] for _, row in usda.iterrows()}
        
        with open(__make_fp__('intermediate_save/usda_ids.pkl'), 'wb') as f:
            pickle.dump(input_dict, f)
    
    df.foodb_id = df.apply(__check_key__, id_col='foodb_id', str_col=chem_key, input_dict=input_dict, axis=1)
                
    return df


def append_pubchem_id(df, chem_key):
    """
        label chemicals with pubchem id by searching the pubchem synonyms

        Input
        ----------------------------------------------------------------
        df : pd.DataFrame
            dataframe containing chemical names
        chem_key : str
            column name with the chemical strings

        Returns
        ----------------------------------------------------------------
        df : pd.DataFrame
            input dataframe with a pubchem_id column
    """
    start = time.time()

    df['pubchem_id'] = np.nan
    df['pubchem_name'] = np.nan

    df['pubchem_id'] = df['pubchem_id'].astype(object)
    df['pubchem_name'] = df['pubchem_name'].astype(object)

    i = 0
    for idx, row in df.iterrows():
        ID, name = get_compound_pubchem_info(row[chem_key])
        df.at[idx, 'pubchem_id'] = ID
        df.at[idx, 'pubchem_name'] = name
        
        if not i % 1000:
            print(idx, 'chems searched in', (time.time() - start) / 60, "min")

        i += 1

    print("Pubchem ids added in", (time.time() - start) / 60, "min")

    return df


def id_searcher(df, chem_key, fdb = True, pubchem = True, use_prefix=True):
    """
        main function to assign chemical keys in pubchem and foodb

        Input
        ----------------------------------------------------------------
        df : pd.DataFrame
            dataframe containing chemical names
        chem_key : str
            column name with the chemical strings
        fdb : bool (default True)
            assign foodb ids
        pubchem : bool (default True)
            assign pubchem ids
        use_prefix : bool (default True)
            only return the prefix of inchikeys (before the first -), which contains the structural information
            (to find out more see https://www.inchi-trust.org/technical-faq-2/)

        Returns
        ----------------------------------------------------------------
        df : pd.DataFrame
            input dataframe with a pubchem_id + foodb_id columns, and composite column chem_id
    """
    if pubchem:
        df = append_pubchem_id(df, chem_key)
        print(len(df))
    if fdb:
        df = append_foodb_id(df, chem_key)
        print(len(df))
    
    total = len(df[chem_key].drop_duplicates())
    
    df.at[df[df.pubchem_id.notnull()].index, 'inchikey'] = cids2inchis(df[df.pubchem_id.notnull()].pubchem_id.tolist(), use_prefix=use_prefix)
    
    # Manually looked up the maximum pubchem index to make sure id's don't overlap
    max_p_index = 134825000

    # Ensure that column type is object
    df['chem_id'] = np.nan
    df['chem_id'] = df['chem_id'].astype(object)
    
    # If there is no pubchem_id, make chem_id the foodb_id + the maximum pubchem id
    for idx, row in df.iterrows():
        if isinstance(row['inchikey'], str):
            df.at[idx, 'chem_id'] = row['inchikey']
        else:
            df.at[idx, 'chem_id'] = row['foodb_id'] + max_p_index
    
    return df
