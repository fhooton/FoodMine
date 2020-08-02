import pandas as pd
import numpy as np
import os
import time

import config
import src.tools.chemidr.labeler as lbr
import src.tools.chemidr.id_map as id_map
import src.collected_data_handling as cdh

report = False

# Function to quickly record given statistics
def report_stat(text, filename, varname=None):
    if not os.path.exists('stats'):
        os.mkdir('stats')
        
    if varname is not None:
        text = text + '\n\tVar: ' + varname
    
    text = text + '\t' + time.strftime("%m/%d/%Y", time.localtime())
    
    with open(config.mfp('stats/' + filename), 'w') as f:
        f.write(text)


# Problem inputing letters into csv, so created system to convert them here
def greek_letter_converter(chem, convert_letter = True):
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

# Clean terms for various file applications
def clean_term(term, convert_letter = True, w_space = True, is_url=True):
    term = term.lower().strip()
    
    if convert_letter:
        term = greek_letter_converter(term)
    else:
        term = greek_letter_converter(term, convert_letter=False)
    
    if w_space:
        if is_url:
            term = term.replace(' ', '%20') # To replace ' ' in request
        else:
            pass
    else:
        term = term.replace(' ', '')
    return term

def id_loader(df, chem_key, load, file, fdb=True, pubchem=True):
    
    if load:
        df = pd.read_pickle(config.mfp(f'data/{file}'))
    else:
        df = lbr.id_searcher(df, chem_key, fdb=fdb, pubchem=pubchem)
        df.to_pickle(config.mfp(f'misc_save/{file}'))
    
    df.rename(columns={'pubchem_id' : 'chem_id_p', 'foodb_id' : 'chem_id_f'}, inplace=True)
    
    return df

def load_raw_data(food, load):
    food_data = pd.read_csv(config.mfp(f'data/{food}_data.csv'), encoding='latin1')
    
    food_scoring = pd.read_csv(config.mfp(f'data/{food}_scoring.csv'), encoding='latin1')

    # Need to remove phenol explorer ids that were manually put into data (for garlic only)
    food_data = food_data[food_data.PMID.isin(food_scoring.PMID.tolist())]

    food_data.chemical = food_data.chemical.str.lower()
    food_data.amount = food_data.amount.str.replace(',', '')

    food_data = food_data.merge(food_scoring[['PMID','is_useful']], how = 'left', on = 'PMID')
    
    if report:
        report_stat(f'Number of papers in search {food}: ' + str(len(food_scoring)), f'num_papers_srch_{food}.txt')
        report_stat(f'Number of papers reviewed {food}: ' + str(len(food_scoring[food_scoring.is_useful.notnull()])), f'num_reviewed_papers_{food}.txt')
        report_stat(f'Number of unique papers {food}: ' + str(len(food_data['PMID'].drop_duplicates())), f'num_unique_papers_{food}.txt')
        report_stat(f'Total number of records {food}: ' + str(len(food_data)), f'num_records_{food}.txt')

    return food_data, food_scoring

def append_keys_raw_data(food_data, food, load):
    if load:
        food_data = pd.read_pickle(config.mfp(f'data/{food}_food_data.pkl'))
    else:
        food_data = id_loader(food_data, 'chemical', load, config.mfp(f'{food}_food_data.pkl'))
    
    return food_data


# Make units processable
# Not all units were entered cleanly, so this function ensures that the data works for collected_data_handling
def unit_clean(df):
    df.units = df.units.replace('ug/ml', 'ug/L')
    df.units = df.units.replace('U/ml', 'ug/L')
    df.units = df.units.replace('mg/ml', 'ug/L')
    df.units = df.units.replace('mg/l', 'ug/L')
    df.units = df.units.replace('g/l', 'ug/L')
    df.units = df.units.replace('mM', 'ug/L')
    df.units = df.units.replace('mg/ ECE g', 'ug/L')
    df.units = df.units.replace('m/z', '%')
    df.units = df.units.replace('mg/lg', 'ug/L')
    df.units = df.units.replace('mmol/g', 'ug/L')
    df.units = df.units.replace('-', 'ug/L')

    df.amount = df.amount.str.strip('~')
    df.amount = df.amount.str.strip('<')

    # Replace with value of 0 when experiments measured quantified contents, but had no detection 
    df.amount = df.amount.replace('not quantified', '0')
    df.amount = df.amount.replace('LOD', '0')
    df.amount = df.amount.replace('ND', '0')
    df.amount = df.amount.replace('NQ', '0')
    df.amount = df.amount.replace('nd', '0')
    df.amount = df.amount.replace('n.q', '0')
    df.amount = df.amount.replace('n.q.', '0')
    df.amount = df.amount.replace('no analyzed', '0')
    df.amount = df.amount.replace('-', '0')
    
    return df

def clean_raw_data_strings(food_data):
    food_data.chemical = food_data.chemical.apply(clean_term, is_url=False)
    food_data = unit_clean(food_data)
    
    return food_data


# Partition data into quantified and unquantified
def partition_raw_data(food_data, food_scoring):

    for idx, row in food_data.iterrows():
        try:
            row['units'].count('g')
        except:
            print(row['units'])

        if row['units'].count('g') > 1:
            food_data.at[idx, 'is_quant'] = 1
        else:
            food_data.at[idx, 'is_quant'] = 0

    # Remove all rows where the units are %'s
    food_data_q = food_data[food_data['is_quant'] == 1].reset_index(drop=True)

    # Have a seperate dataframe for all chemicals that we would put in the category of 'detected but not quantified'
    food_data_dnq = food_data[food_data['is_quant'] == 0].reset_index(drop=True)

     # The quantified dataframe for values that are both quantified and unquantified
    unq_chems = list(set( food_data_dnq['chemical'].tolist() ))
    food_data_both = food_data_q.iloc[[idx for idx, row in food_data_q.fillna('placeholder').iterrows() if row['chemical'] in unq_chems]]

    # Remove occurances of overlaping chemicals from the unquantified garlic data
    q_chems = list(set( food_data_q['chemical'].tolist() ))
    food_data_dnq = food_data_dnq.iloc[[idx for idx, row in food_data_dnq.fillna('placeholder').iterrows() if row['chemical'] not in q_chems]]
    
    return food_data_q, food_data_dnq


# Creates food_mine database data from raw collected data
def build_food_mine(food_data, food_data_q, food_data_dnq):
    food_mine = cdh.build_data_dict(food_data)

    quant_food_mine = cdh.build_data_dict(food_data_q)

    unquant_food_mine = cdh.build_data_dict(food_data_dnq)
    
    # Need to recompare quantified chems and unquantified chems with synonym key to do one last removal
    q_chems = list(set( quant_food_mine['chem_id'].dropna().tolist() ))
    unquant_food_mine = unquant_food_mine[~unquant_food_mine.chem_id.isin(q_chems)].reset_index()
    
    if report:
        report_stat(f'FM size {food}: ' + str(len(food_mine)), f'fm_size_{food}.txt')
        report_stat(f'QFM size {food}: ' + str(len(quant_food_mine)), f'qfm_size_{food}.txt')
        report_stat(f'UQFM size {food}: ' + str(len(unquant_food_mine)), f'uqfm_size_{food}.txt')

    return food_mine, quant_food_mine, unquant_food_mine


# Loads data from FooDB
def load_foodb_data(food, load):
    # Dataframe with contents of foodb
    
    if not load:
        foodb = pd.read_csv(config.mfp('data/contentssql.csv'))
        foodb = foodb[(foodb.source_type != 'Nutrient') & (foodb.source_id != 0) & (foodb.standard_content != 0)]

        compounds = pd.read_csv(config.mfp('data/compounds.csv'), encoding='latin1')

        foodb = foodb.merge(compounds[['id', 'name']], how='left', left_on='source_id', right_on='id')

        if food == 'garlic':
            # Garlic - ["Garlic", "Soft-necked Garlic"]
            target_foodb_food_id = [8, 880]

        if food == 'cocoa':
            # Cocoa - ["cocoa bean", "cocoa butter", "Cocoa powder", "Cocoa Liquor"]
            target_foodb_food_id = [182, 706, 707,708]

        # Gets the subset of the database pertaining to food
        foodb_food = foodb[foodb.food_id.isin(target_foodb_food_id)].reset_index(drop=True)

        # Transforms all the chemical names to lowercase for syncing
        foodb_food.name = foodb_food.name.str.lower()

        foodb_food = foodb_food.rename(index=str, columns={"source_id": "foodb_id"})
    
    if load:
        foodb_food = pd.read_pickle(config.mfp(f'data/{food}_foodb_food.pkl'))
        foodb_food.rename(columns={'orig_source_name' : 'name'}, inplace=True)
                
    foodb_food = id_loader(foodb_food, 'name', load, f'{food}_foodb_food.pkl',fdb=False)

    # Creates a list of the unique chemicals in garlic from foodb
    foodb_food_lower = list(set( foodb_food.chem_id.tolist() ))

    # Creates a seperate dataframe that holds chemicals for garlic in foodb with a real quantification
    quant_foodb_food = foodb_food[foodb_food.standard_content.notnull()][['chem_id', 'chem_id_f', 'orig_source_id','name', 'standard_content']].drop_duplicates()

    # Creates a seperate dataframe that holds chemicals for garlic in foodb without a real quantification
    unquant_foodb_food = foodb_food[foodb_food.standard_content.isnull()][['chem_id', 'chem_id_f', 'orig_source_id','name', 'standard_content']].reset_index()
    
    q_ids = list(set( quant_foodb_food.chem_id.tolist() ))
    q_names = list(set( quant_foodb_food.chem_id_f.tolist() ))
    unquant_foodb_food = unquant_foodb_food[(~unquant_foodb_food.chem_id.fillna('-').isin(q_ids))
                                           & (~unquant_foodb_food.chem_id_f.fillna('-').isin(q_names))]
    
    if report:
        report_stat(f'FDB size {food}: ' + str(len(foodb_food.chem_id.drop_duplicates())), f'fdb_size_{food}.txt')
        report_stat(f'QFDB size {food}: ' + str(len(quant_foodb_food.chem_id.drop_duplicates())), f'qfdb_size_{food}.txt')
        report_stat(f'UQFDB size {food}: ' + str(len(unquant_foodb_food.chem_id.drop_duplicates())), f'uqfdb_size_{food}.txt')

    return foodb_food, quant_foodb_food, unquant_foodb_food


###### Loads USDA data
def load_usda_data(food, load):
    if not load:
        if food == 'garlic':
            # Garlic, 'Garlic, raw', 'Spices, garlic powder'
            NDB_id = [11215, 2020]

        if food == 'cocoa':
            # Cocoa, 'Oil, cocoa butter', 'Cocoa, dry powder, Hershey's European style cocoa', 
            # 'Cocoa, dry powder, unsweetened', 'Cocoa, dry powder, unsweetend, processed with alkali',
            # 'Cocoa, dry powder, hi-fat or breakfast, processed with alkali'
            NDB_id = [4501, 19171, 19165, 19166, 19860]

        # Reads in USDA database
        usda = pd.read_csv(config.mfp('data/SR28_plus_flav.csv'), encoding = 'latin1')

        # Filters out rows not apart of NDB_id
        usda = usda[usda.NDB_No.isin(NDB_id)][['NDB_No','food_name', 'Nutr_No_new', 'nut_desc', 'Nutr_Val', 'unit']]
        usda['num_measures'] = 1

        # Average chemicals that appear in multiple USDA food catagoriess
        for nutr in usda.Nutr_No_new.drop_duplicates().tolist():
            temp = usda[usda.Nutr_No_new == nutr]
            if len(temp) > 1:
                if len(temp.unit.drop_duplicates()) > 1:
                    print(nutr, 'has different units for same nutrient')
                new_row = temp.copy().reset_index(drop=True).loc[0,:]
                new_row['Nutr_Val'] = temp.Nutr_Val.mean()
                new_row['num_measures'] = len(temp)

                usda = usda.drop(temp.index)
                usda = usda.append(new_row)

        usda = usda.reset_index(drop=True)

    # Append chemical key matcher to USDA chemicals
    if load:
        usda = pd.read_pickle(config.mfp(f'data/{food}_usda.pkl'))
    else:
        usda = id_loader(usda, 'nut_desc', load, f'{food}_usda.pkl').reset_index(drop=True)
        
    usda.rename(columns={'foodb_id' : 'chem_id_f'}, inplace=True)
        
    usda = usda[~usda.unit.isin(['IU', 'kcal', 'kJ'])].reset_index(drop=True)

    if report:
        report_stat(f'USDA size {food}: ' + str(len(usda)), f'usda_size_{food}.txt')
        
    return usda


def load_ctd():
    skip = list(range(26)) # First few lines are empty / not useful info
    hdata = pd.read_csv('data/CTD_chemicals_diseases.csv', skiprows=skip).reset_index()
    hdata.columns = ['ChemicalName', 'ChemicalID', 'CasRN', 'DiseaseName', 'DiseaseID', 'DirectEvidence', 'InferenceGeneSymbol', 'InferenceScore', 'OmimIDs', 'PubMedIDs']
    hdata = hdata.drop([0,1], axis = 0).reset_index(drop=True)

    health_pubchem_ids = pd.read_pickle('misc_save/health_chem_pubchem_ids.pickle').drop('pubchem_name', axis=1)

    hdata = hdata.merge(health_pubchem_ids, how='left', on='ChemicalName')

    return hdata



def load_health():
    hdata = load_ctd()

    # Count the number of 'Direct Evidence' listings per chemical with a pubchem id
    de_health = pd.DataFrame(
        hdata[hdata.pubchem_id.notnull() & hdata.DirectEvidence.notnull()][['pubchem_id', 'DirectEvidence', 'ChemicalName']]
        .groupby(['pubchem_id','DirectEvidence']).count()).reset_index()
    
    # The 'ChemicalName' now holds the count of 'Direct Evidence' listings
    return de_health[['pubchem_id', 'ChemicalName']].groupby('pubchem_id').sum()