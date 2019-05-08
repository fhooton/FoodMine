import pandas as pd
import numpy as np
import re
from unit_converter.converter import convert, converts


# Extracts digits from start of string
def __extract_digits__(string):
    i = 0
    
    while string[i].isdigit():
        i += 1
    
    return int(string[0:i])


# recieves a unit in string form, even containing a number, and returns scales and units
def __splice_unit__(unit):
    
    unit = unit.strip()

    if unit[0].isdigit():
            scale = __extract_digits__(unit) # I.e. 10 in 10mg
            num_len = sum(c.isdigit() for c in str(scale)) # gets string length of number, i.e. 100 is 3

            units = unit[num_len:]	# I.e. mg in 10mg
    else:
        scale = 1
        units = unit

    return scale, units


def __converter__(unit, target):
    
    # Package doesn't like u for nano
    if unit[0] == 'u':
        unit = 'µg'

    if target[0] == 'u':
    	target = 'µg'

    unit = '1 ' + unit

    return float(converts(unit, target))



# Handles unit conversions to specified quantity, and makes adjustments as necessary
def __unit_handler__(value, unit, target_unit):

    if re.search(r'to', str(value)) is not None:
        value = np.mean([float(n) for n in value.split(' to ')])

    if unit == target_unit:
        return float(value)

    value = float(value)

    # Handle if unit is percentage
    if unit == '%':
        target_denom_scale, target_denom_unit = __splice_unit__(target_unit.split('/')[1])
        unit_conversion = __converter__(target_unit.split('/')[1], target_unit.split('/')[0])
        converted_value = float(value) * target_denom_scale * unit_conversion / 100

        #return converted_value
        return '% measurement'

    # Handle if unit is based on volume rather than mass
    # Some units give errors; try and except prints those units
    try:
        if unit.split('/')[1][0] == 'L':
            return 'measured in liters'
    except:
        return '% measurement'

    if unit.count('/') is not target_unit.count('/'):
        print("Unit conversion error, one is not fraction")
        return
    
    if unit.count('/') > 0:

        # Check if the first charechter of digit is string, assumed conversion value if it is
        num_scale, num_unit = __splice_unit__(unit.split('/')[0])
        denom_scale, denom_unit = __splice_unit__(unit.split('/')[1])

        target_num_scale, target_num_unit = __splice_unit__(target_unit.split('/')[0])
        target_denom_scale, target_denom_unit = __splice_unit__(target_unit.split('/')[1])

        # Handles the scale conversion not relating to units
        if denom_scale / target_denom_scale is not 1:
            if num_scale / target_num_scale is not 1:
                scale_conversion = (num_scale / target_num_scale) / (denom_scale / target_denom_scale)
            else:
                scale_conversion = 1 / (denom_scale / target_denom_scale)
        else:
            scale_conversion = num_scale / target_num_scal


        # Recieves numerical scales from converter. For __converter__('mg', 'g'), returns 1000
        num_unit_conversion = __converter__(num_unit, target_num_unit)
        denom_unit_conversion = __converter__(denom_unit, target_denom_unit)

        unit_conversion = num_unit_conversion / denom_unit_conversion

        value_conversion = value * scale_conversion * unit_conversion

        return value_conversion
    

# Recieves DataFrame with at least a numerical 'amount' column and calculates summary values to return
def __quant_handler__(df):
    
    target_unit = 'mg/100g'

    absolute_min = None
    absolute_max = None
    
    means =[]
    
    # Calculate avg mean
    for paper_id in df.PMID.drop_duplicates().tolist():
        
        paper_values = []
        
        for _, row in df[df['PMID'] == paper_id].iterrows():
            
            #try:
            conversion = __unit_handler__(row['amount'], row['units'], target_unit)	# Sets the units for the whole output
            #except:
            #    print('\'' + row['units'] + '\'', "broke the conversion")
            
            # If there is an error or issiue in __unit_handler__, will return some sort of string
            if type(conversion) is not str:
            
                # Keeps track of the global min and max for a chemical
                if absolute_min is None:
                    absolute_min = conversion
                    absolute_max = conversion
                else:
                    if absolute_min > conversion:
                        absolute_min = conversion
                    if absolute_max < conversion:
                        absolute_max = conversion

                paper_values.append(conversion)
            
        
        if len(paper_values) > 0:
            means.append(sum(paper_values) / len(paper_values))
    
    if len(means) is not 0:
        avg_mean = np.mean(means)
        min_mean = np.min(means)
        max_mean = np.max(means)
        mean_var = np.var(means)
        median = np.median(means)   
    else:
        avg_mean = None
        min_mean = None
        max_mean = None
        mean_var = None
        median = None 
    
    quant_pack = {
        'average_mean' : avg_mean,
        'absolute_min' : absolute_min,
        'absolute_max' : absolute_max,
        'min_mean' : min_mean,
        'max_mean' : max_mean,
        'mean_variance' : mean_var,
        'median' : median,
        'units' : target_unit
    }

    return quant_pack

# Recieves DataFrame of specific chemical and its various values
def __build_chem_subdict__(df):
    
    # Just take first chemical for general proxy
    chemical = df.chemical.tolist()[0]

    # Number of unique papers
    num_papers = len(df.PMID.drop_duplicates())
    papers = df.PMID.drop_duplicates().tolist()
    
    # Number of terms quantified
    num_terms_quant = df.amount.notnull().sum()
    
    # avg mean, total min, total max, min mean, max mean, mean variance
    quant_pack = __quant_handler__(df[df['amount'].notnull()])
    
    # total # samples
    total_num_samples = df.num_samples.notnull().sum()
    
    # Number without info on total # samples
    num_without_sample_info = df.num_samples.isnull().sum()
    
    # Number without quantity or sample info
    num_lacking_info = len(df[(df['num_samples'].isnull()) & (df['amount'].isnull())])
    
    if len(df[df.chem_id_p.notnull()]) > 0:
        compound_name = df['pubchem_name'].tolist()[0]
    else:
        compound_name = df.chemical.tolist()[0]

    chem_subdict = {
        'num_papers' : num_papers,
        'papers' : papers,
        'num_terms_quantified' : num_terms_quant,
        'chemical' : compound_name,
        'chem_id' : df.chem_id.tolist()[0],
        'chem_id_p' : df.chem_id_p.tolist()[0]
    }

    chem_subdict.update(quant_pack)

    return chem_subdict


# Recieves of DataFrame of compiled data and appropriately calculates/reformats data into final structure
# df must be standardized for project
def build_data_dict(df):
    
    data_dict = {}  # Dict that holds all foods and chemicals, and their calculations

    df['chemical'] = df['chemical'].str.strip()
    df['chemical'] = df['chemical'].apply(greek_letter_converter)
    
    # Iterate over each unique food, and each unique chemical per food to create food, chem dict
    for food in df['food'].drop_duplicates().tolist():
        
        chem_dict = {}   # Holds chemical calculations per food
        
        # need to do one loop for things with a chem id
        for chem_id in df[(df['food'] == food) & (df['chem_id'].notnull())].chem_id.drop_duplicates().tolist():

            #chem = greek_letter_converter(chem)
            
            f_chem_df = df[(df['food'] == food) & (df['chem_id'] == chem_id)]   # Filters DataFrame for rows with specific food-chemical combination
            
            new_entry = __build_chem_subdict__(f_chem_df)  # Reformats and caluclates values from seperate resources as dict

            chem_dict[new_entry['chemical']] = new_entry
        
        for chem in df[(df['food'] == food) & (df['chem_id'].isnull())].chemical.drop_duplicates().tolist():

            f_chem_df = df[(df['food'] == food) & (df['chemical'] == chem)]   # Filters DataFrame for rows with specific food-chemical combination
            
            new_entry = __build_chem_subdict__(f_chem_df)  # Reformats and caluclates values from seperate resources as dict
            
            chem_dict[new_entry['chemical']] = new_entry

        data_dict[food] = chem_dict
    
    return data_dict

# Converts a constructed data dictionary to a pandas dataframe
def dict_to_df(data_dict):
	df = pd.DataFrame()

	for food, food_dict in data_dict.items():
		for chemical, chem_dict in food_dict.items():
			
			chem_dict['food'] = food
			chem_dict['chemical'] = chemical

			df = df.append(chem_dict, ignore_index = True)

	return df

# Converts greek letter fillers to actual greek letters
def greek_letter_converter(chem):
    chem = chem.replace('*alpha*', 'α')
    chem = chem.replace('*beta*', 'β')
    chem = chem.replace('*gamma*', 'γ')
    chem = chem.replace('*rho*', 'ρ')
    chem = chem.replace('*delta*', 'δ')

    return chem






