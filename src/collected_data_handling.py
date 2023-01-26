# Author: Forrest Hooton


import pandas as pd
import numpy as np
import re
from unit_converter.converter import convert, converts


def build_data_dict(df):
    """
        Receives of DataFrame of compiled data and appropriately calculates/reformats data into final structure.

        Parameters
        -----------------
        df : pd.DataFrame
            Dataframe of paper chemicals and metadata. df must be standardized for project. See Data_Statistics for example.


        Returns
        -----------------
        data_df : pd.DataFrame
            Dataframe containing compounds and corresponding metadata
    """

    data_dict = {}  # Dict that holds all foods and chemicals, and their calculations

    df['chemical'] = df['chemical'].str.strip()
    df['chemical'] = df['chemical'].apply(greek_letter_converter, convert_letter = True)

    # Iterate over each unique food, and each unique chemical per food to create food, chem dict
    for food in df['food'].drop_duplicates().tolist():

        chem_dict = {}   # Holds chemical calculations per food

        # need to do one loop for things with a chem id
        for chem_id in df[(df['food'] == food) & (df['chem_id'].notnull())].chem_id.drop_duplicates().tolist():

            f_chem_df = df[(df['food'] == food) & (df['chem_id'] == chem_id)]   # Filters DataFrame for rows with specific food-chemical combination

            new_entry = __build_chem_subdict__(f_chem_df)  # Reformats and calculates values from separate resources as dict

            chem_dict[new_entry['chemical']] = new_entry

        for chem in df[(df['food'] == food) & (df['chem_id'].isnull())].chemical.drop_duplicates().tolist():

            f_chem_df = df[(df['food'] == food) & (df['chemical'] == chem)]   # Filters DataFrame for rows with specific food-chemical combination

            new_entry = __build_chem_subdict__(f_chem_df)  # Reformats and calculates values from separate resources as dict

            chem_dict[new_entry['chemical']] = new_entry

        data_dict[food] = chem_dict

    data_df = dict_to_df(data_dict)

    return data_df


def dict_to_df(data_dict):
    """
        Converts a constructed data dictionary to a pandas dataframe

        Parameters
        -----------------
        df : pd.DataFrame
            Dataframe of paper chemicals and metadata. df must be standardized for project. See Data_Statistics for example.


        Returns
        -----------------
        data_df : pd.DataFrame
            Dataframe containing compounds and corresponding metadata gathered from inputs
    """

    df = pd.DataFrame()

    # Creates dataframe for each food and compound
    for food, food_dict in data_dict.items():
        for chemical, chem_dict in food_dict.items():

            chem_dict['food'] = food
            chem_dict['chemical'] = chemical

            df = df.append(chem_dict, ignore_index = True)

    return df


def __extract_digits__(string):
    """
        Extracts digits from beginning of string up until first non-digit character

        Parameters
        -----------------
        string : string
            Measurement string contain some digits and units


        Returns
        -----------------
        digits : int
            Digits at start of string
    """

    i = 0

    while string[i].isdigit():
        i += 1

    return int(string[0:i])


def __splice_unit__(unit):
    """
        Receives a unit in string form, even containing a number, and returns scales and units.

        Parameters
        -----------------
        unit : string
            Units of measurement


        Returns
        -----------------
        scale : int
            The scale of units, i.e. 10 in 10mg

        units : string
            The units used for a measurement without scale, i.e. mg
    """

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
    """
        Converts entry units to target units using unit_converter package

        Parameters
        -----------------
        unit : string
            Original unit of input

        target : string
            Target unit in conversion


        Returns
        -----------------
        converted unit : float
            Scale of conversion
    """

    # Package doesn't read u as nano
    if unit[0] == 'u':
        unit = 'µg'

    if target[0] == 'u':
    	target = 'µg'

    unit = '1 ' + unit

    return float(converts(unit, target))


def __unit_handler__(value, unit, target_unit):
    """
        Handles unit conversions to specified quantity, and makes adjustments as necessary.

        Parameters
        -----------------
        value : float
            Value of measurement with associated unit

        unit : string
            Original unit of input

        target_unit : string
            Target unit for conversion


        Returns
        -----------------
        value_conversion : float
            Converted value to value for target units
    """

    # Handles specific input of content ranges: "__ to __"
    if re.search(r'to', str(value)) is not None:
        value = np.mean([float(n) for n in value.split(' to ')])

    if unit == target_unit:
        return float(value)

    value = float(value)

    # Handle if unit is percentage
    if unit == '%':
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

        # Check if the first character of digit is string, assumed conversion value if it is
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
            scale_conversion = num_scale / target_num_scale


        # Receives numerical scales from converter. For __converter__('mg', 'g'), returns 1000
        num_unit_conversion = __converter__(num_unit, target_num_unit)
        denom_unit_conversion = __converter__(denom_unit, target_denom_unit)

        unit_conversion = num_unit_conversion / denom_unit_conversion

        value_conversion = value * scale_conversion * unit_conversion

        return value_conversion


def __quant_handler__(df):
    """
        Receives DataFrame for a single chemical with at least a numerical 'amount' column and calculates summary values to return.

        Parameters
        -----------------
        df : pd.DataFrame
            Dataframe of specific chemical, quantification amount, and associated metadata from papers


        Returns
        -----------------
        quant_pack : dict
            Dictionary of select quantitative information
    """

    target_unit = 'mg/100g'

    absolute_min = None
    absolute_max = None

    means =[]

    # Calculate avg mean
    for paper_id in df.PMID.drop_duplicates().tolist():

        paper_values = []

        for _, row in df[df['PMID'] == paper_id].iterrows():

            conversion = __unit_handler__(row['amount'], row['units'], target_unit)	# Sets the units for the whole output

            # If there is an error or issue in __unit_handler__, will return some sort of string
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

        # Calculates average of means from papers
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


def __build_chem_subdict__(df):
    """
        Receives DataFrame of specific chemical and its various values.

        Parameters
        -----------------
        df : pd.DataFrame
            Dataframe of specific chemical, quantification amount, and associated metadata from papers


        Returns
        -----------------
        chem_subdict : dict
            Dictionary with a chemical, amount, and metadata
    """

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

    if len(df[df.pubchem_name.notnull()]) > 0:
        compound_name = df['pubchem_name'].tolist()[0]
    else:
        compound_name = df.chemical.tolist()[0]

    chem_subdict = {
        'num_papers' : num_papers,
        'papers' : papers,
        'num_terms_quantified' : num_terms_quant,
        'chemical' : compound_name,
        'chem_id' : df.chem_id.tolist()[0],
        'pubchem_id' : df.pubchem_id.tolist()[0]
    }

    chem_subdict.update(quant_pack)

    return chem_subdict


# Converts greek letter fillers to actual greek letters
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






