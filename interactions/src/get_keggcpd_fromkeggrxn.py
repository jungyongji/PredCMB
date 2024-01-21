#!/usr/bin/env python

import re



def count_df_columns(dataframe):
    return [data for data in range(len(dataframe))]

def get_dg_product(num, dataframe):
    non_dg = 0
    non_fdg = 0
    formula = dataframe.iloc[num,2]

    try:
        rxn_form = cc.parse_reaction_formula(formula)
        try:
            dg = cc.standard_dg_prime(rxn_form)
            dg = str(dg)
            dg = re.search("(?<=\()(.+)(?=\s\+)",dg).group()
            dg = float(dg)
        except Exception:
            dg = float(0)
            non_dg += 1
    except Exception:
        dg = float(0)
        non_fdg += 1

    if dg == 0 or dg == None:
        tem_form = dataframe.iloc[num,2]
        product = re.sub('\s((\+)|(\<\=\>))\s', ',', tem_form)
        
    elif dg > 0 :
        tem_form = dataframe.iloc[num,2]
        tem_form = re.search('(?<=\>\s)(.+)', tem_form).group()
        product = re.sub('\s\+\s', ',', tem_form)
                         
    elif dg < 0:
        tem_form = dataframe.iloc[num,2]
        tem_form = re.search('(.+)(?=\s\<)', tem_form).group()
        product = re.sub('\s\+\s', ',', tem_form)
        
    return { 'dG' : dg , 'cpd' : product }