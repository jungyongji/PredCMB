#!/usr/bin/env python

import sys
import pandas as pd
from Bio.KEGG import REST
from urllib.error import HTTPError
from pandas import Series, DataFrame



def get_interaction(target, source):
    try:
        interaction = REST.kegg_link(target, source).readlines()
        source_db = []
        target_db = []

        for line in interaction:
            individual = line.strip().split('\t')
            source_db.append(individual[0])
            target_db.append(individual[1])
    
        df_interaction = pd.DataFrame({source: source_db, target: target_db})
        df_interaction.replace('.+(?<=\:)','', regex=True, inplace=True)
        df_interaction.to_csv("./database/"+source+'_'+target+'.tsv', sep = '\t', index = None)
        
    except HTTPError:
        sys.exit("ERROR: Please select the available database according to the following URL. (http://www.kegg.jp/kegg/rest/keggapi.html) ")
        
    return df_interaction



def get_individual(target):
    try:
        individual = REST.kegg_list(target).readlines()
        target_db = []
        name = []

        for line in individual:
            entity = line.strip().split('\t')
            target_db.append(entity[0])
            name.append(entity[1])
            
        df_individual = pd.DataFrame({target : target_db, "name" : name})
        df_individual.replace('.+(?<=\:)','', regex=True, inplace=True)
        df_individual.to_csv("./database/"+'kegg_'+target+'.tsv', sep = '\t', index = None)

    except HTTPError:
        sys.exit("ERROR: Please select the available database according to the following URL. (http://www.kegg.jp/kegg/rest/keggapi.html) ")

    return df_individual