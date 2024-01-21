#!/usr/bin/env python

import os
import sys
import glob
import gzip
import tarfile
import pandas as pd
from urllib.request import urlretrieve


def make_database():
    if not os.path.isdir('database'):
        os.mkdir("database")
        
        
def check_humannapi_input(function, genefamily):
    if ( function =='ec' or function =='ko' ) and ( genefamily == 'uniref50' or genefamily == 'uniref90' ):
        return True
    else :
        sys.exit("ERROR: Please select an available group and uniref. ( group : ko or ec , uniref : uniref50 or uniref90 )")


def gz_preprocessing(gz_file, group, uniref):
    gz_df = pd.read_csv(gz_file, sep = 'delimiter', header=None, engine='python')

    gz_df.columns = [group]

    # Divide multiple values in one column into a separator and add them to a new column
    gz_df = gz_df.join(gz_df[group].str.split('\t', expand=True).add_prefix('ID'))

    # Combine values of multiple columns into one column
    gz_df[uniref] = gz_df[gz_df.columns[2:]].apply(
        lambda x : ','.join(x.dropna().astype(str)),
        axis=1)

    gz_df = gz_df.iloc[:,[-1,1]]
    gz_df.columns = [uniref, group]

    gz_df.to_csv("./database/" + uniref + '_' + group + ".tsv", sep = '\t', index = None)
    return gz_df
    
    
def humann_api(group, uniref):
    """
    Mapping files between databases were downloaded from the api server of the laboratory that developed humann.
    Even though there are other faster and simpler methods, the reason why I used them is because,

    1. When referring to the mapping file of the human library installed by the user,
    it is disadvantageous to users who downloaded multiple versions of human
    (because they have multiple versions of mapping files, so mapping files cannot be distinguished by version).
    (Rather, if you make good use of this, we expect that a specific version of the mapping file can be used.)

    2. Although the mapping file provided in humann github may be available,
    it may become unusable if the repertoire changes or url changes.
    So it takes a long time, but it is recommended to download it directly from the api server and use only the mapping file you want.
    """

    if check_humannapi_input(group, uniref):
        target = "./database/" + uniref + '_' + group + ".tsv"
        if glob.glob(target):
            df_unirefec = pd.read_csv(target, sep = '\t')
            return df_unirefec
        
        if glob.glob('./database/map_level4ec_uniref90.txt.gz'):
            return gz_preprocessing('./database/map_level4ec_uniref90.txt.gz', group, uniref)
        
        if not glob.glob('./database/full_mapping_v201901b.tar.gz'):
            url = "http://huttenhower.sph.harvard.edu/humann_data/full_mapping_v201901b.tar.gz"
            url_file = url.split('/')[-1]
            urlretrieve(url, './database/' + url_file)
            
        # Decompress mapping files
        taropen = tarfile.open('./database/full_mapping_v201901b.tar.gz', "r:gz")
        for tarlist in taropen:
            tarinfo = tarlist.name
            if group in tarinfo and uniref in tarinfo:
                filename = tarinfo
                taropen.extract(filename, "./database/")

        ## Processing
        gz_file = './database/' + filename
        taropen.close()
        return gz_preprocessing(gz_file, group, uniref)
        
    else :
        sys.exit("ERROR: Please select an available group and uniref.")