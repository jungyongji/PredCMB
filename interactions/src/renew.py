#!/usr/bin/env python

import glob
import kegg_api
import humann_api as humann
import pandas as pd

def renew_humann_ec():
    df_ec = kegg_api.get_individual("ec")
    df_unirefec = humann.humann_api("ec","uniref90")

    temp = df_ec[df_ec.name.str.contains('Transferred')]
    temp = temp.replace(('Transferred to ',' and '),('',','),regex = True)

    for i, j in zip(temp.ec, temp.name):
        df_unirefec.replace(i, j, inplace = True)
    
    df_unirefec = \
        (df_unirefec.set_index(df_unirefec.columns.drop('ec',1).tolist())
        .ec.str.split(',', expand=True)
        .stack()
        .reset_index()
        .rename(columns={0:'ec'})
        .loc[:, df_unirefec.columns]
        )

    df_unirefec.to_csv("./database/uniref90_ec.tsv", sep = '\t', index = None)
    return df_unirefec