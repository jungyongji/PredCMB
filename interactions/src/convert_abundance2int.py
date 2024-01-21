import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input', required=True, help='input as genefamilies')
parser.add_argument('-o','--output', required=True, help='ouput as genefamilies_cpm')

args = parser.parse_args()

path_genefamilies = args.input
path_int_genefamilies = args.output


df = pd.read_csv(path_genefamilies, sep = '\t')
df = df.rename(columns={'# Gene Family': 'gene_family'})
df.columns = df.columns.str.replace('_Abundance-RPKs','',regex=True)
df = df[~df.gene_family.str.contains("UNMAPPED")]


df.iloc[:,1:] = df.iloc[:,1:].round().astype("int")
df = df.reset_index(drop=True)
df.to_csv(path_int_genefamilies, sep='\t', index= None)


#python asint_genefamilies.py -i ./table/genefamilies_cpm_unstratified.tsv -o test.tsv
