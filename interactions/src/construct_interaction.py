#!/usr/bin/env python

import renew
import pandas as pd
import kegg_api as kegg
import humann_api as humann

import get_keggrxn_equation as gre
import get_keggcpd_fromkeggrxn as gcr

from functools import partial
from multiprocessing import Pool
from equilibrator_api import ComponentContribution




def main():

    ## Obtain an equation from KEGG reaction ID
    rxn_df = kegg.get_individual('reaction')
    equation_df = pd.DataFrame(map(gre.get_keggrxn_equation, gre.queue_pop_column_from_df(rxn_df, 'reaction')))
    equation_df.to_csv("./database/keggrxn_equation.tsv", sep = '\t', index = None)




    ## Get compounds from reaction equations
    cc = ComponentContribution()
    non_dg = 0
    non_fdg = 0

    coef_list = ['\d\d C', '\d C']
    n_list = [' \w ', ' ', ' \W \w ',' \d \w ',' \w\-\d \w ',' \(side \d\) ', ' \(\w\) ', ' \d\w ', ' \(\w\-\w\) ', ' \(\w\+\w\) ', ' \(\w\-\d\) ', ' \(\w\+\d\) ']

    for coef in coef_list:
        equation_df.replace(coef, 'C', regex = True, inplace = True)
    for n in n_list:
        equation_df.replace(n, ' ', regex = True, inplace = True)

    pool = Pool(processes = 4)
    result = pool.map(partial(gcr.get_dg_product, dataframe = equation_df), gcr.count_df_columns(equation_df))
    print(f" # of Non-ΔG'0 = {non_dg}")
    print(f" # of Non-ΔfG' = {non_fdg}")
    pool.close()
    pool.join()
    cpd_df = pd.DataFrame(result)
    total_df = pd.concat([equation_df, cpd_df], axis = 1 )
    '''
    The duplicated product in reaction is going to be omitted
    ex) Reaction : C04150 + G10509 = C02970 + G10509
    Reaction's product : C04150,C02970,G10509
    '''
    temp = []
    for individual_cpds in total_df.cpd:
        cpd_propr = ','.join(set(individual_cpds.split(',')))
        temp.append(cpd_propr)
    total_df.cpd = temp
    total_df.to_csv("./database/kegg_rxn_cpd_equilibrator.tsv", sep = '\t', index = None)




    ## Assign glycans to compounds
    rxn_mixcpd = total_df[['ENTRY', 'cpd']]
    rxn_mixcpd.columns = ['reaction','cpd']

    rxn_mixcpd = \
    (rxn_mixcpd.set_index(rxn_mixcpd.columns.drop('cpd',1).tolist())
    .cpd.str.split(',', expand=True)
    .stack()
    .reset_index()
    .rename(columns={0:'cpd'})
    .loc[:, rxn_mixcpd.columns]
    )

    rxn_cpd_nog = rxn_mixcpd[~rxn_mixcpd.cpd.str.contains('G')]
    rxn_cpd_g = rxn_mixcpd[rxn_mixcpd.cpd.str.contains('G')]
    rxn_cpd_g.columns = ['reaction','glycan']

    glycan_cpd = kegg.get_interaction("cpd","glycan")
    rxn_glycan_cpd = pd.merge(rxn_cpd_g, glycan_cpd, on = 'glycan')

    rxn_cpd = rxn_glycan_cpd[['reaction','cpd']]
    concat_cpd = pd.concat([rxn_cpd_nog, rxn_cpd]).drop_duplicates()
    concat_cpd.to_csv("./database/kegg_rxn_cpd_nog.tsv", sep = '\t', index = None)




    ## Removal of generic compound
    generic_cpd = pd.read_csv("./database/kegg_generic_cpd.tsv", sep = '\t')
    kegg_rxn_cpd_noggeneric = concat_cpd[~concat_cpd.cpd.isin(generic_cpd.kegg)]
    kegg_rxn_cpd_noggeneric.to_csv("./database/kegg_rxn_cpd_noggeneric.tsv", sep = '\t', index = None)





    ## Map uniref90 to compound


    # EC - Reaction - Cpd
    uniref90_ec = renew.renew_humann_ec()
    ec_rn = kegg.get_interaction("reaction","ec")

    uni_ec_rn = pd.merge(uniref90_ec, ec_rn, on = 'ec')
    uni_ec_rn_cpd = pd.merge(uni_ec_rn, kegg_rxn_cpd_noggeneric, on = 'reaction')
    erc = uni_ec_rn_cpd[['uniref90','cpd']].drop_duplicates()


    # KO - EC - Reaction - Cpd
    uniref90_ko = humann.humann_api("ko", "uniref90")
    ko_ec = kegg.get_interaction("ec","ko")

    uni_ko_ec = pd.merge(uniref90_ko, ko_ec, on = 'ko')
    uni_ko_ec_rn = pd.merge(uni_ko_ec, ec_rn, on = 'ec')
    uni_ko_ec_rn_cpd = pd.merge(uni_ko_ec_rn, kegg_rxn_cpd_noggeneric, on = 'reaction')
    kerc = uni_ko_ec_rn_cpd[['uniref90','cpd']].drop_duplicates()


    # KO - Reaction - Cpd
    ko_rn = kegg.get_interaction("reaction","ko")

    uni_ko_rn = pd.merge(uniref90_ko, ko_rn, on = 'ko')
    uni_ko_rn_cpd = pd.merge(uni_ko_rn, kegg_rxn_cpd_noggeneric, on = 'reaction')
    krc = uni_ko_rn_cpd[['uniref90','cpd']].drop_duplicates()
    

    # # Merge
    kerkc = pd.concat([erc, kerc, krc]).drop_duplicates()

    kerkc = \
    (kerkc.set_index(kerkc.columns.drop('uniref90',1).tolist())
    .uniref90.str.split(',', expand=True)
    .stack()
    .reset_index()
    .rename(columns={0:'uniref90'})
    .loc[:, kerkc.columns]
    )

    kerkc.to_csv("./database/uniref_com_kerk_oc.tsv", sep = '\t', index = None)


if __name__ == '__main__':
    main()