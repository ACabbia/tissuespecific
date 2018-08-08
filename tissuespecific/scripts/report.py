#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 16:21:34 2018

@author: acabbia
"""

#model library report and figures

from cobra.io import read_sbml_model
from tissuespecific.reconstruction.analysis import Summary
import pandas as pd
from scipy.stats import ttest_ind 
from statsmodels.stats import multitest

path = '/home/acabbia/Documents/Muscle_Model/models'
ref_model_path = path + "/recon2.2.xml"
library_folder = path + "/library_GEOmerge/"
outfolder = '/home/acabbia/out/mergione'

def report(ref_model_path,library_folder,outfolder):
    
    ref_model = read_sbml_model(ref_model_path)
    
    table, reactions_matrix, metabolite_matrix, gene_matrix = Summary.report_make_table(library_folder, ref_model)
    table, R_pw_diff, M_pw_diff, G_pw_diff = Summary.report_clustering_plots(table , 
                                                                             reactions_matrix,
                                                                             metabolite_matrix,
                                                                             gene_matrix,
                                                                             outfolder)
    
    '''
    R_pw_sim = Summary.PW_similarity(reactions_matrix)
    M_pw_sim = Summary.PW_similarity(metabolite_matrix)
    G_pw_sim = Summary.PW_similarity(gene_matrix)
    
    R_PW = Summary.merge_PW_df(R_pw_sim, R_pw_diff)
    M_PW = Summary.merge_PW_df(M_pw_sim, M_pw_diff)
    G_PW = Summary.merge_PW_df(G_pw_sim, G_pw_diff)
    
    Summary.heat(R_PW, outfolder +'/reactions_PW.png')
    Summary.heat(M_PW ,outfolder +'/metabolites_PW.png')
    Summary.heat(G_PW, outfolder +'/genes_PW.png'  )
    '''
    
    Summary.table_as_png(table, outfolder +'/library_sumary.png')
    
    return table

def split_classes(df):
    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    
    for c in df.columns:
        if c.split('_')[2] in ['old','Old']:
            df1[c] = df[c]
        elif c.split('_')[2] in ['young', 'Young']:
            df2[c] = df[c]
    return df1 , df2


def t_test_FDR(matrix):
    '''
    Returns variables significantly (0.05+FDR) present in one of the two classes
    '''
    tt = []
    pp = []
    result = pd.DataFrame(index = matrix.index)
    
    old , young = split_classes(matrix)
    
    for i in range(len(matrix.index)):
        t,p = ttest_ind(old.values[i],young.values[i], equal_var= False)
        tt.append(t)
        pp.append(p)
        
    result['t statistic'] = tt
    result['P value'] = pp
    result = result[-result['t statistic'].isna()]
    result = result.sort_values(by='P value')
    result['Adj. P value'] = multitest.fdrcorrection(result['P value'])[1]
    result['Accepted'] = multitest.fdrcorrection(result['P value'])[0]
    result = result[result['Accepted'] == True]
    
    return result            
