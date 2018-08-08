#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 11:12:16 2018

@author: acabbia
"""
######## MAIN ##########
'''
1-model building
2-report
3-distance/clustering

'''
#########################

from build_models_library import build_models_library
from report import report

#set paths
path = '/home/acabbia/Documents/Muscle_Model/models'
ref_model_path = path + '/recon2.2.xml'
library_folder =  path + "/library_GEOmerge_batch/"
outfolder = '/home/acabbia/out/mergione/batch'
#Gene expression data GEO ID (Raue 2012)
GEO_accession_nr = "GSE25941"

# required reactions  
requirements = ['ATPS4m','ENO','PDHm','PYK','G3PD1','G6PDH2r','AKGDm','CYOOm3','r0913','FAOXC4020m','FAOXC80','EX_glc(e)',
                 'EX_hdca(e)','EX_ocdca(e)','EX_arach(e)','EX_doco13ac_','EX_lgnc(e)','EX_ala_L(e)',
                 'EX_arg_L(e)','EX_asn_L(e)','EX_asp_L(e)','EX_cys_L(e)','EX_gln_L(e)','EX_glu_L(e)',
                 'EX_gly(e)','EX_his_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)','EX_met_L(e)',
                 'EX_phe_L(e)','EX_pro_L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)',
                 'EX_val_L(e)','EX_fru(e)','EX_ppa(e)','EX_but(e)','EX_hx(e)','EX_octa(e)','EX_ttdca(e)',
                 'EX_h2o(e)','EX_hco3(e)','EX_co2(e)','EX_h(e)']

# expressed/not expressed cutoff (quantile)
expr_cutoff = 0.65

#%% Step 1: Build library of PD-GSM from GEO expression data
#1)
build_models_library(ref_model_path, GEO_accession_nr, requirements, expr_cutoff, library_folder)

#%% Step 2: Generate report models content and summary figures 
#2)
table = report(ref_model_path,library_folder,outfolder)

print('Average number of reactions:', table['nr. of Reactions'].mean(), ', SD:', table['nr. of Reactions'].std())
print('Average number of metabolites:', table['nr. of Metabolites'].mean(), ', SD:', table['nr. of Metabolites'].std())
print('Average number of genes:', table['nr. of Genes'].mean(), ', SD:', table['nr. of Genes'].std())
#%%