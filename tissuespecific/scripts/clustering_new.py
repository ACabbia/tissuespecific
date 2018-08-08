#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 14:53:00 2018

@author: acabbia
"""

from cobra.io import read_sbml_model
import GEOparse
from sklearn import metrics
from sklearn.cluster import k_means 
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
from tissuespecific.reconstruction.analysis import Summary

outfolder = '/home/acabbia/out/mergione/dist/'
library_folder = '/home/acabbia/Documents/Muscle_Model/models/library_GEOmerge/'
ref_model_file = '/home/acabbia/Documents/Muscle_Model/models/recon2.2.xml'
ref_model = read_sbml_model(ref_model_file)

GEO_accession = ["GSE25941" , "GSE28422"]

def make_exp_mat(geo_accession):
    
    expression_matrix = pd.DataFrame()
    for geo in GEO_accession:
        serie = GEOparse.get_GEO(geo)
        for gsm in serie.gsms:
            if gsm in model_list:
                index = serie.gsms[gsm].table['ID_REF']
                
                lbl = str(serie.gsms[gsm].metadata['title']).split('_')
                
                if geo == "GSE28422":
                    lbl = lbl[3].split(' ')
                    label ='ASK_'+gsm+'_'+lbl[0]
                    
                elif geo == "GSE25941":
                    label ='ASK_'+gsm+'_'+lbl[1]
                    
                expression_matrix[label] = serie.gsms[gsm].table['VALUE']
        
    expression_matrix.index = index
    return expression_matrix

# make expression matrix
model_list = [f.split('_')[1] for f in listdir(library_folder) if isfile(join(library_folder, f))]
expression_matrix = make_exp_mat(GEO_accession)

expression_matrix = expression_matrix.transform(np.log2)
# delete Affy control probes ('AFFX')
expression_matrix= expression_matrix[~expression_matrix.index.str.contains('AFFX')]

# make model matrices
ref_model = read_sbml_model(ref_model_file)
table , reactions_matrix, metabolite_matrix, gene_matrix = Summary.report_make_table(library_folder, ref_model )

#%%

def clust_kMeans_silhouette(x, n_clust, metric):
    centr, km_class , i = k_means(x.T.values, n_clust, random_state = 0)
    s = metrics.silhouette_score(x.T.values, km_class, metric)
    return s

def n_clust_silhouette(df, metric):
    sil = []
    for n_clust in range(2,9):
        sil.append(clust_kMeans_silhouette(df, n_clust ,metric))
    return sil

silhouette = pd.DataFrame(index = range(2,9))
silhouette['Expression'] = n_clust_silhouette(expression_matrix, 'cosine')
silhouette['Reactions']= n_clust_silhouette(reactions_matrix, 'jaccard')
silhouette['Metabolites']= n_clust_silhouette(metabolite_matrix, 'jaccard')
silhouette['Genes']= n_clust_silhouette(gene_matrix, 'jaccard')

ax = silhouette.plot.line()
ax.set_xlabel('Number of K-means clusters')
ax.set_ylabel('Average silhouette value')
    