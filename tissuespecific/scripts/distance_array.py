#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 16:54:27 2018

@author: acabbia
"""

# distance in microarray data 

from cobra.io import read_sbml_model
import GEOparse
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
from scipy.spatial.distance import pdist , correlation , cosine,  squareform
from numpy import median

outfolder = '/home/acabbia/out/mergione/dist/'

library_folder = '/home/acabbia/Documents/Muscle_Model/models/library_GEOmerge/'
model_file = '/home/acabbia/Documents/Muscle_Model/models/recon2.2.xml'
ref_model = read_sbml_model(model_file)

GEO_accession = ["GSE25941" , "GSE28422"]

###

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

def vec2matDist(matrix, median):
    matrix = matrix.T
    d = []
    d_df = pd.DataFrame(index = matrix.columns)
    for c in matrix.columns:
        d.append(cosine(median, matrix[c].values))
    d_df['Cosine dist from the median'] = d
    return d_df

def split_classes(df):
    df = df.T
    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    
    for c in df.columns:
        if c.split('_')[2] == 'Old':
            df1[c] = df[c]
        elif c.split('_')[2] == 'Young':
            df2[c] = df[c]
    return df1 , df2

def savefig(plot, outfolder , title):
    fig = plot
    sav = fig.get_figure()
    sav.savefig(outfolder+title, dpi=300, bbox_inches='tight')

###

# Make expression matrix only for gsm in model folder
model_list = [f.split('_')[1] for f in listdir(library_folder) if isfile(join(library_folder, f))]
e_mat = make_exp_mat(GEO_accession)

# log2 transform
e_mat = e_mat.transform(np.log2)
# delete Affy control probes ('AFFX')
e_mat = e_mat[~e_mat.index.str.contains('AFFX')]
e_mat = e_mat.T
# pairwise distance
PW_e_mat = squareform(pdist(e_mat, metric = cosine))
PW_distance_E = pd.DataFrame(PW_e_mat, index=e_mat.index, columns= e_mat.index)
# distance from median
median_E = median(e_mat, axis = 0)
MD_distance_E = vec2matDist(e_mat, median_E)

# SD dist from median
MD_distance_E.std()


#### BY CLASS
old_E , young_E = split_classes(e_mat)

median_old_E = median(old_E.T, axis = 0)
median_young_E = median(young_E.T, axis = 0)

MD_distance_old_E = vec2matDist(old_E.T, median_old_E)
MD_distance_young_E = vec2matDist(young_E.T, median_young_E)

# SD dist from median
o_sd_E = MD_distance_old_E.std()
y_sd_E = MD_distance_young_E.std()

# dist between class medians

E_class_dist = cosine(median_old_E, median_young_E)
print(E_class_dist)

#plot
savefig(plot=(MD_distance_old_E.sort_values(by= 'Cosine dist from the median').plot.bar(figsize = (10,10),fontsize =7.5,title='class: old, Distance by Expression')),
        outfolder=outfolder, title= 'dist_old_Exp.png')

savefig(plot=(MD_distance_young_E.sort_values(by= 'Cosine dist from the median').plot.bar(figsize = (10,10),fontsize =7.5,title='class: old, Distance by Expression')),
        outfolder=outfolder, title= 'dist_young_Exp.png')