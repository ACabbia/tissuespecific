#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 13:41:37 2018

@author: acabbia
"""

## distance between models 

import pandas as pd 
from cobra.io import read_sbml_model
from tissuespecific.reconstruction.analysis import Summary
from scipy.spatial.distance import pdist , jaccard , squareform
from numpy import median


outfolder = '/home/acabbia/out/mergione/dist/'
path = '/home/acabbia/Documents/Muscle_Model/models'
ref_model_file = path + "/recon2.2.xml"
library_folder = path + "/library_GEOmerge/"

###

def vec2matDist(matrix, median):
    matrix = matrix.T
    d = []
    d_df = pd.DataFrame(index = matrix.columns)
    for c in matrix.columns:
        d.append(jaccard(median, matrix[c].values))
    d_df['Jaccard dist from the median'] = d
    return d_df

def split_classes(df):
    df = df.T
    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    
    for c in df.columns:
        if c.split('_')[2] == 'old':
            df1[c] = df[c]
        elif c.split('_')[2] == 'young':
            df2[c] = df[c]
    return df1 , df2

def savefig(plot, outfolder , title):
    fig = plot
    sav = fig.get_figure()
    sav.savefig(outfolder+title, dpi=300, bbox_inches='tight')

####

ref_model = read_sbml_model(ref_model_file)
table , reactions_matrix, metabolite_matrix, gene_matrix = Summary.report_make_table(library_folder, ref_model )

reactions_matrix = reactions_matrix.T
metabolite_matrix = metabolite_matrix.T
gene_matrix = gene_matrix.T

# pairwise distance

jacc_R = squareform(pdist(reactions_matrix, metric = jaccard))
jacc_M = squareform(pdist(metabolite_matrix, metric = jaccard))
jacc_G = squareform(pdist(gene_matrix, metric=jaccard))

PW_distance_R = pd.DataFrame(jacc_R, index=reactions_matrix.index, columns= reactions_matrix.index)
PW_distance_M = pd.DataFrame(jacc_M, index=metabolite_matrix.index, columns= reactions_matrix.index)
PW_distance_G = pd.DataFrame(jacc_G, index=gene_matrix.index, columns= reactions_matrix.index)

# distance from MEDIAN

median_R = median(reactions_matrix, axis = 0)
median_M = median(metabolite_matrix, axis = 0)
median_G = median(gene_matrix, axis = 0)
    
MD_distance_R = vec2matDist(reactions_matrix , median_R)
MD_distance_M = vec2matDist(metabolite_matrix , median_M )
MD_distance_G = vec2matDist(gene_matrix, median_G)

# plot dist from median 

savefig(plot = (MD_distance_R.sort_values(by= 'Jaccard dist from the median').plot.bar(figsize = (10,10),fontsize =7.5,title = 'Distance by Reactions')),
        outfolder=outfolder, title= 'dist_Rxns.png')

savefig(plot = (MD_distance_M.sort_values(by= 'Jaccard dist from the median').plot.bar(figsize = (10,10),fontsize =7.5, title = 'Distance by Metabolites')),
        outfolder=outfolder, title= 'dist_Mets.png')

savefig(plot = (MD_distance_G.sort_values(by= 'Jaccard dist from the median').plot.bar(figsize = (10,10),fontsize =7.5, title = 'Distance by Genes')),
        outfolder=outfolder, title= 'dist_Genes.png')

# BY CLASS

old_R , young_R = split_classes(reactions_matrix)
old_M , young_M = split_classes(metabolite_matrix)
old_G , young_G = split_classes(gene_matrix)

median_old_R = median(old_R.T, axis = 0)
median_old_M = median(old_M.T, axis = 0)
median_old_G = median(old_G.T, axis = 0)
    
MD_distance_old_R = vec2matDist(old_R.T, median_old_R)
MD_distance_old_M = vec2matDist(old_M.T, median_old_M)
MD_distance_old_G = vec2matDist(old_G.T, median_old_G)

median_young_R = median(young_R.T, axis = 0)
median_young_M = median(young_M.T, axis = 0)
median_young_G = median(young_G.T, axis = 0)
    
MD_distance_young_R = vec2matDist(young_R.T, median_young_R)
MD_distance_young_M = vec2matDist(young_M.T, median_young_M)
MD_distance_young_G = vec2matDist(young_G.T, median_young_G)        


#plot CLASS dist from median

savefig(plot=(MD_distance_old_R.sort_values(by= 'Jaccard dist from the median').plot.bar(figsize = (10,10),fontsize =7.5,title='class: old, Distance by Reactions')),
        outfolder=outfolder, title= 'dist_old_Rxns.png')

savefig(plot= (MD_distance_old_M.sort_values(by= 'Jaccard dist from the median').plot.bar(figsize = (10,10),fontsize =7.5,title='class: old, Distance by  Metabolites')),
        outfolder=outfolder, title= 'dist_old_Mets.png')

savefig(plot= (MD_distance_old_G.sort_values(by= 'Jaccard dist from the median').plot.bar(figsize = (10,10),fontsize =7.5,title='class: old, Distance by Genes')),
        outfolder=outfolder, title= 'dist_old_Genes.png')

savefig(plot= (MD_distance_young_R.sort_values(by= 'Jaccard dist from the median').plot.bar(figsize = (10,10),fontsize =7.5,title='class: young, Distance by Reactions')),
        outfolder=outfolder, title= 'dist_young_Rxns.png')

savefig(plot= (MD_distance_young_M.sort_values(by= 'Jaccard dist from the median').plot.bar(figsize = (10,10),fontsize =7.5,title='class: young, Distance by Metabolites')),
        outfolder=outfolder, title= 'dist_young_Mets.png')

savefig(plot=(MD_distance_young_G.sort_values(by= 'Jaccard dist from the median').plot.bar(figsize = (10,10),fontsize =7.5,title='class: young, Distance by Genes')),
        outfolder=outfolder, title= 'dist_young_Genes.png')
        
#rxn
# SD dist from median
o_sd_R = MD_distance_old_R.std()
y_sd_R = MD_distance_young_R.std()
print(o_sd_R)
print(y_sd_R)
# dist between class medians
R_class_dist = jaccard(median_old_R,median_young_R)
print(R_class_dist)

# mets
# SD dist from median
o_sd_M = MD_distance_old_M.std()
y_sd_M = MD_distance_young_M.std()
print(o_sd_M)
print(y_sd_M)
# dist between class medians
M_class_dist = jaccard(median_old_M,median_young_M)
print(M_class_dist)

#genes
# SD dist from median
o_sd_G = MD_distance_old_G.std()
y_sd_G = MD_distance_young_G.std()
print(o_sd_G)
print(y_sd_G)
# dist between class medians
G_class_dist = jaccard(median_old_G,median_young_G)
print(G_class_dist)










        