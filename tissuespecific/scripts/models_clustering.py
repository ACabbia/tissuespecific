#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 16:34:46 2018

@author: acabbia
"""
import os
import cobra as cb
import pandas as pd 
import seaborn as sns

from scipy.spatial.distance import pdist , jaccard , squareform
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage


#set paths
path = '/home/acabbia/Documents/Muscle_Model/models'
ref_model = path + "/recon2.2.xml"
library_folder = path + "/library_GEO_GSE25941_v2/"

# import reference model (RECON2.2)
recon22 = cb.io.read_sbml_model(ref_model)

#%% cluster (dendro) 

reactions_matrix = pd.DataFrame(index=[r.id for r in recon22.reactions])
metabolite_matrix = pd.DataFrame(index=[m.id for m in recon22.metabolites])
gene_matrix = pd.DataFrame(index=[g.id for g in recon22.genes])

for model_name in os.listdir(library_folder):
    model = cb.io.read_sbml_model(library_folder+model_name)
    rxns = []
    mets = []
    genes = []
    
    for r in recon22.reactions:
        if r in model.reactions:
            rxns.append(1)
        else:
            rxns.append(0)
            
    for m in recon22.metabolites:
        if m in model.metabolites:
            mets.append(1)
        else:
            mets.append(0)
            
    for g in recon22.genes:
        if g in model.genes:
            genes.append(1)
        else:
            genes.append(0)
    
    label = str(model_name).split('_')[1].split('.')[0]
    
    sr = pd.Series(rxns)
    sm = pd.Series(mets)
    sg = pd.Series(genes)
    
    gene_matrix[label] = sg.values
    reactions_matrix[label] = sr.values
    metabolite_matrix[label] = sm.values
    
reactions_matrix = reactions_matrix.T
metabolite_matrix = metabolite_matrix.T
gene_matrix = gene_matrix.T

#%%
Z_r = linkage(reactions_matrix, 'ward')
Z_m = linkage(metabolite_matrix, 'ward')
Z_g = linkage(gene_matrix, 'ward')

#plots

plt.subplot()
plt.figure(figsize=(20, 10))
plt.title('Hierarchical Clustering by Reactions')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(
    Z_r,
    leaf_rotation=90,  
    leaf_font_size=8,  
    labels=reactions_matrix.index
)

plt.subplot()
plt.figure(figsize=(20, 10))
plt.title('Hierarchical Clustering by Metabolites')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(
    Z_m,
    leaf_rotation=90,  
    leaf_font_size=8,
    labels = metabolite_matrix.index
)

plt.subplot()
plt.figure(figsize=(20, 10))
plt.title('Hierarchical Clustering by Genes')
plt.xlabel('sample index')
plt.ylabel('distance')

dendrogram(
    Z_g,
    leaf_rotation=90,  
    leaf_font_size=8,
    labels = gene_matrix.index
)

plt.show()

#%%
# dendro + heatmap plot

#pairwise jaccard distance between the 36 models
jacc_R = squareform(pdist(reactions_matrix, metric = jaccard))
jacc_M = squareform(pdist(metabolite_matrix, metric = jaccard))
jacc_G = squareform(pdist(gene_matrix, metric=jaccard))

distance_R = pd.DataFrame(jacc_R, index=reactions_matrix.index, columns= reactions_matrix.index)
distance_M = pd.DataFrame(jacc_M, index=metabolite_matrix.index, columns= reactions_matrix.index)
distance_G = pd.DataFrame(jacc_G, index=gene_matrix.index, columns= reactions_matrix.index)

#TODO: plot title
sns.clustermap(distance_R)
----ax.set_title('Cluster by Reactions')

sns.clustermap(distance_M)
sns.clustermap(distance_G)
#%%
#number of unique reactions:

R_unique = reactions_matrix.T[reactions_matrix.sum() == 1]
print('Number of unique reactions for each model:',R_unique.sum())

M_unique = metabolite_matrix.T[metabolite_matrix.sum() == 1]
print('Number of unique metabolites for each model:',M_unique.sum())

G_unique = gene_matrix.T[gene_matrix.sum() == 1]
print('Number of unique genes for each model',R_unique.sum())
























