#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:59:11 2018

@author: acabbia
"""
### model metrics & description

import cobra
import os
from pandas import DataFrame
import numpy as np
from tissuespecific.reconstruction.analysis import Summary

path = '/home/acabbia/Documents/Muscle_Model/models'
library_folder = path + "/library_GEOmerge/"


table = DataFrame(index = list(os.listdir(library_folder)))

rxns =[]
mets =[]
gene =[]

##### Number of reactions , metabolites, genes 

for filename in os.listdir(library_folder):
    model = cobra.io.read_sbml_model(library_folder+filename)
    print(filename)
    n_rx, n_mt , n_gn = Summary.model_summary(model)
    rxns.append(n_rx)
    mets.append(n_mt)
    gene.append(n_gn)
    
table['Number of Reactions'] = rxns
table['Number of Metabolites'] = mets
table['Number of Genes'] = gene


##### Pairwise differences in the number of reactions , metabolites, genes 

r = table['Number of Reactions']
R_diff_mat = np.zeros((len(table.index),len(table.index)))
for i in range(len(table.index)):
    for j in range(len(table.index)):
        R_diff_mat[i][j] = abs(r[i]-r[j])

R_diff_df = DataFrame(R_diff_mat, index = table.index, columns = table.index)       

m = table['Number of Metabolites']
M_diff_mat = np.zeros((len(table.index),len(table.index)))
for i in range(len(table.index)):
    for j in range(len(table.index)):
        M_diff_mat[i][j] = abs(m[i]-m[j])

M_diff_df = DataFrame(M_diff_mat, index = table.index, columns = table.index)      

g = table['Number of Genes']
G_diff_mat = np.zeros((len(table.index),len(table.index)))
for i in range(len(table.index)):
    for j in range(len(table.index)):
        G_diff_mat[i][j] = abs(g[i]-g[j])

M_diff_df = DataFrame(G_diff_mat, index = table.index, columns = table.index)