#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 14:09:09 2018

@author: acabbia
"""

import cobra
import pandas as pd
import numpy as np

from sklearn.decomposition import PCA

path = '/home/acabbia/Documents/Muscle_Model/models'
ref_model = path + "/recon2.2.xml"
model = path + "/SKM_v1.1.xml"


# model (SKM v1.1)
skm = cobra.io.read_sbml_model(model)
model= cobra.io.read_sbml_model(model)

#%%
# optimize for a range of oxygen uptake values

flux_matrix = pd.DataFrame(index=[r.id for r in skm.reactions])
for g in range(0, -1001, -1):  # ends at -1000
    skm.reactions.get_by_id('EX_o2(e)').lower_bound = g
    sol = skm.optimize()
    flux_matrix['EX_o2_lb = '+ str(g)]=sol.fluxes.values
        
#PCA of flux matrix
pca = PCA()
pca.fit(flux_matrix.T)
exp_var = pd.Series(pca.explained_variance_ratio_)
exp_var.plot.bar()
#scale loadings by loading of o2 uptake
c1_loadings = pd.Series(pca.components_[0]*np.sqrt(pca.explained_variance_[0]), index= [r.id for r in skm.reactions])
c1_loadings = c1_loadings.T
norm_ld = c1_loadings/c1_loadings['EX_o2(e)']
norm_ld[abs(norm_ld) < 10e-9] = 0 
# plot loadings > 10*e-9 
nnz_ld = norm_ld[abs(norm_ld) > 10e-9]
nnz_ld = nnz_ld*-1 
nnz_ld= nnz_ld.sort_values(ascending = False)
nnz_ld.plot.bar(figsize = (36,18))

#%%
#make recon map overlay
filename='test.csv'  
      
color_code=[]         
for i ,v in nnz_ld.iteritems():
    if nnz_ld.loc[i] > 0:
        color_code.append('#00ff00') #green
    elif nnz_ld.loc[i] < 0:
        color_code.append('#ff0000') #red
    else:
        color_code.append('#ffff00') #yellow
                     
overlay = pd.DataFrame({'color':color_code,
                        'lineWidth':abs(nnz_ld.values)*2,
                        'reactionIdentifier':nnz_ld.index                       
                        })   
 
overlay = overlay[['reactionIdentifier','lineWidth','color']] 
overlay.to_csv('/home/acabbia/'+filename, sep= '\t', header = True, index = False)         
         
###     # reset bounds
for r in skm.exchanges:
    r.bounds = -1000 , 1000
    
skm.reactions.get_by_id('EX_fvs(e)').bounds = 0 , 1000

#%%  
def oneD_perturbation(model, reaction_id , lb_range):
    #### works only for EXchange reactions !!
    
    #initialize flux matrix and range function
    flux_matrix = pd.DataFrame(index=[r.id for r in model.reactions])
    #optimization loop 
    for g in range(lb_range[0],lb_range[1],lb_range[2]):
        model.reactions.get_by_id(reaction_id).lower_bound = g
        sol = model.optimize()         
        flux_matrix[reaction_id +'_lb = '+ str(g)]=sol.fluxes.values
    
    ### Check for and delete reactions (rows) with std = 0
    flux_matrix = flux_matrix[flux_matrix.std(1) !=0]
    #PCA
    pca = PCA(n_components = 5)
    pca.fit(flux_matrix.T)
    exp_var = pd.Series(pca.explained_variance_ratio_)
    exp_var.plot.bar()
    #scale loadings by 'Reaction_id' loading value
    c1_loadings = pd.Series(pca.components_[0]*np.sqrt(pca.explained_variance_[0]),
                            index= [r for r in flux_matrix.index])
    c1_loadings = c1_loadings.T
    try:
        norm_ld = c1_loadings/c1_loadings[reaction_id]   #### warning : division by zero may happen here
    except:
        norm_ld = c1_loadings
        print('Unscaled loadings')
        pass
    
    norm_ld[abs(norm_ld) < 10e-3] = 0 
    
    # plot loadings > 10*e-6 
    nnz_ld = norm_ld[norm_ld != 0]
    nnz_ld = nnz_ld*-1 
    nnz_ld= nnz_ld.sort_values(ascending = False)
    nnz_ld.plot.bar(figsize = (36,18), fontsize = 9)
    return nnz_ld

def oneD_overlay(nnz_ld , filename):
    color_code=[]         
    for i ,v in nnz_ld.iteritems():
        if nnz_ld.loc[i] > 0:
            color_code.append('#00ff00') #green
        elif nnz_ld.loc[i] < 0:
            color_code.append('#ff0000') #red
        else:
            color_code.append('#ffff00') #yellow
                     
    overlay = pd.DataFrame({'color':color_code,
                        'lineWidth':abs(nnz_ld.values)*2,
                        'reactionIdentifier':nnz_ld.index                       
                        })   
 
    overlay = overlay[['reactionIdentifier','lineWidth','color']] 
    overlay.to_csv('/home/acabbia/RMap/'+filename, sep= '\t', header = True, index = False)
    return overlay    