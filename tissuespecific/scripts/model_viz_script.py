#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 15:07:09 2017

@author: acabbia
"""

import cobra 
from tissuespecific.reconstruction.analysis import Builder 
from tissuespecific.sarcoModel import Food
from tissuespecific.reconstruction import Utils
import GEOparse
import pandas as pd 

INF = float('inf')

#paths
path = '/home/acabbia/Documents/Muscle_Model/tissuespecific/tissuespecific'
model_dir = path + '/models'
ref_model_path = model_dir + "/recon2.2.xml"
input_lists_dir = path+'/input_lists/bounds_EU_AVG.tsv'

#Gene expression data GEO ID
GEO_accession_nr = "GSE25941"

# import reference model (RECON2.2)
ref_model = cobra.io.read_sbml_model(ref_model_path)
# import transciptomics data from GEO

'''  
serie=GEOparse.get_GEO(geo=GEO_accession_nr) 

print('Building the tissue specific models ... ')

last 2 digits in GSM6375XX identify gender and age of the subject:

 -13-20 young female
 -21-31 old female
 -32-38 young male
 -39-48 old male
 
table = pd.DataFrame(serie.gsms['GSM637548'].table)

# make reaction confidence score
confidence= Builder.rxn_confidence(ref_model, table)
#make sure that the model can take up (and use?) the nutrients in the inbound list
daily_bounds = Food.load_bounds_file(input_lists_dir)
inbound = list(daily_bounds.keys())

for g in inbound:
     try:
         confidence[g] = 3
         print(confidence[g])
     except:
          print(g, 'not in inbounds list')
     continue
        
# build TS models:
SKM= Builder.build_model(ref_model, confidence, [])
'''
# Write viz file
df = Utils.viz(ref_model,'/home/acabbia/skm_viz.tsv')
