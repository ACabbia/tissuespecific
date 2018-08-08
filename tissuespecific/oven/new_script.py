#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:27:21 2017

@author: acabbia
"""

from cobra.io import read_sbml_model , write_sbml_model

from tissuespecific.reconstruction.protein_metabolism import Protein_metabolism
from tissuespecific.reconstruction.analysis import Builder
from tissuespecific.sarcoModel.individual import Individual, Settings, Food
from tissuespecific.sarcoModel.sim import runSimulation
from tissuespecific.sarcoModel.config_objective import configure_objective

import os
import GEOparse
import pandas as pd
import datetime
import matplotlib.pyplot as plt

INF = float('inf')
#paths
path = os.path.dirname(__file__)
model_dir = path + '/models'
ref_model_path = model_dir + "/recon2.2.xml"
input_lists_dir = path+'/input_lists/bounds_EU_AVG.tsv'

#Gene expression data GEO ID
GEO_accession_nr = "GSE25941"

#%%
serie = GEOparse.get_GEO(GEO_accession_nr) 
# import reference model (RECON2.2)
ref_model = read_sbml_model(ref_model_path)

######################################################
#  Part one: building the tissue specific models
#####################################################

print('Building the tissue specific model ... ')

'''
last 2 digits in GSM6375XX identify gender and age of the subject:

 -13-20 young female
 -21-31 old female
 -32-38 young male
 -39-48 old male
 
'''
#old male 
table = pd.DataFrame(serie.gsms['GSM637548'].table)

# make reaction confidence score:
confidence = Builder.rxn_confidence(ref_model, table)

# assign max confidence to reactions in the inbounds list:

daily_bounds = Food.load_bounds_file(input_lists_dir)
inbound = list(daily_bounds.keys())

## testing
'''
for i in inbound:
     try:
          confidence[i] = 3
          print(i, confidence[i])
     except:
          print(i, ' is not in inbound fluxes list')
     continue

# ensure functionality of reconstructed network (glycolisis, TCA, glycogen storage , FA_oxidation (to check))
to_add = ['ATPS4m','ENO','PDHm','PYK','G3PD1','G6PDH2r']

for a in to_add:
     try:
          confidence[a] = 3
          print(a, confidence[a])
     except:
          print(a, ' is not in inbound fluxes list')
     continue
'''
###

to_include = []
to_include.append('glc_D_e')

for f in inbound:
    g = f.replace('EX_','')
    h = g.replace('(e)','_e')
    try:
         to_include.append(h)
         ref_model.metabolites.get_by_id(h)
    except:
         print(h, 'not in ref_model')
         to_include.pop()
    continue
    
# build TS models:
SKM = Builder.build_model(ref_model, confidence, to_include)

# save
write_sbml_model(SKM,model_dir+'SKM_latest_'+str(datetime.date.today())+'.xml')
#%%

print('Adding protein metabolism ... ')
MW = Protein_metabolism.add_prot_metabolism(SKM)

prot_synth=Protein_metabolism.select_subsystem(SKM, 'protein metabolism')
converted = Protein_metabolism.convert_to_irreversible(SKM, prot_synth)
print('protein metabolism converted to irreversible')

#%%
######################################################
#  Part two: simulating growth
#####################################################

print('Making Individual object to store subjects data')
# male, age 78 , 62 kg
old = Individual(name = 'old',
                 isMale = True, 
                 age = 78,     ## Years <--------- TO DO : linked with uptakeFactor
                 weight = 65)  ## Kg

old.set_fat_ratio (.2)  # 20% body fat
old.set_leanMass_maintenance_cost(60) # kcal/g  <-----TO DO: verify!!! (value taken from stig-met toolbox)

old.set_uptakeFactor(.7) # 70% of nutrients absorbed  <---- TO DO: make it age-dependent and decreasing
old.set_activityfactor(1.3)  ###  1 < x < 2 , ranging from sedentary to very active

print('Configure simulation settings')

#model = skeletal muscle OLD , timestep = 30 days , simulation = 12 steps

settings = Settings ( model = SKM,
                      timeStep = 15, ## days
                      simLength = 300)
settings.proteins_MW = MW

print('Load and normalize uptake bounds values from food intake')

food = Food(daily_bounds, settings, old)
#%%
configure_objective(settings.model, prot_synth, converted)
#%%
########################################################SKM

print( 'Running simulation...')
results, atp_exp, obj_val , TEE = runSimulation(old,settings,food)

#plots
results = pd.Series(results)
results.plot(title='lean mass (Kg)')

plt.figure()
atp_exp = pd.Series(atp_exp)
atp_exp.plot(title= 'ATP_expenditure constraint value')

plt.figure()
obj_val = pd.Series(obj_val)
obj_val.plot(title='objective value')

plt.figure()
TEE = pd.Series(TEE)
TEE.plot(title= 'total energetic expenditure (Kcal/day)')

SKM.summary()






