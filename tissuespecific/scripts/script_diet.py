#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 10:30:19 2017

@author: acabbia
"""
import cobra 
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

# import reference model (RECON2.2)
ref_model = cobra.io.read_sbml_model(ref_model_path)
#%%
# import transciptomics data from GEO
serie=GEOparse.get_GEO(geo=GEO_accession_nr) 

######################################################
#  Part one: building the tissue specific models
#####################################################

print('Building the tissue specific models ... ')

'''
last 2 digits in GSM6375XX identify gender and age of the subject:

 -13-20 young female
 -21-31 old female
 -32-38 young male
 -39-48 old male
 
'''
#old
table = pd.DataFrame(serie.gsms['GSM637548'].table)

# make reaction confidence score
confidence= Builder.rxn_confidence(ref_model, table)
#make sure that the model can take up and use the nutrients in the inbound list
daily_bounds = Food.load_bounds_file(input_lists_dir)
inbound = list(daily_bounds.keys())

for g in inbound:
     try:
         confidence[g] = 3
         print(confidence[g])
     except:
          print(g,'not in inbounds list')
     continue
        
# build TS models:
SKM= Builder.build_model(ref_model, confidence, [])


# save
cobra.io.write_sbml_model(SKM,model_dir+'SKM_latest_'+str(datetime.date.today())+'.xml')

#%%

SKM= cobra.io.read_sbml_model('/home/acabbia/Documents/Muscle_Model/models/SKM_latest_2017-10-26.xml')


print('Adding missing reactions... ')
#old

missing = ['EX_thm_e', 'EX_nac_e', 'EX_na1_e', 'EX_but_e', 'EX_cu2_e', 'EX_zn2_e','EX_pydx_e',
           'EX_lcts_e', 'EX_k_e', 'EX_starch1200_e','EX_mnl_e','EX_caro_e','EX_10fthf_e',
           'EX_pydxn_e', 'EX_cl_e', 'EX_adocbl_e','EX_mn2_e','EX_pydam_e','EX_thf_e',
           'EX_retinol_e','EX_ncam_e','EX_ca2_e','EX_btn_e','EX_sucr_e','EX_phllqne_e',
           'SUCRe','EX_i_e','EX_mg2_e','EX_malt_e','EX_avite1_e','ATPS4m','ENO','PDHm','PYK','G3PD1','G6PDH2r']

for r in missing:
     Builder.add_bigg_reactions(SKM,[r])
     print(r, "added")
     

SKM.reactions.get_by_id('EX_h2o2(e)').bounds = 0,0 
SKM.reactions.get_by_id('EX_o2s(e)').bounds = 0,0
SKM.reactions.get_by_id('EX_atp(e)').bounds = 0,INF

print('Adding protein metabolism ... ')
MW , protein_metab_rxns = Protein_metabolism.add_prot_metabolism(SKM)

converted = Protein_metabolism.convert_to_irreversible(SKM, protein_metab_rxns)
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

old.set_uptakeFactor(1) # 70% of nutrients absorbed  <---- TO DO: make it age-dependent and decreasing
old.set_activityfactor(1.3)  ###  1 < x < 2 , ranging from sedentary to very active

print('Configure simulation settings')

#model = skeletal muscle OLD , timestep = 30 days , simulation = 12 steps

settings = Settings ( model = SKM,
                      timeStep = 15, ## days
                      simLength = 24)
settings.proteins_MW = MW

print('Load and normalize uptake bounds values from food intake')
food = Food(daily_bounds, settings, old)
#%%
prot_synth = Protein_metabolism.select_subsystem(SKM,  'protein metabolism')
configure_objective(settings.model, prot_synth, converted)
#%%
########################################################

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
   