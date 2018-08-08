#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 13:58:41 2018

@author: acabbia
"""
import cobra
import os
import seaborn as sns
import numpy as np
from pandas import DataFrame
from tissuespecific.reconstruction.utils import Utils
from matplotlib.pyplot import show , savefig


def RQ(model):
    model.objective = 'ATPS4m'
    sol = model.optimize()
    rq = np.nan
    try:
        nonzero_fluxes = sol.fluxes[sol.fluxes != 0]
    
        if 'EX_co2(e)' in nonzero_fluxes.index.values and 'EX_hco3(e)' not in nonzero_fluxes.index.values:
            
            rq = abs(nonzero_fluxes['EX_co2(e)'])/abs(nonzero_fluxes['EX_o2(e)'])
        elif 'EX_hco3(e)' in nonzero_fluxes.index.values and 'EX_co2(e)' not in nonzero_fluxes.index.values:
           
            rq = abs(nonzero_fluxes['EX_hco3(e)'])/abs(nonzero_fluxes['EX_o2(e)'])
        elif 'EX_co2(e)' and 'EX_hco3(e)' in nonzero_fluxes.index.values:
          
            rq = abs(nonzero_fluxes['EX_hco3(e)'])+abs(nonzero_fluxes['EX_co2(e)'])/abs(nonzero_fluxes['EX_o2(e)'])
    except:
        print('infeasible solution')
        rq = np.nan
    return rq

def RQ_test(model):
    'Functional Metabolic tests'
    rq_list = []
    for carbon_source in [
                # sugars
                'EX_glc(e)',
                'EX_fru(e)',
                 #fatty acids
                'EX_octa(e)',       # C8:0
                'EX_ttdca(e)',      # C14:0
                'EX_hdca(e)',       # C16:0
                'EX_ocdca(e)',      # C18:0
                'EX_arach(e)',      # C20:0
                'EX_doco13ac_',     # C22:0
                'EX_lgnc(e)',        # C24:0
                # Unsaturated Fatty acids
                'EX_lnlnca(e)',      # C18:3
                'EX_hdcea(e)',       # C16:1
                'EX_lnlc(e)'        # C18:2
                 ]:
        
          #close all inflows , leave outflows open
          ex = Utils.find_exchanges(model)
          for r in ex:
               model.reactions.get_by_id(r).bounds = 0,0
               #allow oxygen and water free exchange, co2 outflow 
          model.reactions.get_by_id('EX_o2(e)').bounds= -1000,1000
          model.reactions.get_by_id('EX_h2o(e)').bounds= -1000,1000
          model.reactions.get_by_id('EX_co2(e)').bounds= 0,1000
          model.reactions.get_by_id('EX_h(e)').bounds= 0,1000 
    
          try:
              model.reactions.get_by_id('sink_octdececoa(c)').bounds = 0,0
              model.reactions.get_by_id('sink_citr(c)').bounds = 0,0
                        
          except:
              pass                 
                      
          # allow inflow of selected carbon source 
          model.reactions.get_by_id(carbon_source).bounds = -1000, 1000     
          # FBA
                   
          rq = RQ(model)
          rq_list.append(rq)
          print(carbon_source)
          print('RQ = ',rq)
          print(model.summary())
          
          print("===============================================================")
                
    return rq_list     

def integrate(y_vals, h):
    i = 1
    total = y_vals[0] + y_vals[-1]
    for y in y_vals[1:-1]:
        if i % 2 == 0:
            total += 2 * y
        else:
            total += 4 * y
        i += 1
    return total * (h / 3.0)             
            
#%%
     
path = '/home/acabbia/Documents/Muscle_Model/models'
library_folder = path + "/library_GEO_GSE25941_v6/"

RQ_prediction = DataFrame( index = ['EX_glc(e)','EX_fru(e)','EX_octa(e)','EX_ttdca(e)','EX_hdca(e)','EX_ocdca(e)',
                                    'EX_arach(e)','EX_doco13ac_','EX_lgnc(e)','EX_lnlnca(e)','EX_hdcea(e)','EX_lnlc(e)'])

for model_name in os.listdir(library_folder):
    model = cobra.io.read_sbml_model(library_folder+model_name)
    rq_list = RQ_test(model)
    RQ_prediction[model_name] = rq_list 
    
##### Swarmplot RQ predictions
g = sns.swarmplot(data = RQ_prediction.T, orient='h')
gg = g.get_figure()
gg.savefig('/home/acabbia/out/RQ_plot.png')

#%%
#### Metflex2 (PARETO OPTIMUM CURVES)
#Objective1 glucose oxidation: 
obj_1 = "HEX1"
#objective2 Fatty acid oxidation
obj_2 = "LNLCCPT2"

auc_list = []       
age_list = []

index = []

names = ["ASK_YF_14.xml","ASK_YM_36.xml"]
for model_name in os.listdir(library_folder):
    model = cobra.io.read_sbml_model(library_folder+model_name)
    
    index.append(model_name)
    if model_name[4] == 'Y':
        age_list.append('Young')
    elif model_name[4] == 'O':
        age_list.append('Old')
        
    ex = Utils.find_exchanges(model)
    #### minimal medium
    #close all exchanges
    for e in ex:
        model.reactions.get_by_id(e).bounds = 0,0
    # carbon sources
    model.reactions.get_by_id('EX_glc(e)').bounds= -100,0
    model.reactions.get_by_id('EX_lnlc(e)').bounds= -100,0
    
    #redundant reactions to be closed (unwanted loops)
    model.reactions.get_by_id('r0355').bounds = 0,0
    model.reactions.get_by_id('r0354').bounds = 0,0
    model.reactions.get_by_id('SBTR').bounds = 0,0
    model.reactions.get_by_id('LINOFATPtc').bounds = 0,0
    model.reactions.get_by_id('r0987').bounds = 0,0
    model.reactions.get_by_id('r1488').bounds = 0,0
    model.reactions.get_by_id('r0408').bounds = 0,0
    model.reactions.get_by_id('FACOAL161').bounds = 0,0
    model.reactions.get_by_id('r2511').bounds = 0,0
    model.reactions.get_by_id('r0438').bounds = 0,0
    model.reactions.get_by_id('r1318').bounds = 0,0
       
    # keep open essential exchanges
    model.reactions.get_by_id('EX_o2(e)').bounds= -1000,1000
    model.reactions.get_by_id('EX_h2o(e)').bounds= -1000,1000
    model.reactions.get_by_id('EX_co2(e)').bounds= 0,1000
    model.reactions.get_by_id('EX_h(e)').bounds= 0,1000 
    
    y_min = []
    y_max = []
    glc_in = []
    lnlc_in = []
    
    flx = DataFrame(index = np.arange(0,1100,100))
    intakes = DataFrame(index = np.arange(0,11,1))
    
    for b in np.arange(0,1100,100):
        model.objective = obj_1
        model.reactions.get_by_id(obj_2).bounds = 0,b
        #model.reactions.get_by_id(obj_3).bounds = 0, (1000 - b)
        sol_max =  model.optimize('maximize')
        model.summary()
        #model.metabolites.atp_c.summary()
        sol_min = model.optimize('minimize')
        y_min.append(abs(sol_min.objective_value))
        y_max.append(abs(sol_max.objective_value))
        glc_in.append(sol_max.fluxes['EX_glc(e)'])
        lnlc_in.append(sol_max.fluxes['EX_lnlc(e)'])
    
    intakes['glc'] = glc_in
    intakes['Lnlc'] = lnlc_in
    
 
    flx['max'] = y_max
    flx['min'] = y_min
    auc_list.append(integrate(flx['max'].values, 100))
    
    ax = flx.plot.area(title = str(model_name), legend = False, xlim = (0,1000), ylim = (0,1000))    
    ax.set(xlabel=obj_2+' flux', ylabel=obj_1+' flux')
    savefig('/home/acabbia/out/flexibility/'+model_name+'_flexi.png')
    
    ax2 = intakes.abs().plot.bar(title = str(model_name),ylim=(0,120))
    ax2.set(xlabel = 'loop nr.',ylabel = 'Intake flux of carbon source (mM/gDw/h)')
    savefig('/home/acabbia/out/flexibility/'+model_name+'_intakes.png')
    
    show()

auc = DataFrame({'AUC':auc_list, 'age': age_list},index = index)
auc['AUC']= auc['AUC']/1000000
xx = auc.groupby('age')['AUC'].mean().plot.bar(title = 'Average AUC of the two classes', ylim = (0,1))
xx.set(ylabel = 'Normalized AUC')
savefig('/home/acabbia/out/flexibility/avg_auc.png')



#%%
#### Metflex3 (individual model demands) 
from cobra.core import Reaction

obj = 'CSm'
obj_3 = "HEX1"
obj_2 = "LNLCCPT2"

auc_list = []       
age_list = []

index = []

names = ["ASK_YF_14.xml","ASK_YM_36.xml"]
for model_name in os.listdir(library_folder):
    model = cobra.io.read_sbml_model(library_folder+model_name)
    
    dm_ac = Reaction('DM_ac_c',name= 'Acetate_demand', lower_bound = 0, upper_bound = 1000)
    dm_ac.add_metabolites({model.metabolites.get_by_id('ac_c'):-1})
    model.add_reactions([dm_ac])
    
    index.append(model_name)
    if model_name[4] == 'Y':
        age_list.append('Young')
    elif model_name[4] == 'O':
        age_list.append('Old')
        
    ex = Utils.find_exchanges(model)
    #### minimal medium
    #close all exchanges
    for e in ex:
        model.reactions.get_by_id(e).bounds = 0,0
    
    #carbon sources
    model.reactions.get_by_id('EX_glc(e)').bounds= -100,0
    model.reactions.get_by_id('EX_lnlc(e)').bounds= -100,0
    
    #redundant reactions to be closed (unwanted loops)
    model.reactions.get_by_id('r0355').bounds = 0,0
    model.reactions.get_by_id('r0354').bounds = 0,0
    model.reactions.get_by_id('SBTR').bounds = 0,0
    model.reactions.get_by_id('LINOFATPtc').bounds = 0,0
    model.reactions.get_by_id('r0987').bounds = 0,0
    model.reactions.get_by_id('r1488').bounds = 0,0
    model.reactions.get_by_id('r0408').bounds = 0,0
    model.reactions.get_by_id('FACOAL161').bounds = 0,0
    model.reactions.get_by_id('r2511').bounds = 0,0
    model.reactions.get_by_id('r0438').bounds = 0,0
    model.reactions.get_by_id('r1318').bounds = 0,0
       
    # keep open essential exchanges
    model.reactions.get_by_id('EX_o2(e)').bounds= -1000,1000
    model.reactions.get_by_id('EX_h2o(e)').bounds= -1000,1000
    model.reactions.get_by_id('EX_co2(e)').bounds= 0,1000
    model.reactions.get_by_id('EX_h(e)').bounds= 0,1000 
    
    y_min = []
    y_max = []
    glc_in = []
    lnlc_in = []
    
    flx = DataFrame(index = np.arange(0,1100,100))
    intakes = DataFrame(index = np.arange(0,11,1))
    
    for b in np.arange(0,1100,100):
        model.objective = obj
        model.reactions.get_by_id(obj_2).bounds = 0,b
        model.reactions.get_by_id(obj_3).bounds = 0, (1000 - b)
        sol_max =  model.optimize('maximize')
        model.summary()
        #model.metabolites.atp_c.summary()
        sol_min = model.optimize('minimize')
        y_min.append(abs(sol_min.objective_value))
        y_max.append(abs(sol_max.objective_value))
        glc_in.append(sol_max.fluxes['EX_glc(e)'])
        lnlc_in.append(sol_max.fluxes['EX_lnlc(e)'])
    
    intakes['glc'] = glc_in
    intakes['Lnlc'] = lnlc_in
    

    ax2 = intakes.abs().plot.bar(title = str(model_name),ylim=(0,120))
    ax2.set(xlabel = 'loop nr.',ylabel = 'Intake flux of carbon source (mM/gDw/h)')
    savefig('/home/acabbia/out/flexibility/in/'+model_name+'_intakes.png')
    
    show()
























































