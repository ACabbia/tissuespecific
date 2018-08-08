#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 14:45:55 2017

@author: acabbia
"""

#script complete model

import cobra
from cobra.io import read_sbml_model
from cobra.core import Reaction 
from tissuespecific.reconstruction import Builder , Utils , Bounds
import pandas as pd

INF = float('inf')
SKM= read_sbml_model('SKM_V1.xml')
#%%
# reactions to be added

missing = ['EX_thm_e', 'EX_nac_e', 'EX_na1_e', 'EX_but_e', 'EX_cu2_e', 'EX_zn2_e','EX_pydx_e',
           'EX_lcts_e', 'EX_k_e', 'EX_starch1200_e','EX_mnl_e','EX_caro_e','EX_10fthf_e',
           'EX_pydxn_e', 'EX_cl_e', 'EX_adocbl_e','EX_mn2_e','EX_pydam_e','EX_thf_e',
           'EX_retinol_e','EX_ncam_e','EX_ca2_e','EX_btn_e','EX_sucr_e','EX_phllqne_e',
           'SUCRe','EX_i_e','EX_mg2_e','EX_malt_e','EX_avite1_e','LEUtec','ACOAD8m','MCCCrm',
           'ATPS4m', # ATP synthase
           'ENO',   # Enolase
           'PDHm',  # Pyruvate dehidrogenase (mito)
           'PYK',   # Pyruvate Kinase
           'G3PD1',     
           'G6PDH2r', # Glucose 6-phosphate dehydrogenase (pentose phosphate pathway)
           'AKGDm',
           'CYOOm3', # Cytochrome c (complex IV)
           'r0913',
           'G6PI',
           'G6PI3',
           'GLCt1r'] # Glucose 6 phosphate isomerase]

for r in missing:
     Builder.add_bigg_reactions(SKM,[r])
     SKM.reactions.get_by_id(r).bounds = -1000, 1000
     print(r, "added")
#%%         
# not present in Bigg, to be added manually :

r0509 = Reaction ( id = 'r0509',
                  name = 'Succinate:ubiquinone oxidoreductase Citrate cycle (TCA cycle) EC:1.3.5.1',
                  subsystem = 'TCA cycle / electron transport chain', 
                  lower_bound = 0 , 
                  upper_bound = 1000)

SKM.add_reaction(r0509)                  
r0509.build_reaction_from_string('succ_m + q10_m --> fum_m+ q10h2_m',verbose = True)
#
r0655 = Reaction ( id = 'r0655',
                  name = '3-Methylbutanoyl-CoA:(acceptor) 2,3-oxidoreductase; EC:1.3.99.10',
                  subsystem = 'Valine, leucine and isoleucine degradation', 
                  lower_bound = -1000, 
                  upper_bound = 1000)

SKM.add_reaction(r0655)                  
r0655.build_reaction_from_string('ivcoa_m + q10_m <=> 3mb2coa_m + q10h2_m ',verbose = True)
# Reactions to be removed    
to_remove = ['r1109','PYRt2r']
SKM.remove_reactions(to_remove, remove_orphans=True)

# reactions to be corrected
SKM.reactions.CYOOm3.add_metabolites({'h_m':-7.9,'h_i':4 ,'h_c':0}, combine=False)

# Fix reaction reversibility/directionality
SKM.reactions.PDHm.bounds = 0, 1000

# Prune unused metabolites

#TO DO fix loop
while len(remov) != 0:
    remov = cobra.manipulation.delete.prune_unused_metabolites(SKM)

#%%


####load bounds
intake_bounds = pd.read_csv(bounds_path, index_col = 0)
daily_bounds_dict= dict(zip(intake_bounds.index, intake_bounds.fluxValue))
ex = Bounds.select_exchanges(SKM)
Bounds.change_exch_bounds(SKM,ex,daily_bounds_dict, ['EX_o2(e)'])

#Viz
SKM.objective = 'ATPS4m'
sol = SKM.optimize()

SKM.summary()

sol.to_frame()
Utils.viz(sol, 'viz5', NZ = True)
