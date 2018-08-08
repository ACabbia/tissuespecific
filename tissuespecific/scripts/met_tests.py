#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 13:16:35 2017

@author: acabbia
"""
import cobra
from cobra.io import read_sbml_model
from cobra.core import Reaction 
from tissuespecific.reconstruction import Builder , Utils , Bounds
import pandas as pd

SKM= cobra.io.read_sbml_model('SKM_v1.1.xml')
recon2 = cobra.io.read_sbml_model('/home/acabbia/Documents/Muscle_Model/models/recon2.2.xml')
# test if the model can grow with minimal inputs (carbon source + oxygen)
INF = float('inf')

#%%
def max_flux(model, carbon_source, objective, normoxic):
          #set infinite bounds
          for r in model.exchanges:
               r.bounds = 0,INF
          #define carbon source
          model.reactions.get_by_id(carbon_source).bounds = -INF , INF     
          #set oxigen exchange
          if normoxic:
              model.reactions.get_by_id('EX_o2(e)').bounds= -INF, INF
          else:
              model.reactions.get_by_id('EX_o2(e)').bounds= 0 , 0
          #change objective
          model.objective = objective
          sol = model.optimize()
          #model.summary()
          return sol.objective_value 

def max_fluxes(model,objective):
          '''
          model = COBRA model
          objective = string, BiGG/COBRA-compliant reaction ID
          '''
    
          if objective not in model.reactions:
               Builder.add_bigg_reactions(model,[objective])
    
          print ('')

          for normoxic in [True, False]:
               for carbon_source in [
                # sugars
                'EX_glc(e)',
                'EX_fru(e)',
                # fatty acids
                'EX_ppa(e)',        # C3:0
                'EX_but(e)',        # C4:0
                'EX_hx(e)',         # C6:0
                'EX_octa(e)',       # C8:0
                'EX_ttdca(e)',      # C14:0
                'EX_hdca(e)',       # C16:0
                'EX_ocdca(e)',      # C18:0
                'EX_arach(e)',      # C20:0
                'EX_doco13ac_',     # C22:0
                'EX_lgnc(e)',       # C24:0
                # amino acids
                'EX_ala_L(e)',
                'EX_arg_L(e)',
                'EX_asn_L(e)',
                'EX_asp_L(e)',
                'EX_cys_L(e)',
                'EX_gln_L(e)',
                'EX_glu_L(e)',
                'EX_gly(e)',
                'EX_his_L(e)',
                'EX_ile_L(e)',
                'EX_leu_L(e)',
                'EX_lys_L(e)',
                'EX_met_L(e)',
                'EX_phe_L(e)',
                'EX_pro_L(e)',
                'EX_ser_L(e)',
                'EX_thr_L(e)',
                'EX_trp_L(e)',
                'EX_tyr_L(e)',
                'EX_val_L(e)',
                ]:
                    f_opt = max_flux(model, carbon_source, objective, normoxic)
                    print ('%s (%s):\t %g' % (carbon_source,'normoxic' if normoxic else 'anaerobic',f_opt))
                    
#%%
print('Number of reactions:', len(SKM.reactions))
print('Number of metabolites:', len(SKM.metabolites))
print('Number of genes:', len(SKM.genes))

from cobra.flux_analysis.sampling import sample

s_skm = sample(SKM, 5000 , processes=4)
s_r22 = sample(recon2, 5000, processes = 4)

to_plot = ['ATPS4m', # ATP synthase
           'ENO',   # Enolase
           'PDHm',  # Pyruvate dehidrogenase (mito)
           'PYK',   # Pyruvate Kinase
           'G3PD1',     
           'G6PDH2r', # Glucose 6-phosphate dehydrogenase (pentose phosphate pathway)
           'AKGDm',
           'CYOOm3', # Cytochrome c (complex IV)
           'G6PI',
           'G6PI3',
           'GLCt1r']

Utils.plot_2(to_plot, s_r22,s_skm)
