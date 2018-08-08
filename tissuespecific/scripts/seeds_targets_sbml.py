#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 16:43:12 2018

@author: acabbia
"""
import pandas as pd
from cobra.core import Model , Metabolite
from cobra.io import read_sbml_model , write_sbml_model

# create model seeds SBML from diet input list
model= read_sbml_model('SKM_v1.1.xml')

bounds_path = '/home/acabbia/Documents/Muscle_Model/tissuespecific/tissuespecific/input_lists/bounds_EU_AVG.tsv'
intake_bounds = pd.read_csv(bounds_path, sep = ' ')
seeds = Model()
for r in intake_bounds.reaction:
    seeds.add_metabolites(model.reactions.get_by_id(r).reactants)
   
write_sbml_model(seeds, 'seeds.xml')
   
# create model targets SBML from diet input list
model = read_sbml_model('skm_prot.xml')
target = Model()
target.add_reaction(model.reactions.myofiber_storage)
write_sbml_model(target, 'target.xml') 