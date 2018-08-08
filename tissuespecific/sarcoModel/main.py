#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:54:03 2017

@author: acabbia
"""
''''
main simulation file
'''
import tissuespecific
from tissuespecific.sarcoModel.individual import Individual , Settings , Food 
from cobra.io import read_sbml_model
from tissuespecific.reconstruction.protein_metabolism import Protein_metabolism
from tissuespecific.sarcoModel.sim import runSimulation 

#load muscle model and add protein metabolism 

model = read_sbml_model('/home/acabbia/Documents/Muscle_Model/models/SKM_old_latest.xml')
#Protein_metabolism.add_prot_metabolism(model)  <- already added

#make 'individual' object to store subject data
# male, age 78 , 62 kg

old = Individual(name = 'old',
                 isMale = True, 
                 age = 78,
                 weight = 62)

old.set_fat_ratio (.26)  # 25% body fat
old.set_leanMass_maintenance_cost(60) # kcal/g
old.set_uptakeFactor(.9) # 90% of nutrients absorbed
old.set_activityfactor(1.1)

#configure simulation settings
#model = skeletal muscle , timestep = 1 year , simulation = 20 years

settings = Settings ( model = model,
                     timeStep = 30,
                     simLength = 12,
                     parsimonious = False)

#load and normalize daily uptake bounds from food intake

daily_bounds = Food.load_bounds_file('/home/acabbia/Documents/Muscle_Model/tissuespecific/tissuespecific/reconstruction/bounds_EU_AVG.tsv')

food = Food(daily_bounds, settings, old)

########################################################

results = runSimulation(old,settings,food)
