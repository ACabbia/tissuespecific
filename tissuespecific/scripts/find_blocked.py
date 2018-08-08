#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 15:41:42 2018

@author: acabbia
"""
import cobra
import os

path = '/home/acabbia/Documents/Muscle_Model/models'
library_folder = path + "/library_GEO_GSE25941_v4/"

blocked = []

for filename in os.listdir(library_folder):
    model = cobra.io.read_sbml_model(library_folder+filename)
    bl = cobra.flux_analysis.find_blocked_reactions(model,open_exchanges= True)
    blocked.append((filename , bl))            
