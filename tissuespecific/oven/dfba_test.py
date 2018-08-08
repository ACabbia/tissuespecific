#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 16:29:09 2018
@author: acabbia
"""
from sbmlutils.dfba.model import DFBAModel
from cobra.io import read_sbml_model

path = '/home/acabbia/Documents/Muscle_Model/models/SKM_v1.1.xml'
skm = read_sbml_model(path)

s = DFBAModel(path)
