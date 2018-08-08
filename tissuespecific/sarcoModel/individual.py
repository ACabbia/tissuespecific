#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 15:01:39 2017

@author: acabbia
"""

class Individual(object):
     
     '''
     class to store input data about a single simulated individual
     
     '''
     proteinLeanFactor = 0.23  ## % proteins in the muscle (by weight)
     ATPconversion = 20 # kcal / Mol ATP
     def __init__(self, name , isMale , age , weight):
          '''
          name = Individual ID
          isMale = bool
          age = age in years
          weight = weight in kg
          
          returns an object of class Individual 
          '''
          self.name = name
          self.isMale = bool(isMale)
          self.age = age
          self.weight = weight
          
     def set_fat_ratio(self, fat_ratio): 
          # sets fat lean ratio 
          self.fat_ratio = fat_ratio
                    
     def set_uptakeFactor(self, uptakeFactor):
          # % of nutrients adsorbed from diet (function of age)
          self.uptakeFactor = uptakeFactor
          
     def set_leanMass_maintenance_cost(self, lean_maintenance):
          # set lean mass maintenance cost (kcal/g)
          self.leanMass_maintenance = lean_maintenance
          
     def set_activityfactor(self, activityFactor): ## TO DO : accept input as category (i.e sedentary/ active...) or directly as factor
          '''
          Sedentary (little or no exercise, desk job).

          BMR x 1.2

          Lightly Active (light exercise/sports 3-5 days/week).

          BMR x 1.3-1.4

          Moderately Active (moderate exercise/sports 3-5 days/week).

          BMR x 1.5-1.6

          Very Active (hard exercise/sports 6-7 days per week).

          BMR x 1.7-1.8

          Extremely Active (very hard daily exercise/sports and physical job or 2/day training).

          BMR x 1.9-2.0
          '''
          self.activityFactor = activityFactor
          
     
class Settings(object):
      '''
      class to store simulation settings
      '''
      def __init__ (self, model, simLength, timeStep):
          '''
          model = a COBRA model structiure
          simLength = length of the sim in n(timesteps)
          timeStep = simulation increments in days
          proteins_MW = vector of molecular weights of the muscle proteins 
          '''
          self.model = model
          self.simLength = simLength
          self.timeStep = timeStep
          
      def proteins_MW(self,MW):
          self.proteins_MW = sum(MW)
           
class Food(object):
     '''
     class to normalize and store data about nutrients uptake bounds
     '''
       
     def __init__(self, daily_bounds_dict, Settings, Individual):
          '''
          daily_bounds_dict = dict {EX_reaction_name : flux value}
          Settings = Settings object
          Individual = Individual object
          
          returns:
               normalized uptake bounds over the time step
          '''
          
          self.uptake_bounds = { k : Settings.timeStep * (v * (1 - Individual.fat_ratio)) for k , v in daily_bounds_dict.items()}
     
     
     def load_bounds_file(path):
          import pandas as pd
          intake_bounds = pd.read_csv(path, index_col = 0)
          daily_bounds_dict= dict(zip(intake_bounds.reaction, intake_bounds.fluxValue))
          return daily_bounds_dict

