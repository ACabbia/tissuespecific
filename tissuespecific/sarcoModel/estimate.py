#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 15:58:54 2017

@author: acabbia
"""

def estimate_muscle_maintenance_cost(Individual):
     '''
     Individual = object of class 'Individual'
     activityCost = cost of physical activity in kcal/kg
     '''
     FFM = (Individual.weight * (1 - Individual.fat_ratio)) ### Fat-Free Mass
     maintenance_cost = Individual.leanMass_maintenance * FFM
     return maintenance_cost 

def estimate_RMR(Individual):
     #RMR [kJ/d]= 3169 + 50.0 x body weight [kg] - 15.3 x age [y] + 746 x sex [female = 0, male = 1]
     eRMR = 3169 + 50*(Individual.weight) - 15.3 * Individual.age + 746 * int(Individual.isMale)
     return eRMR

def estimate_BMR(FFM): # Grande and Keys eqn. (not elderly-specific)
     '''
     Grande F, Keys A. Body weight, body composition, and calorie status. 
     Modern nutrition in health and disease, 27, 1980.
     '''
     eBMR = 31.2 * FFM
     return eBMR

def estimate_TEE(FFM, Individual):
     tee = estimate_BMR(FFM) * Individual.activityFactor
     return tee

def estimate_daily_synthetic_rate(FFM, Settings):
     eDSR = (0.075 * FFM) # fractional synthetic rate (hourly)
     eDSR = eDSR * 24 * Settings.timeStep # normalize to obtain synthetic rate over the timestep
     return eDSR