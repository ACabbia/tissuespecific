#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:14:25 2017

@author: acabbia
"""
from tissuespecific.reconstruction.analysis import Bounds
from tissuespecific.sarcoModel.estimate import estimate_TEE
                    
INF = float('inf')

def runSimulation(Individual, Settings, Food):
     # init data collection list
     m_mass =   [0] * (Settings.simLength +1)
     TEE =      [0] * (Settings.simLength +1)
     atp_exp =  [0] * (Settings.simLength +1)
     obj_val =  [0] * (Settings.simLength +1)
               
     weight = Individual.weight # weight at t=0
     m_mass[0] = (weight * (1-Individual.fat_ratio)) ## weight of muscle mass (Kg)
     
     model=Settings.model 
     ex_list = Bounds.select_exchanges(model)
     
     #simulation loop
     for i in list(range(Settings.simLength)):
          
          print(i)
          
          #reset uptake bounds
          for r in model.exchanges:
               r.bounds = 0, INF        
          
          #set new inbound fluxes
          
          bounds_dict = {k: (( v *  (1-Individual.fat_ratio)) * Individual.activityFactor * Individual.uptakeFactor * Settings.timeStep)
                             for k,v in Food.uptake_bounds.items()
                             }
          
          Bounds.change_exch_bounds(model, ex_list, bounds_dict, keep_open=['EX_h2o(e)','EX_o2(e)','EX_tRNA'])
          
          #set atp expenditure
          FFM = m_mass[i]
          print("Fat-free mass:", FFM)
          
          total_expenditure = estimate_TEE(FFM, Individual)
          TEE[i] = total_expenditure
          print("estimated total expenditure (Kcal/day)", total_expenditure)
          
          # normalize by dry weight
          atp_expenditure = total_expenditure / Individual.ATPconversion ## from kcal/day to mol ATP/day
          #atp_expenditure = atp_expenditure/ ((Individual.weight * (1-Individual.fat_ratio)) * Individual.proteinLeanFactor )
          atp_flux = atp_expenditure * Settings.timeStep
          atp_exp[i] = atp_flux
          print('ATP_expenditure flux:',atp_flux)
          
          #constrain ATP expenditure reaction
          model.reactions.ATP_expenditure.bounds = atp_flux , atp_flux  
          
          
          model.objective = {model.reactions.get_by_id('contr_unit_1_storage'): 1, 
                             model.reactions.get_by_id('contr_unit_2a_storage'): 1,
                             model.reactions.get_by_id('contr_unit_2b_storage'): 1,
                             model.reactions.get_by_id('contr_unit_2x_storage'): 1}
          
          #optimize
          sol = model.optimize()
          
          print("obj.value:", sol.objective_value)
          print('fwd-bwd flux diff:', model.solver.variables.difference.primal)
          obj_val[i] = sol.objective_value
          
          
          
          m_mass[i+1] = m_mass[i] + (model.solver.variables.difference.primal * 1/Individual.proteinLeanFactor * 1000000) 
          print('--------------------------------')
     
     #print(model.summary())     
     return m_mass , atp_exp, obj_val , TEE
     