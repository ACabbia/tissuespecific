#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:14:14 2017

@author: acabbia
"""

def configure_objective(model, prot_synth, converted):
     ## configure objective function
     c = []     
     p = []
     for r in converted:
          c.append(model.reactions.get_by_id(r).flux_expression)
     for s in prot_synth:
          p.append(model.reactions.get_by_id(s).flux_expression)
     
     difference = model.problem.Variable('difference')
     constraint = model.problem.Constraint(sum(p)-sum(c)-difference,lb =0, ub=0)
     model.add_cons_vars([difference, constraint])

