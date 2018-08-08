#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 11:48:35 2017

@author: acabbia
"""
import cobra
import csv
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
import os
import matplotlib.pyplot as plt
import re

class Utils:
    #load bound values for uptake fluxes from .tsv
        
     def load_bounds(path):
         intake_bounds = pd.read_csv(path, sep = '\t', header = None)
         intake_bounds.columns = ['Reaction','Flux Value']
         r_id = [s.replace('[e]','_e')for s in intake_bounds.Reaction.values]
         norm = [((n - intake_bounds['Flux Value'].min())/ (intake_bounds['Flux Value'].max() - intake_bounds['Flux Value'].min())) 
                 for n in intake_bounds['Flux Value'].values]
         bounds_dict = dict(zip(r_id, norm))
         return bounds_dict
     
     def apply_bounds(model,bounds_dict):
        for r in model.exchanges:
            r.bounds = 0,0
        for k, v in bounds_dict.items():
            try:
                model.reactions.get_by_id(k).bounds = -v,0
            except:
                print(k, 'not in model')
                
    # list of model exchange reactions (no sink/demands)
     def find_exchanges(model):
         exchanges = []
         pattern = re.compile('EX_*')
         for r in model.exchanges:
             if pattern.match(r.id):
                 exchanges.append(r.id)
         return exchanges
                 
     ## Tries to optimize all the reactions in the model 
     def optimize_all_rxns(model, print_summary):
          opt=0
          for r in model.reactions:
               model.objective = str(r.id)
               sol=model.optimize()
          if sol.objective_value != 0:
               opt += 1
               print("optimizing for reaction: ", r.id , " Obj. value:", sol.objective_value )
               if print_summary == True:
                    print(model.summary())
                    print('-------------')
                    print('# of non-zero solutions:', opt, '/', len(model.reactions))
          
     ## finds internal blocked reactions 
     def find_blocked(model):
          blocked=cobra.flux_analysis.variability.find_blocked_reactions(model, zero_cutoff=1e-9)
          x=[]
          for b in blocked:
               r=model.reactions.get_by_id(b)
          if r not in model.exchanges:
               x.append(b)    
               included=[r.id for r in model.reactions]
         
          not_blocked=set(included)-set(blocked)    
          return blocked , not_blocked

      ## write file for MINERVA RECON MAP visualization
     def write_viz_file(filename, blocked , not_blocked):
          with open(filename, 'w') as csv_file:
               writer = csv.writer(csv_file, delimiter='\t')
               writer.writerow(['reactionIdentifier','lineWidth','color'])
               for b in blocked:
                    writer.writerow([b,4,'#FF0000']) # red
               for n in not_blocked:
                    writer.writerow([n,3,'#800080']) # green
   
     def nonzero_fluxes(fba_solution):
         '''
         returns pandas.dataframe containing fluxes > 1e-9 in the solution
         '''
         sol = fba_solution.fluxes[fba_solution.fluxes > 1e-9]
         return sol
     
     def normalize_fluxes(fba_solution):
         min_flux = min(fba_solution)
         max_flux = max(fba_solution)
         norm_fluxes = (fba_solution.fluxes - min_flux) / (max_flux - min_flux)
         return norm_fluxes
             
     def viz(solution, filename):
         df = pd.DataFrame(columns = ['reactionIdentifier','lineWidth','color'])
         sol = Utils.nonzero_fluxes(solution)
         
         df['reactionIdentifier'] = sol.index
         df['lineWidth'] = sol.values
         df['color'] = '#FF0000'# red
         
         df.to_csv(filename, sep= '\t', header = True, index = False)
         return df
            
     def plot_sampl(to_plot, *sampl):
          # plot result of MC Sampling experiments
          for r in to_plot:
               join= pd.concat([sampl[r]], axis=1)
               ax = join.plot.hist(bins=50, title = "Reaction "+r )
               ax.get_figure()
               plt.show()
               ax2 = join.plot.box(title = "Reaction "+r )
               ax2.get_figure()
               plt.show()

     def plot_2(to_plot, samplingr2, samplingSKM):
         #plot results of MC sampling (R2 and SKM specific)
         for r in to_plot:
             join = pd.concat([samplingr2[r], samplingSKM[r]], axis=1)
             join.columns = ['Recon2.2', 'SKM_old']
             ax = join.plot.hist(bins=100, title = "Reaction "+r )
             ax.get_figure()
             plt.show()
             ax2 = join.plot.box(title = "Reaction "+r )
             ax2.get_figure()
             plt.show()
                                     
class PCA_utils:
     
     
     def load_sample(filepath, label): 
          ##load saved CSVs & add rename label var
          sample=pd.read_csv(filepath, sep='\t')
          sample.rename(columns={"Unnamed: 0": "label"},inplace = True)
          sample['label']= label
          return sample

     def load_all_csv(path):
          result = pd.concat((PCA_utils.load_sample(os.path.join(path,f),str(f)) for f in os.listdir(path)), axis=0, join='outer')
          return result
     
     def PCA_preprocess(result): 
          ## join dataframes, replace NaN with 0, and scales data
          result= result.fillna(0)
          labels = result['label']
          rx=result.loc[:, result.columns != 'label']
          rx=scale(rx,axis=1)
          return labels, rx

     def samples_PCA(data):
          pca=PCA(n_components=3)
          pca.fit(data)
          print("Percentage of variance explained by the first 3 components:", sum(pca.explained_variance_ratio_))
          transformed=pca.transform(data)
          return pd.DataFrame(transformed)

     def load_and_merge(path):
          r=PCA_utils.load_all_csv(path)
          lbl,rx=PCA_utils.PCA_preprocess(r)
          merged=pd.concat((lbl,r),axis=1)
          merged.to_csv('sampling_merged.csv',sep='\t')
          return merged

