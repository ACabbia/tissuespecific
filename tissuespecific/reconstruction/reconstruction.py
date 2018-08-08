#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 11:25:40 2017

@author: acabbia
"""
import corda
import pandas as pd
import biomart
from cobra import Reaction, Metabolite
import cobrababel
import numpy as np
from interruptingcow import timeout


class Builder:
    
     def affyprobe_translator(table, geneID):
         print('Building translation dictionary ... ')
         # Creates mapping affyprobe ID <-> geneID
         # returns dict{affy_probe ID: geneID}
         # allowed geneID formats: hgnc_id (recon2.2), ENSG/ensemble_gene_id (HMR2), CCDS (Recon3D)
         
         valid_translators = {'hgnc_id', 'hgnc_symbol','ccds'}
         if geneID not in valid_translators:
             raise ValueError("translate_to must be one of: ", valid_translators)
         
         table=table.rename(columns={ table.columns[0]: "probes" })
         table=table.rename(columns={ table.columns[1]: "value" })
         
         probes=[]
         for i, p in table.iterrows():
               probes.append(table.at[i,"probes"])
               
         server=biomart.BiomartServer('http://www.ensembl.org/biomart')  
         dataset=server.datasets['hsapiens_gene_ensembl']
         response=dataset.search({'attributes':['affy_hg_u133_plus_2',geneID]})
         translator_dict={}
         for line in response.iter_lines():
               line = line.decode('utf-8')
               (a,b)=(line.split("\t"))
               if geneID=='ccds':
                   translator_dict[a] = b[4:]
               else:
                   translator_dict[a]=b
         return translator_dict
               
     def rxn_confidence(model, table, translator , gene_id, q):
          # Creates confidence score for each reaction of the model           
          # Define only high confidence reaction (score = 3)
          confidence = {r.id: r.gene_reaction_rule for r in model.reactions}
          # Translate probe id->gene ID
          table=table.rename(columns={ table.columns[0]: "probes" })
          table=table.rename(columns={ table.columns[1]: "value" })
          translator_df=pd.DataFrame.from_dict(translator, orient='index').reset_index()
          translator_df.columns = ['probes', 'hgnc_id']
          #print(table.columns, translator_df.columns)
          translator_df=pd.merge(translator_df,table,on='probes',how='inner')
          tr = translator_df[['hgnc_id','value']]
          # Average duplicated genes (avereps)
          trr = tr.groupby('hgnc_id').mean().reset_index()
          # log normalize and assign confidence score
          trr['value'] = trr['value'].transform(np.log1p)
          trr['confidence']=-1
          idx = [trr['value'] > trr['value'].quantile(q)]
          trr.loc[idx[0],'confidence'] = 3
         
          gene_confidence=pd.Series(trr.confidence.values,index=trr.iloc[:,0].values).to_dict()
         
          ## Evaluate gprs
          for k , v  in confidence.items():
               if isinstance(v, str):
                    confidence[k]=(corda.util.reaction_confidence(v,gene_confidence))
               else:
                    confidence[k]=-1  ## Reaction with no associated GPR rule have confidence score = -1
          for k,v in confidence.items():
              if v == 0:
                  confidence[k]=-1  ### all reactions except high confidence ones have score = -1 
          return confidence 

     def build_model(model, confidence):
          with timeout(1800, exception= TimeoutError):
            try:
              opt = corda.CORDA(model, confidence , n=3)
              opt.tflux = 1e-3
              opt.build()
              print(opt)
              newmodel=opt.cobra_model()        
              return newmodel
            except TimeoutError:
                raise TimeoutError

     def add_bigg_reactions(model,BiGG_id_list):
          for r in BiGG_id_list:
               d=cobrababel.get_bigg_reaction(r, model_bigg_id='universal')
               
               reaction = Reaction(id=d['bigg_id'].replace('_e','(e)'),
                                   name = d['name'],
                                              subsystem = '',
                                              lower_bound = -1000 ,
                                              upper_bound =  1000 ,
                                              objective_coefficient=0.0)
               
          for k, val in enumerate(d['metabolites']):
                    met=Metabolite(val['bigg_id'].replace('__', '_')+'_'+val['compartment_bigg_id'])
                    met.name=val['name']
                    met.compartment=val['compartment_bigg_id']
                    reaction.add_metabolites({met:val['stoichiometry']})

          if reaction not in model.reactions:
                    model.add_reaction(reaction)
          
          return model
      

