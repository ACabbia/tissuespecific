#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import corda
import pandas as pd
import biomart
import numpy as np

class Builder:
    
     def affyprobe_translator(table, geneID):
         print('Building translation dictionary ... ')
         # Creates mapping affyprobe ID <-> geneID
         # returns dict{affy_probe ID: geneID}
         # allowed geneID formats: hgnc_id (recon2.2), ENSG/ensemble_gene_id (HMR2), CCDS (Recon3D)
         
         valid_translators = {'hgnc_id', 'hgnc_symbol','ensembl_id,','ccds'}
         if geneID not in valid_translators:
             raise ValueError("geneID must be one of: ", valid_translators)
         table.columns=['probes', 'value', 'PA_call', 'P-val']
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
               translator_dict[a]=b
               
         return translator_dict
               
     def rxn_confidence(model, table, translator , gene_id):
          # creates confidence score for each reaction of the model           
          # Only high confidence reaction have score = 3, all the remaining reactions have score = -1
          
          confidence = {r.id: r.gene_reaction_rule for r in model.reactions}
          # Normalize expression value (log normalization)
          table.columns=['probes', 'value', 'PA_call', 'P-val']
          table['value'] = table['value'].transform(np.log)
          table['confidence']=-1
          # Here I consider the top 25% expressed probes as 'High Confidence' to avoid choosing a specific threshold value
          idx = [table['value'] > table['value'].quantile(.75)]
          table.loc[idx[0],'confidence'] = 3
          
          translator_df=pd.DataFrame.from_dict(translator, orient='index')
          translator_df[gene_id] = translator_df.index
          translator_df.columns=[gene_id,'probes']
          translator_df=translator_df.drop_duplicates(subset=gene_id)
          table=pd.merge(translator_df,table,on='probes',how='inner')
          table= table.dropna(axis = 0, how = 'any')
          gene_confidence=pd.Series(table.confidence.values,index=table.iloc[:,0].values).to_dict()
         
          #map probe-level confidence on GPR rules to get reaction-level confidence
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
          opt = corda.CORDA(model, confidence)
          opt.build() 
          print(opt)
          newmodel=opt.cobra_model()
                   
          return newmodel
      
#%%
import cobra
import GEOparse
          
# import reference model (RECON2.2)
recon22 = cobra.io.read_sbml_model('/home/acabbia/Downloads/MODEL1603150001.xml') 
#Gene expression data GEO ID (Raue 2012)
GEO_accession_nr = "GSE25941"
#get data from GEO
GEOseries=GEOparse.get_GEO(geo=GEO_accession_nr)
gsm = 'GSM637513'
# required reactions  
requirements = ['ATPS4m','ENO','PDHm','PYK','G3PD1','G6PDH2r','AKGDm','CYOOm3','r0913','GLCt2_2','EX_glc(e)',
                         'EX_hdca(e)','EX_ocdca(e)','EX_arach(e)','EX_doco13ac_','EX_lgnc(e)','EX_ala_L(e)',
                         'EX_arg_L(e)','EX_asn_L(e)','EX_asp_L(e)','EX_cys_L(e)','EX_gln_L(e)','EX_glu_L(e)',
                         'EX_gly(e)','EX_his_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)','EX_met_L(e)',
                         'EX_phe_L(e)','EX_pro_L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)',
                         'EX_val_L(e)','EX_fru(e)','EX_ppa(e)','EX_but(e)','EX_hx(e)','EX_octa(e)','EX_ttdca(e)']

#build translator dict (affy<->geneID)
table =  GEOseries.gsms[gsm].table
translator = Builder.affyprobe_translator(table, 'hgnc_id') #GPR rules in recon2.2 are in HGNC ID format
#build confidence dict
confidence = Builder.rxn_confidence(recon22, table, translator,'hgnc_id')
# to include required reactions, change their confidence score to 3
for r_id in requirements:
    confidence[r_id] = 3

print('building CORDA model...')    
newmodel = Builder.build_model(recon22,confidence) # <- doesn't end
        
        
        
        
        
        
        
        
        
        
        