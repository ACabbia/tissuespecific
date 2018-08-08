#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 14:34:49 2018

@author: acabbia
"""
import cobra
import GEOparse
from tissuespecific.reconstruction import Builder 
import pandas as pd
from interruptingcow import timeout
from os import listdir
from os.path import isfile, join


def build_models(GEO_accession_nr):
    #get data from GEO
    serie=GEOparse.get_GEO(geo=GEO_accession_nr)
    #build translator dict
    table =  serie.gsms[list(serie.gsms.keys())[0]].table
    translator = Builder.affyprobe_translator(table, 'hgnc_id')
    return serie, translator

def build_loop(model, serie, translator, names, label_dict, requirements, expr_cutoff, library_folder):
    # main loop builds model for each subject in serie  (multiprocessing target)
    skipped = []
    print('===============================')
    for gsm in serie.gsms:
        with timeout(1800, exception= TimeoutError):
            try:
                if gsm in names:
                    lll = label_dict[gsm]
                    print('Building model '+ 'ASK_'+gsm+'_'+str(lll))
                    eset = serie.gsms[gsm].table
                    confidence = Builder.rxn_confidence(model, eset, translator,'hgnc_id', expr_cutoff )
                    for r_id in requirements:
                        confidence[r_id] = 3
                    newmodel = Builder.build_model(model,confidence)
                    newmodel.name = 'Skeletal Muscle '+ gsm + lll
                     
                    # Reactions to be corrected (h_i)
                    newmodel.reactions.CYOOm3.add_metabolites({'h_m':-7.9,'h_i':4 ,'h_c':0}, combine=False)
                    # Fix PDH reversibility/directionality
                    newmodel.reactions.PDHm.bounds = 0, 1000
                    #prune unused metabolites
                    remov = cobra.manipulation.delete.prune_unused_metabolites(newmodel)    
                    while len(remov) != 0:
                        remov = cobra.manipulation.delete.prune_unused_metabolites(newmodel)
                    #write newmodel in SBML
                    cobra.io.write_sbml_model(newmodel,filename= library_folder+'ASK_'+gsm+'_'+lll+'.xml')
                    print('===============================')
            except TimeoutError:
                skipped.append(gsm)
                print('Time out occurred:', gsm, 'skipped')
                continue
    print(len(skipped), 'models not built (timeout):' )
    print(skipped)
            
###############################################################################################################
#set paths
path = '/home/acabbia/Documents/Muscle_Model/models'
ref_model_path = path + '/recon2.2.xml'
library_folder =  path + '/library_GEOmerge/'

#load and prepare df
df = pd.read_csv('/home/acabbia/Documents/Muscle_Model/expressionData/mergione_batch_corrected.csv', sep = ',')
label = list(df['label'].values)
batch = list(df['batch'].values)
df.index = df['Unnamed: 0']
df = df.drop(['label','batch','Unnamed: 0'],axis= 1).T

# list of GSEs to be converted in models
names= list(df.columns)
names = [ n.split('_')[0] for n in names]
names = [ n.split('.')[0] for n in names]

names.remove('GSM702445') # <--- remove this GSE from the list of models to make (it hangs)
lbl_dict = dict(zip(names,label))

#Avoid build models who have already been built
file_list = [f.split('_')[1] for f in listdir(library_folder) if isfile(join(library_folder, f))]
names = list(set(names) - set(file_list))

#Gene expression data GEO ID (Raue 2012)
GEO_accession_list = ['GSE28422']
# ['GSE9103','GSE47881',,'GSE25941']

# import reference model (RECON2.2)
model = cobra.io.read_sbml_model(ref_model_path)
# required reactions  (mostly exchanges)
requirements = ['ATPS4m','ENO','PDHm','PYK','G3PD1','G6PDH2r','AKGDm','CYOOm3','r0913','FAOXC4020m','FAOXC80','EX_glc(e)',
                 'EX_hdca(e)','EX_ocdca(e)','EX_arach(e)','EX_doco13ac_','EX_lgnc(e)','EX_ala_L(e)',
                 'EX_arg_L(e)','EX_asn_L(e)','EX_asp_L(e)','EX_cys_L(e)','EX_gln_L(e)','EX_glu_L(e)',
                 'EX_gly(e)','EX_his_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)','EX_met_L(e)',
                 'EX_phe_L(e)','EX_pro_L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)',
                 'EX_val_L(e)','EX_fru(e)','EX_ppa(e)','EX_but(e)','EX_hx(e)','EX_octa(e)','EX_ttdca(e)',
                 'EX_h2o(e)','EX_hco3(e)','EX_co2(e)','EX_h(e)']

#%%
# model building loop
for GEO in GEO_accession_list:
   serie, translator = build_models(GEO)
   build_loop(model, serie, translator, names, lbl_dict, requirements, 0.7, library_folder)
    

    
    
    
    