#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:52:40 2018
@author: acabbia
"""

import cobra
import GEOparse
from tissuespecific.reconstruction import Builder 

#set paths
path = '/home/acabbia/Documents/Muscle_Model/models'
ref_model_path = path + '/recon2.2.xml'
library_folder =  path + '/library_GEOmerge/'
#Gene expression data GEO ID (Raue 2012)
GEO_accession_list = ['GSE9103','GSE25941','GSE28422','GSE47881']

def  build_models_library(ref_model_path, GEO_accession_nr,requirements, expr_cutoff, library_folder):
    # import reference model (RECON2.2)
    recon22 = cobra.io.read_sbml_model(ref_model_path)
    #get data from GEO
    serie=GEOparse.get_GEO(geo=GEO_accession_nr)
    #build translator dict
    table =  serie.gsms[list(serie.gsms.keys())[0]].table
    translator = Builder.affyprobe_translator(table, 'hgnc_id')
      

    #################################################################################################################
    # main loop builds model for each subject in serie
    print('===============================')
    for gsm in serie.gsms:
        label = str(serie.gsms[gsm].metadata['title']).split('_')
        print('Building model '+ 'ASK_'+label[1][0]+label[2][0]+'_'+gsm[-2:])
        eset = serie.gsms[gsm].table
        confidence = Builder.rxn_confidence(recon22, eset, translator,'hgnc_id', expr_cutoff )
        for r_id in requirements:
            confidence[r_id] = 3
        newmodel = Builder.build_model(recon22,confidence)
        newmodel.name = 'Aging Skeletal Muscle '+ gsm + label[1]+label[2]
         
        # Reactions to be corrected (h_i)
        newmodel.reactions.CYOOm3.add_metabolites({'h_m':-7.9,'h_i':4 ,'h_c':0}, combine=False)
        # Fix PDH reversibility/directionality
        newmodel.reactions.PDHm.bounds = 0, 1000
        #prune unused metabolites
        remov = cobra.manipulation.delete.prune_unused_metabolites(newmodel)    
        while len(remov) != 0:
            remov = cobra.manipulation.delete.prune_unused_metabolites(newmodel)
        #write newmodel in SBML
        cobra.io.write_sbml_model(newmodel,filename= library_folder+'ASK_'+label[1][0]+label[2][0]+'_'+gsm[-2:]+'.xml')
        print('===============================')
    
    ###############################################################################################################
