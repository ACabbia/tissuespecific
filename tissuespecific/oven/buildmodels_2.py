#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 15:22:42 2018

@author: acabbia

"""

import cobra
from cobra.core import Reaction 
import GEOparse
from tissuespecific.reconstruction import Builder 

#set paths
path = '/home/acabbia/Documents/Muscle_Model/models'
ref_model = path + "/recon2.2.xml"
output_folder = path + "/library_GEO_GSE25941_v2/"

# import reference model (RECON2.2)
recon22 = cobra.io.read_sbml_model(ref_model)

#Gene expression data GEO ID (Raue 2012)
GEO_accession_nr = "GSE25941"

#get data from GEO
serie=GEOparse.get_GEO(geo=GEO_accession_nr)

gsm = 'GSM637527'
#build translator dict
table =  serie.gsms[gsm].table
translator = Builder.affyprobe_translator(table, 'hgnc_id')
#Build confidence dict
confidence = Builder.rxn_confidence_2(recon22,table, translator, 'hgnc_id')
##########################################################################################################################
#%%

# reactionsto be added
add = ['ATPS4m','ENO','PDHm','PYK','G3PD1','G6PDH2r','AKGDm','CYOOm3','r0913','GLCt2_2',
       'EX_glc(e)','EX_fru(e)','EX_ppa(e)','EX_but(e)','EX_hdca(e)','EX_ocdca(e)',
       'EX_arach(e)','EX_doco13ac_','EX_lgnc(e)','EX_ala_L(e)','EX_arg_L(e)','EX_asn_L(e)',
       'EX_asp_L(e)','EX_cys_L(e)','EX_gln_L(e)','EX_glu_L(e)','EX_gly(e)','EX_his_L(e)',
       'EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)','EX_met_L(e)','EX_phe_L(e)','EX_pro_L(e)',
       'EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)','EX_hx(e)','EX_octa(e)','EX_ttdca(e)']

for r in add:
    confidence[r]= 3
    
label = str(serie.gsms[gsm].metadata['title']).split('_')

print('Building model '+ gsm + label[1]+label[2])
newmodel = Builder.build_model(recon22,confidence,[])
newmodel.name = 'Aging Skeletal Muscle '+ gsm + label[1]+label[2]
       
    ##############################
    # Model Curation :
    ############################
    
    # Reactions to be added
r0509 = Reaction ( id = 'r0509',
                  name = 'Succinate:ubiquinone oxidoreductase Citrate cycle (TCA cycle) EC:1.3.5.1',
                  subsystem = 'TCA cycle / electron transport chain', 
                  lower_bound = 0 , 
                  upper_bound = 1000)

newmodel.add_reaction(r0509)                  
r0509.build_reaction_from_string('succ_m + q10_m --> fum_m+ q10h2_m',verbose = True)

r0655 = Reaction ( id = 'r0655',
                  name = '3-Methylbutanoyl-CoA:(acceptor) 2,3-oxidoreductase; EC:1.3.99.10',
                  subsystem = 'Valine, leucine and isoleucine degradation', 
                  lower_bound = -1000, 
                  upper_bound = 1000)

newmodel.add_reaction(r0655)                  
r0655.build_reaction_from_string('ivcoa_m + q10_m <=> 3mb2coa_m + q10h2_m ',verbose = True)
    # Reactions to be removed    
to_remove = ['r1109','PYRt2r']
try:
    newmodel.remove_reactions(to_remove, remove_orphans=True)
except:
    print('not present')
    pass
    
    # reactions to be corrected
newmodel.reactions.CYOOm3.add_metabolites({'h_m':-7.9,'h_i':4 ,'h_c':0}, combine=False)

    # Fix reaction reversibility/directionality
newmodel.reactions.PDHm.bounds = 0, 1000
    
    ############################
    
    #prune unused metabolites
remov = cobra.manipulation.delete.prune_unused_metabolites(newmodel)    
while len(remov) != 0:
    remov = cobra.manipulation.delete.prune_unused_metabolites(newmodel)
    
    #write newmodel in SBML
cobra.io.write_sbml_model(newmodel, filename= output_folder+'ASK_'+gsm+label[1]+label[2]+'.xml')
print('===============================')

###############################################################################################################
    




