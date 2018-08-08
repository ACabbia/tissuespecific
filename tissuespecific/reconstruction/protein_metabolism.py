#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 10:28:10 2017

@author: acabbia
"""
     
from cobra import Reaction , Metabolite 
from Bio import ExPASy , SeqIO , SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis

INF = float('inf')

class Protein_metabolism:
              
     def add_prot_metabolism(model):
          prot_metab_rxns = []
          added_mets = []
          print('Model name:', model.name)
          print('1. Adding tRNAs ...')
          
          tRNA = Metabolite(id='tRNA_c',
                            formula = '',
                            name= 'tRNA',
                            compartment='c')
          
          tRNA_synthesis = Reaction(id= 'EX_tRNA', 
                                    name= 'tRNA synthesis(exchange)', 
                                    subsystem= 'tRNA metabolism',
                                    lower_bound=-INF,
                                    upper_bound= INF)
          
          tRNA_synthesis.add_metabolites({tRNA: 1})
          
          model.add_reaction(tRNA_synthesis)
          prot_metab_rxns.append(tRNA_synthesis.id)
          added_mets.append(tRNA)
     
          print('2. Adding Aminoacyl-tRNAs synthesis pathways...')
          
          ppi_c = model.metabolites.ppi_c
          pi_c = model.metabolites.pi_c
          atp_c = model.metabolites.atp_c
          adp_c = model.metabolites.adp_c
          amp_c = model.metabolites.adp_c
          
          aa_list = ['tyr_L_c','trp_L_c','thr_L_c','leu_L_c','ile_L_c','lys_L_c','ala_L_c','val_L_c','met_L_c','ser_L_c',
                     'asp_L_c','gly_c','pro_L_c','cys_L_c','glu_L_c','gln_L_c','arg_L_c','phe_L_c','his_L_c','asn_L_c']
          
          for aa in aa_list:
               
               Aminoacyl_tRNA = Metabolite(id= aa[0:3]+"_tRNA_c",
                                           formula='',
                                           name= "amminoacyl_"+aa+"_tRNA",
                                           compartment='c')
               
               aa_tRNA_synth = Reaction(id= aa[0:3]+'_tRNA_synth',
                                        name=aa[0:3]+' tRNA synthesis',
                                        subsystem='protein synthesis',
                                        lower_bound= -INF,
                                        upper_bound=  INF)
               
               aa_tRNA_synth.add_metabolites({               
                         model.metabolites.get_by_id(aa): -1,
                                                  atp_c : -1,
                                                   tRNA : -1,
                                         Aminoacyl_tRNA : 1,
                                                  ppi_c : 1,
                                                  amp_c : 1                                        
                                                  })
    
               model.add_reaction(aa_tRNA_synth)
               prot_metab_rxns.append(aa_tRNA_synth.id)
               added_mets.append(Aminoacyl_tRNA)
               
          print('3. Adding muscle proteins synthesis pathways...')
          
          uniprot_id = ['P68133', # actin
                        'P12882', # myosin heavy chain I and IIx
                        'Q9UKX2', # myosin heavy chain IIa
                        'Q9Y623', # myosin heavy chain IIb
                        'Q96A32', # myosin light chain type 2
                        'P05976', # myosin light chain type 1
                        'Q9H1R3', # myosin light chain kinase
                        'P09493', # tropomyosin type 1 (slow)
                        'P67936', # tropomyosin type 2 (fast)
                        'P63316', # troponin C type 1 (slow)
                        'P02585', # troponin C type 2 (fast)
                        'P19237', # troponin I type 1 (slow)
                        'P48788', # troponin I type 2 (fast)
                        'P13805', # troponin T type 1 (slow)
                        'P45378', # troponin T type 1 (fast)
                        ]
          
          translate_1to3 = {'A': 'ala',
                            'R': 'arg', 
                            'N': 'asn',
                            'D': 'asp',
                            'C': 'cys',
                            'Q': 'gln',
                            'E': 'glu',
                            'G': 'gly',
                            'H': 'his',
                            'I': 'ile',
                            'L': 'leu',
                            'K': 'lys',
                            'M': 'met',
                            'F': 'phe',
                            'P': 'pro',
                            'S': 'ser',
                            'T': 'thr',
                            'W': 'trp',
                            'Y': 'tyr',
                            'V': 'val' } 
          MW = []
          for prot_id in uniprot_id:
               with ExPASy.get_sprot_raw(prot_id) as handle:  #<--- why so slow???
                    Seq_record = SeqIO.read(handle, 'swiss')
                    handle.close()
                    print('Adding: '+ Seq_record.description.split('=')[1].split(';')[0] )
                    MW.append(SeqUtils.molecular_weight(Seq_record.seq, 'protein'))
                    
                    aa_count = ProteinAnalysis((str(Seq_record.seq)),True).count_amino_acids()
                    
                    protein = Metabolite(id=Seq_record.name[0:5]+'_protein_c',
                                         formula='',
                                         name = Seq_record.description.split()[1].split('=')[1],
                                         compartment='c')
                    
                    protein_synth = Reaction(id=Seq_record.name[0:5]+'_protein_synth',
                                             name=Seq_record.description.split()[1].split('=')[1]+' synthesis',
                                             subsystem='protein metabolism',
                                             lower_bound= -INF,
                                             upper_bound=  INF)
                    
                    model.add_metabolites([protein])
                    model.add_reaction(protein_synth)
                    
                    stoich = {translate_1to3[k]+'_tRNA_c':-v for k,v in aa_count.items()}
                    stoich2 = {atp_c.id: -(len(str(Seq_record.seq))),
                              pi_c.id:  (len(str(Seq_record.seq))),
                              adp_c.id:  (len(str(Seq_record.seq))),
                              tRNA.id:  (len(str(Seq_record.seq))),
                              protein.id:  1}
                    
                    stoich.update(stoich2)
                    
                    protein_synth.add_metabolites(stoich)
          
                    prot_metab_rxns.append(protein_synth.id)
                    added_mets.append(protein)
                              
          print('4. Adding muscle fibers synthesis pathways...')
          
          actin = model.metabolites.ACTS__protein_c
          myosinHC1 = model.metabolites.MYH1__protein_c
          myosinHC2a = model.metabolites.MYH2__protein_c
          myosinHC2b = model.metabolites.MYH4__protein_c
          myosinLC2 = model.metabolites.MLRS__protein_c
          myosinLC1 = model.metabolites.MYL1__protein_c
          myosinK = model.metabolites.MYLK2_protein_c
          tropomyosin1 = model.metabolites.TPM1__protein_c
          tropomyosin2 = model.metabolites.TPM4__protein_c
          troponinC1 = model.metabolites.TNNC1_protein_c
          troponinC2 = model.metabolites.TNNC2_protein_c
          troponinI1 = model.metabolites.TNNI1_protein_c
          troponinI2 = model.metabolites.TNNI2_protein_c
          troponinT1 = model.metabolites.TNNT1_protein_c
          troponinT2 = model.metabolites.TNNT3_protein_c
               
          # Myosin Type 1 muscle fiber : 2 myosinHC1 + 2 myosinK + 2 myosinLC1 -> 1 myo_fiber_1
               
          myo_fiber_1 = Metabolite(id='myo_fiber_1_c',
                                  formula= '',
                                  name='Myosin fiber Type I',
                                  compartment = 'c')
         
          myo_fiber_1_synth = Reaction(id='myo_fiber_1_synth',
                                      name='Myosin fiber Type 1 synthesis',
                                      subsystem = 'protein metabolism',
                                      lower_bound = -INF,
                                      upper_bound =  INF)
         
          myo_fiber_1_synth.add_metabolites({
                   myosinHC1: -2,
                   myosinK: -2,
                   myosinLC1:-2,
                   myo_fiber_1: 1})
         
          model.add_reaction(myo_fiber_1_synth)
          
          prot_metab_rxns.append(myo_fiber_1_synth.id)
          added_mets.append(myo_fiber_1)
         
          # Myosin Type 2a muscle fiber: 2 myosinHC2a + 2 myosinK + 2 myosinLC2 -> 1 myo_fiber_2a
          
          myo_fiber_2a = Metabolite(id='myo_fiber_2a_c',
                                    formula= '',
                                    name='Myosin fiber Type IIa',
                                    compartment = 'c')
         
          myo_fiber_2a_synth = Reaction(id='myo_fiber_2a_synth',
                                        name='Myosin fiber Type 2a synthesis',
                                        subsystem = 'protein metabolism',
                                        lower_bound = -INF,
                                        upper_bound =  INF)
         
          myo_fiber_2a_synth.add_metabolites({
                   myosinHC2a: -2,
                   myosinK: -2,
                   myosinLC2:-2,
                   myo_fiber_2a: 1})
         
          model.add_reaction(myo_fiber_2a_synth)
          
          prot_metab_rxns.append(myo_fiber_2a_synth.id)
          added_mets.append(myo_fiber_2a)
              
         # Myosin Type 2b muscle fiber:  2 myosinHC2b + 2 myosinK + 2 myosinLC2 -> 1 myo_fiber_2b
          
         
          
          myo_fiber_2b = Metabolite(id='myo_fiber_2b_c',
                                    formula= '',
                                    name='Myosin fiber Type IIb',
                                    compartment = 'c')
         
          myo_fiber_2b_synth = Reaction(id='myo_fiber_2b_synth',
                                        name='Myosin fiber Type 2b synthesis',
                                        subsystem = 'protein metabolism',
                                        lower_bound = -INF,
                                        upper_bound =  INF)
         
          myo_fiber_2b_synth.add_metabolites({
                    myosinHC2b: -2,
                    myosinK: -2,
                    myosinLC2:-2,
                    myo_fiber_2b: 1})
         
          model.add_reaction(myo_fiber_2b_synth)
          
          prot_metab_rxns.append(myo_fiber_2b_synth.id)
          added_mets.append(myo_fiber_2b)
         
         # Myosin Type 2x muscle fiber : 2 myosinHC1 + 2 myosinK + 2 myosinLC2 -> 1 myo_fiber_2x
          
         
          myo_fiber_2x = Metabolite(id='myo_fiber_2x_c',
                                    formula= '',
                                    name='Myosin fiber Type IIx',
                                    compartment = 'c')
         
          myo_fiber_2x_synth = Reaction(id='myo_fiber_2x_synth',
                                        name='Myosin fiber Type 2x synthesis',
                                        subsystem = 'protein metabolism',
                                        lower_bound = -INF,
                                        upper_bound =  INF)
         
          myo_fiber_2x_synth.add_metabolites({
                    myosinHC1: -2,
                    myosinK: -2,
                    myosinLC2:-2,
                    myo_fiber_2x: 1})
         
          model.add_reaction(myo_fiber_2x_synth)
          
          prot_metab_rxns.append(myo_fiber_2x_synth.id)
          added_mets.append(myo_fiber_2x)
             
          # Tropomyosin Type 1 complex: 2 tropomyosin1 ->  1 tropo_myo_fiber_1
             
          tropo_myo_fiber_1 = Metabolite(id="tropo_myo_fiber_1_c",
                                          formula='',
                                          name="Tropomyosin Type 1 complex",
                                          compartment= 'c')
           
          tropo_myo_fiber_1_synth = Reaction(id="tropo_myo_fiber_1_synth",
                                              name="Troponin type I complex synthesis",
                                              subsystem='protein metabolism',
                                              upper_bound=  INF,
                                              lower_bound= -INF)
           
          tropo_myo_fiber_1_synth.add_metabolites({
                     tropomyosin1: -2,
                     tropo_myo_fiber_1: 1})
                    
          model.add_reaction(tropo_myo_fiber_1_synth)    
         
          prot_metab_rxns.append(tropo_myo_fiber_1_synth.id)
          added_mets.append(tropo_myo_fiber_1)
          
         # Tropomyosin Type 2 complex: 2 tropomyosin2 -> 1 tr_myo_fiber_2
          
          tropo_myo_fiber_2 = Metabolite(id="tropo_myo_fiber_2_c",
                                          formula='',
                                          name="Tropomyosin Type 2 complex",
                                          compartment= 'c')
           
          tropo_myo_fiber_2_synth = Reaction(id="tropo_myo_fiber_2_synth",
                                              name="Troponin type II complex synthesis",
                                              subsystem='protein metabolism',
                                              upper_bound=  INF,
                                              lower_bound= -INF)
           
          tropo_myo_fiber_2_synth.add_metabolites({
                     tropomyosin2: -2,
                     tropo_myo_fiber_2: 1})
                    
          model.add_reaction(tropo_myo_fiber_2_synth)
          
          prot_metab_rxns.append(tropo_myo_fiber_2_synth.id)
          added_mets.append(tropo_myo_fiber_2)
          
         # Troponin complex Type 1: 1 troponinC1 + 1 troponinI1 + 1 troponinT1 -> 1 tropo_complex_1
          tropo_complex_1= Metabolite(id="tropo_complex_1_c",
                                          formula='',
                                          name="Troponin Type 1 complex",
                                          compartment= 'c')
           
          troponin_1_synth = Reaction(id="troponin_1_synth",
                                      name="Troponin I complex synthesis", 
                                      subsystem='protein metabolism',
                                      upper_bound=  INF,
                                      lower_bound= -INF)
           
          troponin_1_synth.add_metabolites({
                    troponinC1: -1,
                    troponinI1: -1,
                    troponinT1: -1,
                    tropo_complex_1:1})
                    
          model.add_reaction(troponin_1_synth)
          
          prot_metab_rxns.append(troponin_1_synth.id)
          added_mets.append(tropo_complex_1)
     
         # Troponin complex Type 2: 1 troponinC2 + 1 troponinI2 + 1 troponinT2 -> 1 tropo_complex_2
         
         
          tropo_complex_2= Metabolite(id="tropo_complex_2_c",
                                          formula='',
                                          name="Troponin Type 2 complex",
                                          compartment= 'c')
           
          troponin_2_synth = Reaction(id="troponin_2_synth",
                                      name="Troponin II complex synthesis", 
                                      subsystem='protein metabolism',
                                      upper_bound=  INF,
                                      lower_bound= -INF)
           
          troponin_2_synth.add_metabolites({
                    troponinC2: -1,
                    troponinI2: -1,
                    troponinT2: -1,
                    tropo_complex_2:1})
                    
          model.add_reaction(troponin_2_synth)
          
          prot_metab_rxns.append(troponin_2_synth.id)
          added_mets.append(tropo_complex_2)
          
         # Contractile unit type 1 : 7 actin + 7 myo_fiber_1 + 1 tropo_myo_fiber_1 + 1 tropo_complex_1 -> 1 contr_unit_1
          
          contr_unit_1 = Metabolite(id='contr_unit_1_c',
                                   formula = '',
                                   name='myofiber contractile unit type I',
                                   compartment='c')
         
          contr_unit_1_synth = Reaction(id='contr_unit_I_synth',
                                       name='Synthesis of myofiber contractile unit type I',
                                       subsystem= 'protein metabolism',
                                       upper_bound=  INF,
                                       lower_bound= -INF)
         
          contr_unit_1_synth.add_metabolites({
                    actin:-7,
                    myo_fiber_1: -7,
                    tropo_myo_fiber_1:-1,
                    tropo_complex_1:-1,
                    contr_unit_1:1})
                    
          model.add_reaction(contr_unit_1_synth)
          
          ###prot_metab_rxns.append(contr_unit_1_synth.id)
          added_mets.append(contr_unit_1)
          
         # Contractile unit type 2a: 7 actin + 7 myo_fiber_2a + 1 tropo_myo_fiber_2 + 1 tropo_complex_2 -> 1 contr_unit_2a
           
          contr_unit_2a = Metabolite(id='contr_unit_2a_c',
                                   formula = '',
                                   name='myofiber contractile unit type IIa',
                                   compartment='c')
         
          contr_unit_2a_synth = Reaction(id='contr_unit_IIa_synth',
                                       name='Synthesis of myofiber contractile unit type IIa',
                                       subsystem= 'protein metabolism',
                                       upper_bound=  INF,
                                       lower_bound= -INF)
         
          contr_unit_2a_synth.add_metabolites({
                    actin:-7,
                    myo_fiber_2a: -7,
                    tropo_myo_fiber_2:-1,
                    tropo_complex_2:-1,
                    contr_unit_2a:1})
                    
          model.add_reaction(contr_unit_2a_synth)
          
          ###prot_metab_rxns.append(contr_unit_2a_synth.id)
          added_mets.append(contr_unit_2a)
          
         # Contractile unit type 2b: 7 actin + 7 myo_fiber_2b + 1 tropo_myo_fiber_2 + 1 tropo_complex_2 -> 1 contr_unit_2b
         
          contr_unit_2b = Metabolite(id='contr_unit_2b_c',
                                   formula = '',
                                   name='myofiber contractile unit type IIb',
                                   compartment='c')
         
          contr_unit_2b_synth = Reaction(id='contr_unit_IIb_synth',
                                       name='Synthesis of myofiber contractile unit type IIb',
                                       subsystem= 'protein metabolism',
                                       upper_bound=  INF,
                                       lower_bound= -INF)
         
          contr_unit_2b_synth.add_metabolites({
                    actin:-7,
                    myo_fiber_2b: -7,
                    tropo_myo_fiber_2:-1,
                    tropo_complex_2:-1,
                    contr_unit_2b:1})
                    
          model.add_reaction(contr_unit_2b_synth)
          
          ###prot_metab_rxns.append(contr_unit_2b_synth.id)
          added_mets.append(contr_unit_2b)
          
         # Contractile unit type 2x: 7 actin + 7 myo_fiber_2x + 1 tropo_myo_fiber_2 + 1 tropo_complex_2 -> 1 contr_unit_2x
          
          contr_unit_2x = Metabolite(id='contr_unit_2x_c',
                                   formula = '',
                                   name='myofiber contractile unit type IIx',
                                   compartment='c')
         
          contr_unit_2x_synth = Reaction(id='contr_unit_IIx_synth',
                                       name='Synthesis of myofiber contractile unit type IIx',
                                       subsystem= 'protein metabolism',
                                       upper_bound=  INF,
                                       lower_bound= -INF)
         
          contr_unit_2x_synth.add_metabolites({
                    actin:-7,
                    myo_fiber_2x: -7,
                    tropo_myo_fiber_2:-1,
                    tropo_complex_2:-1,
                    contr_unit_2x:1})
                    
          model.add_reaction(contr_unit_2x_synth)
          
          
          ###prot_metab_rxns.append(contr_unit_2x_synth.id)
          added_mets.append(contr_unit_2x)
          ''' 
          print('5. Adding contractile units accumulation... (sink reactions)')
    
          # Contractile unit 1 sink: contr_unit_1 ->
          
          contr_unit_1_storage = Reaction(id='contr_unit_1_storage',
                                          name='contractile unit type 1 storage (sink)',
                                          subsystem='protein metabolism_',
                                          upper_bound= INF,
                                          lower_bound=-INF)
          
          contr_unit_1_storage.add_metabolites({
                    contr_unit_1: -1
                    })
                    
          model.add_reaction(contr_unit_1_storage)     
          
          prot_metab_rxns.append(contr_unit_1_storage.id)
          
           # Contractile unit 2a sink:  contr_unit_2a ->
           
          contr_unit_2a_storage = Reaction(id='contr_unit_2a_storage',
                                          name='contractile unit type 2a storage (sink)',
                                          subsystem='protein metabolism_',
                                          upper_bound= INF,
                                          lower_bound=-INF)
          
          contr_unit_2a_storage.add_metabolites({
                    contr_unit_2a: -1
                    })
                    
          model.add_reaction(contr_unit_2a_storage)
          
          prot_metab_rxns.append(contr_unit_2a_storage.id)
         
           # Contractile unit 2b sink: contr_unit_2b ->
          
          contr_unit_2b_storage = Reaction(id='contr_unit_2b_storage',
                                          name='contractile unit type 2b storage (sink)',
                                          subsystem='protein metabolism_',
                                          upper_bound= INF,
                                          lower_bound=-INF)
          
          contr_unit_2b_storage.add_metabolites({
                    contr_unit_2b: -1
                    })
                    
          model.add_reaction(contr_unit_2b_storage)
          
          prot_metab_rxns.append(contr_unit_2b_storage.id)
          
           # Contractile unit 2x sink: contr_unit_2x ->
          
          contr_unit_2x_storage = Reaction(id='contr_unit_2x_storage',
                                          name='contractile unit type 2x storage (sink)',
                                          subsystem='protein metabolism_',
                                          upper_bound= INF,
                                          lower_bound=-INF)
          
          contr_unit_2x_storage.add_metabolites({
                    contr_unit_2x: -1
                    })
                    
          model.add_reaction(contr_unit_2x_storage)
          
          prot_metab_rxns.append(contr_unit_2x_storage.id)
          '''
          
          print('5. Add contractile fibers accumulation (storage)*')
          
          myofibers_storage = Reaction(id='myofibers_storage',
                                          name='contractile units storage (sink)',
                                          subsystem='protein metabolism_',
                                          upper_bound= INF,
                                          lower_bound=-INF)
          
          myofibers_storage.add_metabolites({
                  contr_unit_1: -1,
                  contr_unit_2a: -1,
                  contr_unit_2b: -1,
                  contr_unit_2x: -1
                  })
          
          model.add_reaction(myofibers_storage)
          prot_metab_rxns.append(myofibers_storage.id)
          
          print('Adding ATP_expenditure reaction...')
         
          pi = model.metabolites.pi_c
          atp = model.metabolites.atp_c
          adp = model.metabolites.adp_c
                   
          ATP_expenditure = Reaction(id= 'ATP_expenditure',
                                     name= 'Total energy expenditure',
                                     subsystem = 'energy metabolism',
                                     upper_bound = INF,
                                     lower_bound = -INF)
                                     
          ATP_expenditure.add_metabolites({
                    atp: -1,
                    adp: 1,
                    pi: 1})

          model.add_reaction(ATP_expenditure)
          
          print('----------------------------------------')
          print(len(prot_metab_rxns)+4, "reactions have been added to the model")
          print(len(added_mets), "metabolites have been added to the model")
          print( 'All done!')
          print('----------------------------------------')
     
          return MW , prot_metab_rxns
          
     def select_subsystem(model, subsystem):
          selected = []
          for r in model.reactions:
               if r.subsystem == subsystem:
                    selected.append(r.id)
          return selected
          
     def convert_to_irreversible(model,reactions_to_convert): ##### TO DO : implement ATP cost in reverse reaction (1 ATP)
         """
         Split reversible reactions into two irreversible reactions.
         
         cobra_model: A Cobra model object which will be modified in place.
         reactions_to convert: list of reaction_id to convert 
         
         """
         reactions_to_add = []
         convertedlist = []
         
         atp_metabolites = [
         model.metabolites.ppi_c,
         model.metabolites.pi_c,
         model.metabolites.atp_c,
         model.metabolites.adp_c,
         model.metabolites.adp_c] 
         
         for r in reactions_to_convert:
              
              reaction = model.reactions.get_by_id(r)
              if reaction.lower_bound < 0 and reaction.upper_bound > 0:
                   reverse_reaction = Reaction(id=reaction.id + '_reverse',
                                               lower_bound = reaction.lower_bound,
                                               upper_bound = 0)
                   
                   reaction.lower_bound = 0
                   reaction.upper_bound = reaction.upper_bound
                   
                 
                   reverse_reaction_dict = {k: v * -1 for k, v in reaction.metabolites.items() 
                        if k not in atp_metabolites} # <- disallows recovery of ATP metabolites from protein degradation
                   
                   reverse_reaction.add_metabolites(reverse_reaction_dict)
                                      
                   reverse_reaction.subsystem = reaction.subsystem
                   reactions_to_add.append(reverse_reaction)
                   convertedlist.append(reverse_reaction.id) 
                   
         model.add_reactions(reactions_to_add)
         
         return convertedlist



     
     
     
     
