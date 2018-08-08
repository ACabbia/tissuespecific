#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 15:38:08 2017

@author: acabbia
"""

import cobra.test
import biomart
import pandas as pd
import numpy as np

model = cobra.test.create_test_model('salmonella')

class GSEA:
    
    def GSet_from_CBmodel(model):
        ### define metabolite-associated gene sets from SBML model  {ref: reporter metabolites  / Patil and Nielsen 2005}
        GS = {met.name: [[g.id for g in model.reactions.get_by_id(r.id).genes] for r in met.reactions] for met in model.metabolites}
        GSet = {met : [geneID for sublist in genes for geneID in sublist] for met, genes in GS.items() }
        return GSet
    
    
    def hgnc_to_affyProbes(GSet):  # TO DO: add possibility to skip building of translation dict if it already exists
        # download translation info
        server=biomart.BiomartServer('http://www.ensembl.org/biomart')  
        dataset=server.datasets['hsapiens_gene_ensembl']
        response=dataset.search({'attributes':['hgnc_id','affy_hg_u133_plus_2']})
        # build translation dict
        d= {}
        for line in response.iter_lines():
            line = line.decode('utf-8')
            (a,b)=(line.split("\t"))
            d[a]=b
        # translate hgnc ID into affyProbe ID
        GSaffy = pd.DataFrame([[met,d[gene.replace('HGNC:HGNC:','HGNC:')]] for met, gene_list in GSet.items() for gene in gene_list],columns=['metabolites', 'probes'])
        #drop rows with empty values
        GSaffy['probes'].replace('', np.nan, inplace = True)
        GSaffy.dropna(subset= ['probes'], inplace = True)
        GSaffy = GSaffy[GSaffy.columns[::-1]]
        return GSaffy

    def write_GSet_affy(GSaffy, filename):
        GSaffy.to_csv(filename, sep='\t', header=True, index=False)
        
        
