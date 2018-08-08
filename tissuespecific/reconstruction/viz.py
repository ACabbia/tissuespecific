#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 11:54:09 2017

@author: acabbia
"""

import networkx as nx
import pandas as pd


class visualization:
    
    # graph view
    def save_flow_graph(model, solution, filename):
        
        nnz_sol = solution.fluxes.iloc[solution.fluxes.nonzero()]
        RXN_Nodes = [r for r , v in nnz_sol.iteritems()]
        met_edges_in = { x : [m.id for m in model.reactions.get_by_id(x).reactants] for x in RXN_Nodes}
        met_edges_out = { x : [m.id for m in model.reactions.get_by_id(x).products] for x in RXN_Nodes}
        
        nodef=open(filename, 'w')
        for rID , v in met_edges_in.items():
            for met in v: 
                nodef.write(rID+'\t'+met)
                nodef.write("\n")
        for rID , v in met_edges_out.items():
            for met in v: 
                nodef.write(rID+'\t'+met)
                nodef.write("\n")
        nodef.close()
    
    def read_flow_graph(filename):
        G = nx.readwrite.edgelist.read_edgelist(filename)
        nx.drawing.draw_networkx(G)
        return G
        
    # ReconMap View
    def viz(solution, filename):
        nnz_sol = solution.fluxes.iloc[solution.fluxes.nonzero()]
        df = pd.DataFrame(columns = ['reactionIdentifier','lineWidth','color'])
        df['reactionIdentifier'] = nnz_sol.index
        df['lineWidth'] = (nnz_sol - (min(nnz_sol)) / (max(nnz_sol) - min(nnz_sol))) #<- Normalized fluxes 
        df['color'] = '#FF0000'# red 
        df.to_csv(filename, sep= '\t', header = True, index = False)