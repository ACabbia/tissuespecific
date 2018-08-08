    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 11:32:19 2017

@author: acabbia
"""
import GEOparse
import cobra
import numpy as np
import os
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import six
from scipy.spatial.distance import pdist , jaccard , squareform

INF = float('inf')


#%% 
class Summary:
    
    def model_summary(model):
          '''Statistics on the number of species, reactions and genes in a model'''
          print("This model contains" , len(model.reactions) , "reactions", len(model.metabolites), "metabolites", len(model.genes), "genes")
          return len(model.reactions) , len(model.metabolites), len(model.genes)

    def report_make_table(library_folder,ref_model):
    #returns a table containing statistics about a collection of models:
    #number of total and unique genes, reactions, metabolites
    #matrix info about wheter reaction/metabolite/gene [i] from parent model has beeen added
    
        table = pd.DataFrame(index = list(os.listdir(library_folder)))
        reactions_matrix = pd.DataFrame(index=[r.id for r in ref_model.reactions])
        metabolite_matrix = pd.DataFrame(index=[m.id for m in ref_model.metabolites])
        gene_matrix = pd.DataFrame(index=[g.id for g in ref_model.genes])
        
        n_rxns =[]
        n_mets =[]
        n_gene =[]
        
        for filename in os.listdir(library_folder):
            model = cobra.io.read_sbml_model(library_folder+filename)
            rxns = []
            mets = []
            genes = []
            
            print(filename)
            n_rx, n_mt , n_gn = Summary.model_summary(model)
            n_rxns.append(n_rx)
            n_mets.append(n_mt)
            n_gene.append(n_gn)
            
            label = str(filename).split('.')[0]
            
            for r in ref_model.reactions:
                if r in model.reactions:
                    rxns.append(1)
                else:
                    rxns.append(0)
                    
            for m in ref_model.metabolites:
                if m in model.metabolites:
                    mets.append(1)
                else:
                    mets.append(0)
                    
            for g in ref_model.genes:
                if g in model.genes:
                    genes.append(1)
                else:
                    genes.append(0)
             
            reactions_matrix[label] = pd.Series(rxns).values
            metabolite_matrix[label] = pd.Series(mets).values
            gene_matrix[label] = pd.Series(genes).values
        
        table['nr. of Reactions'] = n_rxns
        table['nr. of Metabolites'] = n_mets
        table['nr. of Genes'] = n_gene
        return table , reactions_matrix, metabolite_matrix, gene_matrix
        

    def report_clustering_plots(table , reactions_matrix, metabolite_matrix, gene_matrix, outfolder):
        ## Generate matrices of pairwise differences in number of reactions metabolites genes between models
        ## Plots hierarchical clustering results (heatmap)
        diff_array = np.zeros((len(table.index),len(table.index), len(table.columns)))
        for c in range(len(table.columns)):
            serie = table[table.columns[c]]
            for i in range(len(table.index)):
                for j in range(len(table.index)):
                    diff_array[i][j][c] = abs(serie[i]-serie[j])
            
        R_pw_diff = pd.DataFrame(data = diff_array[:,:,0],index=table.index, columns= table.index)
        M_pw_diff = pd.DataFrame(data = diff_array[:,:,1],index=table.index, columns= table.index)
        G_pw_diff = pd.DataFrame(data = diff_array[:,:,2],index=table.index, columns= table.index)
        
        reactions_matrix = reactions_matrix.T
        metabolite_matrix = metabolite_matrix.T
        gene_matrix = gene_matrix.T
                
        R_unique = reactions_matrix.T[reactions_matrix.sum() == 1]
        M_unique = metabolite_matrix.T[metabolite_matrix.sum() == 1]
        G_unique = gene_matrix.T[gene_matrix.sum() == 1]
                
        table['Unique Reactions'] = R_unique.sum().values
        table['Unique Metabolites'] = M_unique.sum().values
        table['Unique Genes'] = G_unique.sum().values
        
        #model clustering and plots
        jacc_R = squareform(pdist(reactions_matrix, metric = jaccard))
        jacc_M = squareform(pdist(metabolite_matrix, metric = jaccard))
        jacc_G = squareform(pdist(gene_matrix, metric=jaccard))
        
        distance_R = pd.DataFrame(jacc_R, index=reactions_matrix.index, columns= reactions_matrix.index)
        distance_M = pd.DataFrame(jacc_M, index=metabolite_matrix.index, columns= reactions_matrix.index)
        distance_G = pd.DataFrame(jacc_G, index=gene_matrix.index, columns= reactions_matrix.index)
        
        rg = sns.clustermap(distance_R)
        rg.fig.suptitle('Cluster by reactions', fontsize = 34)
        rg.savefig(outfolder+'/rxn_cluster.png' )
        
        mg = sns.clustermap(distance_M)
        mg.fig.suptitle('Cluster by metabolites', fontsize = 34)
        mg.savefig(outfolder+'/met_cluster.png' )
        
        gg = sns.clustermap(distance_G)
        gg.fig.suptitle('Cluster by genes', fontsize = 34)
        gg.savefig(outfolder+'/gene_cluster.png' )
        
        return table, R_pw_diff , M_pw_diff, G_pw_diff 

    def report(library_folder,ref_model):
        table, reactions_matrix, metabolite_matrix, gene_matrix = Summary.report_make_table(library_folder,ref_model)
        table,R_pw_diff , M_pw_diff, G_pw_diff  = Summary.report_clustering_plots_plots(table, reactions_matrix, metabolite_matrix, gene_matrix)
        return  table, R_pw_diff , M_pw_diff, G_pw_diff
    
    def series_similarity(S1,S2):
        try:
            len(S1)==len(S2)
        except:
            ValueError('The vectors must have the same length')
            raise
            
        v1 = S1.values
        v2 = S2.values
            
        sim = []
        for i in range(len(v1)):
            if v1[i]==1:
                if v2[i]==1:
                    sim.append(1)
        return abs((S1.sum())-len(sim))    

    def PW_similarity(PA_df):
        
        PW_array = np.zeros((len(PA_df.columns),len(PA_df.columns)))
        
        for i in range(len(PA_df.columns)):
            for j in range(len(PA_df.columns)):
                PW_array[i][j] = Summary.series_similarity(PA_df.iloc[:,i], PA_df.iloc[:,j])
        
        PW_df = pd.DataFrame(data = PW_array, columns = PA_df.columns, index = PA_df.columns)
        return PW_df

    def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                         header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                         bbox=[0, 0, 1, 1], header_columns=0,
                         ax=None, **kwargs):
        if ax is None:
            size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
            fig, ax = plt.subplots(figsize=size)
            ax.axis('off')
            
        mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns,rowLabels = data.index, **kwargs)
    
        mpl_table.auto_set_font_size(False)
        mpl_table.set_fontsize(font_size)
    
        for k, cell in  six.iteritems(mpl_table._cells):
            cell.set_edgecolor(edge_color)
            if k[0] == 0 or k[1] < header_columns:
                cell.set_text_props(weight='bold', color='w')
                cell.set_facecolor(header_color)
            else:
                cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
        return ax.get_figure()
    
    def table_as_png(dataframe, outname):
        tab = Summary.render_mpl_table(dataframe)
        tab.savefig(outname, dpi = 300 , bbox = 'tight', pad_inches = 0.1)


    def merge_PW_df(PW_sim_df, PW_diff_df):
        up = np.triu(PW_diff_df.values)
        dn = np.tril(PW_sim_df.values)
        new = up+dn
        new_df = pd.DataFrame(new, columns = PW_diff_df.columns, index = PW_sim_df.index)
        return new_df


    def heat(new_df, outname):
        fig, ax = plt.subplots(figsize=(30,30))
        s = sns.heatmap(new_df, annot=True, linewidths=.5, ax=ax, square = True)
        fig = s.get_figure()    
        fig.savefig(outname)
        
    def find_unique(matrix,model_id):
        '''
        finds unique reactions/metabolites/genes ID for a model
        '''
        um = matrix[matrix.T.sum() == 1]
        _id = um[model_id]
        unique = _id[_id==1]
        unique_list = list(unique.index)
        
        return unique_list    

class GDS_summary:
    
    accession = "GSE25941"
    def series_summary(accession):
        serie=GEOparse.get_GEO(accession)
        PA = pd.DataFrame()
        for gsm in serie.gsms:
             eset = serie.gsms[gsm].table
             call = eset['ABS_CALL']
             call.replace('P', 1 , True)
             call.replace('A', 0,  True)
             call.replace('M', 0,  True)
             PA[gsm] = call.values
        palette = ['#FBA90A','#1430B8']     
        fig, ax = plt.subplots(figsize=(20,20)) 
        sns.heatmap(PA, vmin=0, vmax=1, center = 0.5, cmap = palette)
             
            










                

     
                  
             