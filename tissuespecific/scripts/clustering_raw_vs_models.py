#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 15:21:50 2018

@author: acabbia
"""

#raw data vs models clustering comparison

from cobra.io import read_sbml_model
import GEOparse
from sklearn import metrics
from sklearn.cluster import k_means , KMeans
from sklearn import linear_model, ensemble
from sklearn.model_selection import LeaveOneOut
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, cut_tree
from tissuespecific.reconstruction.analysis import Summary
from os import listdir
from os.path import isfile, join

outfolder = '/home/acabbia/out/mergione/dist/'

library_folder = '/home/acabbia/Documents/Muscle_Model/models/library_GEOmerge/'
model_file = '/home/acabbia/Documents/Muscle_Model/models/recon2.2.xml'
ref_model = read_sbml_model(model_file)

GEO_accession = ["GSE25941" , "GSE28422"]

def make_exp_mat(geo_accession):
    
    expression_matrix = pd.DataFrame()
    for geo in GEO_accession:
        serie = GEOparse.get_GEO(geo)
        for gsm in serie.gsms:
            if gsm in model_list:
                index = serie.gsms[gsm].table['ID_REF']
                
                lbl = str(serie.gsms[gsm].metadata['title']).split('_')
                
                if geo == "GSE28422":
                    lbl = lbl[3].split(' ')
                    label ='ASK_'+gsm+'_'+lbl[0]
                    
                elif geo == "GSE25941":
                    label ='ASK_'+gsm+'_'+lbl[1]
                    
                expression_matrix[label] = serie.gsms[gsm].table['VALUE']
        
    expression_matrix.index = index
    return expression_matrix

model_list = [f.split('_')[1] for f in listdir(library_folder) if isfile(join(library_folder, f))]
expression_matrix = make_exp_mat(GEO_accession)

# 1) Hierarchical clustering

## Raw data clustering
# Hierarchical clustering (metric = cosine similarity)
Z = linkage(expression_matrix,'complete','cosine')
hc_class = pd.DataFrame(cut_tree(Z, n_clusters=2),index = labels, columns = ['hc_raw_class'])

## Models clustering
# Prepare reactions , metabolites and gene matrices
table, reactions_matrix, metabolite_matrix, gene_matrix = Summary.report_make_table(library_folder, ref_model)
# Hierarchical clustering (metric = jaccard)
Z_r = linkage(reactions_matrix.T, 'complete', 'jaccard')
Z_m = linkage(metabolite_matrix.T, 'complete', 'jaccard')
Z_g = linkage(gene_matrix.T, 'complete','jaccard')
modelR_class = pd.DataFrame(cut_tree(Z_r, n_clusters=2),index = labels, columns = ['hc_model_R_class'])
modelM_class = pd.DataFrame(cut_tree(Z_m, n_clusters=2),index = labels, columns = ['hc_model_M_class'])
modelG_class = pd.DataFrame(cut_tree(Z_g, n_clusters=2),index = labels, columns = ['hc_model_G_class'])

# Final class DF
class_df = pd.concat([ground_truth, hc_class, modelR_class, modelM_class, modelG_class], axis = 1)

# Evaluation 1: Adj. Rand Score
metrics.adjusted_rand_score(class_df['ground_truth'].values,class_df['hc_raw_class'].values)
metrics.adjusted_rand_score(class_df['ground_truth'].values,class_df['hc_model_R_class'].values)
metrics.adjusted_rand_score(class_df['ground_truth'].values,class_df['hc_model_M_class'].values)
metrics.adjusted_rand_score(class_df['ground_truth'].values,class_df['hc_model_G_class'].values)
# Evaluation 2: Adj. Mutual Info Score
metrics.adjusted_mutual_info_score(class_df['ground_truth'].values,class_df['hc_raw_class'].values)
metrics.adjusted_mutual_info_score(class_df['ground_truth'].values,class_df['hc_model_R_class'].values)
metrics.adjusted_mutual_info_score(class_df['ground_truth'].values,class_df['hc_model_M_class'].values)
metrics.adjusted_mutual_info_score(class_df['ground_truth'].values,class_df['hc_model_G_class'].values)
# Evaluation 3; homogeneity, completeness , V-score
metrics.homogeneity_completeness_v_measure(class_df['ground_truth'].values,class_df['hc_raw_class'].values)
metrics.homogeneity_completeness_v_measure(class_df['ground_truth'].values,class_df['hc_model_R_class'].values)
metrics.homogeneity_completeness_v_measure(class_df['ground_truth'].values,class_df['hc_model_M_class'].values)
metrics.homogeneity_completeness_v_measure(class_df['ground_truth'].values,class_df['hc_model_G_class'].values)
# Evaluation 4: Silhouette
metrics.silhouette_score(expression_matrix.values,class_df['hc_raw_class'].values, metric = 'cosine')
metrics.silhouette_score(reactions_matrix.T.values,class_df['hc_model_R_class'], metric = 'hamming')
metrics.silhouette_score(metabolite_matrix.T.values,class_df['hc_model_M_class'], metric = 'hamming')
metrics.silhouette_score(gene_matrix.T.values,class_df['hc_model_G_class'], metric = 'hamming')
# silhouette plots
pd.DataFrame(metrics.silhouette_samples(expression_matrix.values,class_df['hc_raw_class'].values, metric = 'cosine')).sort_values(by = 0)
pd.DataFrame(metrics.silhouette_samples(reactions_matrix.T.values,class_df['hc_model_R_class'], metric = 'hamming')).sort_values(by = 0)
pd.DataFrame(metrics.silhouette_samples(metabolite_matrix.T.values,class_df['hc_model_M_class'], metric = 'hamming')).sort_values(by = 0)
pd.DataFrame(metrics.silhouette_samples(gene_matrix.T.values,class_df['hc_model_G_class'], metric = 'hamming')).sort_values(by = 0)

# 2) K-means
centr, km_raw_class , i = k_means(expression_matrix.values, n_clusters=2, random_state = 0)
centr, km_ModR_class , i = k_means(reactions_matrix.T.values, n_clusters=2, random_state = 0)
centr, km_ModM_class , i = k_means(metabolite_matrix.T.values, n_clusters=2, random_state = 0)
centr, km_ModG_class , i = k_means(gene_matrix.T.values, n_clusters=2, random_state = 0)
# Evaluation 1: Adj. Rand Score
metrics.adjusted_rand_score(class_df['ground_truth'].values,km_raw_class)
metrics.adjusted_rand_score(class_df['ground_truth'].values,km_ModR_class)
metrics.adjusted_rand_score(class_df['ground_truth'].values,km_ModM_class)
metrics.adjusted_rand_score(class_df['ground_truth'].values,km_ModG_class)

metrics.adjusted_mutual_info_score(class_df['ground_truth'].values,km_raw_class)
metrics.adjusted_mutual_info_score(class_df['ground_truth'].values,km_ModR_class)
metrics.adjusted_mutual_info_score(class_df['ground_truth'].values,km_ModM_class)
metrics.adjusted_mutual_info_score(class_df['ground_truth'].values,km_ModG_class)

metrics.homogeneity_completeness_v_measure(class_df['ground_truth'].values,km_raw_class)
metrics.homogeneity_completeness_v_measure(class_df['ground_truth'].values,km_ModR_class)
metrics.homogeneity_completeness_v_measure(class_df['ground_truth'].values,km_ModM_class)
metrics.homogeneity_completeness_v_measure(class_df['ground_truth'].values,km_ModG_class)

# silhouette plots
metrics.silhouette_score(expression_matrix.values,km_raw_class, metric = 'cosine')
metrics.silhouette_score(reactions_matrix.T.values,km_ModR_class,metric = 'hamming')
metrics.silhouette_score(metabolite_matrix.T.values,km_ModM_class, metric = 'hamming')
metrics.silhouette_score(gene_matrix.T.values, km_ModG_class, metric = 'hamming')

pd.DataFrame(metrics.silhouette_samples(expression_matrix.values,km_raw_class, metric = 'cosine')).sort_values(by = 0)
pd.DataFrame(metrics.silhouette_samples(reactions_matrix.T.values,km_ModR_class,metric = 'hamming')).sort_values(by = 0)
pd.DataFrame(metrics.silhouette_samples(metabolite_matrix.T.values,km_ModM_class, metric = 'hamming')).sort_values(by = 0)
pd.DataFrame(metrics.silhouette_samples(gene_matrix.T.values, km_ModG_class, metric = 'hamming')).sort_values(by = 0)

#%%

# 3) supervised learning

#TO DO: make sure label matches correct ID (dict?)
def logistic_LOO(X, Y):
    print('LogisticRegression')
    loo = LeaveOneOut()
    TT_split= loo.split(Xm)
    
    counter = 0
    feature_importance = pd.DataFrame(index = X.columns)
    log_classifier = linear_model.LogisticRegression(penalty = 'l1', dual = False)
    
    for train_index, test_index in TT_split:
        log_classifier.fit(X.values[train_index], Y[train_index].ravel())
        feature_importance[X.index[test_index][0]] = (log_classifier.coef_.T)
        y_pred = log_classifier.predict(X.values[test_index])
        if (y_pred == Y[test_index]):
            counter += 1
        else:
            print('Misclassified:',X.index[test_index][0])
    print('prediction accuracy:', counter/len(X))
    return feature_importance

def RF_LOO(X,Y):
    print('RandomForest')
    loo = LeaveOneOut()
    TT_split = loo.split(X)
    
    counter = 0
    feature_importance = pd.DataFrame(index = X.columns)
    rf_classifier = ensemble.RandomForestClassifier(n_estimators = 1000, n_jobs=-1,random_state=0)
    
    for train_index, test_index in TT_split:
        rf_classifier.fit(X.values[train_index], Y[train_index].ravel())
        feature_importance[X.index[test_index][0]] = (rf_classifier.feature_importances_)
        y_pred = rf_classifier.predict(X.values[test_index])
        if (y_pred == Y[test_index]):
            counter += 1
        else:
            print('Misclassified:',X.index[test_index][0])
    print('prediction accuracy:', counter/len(X))
    return feature_importance

def ET_LOO(X,Y):
    print('ExtraTrees')
    loo = LeaveOneOut()
    TT_split = loo.split(X)
    
    counter = 0
    feature_importance = pd.DataFrame(index = X.columns)
    rf_classifier = ensemble.ExtraTreesClassifier(n_estimators = 1000, n_jobs=-1, random_state=0)
    
    for train_index, test_index in TT_split:
        rf_classifier.fit(X.values[train_index], Y[train_index].ravel())
        feature_importance[X.index[test_index][0]] = (rf_classifier.feature_importances_)
        y_pred = rf_classifier.predict(X.values[test_index])
        if (y_pred == Y[test_index]):
            counter += 1
        else:
            print('Misclassified:',X.index[test_index][0])
    print('prediction accuracy:', counter/len(X))
    return feature_importance
#%%    
def plot_clusters_silhouette(X,  n_clusters, metric, cluster_labels = None):
    """Plot the silhouette score for each cluster, given the distance matrix X.
    Parameters
    ----------
    X : array_like, shape [n_samples_a, n_samples_a]
        Distance matrix.
    cluster_labels : array_like
        List of integers which represents the cluster of the corresponding
        point in X. The size must be the same has a dimension of X.
    n_clusters : int
        The number of clusters.
    metric: metric to compute distance between data points
    """
    # if no labels are supplied (default), use K-means to predict them
    # fix random_state for reproducibility
    if cluster_labels is None:
        clusterer = KMeans(n_clusters=n_clusters, random_state=10)
        cluster_labels = clusterer.fit_predict(X)
            
    fig, ax = plt.subplots(1,1)
    fig.set_size_inches(20, 15)
    # Compute the silhouette scores for each sample
    sample_silhouette_values = metrics.silhouette_samples(X, cluster_labels, metric)
    silhouette_avg = metrics.silhouette_score(X, cluster_labels,metric)
    
    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = plt.get_cmap('inferno') (float(i) / n_clusters)
        ax.fill_betweenx(np.arange(y_lower, y_upper), 0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10 
    
    ax.set_xlabel("silhouette coefficient values")
    ax.set_ylabel("cluster label")

    # The vertical line for average silhoutte score of all the values
    ax.axvline(x = 0, color= 'black')
    ax.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax.set_yticks([])  # Clear the yaxis labels / ticks
    ax.set_xticks([-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])

    plt.suptitle("Silhouette analysis (n_clusters = {}), avg score {:.4f}".format(n_clusters, silhouette_avg),fontsize=16, fontweight='bold')    

#%%
# raw data:
print('Raw data')
Xm = expression_matrix
Y = np.array(ground_truth)
imp_raw_log =logistic_LOO(Xm, Y)
imp_raw_rf = RF_LOO(Xm, Y)
imp1 = ET_LOO(Xm, Y)
print('===========================================')
# model-based: 
# Rxns 
print('Reactions')
X_R = reactions_matrix.T
Y = np.array(ground_truth)
imp_rxn_log = logistic_LOO(X_R, Y)
imp_rxn_rf = RF_LOO(X_R, Y)
imp2 = ET_LOO(X_R, Y)
print('===========================================')
# Mets
print('Mets')
X_M = metabolite_matrix.T
Y = np.array(ground_truth)
imp_met_log = logistic_LOO(X_M, Y)
imp_met_rf = RF_LOO(X_M, Y)
print('===========================================')
# Genes 
print('Genes')
X_G = gene_matrix.T
Y = np.array(ground_truth)
imp_gen_log = logistic_LOO(X_G, Y)
imp_gen_rf = RF_LOO(X_G, Y)
print('===========================================')












































































