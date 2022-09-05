


import matplotlib.pyplot as plt
import numpy as np
from graphic_tools.violin_plot import plot_violin
from ML_tools.hierarchical_clustering import plot_hierarchical_clustering
from sklearn.neighbors import kneighbors_graph
from scipy import sparse
from sklearn.manifold import spectral_embedding
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, GridSearchCV, validation_curve, learning_curve
from sklearn.metrics import roc_auc_score
from scipy.sparse.linalg import eigsh, lobpcg
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
import pandas
from seaborn import heatmap
from graphic_tools.plot_decision_tree import plot_tree
from data_infos import protocols, experiments, table_columns, cell_types, dic_protocol_corresp
from pathlib import Path
from paths import figpath
import itertools
import sys
sys.setrecursionlimit(10000)
import matplotlib.cm as cm
from collections import OrderedDict
from itertools import chain
from collections import Counter
from load_data import ProteoDataset
from scipy.stats import pearsonr
from scipy.spatial.distance import pdist, squareform
import copy
from matplotlib.colors import ListedColormap, to_rgba
from matplotlib.colors import TwoSlopeNorm, LogNorm, ListedColormap
import matplotlib as mpl
from enrichments_plot import save_enrichments_results

from matplotlib.colors import LinearSegmentedColormap
from utils.heatmap_utils import create_mask_for_heatmap, make_rgb_transparent
import matplotlib.patches as mpatches

from pyensembl import EnsemblRelease
data = EnsemblRelease(77)
from utils.bioinf_conversion import convert_entrez_gene_id_list_to_gene_name
from graphic_tools.venn_diagrams import plot_venn_diagram, get_venn_sections


prodb = ProteoDataset(reset=False)

cell_type = "Fibro"


log = False
if log:
    text_log="log_"
else:
    text_log=""



list_protocols = protocols[cell_type]

list_replicates = experiments[cell_type]






""""""""""""""""""""""""""""""""""""
"""""""COUNT DIFF PROTEINS EV ProDA"""""""
""""""""""""""""""""""""""""""""""""

table = prodb.load_table(cell_type)

counts = np.zeros((len(list_protocols), len(list_protocols)))
for i, protocol1 in enumerate(list_protocols):
    for j, protocol2 in enumerate(list_protocols):

        if protocol1 == protocol2:
            distance = 0
            continue         
        col = "test_diff_ProDA_"+protocol1 + "-" + protocol2
        
        if col not in table:
            col = "test_diff_ProDA_"+protocol2 + "-" + protocol1

        mask = ~(table[col].isna())
        
        nb_diff = np.sum(np.abs(table[mask][col]))
        
        counts[i,j] = nb_diff

counts = pandas.DataFrame(counts, columns=list_protocols, index=list_protocols)    
print("count differential proteins ProDA")        

print(counts.astype(int).style.to_latex(hrules=True))
               
counts.to_csv(Path(figpath,"similarity_analysis", cell_type+"_counts_proteins_diff_ProDA_EV.csv"))

            
index, labels, silhouette = plot_hierarchical_clustering(counts, \
                        metric="precomputed", linkage_method="average", threshold=0.9, \
                        distance_matrix=False, criterion="threshold", labelsize=20, path=Path(figpath, "similarity_analysis", "HC", cell_type+"_HC_count_differential_proteins_ProDA.pdf"))





""""""""""""""""""""""""""""""""""""
"""""""COUNT DIFF PROTEINS EV DEP"""""""
""""""""""""""""""""""""""""""""""""

table = prodb.load_table(cell_type)

counts = np.zeros((len(list_protocols), len(list_protocols)))
for i, protocol1 in enumerate(list_protocols):
    for j, protocol2 in enumerate(list_protocols):

        if protocol1 == protocol2:
            distance = 0
            continue         
        col = "test_diff_DEP_"+protocol1 + "-" + protocol2
        
        if col not in table:
            col = "test_diff_DEP_"+protocol2 + "-" + protocol1
            
        mask = ~(table[col].isna())

        nb_diff = np.sum(np.abs(table[mask][col]))
        
        counts[i,j] = nb_diff

counts = pandas.DataFrame(counts, columns=list_protocols, index=list_protocols)    
print("count differential proteins DEP")        

print(counts.astype(int).style.to_latex(hrules=True))
               
counts.to_csv(Path(figpath,"similarity_analysis", cell_type+"_counts_proteins_diff_DEP_EV.csv"))

            
index, labels, silhouette = plot_hierarchical_clustering(counts, \
                        metric="precomputed", linkage_method="average", threshold=0.9, \
                        distance_matrix=False, criterion="threshold", labelsize=20, path=Path(figpath, "similarity_analysis", "HC", cell_type+"_HC_count_differential_proteins_DEP.pdf"))








""""""""""""""""""""""""""""""""""""
"""""""COUNT DISTANCE VECTORS EV EXPECTATION"""""""
""""""""""""""""""""""""""""""""""""
table = prodb.load_table(cell_type)

distances = np.zeros((len(list_protocols), len(list_protocols)))
for i, protocol1 in enumerate(list_protocols):
    for j, protocol2 in enumerate(list_protocols):

        
        if protocol1 == protocol2:
            distance = 0
            continue        
        col_protocol1 = "coeffs_"+protocol1
        
        col_protocol2 = "coeffs_"+protocol2
        
        if log:

            distance = np.linalg.norm(table[col_protocol1] - table[col_protocol2])
            
            # distance = np.sqrt(np.sum((table[col_protocol1] - table[col_protocol2])**2))
            
        else:
            distance = np.linalg.norm(2**table[col_protocol1] - 2**table[col_protocol2])

        
        distances[i,j] = distance
        
distances = pandas.DataFrame(distances, columns=list_protocols, index=list_protocols)    
print("distances EV expectation")        
print(distances.round(1).to_latex())
# distances.to_csv(Path(figpath,"similarity_analysis", cell_type_title+"_distances_EV_expectations.csv"))

            
index, labels, silhouette = plot_hierarchical_clustering(distances, \
                        metric="precomputed", linkage_method="average", threshold=0.9, \
                        distance_matrix=False, criterion="threshold", labelsize=20, path=Path(figpath, "similarity_analysis", "HC", text_log+cell_type+"_HC_EV_expectations.pdf"))

    
    

""""""""""""""""""""""""""""""""""""
"""""""COUNT DISTANCE VECTORS EV AVERAGE IMPUTATION"""""""
""""""""""""""""""""""""""""""""""""
table = prodb.load_table(cell_type)


table = prodb.load_table(cell_type)

average_cols = [col for col in table.columns if "average" in col and "imputation" in col]

mask = ~(table[average_cols].isna()).product(axis=1).astype(bool)



distances = np.zeros((len(list_protocols), len(list_protocols)))
for i, protocol1 in enumerate(list_protocols):
    for j, protocol2 in enumerate(list_protocols):

        
        if protocol1 == protocol2:
            distance = 0
            continue        
        col_protocol1 = "average_imputation_DEP_"+protocol1
        
        col_protocol2 = "average_imputation_DEP_"+protocol2
        
        if log:

            distance = np.linalg.norm(table[mask][col_protocol1] - table[mask][col_protocol2])
            
            # distance = np.sqrt(np.sum((table[col_protocol1] - table[col_protocol2])**2))
            
        else:
            distance = np.linalg.norm(2**table[mask][col_protocol1] - 2**table[mask][col_protocol2])

        
        distances[i,j] = distance
        
distances = pandas.DataFrame(distances, columns=list_protocols, index=list_protocols)    
print("distances EV average imputation")        
print(distances.round(1).to_latex())
# distances.to_csv(Path(figpath,"similarity_analysis", cell_type_title+"_distances_EV_expectations.csv"))

            
index, labels, silhouette = plot_hierarchical_clustering(distances, \
                        metric="precomputed", linkage_method="average", threshold=0.9, \
                        distance_matrix=False, criterion="threshold", labelsize=20, path=Path(figpath, "similarity_analysis", "HC", text_log+cell_type+"_HC_average_imputation.pdf"))





""""""""""""""""""""""""""""""""""""
"""""""COUNT DISTANCE APPARTITION DISPARITION EV"""""""
""""""""""""""""""""""""""""""""""""
 
table = prodb.load_table(cell_type)
 
distances = np.zeros((len(list_protocols), len(list_protocols)))
for i, protocol1 in enumerate(list_protocols):
    for j, protocol2 in enumerate(list_protocols):

        
        if protocol1 == protocol2:
            distance = 0
            continue        

        cols_protocol1 = table_columns[cell_type][protocol1]
        
        cols_protocol2 = table_columns[cell_type][protocol2]
        
        is_valid_protocol1 = (table["nb_valid_"+protocol1] >=1)
        is_valid_protocol2 = (table["nb_valid_"+protocol2] >=1)

        distance =  np.sum(is_valid_protocol1 != is_valid_protocol2)
        
        distances[i,j] = distance
        
distances = pandas.DataFrame(distances, columns=list_protocols, index=list_protocols)
print("distances EV hamming apparition disparition replicates")        
print(distances.astype(int).style.to_latex(hrules=True))

distances.to_csv(Path(figpath,"similarity_analysis", cell_type+"_distances_EV_hamming_distance_apparition_disparition.csv"))

index, labels, silhouette = plot_hierarchical_clustering(distances, \
                        metric="precomputed", linkage_method="average", threshold=0.9, \
                        distance_matrix=False, criterion="threshold", labelsize=20, path=Path(figpath, "similarity_analysis", "HC", cell_type+"_HC_EV_hamming_distance_apparition_disparition.pdf"))




""""""""""""""""""""""""""""""""""""
"""""""COUNT DISTANCE VECTORS EV REPLICATES"""""""
""""""""""""""""""""""""""""""""""""

table = prodb.load_table(cell_type)

distances = np.zeros((len(list_replicates), len(list_replicates)))
for i, col1 in enumerate(list_replicates):
    for j, col2 in enumerate(list_replicates):

        
        if col1 == col2:
            distance = 0
            continue        
        
        if log:
            distance =  np.linalg.norm(table[col1].fillna(0) - table[col2].fillna(0), axis=0)
            
        else:
            distance = np.linalg.norm(2**table[col1].fillna(0) - 2**table[col2].fillna(0), axis=0)
        
        distances[i,j] = distance
        
distances = pandas.DataFrame(distances, columns=list_replicates, index=list_replicates)
distances.columns = [col.replace("LFQI ","") for col in distances.columns]  
distances.index = [col.replace("LFQI ","") for col in distances.index]  

print("replicates_missing_values_at_0")
print(distances.round(1).to_latex())
# distances.to_csv(Path(figpath,"similarity_analysis", cell_type_title+"_distances_EV_replicates.csv"))
  
            
index, labels, silhouette = plot_hierarchical_clustering(distances, \
                        metric="precomputed", linkage_method="average", threshold=0.9, \
                        distance_matrix=False, criterion="threshold", labelsize=20, path=Path(figpath, "similarity_analysis", "HC", text_log+cell_type+"_HC_EV_replicates_missing_values_at_0.pdf"))


    
    
    
""""""""""""""""""""""""""""""""""""
"""""""COUNT DISTANCE VECTORS EV REPLICATES IMPUTATION"""""""
""""""""""""""""""""""""""""""""""""

table = prodb.load_table(cell_type)

imputed_cols = [col for col in table.columns if "LFQI" in col and "imputation" in col]

mask = ~(table[imputed_cols].isna()).product(axis=1).astype(bool)

print(len(table), len(table[mask]), "after imputation")



distances = np.zeros((len(imputed_cols), len(imputed_cols)))
for i, col1 in enumerate(imputed_cols):
    for j, col2 in enumerate(imputed_cols):

        
        if col1 == col2:
            distance = 0
            continue        
        
        if log:
            distance =  np.linalg.norm(table[mask][col1] - table[mask][col2],  axis=0)
            
        else:
            distance = np.linalg.norm(2**table[mask][col1] - 2**table[mask][col2], axis=0)
        
        distances[i,j] = distance
        
distances = pandas.DataFrame(distances, columns=list_replicates, index=list_replicates)
distances.columns = [col.replace("LFQI ","") for col in distances.columns]  
distances.index = [col.replace("LFQI ","") for col in distances.index]  

print("replicates imputation")
print(distances.round(1).to_latex())
  
            
index, labels, silhouette = plot_hierarchical_clustering(distances, \
                        metric="precomputed", linkage_method="average", threshold=0.9, \
                        distance_matrix=False, criterion="threshold", labelsize=20, path=Path(figpath, "similarity_analysis", "HC", text_log+cell_type+"_HC_EV_replicates_imputation.pdf"))


    
    
    

    





""""""""""""""""""""""""""""""""""""
"""""""COUNT DISTANCE VECTORS EV REPLICATES DISCARDING"""""""
""""""""""""""""""""""""""""""""""""

table = prodb.load_table(cell_type)

filter_ = table[list_replicates].isna().sum(axis=1)

reduced_table = table[(filter_==0)]

distances = np.zeros((len(list_replicates), len(list_replicates)))
for i, col1 in enumerate(list_replicates):
    for j, col2 in enumerate(list_replicates):
        
        if col1 == col2:
            distance = 0
            continue        
        
        if log:
            distance =  np.linalg.norm(reduced_table[col1].fillna(0) - reduced_table[col2].fillna(0), axis=0)
            
        else:
            distance = np.linalg.norm(2**reduced_table[col1].fillna(0) - 2**reduced_table[col2].fillna(0), axis=0)
        
        distances[i,j] = distance
        
distances = pandas.DataFrame(distances, columns=list_replicates, index=list_replicates)
distances.columns = [col.replace("LFQI ","") for col in distances.columns]  
distances.index = [col.replace("LFQI ","") for col in distances.index]  

print(distances.round(1).to_latex())
# distances.to_csv(Path(figpath,"similarity_analysis", cell_type_title+"_distances_EV_replicates.csv"))
  
            
index, labels, silhouette = plot_hierarchical_clustering(distances, \
                        metric="precomputed", linkage_method="average", threshold=0.9, \
                        distance_matrix=False, criterion="threshold", labelsize=20, path=Path(figpath, "similarity_analysis", "HC", text_log+cell_type+"_HC_EV_replicates_discarding.pdf"))

print("nb before discarding", len(table), "nb after discarding", len(reduced_table))







# # distances.to_csv(Path(figpath,"similarity_analysis", cell_type_title+"_distances_EV_expectations.csv"))

            
# index, labels, silhouette = plot_hierarchical_clustering(distances, \
#                         metric="precomputed", linkage_method="complete", threshold=0.9, \
#                         distance_matrix=False, criterion="threshold", labelsize=15, path=Path(figpath, "similarity_analysis", "HC", text_log+"test_diff_comparison_"+cell_type_title+".pdf"))


   
    
    
    


# """"""""""""""""""""""""""""""""""""
# """""""COUNT DISTANCE VECTORS LYSAT EXPECTATION"""""""
# """"""""""""""""""""""""""""""""""""


# conds_dic = {"Lysat 300 RPM 3h": "Lysat MS", 
#          "Lysat 300 RPM 6h": "Lysat MS.6h",
#          "Lysat 500 RPM 3h": "Lysat HS",
#          "Lysat 2D": "Lysat 2D"}


# list_conds = sorted(list(conds_dic.keys()))
# distances = np.zeros((len(list_conds), len(list_conds)))
# for i, cond1 in enumerate(list_conds):
#     for j, cond2 in enumerate(list_conds):
        
#         if cond1 == cond2:
#             distance = 0
#         continue
        
#         col_cond1 = "coeffs_"+conds_dic[cond1]  + "_" + cell_type
        
#         col_cond2 = "coeffs_"+conds_dic[cond2]  + "_" + cell_type

#         distance = np.linalg.norm(table[col_cond1] - table[col_cond2])
        
#         distances[i,j] = distance
        
# distances = pandas.DataFrame(distances, columns=list_conds, index=list_conds)    
# print("distances Lysat expectation")        
# print(distances)

