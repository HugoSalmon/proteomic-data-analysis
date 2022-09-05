





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
from load_data import ProteoDataset, list_comparisons
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

pandas.set_option('display.expand_frame_repr', False)

prodb = ProteoDataset()



list_cell_types = ["THP1", "MSCH", "Fibro"]

load_lysat=False


log = False
if log:
    text_log="log_"
else:
    text_log=""


table = prodb.load_human_table()


cell_type_title = "_".join(sorted(list_cell_types))




""""""""""""""""""""""""""""""""""""
"""""""COUNT DIFF PROTEINS EV ALL CELL TYPES"""""""
""""""""""""""""""""""""""""""""""""


list_protocols = sum([[protocol+"_"+cell_type for protocol in protocols[cell_type]] for cell_type in cell_types],[])
        
counts = np.zeros((len(list_protocols), len(list_protocols)))
for i, protocol1 in enumerate(list_protocols):
    for j, protocol2 in enumerate(list_protocols):

        if protocol1 == protocol2:
            distance = 0
            continue  
        
        cell_type1 = protocol1.split("_")[1]
            
        cell_type2 = protocol2.split("_")[1]

        if cell_type1!=cell_type2:
            nb_diff = np.nan
            
        else:
            
            table_cell_type = prodb.load_table(cell_type1)
          
            col = "test_diff_"+protocol1.split("_")[0] + "-" + protocol2.split("_")[0]
            
            
            
            if col not in table_cell_type:
                col = "test_diff_"+protocol2.split("_")[0] + "-" + protocol1.split("_")[0]
            
            nb_diff = np.sum(np.abs(table_cell_type[col]))
        
        counts[i,j] = nb_diff

counts = pandas.DataFrame(counts, columns=list_protocols, index=list_protocols)    
print("all diff counts")        

counts.to_csv(Path(figpath,"similarity_analysis", cell_type_title+"_comparison_cell_types_counts_proteins_diff_EV.csv"))

counts.index = [ind.replace("_"," ") for ind in counts.index]
counts.columns = [col.replace("_"," ") for col in counts.columns]
counts = counts.applymap(lambda s: "" if np.isnan(s) else str(int(s)))
print(counts.to_latex())
               



  





""""""""""""""""""""""""""""""""""""
"""""""COUNT DISTANCE FOLD CHANGE"""""""
""""""""""""""""""""""""""""""""""""

list_FC_cols = [col for col in table.columns if "FC" in col]
    
distances = np.zeros((len(list_FC_cols), len(list_FC_cols)))

for i, fold_change1 in enumerate(list_FC_cols):
    for j, fold_change2 in enumerate(list_FC_cols):

            
        if fold_change1 == fold_change2:
            distance = 0
            continue   

        if log:
            prepared_fold_change1 = np.log2(table[fold_change1].apply(lambda s: s if not np.isnan(s) else 1))
            prepared_fold_change2 = np.log2(table[fold_change2].apply(lambda s: s if not np.isnan(s) else 1))

        else:
            prepared_fold_change1 = table[fold_change1].apply(lambda s: s if not np.isnan(s) else 1)
            prepared_fold_change2 = table[fold_change2].apply(lambda s: s if not np.isnan(s) else 1)
        
        distance = np.linalg.norm(prepared_fold_change1 - prepared_fold_change2)
            
    
        
        distances[i,j] = distance
        
distances = pandas.DataFrame(distances, columns=list_FC_cols, index=list_FC_cols)    
print("distances Fold changes")  
distances.columns = [col.replace("FC_","").replace("_"," ") for col in distances.columns]    
distances.index = [col.replace("FC_","").replace("_"," ") for col in distances.index]      
  
print(np.round(distances, 1).to_latex())


import seaborn as sns
fig, ax = plt.subplots(1, figsize=(15,15))
sns.heatmap(distances, ax=ax)
fig.tight_layout()
# distances.to_csv(Path(figpath,"similarity_analysis", cell_type_title+"_distances_EV_expectations.csv"))

            
index, labels, silhouette = plot_hierarchical_clustering(distances, \
                        metric="precomputed", linkage_method="complete", threshold=0.9, \
                        distance_matrix=False, criterion="threshold", labelsize=15, path=Path(figpath, "similarity_analysis", "HC", text_log+"Fold_change_comparison_"+cell_type_title+".pdf"))


for comparison in ["HS/3D", "HS/2D"]:
    
    indexes = [index for index in distances.index if comparison in index]
    columns = [col for col in distances.columns if comparison in col]
    
    subtable = distances.loc[indexes, columns]
    
    print(subtable.round(1).to_latex())
    

print("normalized with 2D")
list_comparisons = ["HS/2D", "MS/2D", "3D/2D"]

indexes = [index for index in distances.index if np.sum([comparison in index for comparison in list_comparisons])==1]
columns = [col for col in distances.columns if np.sum([comparison in col for comparison in list_comparisons])==1]

subtable = distances.loc[indexes, columns]

print(subtable.round(2).to_latex())












""""""""""""""""""""""""""""""""""""
"""""""DISTANCE BETWEEN TEST DIFF PROTEINS"""""""
""""""""""""""""""""""""""""""""""""

list_col_diffs = [col for col in table.columns if "test_diff" in col] 

distances = np.zeros((len(list_col_diffs), len(list_col_diffs)))

for i, col1_diff in enumerate(list_col_diffs):
    for j, col2_diff in enumerate(list_col_diffs):

      

        # distance = np.linalg.norm(table[col1_diff] - table[col2_diff])
        distance = np.sum(table[col1_diff]!=table[col2_diff])
            
    
        
        distances[i,j] = distance
        
        
        
distances = pandas.DataFrame(distances, columns=list_col_diffs, index=list_col_diffs)    
print("distances test diff", "total:", len(table))  
distances.columns = [col.replace("test_diff_","").replace("_"," ") for col in distances.columns]    
distances.index = [col.replace("test_diff_","").replace("_"," ") for col in distances.index]      
distances = distances.astype(int)
# print(distances)

import seaborn as sns
fig, ax = plt.subplots(1, figsize=(15,15))
sns.heatmap(distances, ax=ax, cmap="Blues")
fig.tight_layout()
# # distances.to_csv(Path(figpath,"similarity_analysis", cell_type_title+"_distances_EV_expectations.csv"))

            
# index, labels, silhouette = plot_hierarchical_clustering(distances, \
#                         metric="precomputed", linkage_method="complete", threshold=0.9, \
#                         distance_matrix=False, criterion="threshold", labelsize=15, path=Path(figpath, "similarity_analysis", "HC", text_log+"test_diff_comparison_"+cell_type_title+".pdf"))







for comparison in ["HS-3D", "HS-2D"]:
    
    indexes = [index for index in distances.index if comparison in index]
    columns = [col for col in distances.columns if comparison in col]
    
    subtable = distances.loc[indexes, columns]
    
    print(subtable.astype(int).to_latex())
    
    
    
    
    
    
    
    
    
    

""""""""""""""""""""""""""""""""""""
"""""""NB COMMON PATTERNS TEST DIFF PROTEINS"""""""
""""""""""""""""""""""""""""""""""""

list_col_diffs = [col for col in table.columns if "test_diff" in col] 

distances_increase = np.zeros((len(list_col_diffs), len(list_col_diffs)))
distances_decrease = np.zeros((len(list_col_diffs), len(list_col_diffs)))

for i, col1_diff in enumerate(list_col_diffs):
    for j, col2_diff in enumerate(list_col_diffs):

      

        # distance = np.linalg.norm(table[col1_diff] - table[col2_diff])
        distance = np.sum(table[col1_diff]!=table[col2_diff])
            
    
        
        distances_increase[i,j] = np.sum( (table[col1_diff]==1) & (table[col2_diff]==1) )
        distances_decrease[i,j] = np.sum( (table[col1_diff]==-1) & (table[col2_diff]==-1) )
        
        
        
distances_increase = pandas.DataFrame(distances_increase, columns=list_col_diffs, index=list_col_diffs)    
distances_decrease = pandas.DataFrame(distances_decrease, columns=list_col_diffs, index=list_col_diffs)    

print("distances test diff")  
distances_increase.columns = [col.replace("test_diff_","").replace("_"," ") for col in distances_increase.columns]    
distances_increase.index = [col.replace("test_diff_","").replace("_"," ") for col in distances_increase.index]      

distances_decrease.columns = [col.replace("test_diff_","").replace("_"," ") for col in distances_decrease.columns]    
distances_decrease.index = [col.replace("test_diff_","").replace("_"," ") for col in distances_decrease.index]      
  
  

# import seaborn as sns
# fig, ax = plt.subplots(1, figsize=(15,15))
# sns.heatmap(distances, ax=ax)
# fig.tight_layout()
# # distances.to_csv(Path(figpath,"similarity_analysis", cell_type_title+"_distances_EV_expectations.csv"))

            
# index, labels, silhouette = plot_hierarchical_clustering(distances, \
#                         metric="precomputed", linkage_method="complete", threshold=0.9, \
#                         distance_matrix=False, criterion="threshold", labelsize=15, path=Path(figpath, "similarity_analysis", "HC", text_log+"test_diff_comparison_"+cell_type_title+".pdf"))


for comparison in ["HS-3D", "HS-2D"]:
    
    indexes = [index for index in distances_increase.index if comparison in index]
    columns = [col for col in distances_increase.columns if comparison in col]
    
    subtable = distances_increase.loc[indexes, columns]
    
    print("common increase")
    print(subtable.astype(int).to_latex())
    


    indexes = [index for index in distances_decrease.index if comparison in index]
    columns = [col for col in distances_decrease.columns if comparison in col]
    
    subtable = distances_decrease.loc[indexes, columns]
    
    print("common decrease")
    print(subtable.astype(int).to_latex())
    
    
    print("common increase three types")
    print(np.sum( (table["test_diff_"+comparison+"_THP1"]==1) & (table["test_diff_"+comparison+"_MSC_H"]==1) & (table["test_diff_"+comparison+"_Fibro"]==1) ))
    
    print("common decrease three types")
    print(np.sum( (table["test_diff_"+comparison+"_THP1"]==-1) & (table["test_diff_"+comparison+"_MSC_H"]==-1) & (table["test_diff_"+comparison+"_Fibro"]==-1) ))
            


