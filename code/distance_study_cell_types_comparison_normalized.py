








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
from load_data import ProteoDataset, human_cell_types
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





log = True
if log:
    text_log="log_"
else:
    text_log=""




cell_type_title = "_".join(sorted(human_cell_types))






""""""""""""""""""""""""""""""""""""
"""""""COUNT DISTANCE VECTORS EV EXPECTATION"""""""
""""""""""""""""""""""""""""""""""""

table = prodb.load_human_table(load_lysat=False)
table = prodb.normalize_GAPDH(table)


list_protocols = sum([[protocol+"_"+cell_type for protocol in protocols[cell_type]] for cell_type in human_cell_types],[])
    

distances = np.zeros((len(list_protocols), len(list_protocols)))
for i, protocol1 in enumerate(list_protocols):
    for j, protocol2 in enumerate(list_protocols):
        
        if protocol1 == protocol2:
            distance = 0
            continue        
        col_protocol1 = "normalized_coeffs_"+protocol1
        
        col_protocol2 = "normalized_coeffs_"+protocol2
        
        if log:

            distance = np.linalg.norm(table[col_protocol1] - table[col_protocol2])
            
        else:
            distance = np.linalg.norm(2**table[col_protocol1] - 2**table[col_protocol2])

        
        distances[i,j] = distance
        
distances = pandas.DataFrame(distances, columns=list_protocols, index=list_protocols)    
print("distances EV expectation")   
distances.index = [col.replace("_"," ") for col in distances.index]    
distances.columns = [col.replace("_"," ") for col in distances.columns]    
print(distances.round(2).to_latex())

# distances.to_csv(Path(figpath,"similarity_analysis", cell_type_title+"_distances_EV_expectations.csv"))
            
index, labels, silhouette = plot_hierarchical_clustering(distances, \
                        metric="precomputed", linkage_method="complete", threshold=0.9, \
                        distance_matrix=False, criterion="threshold", labelsize=15, path=Path(figpath, "similarity_analysis", "HC", text_log+cell_type_title+"_HC_EV_expectations_normalized.pdf"))

    
    
    
    


""""""""""""""""""""""""""""""""""""
"""""""COUNT DISTANCE VECTORS EV MEANS WITH DISCARDING"""""""
""""""""""""""""""""""""""""""""""""
    
table = prodb.load_human_table(load_lysat=False)
table = prodb.normalize_GAPDH(table)

table = prodb.discard_missing_values(table)


cell_type_title = "_".join(sorted(human_cell_types))



list_protocols = sum([[protocol+"_"+cell_type for protocol in protocols[cell_type]] for cell_type in human_cell_types],[])
    

distances = np.zeros((len(list_protocols), len(list_protocols)))
for i, protocol1 in enumerate(list_protocols):
    for j, protocol2 in enumerate(list_protocols):
        
        if protocol1 == protocol2:
            distance = 0
            continue        
        col_protocol1 = "normalized_empirical_mean_"+protocol1
        
        col_protocol2 = "normalized_empirical_mean_"+protocol2
        
        if log:

            distance = np.linalg.norm(table[col_protocol1] - table[col_protocol2])
            
        else:
            distance = np.linalg.norm(2**table[col_protocol1] - 2**table[col_protocol2])

        
        distances[i,j] = distance
        
distances = pandas.DataFrame(distances, columns=list_protocols, index=list_protocols)    
print("distances EV expectation")   
distances.index = [col.replace("_"," ") for col in distances.index]    
distances.columns = [col.replace("_"," ") for col in distances.columns]    
print(distances.round(2).to_latex())

# distances.to_csv(Path(figpath,"similarity_analysis", cell_type_title+"_distances_EV_expectations.csv"))
            
index, labels, silhouette = plot_hierarchical_clustering(distances, \
                        metric="precomputed", linkage_method="complete", threshold=0.9, \
                        distance_matrix=False, criterion="threshold", labelsize=15, path=Path(figpath, "similarity_analysis", "HC", text_log+cell_type_title+"_HC_EV_empirical_mean_normalized.pdf"))



