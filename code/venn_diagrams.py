













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
from data_infos import protocols, experiments, table_columns, cell_types, dic_protocol_corresp, protocols_lysat
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
from test_proteins_differences import plot_test_diff_all_methods, plot_test_diff
from matplotlib.colors import ListedColormap, to_rgba
from matplotlib.colors import TwoSlopeNorm, LogNorm, ListedColormap
import matplotlib as mpl
from enrichments_plot import save_enrichments_results

from matplotlib.colors import LinearSegmentedColormap
from utils.heatmap_utils import create_mask_for_heatmap, make_rgb_transparent, make_color_transparent
import matplotlib.patches as mpatches

from pyensembl import EnsemblRelease
data = EnsemblRelease(77)
from utils.bioinf_conversion import convert_entrez_gene_id_list_to_gene_name
from graphic_tools.venn_diagrams import plot_venn_diagram, get_venn_sections

pandas.set_option('display.expand_frame_repr', False)
#    - Un diagramme de Venn avec les conditions 300RPM 3h, 300RPM 6h, 500RPM 3h, 2D starvation pour les EVs, et aussi pour les lysats d'un autre côté (par curiosité). 

cell_type = "THP1"

protocol_dic_colors = {"MS": "mediumspringgreen", 
         "MS.6h": "mediumseagreen",
         "HS": "tomato",
         "2D": "deepskyblue",
         "3D":"royalblue"}


def get_list_proteins_from_protocols_list(table, list_protocols):
    
    lists_proteins = []
    
    for protocol in list_protocols:
        
        where_valid = (table["nb_valid_"+protocol] >= 1)
        
        list_proteins = table[where_valid]["ID"].values.tolist()
        
        lists_proteins.append(list_proteins)
        
    return lists_proteins


prodb = ProteoDataset()
table = prodb.load_table(cell_type, load_lysat=True)


list_protocols = protocols[cell_type]

lists_proteins = get_list_proteins_from_protocols_list(table, list_protocols)

plot_venn_diagram(lists_proteins, list_protocols, list_colors = [make_color_transparent(protocol_dic_colors[prot],0.4) for prot in list_protocols], savepath=Path(figpath, "venn_diagrams", "venn_diagram_EV_"+cell_type+".png"))

sections = get_venn_sections(lists_proteins, list_protocols)  


#verif   

# where = ((table["nb_valid_MS"]>=1) & (table["nb_valid_2D"]>=1))


counts = np.zeros((len(list_protocols), len(list_protocols)))
for i, protocol1 in enumerate(list_protocols):
    for j, protocol2 in enumerate(list_protocols):
        sum_common = np.sum([tup[1] for tup in sections if protocol1 in tup[0] and protocol2 in tup[0]])
        
        counts[i,j] = sum_common
        
counts = pandas.DataFrame(counts, columns=list_protocols, index=list_protocols)    
print("counts")        
print(counts)
  

########"" Lysat #########"

if cell_type=="THP1":
    
    protocol_dic_colors = {"Lysat MS": "mediumspringgreen", 
              "Lysat MS.6h": "mediumseagreen",
              "Lysat HS": "tomato",
              "Lysat 2D": "deepskyblue",
              "Lysat 3D":"royalblue"}
    
    
    list_protocols = protocols_lysat[cell_type]
    
    lists_proteins = get_list_proteins_from_protocols_list(table, list_protocols)
    
    plot_venn_diagram(lists_proteins, list_protocols, list_colors=[make_color_transparent(protocol_dic_colors[prot],0.4) for prot in list_protocols], savepath=Path(figpath, "venn_diagrams", "venn_diagram_Lysat_THP1.png"))
    
    lysat_sections = get_venn_sections(lists_proteins, list_protocols)
    
    
    lysat_counts = np.zeros((len(list_protocols), len(list_protocols)))
    for i, protocol1 in enumerate(list_protocols):
        for j, protocol2 in enumerate(list_protocols):
            sum_common = np.sum([tup[1] for tup in lysat_sections if protocol1 in tup[0] and protocol2 in tup[0]])
             
            lysat_counts[i,j] = sum_common
             
    lysat_counts = pandas.DataFrame(lysat_counts, columns=list_protocols, index=list_protocols) 
    print("lysat counts")           
    print(lysat_counts)     
    
    where = ((table["nb_valid_Lysat MS"]>=1) & (table["nb_valid_Lysat 2D"]>=1))
    print(where.sum())
    
        
        
        
    
    
    ########"" Compare Lysat #########"
    
    protocol_dic_colors = {"Lysat MS": "mediumseagreen", 
             "Lysat MS.6h": "mediumseagreen",
             "Lysat HS": "mediumseagreen",
             "Lysat 2D": "mediumseagreen",
             "Lysat 3D":"mediumseagreen",
             
             "MS": "tomato", 
                      "MS.6h": "tomato",
                      "HS": "tomato",
                      "2D": "tomato",
                      "3D":"tomato"}
        
    for protocol in protocols[cell_type]:
    
    
    
        list_protocols = [protocol, "Lysat "+protocol]
    
        lists_proteins = get_list_proteins_from_protocols_list(table, list_protocols)
    
        plot_venn_diagram(lists_proteins, list_protocols, list_colors=[make_color_transparent(protocol_dic_colors[prot],0.4) for prot in list_protocols], savepath=Path(figpath, "venn_diagrams", "venn_diagram_comparison_EV_Lysat_"+protocol+".png"))
    
        lysat_sections = get_venn_sections(lists_proteins, list_protocols)
    
        
        lysat_counts = np.zeros((len(list_protocols), len(list_protocols)))
        for i, protocol1 in enumerate(list_protocols):
            for j, protocol2 in enumerate(list_protocols):
                sum_common = np.sum([tup[1] for tup in lysat_sections if protocol1 in tup[0] and protocol2 in tup[0]])
                 
                lysat_counts[i,j] = sum_common
                 
        lysat_counts = pandas.DataFrame(lysat_counts, columns=list_protocols, index=list_protocols) 
        print("lysat counts")           
        print(lysat_counts)     
    
