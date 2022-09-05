





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
from data_infos import protocols, experiments, table_columns, dic_species
from pathlib import Path
from paths import figpath, savepath
import itertools
import sys
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
from load_data import ProteoDataset
from scipy.stats import pearsonr
from scipy.spatial.distance import pdist, squareform
import copy
from matplotlib.colors import ListedColormap, to_rgba
from plot_functions.symptoms_figures import barplot_binary_table_prevalence, barplot_results_symptoms, \
                            heatmap_results_symptoms, decision_tree_symptoms, pie_symp_asymp, \
                            barplot_results_symptoms_several_methods
from enrichments_plot import save_enrichments_results

#from ML_tools.pattern_analysis import extract_support
from matplotlib.colors import TwoSlopeNorm, LogNorm, ListedColormap

from statsmodels.stats.multitest import fdrcorrection, fdrcorrection_twostage
from matplotlib.colors import LinearSegmentedColormap
from utils.heatmap_utils import create_mask_for_heatmap, make_rgb_transparent
import matplotlib.patches as mpatches
from filters import get_masks_test_diff, get_masks_on
from load_data import list_comparisons
from enrichments_plot import save_enrichments_results
import os

plot_enrichments=True

figpath = Path(figpath, "grouped_cell_types")



prodb = ProteoDataset(reset=False)
table = prodb.load_human_table(load_lysat=False)



#################### For STRING DB
masks, names, index_ = get_masks_on(table, ["THP1", "Fibro", "MSCH", "MSCSmapH"])


""" Save presence list """

for k, name in enumerate(masks):
    
    mask = masks[name]
    
    
    
    subtable = table[mask].copy()
    subtable.reset_index(inplace=True)
    
    protocol = name.split("_")[1]
    
    cell_types = name.split("_")[-1]
    
    coef_cols = []
    name_cols = []
    
    for cell_type in cell_types.split("-"):
                
        coef_cols += [col for col in subtable.columns if protocol in col and "LFQI" in col and "imputation" not in col]
        
    
    # coef_cols = [names[k].replace("on_","coeffs_")
    # name_cols = 
    
    

    subtable[["ENTREZ_GENE_ID", "GENE_NAME"]+coef_cols].to_csv(Path(figpath, "5. lists_present_proteins_in_at_least_one_replicate_grouped_cell_types", "list_present_proteins_"+name+".txt"), header=["Gene_ID", "Gene_name"]+coef_cols, index=None, sep=' ')
    
####################################"


    

species = "homo_sapiens"
    
if plot_enrichments:
    

    
    save_enrichments_results(prodb, table, masks, index_, names, names_to_plot=names, species=species, savepath=Path(figpath, "6. enrichments_present_proteins_in_at_least_one_replicate_grouped_cell_types"))
