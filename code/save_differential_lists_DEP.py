





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
from filters import get_masks_test_diff, get_masks_test_diff_only_pair
from load_data import list_comparisons
from enrichments_plot import save_enrichments_results
import os
from data_infos import cell_types
import plotly.graph_objects as go


figpath = Path(figpath, "separated_cell_types")



plot_enrichments=True

prodb = ProteoDataset(reset=False)


test="DEP"


for cell_type in cell_types:


    if "mapH" in cell_type:
        print(cell_type)
        
        continue

    table = prodb.load_table(cell_type, load_lysat=True)
    
    
    #################### For STRING DB
    masks, names, index_, cols = get_masks_test_diff_only_pair(table, [cell_type], test=test)
    
    
    """ Save diff list """
    
    for k, name in enumerate(masks):
        
        mask = masks[name]
        
        subtable = table[mask].copy()
        subtable.reset_index(inplace=True)
            
        col_FC = cols[k].replace("test_diff", "FC")
        col_pval = cols[k].replace("test_diff", "pval")
    
        
        cond1 = name.split("_")[1]
        cond2 = name.split("_")[3]
    
        
        LFQI_cols = [col for col in table.columns if "LFQI" in col and cond1 in col and "imputation" not in col] + \
            [col for col in table.columns if "LFQI" in col and cond2 in col and "imputation" not in col]

        LFQI_imputed_cols = [col for col in table.columns if "LFQI" in col and cond1 in col and "imputation" in col] + \
            [col for col in table.columns if "LFQI" in col and cond2 in col and "imputation" in col]

        if not os.path.isdir(Path(figpath, "3. lists_differential_proteins_DEP", cell_type)):
            os.mkdir(Path(figpath, "3. lists_differential_proteins_DEP", cell_type))
            
            
        to_save = subtable[["ENTREZ_GENE_ID", "GENE_NAME", "UNIPROT_ACCESSION", col_pval, col_FC]+LFQI_cols+LFQI_imputed_cols]
        to_save = to_save.sort_values(by=col_FC, ascending=False)
        
        header = ["Gene_ID", "Gene_name", "Protein_ID", "pval", "log2_Fold change"]+ \
            ["log2_"+col for col in LFQI_cols]+["log2_"+col for col in LFQI_imputed_cols]
    
        to_save.to_csv(Path(figpath, "3. lists_differential_proteins_DEP", cell_type, "list_diff_"+cell_type+"_"+test+"_"+name+".csv"), 
                       header=header, index=None, sep=' ')
        
        #### Volcano plot ####
        
        df = subtable.copy()
        df["Gene_ID"] = df["ENTREZ_GENE_ID"]
        df["log2_FC"] = np.log2(df[col_FC])
        df["pvalue"] = np.log2(df[col_pval])
        df = df[["Gene_ID", "log2_FC", "pvalue"]]
        
        fig = go.Figure()
        trace1 = go.Scatter(
         x=df['log2_FC'],
         y=df['pvalue'],
         mode='markers',
         hovertext=list(df["Gene_ID"])
        )
        
        fig.add_trace(trace1)
        fig.update_layout(title=name.replace("_", " "), xaxis={'title':'Log 2 Fold change'}, yaxis={"title":"Log 2 p-value"})        
        fig.write_html(Path(figpath, "3. lists_differential_proteins_DEP", cell_type, "volcano_plot_"+cell_type+"_"+test+"_"+name+".html"))
        
       
        
        
        
    ####################################"
    
    
        
    
    species = dic_species[cell_type]
        
    if plot_enrichments:
        
    
        if not os.path.isdir(Path(figpath, "4. enrichments_differential_proteins_DEP", cell_type)):
            os.mkdir(Path(figpath, "4. enrichments_differential_proteins_DEP", cell_type))
    
        
        save_enrichments_results(prodb, table, masks, index_, names, names_to_plot=names, species=species, savepath=Path(figpath, "4. enrichments_differential_proteins_DEP", cell_type), test=test)
