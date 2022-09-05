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
sys.setrecursionlimit(10000)
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
from load_data import ProteoDataset
from scipy.stats import pearsonr
from scipy.spatial.distance import pdist, squareform
import copy
# from test_proteins_differences import plot_test_diff_all_methods, plot_test_diff
from matplotlib.colors import ListedColormap, to_rgba
from plot_functions.symptoms_figures import barplot_binary_table_prevalence, barplot_results_symptoms, \
                            heatmap_results_symptoms, decision_tree_symptoms, pie_symp_asymp, \
                            barplot_results_symptoms_several_methods

#from ML_tools.pattern_analysis import extract_support
from matplotlib.colors import TwoSlopeNorm, LogNorm, ListedColormap

from statsmodels.stats.multitest import fdrcorrection, fdrcorrection_twostage
from matplotlib.colors import LinearSegmentedColormap
from utils.heatmap_utils import create_mask_for_heatmap, make_rgb_transparent
import matplotlib.patches as mpatches
from filters import get_masks_test_diff_for_heatmap
from load_data import list_comparisons

cell_type = "MSCH"
test="ProDA"


prodb = ProteoDataset()
table = prodb.load_table(cell_type)

list_comparisons = [comparison for comparison in list_comparisons if comparison[0] in protocols[cell_type] and comparison[1] in protocols[cell_type]]


for col in [col for col in table.columns if "FC" in col]:
    table[col.replace("FC", "log2_FC")] = np.log2(table[col])

cols_FC = ["log2_FC_"+test+"_"+method2+ "-" + method1 for (method1, method2) in list_comparisons]


cols_pval = ["pval_"+test+"_"+method2+ "-" + method1 for method1, method2 in list_comparisons]
cols_nb_valid = ["nb_valid_"+method for method in protocols[cell_type]]
cols_total = sum([["pval_"+test+"_"+method2+ "-" + method1, "log2_FC_"+test+"_"+method2+ "-" + method1, "nb_valid_"+method1, "nb_valid_"+method2] for (method1, method2) in list_comparisons], [])
cols_test_diff = ["test_diff_"+test+"_"+method2+ "-" + method1 for method1, method2 in list_comparisons]
cols_is_valid = [col for col in table.columns if "is_valid" in col]
cols_discrete = sum([["test_diff_"+test+"_"+method2+ "-" + method1, "is_valid_"+method1, "is_valid_"+method2] for method1, method2 in list_comparisons], [])
cols_discrete2 = ["test_diff_"+test+"_"+method2+ "-" + method1 for (method1, method2) in list_comparisons]

for col in cols_nb_valid:
    table["is_valid_"+col.split("_")[-1]] = table[col].apply(lambda s: 1 if s>0 else 0)



masks, names, index_ = get_masks_test_diff_for_heatmap(table, [cell_type], test=test)



LOGMIN = 10**(-100)

for col in cols_pval:
    table[col] = table[col].apply(lambda f: LOGMIN if f<=LOGMIN else f)
    




names.append("all")
names_to_plot = names

index_.append([])


annot_size=10
fontsize=15


##### PLOT #####

delimitations = []


for i in range(len(index_)):

    if len(index_[i])==0:
        continue
    


    if ("increased" in names[i] or "decreased" in names[i]):

        to_sort = []
        if "HS" in names[i]:
            if "2D" in names[i]:
                to_sort += [col for col in cols_FC if "HS-2D" in col]
            if "3D" in names[i]:
                to_sort += [col for col in cols_FC if "HS-3D" in col]  
                
        if "MS" in names[i]:
            if "2D" in names[i]:
                to_sort += [col for col in cols_FC if "MS-2D" in col]
            if "3D" in names[i]:
                to_sort += [col for col in cols_FC if "MS-3D" in col]

        sum_FC = table[cols_total].iloc[index_[i]][[col for col in cols_FC if col in to_sort]].sum(axis=1)
          
        ordered_index = sum_FC.argsort()
        
        if "increased" in names[i]:
            ordered_index = ordered_index[::-1]
            
            
        ordered_index = np.array(index_[i])[np.array(ordered_index)]
        
    elif names[i]=="all":

        not_in_index = [u for u in np.arange(len(table)).tolist() if u not in index_[i]]
        
        not_in_index_clustered, labels, silhouette = plot_hierarchical_clustering(table[cols_discrete].iloc[not_in_index], \
                                metric="matching", linkage_method="complete", threshold=1, \
                                distance_matrix=False, criterion="threshold")


        ordered_index = np.array(index_[i] + np.array(not_in_index)[not_in_index_clustered].tolist())

    if names[i]!="all":
        index_[-1] += [u for u in ordered_index.tolist() if u not in index_[-1]]
    
#            delimitations.append(len(index_[-1]))

        continue


    if names[i] not in names_to_plot:
        continue
    

    
    nb_by_pdf = 90
    n_win = int(len(ordered_index)/nb_by_pdf) + 1

    win_size = int(len(ordered_index)/n_win)
    
    if names[i]=="all":
        n_win = 1
        win_size = len(ordered_index)


    
    for ind in range(n_win):
        
        
        
        fontsize=20
        annot_size=10

    
        
        if ind!=n_win-1:
            indexes = ordered_index[ind*win_size:(ind+1)*win_size]
            print(ind*win_size, (ind+1)*win_size, len(ordered_index))
        else:
            indexes = ordered_index[ind*win_size:]
            print(ind*win_size, (ind+1)*win_size, "end", len(ordered_index))




    
    
        fig, ax = plt.subplots(1, figsize=(24,20))
        fig2, ax2 = plt.subplots(1, figsize=(15,20))
        
        labels = []
        patches = []

        mask = create_mask_for_heatmap(table[cols_total].iloc[indexes], cols_to_mask_list=cols_FC + cols_nb_valid,
                                       cols_to_mask_with_cond_dict={col: (col, ">= 0.01") for col in cols_pval})

        cmap = copy.copy(cm.get_cmap("Blues_r"))


        
        cmap = mpl.cm.Blues_r(np.linspace(0.00001,1,2000))
        cmap = mpl.colors.ListedColormap(cmap[:-400])
        

                    
        if names[i]!="all":
            annot = table[cols_total].iloc[indexes]
        else:
            annot = False
            
            
        if names[i]!="all":
                            
            heatmap(table[cols_total].iloc[indexes], ax=ax, xticklabels=True, yticklabels=True,
                annot = annot, annot_kws={"size":annot_size}, mask=mask, cbar_kws={'label': 'Posterior probability of error', 'shrink':0.3},#, "location":"top", "use_gridspec":False},
                fmt=".2", 
                cmap=cmap, norm=LogNorm(vmin=max(table[cols_total].iloc[indexes][cols_pval].min().min(), LOGMIN)),
                vmin=max(table[cols_total].iloc[indexes][cols_pval].min().min(), LOGMIN)
                )
        else:

            mask = create_mask_for_heatmap(table[cols_discrete].iloc[indexes], cols_to_mask_list=cols_is_valid)#,
                                          #cols_to_mask_with_cond_dict={col: (col, "= 0") for col in cols_test_diff})        


            red = [1, 0, 0]
            white = [1, 1, 1]
            alpha = 0.5
            blue = to_rgba("blue")
            cmap = ListedColormap([make_rgb_transparent(red, white, alpha), make_rgb_transparent(blue, white, alpha), "mediumseagreen"], N=3)

#                cmap = ListedColormap(["tomato", "mediumseagreen"], N=2)

#                heatmap(table[cols_discrete].iloc[indexes], ax=ax, xticklabels=True, yticklabels=True,
#                    annot = annot, annot_kws={"size":annot_size}, mask=mask, cbar_kws={'label': 'Evidence', 'shrink':0.3},#, "location":"top", "use_gridspec":False},
#                    fmt=".2", cmap=cmap, cbar=False)#, 
#                    #center=0.01, vmin=table[cols_total][cols_quantif].min().min(), vmax=table[cols_total][cols_quantif].max().max(), 
#                          
            vals = [1,0,-1]
            legends = ["increased", "no change",    "decreased"]
            # Creating 8 Patch instances
            colors = ["mediumseagreen", make_rgb_transparent(blue, white, alpha), make_rgb_transparent(red, white, alpha)]
            patches = [mpatches.Patch(color=b) for b in colors]
            labels = ['{}'.format(legends[i]) for i in range(len(vals))]
            
             
            heatmap(table[cols_discrete].iloc[indexes], ax=ax, xticklabels=True, yticklabels=True,
                annot = False, annot_kws={"size":annot_size}, cbar_kws={'label': 'p-value', 'shrink':0.3},#, "location":"top", "use_gridspec":False},
                fmt=".2", cmap=cmap, cbar=False, mask=mask)
            
             
            heatmap(table[cols_test_diff].iloc[indexes], ax=ax2, xticklabels=True, yticklabels=True,
                annot = False, annot_kws={"size":annot_size}, cbar_kws={'label': 'p-value', 'shrink':0.3},#, "location":"top", "use_gridspec":False},
                fmt=".2", cmap=cmap, cbar=False)
            
            ax2.legend(patches, labels, bbox_to_anchor=(1.04,1), loc="upper left", fontsize=20)
            

        mask = create_mask_for_heatmap(table[cols_total].iloc[indexes], cols_to_mask_list=cols_pval + cols_nb_valid,
                                       cols_to_mask_with_cond_dict={col: (col.replace("log2_FC", "pval"), ">= 0.01") for col in cols_FC})

        norm = TwoSlopeNorm(vmin=table[cols_total][cols_FC].min().min(), vcenter=0, vmax=table[cols_total][cols_FC].max().max())

        if names[i]!="all":
            annot = table[cols_total].iloc[indexes]
        else:
            annot = False
                   
            
        if names[i]!="all":
            
            heatmap(table[cols_total].iloc[indexes], ax=ax, cmap="bwr_r", xticklabels=True, yticklabels=True,
                annot = annot, mask=mask, annot_kws={"size":annot_size}, cbar_kws={'label': 'log2_FC', 'shrink':0.3},#, "location":"left","use_gridspec":False},
                fmt=".2", norm=norm 
                )
            
#            else:
#
#                mask = create_mask_for_heatmap(table[cols_discrete].iloc[indexes], cols_to_mask_list=cols_test_diff+cols_nb_valid,
#                                           cols_to_mask_with_cond_dict={col: (col.replace("log2_FC", "test_diff").replace("/","-"), "= 0") for col in cols_log2_FC})
#            
#           
#                _bwr_data = ((1.0, 0.0, 0.0), (1.0, 1.0, 1.0),  (60/256, 179/256, 113/256))
#                cmap = LinearSegmentedColormap.from_list("u", _bwr_data, N=256, gamma=1.0)
##
#                heatmap(table[cols_discrete].iloc[indexes], ax=ax, cmap=cmap, xticklabels=True, yticklabels=True,
#                    annot = False, mask=mask, annot_kws={"size":annot_size}, cbar_kws={'label': 'log2_FC', 'shrink':0.3},#, "location":"left","use_gridspec":False},
#                    fmt=".2", norm=norm)                
#        

    
        if names[i]!="all":
            ax.vlines([0,4,8,12,16,20,24], *ax.get_ylim(), colors="black", linewidth=2)
#                ax.vlines([9], *ax.get_ylim(), colors="black", linewidth=4)


        else:
            ax.vlines([0,3,6, 9, 12, 15, 18,21], *ax.get_ylim(), colors="black", linewidth=2)
            ax2.vlines([0,1,2, 3, 4, 5, 6,7,8,9,10], *ax.get_ylim(), colors="black", linewidth=2)
           
        ax.hlines(delimitations, *ax.get_xlim(), colors="black", linewidth=2)
        ax.hlines([0,len(indexes)], *ax.get_xlim(), colors="black", linewidth=2)


        grey = to_rgba("grey")
        white = to_rgba("white")
        blue = to_rgba("royalblue")

        grey1 = make_rgb_transparent(grey, white, 0.1)
        grey2 = make_rgb_transparent(grey, white, 0.2)
        grey3 = make_rgb_transparent(grey, white, 0.3)

        custom_cmap = ListedColormap([grey3, grey2, grey1, white], N=4, name="custom_cmap")

        if names[i]!="all":
            annot = table[cols_total].iloc[indexes]

            mask = create_mask_for_heatmap(table[cols_total].iloc[indexes], cols_to_mask_list=cols_FC+cols_pval)        
                    
            heatmap(table[cols_total].iloc[indexes], ax=ax, xticklabels=True, yticklabels=True,
            annot = annot, mask=mask, annot_kws={"size":annot_size}, cbar=False, fmt="", cmap=custom_cmap,
            cbar_kws={'label': 'valid values', 'shrink':0.3})

#         
            vals = [0, 1, 2, 3]
            legends = ["0/3 valid", "1/3 valid", "2/3 valid", "3/3 valid"]
            
            patches += [mpatches.Patch(color=custom_cmap(b)) for b in vals]
            labels += ['{}'.format(legends[i]) for i in range(len(vals))]
                     
        else:
            
            custom_cmap = ListedColormap([grey3, blue], N=2, name="custom_cmap")

            annot = False
        
            mask = create_mask_for_heatmap(table[cols_discrete].iloc[indexes], cols_to_mask_list=cols_test_diff)        
#                        
            heatmap(table[cols_discrete].iloc[indexes], ax=ax, xticklabels=True, yticklabels=True,
            annot = annot, mask=mask, annot_kws={"size":annot_size}, cbar=False, fmt="", cmap=custom_cmap,
            cbar_kws={'label': 'valid values', 'shrink':0.3})

            vals = [0, 1]
            legends = ["0/3 valid", "At least 1/3 valid"]


            patches += [mpatches.Patch(color="white")]
            labels += [""]
            
            patches += [mpatches.Patch(color=custom_cmap(b)) for b in vals]
            labels += ['{}'.format(legends[i]) for i in range(len(vals))]
            

        ax.legend(patches, labels, bbox_to_anchor=(1.04,1), loc="upper left", fontsize=20)
        

        ax.set_yticklabels(table["GENE_NAME"].iloc[indexes].values)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=fontsize)
        ax.set_xticklabels([col for col in ax.get_xticklabels()], rotation=45, ha="left", rotation_mode="anchor", fontsize=fontsize)

   
        if names[i]=="all":
            ax.set_yticks([])
            ax.set_yticklabels([])
            ax2.set_yticks([])
            ax2.set_yticklabels([])
            
                
            ax2.xaxis.tick_top()
            ax2.xaxis.set_label_position('top')
            ax2.set_xticklabels([col.get_text().replace("test_diff_","") for col in ax2.get_xticklabels()], rotation=45, ha="left", rotation_mode="anchor", fontsize=fontsize)
             
            
        fig.tight_layout()
        fig2.tight_layout()
        
        if names[i]!="all":
            fig.savefig(Path(figpath, "heatmaps_test_diff", cell_type, test+"_"+cell_type+"_"+names[i]+"_part_"+str(ind+1)+".pdf"))
        else:
            fig.savefig(Path(figpath, "heatmaps_test_diff", cell_type, test+"_"+cell_type+"_"+names[i]+".png"))
            fig2.savefig(Path(figpath, "heatmaps_test_diff", cell_type, test+"_"+cell_type+"_"+names[i]+"_short"+".pdf"))

            plt.close(fig2)

        plt.close(fig)




#
