





















from sklearn.preprocessing import scale, normalize

from mpl_toolkits import mplot3d
    
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
from test_proteins_differences import plot_test_diff_all_methods, plot_test_diff
from matplotlib.colors import ListedColormap, to_rgba
from matplotlib.colors import TwoSlopeNorm, LogNorm, ListedColormap
import matplotlib as mpl
from enrichments_plot import save_enrichments_results

from matplotlib.colors import LinearSegmentedColormap
from utils.heatmap_utils import create_mask_for_heatmap, make_rgb_transparent
import matplotlib.patches as mpatches

from utils.bioinf_conversion import convert_entrez_gene_id_list_to_gene_name
from graphic_tools.venn_diagrams import plot_venn_diagram, get_venn_sections
import seaborn as sns
import plotly.graph_objects as go
from sklearn.decomposition import PCA
import plotly.express as px
# Import the necessaries libraries
import plotly.offline as pyo
from plotly.offline import plot
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "iframe"
import plotly





# >>> from sklearn.preprocessing import MinMaxScaler
# >>> data = [[-1, 2], [-0.5, 6], [0, 10], [1, 18]]
# >>> scaler = MinMaxScaler()
# >>> print(scaler.fit(data))
# MinMaxScaler()
# >>> print(scaler.data_max_)
# [ 1. 18.]
# >>> print(scaler.transform(data))
# [[0.   0.  ]
#  [0.25 0.25]
#  [0.5  0.5 ]
#  [1.   1.  ]]
# >>> print(scaler.transform([[2, 2]]))
# [[1.5 0. ]]

protocol_dic_colors = {"MS": "mediumspringgreen", 
         "MS.6h": "mediumseagreen",
         "HS": "tomato",
         "2D": "deepskyblue",
         "3D":"royalblue"}



prodb = ProteoDataset()




cell_type = "THP1"

log=True

if log:
    text_log="log_"
else:
    text_log=""


table = prodb.load_table(cell_type)




""""""""""""""""""""""""""""""""""""
"""""""EV EXPECTATION"""""""
""""""""""""""""""""""""""""""""""""



columns = ["coeffs_"+prot for prot in protocols[cell_type]]

X = table[columns].transpose()

index = X.index
columns = X.columns

X = X.astype(float)

if log == False:
    X = 2**X


X = pandas.DataFrame(normalize(X), index=X.index, columns=X.columns)


n_components = 3

pca = PCA(n_components=n_components)
components = pca.fit_transform(X)

cumVar = pandas.DataFrame(np.cumsum(pca.explained_variance_ratio_)*100, 
                      columns=["cumVarExpl"])
expVar = pandas.DataFrame(pca.explained_variance_ratio_*100, columns=["VarExpl"])
results = pandas.concat([expVar, cumVar], axis=1).rename(index={0: "PC1", 1: "PC2", 2:"PC3"})
    
print(results)

colors = ["indianred", "royalblue", "mediumseagreen", "orange", "black"]
 
df = pandas.DataFrame(components, columns=["Component "+str(i) for i in range(1, n_components+1)])
df.index = X.index


# fig = go.Figure(data=[go.Scatter3d(df, "Component 1", y="Component 2", z=["Component 3"],
#                                    mode='markers')])

df["protocol"] = [name.replace("coeffs_","").replace("_THP1","").replace("_MSC_H","").replace("_MSC_S","").replace("_Fibro","") for name in df.index]

df["cell_type"] = [name.replace("coeffs_","").replace("3D_","").replace("2D_","").replace("HS_","").replace("MS.6h_","").replace("MS_","") for name in df.index]



if results["cumVarExpl"].loc["PC2"]==100:   

    fig = px.scatter(df, x="Component 1", y="Component 2", color='protocol', color_discrete_map=protocol_dic_colors,
                     hover_data={"protocol":True,
                                 "cell_type":True})
    fig.update_traces(marker={'size': 15})   
    plot(fig)
    
else:

    fig = px.scatter_3d(df, x="Component 1", y="Component 2", z="Component 3", color='protocol', color_discrete_map=protocol_dic_colors,
                        hover_data={"protocol":True,
                                    "cell_type":True})
    
    plot(fig)    

title = text_log+"PCA_expectation_"+cell_type+".html"
fig.write_html(Path(figpath, "similarity_analysis", "PCA", title))


# fig = px.scatter_3d(df, x="Component 1", y="Component 2", z="Component 3", color='cell_type',
#                     hover_data={"protocol":True,
#                                 "cell_type":True})

# plot(fig)

# title = text_log+"PCA_expectation_by_cell_types"+"_"+cell_type+".html"
# fig.write_html(Path(figpath, "similarity_analysis", "PCA", title))
    

""""""""""""""""""""""""""""""""""""
"""""""EV REPLICATEs"""""""
""""""""""""""""""""""""""""""""""""


columns = [col for col in experiments[cell_type] if "Lysat" not in col]

X = table[columns].transpose()

index = X.index
columns = X.columns

X = X.astype(float)

if log==False:
    X = 2**X


X = X.fillna(0)

X = pandas.DataFrame(normalize(X), index=X.index, columns=X.columns)



n_components = 3

pca = PCA(n_components=n_components)
components = pca.fit_transform(np.array(X))

cumVar = pandas.DataFrame(np.cumsum(pca.explained_variance_ratio_)*100, 
                      columns=["cumVarExpl"])
expVar = pandas.DataFrame(pca.explained_variance_ratio_*100, columns=["VarExpl"])
results = pandas.concat([expVar, cumVar], axis=1).rename(index={0: "PC1", 1: "PC2", 2:"PC3"})
    
print(results)

colors = ["indianred", "royalblue", "mediumseagreen", "orange", "black"]

  
df = pandas.DataFrame(components, columns=["Component "+str(i) for i in range(1, n_components+1)])
df.index = X.index
#
import plotly.io as pio
# import plotly.express as px
pio.renderers.default='plotly_mimetype+notebook_connected'
# pio.renderers.default = "png"


from plotly.offline import plot, iplot, init_notebook_mode
import plotly.graph_objs as go
init_notebook_mode(connected=True)
#

# fig = go.Figure(data=[go.Scatter3d(df, "Component 1", y="Component 2", z=["Component 3"],
#                                    mode='markers')])

df["protocol"] = [name.replace("LFQI ","").replace("_THP1","").replace("_MSC_H","").replace("_MSC_S","").replace("_Fibro","")[:-2] for name in df.index]

df["cell_type"] = [name.replace("3D_ ","").replace("2D_ ","").replace("HS_ ","").replace("MS.6h_ ","").replace("MS_ ","").split("_")[-1] for name in df.index]




fig = px.scatter_3d(df, x="Component 1", y="Component 2", z="Component 3", color='protocol', color_discrete_map=protocol_dic_colors,
                        hover_data={"protocol":True,
                                    "cell_type":True})

plot(fig)
title = text_log+"PCA_replicates_"+"_"+cell_type+".html"
fig.write_html(Path(figpath, "similarity_analysis", "PCA", title))


df["name"] = df.index

# if len(cell_types_list)>1:

#     fig = px.scatter_3d(df, x="Component 1", y="Component 2", z="Component 3", color='cell_type',
#                         hover_data={"protocol":True,
#                                     "cell_type":True})
#     plot(fig)
    
#     title = text_log+"PCA_replicates_by_cell_types_"+"_".join(sorted(cell_types_list))+".html"
#     fig.write_html(Path(figpath, "similarity_analysis", "PCA", title))
      
















