
import numpy as np
import pandas
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.font_manager
rcParams['font.family'] = 'serif'


from seaborn import heatmap


from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, dendrogram, set_link_color_palette
from sklearn.metrics import silhouette_score



import matplotlib.patches as mpatches

ordered_colors = ['indianred','royalblue','mediumseagreen','skyblue','darkturquoise','mediumvioletred', 'darksalmon','blue']

def compute_hierarchical_tree(data, metric, linkage_method="ward"):

    n = len(data)

    if metric!="precomputed":
        flattened_distance_matrix = pdist(data, metric = metric)
        square_distance_matrix = squareform(flattened_distance_matrix)
    else:
#        unsquareform = lambda a: a[np.nonzero(np.triu(a))]
#        flattened_distance_matrix = unsquareform(data)
        flattened_distance_matrix = squareform(data) # function works vice versa
        square_distance_matrix = data

    tree = linkage(flattened_distance_matrix, method=linkage_method)

    tree_order = get_tree_order(tree, n + n-2)
    
    
    ordered_dist_matrix = np.zeros((n,n))    

    a,b = np.triu_indices(n,k=1)
    ordered_dist_matrix[a,b] = square_distance_matrix[ [tree_order[i] for i in a], [tree_order[j] for j in b]]
    ordered_dist_matrix[b,a] = ordered_dist_matrix[a,b]
    
    return tree, tree_order, ordered_dist_matrix

def get_tree_order(linkage_matrix, cluster_index):
  
    n = len(linkage_matrix) + 1

    if cluster_index < n:
        return [cluster_index]
    else:
        left = int(linkage_matrix[cluster_index-n, 0])
        right = int(linkage_matrix[cluster_index-n, 1])
        return (get_tree_order(linkage_matrix, left) + get_tree_order(linkage_matrix, right))




#def plot_results(table, ax):
#
#    fig, ax = plt.subplots(2, 1, figsize=(16,13))
#    heat_map = heatmap(table[index_], ax=ax[0])
#    dend = plot_dendrogram(ax[1], res_linkage, table.index, threshold=0.7)
#    ax[0].set_aspect("equal")
#    ax[0].xaxis.set_tick_params(labelsize=2)
#    ax[0].yaxis.set_tick_params(labelsize=2)
#
#
#

def plot_hierarchical_clustering(table, metric, linkage_method, criterion="maxclust", threshold=0.7, max_clust=2, labelsize=20, distance_matrix=False, path=None):

    data = np.array(table)

    if metric!="precomputed":
        flattened_distance_matrix = pdist(data, metric = metric)
        square_distance_matrix = squareform(flattened_distance_matrix)
    else:
#        unsquareform = lambda a: a[np.nonzero(np.triu(a))]
#        flattened_distance_matrix = unsquareform(data) 
        flattened_distance_matrix = squareform(data) # function works vice versa
        square_distance_matrix = data

       
    title = metric+" distance, "+" agglomerative clustering, "+linkage_method
        
    linkage_matrix, order, ordered_dist_matrix = compute_hierarchical_tree(data, metric=metric, linkage_method=linkage_method)
    
    if distance_matrix:
        fig, ax = plt.subplots(2, 1, figsize=(20,10))
    
        plot_heatmap(ordered_dist_matrix, ax=ax[0])
        ax_dend = ax[1]
        
    else:
        fig, ax = plt.subplots(1, figsize=(15,8))
        ax_dend = ax
    
    corresp_index = {str(table.index[i]):i for i in range(len(table.index))}
    
    dend = plot_dendrogram(ax_dend, linkage_matrix, table.index, threshold=threshold)
    
    ax_dend.set_xticklabels(ax_dend.get_xticklabels(), fontsize=labelsize, rotation=45, rotation_mode="anchor", ha="right")

    ordered_names_dendrogram = [i.get_text() for i in ax_dend.get_xticklabels()]
    ordered_index_dendrogram = [corresp_index[name] for name in ordered_names_dendrogram]
    
    
#    
#
#    fig.suptitle(title)
#    fig.savefig(title.replace(" ","_")+".png")
#    
    if criterion=="threshold":
        labels = fcluster(linkage_matrix, t=threshold*np.max(linkage_matrix[:,2]), criterion="distance")
    elif criterion=="maxclust":
        labels = fcluster(linkage_matrix, t=max_clust, criterion=criterion)
        
    labels-=1

    

    if len(np.unique(labels))>1 and len(np.unique(labels))<len(table.index):
        silhouette_hierarchique = silhouette_score(square_distance_matrix, labels, metric='precomputed')
        print("silhouette score : "+str(silhouette_hierarchique))
    else:
        silhouette_hierarchique = None
    
    
    fig.tight_layout()
    
    
    if path is not None:
        fig.savefig(path)


    return ordered_index_dendrogram, labels, silhouette_hierarchique





def plot_heatmap(table, ax):

    heatmap(table, ax=ax)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")




def plot_dendrogram(ax, res_linkage, results, threshold=0.7):
    
    set_link_color_palette(ordered_colors)
    dend = dendrogram(res_linkage, ax=ax, labels=results, above_threshold_color='k', color_threshold=threshold*np.max(res_linkage[:,2]))
    
    return dend





if __name__=="__main__":     
    
    
    cluster1 = 10 + np.random.randn(5)
    cluster2 = np.random.randn(10)
    data = np.concatenate([cluster1, cluster2]).reshape(-1,1)
    individuals = ["I"+str(i) for i in range(15)]
    table = pandas.DataFrame(data, index=individuals)

    

    """ clustering hierarchique """

    index_, labels, silhouette = plot_hierarchical_clustering(table, metric="euclidean", \
                                                              linkage_method="average", \
                                                              distance_matrix=True, threshold=0.7, \
                                                              criterion="threshold")

    