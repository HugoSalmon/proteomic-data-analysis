

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import patches

import seaborn as sns

import pandas



from matplotlib import rcParams
rcParams['font.family'] = "sans-serif"















def plot_violin(ax, list_data, list_labels, list_colors, percentiles=[0,100], ratio=0):
    
    l_data = list_data.copy()
    
    for k in range(len(l_data)):
        
        
        l_data[k] = [d for d in l_data[k] if d <= np.percentile(l_data[k],percentiles[1])]

        l_data[k] = [d for d in l_data[k] if d >= np.percentile(l_data[k],percentiles[0])]
        
        if len(l_data[k])==0:
            l_data[k] = [0]
            print("empty data, filled by [0]")
            

    parts = ax.violinplot(dataset = l_data,  showmeans=False, showmedians=False,
                          showextrema=False)
 

    ipc = 0
    for ipc, pc in enumerate(parts['bodies']):

        pc.set_facecolor(list_colors[ipc])
        pc.set_edgecolor(list_colors[ipc])
        pc.set_alpha(0.7)  
        


    for ix in range(len(l_data)) : 
        data = l_data[ix]
        
        
        

        mean = np.mean(data)
        std = np.std(data)
        
    

        low_val = mean - 1.96*(std)
        high_val = mean + 1.96*(std)
        
        
        ax.scatter([1+ix], mean, marker='o', color='white', s=30, zorder=3)

        ax.plot([1+ix]*2, [low_val, high_val], color='k', linestyle='-', lw=5)
        ax.plot([1+ix-0.1,1+ix+0.1, ], [low_val]*2, color='k', linestyle='-', lw=5)
        ax.plot([1+ix-0.1,1+ix+0.1, ], [high_val]*2, color='k', linestyle='-', lw=5)
                 

#        if np.sum(np.array(data)<=0)==len(data):
#            ax.plot([1+ix]*2, [low_val, mean], color='k', linestyle='-', lw=5)
#            ax.plot([1+ix-0.1,1+ix+0.1, ], [low_val]*2, color='k', linestyle='-', lw=5)
#        elif np.sum(np.array(data)>=0)==len(data):
#            ax.plot([1+ix]*2, [mean, high_val], color='k', linestyle='-', lw=5)
#            ax.plot([1+ix-0.1,1+ix+0.1, ], [high_val]*2, color='k', linestyle='-', lw=5)
#        else:
#            ax.plot([1+ix]*2, [low_val, high_val], color='k', linestyle='-', lw=5)
#            ax.plot([1+ix-0.1,1+ix+0.1, ], [low_val]*2, color='k', linestyle='-', lw=5)
#            ax.plot([1+ix-0.1,1+ix+0.1, ], [high_val]*2, color='k', linestyle='-', lw=5)
#                            
##

    ax.set_xticklabels(np.array(list_labels))
    
    ax.set_xticks(np.arange(len(l_data))+1)
    
    xlims = ax.get_xlim()

    ax.set_xlim([xlims[0]-(ratio*(xlims[1]-xlims[0])), xlims[1]+(ratio*(xlims[1]-xlims[0]))])
   
#    
#    if cut:
#        percentiles_min = [np.percentile(list_data[k],1) for k in range(len(list_data))]
#        percentiles_max = [np.percentile(list_data[k],99) for k in range(len(list_data))]
#       
#        ax.set_ylim([np.max(percentiles_min), np.min(percentiles_max)])
#    
#
#
#def plot_violin2(list_data, list_labels=None, base_color="lightblue", confidence_color="dodgerblue", ax=None):
#
#    
#    if ax is None:
#        fig, ax = plt.subplots(1, figsize=(10,8))
#    
#
#    width = 1
#    
#    step = 2*width
#    
#    xticks = np.arange(len(list_data))*step
#    
#
#    violins = ax.violinplot(list_data, showmedians=False, showextrema=False, points=300, widths=width, positions = xticks)
#
#    
#    
#    
#    if list_labels is not None:
#
#        ax.set_xticklabels(np.array(list_labels))
#        ax.set_xticks(xticks)
#        
#    else:
#        ax.set_xticks([])
#        ax.set_xticklabels([])
#        
#    
#    if len(list_data)==1:
#        ax.set_xlim([xticks[0]-width*2,xticks[0]+width*2])
#
#        
#
#    k=0
#
#    for i, violin in enumerate(violins["bodies"]):
#
#        path = violin.get_paths()[0]
#        
#        data = list_data[i]
#        conf_bounds = [np.mean(data) - 1.96*np.std(data), np.mean(data) + 1.96*np.std(data)]
#
#        r = patches.Rectangle( (k-width,conf_bounds[0]), width=2*width, height=conf_bounds[1]-conf_bounds[0], color=confidence_color)
#        
#        r.set_clip_path(path, transform=ax.transData)
#        ax.add_patch(r)
#        
#        ax.plot([k,k],[conf_bounds[0],conf_bounds[1]], color="darkblue")
#        ax.plot([k-0.1*width, k+0.1*width], [np.mean(data), np.mean(data)], color="darkblue")
#        
#        if np.sum(data>0)>0:
#            ax.plot([k-0.2*width, k+0.2*width], [conf_bounds[1], conf_bounds[1]], color="darkblue")
#        
#        if np.sum(data<0)>0:
#            ax.plot([k-0.2*width, k+0.2*width], [conf_bounds[0], conf_bounds[0]], color="darkblue")
#
#        
#        k+=step 


if __name__=="__main__":
    
    
    u = np.random.randn(300)
    
    v = 5 + np.random.randn(300)
    
    plot_violin([u,v], ["data1","data2"], "lightblue", "dodgerblue")