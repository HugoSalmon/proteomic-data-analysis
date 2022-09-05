
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import numpy as np
from scipy import special
from scipy import stats
from matplotlib import colors

def create_mask_for_heatmap(df, cols_to_mask_list=[], cols_to_mask_with_cond_dict={}):
    
    mask = np.array(df.applymap(lambda f: False))
    
    columns = list(df.columns)
    
    for col in cols_to_mask_list:

        where_col = np.where(np.array(columns)==col)[0]

        mask[:,where_col] = True
        
    for col, cond in cols_to_mask_with_cond_dict.items():
        
        cond_col = cond[0]
        constraint = cond[1]
        
        sign = constraint.split(" ")[0]
        value = float(constraint.split(" ")[1])
        
        if sign=="=":
            to_mask = (df[cond_col]==value)
        elif sign=="<=":
            to_mask = (df[cond_col]<=value)
        elif sign=="<":
            to_mask = (df[cond_col]<value)
        elif sign==">=":
            to_mask = (df[cond_col]>=value) 
        elif sign==">":
            to_mask = (df[cond_col]>value)            
            
        mask[np.array(to_mask), columns.index(col)] = True

    return mask




def make_rgb_transparent(rgb, bg_rgb, alpha):
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(rgb, bg_rgb)]


def make_color_transparent(color, alpha):
     
    rgb=list(colors.to_rgba(color))
    
    rgb[-1] = alpha
    
    return tuple(rgb)
    # bg_rgb = colors.to_rgba(bg_color)
    # return [alpha * c1 + (1 - alpha) * c2
    #         for (c1, c2) in zip(rgb, bg_rgb)]


    
    
def compute_by_chance_threshold_protein_number(n_proteins, n_comparisons, n_found):
    
    p = 0.5
    
    p_all = p**n_comparisons
    
    probas = []
    
    for k in range(n_proteins):
        
#        pi = special.binom(n_proteins, k) * p_all**k * (1 - p_all)**(n_proteins-k)
        
        pi = stats.binom.pmf(k, p=p_all, n=n_proteins)
        
        probas.append(pi)
        
    probas = np.array(probas)
        
    distribution_mode = np.argmax(probas)
    
    if n_found >= distribution_mode:
        
        pval = 1 - stats.binom.cdf(n_found, p=p_all, n=n_proteins)
        
    else:
        pval = stats.binom.cdf(n_found, p=p_all, n=n_proteins)
        
    print("pval", pval)
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(1)
    bins = np.arange(n_proteins)
    ax.plot(stats.binom.pmf(bins, p=p_all, n=n_proteins))

    return np.array(probas), pval




if __name__=="__main__":
    
    compute_by_chance_threshold_protein_number(n_proteins=5000, n_comparisons=1, n_found=40)
    
    