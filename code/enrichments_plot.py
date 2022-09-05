
from pathlib import Path
from utils.bioinf_conversion import convert_entrez_gene_id_list_to_gene_name
from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import pandas
from background_genomes import load_background_genome_table


import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
plt.rcParams["font.family"] = "serif"

pandas.set_option('display.max_columns', None)

def save_enrichments_results(prodb, table, masks, index_, names, names_to_plot, species, savepath, test=None, display_presence=True):


    for i in range(len(names)):

        

    
        if len(index_[i])==0:
            continue
        

                
        if names[i]=="all" or names[i] not in names_to_plot:
            continue
        
        print(names[i], np.sum(masks[names[i]]))
        
        table_mask = table[masks[names[i]]].reset_index(drop=True)

        chart = prodb.get_enrichment_results(table=table_mask, species=species)
        
        chart["representation"] = chart["Fold Enrichment"].apply(lambda f: "Over" if float(f)>1 else "Under")
    
        
        proteins = OrderedDict({})

        human_table = load_background_genome_table(species=species)
        
        for categ in ["GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT"]:
            
            
            chart_to_save = save_chart(chart, category=categ)
            
            if test is not None:
                title = "enrichments_test_"+test+"_"+names[i]+"_"+categ+".csv"
            else:
                title = "enrichments_test_"+names[i]+"_"+categ+".csv"
            chart_to_save.to_csv(Path(savepath, title))
        
            
            # fig, ax = plt.subplots(1, figsize=(38,70))
            # plot_chart(chart, ax=ax, category = categ, fontsize=25)
            # fig.tight_layout()
            # fig.savefig(Path(savepath, "enrichments_test_"+names[i]+"_"+categ+".pdf")) 
            # plt.close(fig)
            
            
            
            
            
            
            
            
            

#             significative_terms = ((chart["PValue"].astype(float)<0.01)&(chart["Fold Enrichment"].astype(float)>=1.5))

#             terms = chart[significative_terms & (chart["Category"]==categ)]["Term"].values

#             over_under_representations = chart[significative_terms & (chart["Category"]==categ)]["representation"].values
            
#             for u, term in enumerate(terms):
                
#                 where_term = table[categ].astype(str).apply(lambda s: term in s)

# #                where_term_in_chart = (chart["Term"]==term)
# #                over_under = chart[where_term_in_chart]["representation"].values[0] #same

#                 over_under = over_under_representations[u]
                
#                 list_proteins = table[(masks[names[i]] & where_term)]["Gene name"].values

#                 proteins[over_under+"_"+categ.split("_")[1]+"_"+term] = list_proteins



                
#                 where_term_complete_genome = human_table[categ].apply(lambda s: term in str(s))
                
# #                human_id = human_table[where_term_complete_genome]["gene_ID"].values
# #                human_names = convert_entrez_gene_id_list_to_gene_name(human_id)
                
                
#                 human_names = human_table[where_term_complete_genome]["Gene name"].values


#                 human_names = [name for name in human_names if name in list_proteins] + [name for name in human_names if name not in list_proteins] 

# #                for name in proteins[categ.split("_")[1]+"_"+term] :
# #                    print(name in human_names)
# #                    
# #                    print(coucou)
                
                
                
# #                                
                
#                 proteins["Human_"+categ.split("_")[1]+"_"+term] = human_names
                
#                 list_protocols = [protocol for protocol in ["3D","2D","HS","MS"] if protocol in names[i]]
#                 list_cell_types = [cell_type for cell_type in ["MSCH", "THP1", "Fibro"] if cell_type in names[i]]

#                 if display_presence:
                    
                    
#                     all_presences = []
#                     for name in human_names:
#                         where_name_in_table = (table["Gene name"]==name)
# #                        
#                         if np.sum(where_name_in_table)==0:
#                             all_presences.append("")
#                             continue
                        
#                         presence = []

#                         for protocol in list_protocols:

#                             if len(list_cell_types)>0:
                            
#                                 for cell_type in list_cell_types:

#                                     presence += [protocol+"-"+cell_type+":"+ str(int(table[where_name_in_table]["nb_valid_"+protocol+"_"+cell_type].values[0]))]
#                             else:
                                
#                                 presence += [protocol +":"+ str(int(table[where_name_in_table]["nb_valid_"+protocol].values[0]))]
   
#                         char_presence = ", ".join(presence)

#                         all_presences.append(char_presence)
                        
                    
                        
#                     proteins["Presence_"+categ.split("_")[1]+"_"+term] = all_presences
                    
                    
                        
                
#         df = pandas.DataFrame.from_dict(proteins, orient="index").transpose()
#         df.to_csv(Path(savepath, "proteins_lists_GOTERM_"+names[i]+".csv"))
        
        





def save_chart(chart, category):

    chart_increased = chart[((chart["Fold Enrichment"].astype(float)>=1.5)&(chart["PValue"].astype(float)<0.01)&(chart["Category"]==category))]  
    chart_decreased = chart[((chart["Fold Enrichment"].astype(float)<=(1/1.5))&(chart["PValue"].astype(float)<0.01)&(chart["Category"]==category))]  

    chart = chart_increased

    chart.reset_index(inplace=True, drop=True)
    
    sorted_index = np.array(chart["Fold Enrichment"].astype(float)).argsort()
    
    sorted_index = sorted_index[::-1]
        
    chart = chart.iloc[sorted_index]
    
    return chart
        
#         # print(chart[["Term", "Fold Enrichment", "Count", "Pop Hits"]])

#         FC = chart["Fold Enrichment"].astype(float).values
    
#         pvals = chart["PValue"].astype(float).values
        
#         count = chart["Count"].astype(float).values
#         genome_count = chart["Pop Hits"].astype(float).values
    
#         labels = chart["Term"].values
        
#         # labels = [labels[i].split("~")[1]+" (p=%.2e)"%pvals[i] for i in range(len(labels))]
#         labels = [labels[i]+" (p=%.2e)"%pvals[i] for i in range(len(labels))]
        
#         df = pandas.DataFrame(np.array([label, FC, pvals, count, genome_count]))
       
        
            
# #        ax.set_title(names[i], fontsize=13)
        
#         if i==0 :
#             all_FCs += FC.tolist()
            
#         else:
            
#             all_FCs += [((1/fc)*(-1)) if fc>0 else -np.inf for fc in FC]

#         all_counts += count.tolist()
#         all_genome_counts += genome_count.tolist()

#         color = "mediumseagreen" if i==0 else "indianred"
#         colors += [color]*len(FC)
#         all_labels += labels
        

#         if i==1 and np.sum(FC==0)==0:

#             max_fold_change = np.max(1/FC)
#             closest_dozen = np.min(np.where(max_fold_change <= np.arange(0,1000,10)))*10

#         if i==0:
#             legend_patches = [mpatches.Patch(color="mediumseagreen")]
#             legend_labels += ["Fold enrichment"]

#         else:
#             legend_patches += [mpatches.Patch(color="indianred")]
#             legend_labels += ["Fold depletion"]
        
    
        
#     if len(chart)==0:
#         return 

#     total = chart["List Total"].values[0] 
#     total_genome = chart["Pop Total"].values[0]                    
        
#     bins = np.arange(len(all_FCs))

#     to_plot = [fc if fc!=-np.infty else -np.max(np.abs([fc for fc in all_FCs if fc!=-np.inf]))-5 for fc in all_FCs]

#     ax.barh(bins[::-1], to_plot, color=colors, alpha=0.8)

#     ax.set_yticks(bins[::-1])
#     ax.set_yticklabels(all_labels, fontsize=fontsize)
    
#     if closest_dozen is not None:
#         ax.set_xlim(left=-closest_dozen)

#     list_infty = []
#     for u in range(len(all_FCs)):
#         if all_FCs[u]==-np.infty:
#             list_infty.append(u)
    

#     for i, v in enumerate(all_FCs[::-1]):
        
#         count = int(all_counts[::-1][i])
#         genome_count = int(all_genome_counts[::-1][i])
                
# #        if v==-np.infty:
# #            ax.text(-np.max(np.abs([fc for fc in all_FCs if fc!=-np.inf])) - 12, i, "Absent ("+str(count)+", "+str(genome_count)+")", color="indianred", fontweight='bold')    
# #        else:
# #            if v>0:
# #                ax.text(v + 3, i , "%.2f"%v +"("+str(count)+", "+str(genome_count)+")", color="mediumseagreen", fontweight='bold')
# #            else:
# #                ax.text(v - 3, i , "%.2f"%(-v) +"("+str(count)+", "+str(genome_count)+")", color="indianred", fontweight='bold')
# #                

#         if v==-np.infty:
#             ax.text(-np.max(np.abs([fc for fc in all_FCs if fc!=-np.inf])) - 12, i, "Absent", color="indianred", fontweight='bold')    
#         else:
#             if v>0:
#                 # ax.text(v + 3, i , "%.2f"%v, color="mediumseagreen", fontweight='bold')
#                 ax.text(v + 3, i , "Count: "+str(count)+"/"+str(total)+", Genome : "+str(genome_count)+"/"+str(total_genome), color="black", fontweight='bold', fontsize=15)

#             else:
#                 ax.text(v - 3, i-0.1, "%.2f"%(-v), color="indianred", fontweight='bold')
                

#     ax.set_xlim(right = ax.get_xlim()[1] + 40)                  
#     ax.set_xlim(left = ax.get_xlim()[0] - 2)                  

    
#     xticks = ax.get_xticks()

#     ax.xaxis.set_major_locator(mticker.FixedLocator(xticks))
#     xticks = [tick*(-1) if tick<0 else tick for tick in xticks]


#     ax.xaxis.set_major_formatter(mticker.FixedFormatter(xticks))    

#     ax.legend(legend_patches, legend_labels, fontsize=16, loc="lower right")

# #    ax.set_xticklabels([l.get_text() for l in ax.get_xticklabels()])
#     ax.xaxis.set_tick_params(labelsize=13)
    






def plot_chart(chart, ax, category, maxrows=None, fontsize=15):

    chart_increased = chart[((chart["Fold Enrichment"].astype(float)>=1.5)&(chart["PValue"].astype(float)<0.01)&(chart["Category"]==category))]  
    chart_decreased = chart[((chart["Fold Enrichment"].astype(float)<=(1/1.5))&(chart["PValue"].astype(float)<0.01)&(chart["Category"]==category))]  


    all_FCs = []
    all_labels = []
    colors = []
    legend_labels = []
    legend_patches = []
    all_counts = []
    all_genome_counts = []
    
    closest_dozen = None

    for i, chart in enumerate([chart_increased, chart_decreased]):
        
        if len(chart)==0:
            continue
        
        if i==1:
            continue # Only keep increase
        
        if len(chart)==0:
            continue

        chart.reset_index(inplace=True)
        sorted_index = np.array(chart["Fold Enrichment"].astype(float)).argsort()
        if i==0:
            sorted_index = sorted_index[::-1]
            
        if maxrows is not None:
            sorted_index = sorted_index[:maxrows]
        chart = chart.iloc[sorted_index]
        
        # print(chart[["Term", "Fold Enrichment", "Count", "Pop Hits"]])

        FC = chart["Fold Enrichment"].astype(float).values
    
        pvals = chart["PValue"].astype(float).values
        
        count = chart["Count"].astype(float).values
        genome_count = chart["Pop Hits"].astype(float).values
    
        labels = chart["Term"].values
        
        # labels = [labels[i].split("~")[1]+" (p=%.2e)"%pvals[i] for i in range(len(labels))]
        labels = [labels[i]+" (p=%.2e)"%pvals[i] for i in range(len(labels))]
       
        
            
#        ax.set_title(names[i], fontsize=13)
        
        if i==0 :
            all_FCs += FC.tolist()
            
        else:
            
            all_FCs += [((1/fc)*(-1)) if fc>0 else -np.inf for fc in FC]

        all_counts += count.tolist()
        all_genome_counts += genome_count.tolist()

        color = "mediumseagreen" if i==0 else "indianred"
        colors += [color]*len(FC)
        all_labels += labels
        

        if i==1 and np.sum(FC==0)==0:

            max_fold_change = np.max(1/FC)
            closest_dozen = np.min(np.where(max_fold_change <= np.arange(0,1000,10)))*10

        if i==0:
            legend_patches = [mpatches.Patch(color="mediumseagreen")]
            legend_labels += ["Fold enrichment"]

        else:
            legend_patches += [mpatches.Patch(color="indianred")]
            legend_labels += ["Fold depletion"]
        
    
        
    if len(chart)==0:
        return 

    total = chart["List Total"].values[0] 
    total_genome = chart["Pop Total"].values[0]                    
        
    bins = np.arange(len(all_FCs))

    to_plot = [fc if fc!=-np.infty else -np.max(np.abs([fc for fc in all_FCs if fc!=-np.inf]))-5 for fc in all_FCs]

    ax.barh(bins[::-1], to_plot, color=colors, alpha=0.8)

    ax.set_yticks(bins[::-1])
    ax.set_yticklabels(all_labels, fontsize=fontsize)
    
    if closest_dozen is not None:
        ax.set_xlim(left=-closest_dozen)

    list_infty = []
    for u in range(len(all_FCs)):
        if all_FCs[u]==-np.infty:
            list_infty.append(u)
    

    for i, v in enumerate(all_FCs[::-1]):
        
        count = int(all_counts[::-1][i])
        genome_count = int(all_genome_counts[::-1][i])
                
#        if v==-np.infty:
#            ax.text(-np.max(np.abs([fc for fc in all_FCs if fc!=-np.inf])) - 12, i, "Absent ("+str(count)+", "+str(genome_count)+")", color="indianred", fontweight='bold')    
#        else:
#            if v>0:
#                ax.text(v + 3, i , "%.2f"%v +"("+str(count)+", "+str(genome_count)+")", color="mediumseagreen", fontweight='bold')
#            else:
#                ax.text(v - 3, i , "%.2f"%(-v) +"("+str(count)+", "+str(genome_count)+")", color="indianred", fontweight='bold')
#                

        if v==-np.infty:
            ax.text(-np.max(np.abs([fc for fc in all_FCs if fc!=-np.inf])) - 12, i, "Absent", color="indianred", fontweight='bold')    
        else:
            if v>0:
                # ax.text(v + 3, i , "%.2f"%v, color="mediumseagreen", fontweight='bold')
                ax.text(v + 3, i , "Count: "+str(count)+"/"+str(total)+", Genome : "+str(genome_count)+"/"+str(total_genome), color="black", fontweight='bold', fontsize=15)

            else:
                ax.text(v - 3, i-0.1, "%.2f"%(-v), color="indianred", fontweight='bold')
                

    ax.set_xlim(right = ax.get_xlim()[1] + 40)                  
    ax.set_xlim(left = ax.get_xlim()[0] - 2)                  

    
    xticks = ax.get_xticks()

    ax.xaxis.set_major_locator(mticker.FixedLocator(xticks))
    xticks = [tick*(-1) if tick<0 else tick for tick in xticks]


    ax.xaxis.set_major_formatter(mticker.FixedFormatter(xticks))    

    ax.legend(legend_patches, legend_labels, fontsize=16, loc="lower right")

#    ax.set_xticklabels([l.get_text() for l in ax.get_xticklabels()])
    ax.xaxis.set_tick_params(labelsize=13)
    




     
