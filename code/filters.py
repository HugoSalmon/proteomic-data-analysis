

import numpy as np
from data_infos import protocols, table_columns, protocols_lysat
from collections import OrderedDict
from load_data import list_comparisons
import itertools


def get_mask_valid(table, conds_equal={}, conds_sup={}):
    
    mask = np.array([True]*len(table))
    
    for method, value in conds_equal.items():

        mask = mask & np.array(table["nb_valid_"+method]==value)

    for method, value in conds_sup.items():

        mask = mask & np.array(table["nb_valid_"+method]>=value)
                
    return mask




def get_masks_on(table, list_cell_types):


    masks = OrderedDict({})
        
    if len(list_cell_types)==1:
                
        cell_type = list_cell_types[0]

        for protocol in protocols[cell_type] + protocols_lysat[cell_type]:
        
            masks.update({"present_"+protocol : get_mask_valid(table, conds_sup={protocol:1})})

    else:

        
        combi_cell_types = list(itertools.combinations(list_cell_types, 4)) + list(itertools.combinations(list_cell_types, 3)) + list(itertools.combinations(list_cell_types, 2))
                        
        for protocol in ["HS","2D","3D"]:

            for combi in combi_cell_types:    
                
                title = "-".join(sorted(list(combi)))
            
                masks.update({"present_"+protocol+"_"+title : get_mask_valid(table, conds_sup={protocol+"_"+cell_type:1 for cell_type in combi})})

    names = list(masks.keys())

    index_ = [list(np.where(mask)[0]) for k,mask in masks.items()]    
        
    return masks, names, index_





def get_masks_test_diff_only_pair(table, list_cell_types, test):

    cols = []
    masks = OrderedDict({})


    if len(list_cell_types)==1:
        
        cell_type = list_cell_types[0]
            
                        
        for (method1, method2) in list_comparisons:
            
            if method1 not in protocols[cell_type]+protocols_lysat[cell_type] or method2 not in protocols[cell_type]+protocols_lysat[cell_type]:
                continue 
                        
            masks.update({"increased_"+method2.replace(" ","_")+"_versus_"+method1.replace(" ","_"): np.array(table['test_diff_'+test+"_"+method2+"-"+method1]==1).astype(bool)})
            masks.update({"decreased_"+method2.replace(" ","_")+"_versus_"+method1.replace(" ","_"): np.array(table['test_diff_'+test+"_"+method2+"-"+method1]==-1).astype(bool)})
            
            cols.append("test_diff_"+test+"_"+method2+"-"+method1)
            cols.append("test_diff_"+test+"_"+method2+"-"+method1)
            
    else:
                
        combi_cell_types = list(itertools.combinations(list_cell_types, 4)) + list(itertools.combinations(list_cell_types, 3)) + list(itertools.combinations(list_cell_types, 2))

        for (method1, method2) in [("2D", "3D"), ("3D", "HS"), ("2D", "HS")]:

            for combi in combi_cell_types: 
                
                title = "-".join(sorted(list(combi)))

                masks.update({"increased_"+method2.replace(" ","_")+"_versus_"+method1.replace(" ","_")+"_"+title: np.product(np.array([table['test_diff_'+test+"_"+comparison]  for cell_type in combi for comparison in [method2+"-"+method1+"_"+cell_type]])==1, axis=0).astype(bool)})
                masks.update({"decreased_"+method2.replace(" ","_")+"_versus_"+method1.replace(" ","_")+"_"+title: np.product(np.array([table['test_diff_'+test+"_"+comparison]  for cell_type in combi for comparison in [method2+"-"+method1+"_"+cell_type]])==-1, axis=0).astype(bool)})
                           

    names = list(masks.keys())
    
    index_ = [list(np.where(mask)[0]) for k,mask in masks.items()]    
        
    return masks, names, index_, cols




def get_masks_on_off(table, list_cell_types):
        
    if len(list_cell_types)==1:
        
        cell_type = list_cell_types[0]
            
        if "MS" in protocols[cell_type]:
            
            masks = OrderedDict({"all_valid_+": get_mask_valid(table, conds_equal={"HS":3,"MS":3,"3D":3,"2D":3}),
                     "all_2_valid_+": get_mask_valid(table, conds_equal={"HS":2,"MS":2,"3D":2,"2D":2}),
                     "all_at_least_1_valid_+": get_mask_valid(table, conds_sup={"HS":1,"MS":1,"3D":1,"2D":1}),
                     
                     "valid_HS_3_MS_3_not_valid_2D_3D_+": get_mask_valid(table, conds_equal={"HS":3,"MS":3,"3D":0,"2D":0}),
                     "valid_HS_3_not_valid_2D_3D_+": get_mask_valid(table, conds_equal={"HS":3,"3D":0,"2D":0}),
                     "valid_HS_3_MS_3_not_valid_3D_+": get_mask_valid(table, conds_equal={"HS":3,"MS":3,"3D":0}),
                     "valid_HS_3_not_valid_3D_+": get_mask_valid(table, conds_equal={"HS":3,"3D":0}),
        
        
                    "at_least_one_valid_HS_MS_not_valid_2D_3D_+": get_mask_valid(table, conds_equal={"2D":0,"3D":0}, conds_sup={"HS":1,"MS":1}),
                    "at_least_one_valid_HS_not_valid_2D_3D_+": get_mask_valid(table, conds_equal={"2D":0,"3D":0}, conds_sup={"HS":1}),
                    "at_least_one_valid_HS_MS_not_valid_3D_+": get_mask_valid(table, conds_equal={"3D":0}, conds_sup={"HS":1,"MS":1}),
                    "at_least_one_valid_HS_not_valid_3D_+": get_mask_valid(table, conds_equal={"3D":0}, conds_sup={"HS":1}),
        
        
                     "valid_2D_3_3D_3_not_valid_HS_MS_+": get_mask_valid(table, conds_equal={"HS":0,"MS":0,"3D":3,"2D":3}),
                     "valid_3D_3_not_valid_HS_MS_+": get_mask_valid(table, conds_equal={"MS":0,"HS":0,"3D":3}),
                     "valid_3D_3_not_valid_HS_+": get_mask_valid(table, conds_equal={"HS":0,"3D":3}),
        
                    "at_least_one_not_valid_2D_3D_not_valid_HS_MS_+": get_mask_valid(table, conds_sup={"2D":1, "3D":1}, conds_equal={"HS":0,"MS":0}),
                    "at_least_one_not_valid_3D_not_valid_HS_MS_+": get_mask_valid(table, conds_sup={"3D":1}, conds_equal={"HS":0,"MS":0}),
        
                    "at_least_one_not_valid_3D_not_valid_HS_+": get_mask_valid(table, conds_sup={"3D":1}, conds_equal={"HS":0}),
        
        
            })
        
        
        else:
            
        
            masks = OrderedDict({"all_valid_+": get_mask_valid(table, conds_equal={"HS":3,"3D":3,"2D":3}),
                     "all_2_valid_+": get_mask_valid(table, conds_equal={"HS":2,"3D":2,"2D":2}),
                     "all_at_least_1_valid_+": get_mask_valid(table, conds_sup={"HS":1,"3D":1,"2D":1}),
                     
                     "valid_HS_3_not_valid_2D_3D_+": get_mask_valid(table, conds_equal={"HS":3,"3D":0,"2D":0}),
                     "valid_HS_3_not_valid_3D_+": get_mask_valid(table, conds_equal={"HS":3,"3D":0}),
        
        
                    "at_least_one_valid_HS_not_valid_2D_3D_+": get_mask_valid(table, conds_equal={"2D":0,"3D":0}, conds_sup={"HS":1}),
                    "at_least_one_valid_HS_not_valid_3D_+": get_mask_valid(table, conds_equal={"3D":0}, conds_sup={"HS":1}),
        
        
                     "valid_2D_3_3D_3_not_valid_HS_+": get_mask_valid(table, conds_equal={"HS":0,"3D":3,"2D":3}),
                     "valid_3D_3_not_valid_HS_+": get_mask_valid(table, conds_equal={"HS":0,"3D":3}),
        
                    "at_least_one_not_valid_2D_3D_not_valid_HS_+": get_mask_valid(table, conds_sup={"2D":1, "3D":1}, conds_equal={"HS":0}),
        
                    "at_least_one_not_valid_3D_not_valid_HS_+": get_mask_valid(table, conds_sup={"3D":1}, conds_equal={"HS":0}),
        
        
            }) 
            
            
    else:
        
        
        masks = OrderedDict({
            
                 # "at_least_one_valid_THP1_MSCH_Fibro+": get_mask_valid(table, conds_sup={"HS_THP1":3,"3D_THP1":3,"2D_THP1":3,
                 #                                                              "HS_MSCH":3,"3D_MSCH":3,"2D_MSCH":3}),
                
                 "HS_on_3D_off_THP1_MSCH+": get_mask_valid(table, conds_equal={"3D_THP1":0, "3D_MSCH":0},conds_sup={"HS_THP1":1,"HS_MSCH":1}),
                 "HS_off_3D_on_THP1_MSCH-": get_mask_valid(table, conds_equal={"HS_THP1":0, "HS_MSCH":0},conds_sup={"3D_THP1":1,"3D_MSCH":1}),

                 "HS_on_3D_off_THP1_Fibro+": get_mask_valid(table, conds_equal={"3D_THP1":0, "3D_Fibro":0},conds_sup={"HS_THP1":1,"HS_Fibro":1}),
                 "HS_off_3D_on_THP1_Fibro-": get_mask_valid(table, conds_equal={"HS_THP1":0, "HS_Fibro":0},conds_sup={"3D_THP1":1,"3D_Fibro":1}),
 
                 "HS_on_3D_off_MSCH_Fibro+": get_mask_valid(table, conds_equal={"3D_MSCH":0, "3D_Fibro":0},conds_sup={"HS_MSCH":1,"HS_Fibro":1}),
                 "HS_off_3D_on_MSCH_Fibro-": get_mask_valid(table, conds_equal={"HS_MSCH":0, "HS_Fibro":0},conds_sup={"3D_MSCH":1,"3D_Fibro":1}),
                               
                 "HS_on_3D_off_THP1_MSCH_Fibro+": get_mask_valid(table, conds_equal={"3D_THP1":0, "3D_MSCH":0, "3D_Fibro":0},conds_sup={"HS_THP1":1,"HS_MSCH":1, "HS_Fibro":1}),
                 "HS_off_3D_on_THP1_MSCH_Fibro-": get_mask_valid(table, conds_equal={"HS_THP1":0, "HS_MSCH":0, "HS_Fibro":0},conds_sup={"3D_THP1":1,"3D_MSCH":1, "3D_Fibro":1}),


                 "HS_on_2D_off_THP1_MSCH+": get_mask_valid(table, conds_equal={"2D_THP1":0, "2D_MSCH":0},conds_sup={"HS_THP1":1,"HS_MSCH":1}),
                 "HS_off_2D_on_THP1_MSCH-": get_mask_valid(table, conds_equal={"HS_THP1":0, "HS_MSCH":0},conds_sup={"2D_THP1":1,"2D_MSCH":1}),

                 "HS_on_2D_off_THP1_Fibro+": get_mask_valid(table, conds_equal={"2D_THP1":0, "2D_Fibro":0},conds_sup={"HS_THP1":1,"HS_Fibro":1}),
                 "HS_off_2D_on_THP1_Fibro-": get_mask_valid(table, conds_equal={"HS_THP1":0, "HS_Fibro":0},conds_sup={"2D_THP1":1,"2D_Fibro":1}),
 
                 "HS_on_2D_off_MSCH_Fibro+": get_mask_valid(table, conds_equal={"2D_MSCH":0, "2D_Fibro":0},conds_sup={"HS_MSCH":1,"HS_Fibro":1}),
                 "HS_off_2D_on_MSCH_Fibro-": get_mask_valid(table, conds_equal={"HS_MSCH":0, "HS_Fibro":0},conds_sup={"2D_MSCH":1,"2D_Fibro":1}),
                               
                 "HS_on_2D_off_THP1_MSCH_Fibro+": get_mask_valid(table, conds_equal={"2D_THP1":0, "2D_MSCH":0, "2D_Fibro":0},conds_sup={"HS_THP1":1,"HS_MSCH":1, "HS_Fibro":1}),
                 "HS_off_2D_on_THP1_MSCH_Fibro-": get_mask_valid(table, conds_equal={"HS_THP1":0, "HS_MSCH":0, "HS_Fibro":0},conds_sup={"2D_THP1":1,"2D_MSCH":1, "2D_Fibro":1})


        
        })
    


    names = list(masks.keys())
    
    index_ = [list(np.where(mask)[0]) for k,mask in masks.items()]    
        
    return masks, names, index_





def get_masks_test_diff_thp1_for_heatmap(table, list_cell_types):


    masks = OrderedDict({
            
#            "at_least_one_value_and_similar": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-3D","HS-2D","MS-2D","MS-3D","MS]])==0, axis=0).astype(bool) & 

            "increased_MS.6h_versus_2D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["MS.6h-2D"]])==1, axis=0).astype(bool),
            "decreased_MS.6h_versus_2D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["MS.6h-2D"]])==-1, axis=0).astype(bool),

            
            # "increased_HS_MS_versus_2D_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-3D", "HS-2D", "MS-2D", "MS-3D"]])==1, axis=0).astype(bool),
            # "increased_HS_versus_2D_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-3D", "HS-2D"]])==1, axis=0).astype(bool),
            # "increased_HS_versus_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-3D"]])==1, axis=0).astype(bool),
            # "increased_HS_versus_2D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-2D"]])==1, axis=0).astype(bool),
            # "increased_MS_versus_2D_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["MS-3D", "MS-2D"]])==1, axis=0).astype(bool),
            # "increased_MS_versus_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["MS-3D"]])==1, axis=0).astype(bool),

            # "decreased_HS_MS_versus_2D_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-3D", "HS-2D", "MS-2D", "MS-3D"]])==-1, axis=0).astype(bool),
            # "decreased_HS_versus_2D_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-3D", "HS-2D"]])==-1, axis=0).astype(bool),
            # "decreased_HS_versus_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-3D"]])==-1, axis=0).astype(bool),
            # "decreased_HS_versus_2D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-2D"]])==-1, axis=0).astype(bool),
            # "decreased_MS_versus_2D_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["MS-3D", "MS-2D"]])==-1, axis=0).astype(bool),
            # "decreased_MS_versus_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["MS-3D"]])==-1, axis=0).astype(bool),

            # "no_change": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-3D", "HS-2D", "MS-2D", "MS-3D", "3D-2D", "HS-MS"]])==0, axis=0).astype(bool),

            # "no_change_HS_MS_versus_2D_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-3D", "HS-2D", "MS-2D", "MS-3D"]])==0, axis=0).astype(bool),
            # "no_change_HS_versus_2D_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-3D", "HS-2D"]])==0, axis=0).astype(bool),
            # "no_change_HS_versus_3D": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-3D"]])==0, axis=0).astype(bool),

            })

    names = list(masks.keys())
    
    index_ = [list(np.where(mask)[0]) for k,mask in masks.items()]    
        
    return masks, names, index_




def get_masks_test_diff(table, list_cell_types):
    
    if len(list_cell_types)==1:
        
        cell_type = list_cell_types[0]
            
        masks = OrderedDict({})
        
        for (method1, method2) in interesting_comparisons[cell_type]:
            
            masks["increased_"+method2+"vs"+method1] = (table["test_diff_"+method2+"-"+method1] == 1)
            masks["decreased_"+method2+"vs"+method1] = (table["test_diff_"+method2+"-"+method1] == -1)
            
        return masks
        


def get_masks_test_diff_for_heatmap(table, list_cell_types, test="ProDA"):

    
    if len(list_cell_types)==1:
        
        cell_type = list_cell_types[0]
            
        if "MS" in protocols[cell_type]:
        
            masks = OrderedDict({
                    
        #            "at_least_one_value_and_similar": np.product(np.array([table['test_diff_'+comparison] for comparison in ["HS-3D","HS-2D","MS-2D","MS-3D","MS]])==0, axis=0).astype(bool) & 
                    
                    "increased_HS_MS_versus_2D_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D", "HS-2D", "MS-2D", "MS-3D"]])==1, axis=0).astype(bool),
                    "increased_HS_versus_2D_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D", "HS-2D"]])==1, axis=0).astype(bool),
                    "increased_HS_versus_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D"]])==1, axis=0).astype(bool),
                    "increased_HS_versus_2D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-2D"]])==1, axis=0).astype(bool),
                    "increased_MS_versus_2D_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["MS-3D", "MS-2D"]])==1, axis=0).astype(bool),
                    "increased_MS_versus_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["MS-3D"]])==1, axis=0).astype(bool),
        
                    "decreased_HS_MS_versus_2D_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D", "HS-2D", "MS-2D", "MS-3D"]])==-1, axis=0).astype(bool),
                    "decreased_HS_versus_2D_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D", "HS-2D"]])==-1, axis=0).astype(bool),
                    "decreased_HS_versus_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D"]])==-1, axis=0).astype(bool),
                    "decreased_HS_versus_2D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-2D"]])==-1, axis=0).astype(bool),
                    "decreased_MS_versus_2D_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["MS-3D", "MS-2D"]])==-1, axis=0).astype(bool),
                    "decreased_MS_versus_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["MS-3D"]])==-1, axis=0).astype(bool),
        
                    "no_change": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D", "HS-2D", "MS-2D", "MS-3D", "3D-2D", "HS-MS"]])==0, axis=0).astype(bool),
        
                    "no_change_HS_MS_versus_2D_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D", "HS-2D", "MS-2D", "MS-3D"]])==0, axis=0).astype(bool),
                    "no_change_HS_versus_2D_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D", "HS-2D"]])==0, axis=0).astype(bool),
                    "no_change_HS_versus_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D"]])==0, axis=0).astype(bool),
        
                    })
         

        else:
        
            masks = OrderedDict({
                    "increased_HS_versus_2D_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D", "HS-2D"]])==1, axis=0).astype(bool),
                    "increased_HS_versus_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D"]])==1, axis=0).astype(bool),
        
                    "decreased_HS_versus_2D_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D", "HS-2D"]])==-1, axis=0).astype(bool),
                    "decreased_HS_versus_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D"]])==-1, axis=0).astype(bool),
        
                    "no_change": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D", "HS-2D", "3D-2D"]])==0, axis=0).astype(bool),
        
                    "no_change_HS_versus_2D_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D", "HS-2D"]])==0, axis=0).astype(bool),
                    "no_change_HS_versus_3D": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D"]])==0, axis=0).astype(bool),
        
                    })
        


    else:
        
        
        
        masks = OrderedDict({
        
                "increased_HS_versus_2D_3D_MSCH_THP1_Fibro": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_Fibro", "HS-2D_Fibro", "HS-3D_THP1", "HS-2D_THP1","HS-3D_MSCH", "HS-2D_MSCH"]])==1, axis=0).astype(bool),
                
                "increased_HS_versus_2D_3D_MSCH_THP1": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_THP1", "HS-2D_THP1","HS-3D_MSCH", "HS-2D_MSCH"]])==1, axis=0).astype(bool),
                "increased_HS_versus_2D_3D_MSCH_Fibro": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_Fibro", "HS-2D_Fibro","HS-3D_MSCH", "HS-2D_MSCH"]])==1, axis=0).astype(bool),
                "increased_HS_versus_2D_3D_THP1_Fibro": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_THP1", "HS-2D_THP1","HS-3D_Fibro", "HS-2D_Fibro"]])==1, axis=0).astype(bool),
        
                "increased_HS_versus_3D_MSCH_THP1_Fibro": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_Fibro", "HS-3D_THP1", "HS-3D_MSCH"]])==1, axis=0).astype(bool),
        
                "increased_HS_versus_3D_MSCH_THP1": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_THP1", "HS-3D_MSCH"]])==1, axis=0).astype(bool),
                "increased_HS_versus_3D_MSCH_Fibro": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_Fibro", "HS-3D_MSCH"]])==1, axis=0).astype(bool),
                "increased_HS_versus_3D_THP1_Fibro": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_THP1", "HS-3D_Fibro"]])==1, axis=0).astype(bool),
        
                "decreased_HS_versus_2D_3D_MSCH_THP1_Fibro": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_Fibro", "HS-2D_Fibro", "HS-3D_THP1", "HS-2D_THP1","HS-3D_MSCH", "HS-2D_MSCH"]])==-1, axis=0).astype(bool),
        
                "decreased_HS_versus_2D_3D_MSCH_THP1": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_THP1", "HS-2D_THP1","HS-3D_MSCH", "HS-2D_MSCH"]])==-1, axis=0).astype(bool),
                "decreased_HS_versus_2D_3D_MSCH_Fibro": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_Fibro", "HS-2D_Fibro","HS-3D_MSCH", "HS-2D_MSCH"]])==-1, axis=0).astype(bool),
                "decreased_HS_versus_2D_3D_THP1_Fibro": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_THP1", "HS-2D_THP1","HS-3D_Fibro", "HS-2D_Fibro"]])==-1, axis=0).astype(bool),
        
                "decreased_HS_versus_3D_MSCH_THP1_Fibro": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_Fibro", "HS-3D_THP1", "HS-3D_MSCH"]])==-1, axis=0).astype(bool),
        
                "decreased_HS_versus_3D_MSCH_THP1": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_THP1", "HS-3D_MSCH"]])==-1, axis=0).astype(bool),
                "decreased_HS_versus_3D_MSCH_Fibro": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_Fibro", "HS-3D_MSCH"]])==-1, axis=0).astype(bool),
                "decreased_HS_versus_3D_THP1_Fibro": np.product(np.array([table['test_diff_'+test+"_"+comparison] for comparison in ["HS-3D_THP1", "HS-3D_Fibro"]])==-1, axis=0).astype(bool),
        
        
        
        
        
                })
 

    names = list(masks.keys())
    
    index_ = [list(np.where(mask)[0]) for k,mask in masks.items()]    
        
    return masks, names, index_












        


