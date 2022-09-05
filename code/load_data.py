import os
from pathlib import Path
import os
from paths import datapath, savepath
import pandas
import numpy as np
import pickle
from data_infos import protocols, protocols_lysat
from data_infos import experiments, cell_types, table_columns, protocols_lysat, experiments_lysat
from data_infos import dic_protocol_corresp, dic_species
import itertools
from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection, fdrcorrection_twostage
from itertools import chain
from collections import Counter
from scipy.stats import fisher_exact
from scipy.stats import hypergeom
from constants import quantile
from utils.bioinf_conversion import convert_entrez_gene_id_list_to_gene_name
from scipy.stats import ttest_ind
import scipy
import matplotlib.pyplot as plt
from background_genomes import load_background_genome_table
pandas.set_option('display.max_columns', None)

list_comparisons = [("2D", "3D"), ("2D", "HS"), ("3D", "HS"), ("2D", "MS"), ("3D", "MS"), ("MS", "HS"), ("2D", "MS.6h"), ("3D", "MS.6h"), ("MS", "MS.6h"), ("MS.6h", "HS")]

list_comparisons += [("Lysat "+protocol1, "Lysat "+protocol2) for (protocol1, protocol2) in list_comparisons]

list_comparisons += [(protocol, "Lysat "+protocol) for protocol in ["2D", "3D", "MS", "MS.6h", "HS"]]



list_categories = ["GOTERM_CC_DIRECT", "GOTERM_BP_DIRECT"]


excel_files = os.listdir(datapath)


dic_files = {"THP1": 'proteinGroups - THP1.txt',
             "MSCH": 'proteinGroups_MSC H EV CL.txt',
             "MSCS": 'proteinGroups_MSC S.txt',
             "MSCSmapH": 'proteinGroups_MSC S.txt',
             "Fibro": 'proteinGroups Fibroblast.txt'
             }

human_cell_types = ["THP1", "MSCH", "Fibro"]#, "MSCSmapH"]





class ProteoDataset():
    
    
    def __init__(self, prepare_data=False, reset=False):
        
        
        self.complete_human_df = self._get_complete_human_df()
        

        if prepare_data:
            

           self.save_data_for_R()

       
         
        for cell_type in cell_types:
            
            table = self.load_table(cell_type, reset=reset)

        

    
    def load_raw_protein_group_table(self, cell_type):
        
        filename = dic_files[cell_type]
        table = []
        
        with open(Path(datapath, "raw", filename), "r") as file:
            lines = file.readlines()
            
        for line in lines:
            line = line.replace("\n","").split("\t")
            table.append(line)
        
        # Take majority protein ID
        table = pandas.DataFrame(table[1:], columns=table[0])

        return table


    def _preprocess_table(self, table, cell_type, log2=True):

        table["UNIPROT_ACCESSION"] = table["Majority protein IDs"].apply(lambda s: s.split(";")[0])
        table["GENE_NAME1"] = table["Gene names"].apply(lambda s: s.split(";")[0])
        table["PROTEIN_NAME"] = table["Protein names"].apply(lambda s: s.split(";")[0])

        mask = (table['Only identified by site']!="+") \
                & (table['Reverse']!="+") \
                & (table['Potential contaminant']!="+")
                        
        table = table[mask]
        table.reset_index(inplace=True, drop=True)
        
        david_table = load_background_genome_table(dic_species[cell_type])
        
        david_table = david_table.reset_index(drop=True)

        david_table["UNIPROT_ACCESSION"] = david_table["UNIPROT_ACCESSION"].astype(str).apply(lambda x: x.split(','))
 
        david_table = david_table.explode("UNIPROT_ACCESSION")
        
        david_table.drop_duplicates(subset=["UNIPROT_ACCESSION"], inplace=True, keep="first")
                            
        table = pandas.merge(left=table, \
                              right=david_table.astype(str), how="left", \
                              on=["UNIPROT_ACCESSION"])


        table = table.reset_index(drop=True)
        
                 
        mask = table["ENTREZ_GENE_ID"].isna()
        
        ### Remove proteins for which gene ID has not been found
        
        table = table[~mask]
        
        table = table.reset_index(drop=True)
        
                
        # print(cell_type, mask.sum(), "removed")
        
        # print(table[["GENE_NAME1", "GENE_NAME"]]) Comparison David - excel raw files
        # At the end Keep only David
        


        table = table.rename(columns=dic_protocol_corresp[cell_type])

        LFQI_cols = [col for col in table.columns if "LFQI" in col]
        
        if cell_type!="MSCSmapH":
        
            table = table[["ENTREZ_GENE_ID", "UNIPROT_ACCESSION", "GENE_NAME", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT"] + LFQI_cols]
            
        else:
            
            
                    

            mouse_human_mapping = pandas.read_csv(Path(datapath, "human_mouse_mapping.txt"), sep="\t", low_memory=False)
    
            mouse_human_mapping.columns = [col.replace("human_entrez_gene","HUMAN_ENTREZ_GENE_ID").replace("mouse_entrez_gene", "ENTREZ_GENE_ID") for col in mouse_human_mapping.columns]
    
            table = table.merge(mouse_human_mapping[["ENTREZ_GENE_ID", "HUMAN_ENTREZ_GENE_ID"]], on="ENTREZ_GENE_ID", how="left")
            
            human_genome = load_background_genome_table("homo_sapiens")
            
            human_genome["HUMAN_ENTREZ_GENE_ID"] = human_genome["ENTREZ_GENE_ID"]
            human_genome["HUMAN_GENE_NAME"] = human_genome["GENE_NAME"]
            human_genome["HUMAN_GOTERM_BP_DIRECT"] = human_genome["GOTERM_BP_DIRECT"]
            human_genome["HUMAN_GOTERM_CC_DIRECT"] = human_genome["GOTERM_CC_DIRECT"]

            
            table = pandas.merge(table, human_genome[["HUMAN_ENTREZ_GENE_ID", "HUMAN_GENE_NAME", "HUMAN_GOTERM_BP_DIRECT", "HUMAN_GOTERM_CC_DIRECT"]], on="HUMAN_ENTREZ_GENE_ID", how="left")

            table = table[["HUMAN_ENTREZ_GENE_ID", "HUMAN_GENE_NAME", "HUMAN_GOTERM_BP_DIRECT", "HUMAN_GOTERM_CC_DIRECT", "UNIPROT_ACCESSION", "ENTREZ_GENE_ID", "GENE_NAME"] + LFQI_cols]
            table = table.rename(columns={"ENTREZ_GENE_ID": "MOUSE_ENTREZ_GENE_ID",
                                          "GENE_NAME": "MOUSE_GENE_NAME",
                                          "HUMAN_ENTREZ_GENE_ID": "ENTREZ_GENE_ID",
                                          "HUMAN_GENE_NAME": "GENE_NAME",
                                          "HUMAN_GOTERM_BP_DIRECT": "GOTERM_BP_DIRECT",
                                          "HUMAN_GOTERM_CC_DIRECT": "GOTERM_CC_DIRECT"})
            
    
        table[LFQI_cols] = table[LFQI_cols].astype(float)
    
        # Change 0 to missing values
        # MaxQuant codes missing values as 0. The actual abundance was probably not zero, 
        # but just some value too small to be detected by the mass spectrometer
        table[LFQI_cols] = table[LFQI_cols].applymap(lambda s: np.nan if s==0 else s)


        if log2:
            # Apply log2 transformation
            table[LFQI_cols] = table[LFQI_cols].applymap(lambda s: np.log2(s) \
                                                         if np.isnan(s)==False else s)
                

        table = table.drop_duplicates(subset=["ENTREZ_GENE_ID"])
        
        table[["ENTREZ_GENE_ID", "GENE_NAME", "UNIPROT_ACCESSION"]] = table[["ENTREZ_GENE_ID", "GENE_NAME", "UNIPROT_ACCESSION"]].astype(str)

        for cat in list_categories:
            table[cat] = table[cat].fillna("").astype(str).str.replace(", "," - ")

        mask_valid = ((~table[[col for col in table.columns if "LFQ" in col]].isna()).sum(axis=1)>0)
        table = table[mask_valid]     


        return table



    def _get_complete_human_df(self):
        
        list_id = []

        for i, cell_type in enumerate(human_cell_types+["MSCSmapH"]):

            raw_table = self.load_raw_protein_group_table(cell_type)
            preprocessed_table = self._preprocess_table(raw_table, cell_type)
                         
            if i==0:
                concatenated_table = preprocessed_table[["ENTREZ_GENE_ID", "GENE_NAME", "UNIPROT_ACCESSION"]]
                
            else:
                concatenated_table = pandas.concat([concatenated_table, preprocessed_table[["ENTREZ_GENE_ID", "GENE_NAME", "UNIPROT_ACCESSION"]]], axis=0)

        concatenated_table = concatenated_table.drop_duplicates(subset=["ENTREZ_GENE_ID"])

        return concatenated_table



    def _get_preprocessed_table(self, cell_type, human_mapping=False, log2=True):   

        raw_table = self.load_raw_protein_group_table(cell_type)
        
        preprocessed_table = self._preprocess_table(raw_table, cell_type, log2=log2)
        
        preprocessed_table = preprocessed_table.drop_duplicates(subset=["ENTREZ_GENE_ID"])
        
        if cell_type=="MSCS":

            return preprocessed_table

        
        else:
            
            unique_genes_in_table = preprocessed_table["ENTREZ_GENE_ID"].unique()
            absent_genes = [ID for ID in self.complete_human_df["ENTREZ_GENE_ID"] if ID not in unique_genes_in_table]
            

            mask_absent = self.complete_human_df["ENTREZ_GENE_ID"].apply(lambda s: s in absent_genes)

            preprocessed_table = pandas.concat([preprocessed_table, self.complete_human_df[mask_absent]], axis=0)
            
        preprocessed_table.reset_index(inplace=True, drop=True)   
        
        return preprocessed_table



    

        
    def save_data_for_R(self):

        for cell_type in cell_types:
            
            preprocessed_table = self._get_preprocessed_table(cell_type)

            rows = experiments[cell_type] + experiments_lysat[cell_type]
            columns = protocols[cell_type] + protocols_lysat[cell_type]
    
            table_for_analysis = preprocessed_table[rows]        
            design_matrix = np.zeros((len(rows),len(columns)))
                    
            for protocol in range(len(columns)):
                design_matrix[protocol*3:protocol*3+3,protocol] = 1
                
            design_matrix = pandas.DataFrame(design_matrix, columns=columns, index=rows)

            table_for_analysis.to_csv(Path(savepath,"for_R", cell_type+"_table_for_analysis.csv"), index=False)
                            
            design_matrix.to_csv(Path(savepath,"for_R", cell_type+"_design_matrix_for_analysis.csv"), index=False)
            
            preprocessed_table = self._get_preprocessed_table(cell_type, log2=False)
            
            table_for_analysis = pandas.DataFrame(preprocessed_table[["ENTREZ_GENE_ID", "GENE_NAME"]+rows])

            table_for_analysis.columns = [col.replace("LFQI ","LFQ.intensity.") for col in table_for_analysis.columns]
            
            # table_for_analysis = self.add_classical_estimations_imputation(table_for_analysis, cell_type)

            # imputed_cols = [col for col in table_for_analysis.columns if "imputed" in col]
                        
            # table_for_analysis = table_for_analysis[["ID"]+imputed_cols]
                        
            # table_for_analysis.columns = [col.replace("_imputed","").replace(" 1","_1").replace(" 2","_2").replace(" 3", "_3").replace("Lysat_","Lysat.").replace("2D","X2D").replace("3D","X3D").replace("ID","Protein.IDs") for col in table_for_analysis.columns]

            table_for_analysis.columns = [col.replace(" 1","_1").replace(" 2","_2").replace(" 3", "_3").replace("Lysat_","Lysat.").replace("2D","X2D").replace("3D","X3D").replace("ENTREZ_GENE_ID","Gene.IDs").replace("GENE_NAME","Gene.name") for col in table_for_analysis.columns]

            table_for_analysis.to_csv(Path(savepath,"for_R", cell_type+"_table_for_analysis_DEP.csv"), index=False)

            labels = [col.replace("LFQ.intensity.","") for col in table_for_analysis.columns if "ID" not in col and "name" not in col]
                        
            conditions = [label.split("_")[0] for label in labels]
            replicate = [label.split("_")[1] for label in labels]

            design_matrix = pandas.DataFrame(np.array([labels, conditions, replicate]).T, columns = ["label", "condition", "replicate"])

            # design_matrix = design_matrix.applymap(lambda s: s.replace("2D", "X2D").replace("3D","X3D"))

            design_matrix.to_csv(Path(savepath,"for_R", cell_type+"_design_matrix_for_analysis_DEP.csv"), index=False)


        

#

    def load_table(self, cell_type, load_lysat=False, reset=False):
        
        table = self._load_complete_table(cell_type, reset=reset)

        if not load_lysat:
            table = table[[col for col in table.columns if "Lysat" not in col]]

        mask_valid = (table[[col for col in table.columns if "nb_valid" in col]].sum(axis=1)>0)
        table = table[mask_valid]

        filepath = Path(savepath, "table_"+cell_type+".csv")
        
        table.to_csv(str(filepath))

        
        return table


    def _load_complete_table(self, cell_type, reset=False):

        filepath = Path(savepath, "table_"+cell_type+".pkl")
     
        if reset or not os.path.exists(filepath):

            table = self._get_preprocessed_table(cell_type)

            """ Add R results """
            
            R_results = self.get_R_results(cell_type)
            
            # DEP imputation results
            DEP_imputation_table = R_results["DEP_imputation"]
            DEP_imputation_table["name"] = DEP_imputation_table.index
            DEP_imputation_table.reset_index(inplace=True, drop=True)
            
            DEP_imputation_metadata_table = R_results["DEP_imputation_metadata"]
            DEP_imputation_table["name"] = DEP_imputation_table["name"].astype(str)

            DEP_imputation_table = pandas.merge(DEP_imputation_table, DEP_imputation_metadata_table[["ID", "name"]].astype(str), on="name", how="inner")

            DEP_imputation_table.rename(columns={"ID": "ENTREZ_GENE_ID"}, inplace=True)
            DEP_imputation_table.drop("name", axis=1, inplace=True)
 
            DEP_imputation_table.columns = ["LFQI "+col.replace("_"," ")+"_DEP_imputation" if col!="ENTREZ_GENE_ID" else col for col in DEP_imputation_table.columns ]
            DEP_imputation_table[[col for col in DEP_imputation_table.columns if "imputation" in col]] = DEP_imputation_table[[col for col in DEP_imputation_table.columns if "imputation" in col]].astype(float)
            
            table = pandas.merge(table, DEP_imputation_table, on="ENTREZ_GENE_ID", how="left")
             
            # DEP test diff results
            test_diff_DEP_table = pandas.DataFrame(R_results["test_diff_DEP"])
            test_diff_DEP_table = test_diff_DEP_table[["ID"]+[col for col in test_diff_DEP_table.columns if "p.val" in col or "ratio" in col]]
            test_diff_DEP_table.columns = [col.replace("ID","ENTREZ_GENE_ID").replace("Lysat.","Lysat ") for col in test_diff_DEP_table.columns]

            test_diff_DEP_table["ENTREZ_GENE_ID"] = test_diff_DEP_table["ENTREZ_GENE_ID"].astype(str)

            for (protocol1, protocol2) in list_comparisons:
                
                if protocol1 not in protocols[cell_type]+protocols_lysat[cell_type] or protocol2 not in protocols[cell_type]+protocols_lysat[cell_type]:
                    continue
                
                key1 = protocol1+"_vs_"+protocol2
                key2 = protocol2+"_vs_"+protocol1
                                
                if key1+"_p.val" in test_diff_DEP_table.columns:
                    
                    test_diff_DEP_table["pval_DEP_"+protocol2+"-"+protocol1] = test_diff_DEP_table[key1+"_p.val"]
                    test_diff_DEP_table["FC_DEP_"+protocol2+"-"+protocol1] = 1/(2**test_diff_DEP_table[key1+"_ratio"])
                    
                elif key2+"_p.val" in test_diff_DEP_table.columns:
                    test_diff_DEP_table["pval_DEP_"+protocol2+"-"+protocol1] = test_diff_DEP_table[key2+"_p.val"]
                    test_diff_DEP_table["FC_DEP_"+protocol2+"-"+protocol1] = 2**test_diff_DEP_table[key2+"_ratio"]

                FC = test_diff_DEP_table["FC_DEP_"+protocol2+"-"+protocol1]
                pvalues = test_diff_DEP_table["pval_DEP_"+protocol2+"-"+protocol1]
    
                double_filter = ((np.array(pvalues)<0.01) & ((FC>=1.5)|(FC<=(1/1.5))))
                
                test_diff = double_filter.astype(int).copy()
                test_diff[(FC<1) & (test_diff==1)] = -1
                
                test_diff_DEP_table["test_diff_DEP_"+protocol2+"-"+protocol1] = test_diff       
   
            table = table.merge(test_diff_DEP_table[["ENTREZ_GENE_ID"]+[col for col in test_diff_DEP_table.columns if "DEP" in col]], on="ENTREZ_GENE_ID", how="left")
    
            for protocol in protocols[cell_type] + protocols_lysat[cell_type]:
                
                LFQ_imputed_cols = [col for col in table.columns if "imputation" in col and protocol in col]
                
                
                table["average_imputation_DEP_"+protocol] = table[LFQ_imputed_cols].mean(axis=1)
    
            # mask = ~(table["LFQI HS 1_DEP_imputation"].isna())
            # print(table[mask][[col for col in table.columns if ("LFQI HS" in col or "LFQI 3D" in col) or ("FC" in col and "HS" in col and "3D" in col)]])
            # retrieve FC DEP
            # cols_HS = ["LFQI HS "+str(k)+"_DEP_imputation" for k in range(1,4)]
            # cols_3D = ["LFQI 3D "+str(k)+"_DEP_imputation" for k in range(1,4)]
            # average = (2**table[mask][cols_HS].mean(axis=1))/(2**table[mask][cols_3D].mean(axis=1))

            
            # Add ProDA results

            for key in ["coeffs"]:           
                for col in R_results[key].columns:
                    table[key+"_"+col] = np.array(R_results[key][col])  

            coeffs_covariance_matrix = R_results["coeffs_covariance"]
                    
            pval_dic = R_results["pval_dic"]
            
            nb = len(protocols[cell_type]+protocols_lysat[cell_type])
            
            for (protocol1, protocol2) in list_comparisons:
                
                if protocol1 not in protocols[cell_type]+protocols_lysat[cell_type] or protocol2 not in protocols[cell_type]+protocols_lysat[cell_type]:
                    continue    
                
                pvalues = []

                key = protocol2+"-"+protocol1

                for i in range(len(table)):
                    
                    diff = table.iloc[i]["coeffs_"+protocol2] - table.iloc[i]["coeffs_"+protocol1]
                    covariance_mat = coeffs_covariance_matrix[i]
                    var_diff = covariance_mat.loc[protocol1, protocol1] + covariance_mat.loc[protocol2, protocol2] - 2*covariance_mat.loc[protocol1, protocol2]
                    wald_statistic = diff / np.sqrt(var_diff)
                    df =  nb*3 - nb #nb values - number params
                    calculated_pval = 2*scipy.stats.t(df=df).cdf(-np.abs(wald_statistic))
                    
                    pvalues.append(calculated_pval)   

                    # #In MSCS one covar matrix is NaN
                    # if key in pval_dic:
                    #     if np.round(pval_dic[key]["pval"].iloc[i], 6)!=np.round(calculated_pval, 6):
                    #         print("oups", np.round(pval_dic[key]["pval"].iloc[i], 6), np.round(calculated_pval, 6))

                            # print(protocol2, table.iloc[i]["coeffs_"+protocol2], protocol1, table.iloc[i]["coeffs_"+protocol1])
                            # print(cell_type)
                            # print(table.iloc[i][[cols for cols in table.columns if "LFQI" in cols]])
    
                # table["test_pval_"+key] = pval_dic[key].values
                
                FC = 2**table["coeffs_"+protocol2].values / 2**table["coeffs_"+protocol1].values
                
                # log2_FC = table["coeffs_"+protocol2].values - table["coeffs_"+protocol1].values
                
                double_filter = ((np.array(pvalues)<0.01) & ((FC>=1.5)|(FC<=(1/1.5))))
                
                
                
                test_diff = double_filter.astype(int).copy()
                test_diff[(FC<1) & (test_diff==1)] = -1
                
                test_diff = test_diff.astype(float)
                
                mask = np.isnan(pvalues)
                test_diff[mask] = np.nan
                
                
                
                
                # adjusted_prob = fdrcorrection_twostage(np.array(pval_dic[key].values), alpha=0.01, method="bh")[1]
                
                # table["test_diff_adjusted_"+key] = (adjusted_prob<0.01).astype(int)

                test_diff_table = pandas.DataFrame(np.array([pvalues, FC, test_diff]).T, columns=["pval_ProDA_"+key, "FC_ProDA_"+key, "test_diff_ProDA_"+key])

                table = pandas.concat([table, test_diff_table], axis=1)
                
                
                
                # table["log2_FC_"+key] = log2_FC

            # self.add_test_diff_results(table, cell_type)
             
            for protocol in protocols[cell_type] + protocols_lysat[cell_type]:
    
                columns = table_columns[cell_type][protocol]
                nb_valid = np.sum(~(table[columns].isna()),axis=1)
                table["nb_valid_"+protocol] = nb_valid
                

            
            
 
            for col in [col for col in table.columns if "nb_valid" in col]:
                table["all_valid_"+col.replace("nb_valid","")] = table[col].apply(lambda f: True if f==3 else False)
            table["filter"] = (table[[col for col in table.columns if "all_valid" in col]].sum(axis=1)>0)


            
            with open(filepath, "wb") as file:
                pickle.dump(table, file, protocol=pickle.HIGHEST_PROTOCOL)  


        with open(filepath, "rb") as file:
            table = pickle.load(file)      

                    
        return table
    


    def map_uniprot_to_gene_name(self, uniprot_id_list):
        
        uniprot_id_list = np.array(uniprot_id_list)
        
        genome = self.load_background_genome_table("homo_sapiens")

        gene_names_list = []

        genome["UNIPROT_ACCESSION"] = genome["UNIPROT_ACCESSION"].astype(str).apply(lambda str_: str_.split(","))
        
        for uni in uniprot_id_list:
            where_uniprot = genome["UNIPROT_ACCESSION"].apply(lambda list_: uni in list_)
            if np.sum(where_uniprot)==0:
                gene_name = ""
            else:
                gene_name = genome[where_uniprot]["Gene name"].values[0]

            gene_names_list.append(gene_name)
            
        return gene_names_list

    
    def load_human_table(self, reset=False, load_lysat=False):

        for i, cell_type in enumerate(human_cell_types+["MSCSmapH"]):
            
            table = self._load_complete_table(cell_type, reset=reset)
            
            if cell_type=="MSCSmapH":
                table.drop("MOUSE_ENTREZ_GENE_ID", axis=1, inplace=True)
                table.drop("MOUSE_GENE_NAME", axis=1, inplace=True)
                table.drop("UNIPROT_ACCESSION", axis=1, inplace=True)
            
            table.columns = [col+"_"+cell_type if col not in ["ENTREZ_GENE_ID"] else col for col in table.columns]

            if i==0:
                human_table = table 
            else:
                human_table = pandas.merge(human_table, table, on=["ENTREZ_GENE_ID"], how="outer")
                
        human_table["GOTERM_CC_DIRECT"] = human_table["GOTERM_CC_DIRECT_MSCH"]
        human_table["GOTERM_BP_DIRECT"] = human_table["GOTERM_BP_DIRECT_MSCH"]
        human_table["GENE_NAME"] = human_table["GENE_NAME_MSCH"]
                
        mask_valid = (human_table[[col for col in human_table.columns if "nb_valid" in col]].sum(axis=1)>0)
        
        human_table = human_table[mask_valid]
        
        if not load_lysat:
            human_table = human_table[[col for col in human_table.columns if "Lysat" not in col]]
                    
        return human_table
        



    def normalize_GAPDH(self, table):
   
        LFQI_cols = [col for col in table.columns if "LFQI" in col]
        expectation_cols = [col for col in table.columns if "coeffs" in col]
        mean_cols = [col for col in table.columns if "empirical_mean" in col]
        
        print(mean_cols)

        GAPDH_row = (table["Gene name"]=="GAPDH").copy()

        
        for col in expectation_cols:
            table["normalized_"+col] = table[col] - table[GAPDH_row][col].values[0]
        for col in mean_cols:
            table["normalized_"+col] = table[col] - table[GAPDH_row][col].values[0]
                        
        return table


    def discard_missing_values(self, table):
            

        mask = (table[[col for col in table.columns if "at_least_2_valids" in col]].product(axis=1)>0)
        table = table[mask]
        
        return table
    

        
    def get_R_results(self, cell_type):
   
        coeffs = pandas.read_csv(Path(savepath, "R_results", "coeffs_"+cell_type+".csv"))
        dropout_params = pandas.read_csv(Path(savepath, "R_results", "dropout_params_"+cell_type+".csv")) 
        coeffs_covariance = pandas.read_csv(Path(savepath, "R_results", "coef_var_"+cell_type+".csv"))

        coeffs.columns = [col.replace("X2D","2D").replace("X3D","3D").replace("Lysat.","Lysat ") for col in coeffs.columns]
        coeffs_covariance.columns = [col.replace("X2D","2D").replace("X3D","3D").replace("Lysat.", "Lysat ") for col in coeffs_covariance.columns]
        coeffs_covariance.index = [col.replace("X2D","2D").replace("X3D","3D").replace("Lysat.", "Lysat ") for col in coeffs_covariance.index]
        
        list_protocols = protocols[cell_type] + protocols_lysat[cell_type]
        coeffs_covariance = [coeffs_covariance.iloc[u*len(list_protocols):(u+1)*len(list_protocols)] for u in range(len(coeffs))]
        
        # Reproduce Wald test
        # # for (protocol1, protocol2) in list_comparisons:
        # for (protocol1, protocol2) in [("2D","HS")]:
        #     pval = pandas.read_csv(Path(savepath, "R_results", "test_"+protocol2+"_"+protocol1+"_"+cell_type+".csv"))
        #     pval.columns = ["pval"]
        #     t_stats = pandas.read_csv(Path(savepath, "R_results", "test_"+protocol2+"_"+protocol1+"_t_statistic"+"_"+cell_type+".csv"))
            
        #     for i in range(len(pval)):
                
        #         pvalue = pval.iloc[i].values[0]
                
        #         diff = coeffs.iloc[i][protocol2] - coeffs.iloc[i][protocol1]
        #         covariance_mat = coeffs_covariance[i]
        #         var_diff = covariance_mat.loc[protocol1, protocol1] + covariance_mat.loc[protocol2, protocol2] - 2*covariance_mat.loc[protocol1, protocol2]
        #         wald_statistic = diff / np.sqrt(var_diff)
        #         df =  6 #nb restrictions - number params
        #         calculated_pval = 2*scipy.stats.t(df=df).cdf(-np.abs(wald_statistic))
        #         t_stat = t_stats.iloc[i].values[0]       
                
        #         if np.round(pvalue,6)!=np.round(calculated_pval,6):
        #             print("oup")
                
        #         # if np.round(t_stat, 6) != np.round(wald_statistic, 6):
        #         #     print(t_stat, wald_statistic)
        #         #     print("wrong calculus")
                
        #         # print(pvalue, calculated_pval)
            
        #     print(coucou)

    
        pval_dic = {}
        for (protocol1, protocol2) in list_comparisons:
            
            if protocol1 not in protocols[cell_type]+protocols_lysat[cell_type] or protocol2 not in protocols[cell_type]+protocols_lysat[cell_type]:
                continue 
            
            if not os.path.exists(Path(savepath, "R_results", "test_"+protocol2+"_"+protocol1+"_"+cell_type+".csv")):
                continue
            pval = pandas.read_csv(Path(savepath, "R_results", "test_"+protocol2+"_"+protocol1+"_"+cell_type+".csv"))

            pval_dic[protocol2+"-"+protocol1] = pval


        df = pandas.read_csv(Path(savepath,"R_results", "df_"+cell_type+".csv"))
        s2 = pandas.read_csv(Path(savepath, "R_results", "s2_"+cell_type+".csv"))
        
        sigma2 = (s2*df)/(df+1)
        
        
        test_diff_DEP_results = pandas.read_csv(Path(savepath, "R_results", "test_diff_DEP_"+cell_type+".csv"), low_memory=False)
        imputation_DEP = pandas.read_csv(Path(savepath, "R_results", "DEP_imputation_"+cell_type+".csv"), low_memory=False)

        test_diff_DEP_results.columns = [col.replace("X2D","2D").replace("X3D","3D") for col in test_diff_DEP_results.columns]
        imputation_DEP.columns = [col.replace("X2D","2D").replace("X3D","3D") for col in imputation_DEP.columns]
        
        DEP_imputation_metadata = pandas.read_csv(Path(savepath, "R_results", "DEP_imputation_metadata_"+cell_type+".csv"), low_memory=False)
        
        
        results = {"coeffs":coeffs,
                   "dropout_params":dropout_params,
                   "coeffs_covariance":coeffs_covariance,
                   "variance":sigma2,
                    "pval_dic":pval_dic,
                    "test_diff_DEP": test_diff_DEP_results,
                    "DEP_imputation": imputation_DEP,
                    "DEP_imputation_metadata": DEP_imputation_metadata
                   }
        
        

    
        return results



    
    
    def get_uncertain_proteins(self, table, cell_type):
        
        for i, protocol in enumerate(protocols[cell_type]):
            
            nb_NA = table[table_columns[cell_type][protocol]].isna().sum(axis=1)
            
            if i==0:
                concatenated_nb_NA = nb_NA
                
            else:
                concatenated_nb_NA = pandas.concat([concatenated_nb_NA, nb_NA], axis=1)
        concatenated_nb_NA.columns = protocols[cell_type]
        
        concatenated_nb_NA["combination"] = concatenated_nb_NA[protocols[cell_type]].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
        
        
        def extract_type_combi(s):
            
            a = np.array(s.split("_")).astype(int)
            
            counts = np.unique(a, return_counts=True)
            
            dic_count = {0:0,
                         1:0,
                         2:0,
                         3:0}
            
            for i in range(len(counts[0])):
                dic_count[counts[0][i]] = counts[1][i]
                
            label = ", ".join([str(k)+":"+str(dic_count[k]) for k in [0,1,2,3]])
            
            return label
            
            
        
        concatenated_nb_NA["combination_type"] = concatenated_nb_NA["combination"].apply(lambda s: extract_type_combi(s))
  
        counts = np.unique(concatenated_nb_NA["combination_type"], return_counts=True)

        
        concatenated_nb_NA["sum all NA"] = (concatenated_nb_NA[protocols[cell_type]]==3).sum(axis=1)
        
        concatenated_nb_NA["sum at least 2 NA"] = (concatenated_nb_NA[protocols[cell_type]]>=2).sum(axis=1)
        
        
        mask = (concatenated_nb_NA["sum at least 2 NA"]<=len(protocols[cell_type])-2) 
                
        
        mask = mask.apply(lambda b: True)
                
        return mask
 

         

    # def get_binary_functional_table(self, table, list_characteristics=["KEGG_PATHWAY", "GOTERM_CC_DIRECT", "GOTERM_BP_DIRECT"]):

    #     list_cols = []
    
    #     for characteristic in list_characteristics:
    
    #         CC = table[characteristic].astype(str)
    #         all_CC = CC.apply(lambda s: s.split(", "))
    #         all_CC = sum(list(np.array(all_CC)),[])
            
    #         unique_CC = np.unique(all_CC, return_counts=True)
    
    #         for u, col in enumerate(unique_CC[0]):
     
    #             if col != "" and col!="nan" and col!="None":# and unique_CC[1][u]>0.01*len(CC): 
                
    #                 if "~" in col and "GO" in col:
    #                     name = characteristic + ": "+col.split('~')[1]
    #                 elif ':' in col:
                        
    #                     name = characteristic + ": "+col.split(':')[1]
    
    #                 else:
    #                     name = characteristic + ": "+col
    #                 table[name] = CC.apply(lambda s: 1 if col in s else 0)
    #                 list_cols.append(name)
    
    #     return table[list_cols]
        




    def get_background_count_functional_table(self, species, list_categories=["GOTERM_CC_DIRECT", "GOTERM_BP_DIRECT"], mode="unique", reset=False):
        
        path_to_save = Path(savepath, species+"_functional_count.csv")
        
        if reset or not os.path.isfile(path_to_save):
        
            table = self.load_background_genome_table(species)
            
            if mode=="unique":
    
                count_table = self.get_count_functional_table(table, list_categories = list_categories)
            
            elif mode=="combinations":

                count_table = self.get_count_combinations_functional_table(table, list_categories = list_categories)

            count_table.to_csv(path_to_save, index=False)
            
        count_table = pandas.read_csv(path_to_save)
        
        count_table.columns = ["Category", "Term", "Pop Hits", "Pop Total"]
        
        return count_table
    
    
    def compute_functional_chart(self, table, reference_species, list_categories=["GOTERM_CC_DIRECT", "GOTERM_BP_DIRECT"]):
        
        background_functional_table = self.get_background_count_functional_table(species)
        
        functional_table = self.get_count_functional_table(table, list_categories=list_categories)
        
        chart_table = functional_table.merge(background_functional_table, on="Term", how="outer")
        
        return chart_table
 


    def _get_gene_names_goterm_complete_genome(self, genome_table, species, category, term):
                
        where_term_complete_genome = genome_table[category].apply(lambda s: term in str(s))
        
        names = genome_table[where_term_complete_genome]["Gene name"].values
        
        return names

        
        
    def get_gene_names_goterms_complete_genome(self, species, reset=False):
        
        path_to_save = Path(savepath, "gene_names_goterm_complete_genome_"+species+".pkl")
        
        if reset or not os.path.exists(path_to_save):

            genome_table = self.load_background_genome_table(species=species)
            
            dic_terms = {}
            
            for categ in ["GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT"]:
                
                dic_terms[categ] = {}
                
                terms = genome_table[categ].apply(lambda s: str(s).replace(", "," - ").split(","))
       
                counter = Counter(chain.from_iterable(terms))
    
                terms = list(counter.keys())

                for term in terms:
    
                    names = self._get_gene_names_goterm_complete_genome(genome_table, species, categ, term)
                    
                    dic_terms[categ][term] = names
    
                with open(path_to_save, "wb") as file:
                    pickle.dump(dic_terms, file, protocol=pickle.HIGHEST_PROTOCOL)  

                   
        with open(path_to_save, "rb") as file:
            dic_terms = pickle.load(file)                  
                
        return dic_terms
                
                
            

        

 
    def get_count_functional_table(self, table, list_categories=["GOTERM_CC_DIRECT", "GOTERM_BP_DIRECT"]):

        for i, category in enumerate(list_categories):

            CC = table[category].astype(str)
          
            all_CC = CC.apply(lambda s: s.replace(", "," - ").split(","))

            counter = Counter(chain.from_iterable(all_CC))

            terms = list(counter.keys())

            terms = [term for term in terms if term!="nan" and term!="None"]

            count_terms = [counter[term] for term in terms]

            total_count_category = CC.apply(lambda s: s!="nan" and s!="None" and s!="").sum()

            categories = [category]*len(terms)
            total_counts = [total_count_category]*len(terms)
                
            counts = pandas.DataFrame(np.array([categories, terms, count_terms, total_counts]).T, 
                                      columns=["Category", "Term", "Hits", "Total"])
                
            counts["Hits"] = counts["Hits"].astype(int)
            counts["Total"] = counts["Total"].astype(int)

            
            if i==0:
                count_df = counts
            else:
                count_df = pandas.concat([count_df, counts], axis=0)

        return count_df


    def get_count_combinations_functional_table(self, table, list_categories=["GOTERM_CC_DIRECT", "GOTERM_BP_DIRECT"]):

        check_categories = table[list_categories].astype(str).applymap(lambda s: (s!="nan" and s!="None" and s!=""))
        
        check_all_categories = check_categories.product(axis=1).astype(bool)
        
        total_count = check_all_categories.sum()

        CC = table[list_categories[0]].astype(str).apply(lambda s: list_categories[0]+"_"+s)

        for i in range(1, len(list_categories)):
            CC = CC.str.cat(table[list_categories[i]].astype(str).apply(lambda s: list_categories[0]+"_"+s), sep=",")

        all_CC = CC.apply(lambda s: s.replace(", "," - ").split(","))
        
        def get_all_combinations_without_nan_or_None(row):
            # remove nan, None and double values
            set_without_nan = {value for value in row if isinstance(value, str)}

            # generate all possible combinations of the values in a row
            all_combinations = []
            for i in range(1, len(set_without_nan)+1):
                result = list(itertools.combinations(set_without_nan, i))
                all_combinations.extend(result)
                
            return all_combinations
            
        # get all posssible combinations of values in a row
        all_rows = all_CC.apply(get_all_combinations_without_nan_or_None, 1).values
        all_rows_flatten = list(itertools.chain.from_iterable(all_rows))
        
        # use Counter to count how many there are of each combination
        counter = Counter(all_rows_flatten)

#        counter = Counter(chain.from_iterable(all_combinations))

        terms = list(counter.keys())

        terms = [term for term in terms if "nan" not in term and "None" not in term]

        count_terms = [counter[term] for term in terms]

        categories = ["_".join(list_categories+["combinations"])]*len(terms)
        total_counts = [total_count]*len(terms)
        
        counts = pandas.DataFrame(np.array([categories, terms, count_terms, total_counts]).T, 
                                  columns=["Category", "Term", "Hits", "Total"])
            
        counts["Hits"] = counts["Hits"].astype(int)
        counts["Total"] = counts["Total"].astype(int)

        

        return counts
    



    def test_imputation_matrix(self, table, cell_type):
        
        path_file = Path(savepath, "R_results", "imputed_matrix_"+cell_type+".csv")
        
        DEP_path_file = Path(savepath, "R_results", "DEP_imputation_"+cell_type+".csv")
        DEP_imputation_table = pandas.read_csv(DEP_path_file)
        DEP_imputation_table.columns = [col.replace("X2D","2D").replace("X3D","3D") for col in DEP_imputation_table.columns]

        
        imput_table = pandas.read_csv(path_file).reset_index(drop=True)
        
        fig, ax = plt.subplots(imput_table.shape[1], 2, sharex=True)
        
        LFQI_cols = [col for col in table.columns if "LFQI" in col]

        values = np.array(table[LFQI_cols]).flatten()
        valid_values = values[~np.isnan(values)]
        
        
        # data = table[LFQI_cols]
        
        # mask = (data.isna().sum(axis=1)<int(data.shape[1]/2))
        
        # data = data[mask]
        # std = np.std(data, axis=1)
        # scale = np.median(std[~np.isnan(std)])
        
        
        for i, col in enumerate(imput_table.columns):
            
            values = np.array(table[col.replace(".", " ").replace("MS 6h", "MS.6h")]).flatten()
            valid_values = values[~np.isnan(values)]
            quantile = np.percentile(valid_values, 1)


    
            data = table[[column for column in table.columns if "LFQ" in column and col.split("_")[0] and column]]
            mask = (data.isna().sum(axis=1)<int(data.shape[1]/2))
            data = data[mask]
            std = np.std(data, axis=1)
            scale = np.median(std[~np.isnan(std)])            



            gaussian = norm(loc=quantile, scale=scale)
            
                
            
            nan_values = (table[col.replace(".", " ").replace("MS 6h", "MS.6h")].isna())

            imputed_values = gaussian.rvs(nan_values.sum())

            
            DEP_imputed_values = DEP_imputation_table.reset_index(drop=True)[nan_values][col.replace("LFQI.","").replace(".","_")]
            
         
            # ax[i, 0].hist(imput_table[nan_values][col], bins=100, alpha=0.3, color="red")
            ax[i, 0].hist(imputed_values, bins=100, alpha=0.3, color="blue")
            ax[i,0].hist(DEP_imputed_values, bins=100, alpha=0.3, color="green")
            
            ax[i, 1].hist(table[col.replace(".", " ").replace("MS 6h", "MS.6h")], bins=100)
            
            
            median = table[col.replace(".", " ").replace("MS 6h", "MS.6h")].median()
            
            ax[i,1].plot([median,median],[0, 40], color="white")
 

        

    def get_enrichment_results(self, table, species, threshold_pval=0.1, list_categories=["GOTERM_CC_DIRECT", "GOTERM_BP_DIRECT"],
                               only_enrichment=False, only_depletion=False, keep_significative=True):
        
        
        functional_count = self.get_count_functional_table(table, list_categories=list_categories)
            
        background_functional_count = self.get_background_count_functional_table(species)



#        elif mode=="combinations":
#            
#            background_functional_combinations_count = self.get_background_count_functional_table(species, mode="combinations")
#        
#            functional_combinations_count = self.get_count_functional_combinations_table(table)

        
        terms = []
        pvals = []
        fold_enrichments = []
        pop_total = []
        pop_hits = []
        count = []
        list_total = []
        categories = []

                
        for i in range(len(background_functional_count)):
            
            term = background_functional_count["Term"].iloc[i]
            
            category = background_functional_count["Category"].iloc[i]

            index_term_in_list = np.where((functional_count["Term"]==term)&(functional_count["Category"]==category))[0]

            if len(index_term_in_list)==0:
                index_category = np.where(functional_count["Category"]==category)[0][0]
                nb_total_in_list = functional_count["Total"].iloc[index_category]
                nb_in_list = 0
            else:
                index_term_in_list = index_term_in_list[0]
                nb_total_in_list = functional_count["Total"].iloc[index_term_in_list]
                nb_in_list = functional_count["Hits"].iloc[index_term_in_list]
                
#            if nb_in_list<2:
#                continue

            nb_total_in_genome = background_functional_count["Pop Total"].iloc[i]
            nb_in_genome = background_functional_count["Pop Hits"].iloc[i]

            contingency_table = np.array([[nb_in_list, nb_in_genome - nb_in_list], [nb_total_in_list - nb_in_list, nb_total_in_genome - nb_total_in_list - (nb_in_genome - nb_in_list)]])


            fold_enrichment = (nb_in_list/nb_total_in_list)/(nb_in_genome/nb_total_in_genome)

            if fold_enrichment>1:
                if only_depletion:
                    continue
                alternative="greater"
            else:
                if only_enrichment:
                    continue
                alternative="less"




            oddsr, p = fisher_exact(contingency_table, alternative=alternative)
                
                
            if p>=0.1 and keep_significative:
                continue

#            print(p)

            categories.append(category)
            terms.append(term)
            pvals.append(p)
            fold_enrichments.append(fold_enrichment)
            pop_total.append(nb_total_in_genome)
            pop_hits.append(nb_in_genome)
            list_total.append(nb_total_in_list)
            count.append(nb_in_list)
            
            
            
            if np.isnan(nb_in_list):
                print(contingency_table)
            
        corrected_pvals = fdrcorrection_twostage(np.array(pvals), alpha=0.01, method="bh")[1]

        array = np.array([categories, terms, pvals, fold_enrichments, corrected_pvals, pop_total, pop_hits, count, list_total]).T

        table = pandas.DataFrame(array, columns=["Category", "Term", "PValue", "Fold Enrichment", "Benjamini", "Pop Total", "Pop Hits", "Count", "List Total"])            

        return table
        





if __name__=="__main__":
    

    prodb = ProteoDataset(prepare_data=False, reset=True)
    
    table = prodb.load_table("MSCS")
    
    test = "MS-2D"
    
    mask = ~(table[[col for col in table.columns if "test_diff_DEP" in col]].isna()).product(axis=1).astype(bool) #Evaluates only where test_diff_DEP has been computed

    # print(np.sum(table[mask]["test_diff_ProDA_"+test]==table[mask]["test_diff_DEP_"+test]))
    
    # print(np.sum(table[mask]["test_diff_ProDA_"+test]!=table[mask]["test_diff_DEP_"+test]))
 
    # print(np.sum((np.abs(table[mask]["test_diff_ProDA_"+test])==1) & ((table[mask]["test_diff_DEP_"+test])==0)))
    # print(np.sum((np.abs(table[mask]["test_diff_ProDA_"+test])==0) & (np.abs(table[mask]["test_diff_DEP_"+test])==1)))


    # print(np.sum((table[mask]["test_diff_ProDA_"+test]==-1)&(table[mask]["test_diff_DEP_"+test]==1)))
    # print(np.sum((table[mask]["test_diff_ProDA_"+test]==1)&(table[mask]["test_diff_DEP_"+test]==-1)))
    