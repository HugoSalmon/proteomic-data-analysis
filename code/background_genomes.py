from pathlib import Path
from paths import savepath, datapath
import pandas
import os
import pickle
import numpy as np

from utils.bioinf_conversion import convert_entrez_gene_id_list_to_gene_name




path_background_genomes = {"homo_sapiens": Path(datapath, "background_genome", "DAVIDKnowledgebase_homo_sapiens"),
                            "mus_musculus": Path(datapath, "background_genome", "DAVIDKnowledgebase_mus_musculus")}


def load_background_genome_table(species, reset=False):

    path_to_save = Path(savepath, "background_genome_table_"+species+".pkl")

    if reset or not os.path.exists(path_to_save):

        directory = str(path_background_genomes[species])
        
        files = ["ENTREZ_GENE_ID2GOTERM_BP_DIRECT.txt",
                 "ENTREZ_GENE_ID2GOTERM_CC_DIRECT.txt",
                 "ENTREZ_GENE_ID2UNIPROT_ACCESSION.txt"]

        tables = {}

        for i in range(len(files)):
    
            path = files[i]
            
            with open(Path(directory, path), "r") as file:
                data = file.read()
                
            data = data.split("\n")
            data = [line.split("\t") for line in data]
            
            table = pandas.DataFrame(data)
            
            table.columns = ["ENTREZ_GENE_ID", files[i].replace("ENTREZ_GENE_ID2","").replace(".txt","")]

            tables[file] = table
    
            table = table.groupby('ENTREZ_GENE_ID').agg(lambda x: x.tolist())

            table.reset_index(inplace=True)
   #
            table.iloc[:,1] = table.iloc[:,1].apply(lambda l: ",".join(l).replace(", "," - ") if type(l)==list and l[0]!=None else l[0])

            table.columns = ["ENTREZ_GENE_ID", files[i].replace("ENTREZ_GENE_ID2","").replace(".txt","")]
            
            if i==0:
                concatenated_table = table
            else:
                concatenated_table = concatenated_table.merge(table, how="outer", on="ENTREZ_GENE_ID")


        concatenated_table = concatenated_table[concatenated_table["ENTREZ_GENE_ID"]!=""]
        concatenated_table.reset_index(inplace=True, drop=True)

        list_gene_id = concatenated_table["ENTREZ_GENE_ID"].values.tolist()
        
        first_split = list_gene_id[:10000]
        second_split = list_gene_id[10000:20000]
        third_split = list_gene_id[20000:]

        gene_names_first_split = convert_entrez_gene_id_list_to_gene_name(first_split)
        gene_names_second_split = convert_entrez_gene_id_list_to_gene_name(second_split)
        gene_names_third_split = convert_entrez_gene_id_list_to_gene_name(third_split) #one is not found

        gene_names = gene_names_first_split + gene_names_second_split + gene_names_third_split
        
        concatenated_table["GENE_NAME"] = np.array(gene_names)




        with open(path_to_save, "wb") as file:
            pickle.dump(concatenated_table, file, protocol=pickle.HIGHEST_PROTOCOL)  
            
            

               
    with open(path_to_save, "rb") as file:
        concatenated_table = pickle.load(file)   
        

    return concatenated_table


# def convert_protein_ID_to_gene_infos(array_protein_ID, species):
    
#     background_genome = load_background_genome_table(species)
    
#     background_genome = background_genome.explode("UNIPROT_ACCESSION")
    
    
#     print(background_genome["ENTREZ_GENE_ID"].tolist())
    


        
#     df = pandas.DataFrame(array_protein_ID.reshape(-1,1), columns=["UNIPROT_ACCESSION"])

#     df['ENTREZ_GENE_ID_FROM_GENOME'] = df['UNIPROT_ACCESSION'].map(background_genome.set_index('UNIPROT_ACCESSION')['ENTREZ_GENE_ID'])
#     df['GENE_NAME_FROM_GENOME'] = df['UNIPROT_ACCESSION'].map(background_genome.set_index('UNIPROT_ACCESSION')['GENE_NAME'])


        
#     return df


if __name__=="__main__":
    
    load_background_genome_table("homo_sapiens", reset=False)
    load_background_genome_table("mus_musculus", reset=True)
    

#    