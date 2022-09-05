
#
#import biomart
#
#
#def convert_entrez_gene_id_to_gene_name():
#    
#    
#    server = biomart.BiomartServer('http://uswest.ensembl.org/biomart')         
#    mart = server.datasets['mmusculus_gene_ensembl']  
#
# 
#    
#    converter = {}
#    
#    # List the types of data we want                                            
#    attributes = ['entrez_gene_id', 'gene_name']
#    
#    
#    response = mart.search({'attributes': attributes})                          
#    
#    
#    for line in data.splitlines():                                              
#        line = line.split('\t')                                                 
#        # The entries are in the same order as in the `attributes` variable
#        entrez_gene_id = line[0]                                                 
#        gene_name = line[2]                                                  
#                                                                                
#        
#        converter[entrez_gene_id] = gene_name
#        
#        
        
        
        
        
        
import sys

from Bio import Entrez

Entrez.email = 'alice.nicolai@u-paris.fr'

## *Always* tell NCBI who you are
#Entrez.email = "your email here"


def retrieve_annotation(id_list):

    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    request = Entrez.epost("gene", id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        # FIXME: How generate NAs instead of causing an error with invalid IDs?
        print("An error occurred while retrieving the annotations.")
        print("The error returned was %s" % e)
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)
#
#    print("Retrieved %d annotations for %d genes" % (len(annotations["DocumentSummarySet"]["DocumentSummary"]), len(id_list)))

    return annotations["DocumentSummarySet"]["DocumentSummary"]

 


def convert_(annotation):
    pass
#    for gene_data in annotation:
#        gene_id = gene_data["Id"]
#        gene_symbol = gene_data["NomenclatureSymbol"]
#        gene_name = gene_data["Description"]
#        print "ID: %s - Gene Symbol: %s - Gene Name: %s" % (
#            gene_id,
#            gene_symbol,
#            gene_name,
#        )


def convert_entrez_gene_id_list_to_gene_name(entrez_gene_id_list):

    if len(entrez_gene_id_list)==0:
        return []
    
    annotations = retrieve_annotation(entrez_gene_id_list)
    

    
    gene_names_list = []
    for gene_data in annotations:
        gene_names_list.append(gene_data["Name"])
        
    return gene_names_list




def convert_entrez_gene_id_list_to_gene_id(entrez_gene_id_list):
    
    if len(entrez_gene_id_list)==0:
        return []
    
    annotations = retrieve_annotation(entrez_gene_id_list)
    
    print(annotations)
    

    
    gene_id_list = []
    for gene_data in annotations:
        gene_id_list.append(gene_data["Id"])
        
    return gene_id_list


       
        
if __name__=="__main__":
    
    
    gene_names_list = convert_entrez_gene_id_list_to_gene_name(entrez_gene_id_list=['148014'])
    
    gene_id_list = convert_entrez_gene_id_list_to_gene_id(entrez_gene_id_list=['148014'])