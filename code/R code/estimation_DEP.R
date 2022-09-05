
library("dplyr")

library("DEP")

library("base")

library(SummarizedExperiment)
library(extraDistr)



cellular_type = "MSCSmapH"

datapath = "/home/nicolai/Documents/MSC-Med/proteomic-data-analysis/save/for_R/"

savepath = "/home/nicolai/Documents/MSC-Med/proteomic-data-analysis/save/R_results/"

file <- paste(datapath, cellular_type, "_table_for_analysis_DEP.csv", sep="")
dataset <- read.csv(file)
#dataset <- as.matrix(dataset, rownames = TRUE)

file_design_matrix = paste(datapath, cellular_type, "_design_matrix_for_analysis_DEP.csv", sep="")
design_matrix <- read.csv(file_design_matrix)
#design_matrix <- as.matrix(design_matrix, colnames = TRUE)

data_unique <- make_unique(dataset, "Gene.name", "Gene.IDs", delim = ";")

LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_se <- make_se(data_unique, LFQ_columns, design_matrix)

data_filt <- filter_missval(data_se, thr = 1)



#plot_missval(data_se)
#plot_frequency(data_se)
#plot_detect(data_se)




data_imp <- impute(data_filt, fun = "MinProb", q = 0.01)

imputed_data = assays(data_imp)[[1]]
write.table(imputed_data, paste(savepath, "DEP_imputation_",cellular_type,".csv", sep=""), sep=",")

metadata = attributes(data_imp)$elementMetadata
write.table(metadata, paste(savepath, "DEP_imputation_metadata_",cellular_type,".csv", sep=""), sep=",")


data_diff <- DEP::test_diff(data_imp, type = "all")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))
data_results <- get_results(dep)

write.table(data_results, paste(savepath, "test_diff_DEP_",cellular_type,".csv", sep=""), sep=",")

