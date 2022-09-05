

library("base")

library("proDA")



cellular_type = "THP1"

datapath = "/home/nicolai/Documents/MSC-Med/proteomic-data-analysis/save/for_R/"

savepath = "/home/nicolai/Documents/MSC-Med/proteomic-data-analysis/save/R_results/"

file <- paste(datapath, cellular_type, "_table_for_analysis.csv", sep="")
dataset <- read.csv(file)
dataset <- as.matrix(dataset, rownames = TRUE)

file_design_matrix = paste(datapath, cellular_type, "_design_matrix_for_analysis.csv", sep="")
design_matrix <- read.csv(file_design_matrix)
design_matrix <- as.matrix(design_matrix, rownames = TRUE, colnames = TRUE)

# Generate some dataset with known structure
#syn_dataset <- generate_synthetic_data(n_proteins = 100, n_conditions = 2)   # abundance matrix

#syn_dataset$groups

#Y <- syn_dataset$Y
#fit <- proDA(Y, design = syn_dataset$groups)

row.names(dataset) <- NULL



fit <- proDA::proDA(dataset, design=design_matrix)

write.table(fit$coefficients, paste(savepath, "coeffs_", cellular_type, ".csv", sep=""), sep=",") 

hyperparams = fit$hyper_parameters
priors = c(hyperparams$location_prior_mean, hyperparams$location_prior_scale, hyperparams$location_prior_df, hyperparams$variance_prior_scale)

write.table(priors, paste(savepath, "priors_",cellular_type,".csv",sep=""), sep=",") 


dropout_params = fit$colData

write.table(dropout_params, paste(savepath, "dropout_params_",cellular_type,".csv",sep=""), sep=",")


#coef_variance = attributes(fit)$elementMetadata$coef_var

coef_variance = fit$coefficient_variance_matrices
#coef_variance_processed = lapply(coef_variance, diag)
#coef_variance_processed = do.call(rbind, coef_variance_processed)

for (i in 1:length(coef_variance)) {
  print(i)
  if (any(is.na(coef_variance[i]))){
    print(coef_variance[i],i)}
}
k=1
for (covar in coef_variance) {
  
  if (any(is.na(covar))){
    print(k)
    print(covar)}
k=k+1
  }

write.table(do.call(rbind, coef_variance), paste(savepath, "coef_var_", cellular_type, ".csv", sep=""), sep=",")


df = attributes(fit)$elementMetadata$df
s2 = attributes(fit)$elementMetadata$s2

write.table(df, paste(savepath, "df_",cellular_type,".csv", sep=""), sep=",")
write.table(s2, paste(savepath, "s2_",cellular_type,".csv", sep=""), sep=",")


#distances = dist_approx(fit, by_sample = FALSE)
#write.table(as.matrix(distances$mean), paste(savepath, "distances_protein_",cellular_type,".csv",""), sep=",")


test_res <- proDA::test_diff(fit, "HS-X2D")
write.table(test_res, paste(savepath, "test_HS_2D_", cellular_type, ".csv", sep=""), sep=",")

#write.table(test_res$t_statistic, paste(savepath, "test_HS_2D_t_statistic_", cellular_type, ".csv", sep=""), sep=",")

test_res <- proDA::test_diff(fit, "HS-X2D")
write.table(test_res, paste(savepath, "test_HS_2D_", cellular_type, ".csv", sep=""), sep=",")

test_res <- proDA::test_diff(fit, "HS-X3D")
write.table(test_res, paste(savepath, "test_HS_3D_", cellular_type, ".csv", sep=""), sep=",")

test_res <- proDA::test_diff(fit, "X3D-X2D")
write.table(test_res, paste(savepath, "test_3D_2D_", cellular_type, ".csv", sep=""), sep=",")

test_res <- proDA::test_diff(fit, "HS-MS")
write.table(test_res, paste(savepath, "test_HS_MS_", cellular_type, ".csv", sep=""), sep=",")

test_res <- proDA::test_diff(fit, "MS-X3D")
write.table(test_res, paste(savepath, "test_MS_3D_", cellular_type, ".csv", sep=""), sep=",")

test_res <- proDA::test_diff(fit, "MS-X2D")
write.table(test_res, paste(savepath, "test_MS_2D_", cellular_type, ".csv", sep=""), sep=",")

test_res <- proDA::test_diff(fit, "MS.6h-X2D")
write.table(test_res, paste(savepath, "test_MS.6h_2D_", cellular_type, ".csv", sep=""), sep=",")

test_res <- proDA::test_diff(fit, "MS.6h-X3D")
write.table(test_res, paste(savepath, "test_MS.6h_3D_", cellular_type, ".csv", sep=""), sep=",")

test_res <- proDA::test_diff(fit, "MS.6h-MS")
write.table(test_res, paste(savepath, "test_MS.6h_MS_", cellular_type, ".csv", sep=""), sep=",")

test_res <- proDA::test_diff(fit, "HS-MS.6h")
write.table(test_res, paste(savepath, "test_HS_MS.6h_", cellular_type, ".csv", sep=""), sep=",")


