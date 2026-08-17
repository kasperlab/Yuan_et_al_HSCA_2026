######
py_run_string("adata.var_names_make_unique()")
py_run_string("adata.var['mt'] = adata.var_names.str.startswith('MT-')")
py_run_string("sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)")
######
py$adata$obs["StudyID"] = as.character(py$adata$obs[["ExperimentID"]])
# Make a shorter DonorID
py$adata$obs["DonorID_short"] = sapply(as.character(py$adata$obs[["DonorID"]]), function(x){
  unlist(strsplit(x, split = "___", fixed = T))[2]
})
######
py_run_string("adata.obs['merge_by'] = adata.obs[merge_by].agg('___'.join, axis=1).values") # Only use for record.
######
na_value_target = c("NA", "nan", "unknown", "True")
for(ii in c(py$factors, py$annotation)){
  this_factor = as.character(py$adata$obs[[ii]])
  this_factor[tolower(this_factor) %in% tolower(na_value_target)] = NA
  py_run_string(paste0("del adata.obs['", ii, "']"))
  py$adata$obs[ii] = as.factor(this_factor)
  rm(this_factor)
}
######
py_run_string("qc_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']")
py_run_string("sc.pl.violin(adata, qc_metrics, jitter=0.4, multi_panel=True, save='_before_QC.pdf')")
qc_state = "before_QC"
summary_list[[qc_state]] = list()
summary_list[[qc_state]][["dim"]] = c(py$adata$n_obs, py$adata$n_vars)
py_run_string("tmp_expressed_gene = adata.obs['n_genes_by_counts'].values.copy()")
py_run_string("tmp_expression = adata.obs['total_counts'].values.copy()")
summary_list[[qc_state]][["expressed_gene"]] = summary(py$tmp_expressed_gene)
summary_list[[qc_state]][["expression"]] = summary(py$tmp_expression)
py_run_string("del tmp_expressed_gene")
py_run_string("del tmp_expression")
for(ii in unlist(py$factors)){
  py_run_string(paste0("tmp_factor = adata.obs['", ii, "'].values.astype('str').copy()"))
  summary_list[[qc_state]][[ii]] = table(py$tmp_factor)
  if(length(summary_list[[qc_state]][[ii]]) == 0)(
    stop(Paste0("Factor ", ii, " consist with only NA!"))
  )
  py_run_string("del tmp_factor")
}
# Already filtered in pre-processing
######
# QC #
######
comparison_mat = matrix(ncol = 5, nrow = py$adata$n_obs,
                        dimnames = list(py$adata$obs_names$values,
                                        c("nUMI", "nFeatures", "MitoPrecentage", "XPrecentage", "YPrecentage")))
py_run_string("M_Genes = adata.var_names.str.startswith('MT-')")
py_run_string("tmp_numi = adata.X.sum(1)")
py_run_string("tmp_nfeatures = (adata.X > 0).sum(1)")
py_run_string("tmp_mito_pct = adata[:, M_Genes].X.sum(1) / tmp_numi")
py_run_string("tmp_x_pct = adata[:, np.intersect1d(X_Genes, adata.var_names.values)].X.sum(1) / tmp_numi")
py_run_string("tmp_y_pct = adata[:, np.intersect1d(Y_Genes, adata.var_names.values)].X.sum(1) / tmp_numi")
comparison_mat[, "nUMI"] = py$tmp_numi
comparison_mat[, "nFeatures"] = py$tmp_nfeatures
comparison_mat[, "MitoPrecentage"] = py$tmp_mito_pct
comparison_mat[, "XPrecentage"] = py$tmp_x_pct
comparison_mat[, "YPrecentage"] = py$tmp_y_pct
py_run_string("del tmp_numi")
py_run_string("del tmp_nfeatures")
py_run_string("del tmp_mito_pct")
py_run_string("del tmp_x_pct")
py_run_string("del tmp_y_pct")
comparison_df = cbind(as.data.frame(comparison_mat), cbind(as.character(py$adata$obs$SampleID), as.character(py$adata$obs$DonorID), as.character(py$adata$obs$ExperimentID)))
colnames(comparison_df) = c(colnames(comparison_mat), "SampleID", "DonorID", "ExperimentID")
p = StackedVlnPlot(comparison_df, "SampleID", colnames(comparison_mat), "ExperimentID", ymaxs = c(200000, max(comparison_df$nFeatures), 0.2, 0.1, 0.001))
ggsave(paste0(py$figure_dir, "/qc.pdf"), p, width = 1.5 * length(unique(comparison_df[, "SampleID"])) + 3, height = 12, limitsize = FALSE)
saveRDS(comparison_df, paste0(py$output_dir, "/comparison_df.rds"))
write.table(comparison_df, paste0(py$output_dir, "/comparison_df.tsv"), sep = "\t", row.names = T, col.names = T, quote = F)
######
cell_number = py$adata$n_obs
py$cell_sorting = py$adata$obs_names$values
saveRDS(py$cell_sorting, paste0(py$output_dir, "/cell_sorting.rds"))
print(summary_list)
