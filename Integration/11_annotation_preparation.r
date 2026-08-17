######
py$Keratinocyte_markers = markers_list[["Keratinocyte"]][markers_list[["Keratinocyte"]] %in% py$adata_imbalance$var_names$values]
py$Fibroblast_MuralCell_markers = markers_list[["Fibroblast_MuralCell"]][markers_list[["Fibroblast_MuralCell"]] %in% py$adata_imbalance$var_names$values]
py$NeuralCrestderivedCell_markers = markers_list[["NeuralCrestderivedCell"]][markers_list[["NeuralCrestderivedCell"]] %in% py$adata_imbalance$var_names$values]
py$EndothelialCell_markers = markers_list[["EndothelialCell"]][markers_list[["EndothelialCell"]] %in% py$adata_imbalance$var_names$values]
py$ImmuneCell_markers = markers_list[["ImmuneCell"]][markers_list[["ImmuneCell"]] %in% py$adata_imbalance$var_names$values]
py$Plasma_Erythrocyte_markers = markers_list[["Plasma_Erythrocyte"]][markers_list[["Plasma_Erythrocyte"]] %in% py$adata_imbalance$var_names$values]
######
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=Keratinocyte_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=6, frameon=False, save='_Signatures_Keratinocyte.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=Fibroblast_MuralCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Signatures_Fibroblast_MuralCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=NeuralCrestderivedCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Signatures_NeuralCrestderivedCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=EndothelialCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Signatures_EndothelialCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=ImmuneCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=6, frameon=False, save='_Signatures_ImmuneCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=Plasma_Erythrocyte_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Signatures_Plasma_Erythrocyte.pdf')"))
######
##############
# DEG Analysis
##############
if(do_DEG){
  DEG_leiden_path_list[["All"]] = list()
  for(ii in c("leiden_1_015")){
    this_path = paste0(py$output_dir, "/DEG_", ii, ".rds")
    DEG_leiden_path_list[["All"]][[ii]] = this_path
    if(file.exists(this_path)){
      DEG_leiden_list[["All"]][[ii]] = readRDS(this_path)
    }else{
      DEG_leiden_list[["All"]][[ii]] = DEG_Analysis(py$adata_imbalance, ii, py$output_dir, prefix = paste0(ii, "_"), top_DE_range = top_DE_range)
      saveRDS(DEG_leiden_list[["All"]][[ii]], this_path)
    }
  }
}
