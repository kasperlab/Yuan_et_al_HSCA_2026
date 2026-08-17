#####################
# Cell identification
#####################
adata_imbalance_annotated_1_path = paste0(py$output_dir_cache, "/adata_imbalance_annotated_1.h5ad")
if(file.exists(adata_imbalance_annotated_1_path)){
  py_run_string(paste0("adata_imbalance = ad.read_h5ad('", adata_imbalance_annotated_1_path, "')"))
  py_run_string("annotation = annotation if np.isin('IntegratedAnnotation_1', annotation) else annotation + ['IntegratedAnnotation_1']")
  ###
  py$cell_sorting = readRDS(paste0(py$output_dir_cache, "/cell_sorting.rds"))
  py_run_string("tmp_IntegratedAnnotation = adata_imbalance.obs['IntegratedAnnotation_1'].astype(str).values")
  summary_list[[qc_state]][["IntegratedAnnotation_1"]] = list("AllCells" = table(py$tmp_IntegratedAnnotation))
  py_run_string("del tmp_IntegratedAnnotation")
}else{
  py_run_string("tmp_leiden_1_015 = adata_imbalance.obs['leiden_1_015'][cell_sorting].astype(str).values")
  leiden_1_015 = py$tmp_leiden_1_015
  py_run_string("del tmp_leiden_1_015")
  cell_annotation = rep(NA, py$adata_imbalance$n_obs)
  names(cell_annotation) = py$cell_sorting
  ######
  cell_annotation[leiden_1_015 == "0"] = "Keratinocyte" # KC
  cell_annotation[leiden_1_015 == "1"] = "Keratinocyte" # KC
  cell_annotation[leiden_1_015 == "2"] = "Fibroblast" # FIB
  cell_annotation[leiden_1_015 == "3"] = "Keratinocyte" # KC
  cell_annotation[leiden_1_015 == "4"] = "Mural Cell" # MUR
  cell_annotation[leiden_1_015 == "5"] = "Lymphocyte" # LYM
  cell_annotation[leiden_1_015 == "6"] = "Vascular Endothelial Cell" # VEC
  cell_annotation[leiden_1_015 == "7"] = "Myeloid Cell" # MYE
  cell_annotation[leiden_1_015 == "8"] = "Melanocyte" # MEL
  cell_annotation[leiden_1_015 == "9"] = "Keratinocyte" # KC
  cell_annotation[leiden_1_015 == "10"] = "Lymphatic Endothelial Cell" # LEC
  cell_annotation[leiden_1_015 == "11"] = "Keratinocyte" # KC
  cell_annotation[leiden_1_015 == "12"] = "Mast Cell" # MAST
  cell_annotation[leiden_1_015 == "13"] = "Schwann Cell" # SC
  #
  print(sum(is.na(cell_annotation)))
  py$tmp_clustering_result = cell_annotation
  ###
  py_run_string("adata_imbalance.obs['IntegratedAnnotation_1'] = tmp_clustering_result")
  py_run_string("del tmp_clustering_result")
  ###
  py_run_string("tmp_IntegratedAnnotation = adata_imbalance.obs['IntegratedAnnotation_1'].astype(str).values")
  summary_list[[qc_state]][["IntegratedAnnotation_1"]] = list("AllCells" = table(py$tmp_IntegratedAnnotation))
  py_run_string("del tmp_IntegratedAnnotation")
  ######
  for(ii in unlist(py$annotation)){
    py_run_string(paste0("tmp_this_annotation = adata_imbalance.obs['", ii, "'].astype(str).values"))
    summary_list[[qc_state]][[ii]] = table(py$tmp_this_annotation)
    py_run_string("del tmp_this_annotation")
  }
  py_run_string("annotation = annotation if np.isin('IntegratedAnnotation_1', annotation) else annotation + ['IntegratedAnnotation_1']")
  ######################
  # Cache - Querying - 2
  ######################
  if(use_cache){
    py_run_string(paste0("adata_imbalance.write('", adata_imbalance_annotated_1_path, "')"))
  }
}
###
# Ploting will tickle .strings_to_categoricals() automatically, cause error when accessing adata attributes in r
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['", paste(clustering_keys_1, collapse = "', '"), "'] + annotation, legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=7, save='_annotation.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['", paste(clustering_keys_1, collapse = "', '"), "'] + annotation, legend_loc='right margin', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=1, save='_annotation_1.pdf')"))
for(ii in c("select", "phase", unlist(py$annotation))){
  py_run_string(paste0("clustering_plot(adata_imbalance, '", ii, "', basis='X_umap_imbalance', size=", get_dot_size(cell_number), ", colorbar_loc=None, ncols=7, save='_detail_", ii, ".pdf')"))
}
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['SeparateAnnotation_1', 'KRT5', 'MGST1', 'COL1A2', 'RGS5', 'DCT', 'VWF', 'LYVE1', 'TPSAB1', 'IL1B', 'CD3D', 'HBA1'], layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Annotation_Signatures.pdf')"))
density_visualization("adata_imbalance", "umap_imbalance", "StudyID", "imbalance.pdf", "density_color", "6", as.character(get_dot_size(cell_number)))
density_visualization("adata_imbalance", "umap_imbalance", "Sex", "imbalance.pdf", "density_color", "5", as.character(get_dot_size(cell_number)))
density_visualization("adata_imbalance", "umap_imbalance", "IntegratedAnnotation_1", "imbalance.pdf", "density_color", "5", as.character(get_dot_size(cell_number)))
density_visualization("adata_imbalance", "umap_imbalance", "AnatomicalRegionLevel2", "imbalance.pdf", "density_color", "7", as.character(get_dot_size(cell_number)))
py_run_string("factor_distribution(adata_imbalance, factors + ['select', 'phase'], 'IntegratedAnnotation_1', figure_dir + '/Distribution_IntegratedAnnotation.pdf')")
py_run_string("factor_distribution_per_method(adata_imbalance, factors + ['select', 'phase'], 'IntegratedAnnotation_1', precent='y', path=figure_dir + '/Distribution_IntegratedAnnotation_divided.pdf')")
saveRDS(summary_list, paste0(py$output_dir, "/summary_list.rds"))
sink(paste0(py$output_dir, "/summary_list.txt"))
print(summary_list)
sink()
neighbor_list = list()
for(ii in c("clustering", "visualization")){
  this_neighbor_uns = py$adata_imbalance_hvg$uns[[ii]]
  neighbor_list[[ii]] = list("uns" = this_neighbor_uns,
                             "connectivities" = py$adata_imbalance_hvg$obsp[[this_neighbor_uns[["connectivities_key"]]]],
                             "distances" = py$adata_imbalance_hvg$obsp[[this_neighbor_uns[["distances_key"]]]])
}
saveRDS(neighbor_list, paste0(py$output_dir, "/all_neighbor_list.rds"))
############
py$loom_path_all = paste0(py$output_dir, "/all.loom")
py$loom_path_all_hvg = paste0(py$output_dir, "/all_hvg.loom")
metadata_path_all = paste0(py$output_dir, "/metadata_all.tsv")
if(file.exists(py$loom_path_all)){
  unlink(py$loom_path_all)
}
if(file.exists(py$loom_path_all_hvg)){
  unlink(py$loom_path_all_hvg)
}
py_run_string("loompy.create(loom_path_all, adata_imbalance.layers['raw_counts'].T, {'Gene': adata_imbalance.var_names.values}, {'CellID': adata_imbalance.obs_names.values})")
py_run_string("loompy.create(loom_path_all_hvg, adata_imbalance_hvg.layers['raw_counts'].T, {'Gene': adata_imbalance_hvg.var_names.values}, {'CellID': adata_imbalance_hvg.obs_names.values})")
py_run_string("tmp_metadata_all = adata_imbalance.obs.astype(str)")
py_run_string("del loom_path_all")
py_run_string("del loom_path_all_hvg")
write.table(py$tmp_metadata_all, metadata_path_all, sep = "\t", row.names = T, col.names = T, quote = F)
saveRDS(py$tmp_metadata_all, paste0(py$output_dir, "/metadata.rds"))
py_run_string("del tmp_metadata_all")
#############################
# End of Cell Type Clustering
#############################
py_run_string("annotation_1 = copy.deepcopy(annotation)")
annotation_1_path = paste0(py$output_dir, "/annotation_1.rds")
if(use_cache){
  saveRDS(py$annotation_1, annotation_1_path)
}
#####################
# Representative Cell
#####################
py_run_string("adata_representative = adata_imbalance[adata_imbalance.obs['select'] == 'Y', :].copy()")
py_run_string("adata_representative.write(output_dir + '/adata_representative_annotated.h5ad')")
