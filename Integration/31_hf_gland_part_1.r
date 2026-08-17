adata_imbalance_annotated_path = paste0("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss/IntegratedAnnotation_2.h5ad")
tryCatch({
  py_run_string(paste0("print(adata_imbalance)"))
}, error = function(x){
  py_run_string(paste0("adata_imbalance = ad.read_h5ad('", adata_imbalance_annotated_path, "')"))
  py_run_string("adata_imbalance.uns['log1p']['base'] = None")
  print("Loaded adata_imbalance from cache.")
})

HF_Gland_Karl_dt = fread("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss/HF_Gland_Cells_and_Metadata.csv", header = T)
HF_Gland_Karl_df = as.data.frame(HF_Gland_Karl_dt)
rownames(HF_Gland_Karl_df) = HF_Gland_Karl_df[, 1]
py$HF_Gland_Karl_cells = HF_Gland_Karl_df[, 1]

this_celltype = "HF_Gland"
this_celltype_name = gsub("_", "", gsub("[ -]", "", this_celltype), fixed = T)
adata_str = paste0("adata_", this_celltype_name, "_imbalance")
##############
#Configuration
##############
py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs_Karl_20250603")
dir.create(py$this_output_dir, showWarnings = F, recursive = T)
py_run_string("this_figure_dir = this_output_dir + '/figures'")
py_run_string("sc.settings.figdir = this_figure_dir")
############
# Clustering
############
path_list = list()
path_list[["clustering"]] = paste0(py$this_output_dir, "/clustering.h5ad")
path_list[["clustering_hvg"]] = paste0(py$this_output_dir, "/clustering_hvg.h5ad")
path_list[["clustering_keys"]] = paste0(py$this_output_dir, "/clustering_keys.rds")
if(file.exists(path_list[["clustering"]]) && file.exists(path_list[["clustering_hvg"]])){
  py_run_string(paste0(adata_str, " = ad.read_h5ad('", path_list[["clustering"]], "')"))
  py_run_string(paste0(adata_str, "_hvg = ad.read_h5ad('", path_list[["clustering_hvg"]], "')"))
  py_run_string(paste0(adata_str, "_hvg.uns['log1p']['base'] = None"))
  clustering_keys = readRDS(path_list[["clustering_keys"]])
  this_pca_dim = min(c(pca_dim, eval(parse(text = paste0(adata_str, "$n_obs"))) - 1))
}else{
  py_run_string(paste0("adata_", this_celltype_name, " = adata_imbalance[HF_Gland_Karl_cells].copy()"))
  py_run_string(paste0("adata_", this_celltype_name, ".X = adata_", this_celltype_name, ".layers['raw_counts'].copy()"))
  tryCatch({
    py_run_string(paste0("del adata_", this_celltype_name, ".uns['log1p']"))
  }, error = function(x){
    print(x)
  })
  eval(parse(text = paste0("py$adata_", this_celltype_name, "$obs['karl_clustering'] = as.character(HF_Gland_Karl_df[, 'leiden_res2_harmonyDonorID_bbknnLibraryPlatform_hvg2000-StudyID'])")))
  eval(parse(text = paste0("py$adata_", this_celltype_name, "$obs['karl_annotation'] = HF_Gland_Karl_df[, 'leiden_new_named']")))
  ###
  pca_cal_list = PCA_calculation_python(adata = eval(parse(text = paste0("py$adata_", this_celltype_name))), random_state = random_state, merge_by = unlist(py$merge_by), hvg_number = hvg_number, pca_dim = pca_dim, return_all = T)
  eval(parse(text = paste0("py$", adata_str, " = pca_cal_list[['complete_matrix']]")))
  eval(parse(text = paste0("py$", adata_str, "$obsm['PCA_use'] = pca_cal_list[['hvg_matrix']]$obsm['PCA_use']")))
  summary_list[["merge_by"]] = list("merge_by" = unlist(py$merge_by))
  rm(pca_cal_list)
  # Clustering
  clustering_embedding = "PCA_use"
  visualization_embedding = "X_umap_3_imbalance"
  clustering_keys = pipeline_post_PCA(adata_str, clustering_embedding, visualization_embedding,
                                      n_neighbors_clustering, n_neighbors_visualization, leiden_resolutions,
                                      UMAP_min_dist=0.5,
                                      clustering_prefix = "_leiden_",
                                      random_state = random_state)
  ####################
  # Cache - Clustering
  ####################
  if(use_cache){
    py_run_string(paste0(adata_str, ".write('", path_list[["clustering"]], "')"))
    saveRDS(clustering_keys, path_list[["clustering_keys"]])
  }
}
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", paste(clustering_keys, collapse = "', '"), "'] + ['IntegratedAnnotation_2', 'IntegratedAnnotation_3_v65', 'karl_clustering', 'karl_annotation'], legend_loc='on data', legend_fontsize=6, legend_fontoutline=0.5, size=8, ncols=6, frameon=False, save='_", adata_str, "_clustering_", clustering_embedding, "_OD.pdf')"))
for(ii in c('IntegratedAnnotation_2', 'IntegratedAnnotation_3_v65', 'karl_clustering', 'karl_annotation')){
  py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=8, ncols=5, save='_", adata_str, "_detail_", ii, ".pdf')"))
}





py_run_string("fig = sc.pl.embedding(adata_HFGland_imbalance,
                                     basis='X_umap_3_imbalance',
                                     legend_loc='right margin',
                                     legend_fontoutline=0.5,
                                     legend_fontsize=6,
                                     show=False,
                                     ncols=5,
                                     size=8,
                                     frameon=False,
                                     add_outline=False)")
py_run_string(paste0("sc.pl.embedding(adata_HFGland_imbalance[adata_HFGland_imbalance.obs['IntegratedAnnotation_3_v65'] == 'IRS/Cortex - KRT25/KRT27/KRT28/KRT35/KRT85'], '", visualization_embedding, "', color=['StudyID'], legend_loc='right margin', legend_fontsize=6, legend_fontoutline=0.5, size=8, ncols=6, frameon=False, ax=fig, save='_", adata_str, "_StudyID.pdf')"))
py_run_string("fig = sc.pl.embedding(adata_HFGland_imbalance,
                                     basis='X_umap_3_imbalance',
                                     legend_loc='right margin',
                                     legend_fontoutline=0.5,
                                     legend_fontsize=6,
                                     show=False,
                                     ncols=5,
                                     size=8,
                                     frameon=False,
                                     add_outline=False)")
py_run_string(paste0("sc.pl.embedding(adata_HFGland_imbalance[adata_HFGland_imbalance.obs['IntegratedAnnotation_3_v65'] == 'IRS/Cortex - KRT25/KRT27/KRT28/KRT35/KRT85'], '", visualization_embedding, "', color=['AnatomicalRegionLevel2'], legend_loc='right margin', legend_fontsize=6, legend_fontoutline=0.5, size=8, ncols=6, frameon=False, ax=fig, save='_", adata_str, "_AnatomicalRegionLevel2.pdf')"))
py_run_string("fig = sc.pl.embedding(adata_HFGland_imbalance,
                                     basis='X_umap_3_imbalance',
                                     legend_loc='right margin',
                                     legend_fontoutline=0.5,
                                     legend_fontsize=6,
                                     show=False,
                                     ncols=5,
                                     size=8,
                                     frameon=False,
                                     add_outline=False)")
py_run_string(paste0("sc.pl.embedding(adata_HFGland_imbalance[adata_HFGland_imbalance.obs['IntegratedAnnotation_3_v65'] == 'IRS/Cortex - KRT25/KRT27/KRT28/KRT35/KRT85'], '", visualization_embedding, "', color=['AnatomicalRegionLevel3'], legend_loc='right margin', legend_fontsize=6, legend_fontoutline=0.5, size=8, ncols=6, frameon=False, ax=fig, save='_", adata_str, "_AnatomicalRegionLevel3.pdf')"))


this_celltype = "HF_Gland_20_23"
this_celltype_name = gsub("_", "", gsub("[ -]", "", this_celltype), fixed = T)
adata_str = paste0("adata_", this_celltype_name, "_imbalance")
##############
#Configuration
##############
py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs_20250603")
dir.create(py$this_output_dir, showWarnings = F, recursive = T)
py_run_string("this_figure_dir = this_output_dir + '/figures'")
py_run_string("sc.settings.figdir = this_figure_dir")
############
# Clustering
############
path_list = list()
path_list[["clustering"]] = paste0(py$this_output_dir, "/clustering.h5ad")
path_list[["clustering_hvg"]] = paste0(py$this_output_dir, "/clustering_hvg.h5ad")
path_list[["clustering_keys"]] = paste0(py$this_output_dir, "/clustering_keys.rds")
if(file.exists(path_list[["clustering"]]) && file.exists(path_list[["clustering_hvg"]])){
  py_run_string(paste0(adata_str, " = ad.read_h5ad('", path_list[["clustering"]], "')"))
  py_run_string(paste0(adata_str, "_hvg = ad.read_h5ad('", path_list[["clustering_hvg"]], "')"))
  py_run_string(paste0(adata_str, "_hvg.uns['log1p']['base'] = None"))
  clustering_keys = readRDS(path_list[["clustering_keys"]])
  this_pca_dim = min(c(pca_dim, eval(parse(text = paste0(adata_str, "$n_obs"))) - 1))
}else{
  py_run_string("adata_HFGland_imbalance.obs['HFGland_leiden_080'] = adata_HFGland_imbalance.obs['PCA_use_leiden_080'].copy()")
  py_run_string(paste0("adata_", this_celltype_name, " = adata_HFGland_imbalance[adata_HFGland_imbalance.obs['PCA_use_leiden_080'].isin(['20', '23'])].copy()"))
  py_run_string(paste0("adata_", this_celltype_name, ".X = adata_", this_celltype_name, ".layers['raw_counts'].copy()"))
  tryCatch({
    py_run_string(paste0("del adata_", this_celltype_name, ".uns['log1p']"))
  }, error = function(x){
    print(x)
  })
  ###
  pca_cal_list = PCA_calculation_python(adata = eval(parse(text = paste0("py$adata_", this_celltype_name))), random_state = random_state, merge_by = NULL, hvg_number = hvg_number, pca_dim = pca_dim, return_all = T)
  eval(parse(text = paste0("py$", adata_str, " = pca_cal_list[['complete_matrix']]")))
  eval(parse(text = paste0("py$", adata_str, "$obsm['PCA_use'] = pca_cal_list[['hvg_matrix']]$obsm['PCA_use']")))
  summary_list[["merge_by"]] = list("merge_by" = unlist(py$merge_by))
  rm(pca_cal_list)
  # Clustering
  clustering_embedding = "PCA_use"
  visualization_embedding = "X_umap_3_imbalance"
  clustering_keys = pipeline_post_PCA(adata_str, clustering_embedding, visualization_embedding,
                                      n_neighbors_clustering, n_neighbors_visualization, leiden_resolutions,
                                      UMAP_min_dist=0.5,
                                      clustering_prefix = "_leiden_",
                                      random_state = random_state)
  ####################
  # Cache - Clustering
  ####################
  if(use_cache){
    py_run_string(paste0(adata_str, ".write('", path_list[["clustering"]], "')"))
    saveRDS(clustering_keys, path_list[["clustering_keys"]])
  }
}
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", paste(clustering_keys, collapse = "', '"), "'] + ['IntegratedAnnotation_2', 'IntegratedAnnotation_3_v65', 'karl_clustering', 'karl_annotation', 'HFGland_leiden_080', 'AnatomicalRegionLevel2', 'StudyID'], legend_loc='right margin', legend_fontsize=6, legend_fontoutline=0.5, size=8, ncols=1, frameon=False, save='_", adata_str, "_clustering_", clustering_embedding, ".pdf')"))
for(ii in c('IntegratedAnnotation_2', 'IntegratedAnnotation_3_v65', 'karl_clustering', 'karl_annotation')){
  py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=8, ncols=5, save='_", adata_str, "_detail_", ii, ".pdf')"))
}
