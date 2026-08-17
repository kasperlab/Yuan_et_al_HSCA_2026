############
# Cell cycle
############
py_run_string(paste0("sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes, random_state=", random_state, ")"))
############
# Clustering
############
adata_imbalance_clustering_1_path = paste0(py$output_dir_cache, "/adata_imbalance_clustering_1.h5ad")
adata_imbalance_hvg_clustering_1_path = paste0(py$output_dir_cache, "/adata_imbalance_hvg_clustering_1.h5ad")
clustering_keys_1_path = paste0(py$output_dir_cache, "/clustering_keys_1.rds")
if(file.exists(adata_imbalance_clustering_1_path) && file.exists(adata_imbalance_hvg_clustering_1_path)){
  py_run_string(paste0("adata_imbalance = ad.read_h5ad('", adata_imbalance_clustering_1_path, "')"))
  py_run_string(paste0("adata_imbalance_hvg = ad.read_h5ad('", adata_imbalance_hvg_clustering_1_path, "')"))
  py_run_string("adata_imbalance.uns['log1p']['base'] = None")
  clustering_keys_1 = readRDS(clustering_keys_1_path)
}else{
  ############
  # Imbalanced
  ############
  # Clustering
  pca_cal_list[["All"]] = PCA_calculation_python(adata = py$adata, random_state = random_state, merge_by = unlist(py$merge_by), hvg_number = hvg_number, pca_dim = pca_dim, return_all = T)
  py$adata_imbalance_hvg = pca_cal_list[["All"]][["hvg_matrix"]]
  py$adata_imbalance = pca_cal_list[["All"]][["complete_matrix"]]
  py_run_string(paste0("sc.pl.pca_variance_ratio(adata_imbalance_hvg, n_pcs=", pca_dim, ", log=True, save='_imbalance.pdf')"))
  py_run_string(paste0("sc.pp.neighbors(adata_imbalance_hvg, n_neighbors=", n_neighbors_clustering, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='clustering')"))
  # Leiden Clustering
  leiden_clustering_keys_1 = c()
  for(ii in leiden_resolutions){
    this_key = paste0("leiden_1_", gsub(pattern = ".", replacement = "", sprintf("%.2f", ii), fixed = T))
    py_run_string(paste0("sc.tl.leiden(adata_imbalance_hvg, resolution=", ii, ", key_added='", this_key, "', neighbors_key='clustering')"))
    py_run_string(paste0("adata_imbalance.obs['", this_key, "'] = adata_imbalance_hvg.obs['", this_key, "']"))
    leiden_clustering_keys_1 = c(leiden_clustering_keys_1, this_key)
  }
  # Visualization
  py_run_string(paste0("sc.pp.neighbors(adata_imbalance_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_imbalance_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string("adata_imbalance.obsm['X_umap_imbalance'] = adata_imbalance_hvg.obsm['X_umap']")
  clustering_keys_1 = c(leiden_clustering_keys_1)
  ####################
  # Cache - Clustering
  ####################
  if(use_cache){
    py_run_string(paste0("adata_imbalance.write('", adata_imbalance_clustering_1_path, "')"))
    py_run_string(paste0("adata_imbalance_hvg.write('", adata_imbalance_hvg_clustering_1_path, "')"))
    saveRDS(clustering_keys_1, clustering_keys_1_path)
  }
}
py_run_string(paste0("sc.pl.pca_loadings(adata_imbalance_hvg, components = '", paste(seq(pca_dim), collapse = ", "), "', save='.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=np.setdiff1d(factors, factors_only_binary).tolist() + ['phase'], legend_loc='right margin', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=1, save='_factors.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['", paste(clustering_keys_1, collapse = "', '"), "'] + annotation, legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=7, save='_clustering.pdf')"))
for(ii in c(clustering_keys_1, unlist(py$factors))){
  py_run_string(paste0("clustering_plot(adata_imbalance, '", ii, "', basis='X_umap_imbalance', size=", get_dot_size(cell_number), ", colorbar_loc=None, ncols=7, save='_detail_", ii, ".pdf')"))
}
