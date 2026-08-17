#####################
# Secondly Clustering
#####################


################
# Sub-clustering
################
#celltype_1 = c("Melanocyte", "Lymphocyte", "Keratinocyte", "Mural Cell", "Lymphatic Endothelial Cell", "Mast Cell", "Schwann Cell", "Myeloid Cell", "Vascular Endothelial Cell", "Fibroblast")
celltype_1 = as.character(unique(py$adata_imbalance$obs[["IntegratedAnnotation_1"]]))

#this_celltype = "Keratinocyte"
#this_celltype = "Fibroblast"
#this_celltype = "Melanocyte"
#this_celltype = "Myeloid Cell"
#this_celltype = "Vascular Endothelial Cell"
#this_celltype = "Lymphocyte"
#this_celltype = "Mural Cell"
#this_celltype = "Lymphatic Endothelial Cell"
#this_celltype = "Mast Cell"
upper_level_name = "IntegratedAnnotation_1"
for(this_celltype in celltype_1){
  this_celltype_name = gsub("[ -]", "", this_celltype)
  ##############
  #Configuration
  ##############
  #py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name)
  py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs")
  dir.create(py$this_output_dir, showWarnings = F, recursive = T)
  py_run_string("this_figure_dir = this_output_dir + '/figures'")
  py_run_string("sc.settings.figdir = this_figure_dir")
  py_run_string("annotation_tmp = copy.deepcopy(annotation_1)")
  ############
  # Clustering
  ############
  path_list[[this_celltype_name]] = list()
  path_list[[this_celltype_name]][["clustering"]] = paste0(py$this_output_dir, "/clustering.h5ad")
  path_list[[this_celltype_name]][["clustering_hvg"]] = paste0(py$this_output_dir, "/clustering_hvg.h5ad")
  path_list[[this_celltype_name]][["clustering_keys"]] = paste0(py$this_output_dir, "/clustering_keys.rds")
  if(file.exists(path_list[[this_celltype_name]][["clustering"]]) && file.exists(path_list[[this_celltype_name]][["clustering_hvg"]])){
    py_run_string(paste0("adata_", this_celltype_name, "_imbalance = ad.read_h5ad('", path_list[[this_celltype_name]][["clustering"]], "')"))
    py_run_string(paste0("adata_", this_celltype_name, "_imbalance_hvg = ad.read_h5ad('", path_list[[this_celltype_name]][["clustering_hvg"]], "')"))
    py_run_string(paste0("adata_", this_celltype_name, "_imbalance_hvg.uns['log1p']['base'] = None"))
    clustering_keys_2 = readRDS(path_list[[this_celltype_name]][["clustering_keys"]])
    this_pca_dim = min(c(pca_dim, eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance$n_obs"))) - 1))
  }else{
    py_run_string(paste0("adata_", this_celltype_name, " = adata_imbalance[adata_imbalance.obs['", upper_level_name, "'] == '", this_celltype, "', :].copy()"))
    py_run_string(paste0("adata_", this_celltype_name, ".X = adata_", this_celltype_name, ".layers['raw_counts'].copy()"))
    py_run_string(paste0("del adata_", this_celltype_name, ".uns['log1p']"))
    py_run_string(paste0("cell_sorting_", this_celltype_name, " = adata_", this_celltype_name, ".obs_names.values"))
    ############
    # Cell cycle
    ############
    py_run_string(paste0("sc.tl.score_genes_cell_cycle(adata_", this_celltype_name, ", s_genes=s_genes, g2m_genes=g2m_genes, random_state=", random_state, ")"))
    ############
    # Imbalanced
    ############
    # Clustering
    this_pca_dim = min(c(pca_dim, eval(parse(text = paste0("py$adata_", this_celltype_name, "$n_obs"))) - 1))
    pca_cal_list[[this_celltype_name]] = PCA_calculation_python(adata = eval(parse(text = paste0("py$adata_", this_celltype_name))), random_state = random_state, merge_by = unlist(py$merge_by), hvg_number = hvg_number, pca_dim = pca_dim, return_all = T)
    eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance_hvg = pca_cal_list[[this_celltype_name]][['hvg_matrix']]")))
    eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance = pca_cal_list[[this_celltype_name]][['complete_matrix']]")))
    py_run_string(paste0("sc.pl.pca_variance_ratio(adata_", this_celltype_name, "_imbalance_hvg, n_pcs=", pca_dim, ", log=True, save='_imbalance.pdf')"))
    py_run_string(paste0("sc.pp.neighbors(adata_", this_celltype_name, "_imbalance_hvg, n_neighbors=", n_neighbors_clustering, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='clustering_2')"))
    # Leiden Clustering
    leiden_clustering_keys_2 = c()
    for(ii in leiden_resolutions){
      this_key = paste0("leiden_2_", gsub(pattern = ".", replacement = "", sprintf("%.2f", ii), fixed = T))
      py_run_string(paste0("sc.tl.leiden(adata_", this_celltype_name, "_imbalance_hvg, resolution=", ii, ", key_added='", this_key, "', neighbors_key='clustering_2')"))
      py_run_string(paste0("adata_", this_celltype_name, "_imbalance.obs['", this_key, "'] = adata_", this_celltype_name, "_imbalance_hvg.obs['", this_key, "']"))
      leiden_clustering_keys_2 = c(leiden_clustering_keys_2, this_key)
    }
    # Visualization
    py_run_string(paste0("sc.pp.neighbors(adata_", this_celltype_name, "_imbalance_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
    py_run_string(paste0("sc.tl.umap(adata_", this_celltype_name, "_imbalance_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
    py_run_string(paste0("adata_", this_celltype_name, "_imbalance.obsm['X_umap_2_imbalance'] = adata_", this_celltype_name, "_imbalance_hvg.obsm['X_umap']"))
    clustering_keys_2 = c(leiden_clustering_keys_2)
    ####################
    # Cache - Clustering
    ####################
    if(use_cache){
      py_run_string(paste0("adata_", this_celltype_name, "_imbalance.write('", path_list[[this_celltype_name]][["clustering"]], "')"))
      py_run_string(paste0("adata_", this_celltype_name, "_imbalance_hvg.write('", path_list[[this_celltype_name]][["clustering_hvg"]], "')"))
      saveRDS(clustering_keys_2, path_list[[this_celltype_name]][["clustering_keys"]])
    }
  }
  ######
  cell_number = eval(parse(text = paste0("py$adata_", this_celltype_name, "$n_obs")))
  ##########
  # Bodysite
  ##########
  bodysite = eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance$obs[['AnatomicalRegionLevel3']]")))
  names(bodysite) = eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance$obs_names$values")))
  study_id = eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance$obs[['StudyID']]")))
  bodysite[study_id != "SS3_HumanSkin_20K"] = NA
  eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance$obs[['SS3_HumanSkin_20K_AnatomicalRegionLevel3']] = bodysite")))
  py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'X_umap_2_imbalance', color=['SS3_HumanSkin_20K_AnatomicalRegionLevel3'], size=", get_dot_size(cell_number), ", ncols=1, legend_loc='right margin', legend_fontsize=6, frameon=False, save='_SS3_HumanSkin_20K_AnatomicalRegionLevel3.pdf')"))
  ######
  py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'X_umap_2_imbalance', color=np.setdiff1d(factors, factors_only_binary).tolist() + ['phase'], size=", get_dot_size(cell_number), ", ncols=1, legend_loc='right margin', legend_fontsize=6, frameon=False, save='_factors.pdf')"))
  py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'X_umap_2_imbalance', color=['", paste(clustering_keys_2, collapse = "', '"), "'] + annotation, legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=8, frameon=False, save='_clustering.pdf')"))
  for(ii in c(py$factors, "phase", clustering_keys_2)){
    py_run_string(paste0("clustering_plot(adata_", this_celltype_name, "_imbalance, '", ii, "', basis='X_umap_2_imbalance', size=", get_dot_size(cell_number), ", colorbar_loc=None, ncols=7, save='_detail_", ii, ".pdf')"))
  }
  density_visualization(paste0("adata_", this_celltype_name, "_imbalance"), "umap_2_imbalance", "StudyID", ".pdf", "density_color", "7", as.character(get_dot_size(cell_number)))
  density_visualization(paste0("adata_", this_celltype_name, "_imbalance"), "umap_2_imbalance", "Sex", ".pdf", "density_color", "2", as.character(get_dot_size(cell_number)))
  density_visualization(paste0("adata_", this_celltype_name, "_imbalance"), "umap_2_imbalance", "AnatomicalRegionLevel1", ".pdf", "density_color", "4", as.character(get_dot_size(cell_number)))
  density_visualization(paste0("adata_", this_celltype_name, "_imbalance"), "umap_2_imbalance", "AnatomicalRegionLevel2", ".pdf", "density_color", "6", as.character(get_dot_size(cell_number)))
  density_visualization(paste0("adata_", this_celltype_name, "_imbalance"), "umap_2_imbalance", "AnatomicalRegionLevel3", ".pdf", "density_color", "7", as.character(get_dot_size(cell_number)))
  # Signatures
  py$Keratinocyte_markers = markers_list[["Keratinocyte"]][markers_list[["Keratinocyte"]] %in% eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance")))$var_names$values]
  py$Fibroblast_MuralCell_markers = markers_list[["Fibroblast_MuralCell"]][markers_list[["Fibroblast_MuralCell"]] %in% eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance")))$var_names$values]
  py$NeuralCrestderivedCell_markers = markers_list[["NeuralCrestderivedCell"]][markers_list[["NeuralCrestderivedCell"]] %in% eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance")))$var_names$values]
  py$EndothelialCell_markers = markers_list[["EndothelialCell"]][markers_list[["EndothelialCell"]] %in% eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance")))$var_names$values]
  py$ImmuneCell_markers = markers_list[["ImmuneCell"]][markers_list[["ImmuneCell"]] %in% eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance")))$var_names$values]
  py$Plasma_Erythrocyte_markers = markers_list[["Plasma_Erythrocyte"]][markers_list[["Plasma_Erythrocyte"]] %in% eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance")))$var_names$values]
  py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_2_imbalance', color=Keratinocyte_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=6, frameon=False, save='_Signatures_Keratinocyte.pdf')"))
  py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_2_imbalance', color=Fibroblast_MuralCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Signatures_Fibroblast_MuralCell.pdf')"))
  py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_2_imbalance', color=NeuralCrestderivedCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Signatures_NeuralCrestderivedCell.pdf')"))
  py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_2_imbalance', color=EndothelialCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Signatures_EndothelialCell.pdf')"))
  py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_2_imbalance', color=ImmuneCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=6, frameon=False, save='_Signatures_ImmuneCell.pdf')"))
  py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_2_imbalance', color=Plasma_Erythrocyte_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=4, frameon=False, save='_Signatures_Plasma_Erythrocyte.pdf')"))
}


celltype_1 = as.character(unique(py$adata_imbalance$obs[["IntegratedAnnotation_1"]]))
upper_level_name = "IntegratedAnnotation_1"
for(this_celltype in celltype_1){
  this_celltype_name = gsub("[ -]", "", this_celltype)
  ##############
  #Configuration
  ##############
  if(this_celltype == "Mast Cell"){
    py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs_supported_using_ImmuneCell")
    dir.create(py$this_output_dir, showWarnings = F, recursive = T)
    py_run_string("this_figure_dir = this_output_dir + '/figures'")
    py_run_string("sc.settings.figdir = this_figure_dir")
    py_run_string("annotation_tmp = copy.deepcopy(annotation_1)")
    #
    adata_str = "adata_MastCell"
    visualization_embedding = "UMAP_ImmuneCell_MastCell"
  }else{
    py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs")
    dir.create(py$this_output_dir, showWarnings = F, recursive = T)
    py_run_string("this_figure_dir = this_output_dir + '/figures'")
    py_run_string("sc.settings.figdir = this_figure_dir")
    py_run_string("annotation_tmp = copy.deepcopy(annotation_1)")
    #
    adata_str = paste0("adata_", this_celltype_name, "_imbalance")
    visualization_embedding = "X_umap_2_imbalance"
  }
  #
  cell_number = eval(parse(text = paste0("py$", adata_str, "$n_obs")))
  ##########
  # Bodysite
  ##########
  bodysite = as.character(eval(parse(text = paste0("py$", adata_str, "$obs[['AnatomicalRegionLevel3']]"))))
  names(bodysite) = eval(parse(text = paste0("py$", adata_str, "$obs_names$values")))
  study_id = eval(parse(text = paste0("py$", adata_str, "$obs[['StudyID']]")))
  bodysite[study_id != "SS3_HumanSkin_20K"] = NA
  bodysite_slot = "SS3_HumanSkin_20K_AnatomicalRegionLevel3"
  eval(parse(text = paste0("py$", adata_str, "$obs[['", bodysite_slot, "']] = bodysite")))
  py_run_string(paste0("tmp_this_factor = pd.Series(", adata_str, ".obs['", bodysite_slot, "'].astype(str).values)"))
  for(jj in tolower(c("NA", "nan", "unknown", "True"))){
    py_run_string(paste0("tmp_this_factor[tmp_this_factor.str.lower() == '", jj, "'] = pd.NA")) # Only pd.NA can make <NA> in Pandas (NA in R), np.nan will make NaN and None will make None.
  }
  py_run_string(paste0("del ", adata_str, ".obs['", bodysite_slot, "']"))
  py_run_string(paste0(adata_str, ".obs['", bodysite_slot, "'] = tmp_this_factor.values"))
  py_run_string("del tmp_this_factor")
  py$dotsize = max(c(40, get_dot_size(cell_number)))
  #py$dotsize = max(c(8, get_dot_size(cell_number)))
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", bodysite_slot, "'], size=dotsize, ncols=1, legend_loc='right margin', legend_fontsize=6, frameon=False, save='_SS3_HumanSkin_20K_AnatomicalRegionLevel3.pdf')"))
  py_run_string(paste0("clustering_plot(", adata_str, ", '", bodysite_slot, "', basis='", visualization_embedding, "', size=dotsize, colorbar_loc=None, ncols=5, save='_detail_", bodysite_slot, ".pdf')"))
  density_visualization(adata_str, visualization_embedding, bodysite_slot, ".pdf", "density_color", "5", as.character(py$dotsize))
  #
  eval(parse(text = paste0("StudyID = py$", adata_str, "$obs[['StudyID']]")))
  eval(parse(text = paste0("py$", adata_str, "$obs[['StudyID_ordered']] = factor(StudyID, levels = StudyID[!duplicated(StudyID)])")))
  py_run_string(paste0("clustering_plot(", adata_str, ", 'StudyID_ordered', basis='", visualization_embedding, "', size=dotsize, colorbar_loc=None, ncols=6, save='_detail_StudyID_ordered.pdf')"))
  density_visualization(adata_str, visualization_embedding, "StudyID_ordered", ".pdf", "density_color", "4", as.character(py$dotsize))
}
#############################################################################
#############################################################################



