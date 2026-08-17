# SS3_HumanSkin_20K
adata_2_path = paste0(py$output_dir_cache, "/adata_2.h5ad")
adata_2_hvg_path = paste0(py$output_dir_cache, "/adata_2_hvg.h5ad")
if(file.exists(adata_2_path) && file.exists(adata_2_hvg_path)){
  py_run_string(paste0("adata_2 = ad.read_h5ad('", adata_2_path, "')"))
  py_run_string(paste0("adata_2_hvg = ad.read_h5ad('", adata_2_hvg_path, "')"))
}else{
  py_run_string("merge_by_2 = ['DonorID']")
  adata_2_list = load_data_wrap(loom_path_2, metadata_path_2, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_2), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_2 = adata_2_list[["adata"]]
  py$adata_2_hvg = adata_2_list[["hvg"]]
  rm(adata_2_list)
  py_run_string(paste0("sc.pp.neighbors(adata_2_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_2_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_2_hvg, 'select', size=8, save='_data_2_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_2_hvg, 'select', 'PCA_use', size=8, save='_data_2_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_2_hvg.obsm['X_pca_use'] = adata_2_hvg.obsm['PCA_use']")
  density_visualization("adata_2_hvg", "pca_use", "select", paste0("_data_2_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_2).write('", adata_2_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_2_hvg).write('", adata_2_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["SS3_HumanSkin_20K"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_2$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["SS3_HumanSkin_20K"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_2$obs[["SeparateAnnotation_1"]][py$adata_2$obs[['select']] == 'Y']))
# Ning
adata_3_path = paste0(py$output_dir_cache, "/adata_3.h5ad")
adata_3_hvg_path = paste0(py$output_dir_cache, "/adata_3_hvg.h5ad")
if(file.exists(adata_3_path) && file.exists(adata_3_hvg_path)){
  py_run_string(paste0("adata_3 = ad.read_h5ad('", adata_3_path, "')"))
  py_run_string(paste0("adata_3_hvg = ad.read_h5ad('", adata_3_hvg_path, "')"))
}else{
  py_run_string("merge_by_3 = ['DonorID']")
  adata_3_list = load_data_wrap(loom_path_3, metadata_path_3, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_3), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_3 = adata_3_list[["adata"]]
  py$adata_3_hvg = adata_3_list[["hvg"]]
  rm(adata_3_list)
  py_run_string(paste0("sc.pp.neighbors(adata_3_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_3_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_3_hvg, 'select', size=8, save='_data_3_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_3_hvg, 'select', 'PCA_use', size=8, save='_data_3_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_3_hvg.obsm['X_pca_use'] = adata_3_hvg.obsm['PCA_use']")
  density_visualization("adata_3_hvg", "pca_use", "select", paste0("_data_3_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_3).write('", adata_3_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_3_hvg).write('", adata_3_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Ning"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_3$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Ning"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_3$obs[["SeparateAnnotation_1"]][py$adata_3$obs[['select']] == 'Y']))
# Gaydosik_Fuschiotti_ClinicalCancerResearch_2019
adata_4_path = paste0(py$output_dir_cache, "/adata_4.h5ad")
adata_4_hvg_path = paste0(py$output_dir_cache, "/adata_4_hvg.h5ad")
if(file.exists(adata_4_path) && file.exists(adata_4_hvg_path)){
  py_run_string(paste0("adata_4 = ad.read_h5ad('", adata_4_path, "')"))
  py_run_string(paste0("adata_4_hvg = ad.read_h5ad('", adata_4_hvg_path, "')"))
}else{
  py_run_string("merge_by_4 = ['DonorID']")
  adata_4_list = load_data_wrap(loom_path_4, metadata_path_4, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_4), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_4 = adata_4_list[["adata"]]
  py$adata_4_hvg = adata_4_list[["hvg"]]
  rm(adata_4_list)
  py_run_string(paste0("sc.pp.neighbors(adata_4_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_4_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_4_hvg, 'select', size=8, save='_data_4_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_4_hvg, 'select', 'PCA_use', size=8, save='_data_4_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_4_hvg.obsm['X_pca_use'] = adata_4_hvg.obsm['PCA_use']")
  density_visualization("adata_4_hvg", "pca_use", "select", paste0("_data_4_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_4).write('", adata_4_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_4_hvg).write('", adata_4_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Gaydosik_Fuschiotti_ClinicalCancerResearch_2019"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_4$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Gaydosik_Fuschiotti_ClinicalCancerResearch_2019"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_4$obs[["SeparateAnnotation_1"]][py$adata_4$obs[['select']] == 'Y']))
# Wang_Atwood_NatureCommunications_2020
adata_5_path = paste0(py$output_dir_cache, "/adata_5.h5ad")
adata_5_hvg_path = paste0(py$output_dir_cache, "/adata_5_hvg.h5ad")
if(file.exists(adata_5_path) && file.exists(adata_5_hvg_path)){
  py_run_string(paste0("adata_5 = ad.read_h5ad('", adata_5_path, "')"))
  py_run_string(paste0("adata_5_hvg = ad.read_h5ad('", adata_5_hvg_path, "')"))
}else{
  py_run_string("merge_by_5 = ['DonorID']")
  adata_5_list = load_data_wrap(loom_path_5, metadata_path_5, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_5), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_5 = adata_5_list[["adata"]]
  py$adata_5_hvg = adata_5_list[["hvg"]]
  rm(adata_5_list)
  py_run_string(paste0("sc.pp.neighbors(adata_5_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_5_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_5_hvg, 'select', size=8, save='_data_5_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_5_hvg, 'select', 'PCA_use', size=8, save='_data_5_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_5_hvg.obsm['X_pca_use'] = adata_5_hvg.obsm['PCA_use']")
  density_visualization("adata_5_hvg", "pca_use", "select", paste0("_data_5_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_5).write('", adata_5_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_5_hvg).write('", adata_5_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Wang_Atwood_NatureCommunications_2020"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_5$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Wang_Atwood_NatureCommunications_2020"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_5$obs[["SeparateAnnotation_1"]][py$adata_5$obs[['select']] == 'Y']))
# Takahashi_Lowry_JournalofInvestigativeDermatology_2020
adata_6_path = paste0(py$output_dir_cache, "/adata_6.h5ad")
adata_6_hvg_path = paste0(py$output_dir_cache, "/adata_6_hvg.h5ad")
if(file.exists(adata_6_path) && file.exists(adata_6_hvg_path)){
  py_run_string(paste0("adata_6 = ad.read_h5ad('", adata_6_path, "')"))
  py_run_string(paste0("adata_6_hvg = ad.read_h5ad('", adata_6_hvg_path, "')"))
}else{
  py_run_string("merge_by_6 = ['DonorID', 'LibraryPlatform']")
  adata_6_list = load_data_wrap(loom_path_6, metadata_path_6, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_6), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_6 = adata_6_list[["adata"]]
  py$adata_6_hvg = adata_6_list[["hvg"]]
  rm(adata_6_list)
  py_run_string(paste0("sc.pp.neighbors(adata_6_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_6_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_6_hvg, 'select', size=8, save='_data_6_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_6_hvg, 'select', 'PCA_use', size=8, save='_data_6_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_6_hvg.obsm['X_pca_use'] = adata_6_hvg.obsm['PCA_use']")
  density_visualization("adata_6_hvg", "pca_use", "select", paste0("_data_6_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_6).write('", adata_6_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_6_hvg).write('", adata_6_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Takahashi_Lowry_JournalofInvestigativeDermatology_2020"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_6$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Takahashi_Lowry_JournalofInvestigativeDermatology_2020"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_6$obs[["SeparateAnnotation_1"]][py$adata_6$obs[['select']] == 'Y']))
# Belote_Torres_NatureCellBiology_2021
adata_7_path = paste0(py$output_dir_cache, "/adata_7.h5ad")
adata_7_hvg_path = paste0(py$output_dir_cache, "/adata_7_hvg.h5ad")
if(file.exists(adata_7_path) && file.exists(adata_7_hvg_path)){
  py_run_string(paste0("adata_7 = ad.read_h5ad('", adata_7_path, "')"))
  py_run_string(paste0("adata_7_hvg = ad.read_h5ad('", adata_7_hvg_path, "')"))
}else{
  py_run_string("merge_by_7 = ['DonorID']")
  adata_7_list = load_data_wrap(loom_path_7, metadata_path_7, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_7), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_7 = adata_7_list[["adata"]]
  py$adata_7_hvg = adata_7_list[["hvg"]]
  rm(adata_7_list)
  py_run_string(paste0("sc.pp.neighbors(adata_7_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_7_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_7_hvg, 'select', size=8, save='_data_7_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_7_hvg, 'select', 'PCA_use', size=8, save='_data_7_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_7_hvg.obsm['X_pca_use'] = adata_7_hvg.obsm['PCA_use']")
  density_visualization("adata_7_hvg", "pca_use", "select", paste0("_data_7_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_7).write('", adata_7_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_7_hvg).write('", adata_7_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Belote_Torres_NatureCellBiology_2021"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_7$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Belote_Torres_NatureCellBiology_2021"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_7$obs[["SeparateAnnotation_1"]][py$adata_7$obs[['select']] == 'Y']))
# Boldo_Lyko_CommunicationsBiology_2020
adata_8_path = paste0(py$output_dir_cache, "/adata_8.h5ad")
adata_8_hvg_path = paste0(py$output_dir_cache, "/adata_8_hvg.h5ad")
if(file.exists(adata_8_path) && file.exists(adata_8_hvg_path)){
  py_run_string(paste0("adata_8 = ad.read_h5ad('", adata_8_path, "')"))
  py_run_string(paste0("adata_8_hvg = ad.read_h5ad('", adata_8_hvg_path, "')"))
}else{
  py_run_string("merge_by_8 = ['DonorID']")
  adata_8_list = load_data_wrap(loom_path_8, metadata_path_8, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_8), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_8 = adata_8_list[["adata"]]
  py$adata_8_hvg = adata_8_list[["hvg"]]
  rm(adata_8_list)
  py_run_string(paste0("sc.pp.neighbors(adata_8_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_8_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_8_hvg, 'select', size=8, save='_data_8_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_8_hvg, 'select', 'PCA_use', size=8, save='_data_8_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_8_hvg.obsm['X_pca_use'] = adata_8_hvg.obsm['PCA_use']")
  density_visualization("adata_8_hvg", "pca_use", "select", paste0("_data_8_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_8).write('", adata_8_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_8_hvg).write('", adata_8_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Boldo_Lyko_CommunicationsBiology_2020"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_8$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Boldo_Lyko_CommunicationsBiology_2020"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_8$obs[["SeparateAnnotation_1"]][py$adata_8$obs[['select']] == 'Y']))
# Cheng_Cho_CellReports_2018
adata_9_path = paste0(py$output_dir_cache, "/adata_9.h5ad")
adata_9_hvg_path = paste0(py$output_dir_cache, "/adata_9_hvg.h5ad")
if(file.exists(adata_9_path) && file.exists(adata_9_hvg_path)){
  py_run_string(paste0("adata_9 = ad.read_h5ad('", adata_9_path, "')"))
  py_run_string(paste0("adata_9_hvg = ad.read_h5ad('", adata_9_hvg_path, "')"))
}else{
  py_run_string("merge_by_9 = ['DonorID']")
  adata_9_list = load_data_wrap(loom_path_9, metadata_path_9, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_9), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_9 = adata_9_list[["adata"]]
  py$adata_9_hvg = adata_9_list[["hvg"]]
  rm(adata_9_list)
  py_run_string(paste0("sc.pp.neighbors(adata_9_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_9_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_9_hvg, 'select', size=8, save='_data_9_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_9_hvg, 'select', 'PCA_use', size=8, save='_data_9_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_9_hvg.obsm['X_pca_use'] = adata_9_hvg.obsm['PCA_use']")
  density_visualization("adata_9_hvg", "pca_use", "select", paste0("_data_9_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_9).write('", adata_9_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_9_hvg).write('", adata_9_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Cheng_Cho_CellReports_2018"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_9$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Cheng_Cho_CellReports_2018"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_9$obs[["SeparateAnnotation_1"]][py$adata_9$obs[['select']] == 'Y']))
# Ji_Khavari_Cell_2020
adata_10_path = paste0(py$output_dir_cache, "/adata_10.h5ad")
adata_10_hvg_path = paste0(py$output_dir_cache, "/adata_10_hvg.h5ad")
if(file.exists(adata_10_path) && file.exists(adata_10_hvg_path)){
  py_run_string(paste0("adata_10 = ad.read_h5ad('", adata_10_path, "')"))
  py_run_string(paste0("adata_10_hvg = ad.read_h5ad('", adata_10_hvg_path, "')"))
}else{
  py_run_string("merge_by_10 = ['DonorID']")
  adata_10_list = load_data_wrap(loom_path_10, metadata_path_10, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_10), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_10 = adata_10_list[["adata"]]
  py$adata_10_hvg = adata_10_list[["hvg"]]
  rm(adata_10_list)
  py_run_string(paste0("sc.pp.neighbors(adata_10_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_10_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_10_hvg, 'select', size=8, save='_data_10_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_10_hvg, 'select', 'PCA_use', size=8, save='_data_10_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_10_hvg.obsm['X_pca_use'] = adata_10_hvg.obsm['PCA_use']")
  density_visualization("adata_10_hvg", "pca_use", "select", paste0("_data_10_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_10).write('", adata_10_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_10_hvg).write('", adata_10_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Ji_Khavari_Cell_2020"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_10$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Ji_Khavari_Cell_2020"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_10$obs[["SeparateAnnotation_1"]][py$adata_10$obs[['select']] == 'Y']))
# Tabib_Lafyatis_NatureCommunications_2021
adata_11_path = paste0(py$output_dir_cache, "/adata_11.h5ad")
adata_11_hvg_path = paste0(py$output_dir_cache, "/adata_11_hvg.h5ad")
if(file.exists(adata_11_path) && file.exists(adata_11_hvg_path)){
  py_run_string(paste0("adata_11 = ad.read_h5ad('", adata_11_path, "')"))
  py_run_string(paste0("adata_11_hvg = ad.read_h5ad('", adata_11_hvg_path, "')"))
}else{
  py_run_string("merge_by_11 = ['DonorID', 'LibraryPlatform']")
  adata_11_list = load_data_wrap(loom_path_11, metadata_path_11, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_11), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_11 = adata_11_list[["adata"]]
  py$adata_11_hvg = adata_11_list[["hvg"]]
  rm(adata_11_list)
  py_run_string(paste0("sc.pp.neighbors(adata_11_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_11_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_11_hvg, 'select', size=8, save='_data_11_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_11_hvg, 'select', 'PCA_use', size=8, save='_data_11_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_11_hvg.obsm['X_pca_use'] = adata_11_hvg.obsm['PCA_use']")
  density_visualization("adata_11_hvg", "pca_use", "select", paste0("_data_11_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_11).write('", adata_11_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_11_hvg).write('", adata_11_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Tabib_Lafyatis_NatureCommunications_2021"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_11$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Tabib_Lafyatis_NatureCommunications_2021"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_11$obs[["SeparateAnnotation_1"]][py$adata_11$obs[['select']] == 'Y']))
# He_Guttman_JournalofAllergyandClinicalImmunology_2020
adata_12_path = paste0(py$output_dir_cache, "/adata_12.h5ad")
adata_12_hvg_path = paste0(py$output_dir_cache, "/adata_12_hvg.h5ad")
if(file.exists(adata_12_path) && file.exists(adata_12_hvg_path)){
  py_run_string(paste0("adata_12 = ad.read_h5ad('", adata_12_path, "')"))
  py_run_string(paste0("adata_12_hvg = ad.read_h5ad('", adata_12_hvg_path, "')"))
}else{
  py_run_string("merge_by_12 = ['DonorID']")
  adata_12_list = load_data_wrap(loom_path_12, metadata_path_12, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_12), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_12 = adata_12_list[["adata"]]
  py$adata_12_hvg = adata_12_list[["hvg"]]
  rm(adata_12_list)
  py_run_string(paste0("sc.pp.neighbors(adata_12_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_12_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_12_hvg, 'select', size=8, save='_data_12_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_12_hvg, 'select', 'PCA_use', size=8, save='_data_12_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_12_hvg.obsm['X_pca_use'] = adata_12_hvg.obsm['PCA_use']")
  density_visualization("adata_12_hvg", "pca_use", "select", paste0("_data_12_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_12).write('", adata_12_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_12_hvg).write('", adata_12_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["He_Guttman_JournalofAllergyandClinicalImmunology_2020"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_12$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["He_Guttman_JournalofAllergyandClinicalImmunology_2020"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_12$obs[["SeparateAnnotation_1"]][py$adata_12$obs[['select']] == 'Y']))
# Zou_Liu_DevelopmentalCell_2020
adata_13_path = paste0(py$output_dir_cache, "/adata_13.h5ad")
adata_13_hvg_path = paste0(py$output_dir_cache, "/adata_13_hvg.h5ad")
if(file.exists(adata_13_path) && file.exists(adata_13_hvg_path)){
  py_run_string(paste0("adata_13 = ad.read_h5ad('", adata_13_path, "')"))
  py_run_string(paste0("adata_13_hvg = ad.read_h5ad('", adata_13_hvg_path, "')"))
}else{
  py_run_string("merge_by_13 = ['DonorID']")
  adata_13_list = load_data_wrap(loom_path_13, metadata_path_13, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_13), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_13 = adata_13_list[["adata"]]
  py$adata_13_hvg = adata_13_list[["hvg"]]
  rm(adata_13_list)
  py_run_string(paste0("sc.pp.neighbors(adata_13_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_13_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_13_hvg, 'select', size=8, save='_data_13_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_13_hvg, 'select', 'PCA_use', size=8, save='_data_13_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_13_hvg.obsm['X_pca_use'] = adata_13_hvg.obsm['PCA_use']")
  density_visualization("adata_13_hvg", "pca_use", "select", paste0("_data_13_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_13).write('", adata_13_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_13_hvg).write('", adata_13_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Zou_Liu_DevelopmentalCell_2020"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_13$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Zou_Liu_DevelopmentalCell_2020"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_13$obs[["SeparateAnnotation_1"]][py$adata_13$obs[['select']] == 'Y']))
# Wiedemann_Andersen_CellReports_2023
adata_14_path = paste0(py$output_dir_cache, "/adata_14.h5ad")
adata_14_hvg_path = paste0(py$output_dir_cache, "/adata_14_hvg.h5ad")
if(file.exists(adata_14_path) && file.exists(adata_14_hvg_path)){
  py_run_string(paste0("adata_14 = ad.read_h5ad('", adata_14_path, "')"))
  py_run_string(paste0("adata_14_hvg = ad.read_h5ad('", adata_14_hvg_path, "')"))
}else{
  py_run_string("merge_by_14 = ['DonorID']")
  adata_14_list = load_data_wrap(loom_path_14, metadata_path_14, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_14), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_14 = adata_14_list[["adata"]]
  py$adata_14_hvg = adata_14_list[["hvg"]]
  rm(adata_14_list)
  py_run_string(paste0("sc.pp.neighbors(adata_14_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_14_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_14_hvg, 'select', size=8, save='_data_14_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_14_hvg, 'select', 'PCA_use', size=8, save='_data_14_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_14_hvg.obsm['X_pca_use'] = adata_14_hvg.obsm['PCA_use']")
  density_visualization("adata_14_hvg", "pca_use", "select", paste0("_data_14_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_14).write('", adata_14_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_14_hvg).write('", adata_14_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Wiedemann_Andersen_CellReports_2023"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_14$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Wiedemann_Andersen_CellReports_2023"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_14$obs[["SeparateAnnotation_1"]][py$adata_14$obs[['select']] == 'Y']))
# Alkon_Stingl_JournalofAllergyandClinicalImmunology_2022
adata_15_path = paste0(py$output_dir_cache, "/adata_15.h5ad")
adata_15_hvg_path = paste0(py$output_dir_cache, "/adata_15_hvg.h5ad")
if(file.exists(adata_15_path) && file.exists(adata_15_hvg_path)){
  py_run_string(paste0("adata_15 = ad.read_h5ad('", adata_15_path, "')"))
  py_run_string(paste0("adata_15_hvg = ad.read_h5ad('", adata_15_hvg_path, "')"))
}else{
  py_run_string("merge_by_15 = ['DonorID']")
  adata_15_list = load_data_wrap(loom_path_15, metadata_path_15, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_15), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_15 = adata_15_list[["adata"]]
  py$adata_15_hvg = adata_15_list[["hvg"]]
  rm(adata_15_list)
  py_run_string(paste0("sc.pp.neighbors(adata_15_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_15_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_15_hvg, 'select', size=8, save='_data_15_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_15_hvg, 'select', 'PCA_use', size=8, save='_data_15_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_15_hvg.obsm['X_pca_use'] = adata_15_hvg.obsm['PCA_use']")
  density_visualization("adata_15_hvg", "pca_use", "select", paste0("_data_15_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_15).write('", adata_15_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_15_hvg).write('", adata_15_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Alkon_Stingl_JournalofAllergyandClinicalImmunology_2022"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_15$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Alkon_Stingl_JournalofAllergyandClinicalImmunology_2022"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_15$obs[["SeparateAnnotation_1"]][py$adata_15$obs[['select']] == 'Y']))
# Dunlap_Rao_JCIInsight_2022
adata_16_path = paste0(py$output_dir_cache, "/adata_16.h5ad")
adata_16_hvg_path = paste0(py$output_dir_cache, "/adata_16_hvg.h5ad")
if(file.exists(adata_16_path) && file.exists(adata_16_hvg_path)){
  py_run_string(paste0("adata_16 = ad.read_h5ad('", adata_16_path, "')"))
  py_run_string(paste0("adata_16_hvg = ad.read_h5ad('", adata_16_hvg_path, "')"))
}else{
  py_run_string("merge_by_16 = ['DonorID']")
  adata_16_list = load_data_wrap(loom_path_16, metadata_path_16, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_16), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_16 = adata_16_list[["adata"]]
  py$adata_16_hvg = adata_16_list[["hvg"]]
  rm(adata_16_list)
  py_run_string(paste0("sc.pp.neighbors(adata_16_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_16_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_16_hvg, 'select', size=8, save='_data_16_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_16_hvg, 'select', 'PCA_use', size=8, save='_data_16_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_16_hvg.obsm['X_pca_use'] = adata_16_hvg.obsm['PCA_use']")
  density_visualization("adata_16_hvg", "pca_use", "select", paste0("_data_16_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_16).write('", adata_16_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_16_hvg).write('", adata_16_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Dunlap_Rao_JCIInsight_2022"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_16$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Dunlap_Rao_JCIInsight_2022"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_16$obs[["SeparateAnnotation_1"]][py$adata_16$obs[['select']] == 'Y']))
# Boothby_Rosenblum_Nature_2021
adata_17_path = paste0(py$output_dir_cache, "/adata_17.h5ad")
adata_17_hvg_path = paste0(py$output_dir_cache, "/adata_17_hvg.h5ad")
if(file.exists(adata_17_path) && file.exists(adata_17_hvg_path)){
  py_run_string(paste0("adata_17 = ad.read_h5ad('", adata_17_path, "')"))
  py_run_string(paste0("adata_17_hvg = ad.read_h5ad('", adata_17_hvg_path, "')"))
}else{
  py_run_string("merge_by_17 = ['DonorID']")
  adata_17_list = load_data_wrap(loom_path_17, metadata_path_17, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_17), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_17 = adata_17_list[["adata"]]
  py$adata_17_hvg = adata_17_list[["hvg"]]
  rm(adata_17_list)
  py_run_string(paste0("sc.pp.neighbors(adata_17_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_17_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_17_hvg, 'select', size=8, save='_data_17_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_17_hvg, 'select', 'PCA_use', size=8, save='_data_17_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_17_hvg.obsm['X_pca_use'] = adata_17_hvg.obsm['PCA_use']")
  density_visualization("adata_17_hvg", "pca_use", "select", paste0("_data_17_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_17).write('", adata_17_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_17_hvg).write('", adata_17_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Boothby_Rosenblum_Nature_2021"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_17$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Boothby_Rosenblum_Nature_2021"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_17$obs[["SeparateAnnotation_1"]][py$adata_17$obs[['select']] == 'Y']))
# Gao_Hu_CellDeathDisease_2021
adata_18_path = paste0(py$output_dir_cache, "/adata_18.h5ad")
adata_18_hvg_path = paste0(py$output_dir_cache, "/adata_18_hvg.h5ad")
if(file.exists(adata_18_path) && file.exists(adata_18_hvg_path)){
  py_run_string(paste0("adata_18 = ad.read_h5ad('", adata_18_path, "')"))
  py_run_string(paste0("adata_18_hvg = ad.read_h5ad('", adata_18_hvg_path, "')"))
}else{
  py_run_string("merge_by_18 = ['DonorID']")
  adata_18_list = load_data_wrap(loom_path_18, metadata_path_18, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_18), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_18 = adata_18_list[["adata"]]
  py$adata_18_hvg = adata_18_list[["hvg"]]
  rm(adata_18_list)
  py_run_string(paste0("sc.pp.neighbors(adata_18_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_18_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_18_hvg, 'select', size=8, save='_data_18_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_18_hvg, 'select', 'PCA_use', size=8, save='_data_18_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_18_hvg.obsm['X_pca_use'] = adata_18_hvg.obsm['PCA_use']")
  density_visualization("adata_18_hvg", "pca_use", "select", paste0("_data_18_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_18).write('", adata_18_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_18_hvg).write('", adata_18_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Gao_Hu_CellDeathDisease_2021"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_18$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Gao_Hu_CellDeathDisease_2021"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_18$obs[["SeparateAnnotation_1"]][py$adata_18$obs[['select']] == 'Y']))
# Hughes_Shalek_Immunity_2020
adata_19_path = paste0(py$output_dir_cache, "/adata_19.h5ad")
adata_19_hvg_path = paste0(py$output_dir_cache, "/adata_19_hvg.h5ad")
if(file.exists(adata_19_path) && file.exists(adata_19_hvg_path)){
  py_run_string(paste0("adata_19 = ad.read_h5ad('", adata_19_path, "')"))
  py_run_string(paste0("adata_19_hvg = ad.read_h5ad('", adata_19_hvg_path, "')"))
}else{
  py_run_string("merge_by_19 = ['DonorID']")
  adata_19_list = load_data_wrap(loom_path_19, metadata_path_19, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_19), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_19 = adata_19_list[["adata"]]
  py$adata_19_hvg = adata_19_list[["hvg"]]
  rm(adata_19_list)
  py_run_string(paste0("sc.pp.neighbors(adata_19_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_19_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_19_hvg, 'select', size=8, save='_data_19_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_19_hvg, 'select', 'PCA_use', size=8, save='_data_19_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_19_hvg.obsm['X_pca_use'] = adata_19_hvg.obsm['PCA_use']")
  density_visualization("adata_19_hvg", "pca_use", "select", paste0("_data_19_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_19).write('", adata_19_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_19_hvg).write('", adata_19_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Hughes_Shalek_Immunity_2020"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_19$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Hughes_Shalek_Immunity_2020"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_19$obs[["SeparateAnnotation_1"]][py$adata_19$obs[['select']] == 'Y']))
# Kim_Nagao_NatureMedicine_2020
adata_20_path = paste0(py$output_dir_cache, "/adata_20.h5ad")
adata_20_hvg_path = paste0(py$output_dir_cache, "/adata_20_hvg.h5ad")
if(file.exists(adata_20_path) && file.exists(adata_20_hvg_path)){
  py_run_string(paste0("adata_20 = ad.read_h5ad('", adata_20_path, "')"))
  py_run_string(paste0("adata_20_hvg = ad.read_h5ad('", adata_20_hvg_path, "')"))
}else{
  py_run_string("merge_by_20 = ['SampleID', 'LibraryPlatform']")
  adata_20_list = load_data_wrap(loom_path_20, metadata_path_20, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_20), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_20 = adata_20_list[["adata"]]
  py$adata_20_hvg = adata_20_list[["hvg"]]
  rm(adata_20_list)
  py_run_string(paste0("sc.pp.neighbors(adata_20_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_20_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_20_hvg, 'select', size=8, save='_data_20_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_20_hvg, 'select', 'PCA_use', size=8, save='_data_20_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_20_hvg.obsm['X_pca_use'] = adata_20_hvg.obsm['PCA_use']")
  density_visualization("adata_20_hvg", "pca_use", "select", paste0("_data_20_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_20).write('", adata_20_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_20_hvg).write('", adata_20_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Kim_Nagao_NatureMedicine_2020"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_20$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Kim_Nagao_NatureMedicine_2020"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_20$obs[["SeparateAnnotation_1"]][py$adata_20$obs[['select']] == 'Y']))
# Rindler_Brunner_MolecularCancer_2021
adata_21_path = paste0(py$output_dir_cache, "/adata_21.h5ad")
adata_21_hvg_path = paste0(py$output_dir_cache, "/adata_21_hvg.h5ad")
if(file.exists(adata_21_path) && file.exists(adata_21_hvg_path)){
  py_run_string(paste0("adata_21 = ad.read_h5ad('", adata_21_path, "')"))
  py_run_string(paste0("adata_21_hvg = ad.read_h5ad('", adata_21_hvg_path, "')"))
}else{
  py_run_string("merge_by_21 = ['DonorID']")
  adata_21_list = load_data_wrap(loom_path_21, metadata_path_21, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_21), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_21 = adata_21_list[["adata"]]
  py$adata_21_hvg = adata_21_list[["hvg"]]
  rm(adata_21_list)
  py_run_string(paste0("sc.pp.neighbors(adata_21_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_21_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_21_hvg, 'select', size=8, save='_data_21_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_21_hvg, 'select', 'PCA_use', size=8, save='_data_21_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_21_hvg.obsm['X_pca_use'] = adata_21_hvg.obsm['PCA_use']")
  density_visualization("adata_21_hvg", "pca_use", "select", paste0("_data_21_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_21).write('", adata_21_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_21_hvg).write('", adata_21_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Rindler_Brunner_MolecularCancer_2021"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_21$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Rindler_Brunner_MolecularCancer_2021"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_21$obs[["SeparateAnnotation_1"]][py$adata_21$obs[['select']] == 'Y']))
# Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020
adata_22_path = paste0(py$output_dir_cache, "/adata_22.h5ad")
adata_22_hvg_path = paste0(py$output_dir_cache, "/adata_22_hvg.h5ad")
if(file.exists(adata_22_path) && file.exists(adata_22_hvg_path)){
  py_run_string(paste0("adata_22 = ad.read_h5ad('", adata_22_path, "')"))
  py_run_string(paste0("adata_22_hvg = ad.read_h5ad('", adata_22_hvg_path, "')"))
}else{
  py_run_string("merge_by_22 = ['DonorID', 'LibraryPlatform']")
  adata_22_list = load_data_wrap(loom_path_22, metadata_path_22, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_22), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_22 = adata_22_list[["adata"]]
  py$adata_22_hvg = adata_22_list[["hvg"]]
  rm(adata_22_list)
  py_run_string(paste0("sc.pp.neighbors(adata_22_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_22_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_22_hvg, 'select', size=8, save='_data_22_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_22_hvg, 'select', 'PCA_use', size=8, save='_data_22_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_22_hvg.obsm['X_pca_use'] = adata_22_hvg.obsm['PCA_use']")
  density_visualization("adata_22_hvg", "pca_use", "select", paste0("_data_22_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_22).write('", adata_22_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_22_hvg).write('", adata_22_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_22$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_22$obs[["SeparateAnnotation_1"]][py$adata_22$obs[['select']] == 'Y']))
# Theocharidis_Bhasin_NatureCommunications_2022
adata_23_path = paste0(py$output_dir_cache, "/adata_23.h5ad")
adata_23_hvg_path = paste0(py$output_dir_cache, "/adata_23_hvg.h5ad")
if(file.exists(adata_23_path) && file.exists(adata_23_hvg_path)){
  py_run_string(paste0("adata_23 = ad.read_h5ad('", adata_23_path, "')"))
  py_run_string(paste0("adata_23_hvg = ad.read_h5ad('", adata_23_hvg_path, "')"))
}else{
  py_run_string("merge_by_23 = ['SampleID']")
  adata_23_list = load_data_wrap(loom_path_23, metadata_path_23, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_23), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_23 = adata_23_list[["adata"]]
  py$adata_23_hvg = adata_23_list[["hvg"]]
  rm(adata_23_list)
  py_run_string(paste0("sc.pp.neighbors(adata_23_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_23_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_23_hvg, 'select', size=8, save='_data_23_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_23_hvg, 'select', 'PCA_use', size=8, save='_data_23_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_23_hvg.obsm['X_pca_use'] = adata_23_hvg.obsm['PCA_use']")
  density_visualization("adata_23_hvg", "pca_use", "select", paste0("_data_23_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_23).write('", adata_23_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_23_hvg).write('", adata_23_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Theocharidis_Bhasin_NatureCommunications_2022"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_23$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Theocharidis_Bhasin_NatureCommunications_2022"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_23$obs[["SeparateAnnotation_1"]][py$adata_23$obs[['select']] == 'Y']))
# Vorstandlechner_Mildner_NatureCommunications_2021
adata_24_path = paste0(py$output_dir_cache, "/adata_24.h5ad")
adata_24_hvg_path = paste0(py$output_dir_cache, "/adata_24_hvg.h5ad")
if(file.exists(adata_24_path) && file.exists(adata_24_hvg_path)){
  py_run_string(paste0("adata_24 = ad.read_h5ad('", adata_24_path, "')"))
  py_run_string(paste0("adata_24_hvg = ad.read_h5ad('", adata_24_hvg_path, "')"))
}else{
  py_run_string("merge_by_24 = ['DonorID']")
  adata_24_list = load_data_wrap(loom_path_24, metadata_path_24, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_24), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_24 = adata_24_list[["adata"]]
  py$adata_24_hvg = adata_24_list[["hvg"]]
  rm(adata_24_list)
  py_run_string(paste0("sc.pp.neighbors(adata_24_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_24_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_24_hvg, 'select', size=8, save='_data_24_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_24_hvg, 'select', 'PCA_use', size=8, save='_data_24_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_24_hvg.obsm['X_pca_use'] = adata_24_hvg.obsm['PCA_use']")
  density_visualization("adata_24_hvg", "pca_use", "select", paste0("_data_24_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_24).write('", adata_24_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_24_hvg).write('", adata_24_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Vorstandlechner_Mildner_NatureCommunications_2021"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_24$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Vorstandlechner_Mildner_NatureCommunications_2021"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_24$obs[["SeparateAnnotation_1"]][py$adata_24$obs[['select']] == 'Y']))
# Xu_Chen_Nature_2022
adata_25_path = paste0(py$output_dir_cache, "/adata_25.h5ad")
adata_25_hvg_path = paste0(py$output_dir_cache, "/adata_25_hvg.h5ad")
if(file.exists(adata_25_path) && file.exists(adata_25_hvg_path)){
  py_run_string(paste0("adata_25 = ad.read_h5ad('", adata_25_path, "')"))
  py_run_string(paste0("adata_25_hvg = ad.read_h5ad('", adata_25_hvg_path, "')"))
}else{
  py_run_string("merge_by_25 = ['DonorID']")
  adata_25_list = load_data_wrap(loom_path_25, metadata_path_25, max_number = max_number, sample_method = "geosketch", random_state = random_state, merge_by = unlist(py$merge_by_25), hvg_number = hvg_number, pca_dim = pca_dim, keep = T, alias_unifying = alias_gene_c_seurat, return_all = T)
  py$adata_25 = adata_25_list[["adata"]]
  py$adata_25_hvg = adata_25_list[["hvg"]]
  rm(adata_25_list)
  py_run_string(paste0("sc.pp.neighbors(adata_25_hvg, n_neighbors=", n_neighbors_visualization, ", n_pcs=", pca_dim, ", use_rep='PCA_use', random_state=", random_state, ", key_added='visualization')"))
  py_run_string(paste0("sc.tl.umap(adata_25_hvg, random_state=", random_state, ", neighbors_key='visualization')"))
  py_run_string(paste0("clustering_plot(adata_25_hvg, 'select', size=8, save='_data_25_select_geosketch_", max_number, ".pdf')"))
  py_run_string(paste0("clustering_plot(adata_25_hvg, 'select', 'PCA_use', size=8, save='_data_25_select_geosketch_", max_number, ".pdf')"))
  py_run_string("adata_25_hvg.obsm['X_pca_use'] = adata_25_hvg.obsm['PCA_use']")
  density_visualization("adata_25_hvg", "pca_use", "select", paste0("_data_25_select_geosketch_", max_number, ".pdf"), "density_color")
  if(use_cache){
    py_run_string(paste0("adata_obs_to_string(adata_25).write('", adata_25_path, "')"))
    py_run_string(paste0("adata_obs_to_string(adata_25_hvg).write('", adata_25_hvg_path, "')"))
  }
}
summary_list[["Selection"]][["Xu_Chen_Nature_2022"]] = list ("all_cells" = list("SeparateAnnotation_1" = table(py$adata_25$obs[["SeparateAnnotation_1"]])))
summary_list[["Selection"]][["Xu_Chen_Nature_2022"]][["selected_cells"]] = list("SeparateAnnotation_1" = table(py$adata_25$obs[["SeparateAnnotation_1"]][py$adata_25$obs[['select']] == 'Y']))
######
dir.create(py$figure_dir, showWarnings = F, recursive = T)
######
# Statistics of the selection
cell_type = c()
datasets = names(summary_list[["Selection"]])
for(ii in datasets){
  all_cells = summary_list[["Selection"]][[ii]][["all_cells"]][["SeparateAnnotation_1"]]
  cell_type = c(cell_type, names(all_cells))
}
cell_type = unique(cell_type)
all_cells_matrix = matrix(0, nrow = length(datasets), ncol = length(cell_type), dimnames = list(datasets, cell_type))
selected_cells_matrix = matrix(0, nrow = length(datasets), ncol = length(cell_type), dimnames = list(datasets, cell_type))
selection_ratio_matrix = matrix(NA, nrow = length(datasets), ncol = length(cell_type), dimnames = list(datasets, cell_type))
for(ii in datasets){
  this_all_cells = summary_list[["Selection"]][[ii]][["all_cells"]][["SeparateAnnotation_1"]]
  this_selected_cells = summary_list[["Selection"]][[ii]][["selected_cells"]][["SeparateAnnotation_1"]]
  this_ratio = (this_selected_cells / sum(this_selected_cells)) / (this_all_cells / sum(this_all_cells))
  all_cells_matrix[ii, names(this_all_cells)] = this_all_cells
  selected_cells_matrix[ii, names(this_selected_cells)] = this_selected_cells
  selection_ratio_matrix[ii, names(this_ratio)] = this_ratio
}
py$all_cells_matrix = as.data.frame(t(all_cells_matrix))
py$selected_cells_matrix = as.data.frame(t(selected_cells_matrix))
py$selection_ratio_matrix = as.data.frame(t(selection_ratio_matrix))
py_run_string("distribution_matrix(all_cells_matrix, path=figure_dir + '/Distribution_All_Cells.pdf')")
py_run_string("distribution_matrix(all_cells_matrix, path=figure_dir + '/Distribution_All_Cells.pdf')")
py_run_string("distribution_matrix(selected_cells_matrix, path=figure_dir + '/Distribution_Selected_Cells.pdf')")
py_run_string("distribution_matrix(selection_ratio_matrix, path=figure_dir + '/Distribution_Selection_Ratio.pdf')")
py_run_string("del all_cells_matrix")
py_run_string("del selected_cells_matrix")
py_run_string("del selection_ratio_matrix")
write.table(all_cells_matrix, paste0(py$output_dir, "/Selection_all_cells.tsv"), sep = "\t", row.names = T, col.names = T, quote = F)
write.table(selected_cells_matrix, paste0(py$output_dir, "/Selection_Select_cells.tsv"), sep = "\t", row.names = T, col.names = T, quote = F)
write.table(selection_ratio_matrix, paste0(py$output_dir, "/Selection_Ratio.tsv"), sep = "\t", row.names = T, col.names = T, quote = F)
rm(all_cells_matrix, selected_cells_matrix, selection_ratio_matrix)
######
# Adata of all cells
py_run_string("adata = ad.concat({
              'SS3_HumanSkin_20K': adata_2,
              'Ning': adata_3,
              'Gaydosik_Fuschiotti_ClinicalCancerResearch_2019': adata_4,
              'Wang_Atwood_NatureCommunications_2020': adata_5, 
              'Takahashi_Lowry_JournalofInvestigativeDermatology_2020': adata_6, 
              'Belote_Torres_NatureCellBiology_2021': adata_7,
              'Boldo_Lyko_CommunicationsBiology_2020': adata_8,
              'Cheng_Cho_CellReports_2018': adata_9,
              'Ji_Khavari_Cell_2020': adata_10,
              'Tabib_Lafyatis_NatureCommunications_2021': adata_11,
              'He_Guttman_JournalofAllergyandClinicalImmunology_2020': adata_12,
              'Zou_Liu_DevelopmentalCell_2020': adata_13,
              'Wiedemann_Andersen_CellReports_2023': adata_14,
              'Alkon_Stingl_JournalofAllergyandClinicalImmunology_2022': adata_15,
              'Dunlap_Rao_JCIInsight_2022': adata_16,
              'Boothby_Rosenblum_Nature_2021': adata_17,
              'Gao_Hu_CellDeathDisease_2021': adata_18,
              'Hughes_Shalek_Immunity_2020': adata_19,
              'Kim_Nagao_NatureMedicine_2020': adata_20,
              'Rindler_Brunner_MolecularCancer_2021': adata_21,
              'Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020': adata_22,
              'Theocharidis_Bhasin_NatureCommunications_2022': adata_23,
              'Vorstandlechner_Mildner_NatureCommunications_2021': adata_24,
              'Xu_Chen_Nature_2022': adata_25,
              }, axis=0, join='outer', index_unique='---')")
#py_run_string(paste0("adata = ad.read_h5ad('", paste0(py$output_dir, "/adata_merged.h5ad"), "')"))
PreviousIntegratedAnnotation_2 = readRDS("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v65_gauss/PreviousIntegratedAnnotation_2.rds")
IntegratedAnnotation_2 = readRDS("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v65_gauss/IntegratedAnnotation_2.rds")
IntegratedAnnotation_3_Keratinocyte = readRDS("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v65_gauss/IntegratedAnnotation_3_Keratinocyte.rds")
py$adata$obs[["IntegratedAnnotation_2_v35"]] = PreviousIntegratedAnnotation_2[py$adata$obs_names$values]
py$adata$obs[["IntegratedAnnotation_2_v65"]] = IntegratedAnnotation_2[py$adata$obs_names$values]
py$adata$obs[["IntegratedAnnotation_3_v65"]] = IntegratedAnnotation_3_Keratinocyte[py$adata$obs_names$values]
py_run_string(paste0("adata_obs_to_string(adata).write('", paste0(py$output_dir, "/adata_merged.h5ad"), "')"))
######
py_run_string("del adata_2")
py_run_string("del adata_3")
py_run_string("del adata_4")
py_run_string("del adata_5")
py_run_string("del adata_6")
py_run_string("del adata_7")
py_run_string("del adata_8")
py_run_string("del adata_9")
py_run_string("del adata_10")
py_run_string("del adata_11")
py_run_string("del adata_12")
py_run_string("del adata_13")
py_run_string("del adata_14")
py_run_string("del adata_15")
py_run_string("del adata_16")
py_run_string("del adata_17")
py_run_string("del adata_18")
py_run_string("del adata_19")
py_run_string("del adata_20")
py_run_string("del adata_21")
py_run_string("del adata_22")
py_run_string("del adata_23")
py_run_string("del adata_24")
py_run_string("del adata_25")
py_run_string("del adata_2_hvg")
py_run_string("del adata_3_hvg")
py_run_string("del adata_4_hvg")
py_run_string("del adata_5_hvg")
py_run_string("del adata_6_hvg")
py_run_string("del adata_7_hvg")
py_run_string("del adata_8_hvg")
py_run_string("del adata_9_hvg")
py_run_string("del adata_10_hvg")
py_run_string("del adata_11_hvg")
py_run_string("del adata_12_hvg")
py_run_string("del adata_13_hvg")
py_run_string("del adata_14_hvg")
py_run_string("del adata_15_hvg")
py_run_string("del adata_16_hvg")
py_run_string("del adata_17_hvg")
py_run_string("del adata_18_hvg")
py_run_string("del adata_19_hvg")
py_run_string("del adata_20_hvg")
py_run_string("del adata_21_hvg")
py_run_string("del adata_22_hvg")
py_run_string("del adata_23_hvg")
py_run_string("del adata_24_hvg")
py_run_string("del adata_25_hvg")
######