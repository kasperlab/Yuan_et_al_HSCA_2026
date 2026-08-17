##########
# HF&Gland
##########
this_celltype = "HF_Gland"
this_celltype_name = gsub("_", "", gsub("[ -]", "", this_celltype), fixed = T)
py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs_Karl_20250603")
dir.create(py$this_output_dir, showWarnings = F, recursive = T)
py_run_string("this_figure_dir = this_output_dir + '/figures'")
py_run_string("sc.settings.figdir = this_figure_dir")
qc_state = "before_QC"
max_gene_number = 200
######
adata_str = paste0("adata_", this_celltype_name, "_imbalance")
clustering_embedding = "PCA_use"
visualization_embedding = "X_umap_3_imbalance"
######
if(!exists("annotation_1", envir = .GlobalEnv)){
  py$annotation_1 = readRDS(paste0(py$output_dir, "/annotation_1.rds"))
}
clustering_keys = readRDS(paste0(py$this_output_dir, "/clustering_keys.rds"))
tryCatch({
  stop()
  py_run_string(paste0("print(", adata_str, ")"))
}, error = function(x){
  py_run_string(paste0(adata_str, " = ad.read_h5ad(this_output_dir + '/clustering.h5ad')"))
  py_run_string(paste0("del ", adata_str, ".uns"))
  print(paste0("Loaded ", this_celltype_name, " from cache."))
})
######
py_run_string("annotation_tmp = copy.deepcopy(annotation_1)")
cell_number = eval(parse(text = paste0("py$", adata_str, "$n_obs")))
py$dotsize = max(c(8, get_dot_size(cell_number)))
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
py_run_string(paste0("", adata_str, ".obs['", bodysite_slot, "'] = tmp_this_factor.values"))
py_run_string("del tmp_this_factor")
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", bodysite_slot, "'], size=dotsize, ncols=1, legend_loc='right margin', legend_fontsize=6, frameon=False, save='_SS3_HumanSkin_20K_AnatomicalRegionLevel3.pdf')"))
py_run_string(paste0("clustering_plot(", adata_str, ", '", bodysite_slot, "', basis='", visualization_embedding, "', size=dotsize, ncols=5, save='_detail_", bodysite_slot, ".pdf')"))
density_visualization(adata_str, visualization_embedding, bodysite_slot, ".pdf", "density_color", "5", "8")
######
py_run_string(paste0("fig = sc.pl.embedding(", adata_str, ",
                                            basis='", visualization_embedding, "',
                                            legend_loc='right margin',
                                            legend_fontoutline=0.5,
                                            legend_fontsize=6,
                                            show=False,
                                            ncols=5,
                                            size=dotsize,
                                            frameon=False,
                                            add_outline=False)"))
py_run_string(paste0("sc.pl.embedding(", adata_str, "[adata_HFGland_imbalance.obs['PCA_use_leiden_200'] == '27'], '", visualization_embedding, "', color=['AnatomicalRegionLevel2'], legend_loc='right margin', legend_fontsize=6, legend_fontoutline=0.5, size=dotsize, ncols=6, frameon=False, ax=fig, save='_", adata_str, "_leiden_200_C27_AnatomicalRegionLevel2.pdf')"))
# DEG
DEG_key = "PCA_use_leiden_300"
for(used_table in c("top_p_by_score", "top_p_by_FC", "top_p_by_P")){
  this_path = paste0(py$this_output_dir, "/DEG_", adata_str, "_" , DEG_key, ".rds")
  if(file.exists(this_path)){
    DEG_result = readRDS(this_path)
  }else{
    DEG_result = DEG_Analysis(eval(parse(text = paste0("py$", adata_str))), DEG_key, py$this_output_dir, prefix = paste0(adata_str, "_" , DEG_key, "_"), top_DE_range = top_DE_range)
    saveRDS(DEG_result, this_path)
  }
  this_table = DEG_result[["rest"]][[used_table]]
  py_run_string(paste0("sc.tl.dendrogram(", adata_str, ", groupby='", DEG_key, "', use_rep='", clustering_embedding, "')"))
  this_clusters = unique(this_table$Cluster)
  top_DEGs = c()
  py_run_string("DEG_dict = {}")
  for(this_cluster in this_clusters){
    top_DEGs = c(top_DEGs, head(this_table[this_table$Cluster == this_cluster, ], 10)$Gene)
    py_run_string(paste0("DEG_dict['", this_cluster, "'] = ['", paste(head(this_table[this_table$Cluster == this_cluster, ], 10)$Gene, collapse = "', '"), "']"))
  }
  top_DEGs_ordered = c()
  for(this_cluster in eval(parse(text = paste0("py$", adata_str, "$uns[['dendrogram_", DEG_key, "']][['categories_ordered']]")))){
    top_DEGs_ordered = c(top_DEGs_ordered, head(this_table[this_table$Cluster == this_cluster, ], 10)$Gene)
  }
  if(length(top_DEGs) <= max_gene_number){
    #py_run_string(paste0("sc.pl.heatmap(", adata_str, ", var_names=DEG_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=True, show_gene_labels=True, save='_", adata_str, "_" , DEG_key, "_ordered.pdf')"))
    py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=DEG_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
    py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", used_table, ".pdf', bbox_inches = 'tight')"))
  }else{
    for(ii in 1:ceiling(length(top_DEGs) / max_gene_number)){
      this_index = seq(max_gene_number) + max_gene_number * (ii - 1)
      this_index = this_index[this_index <= length(top_DEGs)]
      py$tmp_genes = top_DEGs[this_index]
      py$tmp_genes_ordered = top_DEGs_ordered[this_index]
      #py_run_string(paste0("sc.pl.heatmap(", adata_str, ", var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, show_gene_labels=True, save='_", adata_str, "_" , DEG_key, "_ordered_", ii, ".pdf')"))
      py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
      py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", ii, "_", used_table, ".pdf', bbox_inches = 'tight')"))
      py_run_string("del tmp_genes")
    }
  }
}
IntegratedAnnotation_3_Keratinocyte = readRDS("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss/IntegratedAnnotation_3_Keratinocyte.rds")
py$adata_HFGland_imbalance$obs[["IntegratedAnnotation_3_v67"]] = IntegratedAnnotation_3_Keratinocyte[py$adata_HFGland_imbalance$obs_names$values]
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", paste(clustering_keys, collapse = "', '"), "'] + ['IntegratedAnnotation_2', 'IntegratedAnnotation_3_v67', 'IntegratedAnnotation_3_v65', 'karl_clustering', 'karl_annotation', 'AnatomicalRegionLevel2', 'StudyID'], legend_loc='right margin', legend_fontsize=6, legend_fontoutline=0.5, size=dotsize, ncols=1, frameon=False, save='_", adata_str, "_clustering_", clustering_embedding, ".pdf')"))
ii = "IntegratedAnnotation_3_v67"
py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=5, save='_", adata_str, "_detail_", ii, ".pdf')"))
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
py_run_string(paste0("", adata_str, ".obs['", bodysite_slot, "'] = tmp_this_factor.values"))
py_run_string("del tmp_this_factor")
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", bodysite_slot, "'], size=dotsize, ncols=1, legend_loc='right margin', legend_fontsize=6, frameon=False, save='_SS3_HumanSkin_20K_AnatomicalRegionLevel3.pdf')"))
py_run_string(paste0("clustering_plot(", adata_str, ", '", bodysite_slot, "', basis='", visualization_embedding, "', size=dotsize, ncols=2, save='_detail_", bodysite_slot, ".pdf')"))
density_visualization(adata_str, visualization_embedding, bodysite_slot, ".pdf", "density_color", "2", as.character(get_dot_size(cell_number)))
######
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=np.setdiff1d(factors, factors_only_binary).tolist() + ['phase'], size=dotsize, ncols=1, legend_loc='right margin', legend_fontsize=6, frameon=False, save='_factors.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['", paste(clustering_keys, collapse = "', '"), "'] + annotation_tmp, legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=8, frameon=False, save='_clustering.pdf')"))
for(ii in c(setdiff(py$factors, c("SampleID", "DonorID_short", "DonorID")), "phase", clustering_keys)){
  py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=7, save='_detail_", ii, ".pdf')"))
}
for(ii in unlist(py$annotation_tmp)){
  py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=7, save='_detail_", ii, ".pdf')"))
}
density_visualization(adata_str, visualization_embedding, "StudyID", ".pdf", "density_color", "6", as.character(py$dotsize))
density_visualization(adata_str, visualization_embedding, "Sex", ".pdf", "density_color", "2", as.character(get_dot_size(cell_number)))
density_visualization(adata_str, visualization_embedding, "AnatomicalRegionLevel1", ".pdf", "density_color", "4", as.character(get_dot_size(cell_number)))
density_visualization(adata_str, visualization_embedding, "AnatomicalRegionLevel2", ".pdf", "density_color", "6", as.character(get_dot_size(cell_number)))
density_visualization(adata_str, visualization_embedding, "AnatomicalRegionLevel3", ".pdf", "density_color", "7", as.character(get_dot_size(cell_number)))
# Signatures
py$Keratinocyte_markers = markers_list[["Keratinocyte"]][markers_list[["Keratinocyte"]] %in% eval(parse(text = paste0("py$", adata_str, "")))$var_names$values]
py$Fibroblast_MuralCell_markers = markers_list[["Fibroblast_MuralCell"]][markers_list[["Fibroblast_MuralCell"]] %in% eval(parse(text = paste0("py$", adata_str, "")))$var_names$values]
py$NeuralCrestderivedCell_markers = markers_list[["NeuralCrestderivedCell"]][markers_list[["NeuralCrestderivedCell"]] %in% eval(parse(text = paste0("py$", adata_str, "")))$var_names$values]
py$EndothelialCell_markers = markers_list[["EndothelialCell"]][markers_list[["EndothelialCell"]] %in% eval(parse(text = paste0("py$", adata_str, "")))$var_names$values]
py$ImmuneCell_markers = markers_list[["ImmuneCell"]][markers_list[["ImmuneCell"]] %in% eval(parse(text = paste0("py$", adata_str, "")))$var_names$values]
py$Plasma_Erythrocyte_markers = markers_list[["Plasma_Erythrocyte"]][markers_list[["Plasma_Erythrocyte"]] %in% eval(parse(text = paste0("py$", adata_str, "")))$var_names$values]
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=Keratinocyte_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=5, frameon=False, save='_Signatures_Keratinocyte.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=Fibroblast_MuralCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=4, frameon=False, save='_Signatures_Fibroblast_MuralCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=NeuralCrestderivedCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=4, frameon=False, save='_Signatures_NeuralCrestderivedCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=EndothelialCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=4, frameon=False, save='_Signatures_EndothelialCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=ImmuneCell_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=6, frameon=False, save='_Signatures_ImmuneCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=Plasma_Erythrocyte_markers, layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=4, frameon=False, save='_Signatures_Plasma_Erythrocyte.pdf')"))





