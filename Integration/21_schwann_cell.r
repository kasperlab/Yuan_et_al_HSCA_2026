##############
# Schwann Cell
##############
this_celltype = "Schwann Cell"
this_celltype_name = gsub("_", "", gsub("[ -]", "", this_celltype), fixed = T)
py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs")
dir.create(py$this_output_dir, showWarnings = F, recursive = T)
py_run_string("this_figure_dir = this_output_dir + '/figures'")
py_run_string("sc.settings.figdir = this_figure_dir")
qc_state = "before_QC"
max_gene_number = 200
######
adata_str = paste0("adata_", this_celltype_name, "_imbalance")
clustering_embedding = "PCA_use"
visualization_embedding = "X_umap_2_imbalance"
######
if(!exists("annotation_1", envir = .GlobalEnv)){
  py$annotation_1 = readRDS(paste0(py$output_dir, "/annotation_1.rds"))
}
tryCatch({
  stop()
  py_run_string(paste0("print(", adata_str, ")"))
}, error = function(x){
  py_run_string(paste0(adata_str, " = ad.read_h5ad('", py$this_output_dir, "/clustering.h5ad')"))
  py_run_string(paste0(adata_str, ".uns['log1p'] = None"))
  print(paste0("Loaded ", this_celltype_name, " from cache."))
})
tryCatch({
  stop()
  py_run_string(paste0("print(", adata_str, "_hvg)"))
}, error = function(x){
  py_run_string(paste0(adata_str, "_hvg = ad.read_h5ad('", py$this_output_dir, "/clustering_hvg.h5ad')"))
  py_run_string(paste0(adata_str, "_hvg.uns['log1p'] = None"))
  print(paste0("Loaded ", this_celltype_name, "_hvg from cache."))
})
######
py_run_string(paste0(adata_str, ".obsm['PCA_use'] = ", adata_str, "_hvg.obsm['PCA_use'].copy()"))
py_run_string("annotation_tmp = copy.deepcopy(annotation_1)")
cell_number = eval(parse(text = paste0("py$", adata_str, "$n_obs")))
py$dotsize = max(c(8, get_dot_size(cell_number)))
######
leiden_2_150 = as.character(eval(parse(text = paste0("py$", adata_str, "")))$obs[['leiden_2_150']])
leiden_2_300 = as.character(eval(parse(text = paste0("py$", adata_str, "")))$obs[['leiden_2_300']])
#
cell_group = leiden_2_150
names(cell_group) = eval(parse(text = paste0("py$", adata_str, "$obs_names$values")))
#
cell_group[leiden_2_300 %in% c("9", "22")] = "9_or_22"
cell_group[leiden_2_150 == "9" & leiden_2_300 == "25"] = "9_and_25"
cell_group[leiden_2_150 %in% c("9", "10") & leiden_2_300 == "21"] = "9_or_10_and_21"
#
print(sum(is.na(cell_group)))
#
eval(parse(text = paste0("py$", adata_str, "$obs[['cell_group_customized']] = as.factor(cell_group)")))
#
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['cell_group_customized'], legend_loc='on data', legend_fontsize=6, legend_fontoutline=0.5, size=dotsize, frameon=False, ncols=4, save='_", adata_str, "_cell_group_customized_OD.pdf')"))
#
DEG_key = "cell_group_customized"
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
    py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=DEG_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=True, n_genes=10, swap_axes = False, highlight_genes = DEG_dict, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
    py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", used_table, ".pdf', bbox_inches = 'tight')"))
  }else{
    for(ii in 1:ceiling(length(top_DEGs) / max_gene_number)){
      this_index = seq(max_gene_number) + max_gene_number * (ii - 1)
      this_index = this_index[this_index <= length(top_DEGs)]
      py$tmp_genes = top_DEGs[this_index]
      py$tmp_genes_ordered = top_DEGs_ordered[this_index]
      #py_run_string(paste0("sc.pl.heatmap(", adata_str, ", var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, show_gene_labels=True, save='_", adata_str, "_" , DEG_key, "_ordered_", ii, ".pdf')"))
      py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, n_genes=10, swap_axes = False, highlight_genes = DEG_dict, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
      py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", ii, "_", used_table, ".pdf', bbox_inches = 'tight')"))
      py_run_string("del tmp_genes")
    }
  }
}
############
# Signatures
############
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[
              'SOX10', 'NCMAP', 'GFRA3', 'BCHE'
              ], layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ",
                     ncols=4, frameon=False, save='_Signatures_SchwannCell.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[
              'MARCKS', 'SCN7A', 'PRNP', 'DBI', 'PAK1', 'LMNA', 'CLDN14', 'ADAMTSL1', 'PMP2', 'SLIT2', 'COL23A1'
              ], layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ",
                     ncols=4, frameon=False, save='_Signatures_s41593_021_01005_1.pdf')"))
# https://www.nature.com/articles/s12276-022-00874-1
DEG_key = "cell_group_customized"
py_run_string("marker_dicts_12276_022_00874_1 = {
  'General Schwann Cell Marker': ['S100B', 'NGFR'],
  'SC - Myelinating': ['MBP', 'MPZ', 'PLP1', 'PMP22'],
  'SC - Non-myelinating': ['NCAM1', 'L1CAM', 'SCN7A'],
  'SC - Keloidal': ['NES', 'IGFBP5', 'CCN3'],
  'SC - Proliferating': ['TOP2A', 'MKI67']
}")
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_dicts_12276_022_00874_1, groupby='", DEG_key, "', dendrogram=False, layer='normalized', n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = marker_dict, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_var_marker_dicts_12276_022_00874_1.pdf', bbox_inches = 'tight')"))
######
for(ii in unlist(py$annotation_tmp)){
  py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=7, save='_detail_", ii, ".pdf')"))
}
density_visualization(adata_str, visualization_embedding, "StudyID", ".pdf", "density_color", "6", as.character(py$dotsize))
#
cell_annotation = rep(NA, eval(parse(text = paste0("py$", adata_str, "$n_obs"))))
names(cell_annotation) = eval(parse(text = paste0("py$", adata_str, "$obs_names$values")))
#
cell_annotation[cell_group == "0"] = "SC - Remak(NRXN1/SCN7A)"
cell_annotation[cell_group == "1"] = "SC - Myelinated (Stressed)"
cell_annotation[cell_group == "2"] = "SC - Remak(NRXN1/SCN7A/SOCS3/EGR1/S100A16/NGFR/CCN5/IFI27)"
cell_annotation[cell_group == "3"] = "SC - Remak(NRXN1/SCN7A)"
cell_annotation[cell_group == "4"] = "SC - Remak(NRXN1/SCN7A/XKR4/CADM2)"
cell_annotation[cell_group == "5"] = "SC - Remak(NRXN1/SCN7A/SOCS3/EGR1/S100A16/NGFR/ABL2)"
cell_annotation[cell_group == "7"] = "SC - Remak(NRXN1)"
cell_annotation[cell_group == "8"] = "SC - Remak(NRXN1/SCN7A/SOCS3/EGR1)"
cell_annotation[cell_group == "9"] = "SC - Remak(NRXN1/SCN7A/SOCS3/EGR1)"
cell_annotation[cell_group == "9_and_25"] = "SC - Remak(NRXN1/SCN7A/CFD/TMOD2/SBSPON)"
cell_annotation[cell_group == "9_or_22"] = "SC - Remak(NRXN1/SCN7A/SPARCL1)"
cell_annotation[cell_group == "9_or_10_and_21"] = "SC - Remak(NRXN1/SCN7A/CFD/TMOD2/SBSPON/CADM2)"
cell_annotation[cell_group == "10"] = "SC - Remak - MHC II (HLA-DRB1/CD74/HLA-DRB5)" # 
cell_annotation[cell_group == "11"] = "SC - Remak(NRXN1/SCN7A/XKR4/CADM2/FN1/ANK2/NRXN3/COL28A1/SPARCL1)"
cell_annotation[cell_group == "12"] = "SC - Nociceptive(BCHE/AQP1/PALMD)"
cell_annotation[cell_group == "13"] = "SC - Remak(NRXN1/SCN7A/SOCS3/EGR1/S100A16/NGFR/ABL2/PDLIM3)"
cell_annotation[cell_group == "14"] = "SC - Remak(NRXN1/SCN7A)"
cell_annotation[cell_group == "15"] = "SC - MIX(MEL)" # DCT
cell_annotation[cell_group == "16"] = "SC - Myelinated"
cell_annotation[cell_group == "17"] = "SC - MIX(MUR)" # TAGLN, ACTA2
cell_annotation[cell_group == "18"] = "SC - Remak(NRXN1/SCN7A/XKR4/FN1/NGFR)"
cell_annotation[cell_group == "19"] = "SC - Remak(NRXN1/S100A16/NGFR)"
cell_annotation[cell_group == "20"] = "SC - MIX(KC/FIB)" # DCN, CXCL14, APOD
#
print(sum(is.na(cell_annotation)))
######
eval(parse(text = paste0("py$", adata_str, "$obs[['IntegratedAnnotation_2']] = as.factor(cell_annotation)")))
summary_list[[qc_state]][["IntegratedAnnotation_2"]] = list("AllCells" = table(cell_annotation))
py_run_string("annotation_tmp = annotation_tmp if np.isin('IntegratedAnnotation_2', annotation_tmp) else annotation_tmp + ['IntegratedAnnotation_2']")
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=annotation_tmp, legend_loc='on data', legend_fontsize=6, legend_fontoutline=0.5, size=dotsize, frameon=False, ncols=4, save='_", adata_str, "_annotation_OD.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=annotation_tmp, legend_loc='right margin', legend_fontsize=6, legend_fontoutline=0.5, size=dotsize, frameon=False, ncols=1, save='_", adata_str, "_annotation.pdf')"))
######
ii = "IntegratedAnnotation_2"
py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=7, save='_detail_", ii, ".pdf')"))
#########
# Markers
#########
# DEG
DEG_key = "IntegratedAnnotation_2"
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
    py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=DEG_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=True, n_genes=10, swap_axes = False, highlight_genes = DEG_dict, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
    py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", used_table, ".pdf', bbox_inches = 'tight')"))
  }else{
    for(ii in 1:ceiling(length(top_DEGs) / max_gene_number)){
      this_index = seq(max_gene_number) + max_gene_number * (ii - 1)
      this_index = this_index[this_index <= length(top_DEGs)]
      py$tmp_genes = top_DEGs[this_index]
      py$tmp_genes_ordered = top_DEGs_ordered[this_index]
      #py_run_string(paste0("sc.pl.heatmap(", adata_str, ", var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, show_gene_labels=True, save='_", adata_str, "_" , DEG_key, "_ordered_", ii, ".pdf')"))
      py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, n_genes=10, swap_axes = False, highlight_genes = DEG_dict, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
      py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", ii, "_", used_table, ".pdf', bbox_inches = 'tight')"))
      py_run_string("del tmp_genes")
    }
  }
}
py_run_string("marker_dict = {
  'Schwann Cell': ['S100B'],
  'SC - Remak': ['NRXN1', 'SCN7A'],
  'SC - Remak Subtypes (SOCS3-/EGR1-/S100A16-)': ['XKR4', 'FN1', 'ANK2', 'NRXN3', 'COL28A1', 'SPARCL1', 'CFD', 'TMOD2', 'SBSPON', 'CADM2'],
  'SC - Remak - MHC II': ['HLA-DRB1', 'CD74', 'HLA-DRB5'],
  'SOCS3/EGR1/S100A16': ['SOCS3', 'EGR1', 'S100A16'],
  'SC - Remak Subtypes (SOCS3+/EGR1+)': ['CCN5', 'IFI27', 'ABL2', 'PDLIM3'],
  'SC - Remak Repair': ['NGFR'],
  'SC - Nociceptive': ['BCHE', 'AQP1', 'PALMD'],
  'SC - Nociceptive (Extend)': ['PRSS23', 'PTN', 'IGFBP5', 'PTPRZ1'],
  'SC - Myelinated': ['PRX', 'GLDN', 'MLIP', 'NCMAP', 'BCAS1'],
  'SC - Myelinated (Stressed)': ['JUNB', 'ATF3', 'FOS', 'FOSB', 'HSPA1B', 'ZFP36'],
  'MPZ/SOX': ['MPZ', 'SOX2', 'SOX10'],
  'SC - MIX(MUR)': ['TAGLN', 'ACTA2'],
  'SC - MIX(KC/FIB)': ['DCN', 'CXCL14', 'APOD'],
  'SC - MIX(MEL)': ['DCT', 'KIT', 'ANXA2P2', 'QPCT']
}")
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=True, n_genes=10, standard_scale='var', swap_axes=False, highlight_genes=None, highlight_params={'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line=False, cluster_line_highlight_params={'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_var.pdf', bbox_inches='tight')"))
#######
# Write
#######
h5ad_path = paste0(py$this_output_dir, "/IntegratedAnnotation_2.h5ad")
py_run_string(paste0("", adata_str, ".write('", h5ad_path, "')"))
#
# Further analysis
#
Annotation_Level_2 = as.character(eval(parse(text = paste0("py$", adata_str, "$obs[['IntegratedAnnotation_2']]"))))
print(table(Annotation_Level_2))
#
py$index_name_dict = list(
  "0" = "SC - Remak(NRXN1)",
  "1" = "SC - Remak(NRXN1/SCN7A)",
  "2" = "SC - Remak(NRXN1/SCN7A/XKR4/CADM2)",
  "3" = "SC - Remak(NRXN1/SCN7A/XKR4/FN1/NGFR)",
  "4" = "SC - Remak(NRXN1/SCN7A/XKR4/CADM2/FN1/ANK2/NRXN3/COL28A1/SPARCL1)",
  "5" = "SC - Remak(NRXN1/SCN7A/SPARCL1)",
  "6" = "SC - Remak(NRXN1/SCN7A/CFD/TMOD2/SBSPON)",
  "7" = "SC - Remak(NRXN1/SCN7A/CFD/TMOD2/SBSPON/CADM2)",
  "8" = "SC - Remak - MHC II (HLA-DRB1/CD74/HLA-DRB5)",
  "9" = "SC - Remak(NRXN1/SCN7A/SOCS3/EGR1)",
  "10" = "SC - Remak(NRXN1/SCN7A/SOCS3/EGR1/S100A16/NGFR/CCN5/IFI27)",
  "11" = "SC - Remak(NRXN1/SCN7A/SOCS3/EGR1/S100A16/NGFR/ABL2)",
  "12" = "SC - Remak(NRXN1/SCN7A/SOCS3/EGR1/S100A16/NGFR/ABL2/PDLIM3)",
  "13" = "SC - Remak(NRXN1/S100A16/NGFR)",
  "14" = "SC - Nociceptive(BCHE/AQP1/PALMD)",
  "15" = "SC - Myelinated",
  "16" = "SC - Myelinated (Stressed)",
  "17" = "SC - MIX(MUR)",
  "18" = "SC - MIX(KC/FIB)",
  "19" = "SC - MIX(MEL)"
)
py_run_string(paste0(adata_str, ".obs['Annotation_Level_2_renumber'] = ", adata_str, ".obs['IntegratedAnnotation_2'].map({name: f'{index}' for index, name in index_name_dict.items()})"))
py_run_string(paste0(adata_str, ".obs['Annotation_Level_2_FULL'] = ", adata_str, ".obs['IntegratedAnnotation_2'].map({name: f'{index}: {name}' for index, name in index_name_dict.items()})"))
py_run_string("new_order = [f'{index}: {name}' for index, name in index_name_dict.items()]")
#
# Color type:
# 1, 0.7
# 0.9, 0.9
# 0.7, 0.4
# 0.6, 0.6
# 0.4, 1
#
py_run_string("color_dict={'0': colorsys.hsv_to_rgb(0.9, 1, 0.7),
                           '1': colorsys.hsv_to_rgb(0.9, 0.9, 0.9),
                           '2': colorsys.hsv_to_rgb(0.9, 0.4, 1),
                           '3': colorsys.hsv_to_rgb(0, 0.9, 0.9),
                           '4': colorsys.hsv_to_rgb(0, 0.6, 0.6),
                           '5': colorsys.hsv_to_rgb(0, 1, 0.7),
                           '6': colorsys.hsv_to_rgb(0.1, 1, 0.7),
                           '7': colorsys.hsv_to_rgb(0.1, 0.4, 1),
                           '8': colorsys.hsv_to_rgb(0.1, 0.9, 0.9),
                           '9': colorsys.hsv_to_rgb(0.2, 1, 0.7),
                           '10': colorsys.hsv_to_rgb(0.25, 0.9, 0.9),
                           '11': colorsys.hsv_to_rgb(0.3, 0.6, 0.6),
                           '12': colorsys.hsv_to_rgb(0.35, 0.4, 1),
                           '13': colorsys.hsv_to_rgb(0.4, 1, 0.7),
                           '14': colorsys.hsv_to_rgb(0.5, 0.9, 0.9),
                           '15': colorsys.hsv_to_rgb(0.6, 0.4, 1),
                           '16': colorsys.hsv_to_rgb(0.7, 0.9, 0.9),
                           '17': colorsys.hsv_to_rgb(0, 0.05, 0.9),
                           '18': colorsys.hsv_to_rgb(0.3, 0.05, 0.9),
                           '19': colorsys.hsv_to_rgb(0.6, 0.05, 0.9)
                          }")
py_run_string("color_dict_FULL = {f'{index}: {name}': color_dict[index] for index, name in index_name_dict.items()}")
#
py_run_string(paste0(adata_str, ".obs['Annotation_Level_2_FULL'] = ", adata_str, ".obs['Annotation_Level_2_FULL'].cat.reorder_categories(new_order)"))
py_run_string(paste0("print(", adata_str, ".obs['Annotation_Level_2_FULL'].cat.categories)"))
#
py_run_string(paste0("tmp_adata = ad.AnnData(obs=", adata_str, ".obs.copy(), obsm={'", visualization_embedding, "': ", adata_str, ".obsm['", visualization_embedding, "'].copy()})"))
py_run_string("np.random.seed(0)")
py_run_string("tmp_adata_random = tmp_adata[np.random.permutation(list(range(tmp_adata.n_obs))), :]")
py_run_string("del tmp_adata")
py_run_string("gc.collect()")
py_run_string(paste0("sc.pl.embedding(tmp_adata_random,
                                      basis='", visualization_embedding, "',
                                      color=['Annotation_Level_2_FULL'],
                                      legend_loc='none',
                                      legend_fontoutline=0.5,
                                      legend_fontsize=6,
                                      show=False,
                                      ncols=5,
                                      size=8,
                                      frameon=False,
                                      palette=color_dict_FULL,
                                      add_outline=False,
                                      save='_Annotation_Level_2_NL.png'
)"))
py_run_string(paste0("fig = sc.pl.embedding(tmp_adata_random,
                                            basis='", visualization_embedding, "',
                                            color=['Annotation_Level_2_renumber'],
                                            legend_loc='on data',
                                            legend_fontoutline=0.5,
                                            legend_fontsize=6,
                                            show=False,
                                            ncols=1,
                                            size=0,
                                            frameon=False,
                                            add_outline=False,
                                            return_fig=True
)"))
py_run_string(paste0("sc.pl.embedding(tmp_adata_random,
                                      basis='", visualization_embedding, "',
                                      color=['Annotation_Level_2_FULL'],
                                      legend_loc='right margin',
                                      legend_fontsize=6,
                                      show=False,
                                      size=8,
                                      add_outline=False,
                                      palette=color_dict_FULL,
                                      ax=fig.axes[0],
                                      save='_Annotation_Level_2_RM.pdf'
)"))
py_run_string("del tmp_adata_random")
py_run_string("gc.collect()")
#
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict.values() for x in y], layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=6, frameon=False, save='_marker_genes.pdf')"))
#
DEG_key = "Annotation_Level_2_FULL"
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_var.pdf', bbox_inches = 'tight')"))
#
ii = "Annotation_Level_2_FULL"
py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', mask=", adata_str, ".obs['StudyID'] == 'SS3_HumanSkin_20K', size=dotsize, ncols=4, save='_detail_", ii, "_SS3_HumanSkin_20K.pdf')"))
py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=4, save='_detail_", ii, ".pdf')"))
#
py_run_string(paste0("print(", adata_str, "[", adata_str, ".obs['StudyID'] == 'SS3_HumanSkin_20K'].obs['Annotation_Level_2_FULL'].value_counts())"))
#
py$sankey_categories = c("StudyID", "Annotation_Level_2_FULL", "AnatomicalRegionLevel2")
py_run_string(paste0("sankey_fig = plot_sankey_diagram(", adata_str, ", categories=sankey_categories, return_fig=True)"))
py_run_string(paste0("sankey_fig.write_image(f'{this_figure_dir}/sankey_", adata_str, "_", paste(py$sankey_categories, collapse = "_"), ".pdf')"))
py$sankey_categories = c("StudyID", "Annotation_Level_2_FULL", "AnatomicalRegionLevel2", "LibraryPlatform")
py_run_string(paste0("sankey_fig = plot_sankey_diagram(", adata_str, ", categories=sankey_categories, return_fig=True)"))
py_run_string(paste0("sankey_fig.write_image(f'{this_figure_dir}/sankey_", adata_str, "_", paste(py$sankey_categories, collapse = "_"), ".pdf')"))
py$sankey_categories = c("StudyID", "Annotation_Level_2_FULL", "LibraryPlatform")
py_run_string(paste0("sankey_fig = plot_sankey_diagram(", adata_str, ", categories=sankey_categories, return_fig=True)"))
py_run_string(paste0("sankey_fig.write_image(f'{this_figure_dir}/sankey_", adata_str, "_", paste(py$sankey_categories, collapse = "_"), ".pdf')"))
#######
# Write
#######
h5ad_path = paste0(py$this_output_dir, "/IntegratedAnnotation_2_20250531.h5ad")
py_run_string(paste0("", adata_str, ".write('", h5ad_path, "')"))








# 20251002
py_run_string(paste0(adata_str, " = ad.read_h5ad('", py$this_output_dir, "/IntegratedAnnotation_2_20250531.h5ad')"))

DEG_key = "Annotation_Level_2_renumber"
clustering_embedding = "PCA_use"
py_run_string(paste0("print(", adata_str, ".obs['", DEG_key, "'].cat.categories.size)"))
py_run_string(paste0("sc.tl.dendrogram(", adata_str, ", groupby='", DEG_key, "', use_rep='", clustering_embedding, "')"))


py_run_string(paste0("fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(", adata_str, ".obs['", DEG_key, "'].cat.categories.size * 0.3, 20))"))
py_run_string(paste0("sc.pl.dendrogram(", adata_str, ", '", DEG_key, "', ax=axs)"))
py_run_string(paste0("fig.savefig(f'{this_figure_dir}/DendrogramPlot_", adata_str, "_", DEG_key, ".pdf')"))
py_run_string("plt.close('all')")
py_run_string("gc.collect()")







