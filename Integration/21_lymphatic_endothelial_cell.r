############################
# Lymphatic Endothelial Cell
############################
this_celltype = "Lymphatic Endothelial Cell"
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
############
# Signatures
############
DEG_key = "leiden_2_100"
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
# VEGFR3: FLT4
# CD31: PECAM1
# CD73: NT5E
# ITGAB2: ITGA2B
# PDL1: CD274
py_run_string("marker_dict_ijms222111976_F2 = {
  'LEC - Cortical': ['MARCO', 'KCNJ8', 'LYVE1', 'ITIH5', 'ANXA2', 'PTX3'],
  'LEC - Ceiling': ['ACKR4', 'NT5E', 'CAV1', 'CD36'],
  'LEC - Valve': ['PROX1', 'CLDN11'],
  'LEC - Capillary': ['PDPN', 'PECAM1', 'CCL21', 'MMRN1', 'LYVE1', 'PROX1', 'FLT4'],
  'LEC - Collecting Vessel': ['PDPN', 'PECAM1', 'CCL21', 'MMRN1', 'LYVE1', 'PROX1', 'FLT4', 'NT5E', 'CAV1'],
  'LEC - Floor': ['GLYCAM1', 'TNFRSF9', 'MADCAM1', 'IFNGR1', 'IFNGR2', 'CCL20', 'IL33', 'ITGA2B', 'CD274', 'CD44', 'CXCL1'],
  'LEC - Medullary': ['MARCO', 'CLEC4G', 'LYVE1', 'CD209', 'IFNGR1', 'IFNGR2', 'CXCL1', 'CLEC4M'],
  'LEC - Medullary Sinus': ['CXCL2', 'CXCL3'],
  'LEC - Medulla Ceiling': ['MFAP4']
}")
# 10.3390/ijms222111976 Figure 2
py_run_string(paste0("sc.pl.dotplot(", adata_str, ", marker_dict_ijms222111976_F2, groupby='leiden_2_100', layer='normalized', dendrogram='dendrogram_leiden_2_100', save='Signatures_marker_dict_ijms222111976_F2.pdf')"))
######
# CD206: MRC1
py_run_string("marker_dict_fimmu_2023_1235812_F2 = {
  'LEC - SCS Ceiling': ['ACKR4', 'NT5E', 'CAV1', 'NTS', 'NUDT4'],
  'LEC - SCS Floor': ['CCL20', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL5', 'TNFRSF9', 'ACKR1'],
  'LEC - Medullary Capsule-lining': ['MFAP4', 'LYVE1'],
  'LEC - Paracortical Sinus': ['LYVE1', 'CCL21', 'PDPN', 'NRP2', 'FLT4', 'ITIH3', 'SPHK1', 'PTX3', 'CD36'],
  'LEC - Valve': ['CLDN11', 'GJA1', 'ESAM', 'ANGPT2', 'FOXC2', 'CLDN5', 'CD9', 'GJA4', 'CAV1'],
  'LEC - Medullar Sinus': ['LYZ', 'ACKR1', 'CXCL1', 'CXCL2', 'CXCL3', 'MRC1', 'CD209', 'CLEC4M', 'CLEC4G', 'MARCO', 'ICAM1', 'VCAM1', 'CADM3', 'CD44', 'LYVE1']
}")
# 10.3389/fimmu.2023.1235812 Figure 2
py_run_string(paste0("sc.pl.dotplot(", adata_str, ", marker_dict_fimmu_2023_1235812_F2, groupby='leiden_2_100', layer='normalized', dendrogram='dendrogram_leiden_2_100', save='Signatures_marker_dict_fimmu_2023_1235812_F2.pdf')"))
######
for(ii in unlist(py$annotation_tmp)){
  py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=7, save='_detail_", ii, ".pdf')"))
}
density_visualization(adata_str, visualization_embedding, "StudyID", ".pdf", "density_color", "6", as.character(py$dotsize))
#
cell_annotation = rep(NA, eval(parse(text = paste0("py$", adata_str, "$n_obs"))))
names(cell_annotation) = eval(parse(text = paste0("py$", adata_str, "$obs_names$values")))
#
leiden_2_100 = as.character(eval(parse(text = paste0("py$", adata_str, "")))$obs[['leiden_2_100']])
#
cell_annotation[leiden_2_100 == "0"] = "LEC - Common(FLT4/NT5E/NUDT4)"
cell_annotation[leiden_2_100 == "1"] = "LEC - Common"
cell_annotation[leiden_2_100 == "2"] = "LEC - Common"
cell_annotation[leiden_2_100 == "3"] = "LEC - CXCL1/CXCL2/CXCL3/FABP5/FKBP1A"
cell_annotation[leiden_2_100 == "4"] = "LEC - Common(FLT4/NT5E/NUDT4)"
cell_annotation[leiden_2_100 == "5"] = "LEC - Common(FLT4)"
cell_annotation[leiden_2_100 == "6"] = "LEC - Common(FLT4/NT5E/NUDT4)"
cell_annotation[leiden_2_100 == "7"] = "LEC - Lymphatic Valve Downstream (ADM/ANGPT2/SCG3/GJA4/CLDN11/PDLIM1/PROCR/RAMP3)" # 10.1038/s41577-020-0281-x Figure 1
cell_annotation[leiden_2_100 == "8"] = "LEC - Subcapsular Sinus Ceiling (NTS/PLAAT4/IL33)" # 10.3389/fimmu.2023.1235812 Figure 2 
cell_annotation[leiden_2_100 == "9"] = "LEC - SNHG5/GAS5/SNHG29/SNHG6"
cell_annotation[leiden_2_100 == "10"] = "LEC - Lymphatic Valve Downstream (ADM/ANGPT2/SCG3/GJA4/CLDN11)" # 10.1038/s41577-020-0281-x Figure 1
cell_annotation[leiden_2_100 == "11"] = "LEC - ADAMTS6/SLC41A1/NEO1/AGRN/CD200/SBSPON"
cell_annotation[leiden_2_100 == "12"] = "LEC - MIX(MYE)"
cell_annotation[leiden_2_100 == "13"] = "LEC - Common(FLT4)"
cell_annotation[leiden_2_100 == "14"] = "LEC - MIX(KC)"
cell_annotation[leiden_2_100 == "15"] = "LEC - MIX(FIB)"
cell_annotation[leiden_2_100 == "16"] = "LEC - ADAMTS6/SLC41A1/NEO1/AGRN"
cell_annotation[leiden_2_100 == "17"] = "LEC - MIX(KC)"
cell_annotation[leiden_2_100 == "18"] = "LEC - CXCL1/CXCL2/CXCL3"
cell_annotation[leiden_2_100 == "19"] = "LEC - MIX(MUR)"
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
  'Common': ['FLT4', 'NT5E', 'NUDT4'],
  'Lymphatic Valve Downstream': ['ADM', 'ANGPT2', 'SCG3', 'GJA4', 'CLDN11', 'PDLIM1', 'PROCR', 'RAMP3'],
  'Subcapsular Sinus Ceiling': ['NTS', 'PLAAT4', 'IL33'],
  'CXCL1/CXCL2/CXCL3/FABP5/FKBP1A': ['CXCL1', 'CXCL2', 'CXCL3', 'FABP5', 'FKBP1A'],
  'SNHG5/GAS5/SNHG29/SNHG6': ['SNHG5', 'GAS5', 'SNHG29', 'SNHG6'],
  'ADAMTS6/SLC41A1/NEO1/AGRN/CD200/SBSPON': ['ADAMTS6', 'SLC41A1', 'NEO1', 'AGRN', 'CD200', 'SBSPON'],
  'MIX(FIB)': ['FBLN1', 'COL1A1', 'COL1A2'],
  'MIX(MUR)': ['TAGLN', 'MYL9', 'ACTA2'],
  'MIX(MYE)': ['HLA-DRA', 'HLA-DPA1', 'CD74'],
  'MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10']
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
  "0" = "LEC - Common",
  "1" = "LEC - Common(FLT4)",
  "2" = "LEC - Common(FLT4/NT5E/NUDT4)",
  "3" = "LEC - Lymphatic Valve Downstream (ADM/ANGPT2/SCG3/GJA4/CLDN11)",
  "4" = "LEC - Lymphatic Valve Downstream (ADM/ANGPT2/SCG3/GJA4/CLDN11/PDLIM1/PROCR/RAMP3)",
  "5" = "LEC - Subcapsular Sinus Ceiling (NTS/PLAAT4/IL33)",
  "6" = "LEC - CXCL1/CXCL2/CXCL3",
  "7" = "LEC - CXCL1/CXCL2/CXCL3/FABP5/FKBP1A",
  "8" = "LEC - SNHG5/GAS5/SNHG29/SNHG6",
  "9" = "LEC - ADAMTS6/SLC41A1/NEO1/AGRN",
  "10" = "LEC - ADAMTS6/SLC41A1/NEO1/AGRN/CD200/SBSPON",
  "11" = "LEC - MIX(FIB)",
  "12" = "LEC - MIX(MUR)",
  "13" = "LEC - MIX(MYE)",
  "14" = "LEC - MIX(KC)"
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
                           '2': colorsys.hsv_to_rgb(0, 0.7, 0.4),
                           '3': colorsys.hsv_to_rgb(0.1, 0.9, 0.9),
                           '4': colorsys.hsv_to_rgb(0.1, 0.6, 0.6),
                           '5': colorsys.hsv_to_rgb(0.2, 1, 0.7),
                           '6': colorsys.hsv_to_rgb(0.3, 0.9, 0.9),
                           '7': colorsys.hsv_to_rgb(0.3, 0.7, 0.4),
                           '8': colorsys.hsv_to_rgb(0.4, 0.6, 0.6),
                           '9': colorsys.hsv_to_rgb(0.55, 0.7, 0.4),
                           '10': colorsys.hsv_to_rgb(0.7, 0.6, 0.6),
                           '11': colorsys.hsv_to_rgb(0, 0.05, 0.9),
                           '12': colorsys.hsv_to_rgb(0.25, 0.05, 0.9),
                           '13': colorsys.hsv_to_rgb(0.5, 0.05, 0.9),
                           '14': colorsys.hsv_to_rgb(0.75, 0.05, 0.9)
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








