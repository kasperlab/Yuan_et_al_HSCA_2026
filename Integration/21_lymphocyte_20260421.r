############
# Lymphocyte
############
this_celltype = "Lymphocyte"
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
############
# Signatures
############
DEG_key = "leiden_2_200"
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
  #########
  # Markers
  #########
  py_run_string(paste0("sc.pl.dotplot(", adata_str, ", var_names=['IFNG', 'IL2', 'IL9', 'IL4', 'IL5', 'IL13', 'IL17A', 'IL17F', 'IL22', 'IL21'], groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_T_helper.pdf')"))
  py_run_string(paste0("sc.pl.dotplot(", adata_str, ", var_names=['NKG7', 'CD8A', 'CD8B', 'GZMK', 'IL7R', 'XCL1', 'XCL2', 'KIT'], groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_NK_Tc_ILC.pdf')"))
}


this_table = DEG_leiden_list[[this_celltype_name]][[DEG_key]][["rest"]][["table"]]
for(this_cluster in unique(this_table$Cluster)){
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", 'umap_2_imbalance', color=[
            '", paste(head(this_table[this_table$Cluster == this_cluster, ], top_DE_range)$Gene, collapse = "', '"), "'
            ], layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=10, frameon=False, save='_Signatures_", DEG_key, "_", this_cluster, ".pdf')"))
}

py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=[
  'CD79A', 'JCHAIN', 'IGKC', 'CD3D', 'CD3G', 'CD5', 'CD6', 'CD8A', 
  'CD8B', 'CD40LG', 'CD4', 'PTGER4', 'BATF', 'FOXP3', 'TIGIT', 'CTLA4',
  'KLRB1', 'XCL1', 'XCL2', 'TNFRSF18', 'TNFSF11', 'FCER1G', 'KLRC1', 'PIK3R1',
  'KLRD1', 'GNLY', 'PRF1', 'GZMB', 'FCGR3A', 'IL7R', 'IL1RL1', 'PTGDR2',
  'GATA3', 'TNFRSF25', 'KIT', 'RORC', 'NKG7', 'GZMA', 'GZMH', 'GZMK', 'GZMM'
], groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='LC.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_2_imbalance', color=[
  'CD79A', 'JCHAIN', 'IGKC', 'CD3D', 'CD3G', 'CD5', 'CD6', 'CD8A', 
  'CD8B', 'CD40LG', 'CD4', 'PTGER4', 'BATF', 'FOXP3', 'TIGIT', 'CTLA4',
  'KLRB1', 'XCL1', 'XCL2', 'TNFRSF18', 'TNFSF11', 'FCER1G', 'KLRC1', 'PIK3R1',
  'KLRD1', 'GNLY', 'PRF1', 'GZMB', 'FCGR3A', 'IL7R', 'IL1RL1', 'PTGDR2',
  'GATA3', 'TNFRSF25', 'KIT', 'RORC', 'NKG7', 'GZMA', 'GZMH', 'GZMK', 'GZMM'
], layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=10, frameon=False, save='_LC.pdf')"))

py_run_string("marker_dict_JI = {
  'Epithelial': ['KRT5', 'KRT14', 'KRT1', 'KRT10'],
  'HF': ['KRT15', 'LHX2', 'SOX9'],
  'Fibroblast': ['COL1A1', 'COL1A2', 'LUM'],
  'Melanocyte': ['MLANA', 'DCT', 'PMEL'],
  'Endothelial Cell': ['TFF3', 'CLDN5', 'VWF'],
  'B/Plasma': ['IGLL5', 'JCHAIN', 'MS4A1', 'CD79A'],
  'Myeloid Cell': ['LYZ', 'HLA-DRB1', 'HLA-DRA', 'HLA-DQB2'],
  'Monocyte (Myeloid Cell Sub)': ['CD14'],
  'Macrophage (Monocyte Sub)': ['CD163', 'CD68'],
  'MDSC': ['S100A8', 'S100A9', 'TREM1'],
  'CD1C DC (Non-monocyte)': ['CD1C', 'CLEC10A'],
  'CLEC9A (Non-monocyte)': ['CLEC9A', 'CADM1', 'XCR1'],
  'LC (Non-monocyte)': ['CD207', 'CD1A', 'S100B'],
  'AS DC (Non-monocyte)': ['AXL', 'SIGLEC6', 'IGFBP5', 'PPP1R14A'],
  'plasmacytoid DC': ['CLEC4C', 'IL3RA'],
  'migrating DC': ['CCR7', 'CCL19'],
  'plasmacytoid DC': ['CLEC4C', 'IL3RA'],
  'T Cell': ['CD3D', 'CD2', 'CD7'],
  'Pan-T Cell': ['CD3D', 'CD3G'],
  'Naive': ['CCR7'],
  'CD4': ['CD4'],
  'Treg': ['FOXP3'],
  'Cycling': ['MKI67'],
  'CD8': ['CD8A'],
  'NK': ['XCL1'],
  'CD4+ Naive': ['CCR7', 'SELL', 'IL7R'],
  'CD8+ Naive': ['ANXA1', 'NBAS', 'CLU', 'LMNA'],
  'CD4+ RGCC': ['RGCC', 'PLIN2', 'ZNF331', 'CREM', 'SLC2A3', 'KLRB1', 'ZFP36L2'],
  'CD8+ TEM': ['GZMK', 'CCL5', 'CD8B', 'HSPA1B', 'DNAJB1', 'DUSP2', 'CMC1', 'SAMD3', 'GIMAP7', 'GIMAP4'],
  'CD8+ TEMRA': ['FGFBP2', 'GZMH', 'KLRG1', 'NKG7', 'PLAC8', 'FCGR3A', 'PLEK'],
  'CD8+ Exhausted': ['CCL4', 'PRF1', 'GZMB', 'CD8A', 'GZMA', 'AC092580.4', 'KLRC2', 'KLRC1', 'KRT86', 'RGS1', 'CCL3'],
  'CD4+ Pre-Exh': ['IGFL2', 'NMB', 'CXCL13', 'FKBP5', 'CHN1', 'TSHZ2', 'BCAS3', 'NR3C1', 'GRAMD1A'],
  'CD4+ Exhausted': ['LAG3', 'PDCD1', 'MAL', 'CSF2', 'IFNG', 'GADD45G', 'IFI6'],
  'Treg': ['FOXP3', 'IL2RA', 'CD27', 'TNFRSF4', 'BATF', 'AC017002.1', 'CTLA4', 'SAT1', 'LAIR2', 'TNFRSF18'],
  'NK Cell': ['FCER1G', 'TYROBP', 'XCL1', 'XCL2', 'AREG', 'KRT81', 'CTSW', 'KLRD1', 'GNLY'],
}")
# JI Figure S4B & Cell Type Annotation
# IGJ: JCHAIN
py_run_string(paste0("sc.pl.dotplot(", adata_str, ", marker_dict_JI, groupby='leiden_2_200', layer='normalized', dendrogram='dendrogram_leiden_2_200', save='Signatures_marker_dict_JI.pdf')"))


######
for(ii in unlist(py$annotation_tmp)){
  py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=7, save='_detail_", ii, ".pdf')"))
}
density_visualization(adata_str, visualization_embedding, "StudyID", ".pdf", "density_color", "6", as.character(py$dotsize))
#
cell_annotation = rep(NA, eval(parse(text = paste0("py$", adata_str, "$n_obs"))))
names(cell_annotation) = eval(parse(text = paste0("py$", adata_str, "$obs_names$values")))
#
leiden_2_200 = as.character(eval(parse(text = paste0("py$", adata_str, "")))$obs[['leiden_2_200']])
#
cell_annotation[leiden_2_200 == "0"] = "LYM - Th"
cell_annotation[leiden_2_200 == "1"] = "LYM - Th"
cell_annotation[leiden_2_200 == "2"] = "LYM - Th"
cell_annotation[leiden_2_200 == "3"] = "LYM - Tc(NKG7/GZMK/CD8A/CD8B)"
cell_annotation[leiden_2_200 == "4"] = "LYM - Th(IL13)"
cell_annotation[leiden_2_200 == "5"] = "LYM - Th"
cell_annotation[leiden_2_200 == "6"] = "LYM - ILC1/3(XCL1/XCL2/SPINK2)"
cell_annotation[leiden_2_200 == "7"] = "LYM - NK(NKG7/GZMB/FCGR3A)"
cell_annotation[leiden_2_200 == "8"] = "LYM - Treg(FOXP3/TIGIT/CTLA4)"
cell_annotation[leiden_2_200 == "9"] = "LYM - Th"
cell_annotation[leiden_2_200 == "10"] = "LYM - Th"
cell_annotation[leiden_2_200 == "11"] = "LYM - Tc(CD8A/CD8B)"
cell_annotation[leiden_2_200 == "12"] = "LYM - Th"
cell_annotation[leiden_2_200 == "13"] = "LYM - MIX(KC)"
cell_annotation[leiden_2_200 == "14"] = "LYM - Th(IL17A/CCL20/CCR6/PDE4D)"
cell_annotation[leiden_2_200 == "15"] = "LYM - Tc(NKG7/GZMK/CD8A/CD8B)"
cell_annotation[leiden_2_200 == "16"] = "LYM - ILC1/NK(NKG7/GZMK/XCL1/XCL2)"
cell_annotation[leiden_2_200 == "17"] = "LYM - T_Naive(CCR7/SELL)" # Ji et al Fig S4 B
cell_annotation[leiden_2_200 == "18"] = "LYM - Th"
cell_annotation[leiden_2_200 == "19"] = "LYM - Tc(CD8A/CD8B)"
cell_annotation[leiden_2_200 == "20"] = "LYM - T_Cycling(MKI67/STMN1/PCLAF/TYMS)" # Ji et al Fig S4 A
cell_annotation[leiden_2_200 == "21"] = "LYM - ILC1/3(XCL1/XCL2/SPINK2)"
cell_annotation[leiden_2_200 == "22"] = "LYM - Plasma"
cell_annotation[leiden_2_200 == "23"] = "LYM - ADIRF/IGFBP7/CCL2"
cell_annotation[leiden_2_200 == "24"] = "LYM - MIX(MYE)"
cell_annotation[leiden_2_200 == "25"] = "LYM - Th"
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
py_run_string("marker_dict_0 = {
  'Tc': ['CD8A', 'CD8B'],
  'NK': ['NKG7', 'GZMB', 'FCGR3A'],
  'ILC': ['GZMK', 'XCL1', 'XCL2', 'SPINK2'],
  'Th_sub': ['IL13', 'IL17A', 'CCL20', 'CCR6', 'PDE4D'],
  'Naive': ['CCR7', 'SELL'],
  'Treg': ['FOXP3', 'TIGIT', 'CTLA4'],
  'Plasma': ['JCHAIN', 'CD79A'],
  'Cycling': ['MKI67', 'STMN1', 'PCLAF', 'TYMS'],
  'ADIRF/IGFBP7/CCL2': ['ADIRF', 'IGFBP7', 'CCL2'],
  'MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10'],
  'MIX(MYE)': ['HLA-DRA', 'HLA-DPA1', 'CD74']
}")
py_run_string("marker_dict_1 = {
  'Tc': ['CD8A', 'CD8B'],
  'NK': ['NKG7', 'GZMB', 'FCGR3A'],
  'ILC': ['GZMK', 'XCL1', 'XCL2', 'SPINK2'],
  'Th': ['CD40LG'],
  'Th_sub': ['IL13', 'IL17A', 'CCL20', 'CCR6', 'PDE4D'],
  'Naive': ['CCR7', 'SELL'],
  'Treg': ['FOXP3', 'TIGIT', 'CTLA4'],
  'Plasma': ['JCHAIN', 'CD79A'],
  'Cycling': ['MKI67', 'STMN1', 'PCLAF', 'TYMS'],
  'ADIRF/IGFBP7/CCL2': ['ADIRF', 'IGFBP7', 'CCL2'],
  'MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10'],
  'MIX(MYE)': ['HLA-DRA', 'HLA-DPA1', 'CD74']
}")
py_run_string("marker_dict = {
  'Tc': ['CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD8B'],
  'NK': ['NCR1', 'NCAM1', 'NKG7', 'GZMB', 'FCGR3A'],
  'ILC': ['IL7R', 'GZMK', 'XCL1', 'XCL2', 'SPINK2'],
  'Th': ['CD40LG'],
  'Th_sub': ['IL13', 'IL17A', 'CCL20', 'CCR6', 'PDE4D'],
  'Naive': ['CCR7', 'SELL'],
  'Treg': ['FOXP3', 'TIGIT', 'CTLA4'],
  'Plasma': ['JCHAIN', 'CD79A'],
  'Cycling': ['MKI67', 'STMN1', 'PCLAF', 'TYMS'],
  'ADIRF/IGFBP7/CCL2': ['ADIRF', 'IGFBP7', 'CCL2'],
  'MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10'],
  'MIX(MYE)': ['HLA-DRA', 'HLA-DPA1', 'CD74']
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
  "0" = "LYM - Tc(CD8A/CD8B)",
  "1" = "LYM - Tc(NKG7/GZMK/CD8A/CD8B)",
  "2" = "LYM - NK(NKG7/GZMB/FCGR3A)",
  "3" = "LYM - ILC1/NK(NKG7/GZMK/XCL1/XCL2)",
  "4" = "LYM - ILC1/3(XCL1/XCL2/SPINK2)",
  "5" = "LYM - Th",
  "6" = "LYM - Th(IL13)",
  "7" = "LYM - Th(IL17A/CCL20/CCR6/PDE4D)",
  "8" = "LYM - T_Naive(CCR7/SELL)",
  "9" = "LYM - Treg(FOXP3/TIGIT/CTLA4)",
  "10" = "LYM - Plasma",
  "11" = "LYM - T_Cycling(MKI67/STMN1/PCLAF/TYMS)",
  "12" = "LYM - ADIRF/IGFBP7/CCL2",
  "13" = "LYM - MIX(KC)",
  "14" = "LYM - MIX(MYE)"
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
py_run_string("color_dict={'0': colorsys.hsv_to_rgb(0.9, 0.9, 0.9),
                           '1': colorsys.hsv_to_rgb(0.9, 1, 0.7),
                           '2': colorsys.hsv_to_rgb(0, 0.6, 0.6),
                           '3': colorsys.hsv_to_rgb(0.1, 1, 0.7),
                           '4': colorsys.hsv_to_rgb(0.1, 0.9, 0.9),
                           '5': colorsys.hsv_to_rgb(0.2, 0.6, 0.6),
                           '6': colorsys.hsv_to_rgb(0.2, 0.7, 0.4),
                           '7': colorsys.hsv_to_rgb(0.2, 0.9, 0.9),
                           '8': colorsys.hsv_to_rgb(0.3, 1, 0.7),
                           '9': colorsys.hsv_to_rgb(0.5, 0.6, 0.6),
                           '10': colorsys.hsv_to_rgb(0.6, 0.6, 0.6),
                           '11': colorsys.hsv_to_rgb(0.7, 0.4, 1),
                           '12': colorsys.hsv_to_rgb(0.8, 0.7, 0.4),
                           '13': colorsys.hsv_to_rgb(0, 0.075, 0.9),
                           '14': colorsys.hsv_to_rgb(0.5, 0.075, 0.9)
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





# 20260421
py_run_string(paste0(adata_str, " = ad.read_h5ad('", py$this_output_dir, "/IntegratedAnnotation_2_20250531.h5ad')"))
DEG_key = "Annotation_Level_2_FULL"
py_run_string("marker_dict = {
  'Tc': ['CD3D', 'CD8A', 'CD8B'],
  'NK': ['NKG7', 'GZMB', 'FCGR3A'],
  'ILC': ['GZMK', 'XCL1', 'XCL2', 'SPINK2', 'IL7R'],
  'Th': ['CD40LG'],
  'Th-sub': ['IL13', 'IL17A', 'CCL20', 'CCR6', 'PDE4D'],
  'Naive': ['CCR7', 'SELL'],
  'Treg': ['FOXP3', 'TIGIT', 'CTLA4'],
  'Plasma': ['JCHAIN', 'CD79A'],
  'Cycling': ['MKI67', 'STMN1', 'PCLAF', 'TYMS'],
  'ADIRF': ['ADIRF', 'IGFBP7', 'CCL2'],
  'KC\\nMIX': ['KRT14', 'KRT5', 'KRT1', 'KRT10'],
  'MYE\\nMIX': ['HLA-DRA', 'HLA-DPA1', 'CD74']
}")
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale = 'var', swap_axes = False, var_group_rotation=0)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_var_20260421.pdf', bbox_inches = 'tight')"))


