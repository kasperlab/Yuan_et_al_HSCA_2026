##############
# Myeloid Cell
##############
this_celltype = "Myeloid Cell"
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

py_run_string(paste0("sc.pl.dotplot(", adata_str, ", [
              'CD68', 'C1QB', 'C1QC', 'CD163', 'MARCO', 'FCGR2A', 'CTSB', 'F13A1', 
              'NR4A1', 'NR4A2', 'KLF4', 'IL23A', 'CLEC9A', 'CLEC10A', 'CD207', 'CCR7', 
              'LAMP3', 'CD40', 'CD274', 'IDO1', 'CD200', 'PDCD1LG2', 'SOCS1', 'CD80', 
              'CD83', 'CD86'
                     ], groupby='leiden_2_100', layer='normalized', dendrogram='dendrogram_leiden_2_100', save='Signatures_Muzz_F3B.pdf')"))

py_run_string(paste0("sc.pl.dotplot(", adata_str, ", [
              'CD14', 'CD163', 'F13A1', 'CD1A', 'CD207', 'CD1C', 'MRC1', 'IL1B', 
              'IRF8', 'FCER1A', 'LAMP3', 'CD83', 'CCR7', 'LILRA4', 'CLEC4C'
                     ], groupby='leiden_2_100', layer='normalized', dendrogram='dendrogram_leiden_2_100', save='Signatures_Rojahn_F4.pdf')"))

py_run_string(paste0("sc.pl.dotplot(", adata_str, ", [
              'CD68', 'IGFBP5', 'CD1C', 'CLEC9A', 'CD207', 'CD68', 'S100A8', 'CLEC4C'
                     ], groupby='leiden_2_100', layer='normalized', dendrogram='dendrogram_leiden_2_100', save='Signatures_Ji_FS1C.pdf')"))

py_run_string(paste0("sc.pl.dotplot(", adata_str, ", [
              'CD14', 'CD68', 'S100A8', 'MKI67', 'CCR7', 'AXL', 'IGFBP5'
                     ], groupby='leiden_2_100', layer='normalized', dendrogram='dendrogram_leiden_2_100', save='Signatures_Ji_F1C.pdf')"))

py_run_string(paste0("sc.pl.dotplot(", adata_str, ", [
              'MNDA', 'IGSF6'
                     ], groupby='leiden_2_100', layer='normalized', dendrogram='dendrogram_leiden_2_100', save='Signatures_C17.pdf')"))

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
py_run_string(paste0("sc.pl.dotplot(", adata_str, ", marker_dict_JI, groupby='leiden_2_100', layer='normalized', dendrogram='dendrogram_leiden_2_100', save='Signatures_marker_dict_JI.pdf')"))


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
cell_annotation[leiden_2_100 == "0"] = "MYE - LC(CD207)" # Ji_Khavari_Cell_2020 Fig S1
cell_annotation[leiden_2_100 == "1"] = "MYE - DC(IL1B/IL23A/CCR7)" # Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020
cell_annotation[leiden_2_100 == "2"] = "MYE - Macrophage(CD163/MRC1/LYVE1)"
cell_annotation[leiden_2_100 == "3"] = "MYE - DC/Macrophage"
cell_annotation[leiden_2_100 == "4"] = "MYE - DC(CD1C/CLEC10A)" # Ji_Khavari_Cell_2020 Fig S1
cell_annotation[leiden_2_100 == "5"] = "MYE - Macrophage(CD163/CXCL2/CXCL3)"
cell_annotation[leiden_2_100 == "6"] = "MYE - DC(CLEC9A)" # Ji_Khavari_Cell_2020 Fig S1, PDC @ Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020
cell_annotation[leiden_2_100 == "7"] = "MYE - MIX(KC)"
cell_annotation[leiden_2_100 == "8"] = "MYE - MDSC(S100A8/S100A9/TREM1)" # myeloid-derived suppressor cells @ Ji_Khavari_Cell_2020 Fig S1
cell_annotation[leiden_2_100 == "9"] = "MYE - Cycling(PCLAF/CKS1B/MKI67)"
cell_annotation[leiden_2_100 == "10"] = "MYE - MIX(FIB)"
cell_annotation[leiden_2_100 == "11"] = "MYE - DC/Macrophage"
cell_annotation[leiden_2_100 == "12"] = "MYE - DC(CCL17)"
cell_annotation[leiden_2_100 == "13"] = "MYE - DC_migratory(LAMP3/CCR7)" # Reynolds_Haniffa_Science_2021 Fig 3 B & Ji_Khavari_Cell_2020 Fig S4B
cell_annotation[leiden_2_100 == "14"] = "MYE - Monocyte-derived Macrophage(MSR1/FCGR3A/OLR1)"
cell_annotation[leiden_2_100 == "15"] = "MYE - LC(CD207)" # Ji_Khavari_Cell_2020 Fig S1
cell_annotation[leiden_2_100 == "16"] = "MYE - MIX(MUR)"
cell_annotation[leiden_2_100 == "17"] = "MYE - DC(CD1C/CLEC10A)"
cell_annotation[leiden_2_100 == "18"] = "MYE - Macrophage(STAT1/GBP1)"
cell_annotation[leiden_2_100 == "19"] = "MYE - Plasma(JCHAIN/IGKC)"
cell_annotation[leiden_2_100 == "20"] = "MYE - Monocyte-derived Macrophage(MSR1/GPNMB)"
cell_annotation[leiden_2_100 == "21"] = "MYE - MIX(MEL)"
cell_annotation[leiden_2_100 == "22"] = "MYE - MIX(LYM)"
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
  'Monocyte': ['CD14'],
  'MDSC': ['S100A8', 'S100A9', 'TREM1'],
  'Macrophage': ['CD68', 'CD163', 'MARCO'],
  'Macrophage_sub': ['MSR1', 'GPNMB', 'FCGR3A', 'OLR1', 'STAT1', 'GBP1', 'CXCL2', 'CXCL3', 'MRC1', 'LYVE1'],
  'Dendritic Cell': ['CD1C', 'CLEC10A', 'IL1B'],
  'DC_sub': ['IL23A', 'CCL17', 'CLEC9A', 'LAMP3', 'CCR7'],
  'Langerhans Cell': ['CD207', 'CD1A'],
  'Cycling': ['MKI67', 'PCLAF', 'CKS1B'],
  'Plasma': ['JCHAIN', 'IGKC'],
  'MIX(LYM)': ['CD3D', 'TRBC2'],
  'MIX(KC)': ['KRT10', 'KRT14'],
  'MIX(MUR)': ['TAGLN', 'MYL9'],
  'MIX(FIB)': ['COL1A2', 'COL6A2'],
  'MIX(MEL)': ['DCT', 'TYRP1']
}")
py_run_string("marker_dict = {
  'Monocyte/Macrophage': ['CD14', 'CD163', 'MARCO'],
  'Monocyte_sub': ['S100A8', 'S100A9', 'TREM1', 'MSR1', 'GPNMB', 'FCGR3A', 'OLR1', 'STAT1', 'GBP1', 'CXCL2', 'CXCL3', 'MRC1', 'LYVE1'],
  'Dendritic Cell': ['CD1C', 'CLEC10A', 'IL1B'],
  'DC_sub': ['IL23A', 'CCL17', 'CLEC9A', 'LAMP3', 'CCR7'],
  'Langerhans Cell': ['CD207', 'CD1A'],
  'Cycling': ['MKI67', 'PCLAF', 'CKS1B'],
  'Plasma': ['JCHAIN', 'IGKC'],
  'MIX(LYM)': ['CD3D', 'TRBC2'],
  'MIX(KC)': ['KRT10', 'KRT14'],
  'MIX(MUR)': ['TAGLN', 'MYL9'],
  'MIX(FIB)': ['COL1A2', 'COL6A2'],
  'MIX(MEL)': ['DCT', 'TYRP1']
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
  "0" = "MYE - MDSC(S100A8/S100A9/TREM1)",
  "1" = "MYE - Monocyte-derived Macrophage(MSR1/GPNMB)",
  "2" = "MYE - Monocyte-derived Macrophage(MSR1/FCGR3A/OLR1)",
  "3" = "MYE - Macrophage(STAT1/GBP1)",
  "4" = "MYE - Macrophage(CD163/CXCL2/CXCL3)",
  "5" = "MYE - Macrophage(CD163/MRC1/LYVE1)",
  "6" = "MYE - DC/Macrophage",
  "7" = "MYE - DC(CD1C/CLEC10A)",
  "8" = "MYE - DC(IL1B/IL23A/CCR7)",
  "9" = "MYE - DC(CCL17)",
  "10" = "MYE - DC(CLEC9A)",
  "11" = "MYE - DC_migratory(LAMP3/CCR7)",
  "12" = "MYE - LC(CD207)",
  "13" = "MYE - Cycling(PCLAF/CKS1B/MKI67)",
  "14" = "MYE - Plasma(JCHAIN/IGKC)",
  "15" = "MYE - MIX(LYM)",
  "16" = "MYE - MIX(KC)",
  "17" = "MYE - MIX(MUR)",
  "18" = "MYE - MIX(FIB)",
  "19" = "MYE - MIX(MEL)"
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
                           '1': colorsys.hsv_to_rgb(0.95, 0.9, 0.9),
                           '2': colorsys.hsv_to_rgb(0.95, 0.4, 0.95),
                           '3': colorsys.hsv_to_rgb(0, 0.7, 0.4),
                           '4': colorsys.hsv_to_rgb(0.05, 1, 0.7),
                           '5': colorsys.hsv_to_rgb(0.05, 0.9, 0.9),
                           '6': colorsys.hsv_to_rgb(0.15, 0.6, 0.6),
                           '7': colorsys.hsv_to_rgb(0.25, 1, 0.7),
                           '8': colorsys.hsv_to_rgb(0.25, 0.7, 0.4),
                           '9': colorsys.hsv_to_rgb(0.3, 0.4, 0.95),
                           '10': colorsys.hsv_to_rgb(0.35, 0.9, 0.9),
                           '11': colorsys.hsv_to_rgb(0.35, 1, 0.7),
                           '12': colorsys.hsv_to_rgb(0.45, 0.6, 0.6),
                           '13': colorsys.hsv_to_rgb(0.55, 1, 0.7),
                           '14': colorsys.hsv_to_rgb(0.65, 0.6, 0.6),
                           '15': colorsys.hsv_to_rgb(0, 0.075, 0.9),
                           '16': colorsys.hsv_to_rgb(0.2, 0.075, 0.9),
                           '17': colorsys.hsv_to_rgb(0.4, 0.075, 0.9),
                           '18': colorsys.hsv_to_rgb(0.6, 0.075, 0.9),
                           '19': colorsys.hsv_to_rgb(0.8, 0.075, 0.9)
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


# 20260427
py_run_string(paste0(adata_str, " = ad.read_h5ad('", py$this_output_dir, "/IntegratedAnnotation_2_20250531.h5ad')"))


py_run_string("marker_dict_0 = {
  'MD': ['CD14', 'CD68'],
  'Monocyte-derived\\nsubtypes': ['VCAN', 'FCN1', 'FCGR3A', 'FABP4', 'LPL', 'TREM2', 'APOC1', 'C3', 'TREM2', 'FCGR3A', 'S100A8', 'S100A9', 'TREM1', 'MSR1', 'GPNMB', 'FCGR3A', 'OLR1'],
  'M': ['C1QC', 'CD163', 'MARCO'],
  'Macrophage\\nsubtypes': ['CXCL9', 'CXCL10', 'IL10', 'MMP19', 'FOLR2', 'MARCO', 'STAT1', 'GBP1', 'CXCL2', 'CXCL3', 'MRC1', 'LYVE1'],
  'DC': ['CD1C', 'CLEC10A', 'IL1B'],
  'Dendritic Cell\\nsubtypes': ['FCER1A', 'IL7R', 'G0S2', 'CD1B', 'MMP12', 'LIPA', 'XCR1', 'CADM1', 'IL7R', 'ADAM12', 'IL23A', 'CCL17', 'CLEC9A', 'LAMP3', 'CCR7'],
  'LC': ['CD207', 'CD1A', 'FCGBP'],
  'Cycling': ['AURKB', 'MKI67', 'PCLAF', 'CKS1B'],
  'Plasma': ['GZMB', 'IRF7', 'LILRA4'],
  'LYM\\nMIX': ['CD3D', 'TRBC2'],
  'KC\\nMIX': ['KRT10', 'KRT14'],
  'MUR\\nMIX': ['TAGLN', 'MYL9'],
  'FIB\\nMIX': ['COL1A2', 'COL6A2'],
  'MEL\\nMIX': ['DCT', 'TYRP1']
}")

py_run_string("marker_dict = {
  'MD': ['CD14', 'CD68'],
  'Monocyte-derived\\nsubtypes': ['VCAN', 'FCN1', 'S100A8', 'S100A9', 'FABP4', 'LPL', 'TREM2', 'APOC1', 'C3', 'TREM2', 'FCGR3A'],
  'M': ['C1QC', 'CD163'],
  'Macrophage\\nsubtypes': ['STAT1', 'GBP1', 'CXCL9', 'CXCL10', 'CXCL2', 'CXCL3', 'IL10', 'MMP19', 'FOLR2', 'MARCO', 'LYVE1'],
  'DC': ['CD1C', 'CLEC10A', 'IL1B'],
  'Dendritic Cell\\nsubtypes': ['FCER1A', 'IL23A', 'G0S2', 'CCL17', 'CD1B', 'MMP12', 'LIPA', 'CLEC9A', 'XCR1', 'CADM1', 'LAMP3', 'CCR7', 'IL7R', 'ADAM12'],
  'LC': ['CD207', 'CD1A', 'FCGBP'],
  'Cycling': ['AURKB', 'MKI67', 'PCLAF', 'CKS1B'],
  'Plasma': ['GZMB', 'IRF7', 'LILRA4'],
  'LYM\\nMIX': ['CD3D', 'TRBC2'],
  'KC\\nMIX': ['KRT10', 'KRT14'],
  'MUR\\nMIX': ['TAGLN', 'MYL9'],
  'FIB\\nMIX': ['COL1A2', 'COL6A2'],
  'MEL\\nMIX': ['DCT', 'TYRP1']
}")



DEG_key = "Annotation_Level_2_FULL"
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale = 'var', swap_axes = False, var_group_rotation=0)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_var_20260427.pdf', bbox_inches = 'tight')"))



visualization_embedding = "X_umap_2_imbalance"
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=['CLEC9A', 'FCER1A', 'S100A9', 'LYVE1', 'IL7R'], layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=6, frameon=False, save='_marker_genes_new.pdf')"))
















