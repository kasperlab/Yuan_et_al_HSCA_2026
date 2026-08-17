############
# Fibroblast
############
this_celltype = "Fibroblast"
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
  #py_run_string(paste0(adata_str, " = ad.read_h5ad('", py$this_output_dir, "/clustering_20250118.h5ad')"))
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
DEG_keys = c("leiden_2_030", "leiden_2_150", "leiden_2_300")
for(DEG_key in DEG_keys){
  this_path = paste0(py$this_output_dir, "/DEG_", adata_str, "_" , DEG_key, ".rds")
  if(file.exists(this_path)){
    DEG_result = readRDS(this_path)
  }else{
    DEG_result = DEG_Analysis(eval(parse(text = paste0("py$", adata_str))), DEG_key, py$this_output_dir, prefix = paste0(adata_str, "_" , DEG_key, "_"), top_DE_range = top_DE_range)
    saveRDS(DEG_result, this_path)
  }
  this_table = DEG_result[["rest"]][["table"]]
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
    py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, ".pdf', bbox_inches = 'tight')"))
  }else{
    for(ii in 1:ceiling(length(top_DEGs) / max_gene_number)){
      this_index = seq(max_gene_number) + max_gene_number * (ii - 1)
      this_index = this_index[this_index <= length(top_DEGs)]
      py$tmp_genes = top_DEGs[this_index]
      py$tmp_genes_ordered = top_DEGs_ordered[this_index]
      #py_run_string(paste0("sc.pl.heatmap(", adata_str, ", var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, show_gene_labels=True, save='_", adata_str, "_" , DEG_key, "_ordered_", ii, ".pdf')"))
      py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, n_genes=10, swap_axes = False, highlight_genes = DEG_dict, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
      py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", ii, ".pdf', bbox_inches = 'tight')"))
      py_run_string("del tmp_genes")
    }
  }
  #
  # DEGs
  #
  #########
  # Markers
  #########
  py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=[
              'PDGFRA', 'DPT', 'COL1A2', 'COL3A1', 'TWIST2', 'VIM', 'EN1', 'DLK1', 'HIC1',
              'DPP4', 'DLK1', 'ATXN1', 'PRDM1', 'EPHB2', 'LRIG1', 'TRPS1',
              'ACTA2', 'ITGA8', 'ITGA5', 'COL11A1', 'ACAN', 'CD200', 'MYH10', 'MYLK', 'MYL9',
              'SOX2', 'LEF1', 'CRABP1', 'RSPO3', 'CORIN', 'VCAN', 'ALPL',
              'SOX2', 'ACTA2', 'ITGA8', 'MGP', 'ACAN',
              'ACTA2', 'TAGLN', 'FN1', 
              'CRABP1', 'DES', 'PRSS35'
              ], groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_Max_review.pdf')"))
  py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=[
              'CCN5', 'SLPI', 'CTHRC1', 'MFAP5', 'TSPAN8',
              'CCL19', 'APOE', 'CXCL2', 'CXCL3', 'EFEMP1',
              'APCDD1', 'ID1', 'WIF1', 'COL18A1', 'PTGDS',
              'ASPN', 'POSTN', 'GPC3', 'TNN', 'SFRP1'
              ], groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_Signatures_Boldo_t1.pdf')"))
  py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=[
              'CPE', 'SFRP1', 'FGFBP2', 'IGFBP2', 'OLFML2A',
              'ECRG4', 'TM4SF1', 'NR2F2', 'APOD', 'CLDN1'
              ], groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_Signatures_11_12.pdf')"))
  
  py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=[
              'CPE', 'SFRP1', 'FGFBP2', 'IGFBP2', 'OLFML2A',
              'ECRG4', 'TM4SF1', 'NR2F2', 'APOD', 'CLDN1'
              ], groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_Signatures_11_12.pdf')"))
}
py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=[
              'COCH', # leiden_2_300 - #5
              'CXCL12', # leiden_2_300 - #27
              'SLC5A3', 'TFAP2A', # leiden_2_300 - #30
              'POSTN', 'GPC3', 'COL11A1', # leiden_2_300 - #17, #19, #32
              'ASPN', 'TNN', 'SFRP1', # Mesenchymal
              'DSP', 'DSC3', 'KRT14', 'KRT5', 'DMKN', # #13
              'TAGLN', 'TM4SF1', 'ECRG4', # #11
              'SFRP1', 'IGFBP2', # #12
              'ACTA2', 'RGS5', # leiden_2_300 - #39
              'HLA-DRA', 'HLA-DPA1', 'CD74', # #14
              'PTPRC', 'CD3D', # #16
              'CCL19', # #1, #17
              'MEDAG', 'SERPINE2', 'PNPLA8', # #19
              'CCN5', 'SLPI', 'TSPAN8', # Secretory-reticular, #6
              'ISG15', 'MX1', 'IFIT3', 'IFI44L', 'OAS1', # #18
              'COL23A1', # #3
              'APCDD1', 'WIF1', # Secretory-papillary
              ], groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_Signatures_previous.pdf')"))
######
for(ii in unlist(py$annotation_tmp)){
  py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=7, save='_detail_", ii, ".pdf')"))
}
density_visualization(adata_str, visualization_embedding, "StudyID", ".pdf", "density_color", "6", as.character(py$dotsize))
#
cell_annotation = rep(NA, eval(parse(text = paste0("py$", adata_str, "$n_obs"))))
names(cell_annotation) = eval(parse(text = paste0("py$", adata_str, "$obs_names$values")))
#
leiden_2_030 = as.character(eval(parse(text = paste0("py$", adata_str, "")))$obs[['leiden_2_030']])
leiden_2_150 = as.character(eval(parse(text = paste0("py$", adata_str, "")))$obs[['leiden_2_150']])
leiden_2_300 = as.character(eval(parse(text = paste0("py$", adata_str, "")))$obs[['leiden_2_300']])
#
cell_annotation[leiden_2_150 == "0"] = "FIB - Bridge Cell"
cell_annotation[leiden_2_150 == "1"] = "FIB - Secretory-reticular"
cell_annotation[leiden_2_150 == "2"] = "FIB - Secretory-papillary"
cell_annotation[leiden_2_150 == "3"] = "FIB - Pro-inflammatory(MEDAG/HMOX1/PTGS2/IL6/GPC3)"
cell_annotation[leiden_2_150 == "4"] = "FIB - Secretory-reticular"
cell_annotation[leiden_2_150 == "5"] = "FIB - Pro-inflammatory(MEDAG/HMOX1/PTGS2/IL6/CCL19)"
cell_annotation[leiden_2_150 == "6"] = "FIB - Pro-inflammatory(CCL19/IGFBP3)"
cell_annotation[leiden_2_150 == "7"] = "FIB - Secretory-papillary/Secretory-reticular"
cell_annotation[leiden_2_150 == "8"] = "FIB - MYOC" # Gur_Amit_Cell_2022
cell_annotation[leiden_2_150 == "9"] = "FIB - Mesenchymal(POSTN/GPC3/COL11A1/DPEP1)"
cell_annotation[leiden_2_150 == "10"] = "FIB - Mesenchymal(COCH)"
cell_annotation[leiden_2_150 == "11"] = "FIB - Secretory-papillary"
cell_annotation[leiden_2_150 == "12"] = "FIB - Pro-inflammatory/Secretory-papillary"
#cell_annotation[leiden_2_150 == "13"] = "FIB - TAGLN/TM4SF1/ECRG4/CLDN1"
cell_annotation[leiden_2_150 == "14"] = "FIB - Pro-inflammatory"
cell_annotation[leiden_2_150 == "15"] = "FIB - MIX(KC)"
cell_annotation[leiden_2_150 == "16"] = "FIB - IGFBP2/FGFBP2/OLFML2A" # Gur_Amit_Cell_2022
cell_annotation[leiden_2_150 == "17"] = "FIB - Mesenchymal(POSTN/GPC3/PI16/CXCL12)"
cell_annotation[leiden_2_150 == "18"] = "FIB - Mesenchymal(POSTN/COL18A1/BGN/LOXL1/TNC)"
cell_annotation[leiden_2_150 == "19"] = "FIB - Mesenchymal(INHBA/TFAP2A/WNT5A)"
cell_annotation[leiden_2_150 == "20"] = "FIB - ISG15/MX1/IFIT3/IFI44L/OAS1"
cell_annotation[leiden_2_150 == "21"] = "FIB - MIX(LYM)"
cell_annotation[leiden_2_150 == "22"] = "FIB - Secretory-papillary(Stressed)"
cell_annotation[leiden_2_150 == "23"] = "FIB - MIX(MYE)"
cell_annotation[leiden_2_150 == "24"] = "FIB - MIX(MUR)"
cell_annotation[leiden_2_150 == "25"] = "FIB - MYOC(Stressed)"
cell_annotation[leiden_2_150 == "26"] = "FIB - Secretory-reticular(Stressed)"
#cell_annotation[leiden_2_150 == "27"] = "FIB - Bridge Cell(Stressed)"
#
cell_annotation[leiden_2_150 == "13" & leiden_2_300 == "31"] = "FIB - ECRG4/EBF2/CDH19/ANGPTL7"
cell_annotation[leiden_2_150 == "13" & leiden_2_300 == "32"] = "FIB - ECRG4/EBF2/ITGA6"
#
print(leiden_2_300[leiden_2_150 == "13" & !leiden_2_300 %in% c("31", "32")])
print(table(cell_annotation[leiden_2_300 == "21"]))
cell_annotation[leiden_2_150 == "13" & !leiden_2_300 %in% c("31", "32")] = "FIB - Pro-inflammatory"
#
cell_annotation[leiden_2_030 == "3"] = "FIB - Bridge Cell"
cell_annotation[leiden_2_030 == "3" & leiden_2_150 == "27"] = "FIB - Bridge Cell(Stressed)"
print(table(leiden_2_300[is.na(cell_annotation)]))
uncovered_mask = leiden_2_150 %in% c("0", "27") & !leiden_2_300 == "3"
cell_annotation[uncovered_mask & leiden_2_300 == "15"] = names(table(cell_annotation[leiden_2_300 == "15"]))[1]
cell_annotation[uncovered_mask & leiden_2_300 == "18"] = names(table(cell_annotation[leiden_2_300 == "18"]))[1]
cell_annotation[uncovered_mask & leiden_2_300 == "26"] = names(table(cell_annotation[leiden_2_300 == "26"]))[1]
cell_annotation[uncovered_mask & leiden_2_300 == "28"] = names(table(cell_annotation[leiden_2_300 == "28"]))[1]
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
  'Reticular/Pre-adipocytes': ['MYOC', 'FMO1', 'PPARG', 'FABP4'],
  'FRC-like': ['CCL19', 'CD74', 'TNFSF13B', 'VCAM1'],
  'FRC-like subtypes': ['MEDAG', 'PTGS2', 'IL6', 'IGFBP3'],
  'Papillary': ['LEPR', 'WIF1', 'APCDD1', 'NKD2', 'COL23A1', 'COL13A1', 'NPTX2', 'COL6A5', 'HSPB3', 'CLEC2A'],
  'Reticular': ['THY1', 'CCN5', 'APOD'],
  'Reticular (stressed)': ['HSPB2'],
  'Interferon': ['ISG15', 'MX1', 'IFIT3', 'IFI44L', 'OAS1'],
  'HF-associated': ['ASPN', 'TNMD', 'PPP1R14A', 'TNN'],
  'HF-associated subtypes': ['INHBA', 'TFAP2A', 'WNT5A', 'COCH', 'POSTN', 'COL18A1', 'BGN', 'LOXL1', 'TNC', 'GPC3', 'COL11A1', 'DPEP1', 'PI16', 'CXCL12'],
  'Schwann-like (RAMP1+)': ['RAMP1', 'IGFBP2', 'FGFBP2', 'OLFML2A'],
  'Schwann-like (NGFR+)': ['NGFR', 'ECRG4', 'EBF2'],
  'Schwann-like (NGFR+) subtypes': ['CDH19', 'ANGPTL7', 'ITGA6'],
  'Bridge (Reticular-Adipo) (stressed)': ['PTGS2', 'S100A2', 'HMOX1'],
  'MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10'],
  'MIX(LYM)': ['PTPRC', 'CD3D'],
  'MIX(MYE)': ['HLA-DRA', 'HLA-DPA1', 'CD74'],
  'MIX(MUR)': ['MCAM', 'RGS5']
}")
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=True, n_genes=10, standard_scale='var', swap_axes=False, highlight_genes=None, highlight_params={'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line=False, cluster_line_highlight_params={'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_var.pdf', bbox_inches='tight')"))
#######
# Write
#######
#h5ad_path = paste0(py$this_output_dir, "/IntegratedAnnotation_2.h5ad")
#py_run_string(paste0("", adata_str, ".write('", h5ad_path, "')"))
h5ad_path = paste0(py$this_output_dir, "/IntegratedAnnotation_2_20250428.h5ad")
py_run_string(paste0("", adata_str, ".write('", h5ad_path, "')"))
#
# Further analysis
#
Annotation_Level_2 = as.character(eval(parse(text = paste0("py$", adata_str, "$obs[['IntegratedAnnotation_2']]"))))
print(table(Annotation_Level_2))
#
py$index_name_dict = list(
  "0" = "FIB - MYOC",
  "1" = "FIB - MYOC(Stressed)",
  "2" = "FIB - Pro-inflammatory",
  "3" = "FIB - Pro-inflammatory(MEDAG/HMOX1/PTGS2/IL6/GPC3)",
  "4" = "FIB - Pro-inflammatory(MEDAG/HMOX1/PTGS2/IL6/CCL19)",
  "5" = "FIB - Pro-inflammatory(CCL19/IGFBP3)",
  "6" = "FIB - Pro-inflammatory/Secretory-papillary",
  "7" = "FIB - Secretory-papillary",
  "8" = "FIB - Secretory-papillary(Stressed)",
  "9" = "FIB - Secretory-papillary/Secretory-reticular",
  "10" = "FIB - Secretory-reticular",
  "11" = "FIB - Secretory-reticular(Stressed)",
  "12" = "FIB - ISG15/MX1/IFIT3/IFI44L/OAS1",
  "13" = "FIB - Mesenchymal(INHBA/TFAP2A/WNT5A)",
  "14" = "FIB - Mesenchymal(COCH)",
  "15" = "FIB - Mesenchymal(POSTN/COL18A1/BGN/LOXL1/TNC)",
  "16" = "FIB - Mesenchymal(POSTN/GPC3/COL11A1/DPEP1)",
  "17" = "FIB - Mesenchymal(POSTN/GPC3/PI16/CXCL12)",
  "18" = "FIB - IGFBP2/FGFBP2/OLFML2A",
  "19" = "FIB - ECRG4/EBF2/CDH19/ANGPTL7",
  "20" = "FIB - ECRG4/EBF2/ITGA6",
  "21" = "FIB - Bridge Cell",
  "22" = "FIB - Bridge Cell(Stressed)",
  "23" = "FIB - MIX(KC)",
  "24" = "FIB - MIX(LYM)",
  "25" = "FIB - MIX(MYE)",
  "26" = "FIB - MIX(MUR)"
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
py_run_string("color_dict={'0': colorsys.hsv_to_rgb(0.15, 0.4, 1),
                           '1': colorsys.hsv_to_rgb(0.15, 0.9, 0.9),
                           '2': colorsys.hsv_to_rgb(0.05, 1, 1),
                           '3': colorsys.hsv_to_rgb(0.05, 0.4, 1),
                           '4': colorsys.hsv_to_rgb(0.1, 0.9, 0.9),
                           '5': colorsys.hsv_to_rgb(0.1, 0.4, 1),
                           '6': colorsys.hsv_to_rgb(0.0, 1, 0.7),
                           '7': colorsys.hsv_to_rgb(0.9, 0.9, 0.9),
                           '8': colorsys.hsv_to_rgb(0.9, 0.7, 0.4),
                           '9': colorsys.hsv_to_rgb(0.9, 1, 0.7),
                           '10': colorsys.hsv_to_rgb(0.35, 1, 0.7),
                           '11': colorsys.hsv_to_rgb(0.35, 0.9, 0.9),
                           '12': colorsys.hsv_to_rgb(0.9, 0.4, 1),
                           '13': colorsys.hsv_to_rgb(0.6, 1, 0.7),
                           '14': colorsys.hsv_to_rgb(0.55, 1, 0.7),
                           '15': colorsys.hsv_to_rgb(0.6, 0.4, 1),
                           '16': colorsys.hsv_to_rgb(0.6, 0.9, 0.9),
                           '17': colorsys.hsv_to_rgb(0.55, 0.4, 1),
                           '18': colorsys.hsv_to_rgb(0.55, 0.9, 0.9),
                           '19': colorsys.hsv_to_rgb(0.75, 0.4, 1),
                           '20': colorsys.hsv_to_rgb(0.75, 0.9, 0.9),
                           '21': colorsys.hsv_to_rgb(0.7, 1, 0.7),
                           '22': colorsys.hsv_to_rgb(0.15, 1, 0.7),
                           '23': colorsys.hsv_to_rgb(0, 0.05, 0.9),
                           '24': colorsys.hsv_to_rgb(0.25, 0.05, 0.9),
                           '25': colorsys.hsv_to_rgb(0.5, 0.05, 0.9),
                           '26': colorsys.hsv_to_rgb(0.75, 0.05, 0.9)
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
                                      save='_Annotation_Level_2_OC_NL.png'
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
                                      save='_Annotation_Level_2_OC_RM.pdf'
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



# 20260402
py_run_string(paste0(adata_str, " = ad.read_h5ad('", py$this_output_dir, "/IntegratedAnnotation_2_20250531.h5ad')"))
#
py_run_file(py_visualization, convert = F)
py_run_string("from matplotlib import font_manager as fm, rcParams")
py_run_string("fm.fontManager.addfont('/home/haoy/projects/Analysis/fonts/ARIAL.TTF')")
py_run_string("fm.fontManager.addfont('/home/haoy/projects/Analysis/fonts/ARIALBD.TTF')")
py_run_string("plt.rcParams['font.family'] = 'Arial'")
DEG_key = "Annotation_Level_2_FULL"
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_var.pdf', bbox_inches = 'tight')"))
#
for(ii in unlist(py$marker_dict)){
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color='", ii, "', layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=6, frameon=False, save='_", ii, ".pdf')"))

}


# 20260521
py_run_string(paste0(adata_str, " = ad.read_h5ad('", py$this_output_dir, "/IntegratedAnnotation_2_20250531.h5ad')"))
#
py_run_file(py_visualization, convert = F)
py_run_string("from matplotlib import font_manager as fm, rcParams")
py_run_string("fm.fontManager.addfont('/home/haoy/projects/Analysis/fonts/ARIAL.TTF')")
py_run_string("fm.fontManager.addfont('/home/haoy/projects/Analysis/fonts/ARIALBD.TTF')")
py_run_string("plt.rcParams['font.family'] = 'Arial'")
DEG_key = "Annotation_Level_2_FULL"
py_run_string("marker_dict = {
  'Reticular/Pre-adipocytes': ['MYOC', 'PPARG'],
  'FRC-like': ['APOE', 'CCL19', 'CD74', 'TNFSF13B', 'VCAM1', 'C7', 'APOC1', 'IL33'],
  'FRC-like subtypes': ['MEDAG', 'PTGS2', 'IL6', 'IGFBP3'],
  'Papillary': ['LEPR', 'WIF1', 'APCDD1', 'NKD2', 'COL23A1', 'COL13A1', 'NPTX2', 'COL6A5', 'HSPB3', 'CLEC2A'],
  'Reticular': ['SLPI', 'C1QTNF3', 'HSPB2'],
  'Interferon': ['ISG15', 'MX1', 'IFIT3', 'IFI44L', 'OAS1'],
  'HF-associated': ['ASPN', 'TNN', 'CRABP1', 'COL24A1'],
  'HF-associated subtypes': ['INHBA', 'PTCH1', 'TFAP2A', 'WNT5A', 'HHIP', 'COCH', 'POSTN', 'COL18A1', 'BGN', 'LOXL1', 'TNC', 'GPC3', 'COL11A1', 'DPEP1', 'PI16', 'CXCL12'],
  'Schwann-like (RAMP1+)': ['RAMP1', 'IGFBP2', 'FGFBP2', 'OLFML2A'],
  'Schwann-like (NGFR+)': ['NGFR', 'ECRG4', 'EBF2', 'CDH19', 'ANGPTL7', 'ITGA6'],
  'Bridge (Reticular-Adipo) (stressed)': ['PTGS2', 'S100A2', 'HMOX1'],
  'MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10'],
  'MIX(LYM)': ['PTPRC', 'CD3D'],
  'MIX(MYE)': ['HLA-DRA', 'HLA-DPA1', 'CD74'],
  'MIX(MUR)': ['MCAM', 'RGS5']
}")
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_var.pdf', bbox_inches = 'tight')"))
#

