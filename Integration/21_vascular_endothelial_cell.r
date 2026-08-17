###########################
# Vascular Endothelial Cell
###########################
this_celltype = "Vascular Endothelial Cell"
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
DEG_key = "leiden_2_008"
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
######
# 10.1161/CIRCULATIONAHA.120.052318 Figure 3
py_run_string("marker_dict_CIRCULATIONAHA_120_052318_F3 = {
  'VEC - Arterial': ['DKK2', 'IGFBP3', 'SERPINE2', 'CLDN10', 'GJA5', 'CXCL12', 'BMX', 'LTBP4', 'HEY1', 'SOX5', 'SEMA3G'],
  'VEC - Aerocyte': ['SOSTDC1', 'EDNRB', 'HPGD', 'CYP3A5', 'PRKG1', 'TBX2', 'RCSD1', 'EDA', 'B3GALNT1', 'EXPH5', 'NCALD', 'S100A4'],
  'VEC - Aerocyte/General Capillary': ['CA4', 'AFF3', 'ADGRL2', 'BTNL9', 'RGCC', 'ADGRF5', 'KIAA1217'],
  'VEC - General Capillary': ['FCN3', 'IL7R', 'CD36', 'NRXN3', 'SLC6A4', 'GPIHBP1', 'ARHGAP18', 'IL18R1'],
  'VEC - Pulmonary-venous': ['CPE', 'CLU', 'C7', 'PTGS1', 'EFEMP1', 'MMRN1', 'PKHD1L1', 'PDZRN4', 'DKK3', 'PLAT', 'CDH11', 'HDAC9'],
  'VEC - Pulmonary-venous/Systemic-venous': ['ACKR1', 'IGFBP7', 'MCTP1', 'VWF'],
  'VEC - Systemic-venous': ['COL15A1', 'ZNF385D', 'EBF1', 'TSHZ2', 'FLRT2', 'OLFM1', 'CPXM2', 'PLVAP', 'TPD52L1', 'PDE7B', 'VWA1', 'SPRY1']
}")
py_run_string(paste0("sc.pl.dotplot(", adata_str, ", marker_dict_CIRCULATIONAHA_120_052318_F3, groupby='leiden_2_008', layer='normalized', dendrogram='dendrogram_leiden_2_008', save='Signatures_marker_dict_CIRCULATIONAHA_120_052318_F3.pdf')"))
# 10.1093/nsr/nwae231 Figure 1
py_run_string("marker_dict_nsr_nwae231_F1 = {
  'VEC - Arteries': ['GJA5', 'FBLN5', 'GJA4'],
  'VEC - Capillaries': ['CA4', 'CD36', 'RGCC'],
  'VEC - Hypoxia': ['MT1X', 'MT1E', 'MT2A'],
  'VEC - Lymphatics': ['PROX1', 'LYVE1', 'CCL21'],
  'VEC - Tip Cell': ['COL4A1', 'KDR', 'ESM1'],
  'VEC - Veins': ['ACKR1', 'SELP', 'CLU']
}")
py_run_string(paste0("sc.pl.dotplot(", adata_str, ", marker_dict_nsr_nwae231_F1, groupby='leiden_2_008', layer='normalized', dendrogram='dendrogram_leiden_2_008', save='Signatures_marker_dict_nsr_nwae231_F1.pdf')"))
# 10.1093/nsr/nwae231 Figure S1
# ATP5E: ATP5F1E
# ATP5G2: ATP5MC2
# ATP5L: ATP5MG
# H2AFZ: H2AZ1
# PTRF: CAVIN1
py_run_string("marker_dict_nsr_nwae231_FS1 = {
  'VEC - Tip Cell': ['KDR', 'APOE', 'APOA2', 'CXCR4', 'ESM1', 'ACKR3', 'DNAJB1', 'JUND', 'FOSB', 'ATP5F1E', 'ATP5MC2', 'ATP5MG', 'STMN1', 'MKI67', 'H2AZ1'],
  'VEC - Veins': ['CLU', 'COL4A1', 'COL4A2', 'SELE', 'CCL2', 'ICAM1', 'FABP4', 'CD36', 'IFITM1', 'IGFBP5'],
  'VEC - Capillaries': ['TMEM100', 'HPGD', 'FABP4', 'MARCKS', 'RGCC', 'JUN', 'CD320'],
  'VEC - Arteries': ['DNAJB1', 'JUND', 'FOSB', 'CLU', 'MALAT1', 'CAVIN1', 'COL4A1', 'COL4A2', 'HSPG2']
}")
py_run_string(paste0("sc.pl.dotplot(", adata_str, ", marker_dict_nsr_nwae231_FS1, groupby='leiden_2_008', layer='normalized', dendrogram='dendrogram_leiden_2_008', save='Signatures_marker_dict_nsr_nwae231_FS1.pdf')"))
# thno.54917 Figure 1
# SLC9A3R2: NHERF2
py_run_string("marker_dict_thno_54917_F1 = {
  'VEC - A': ['IGFBP3', 'SRGN', 'SAT1', 'SEMA3G', 'NHERF2', 'RHOB', 'HEY1', 'BST2', 'JAG1', 'SOX17'],
  'VEC - C1': ['RBP7', 'RGCC', 'TSC22D1', 'SPRY1', 'SPARC', 'PLVAP', 'ITIH5', 'A2M', 'AQP1', 'ADGRF5'],
  'VEC - C2': ['ACKR1', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DRA', 'CD74', 'HLA-DQA1', 'HLA-DPB1', 'SPARCL1', 'CTSC', 'HLA-DRB5'],
  'VEC - P': ['SELE', 'C2CD4B', 'SOCS3', 'ZFP36', 'HSPA1A', 'ATF3', 'HSPA1B', 'CDKN1A', 'NFKBIA', 'DNAJB1'],
  'VEC - V': ['CYP1B1', 'CLU', 'VWF', 'MGP', 'IER3', 'PERP', 'TSC22D3', 'EDN1', 'FBLN2', 'LYST']
}")
py_run_string(paste0("sc.pl.dotplot(", adata_str, ", marker_dict_thno_54917_F1, groupby='leiden_2_008', layer='normalized', dendrogram='dendrogram_leiden_2_008', save='Signatures_marker_dict_thno_54917_F1.pdf')"))
# PMC8997372 Figure 1
# SLC9A3R2: NHERF2
py_run_string("marker_dict_PMC8997372_F1 = {
  'VEC - A': ['SEMA3G', 'GJA4', 'GJA5', 'CXCL12', 'HEY1', 'NOTCH4'],
  'VEC - PAC': ['ASS1', 'S100A4', 'RGCC'],
  'VEC - PVC': ['EFNB2', 'SOX17', 'PLAUR', 'KLHL21', 'GJA1', 'LITAF'],
  'VEC - PCV': ['ICAM1', 'IRF1', 'SELE', 'SELP', 'LRG1', 'EGR2'],
  'VEC - V': ['ACKR1', 'CCL14', 'CLU', 'VWF'],
}")
py_run_string(paste0("sc.pl.dotplot(", adata_str, ", marker_dict_PMC8997372_F1, groupby='leiden_2_008', layer='normalized', dendrogram='dendrogram_leiden_2_008', save='Signatures_marker_dict_PMC8997372_F1.pdf')"))
######
for(ii in unlist(py$annotation_tmp)){
  py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=7, save='_detail_", ii, ".pdf')"))
}
density_visualization(adata_str, visualization_embedding, "StudyID", ".pdf", "density_color", "6", as.character(py$dotsize))
#
cell_annotation = rep(NA, eval(parse(text = paste0("py$", adata_str, "$n_obs"))))
names(cell_annotation) = eval(parse(text = paste0("py$", adata_str, "$obs_names$values")))
#
leiden_2_008 = as.character(eval(parse(text = paste0("py$", adata_str, "")))$obs[['leiden_2_008']])
#
cell_annotation[leiden_2_008 == "0"] = "VEC - Collecting Venules"
cell_annotation[leiden_2_008 == "1"] = "VEC - Pre-venular Capillaries(S100A4/EFNB2/SOX17/ICAM1/IRF1/SELE)"
cell_annotation[leiden_2_008 == "2"] = "VEC - Arterioles(S100A4/SEMA3G/GJA4/CXCL12/IGFBP3)"
cell_annotation[leiden_2_008 == "3"] = "VEC - MIX(MUR/FIB)"
cell_annotation[leiden_2_008 == "4"] = "VEC - MIX(KC)"
cell_annotation[leiden_2_008 == "5"] = "VEC - MIX(LYM)"
cell_annotation[leiden_2_008 == "6"] = "VEC - MIX(MEL)"
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
  'Collecting Venules': ['S100A4'],
  'Pre-venular Capillaries': ['EFNB2', 'SOX17', 'ICAM1', 'IRF1', 'SELE'],
  'Arterioles': ['SEMA3G', 'GJA4', 'CXCL12', 'IGFBP3'],
  'MIX(MUR/FIB)': ['RGS5', 'MYL9', 'COL1A2', 'COL6A2'],
  'MIX(KC)': ['KRT14', 'KRT5', 'KRT1', 'KRT10'],
  'MIX(LYM)': ['CD3D', 'CD3E'],
  'MIX(MEL)': ['DCT', 'TYRP1', 'MLANA', 'PMEL']
}")
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_dict, groupby='", DEG_key, "', dendrogram=False, layer='normalized', n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = marker_dict, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_var.pdf', bbox_inches = 'tight')"))
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
  "0" = "VEC - Collecting Venules",
  "1" = "VEC - Pre-venular Capillaries(S100A4/EFNB2/SOX17/ICAM1/IRF1/SELE)",
  "2" = "VEC - Arterioles(S100A4/SEMA3G/GJA4/CXCL12/IGFBP3)",
  "3" = "VEC - MIX(MUR/FIB)",
  "4" = "VEC - MIX(KC)",
  "5" = "VEC - MIX(LYM)",
  "6" = "VEC - MIX(MEL)"
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
                           '1': colorsys.hsv_to_rgb(0.1, 0.9, 0.9),
                           '2': colorsys.hsv_to_rgb(0.3, 0.9, 0.9),
                           '3': colorsys.hsv_to_rgb(0, 0.05, 0.9),
                           '4': colorsys.hsv_to_rgb(0.25, 0.05, 0.9),
                           '5': colorsys.hsv_to_rgb(0.5, 0.05, 0.9),
                           '6': colorsys.hsv_to_rgb(0.75, 0.05, 0.9)
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






