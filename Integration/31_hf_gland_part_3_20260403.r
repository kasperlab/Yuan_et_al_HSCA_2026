##########
# HF&Gland
##########
this_celltype = "HF_Gland"
this_celltype_name = gsub("_", "", gsub("[ -]", "", this_celltype), fixed = T)
py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs_Karl_20250603_20241209")
py_run_string("this_figure_dir = this_output_dir + '/figures'")
py_run_string("sc.settings.figdir = this_figure_dir")
dir.create(py$this_figure_dir, showWarnings = F, recursive = T)
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
clustering_keys = readRDS("/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss/HFGland_2000_hvgs_50_pcs_Karl/clustering_keys.rds")
tryCatch({
  stop()
  py_run_string(paste0("print(", adata_str, ")"))
}, error = function(x){
  py_run_string(paste0(adata_str, " = ad.read_h5ad('/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss/HFGland_2000_hvgs_50_pcs_Karl/IntegratedAnnotation_3_20241209.h5ad')"))
  py_run_string(paste0("del ", adata_str, ".uns"))
  print(paste0("Loaded ", this_celltype_name, " from cache."))
})
######
py_run_string("annotation_tmp = copy.deepcopy(annotation_1)")
cell_number = eval(parse(text = paste0("py$", adata_str, "$n_obs")))
py$dotsize = max(c(8, get_dot_size(cell_number)))
#
# Annotation
#
#
# DEG
#
DEG_keys = c("PCA_use_leiden_080", "PCA_use_leiden_150")
for(DEG_key in DEG_keys){
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
      #for(this_study in eval(parse(text = paste0("unique(as.character(py$", adata_str, "$obs[['StudyID']]))")))){
      #  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, "[", adata_str, ".obs['StudyID'] == '", this_study, "'], var_names=DEG_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
      #  py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", this_study, "_", used_table, ".pdf', bbox_inches = 'tight')"))
      #}
    }else{
      for(ii in 1:ceiling(length(top_DEGs) / max_gene_number)){
        this_index = seq(max_gene_number) + max_gene_number * (ii - 1)
        this_index = this_index[this_index <= length(top_DEGs)]
        py$tmp_genes = top_DEGs[this_index]
        py$tmp_genes_ordered = top_DEGs_ordered[this_index]
        #py_run_string(paste0("sc.pl.heatmap(", adata_str, ", var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, show_gene_labels=True, save='_", adata_str, "_" , DEG_key, "_ordered_", ii, ".pdf')"))
        py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
        py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", ii, "_", used_table, ".pdf', bbox_inches = 'tight')"))
        #for(this_study in eval(parse(text = paste0("unique(as.character(py$", adata_str, "$obs[['StudyID']]))")))){
        #  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, "[", adata_str, ".obs['StudyID'] == '", this_study, "'], var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
        #  py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", ii, "_", this_study, "_", used_table, ".pdf', bbox_inches = 'tight')"))
        #}
        py_run_string("del tmp_genes")
      }
    }
  }
}
######
######
cell_annotation = rep(NA, eval(parse(text = paste0("py$", adata_str, "$n_obs"))))
names(cell_annotation) = eval(parse(text = paste0("py$", adata_str, "$obs_names$values")))
#
leiden_050 = as.character(eval(parse(text = paste0("py$", adata_str, "")))$obs[['PCA_use_leiden_050']])
leiden_080 = as.character(eval(parse(text = paste0("py$", adata_str, "")))$obs[['PCA_use_leiden_080']])
leiden_150 = as.character(eval(parse(text = paste0("py$", adata_str, "")))$obs[['PCA_use_leiden_150']])
#
cell_annotation[leiden_080 == "0"] = "KC - Bulge - DIO2"
cell_annotation[leiden_080 == "1"] = "KC - Gland Transition"
cell_annotation[leiden_080 == "2"] = "KC - Upper Hair Follicle - KRT1"
cell_annotation[leiden_080 == "3"] = "KC - Upper Hair Follicle - PTN"
cell_annotation[leiden_080 == "4"] = "KC - Channel - KRT23/CLCA2"
cell_annotation[leiden_080 == "5"] = "KC - ORS/CP - KRT5/KRT17/GJA1"
cell_annotation[leiden_080 == "6"] = "KC - Sweat Gland - KRT19/AQP5/CA6/SNORC"
cell_annotation[leiden_080 == "7"] = "KC - Channel - KRT23/SEMA3C/FGF7/LAMB3/PLAUR/MYC/IER3/KRT17"
cell_annotation[leiden_080 == "8"] = "KC - Upper Hair Follicle - COL17A1/POSTN"
cell_annotation[leiden_080 == "9"] = "KC - Sebaceous Gland - MGST1/IL1R2"
cell_annotation[leiden_080 == "10"] = "KC - Channel - KRT23/SEMA3C/FGF7/LAMB3"
cell_annotation[leiden_080 == "11"] = "KC - Sweat Gland - KRT19/AQP5/PIP/DCD/MUCL1"
cell_annotation[leiden_080 == "12"] = "KC - Upper Hair Follicle - PTN/SPRR1B"
cell_annotation[leiden_080 == "13"] = "KC - Upper Hair Follicle - NCOA7"
cell_annotation[leiden_080 == "14"] = "KC - Sweat Duct - WFDC3/MMP7"
cell_annotation[leiden_080 == "15"] = "KC - Lower Bulge/Hair Germ/Lower Proximal Cup - LGR5"
cell_annotation[leiden_080 == "16"] = "KC - Sebaceous Gland - MGST1"
cell_annotation[leiden_080 == "17"] = "KC - Sebaceous Gland - MGST1/FASN/AWAT2"
cell_annotation[leiden_080 == "18"] = "KC - Isthmus - DIO2/PTHLH"
cell_annotation[leiden_080 == "19"] = "KC - Germinative Layer - TK1/MCM3/HELLS"
cell_annotation[leiden_080 == "20"] = "KC - Matrix/Cortex/Medulla - KRT35/KRT85/MT4" # Takahashi
cell_annotation[leiden_080 == "21"] = "KC - Gland Transition - APOE/KRT1/KRT10/CLEC2B"
cell_annotation[leiden_080 == "22"] = "KC - Sweat Gland - KRT19/AQP5/WFDC2/SLPI/AZGP1"
cell_annotation[leiden_080 == "23"] = "KC - IRS H/H - KRT25/KRT27/KRT28/KRT71/TCHH" # Takahashi
#
cell_annotation[leiden_080 == "1" & leiden_150 == "24"] = "KC - Gland Transition - APOE"
cell_annotation[leiden_080 == "4" & leiden_150 == "7"] = "KC - Channel - KRT23/CLCA2/KRT17"
cell_annotation[leiden_080 == "4" & leiden_150 == "19"] = "KC - Channel - KRT23/CLCA2"
cell_annotation[leiden_080 == "5" & leiden_150 == "5"] = "KC - ORS/CP - KRT5/KRT17/GJA1"
cell_annotation[leiden_080 == "5" & leiden_150 == "26"] = "KC - ORS/CP - KRT5/KRT17/GJA1/KRT6B/CLDN4"
cell_annotation[leiden_080 == "6" & leiden_150 == "9"] = "KC - Sweat Gland - KRT19/AQP5/CA6/SNORC"
cell_annotation[leiden_080 == "6" & leiden_150 == "17"] = "KC - Sweat Gland - KRT19/AQP5/CA6/SNORC/CLDN10"
cell_annotation[leiden_080 == "15" & leiden_150 == "20"] = "KC - Lower Bulge - LGR5/TIMP3/COMP/MGP" # Takahashi
cell_annotation[leiden_080 == "15" & leiden_150 == "28"] = "KC - Hair Germ/Lower Proximal Cup - LGR5/EPCAM/FBLN1"
# Fix
cell_annotation[leiden_050 == "0" & leiden_080 == "23"] = "KC - Gland Transition"
cell_annotation[leiden_080 == "15" & (leiden_150 != "20" & leiden_150 != "28")] = "KC - Lower Bulge - LGR5/TIMP3/COMP/MGP"
#print(table(cell_annotation[leiden_2_300 == "21"]))
#
print(sum(is.na(cell_annotation)))
######
eval(parse(text = paste0("py$", adata_str, "$obs[['IntegratedAnnotation_3']] = as.factor(cell_annotation)")))
summary_list[[qc_state]][["IntegratedAnnotation_3"]] = list("AllCells" = table(cell_annotation))
py_run_string("annotation_tmp = annotation_tmp if np.isin('IntegratedAnnotation_2', annotation_tmp) else annotation_tmp + ['IntegratedAnnotation_2']")
py_run_string("annotation_tmp = annotation_tmp if np.isin('IntegratedAnnotation_3', annotation_tmp) else annotation_tmp + ['IntegratedAnnotation_3']")
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=annotation_tmp, legend_loc='on data', legend_fontsize=6, legend_fontoutline=0.5, size=dotsize, frameon=False, ncols=4, save='_", adata_str, "_annotation_OD.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=annotation_tmp, legend_loc='right margin', legend_fontsize=6, legend_fontoutline=0.5, size=dotsize, frameon=False, ncols=1, save='_", adata_str, "_annotation.pdf')"))
######
ii = "IntegratedAnnotation_3"
py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=7, save='_detail_", ii, ".pdf')"))
#########
# Markers
#########
# DEG
DEG_key = "IntegratedAnnotation_3"
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
    #for(this_study in eval(parse(text = paste0("unique(as.character(py$", adata_str, "$obs[['StudyID']]))")))){
    #  py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, "[", adata_str, ".obs['StudyID'] == '", this_study, "'], var_names=DEG_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
    #  py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", this_study, "_", used_table, ".pdf', bbox_inches = 'tight')"))
    #}
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
py_run_string("marker_dict = {
  'Bulge (Related)': ['DIO2', 'LGR5', 'TIMP3', 'COMP', 'MGP', 'EPCAM', 'FBLN1', 'PTHLH'],
  'IRS/Cortex': ['KRT25', 'KRT27', 'KRT28', 'KRT71', 'TCHH', 'KRT35', 'KRT85', 'MT4'],
  'GL/ORS/CP': ['TK1', 'MCM3', 'HELLS', 'KRT5', 'GJA1', 'KRT6B', 'CLDN4'],
  'uHF': ['PTN', 'SPRR1B', 'NCOA7', 'KRT1', 'COL17A1', 'POSTN'],
  'Gland Transition': ['APOE', 'KRT10', 'CLEC2B'],
  'Sebaceous Gland': ['MGST1', 'FASN', 'AWAT2', 'IL1R2'],
  'Sweat Gland/Duct': ['KRT19', 'AQP5', 'CA6', 'SNORC', 'CLDN10', 'WFDC2', 'SLPI', 'AZGP1', 'PIP', 'DCD', 'MUCL1', 'WFDC3', 'MMP7'],
  'Channel': ['KRT23', 'SEMA3C', 'FGF7', 'LAMB3', 'PLAUR', 'MYC', 'IER3', 'KRT17', 'CLCA2']
}")
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=True, n_genes=10, standard_scale='var', swap_axes=False, highlight_genes=None, highlight_params={'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line=False, cluster_line_highlight_params={'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_var.pdf', bbox_inches='tight')"))
#######
# Write
#######
h5ad_path = paste0(py$this_output_dir, "/IntegratedAnnotation_3.h5ad")
py_run_string(paste0("", adata_str, ".write('", h5ad_path, "')"))
#
# Further analysis
#
Annotation_Level_3 = as.character(eval(parse(text = paste0("py$", adata_str, "$obs[['IntegratedAnnotation_3']]"))))
print(table(Annotation_Level_3))
#
py$index_name_dict = list(
  "0" = "KC - Bulge - DIO2",
  "1" = "KC - Lower Bulge - LGR5/TIMP3/COMP/MGP",
  "2" = "KC - Hair Germ/Lower Proximal Cup - LGR5/EPCAM/FBLN1",
  "3" = "KC - Isthmus - DIO2/PTHLH",
  "4" = "KC - IRS H/H - KRT25/KRT27/KRT28/KRT71/TCHH",
  "5" = "KC - Matrix/Cortex/Medulla - KRT35/KRT85/MT4",
  "6" = "KC - Germinative Layer - TK1/MCM3/HELLS",
  "7" = "KC - ORS/CP - KRT5/KRT17/GJA1",
  "8" = "KC - ORS/CP - KRT5/KRT17/GJA1/KRT6B/CLDN4",
  "9" = "KC - Upper Hair Follicle - PTN",
  "10" = "KC - Upper Hair Follicle - PTN/SPRR1B",
  "11" = "KC - Upper Hair Follicle - NCOA7",
  "12" = "KC - Upper Hair Follicle - KRT1",
  "13" = "KC - Upper Hair Follicle - COL17A1/POSTN",
  "14" = "KC - Gland Transition",
  "15" = "KC - Gland Transition - APOE",
  "16" = "KC - Gland Transition - APOE/KRT1/KRT10/CLEC2B",
  "17" = "KC - Sebaceous Gland - MGST1",
  "18" = "KC - Sebaceous Gland - MGST1/FASN/AWAT2",
  "19" = "KC - Sebaceous Gland - MGST1/IL1R2",
  "20" = "KC - Sweat Gland - KRT19/AQP5/CA6/SNORC",
  "21" = "KC - Sweat Gland - KRT19/AQP5/CA6/SNORC/CLDN10",
  "22" = "KC - Sweat Gland - KRT19/AQP5/WFDC2/SLPI/AZGP1",
  "23" = "KC - Sweat Gland - KRT19/AQP5/PIP/DCD/MUCL1",
  "24" = "KC - Sweat Duct - WFDC3/MMP7",
  "25" = "KC - Channel - KRT23/SEMA3C/FGF7/LAMB3",
  "26" = "KC - Channel - KRT23/SEMA3C/FGF7/LAMB3/PLAUR/MYC/IER3/KRT17",
  "27" = "KC - Channel - KRT23/CLCA2/KRT17",
  "28" = "KC - Channel - KRT23/CLCA2"
)
py_run_string(paste0(adata_str, ".obs['Annotation_Level_3_renumber'] = ", adata_str, ".obs['IntegratedAnnotation_3'].map({name: f'{index}' for index, name in index_name_dict.items()})"))
py_run_string(paste0(adata_str, ".obs['Annotation_Level_3_FULL'] = ", adata_str, ".obs['IntegratedAnnotation_3'].map({name: f'{index}: {name}' for index, name in index_name_dict.items()})"))
py_run_string("new_order = [f'{index}: {name}' for index, name in index_name_dict.items()]")
#
py_run_string("color_dict={'0': colorsys.hsv_to_rgb(0.9, 1, 0.7),
                           '1': colorsys.hsv_to_rgb(0.9, 0.9, 0.9),
                           '2': colorsys.hsv_to_rgb(0.9, 0.4, 1),
                           '3': colorsys.hsv_to_rgb(0.9, 0.7, 0.4),
                           '4': colorsys.hsv_to_rgb(0.0, 1, 0.7),
                           '5': colorsys.hsv_to_rgb(0.0, 0.9, 0.9),
                           '6': colorsys.hsv_to_rgb(0.1, 0.4, 1),
                           '7': colorsys.hsv_to_rgb(0.15, 1, 0.7),
                           '8': colorsys.hsv_to_rgb(0.15, 0.9, 0.9),
                           '9': colorsys.hsv_to_rgb(0.25, 1, 0.7),
                           '10': colorsys.hsv_to_rgb(0.25, 0.9, 0.9),
                           '11': colorsys.hsv_to_rgb(0.25, 0.4, 1),
                           '12': colorsys.hsv_to_rgb(0.25, 0.7, 0.4),
                           '13': colorsys.hsv_to_rgb(0.3, 0.6, 0.6),
                           '14': colorsys.hsv_to_rgb(0.4, 1, 0.7),
                           '15': colorsys.hsv_to_rgb(0.4, 0.9, 0.9),
                           '16': colorsys.hsv_to_rgb(0.4, 0.4, 1),
                           '17': colorsys.hsv_to_rgb(0.5, 1, 0.7),
                           '18': colorsys.hsv_to_rgb(0.5, 0.9, 0.9),
                           '19': colorsys.hsv_to_rgb(0.5, 0.4, 1),
                           '20': colorsys.hsv_to_rgb(0.6, 1, 0.7),
                           '21': colorsys.hsv_to_rgb(0.6, 0.9, 0.9),
                           '22': colorsys.hsv_to_rgb(0.6, 0.4, 1),
                           '23': colorsys.hsv_to_rgb(0.6, 0.7, 0.4),
                           '24': colorsys.hsv_to_rgb(0.65, 0.6, 0.6),
                           '25': colorsys.hsv_to_rgb(0.75, 1, 0.7),
                           '26': colorsys.hsv_to_rgb(0.5, 0.9, 0.9),
                           '27': colorsys.hsv_to_rgb(0.75, 0.4, 1),
                           '28': colorsys.hsv_to_rgb(0.75, 0.7, 0.4)
                          }")
py_run_string("color_dict_FULL = {f'{index}: {name}': color_dict[index] for index, name in index_name_dict.items()}")
#
py_run_string(paste0(adata_str, ".obs['Annotation_Level_3_FULL'] = ", adata_str, ".obs['Annotation_Level_3_FULL'].cat.reorder_categories(new_order)"))
py_run_string(paste0("print(", adata_str, ".obs['Annotation_Level_3_FULL'].cat.categories)"))
#
py_run_string(paste0("tmp_adata = ad.AnnData(obs=", adata_str, ".obs.copy(), obsm={'", visualization_embedding, "': ", adata_str, ".obsm['", visualization_embedding, "'].copy()})"))
py_run_string("np.random.seed(0)")
py_run_string("tmp_adata_random = tmp_adata[np.random.permutation(list(range(tmp_adata.n_obs))), :]")
py_run_string("del tmp_adata")
py_run_string("gc.collect()")
py_run_string(paste0("sc.pl.embedding(tmp_adata_random,
                                      basis='", visualization_embedding, "',
                                      color=['Annotation_Level_3_FULL'],
                                      legend_loc='none',
                                      legend_fontoutline=0.5,
                                      legend_fontsize=6,
                                      show=False,
                                      ncols=5,
                                      size=8,
                                      frameon=False,
                                      palette=color_dict_FULL,
                                      add_outline=False,
                                      save='_Annotation_Level_3_NL.png'
)"))
py_run_string(paste0("fig = sc.pl.embedding(tmp_adata_random,
                                            basis='", visualization_embedding, "',
                                            color=['Annotation_Level_3_renumber'],
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
                                      color=['Annotation_Level_3_FULL'],
                                      legend_loc='right margin',
                                      legend_fontsize=6,
                                      show=False,
                                      size=8,
                                      add_outline=False,
                                      palette=color_dict_FULL,
                                      ax=fig.axes[0],
                                      save='_Annotation_Level_3_RM.pdf'
)"))
py_run_string("del tmp_adata_random")
py_run_string("gc.collect()")
#
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[x for y in marker_dict.values() for x in y], layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=6, frameon=False, save='_marker_genes.pdf')"))
#
DEG_key = "Annotation_Level_3_FULL"
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_var.pdf', bbox_inches = 'tight')"))
#
ii = "Annotation_Level_3_FULL"
py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', mask=", adata_str, ".obs['StudyID'] == 'SS3_HumanSkin_20K', size=dotsize, ncols=4, save='_detail_", ii, "_SS3_HumanSkin_20K.pdf')"))
py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, ncols=4, save='_detail_", ii, ".pdf')"))
#
py_run_string(paste0("print(", adata_str, "[", adata_str, ".obs['StudyID'] == 'SS3_HumanSkin_20K'].obs['Annotation_Level_3_FULL'].value_counts())"))
#
py$sankey_categories = c("StudyID", "Annotation_Level_3_FULL", "AnatomicalRegionLevel2")
py_run_string(paste0("sankey_fig = plot_sankey_diagram(", adata_str, ", categories=sankey_categories, return_fig=True)"))
py_run_string(paste0("sankey_fig.write_image(f'{this_figure_dir}/sankey_", adata_str, "_", paste(py$sankey_categories, collapse = "_"), ".pdf')"))
py$sankey_categories = c("StudyID", "Annotation_Level_3_FULL", "AnatomicalRegionLevel2", "LibraryPlatform")
py_run_string(paste0("sankey_fig = plot_sankey_diagram(", adata_str, ", categories=sankey_categories, return_fig=True)"))
py_run_string(paste0("sankey_fig.write_image(f'{this_figure_dir}/sankey_", adata_str, "_", paste(py$sankey_categories, collapse = "_"), ".pdf')"))
py$sankey_categories = c("StudyID", "Annotation_Level_3_FULL", "LibraryPlatform")
py_run_string(paste0("sankey_fig = plot_sankey_diagram(", adata_str, ", categories=sankey_categories, return_fig=True)"))
py_run_string(paste0("sankey_fig.write_image(f'{this_figure_dir}/sankey_", adata_str, "_", paste(py$sankey_categories, collapse = "_"), ".pdf')"))
#######
# Write
#######
h5ad_path = paste0(py$this_output_dir, "/IntegratedAnnotation_3_20250603_20241209.h5ad")
py_run_string(paste0("", adata_str, ".write('", h5ad_path, "')"))



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
py_run_string(paste0("clustering_plot(", adata_str, ", '", bodysite_slot, "', basis='", visualization_embedding, "', size=dotsize, colorbar_loc=None, ncols=5, save='_detail_", bodysite_slot, ".pdf')"))
#
eval(parse(text = paste0("StudyID = py$", adata_str, "$obs[['StudyID']]")))
eval(parse(text = paste0("py$", adata_str, "$obs[['StudyID_ordered']] = factor(StudyID, levels = StudyID[!duplicated(StudyID)])")))
density_visualization(adata_str, visualization_embedding, "StudyID_ordered", ".pdf", "density_color", "4", as.character(py$dotsize))


py_run_string(paste0("adata_KC = adata[np.logical_and(adata.obs['integrated_annotation_L1'].astype(str) == 'Keratinocyte', adata.obs['integrated_annotation_L2'].astype(str) != 'nan')]"))
py_run_string(paste0("adata_KC.obs['isHFG'] = np.where(adata_KC.obs['integrated_annotation_L3'].astype(str).str.contains('KC - HFGland', regex=False), 'KC - HFGland', 'Others')"))
py_run_string(paste0("clustering_plot(adata_KC, 'isHFG', basis='integrated_UMAP_L2', size=8, ncols=5, save='_detail_isHFG_size_8.pdf')"))







# 20251002
py_run_string(paste0(adata_str, " = ad.read_h5ad('", py$this_output_dir, "/IntegratedAnnotation_3_20250603_20241209.h5ad')"))

DEG_key = "Annotation_Level_3_renumber"
clustering_embedding = "PCA_use"
py_run_string(paste0("print(", adata_str, ".obs['", DEG_key, "'].cat.categories.size)"))
py_run_string(paste0("sc.tl.dendrogram(", adata_str, ", groupby='", DEG_key, "', use_rep='", clustering_embedding, "')"))


py_run_string(paste0("fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(", adata_str, ".obs['", DEG_key, "'].cat.categories.size * 0.3, 20))"))
py_run_string(paste0("sc.pl.dendrogram(", adata_str, ", '", DEG_key, "', ax=axs)"))
py_run_string(paste0("fig.savefig(f'{this_figure_dir}/DendrogramPlot_", adata_str, "_", DEG_key, ".pdf')"))
py_run_string("plt.close('all')")
py_run_string("gc.collect()")



# 20260403
py_run_string(paste0(adata_str, " = ad.read_h5ad('", py$this_output_dir, "/IntegratedAnnotation_3_20250603_20241209.h5ad')"))
visualization_embedding = "X_umap_3_imbalance"
for(ii in c("GJB2", "KRT23", "KRT5")){
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color='", ii, "', layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=6, frameon=False, save='_", ii, ".pdf')"))
}

DEG_key = "Annotation_Level_3_FULL"
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=['GJB2', 'KRT23'], groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes=None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_GJB2&KRT23_var.pdf', bbox_inches = 'tight')"))




# 20260422
py_run_string(paste0(adata_str, " = ad.read_h5ad('", py$this_output_dir, "/IntegratedAnnotation_3_20250603_20241209.h5ad')"))
py_run_string(paste0(adata_str, "_subset = ", adata_str, "[np.logical_not(", adata_str, ".obs['Annotation_Level_3_renumber'].isin(np.array([14, 15, 16, 20, 21, 22, 23, 24]).astype(str)))].copy()"))
py_run_string("color_dict={'0': colorsys.hsv_to_rgb(0.9, 1, 0.7),
                           '1': colorsys.hsv_to_rgb(0.9, 0.9, 0.9),
                           '2': colorsys.hsv_to_rgb(0.9, 0.4, 1),
                           '3': colorsys.hsv_to_rgb(0.9, 0.7, 0.4),
                           '4': colorsys.hsv_to_rgb(0.0, 1, 0.7),
                           '5': colorsys.hsv_to_rgb(0.0, 0.9, 0.9),
                           '6': colorsys.hsv_to_rgb(0.1, 0.4, 1),
                           '7': colorsys.hsv_to_rgb(0.15, 1, 0.7),
                           '8': colorsys.hsv_to_rgb(0.15, 0.9, 0.9),
                           '9': colorsys.hsv_to_rgb(0.25, 1, 0.7),
                           '10': colorsys.hsv_to_rgb(0.25, 0.9, 0.9),
                           '11': colorsys.hsv_to_rgb(0.25, 0.4, 1),
                           '12': colorsys.hsv_to_rgb(0.25, 0.7, 0.4),
                           '13': colorsys.hsv_to_rgb(0.3, 0.6, 0.6),
                           '14': colorsys.hsv_to_rgb(0.4, 1, 0.7),
                           '15': colorsys.hsv_to_rgb(0.4, 0.9, 0.9),
                           '16': colorsys.hsv_to_rgb(0.4, 0.4, 1),
                           '17': colorsys.hsv_to_rgb(0.5, 1, 0.7),
                           '18': colorsys.hsv_to_rgb(0.5, 0.9, 0.9),
                           '19': colorsys.hsv_to_rgb(0.5, 0.4, 1),
                           '20': colorsys.hsv_to_rgb(0.6, 1, 0.7),
                           '21': colorsys.hsv_to_rgb(0.6, 0.9, 0.9),
                           '22': colorsys.hsv_to_rgb(0.6, 0.4, 1),
                           '23': colorsys.hsv_to_rgb(0.6, 0.7, 0.4),
                           '24': colorsys.hsv_to_rgb(0.65, 0.6, 0.6),
                           '25': colorsys.hsv_to_rgb(0.75, 1, 0.7),
                           '26': colorsys.hsv_to_rgb(0.5, 0.9, 0.9),
                           '27': colorsys.hsv_to_rgb(0.75, 0.4, 1),
                           '28': colorsys.hsv_to_rgb(0.75, 0.7, 0.4)
                          }")
visualization_embedding = "X_umap_3_imbalance"
py_run_string(paste0("sc.pl.embedding(", adata_str, "_subset, '", visualization_embedding, "', color='Annotation_Level_3_renumber', layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=6, frameon=False, save='_subset_Annotation_Level_3_renumber.pdf')"))









