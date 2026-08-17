##########
# IFE_like
##########
this_celltype = "IFE_like"
this_celltype_name = gsub("_", "", gsub("[ -]", "", this_celltype), fixed = T)
py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs")
dir.create(py$this_output_dir, showWarnings = F, recursive = T)
py_run_string("this_figure_dir = this_output_dir + '/figures'")
py_run_string("sc.settings.figdir = this_figure_dir")
qc_state = "before_QC"
######
if(!exists("annotation_1", envir = .GlobalEnv)){
  py$annotation_1 = readRDS(paste0(py$output_dir_cache, "/annotation_1.rds"))
}
tryCatch({
  py_run_string(paste0("print(adata_", this_celltype_name, "_imbalance)"))
}, error = function(x){
  py_run_string(paste0("adata_", this_celltype_name, "_imbalance = ad.read_h5ad('", py$output_dir_cache, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs/clustering.h5ad')"))
  py_run_string(paste0("adata_", this_celltype_name, "_imbalance.uns['log1p'] = None"))
  py_run_string(paste0("del adata_", this_celltype_name, "_imbalance.uns['log1p']"))
  print(paste0("Loaded ", this_celltype_name, " from cache."))
})
tryCatch({
  py_run_string(paste0("print(adata_", this_celltype_name, "_imbalance_hvg)"))
}, error = function(x){
  py_run_string(paste0("adata_", this_celltype_name, "_imbalance_hvg = ad.read_h5ad('", py$output_dir_cache, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs/clustering_hvg.h5ad')"))
  py_run_string(paste0("adata_", this_celltype_name, "_imbalance_hvg.uns['log1p'] = None"))
  py_run_string(paste0("del adata_", this_celltype_name, "_imbalance_hvg.uns['log1p']"))
  print(paste0("Loaded ", this_celltype_name, "_hvg from cache."))
})
######
py_run_string("annotation_tmp = copy.deepcopy(annotation_1)")
cell_annotation = rep(NA, eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance$n_obs"))))
names(cell_annotation) = eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance$obs_names$values")))
cell_number = eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance$n_obs")))
######
leiden_3_100 = as.character(eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance")))$obs[['leiden_3_100']])
############
# Signatures
############
genename = eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance")))$var_names$values
krt_genes = genename[grep("^KRT.+", genename)]
# py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_3_imbalance', color=[], layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=3, frameon=False, save='_T_helper.pdf')"))
######
DEG_keys = c("leiden_3_100")
DEG_leiden_path_list[[this_celltype_name]] = list()
for(DEG_key in DEG_keys){
  this_path = paste0(py$this_output_dir, "/DEG_", DEG_key, ".rds")
  DEG_leiden_path_list[[this_celltype_name]][[DEG_key]] = this_path
  if(file.exists(this_path)){
    DEG_leiden_list[[this_celltype_name]][[DEG_key]] = readRDS(this_path)
  }else{
    DEG_leiden_list[[this_celltype_name]][[DEG_key]] = DEG_Analysis(eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance"))), DEG_key, py$this_output_dir, prefix = paste0(DEG_key, "_"), top_DE_range = top_DE_range)
    saveRDS(DEG_leiden_list[[this_celltype_name]][[DEG_key]], this_path)
  }
}

for(DEG_key in DEG_keys){
  this_table = DEG_leiden_list[[this_celltype_name]][[DEG_key]][["rest"]][["table"]]
  for(this_cluster in unique(this_table$Cluster)){
    py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_3_imbalance', color=[
              '", paste(head(this_table[this_table$Cluster == this_cluster, ], top_DE_range)$Gene, collapse = "', '"), "'
              ], layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=10, frameon=False, save='_Signatures_", DEG_key, "_", this_cluster, ".pdf')"))
  }
}
max_gene_number = 200
for(DEG_key in DEG_keys){
  this_table = DEG_leiden_list[[this_celltype_name]][[DEG_key]][["rest"]][["table"]]
  py_run_string(paste0("sc.tl.dendrogram(adata_", this_celltype_name, "_imbalance_hvg, groupby='", DEG_key, "', use_rep='PCA_use')"))
  py_run_string(paste0("adata_", this_celltype_name, "_imbalance.uns['dendrogram_", DEG_key, "'] = adata_", this_celltype_name, "_imbalance_hvg.uns['dendrogram_", DEG_key, "']"))
  this_clusters = unique(this_table$Cluster)
  top_DEGs = c()
  for(this_cluster in this_clusters){
    top_DEGs = c(top_DEGs, head(this_table[this_table$Cluster == this_cluster, ], 10)$Gene)
  }
  top_DEGs_ordered = c()
  for(this_cluster in eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance$uns[['dendrogram_", DEG_key, "']][['categories_ordered']]")))){
    top_DEGs_ordered = c(top_DEGs_ordered, head(this_table[this_table$Cluster == this_cluster, ], 10)$Gene)
  }
  if(length(top_DEGs) <= max_gene_number){
    py_run_string(paste0("sc.pl.heatmap(adata_", this_celltype_name, "_imbalance, var_names=[
              '", paste(top_DEGs, collapse = "', '"), "'
              ], groupby='", DEG_key, "', layer='normalized', show_gene_labels=True, save='_", this_celltype_name, "_", DEG_key, ".pdf')"))
    py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=[
              '", paste(top_DEGs, collapse = "', '"), "'
              ], groupby='", DEG_key, "', layer='normalized', save='", this_celltype_name, "_", DEG_key, ".pdf')"))
    py_run_string(paste0("sc.pl.heatmap(adata_", this_celltype_name, "_imbalance, var_names=[
              '", paste(top_DEGs_ordered, collapse = "', '"), "'
              ], groupby='", DEG_key, "', layer='normalized', dendrogram=True, show_gene_labels=True, save='_", this_celltype_name, "_", DEG_key, "_ordered.pdf')"))
    py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=[
              '", paste(top_DEGs_ordered, collapse = "', '"), "'
              ], groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_ordered.pdf')"))
  }else{
    for(ii in 1:ceiling(length(top_DEGs) / max_gene_number)){
      this_index = seq(max_gene_number) + max_gene_number * (ii - 1)
      this_index = this_index[this_index <= length(top_DEGs)]
      py$tmp_genes = top_DEGs[this_index]
      py$tmp_genes_ordered = top_DEGs_ordered[this_index]
      # py_run_string(paste0("sc.pl.heatmap(adata_", this_celltype_name, "_imbalance, var_names=tmp_genes, groupby='", DEG_key, "', layer='normalized', dendrogram=True, show_gene_labels=True, save='_", this_celltype_name, "_", DEG_key, "_", ii, ".pdf')"))
      py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=tmp_genes, groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_", ii, ".pdf')"))
      # py_run_string(paste0("sc.pl.heatmap(adata_", this_celltype_name, "_imbalance, var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, show_gene_labels=True, save='_", this_celltype_name, "_", DEG_key, "_ordered_", ii, ".pdf')"))
      py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_ordered_", ii, ".pdf')"))
      py_run_string("del tmp_genes")
    }
  }
  #########
  # Markers
  #########
  #py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=[], groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_T_helper.pdf')"))
  
  
}
for(DEG_key in DEG_keys){
  for(ii in 1:ceiling(length(krt_genes) / 80)){
    this_index = seq(80) + 80 * (ii - 1)
    this_index = this_index[this_index <= length(krt_genes)]
    py$tmp_genes = krt_genes[this_index]
    py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=tmp_genes, groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_krt_genes_", ii, ".pdf')"))
    py_run_string("del tmp_genes")
  }
}
######
DEG_key = "leiden_3_100"
basement_membrane = toupper(c("LAMA5", "ITGA3", "ITGB1", "ITGA6", "ITGB4"))
for(ii in c("basement_membrane")){
  py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=['", paste(eval(parse(text = ii)), collapse = "', '"), "'], groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_", ii, ".pdf')"))
  py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_3_imbalance', color=['", paste(eval(parse(text = ii)), collapse = "', '"), "'], layer='normalized', legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", ncols=10, frameon=False, save='_Signatures_", DEG_key, "_", ii, ".pdf')"))
}
######
for(ii in unlist(py$annotation_tmp)){
  py_run_string(paste0("clustering_plot(adata_", this_celltype_name, "_imbalance, '", ii, "', basis='umap_3_imbalance', size=", get_dot_size(cell_number), ", colorbar_loc=None, ncols=7, save='_detail_", ii, ".pdf')"))
}
#
cell_annotation[leiden_3_100 == "0"] = "KC - IFE Basal"
cell_annotation[leiden_3_100 == "1"] = "KC - IFE Spinous II"
cell_annotation[leiden_3_100 == "2"] = "KC - IFE Granular II"
cell_annotation[leiden_3_100 == "3"] = "KC - IFE Granular III"
cell_annotation[leiden_3_100 == "4"] = "KC - IFE Spinous IV"
cell_annotation[leiden_3_100 == "5"] = "KC - IFE Cycle - G2M"
cell_annotation[leiden_3_100 == "6"] = "KC - IFE Cycle - S"
cell_annotation[leiden_3_100 == "7"] = "KC - IFE Granular IV"
cell_annotation[leiden_3_100 == "8"] = "KC - IFE Spinous I"
cell_annotation[leiden_3_100 == "9"] = "KC - IFE Spinous V"
cell_annotation[leiden_3_100 == "10"] = "KC - Infundibulum"
cell_annotation[leiden_3_100 == "11"] = "KC - IFE Granular I"
cell_annotation[leiden_3_100 == "12"] = "KC - IFE Basal"
cell_annotation[leiden_3_100 == "13"] = "KC - IFE Cornified"
cell_annotation[leiden_3_100 == "14"] = "KC - IFE Basal"
cell_annotation[leiden_3_100 == "15"] = "KC - Spinous III"
cell_annotation[leiden_3_100 == "16"] = "KC - IRS/Cortex"
cell_annotation[leiden_3_100 == "17"] = "KC - IFE Granular - SLC35F1/SLC24A3/SOX5/TRAPPC5"
cell_annotation[leiden_3_100 == "18"] = "MIX - KC/Plasma"
cell_annotation[leiden_3_100 == "19"] = "KC - IFE Basal"
#
#
print(sum(is.na(cell_annotation)))
######
eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance$obs[['IntegratedAnnotation_3']] = cell_annotation")))
summary_list[[qc_state]][["IntegratedAnnotation_3"]] = list("AllCells" = table(cell_annotation))
py_run_string("annotation_tmp = annotation_tmp if np.isin('IntegratedAnnotation_2', annotation_tmp) else annotation_tmp + ['IntegratedAnnotation_2']")
py_run_string("annotation_tmp = annotation_tmp if np.isin('IntegratedAnnotation_3', annotation_tmp) else annotation_tmp + ['IntegratedAnnotation_3']")
py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_3_imbalance', color=annotation_tmp, legend_loc='on data', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=4, save='_annotation.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_3_imbalance', color=annotation_tmp, legend_loc='right margin', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=1, save='_annotation_1.pdf')"))
######
ii = "IntegratedAnnotation_3"
py_run_string(paste0("clustering_plot(adata_", this_celltype_name, "_imbalance, '", ii, "', basis='umap_3_imbalance', size=", get_dot_size(cell_number), ", colorbar_loc=None, ncols=7, save='_detail_", ii, ".pdf')"))
#########
# Markers
#########
eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance_hvg$obs[['IntegratedAnnotation_3']] = as.factor(cell_annotation)")))
py_run_string(paste0("sc.tl.dendrogram(adata_", this_celltype_name, "_imbalance_hvg, groupby='IntegratedAnnotation_3', use_rep='PCA_use')"))
py_run_string(paste0("adata_", this_celltype_name, "_imbalance.uns['dendrogram_IntegratedAnnotation_3'] = adata_", this_celltype_name, "_imbalance_hvg.uns['dendrogram_IntegratedAnnotation_3']"))
DEG_key = "IntegratedAnnotation_3"
this_path = paste0(py$this_output_dir, "/DEG_", DEG_key, ".rds")
DEG_leiden_path_list[[this_celltype_name]][[DEG_key]] = this_path
if(file.exists(this_path)){
  DEG_leiden_list[[this_celltype_name]][[DEG_key]] = readRDS(this_path)
}else{
  DEG_leiden_list[[this_celltype_name]][[DEG_key]] = DEG_Analysis(eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance"))), DEG_key, py$this_output_dir, prefix = paste0(DEG_key, "_"), top_DE_range = top_DE_range)
  saveRDS(DEG_leiden_list[[this_celltype_name]][[DEG_key]], this_path)
}
this_table = DEG_leiden_list[[this_celltype_name]][[DEG_key]][["rest"]][["table"]]
py_run_string(paste0("sc.tl.dendrogram(adata_", this_celltype_name, "_imbalance_hvg, groupby='", DEG_key, "', use_rep='PCA_use')"))
py_run_string(paste0("adata_", this_celltype_name, "_imbalance.uns['dendrogram_", DEG_key, "'] = adata_", this_celltype_name, "_imbalance_hvg.uns['dendrogram_", DEG_key, "']"))
this_clusters = unique(this_table$Cluster)
top_DEGs_ordered = c()
for(this_cluster in eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance$uns[['dendrogram_", DEG_key, "']][['categories_ordered']]")))){
  top_DEGs_ordered = c(top_DEGs_ordered, head(this_table[this_table$Cluster == this_cluster, ], 10)$Gene)
}
if(length(top_DEGs) <= max_gene_number){
  py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=[
              '", paste(top_DEGs_ordered, collapse = "', '"), "'
              ], groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_ordered.pdf')"))
}else{
  for(ii in 1:ceiling(length(top_DEGs_ordered) / max_gene_number)){
    this_index = seq(max_gene_number) + max_gene_number * (ii - 1)
    this_index = this_index[this_index <= length(top_DEGs_ordered)]
    py$tmp_genes_ordered = top_DEGs_ordered[this_index]
    py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=tmp_genes_ordered, groupby='", DEG_key, "', layer='normalized', dendrogram=True, save='", this_celltype_name, "_", DEG_key, "_ordered_", ii, ".pdf')"))
    py_run_string("del tmp_genes_ordered")
  }
}
py_run_string(paste0("sc.pl.dotplot(adata_", this_celltype_name, "_imbalance, var_names=[
              'H19', 'LOR', 'FLG', 'SPINK5', 'CALML5', 'KRT2', 'KRTDAP', 'KRT1', 'KRT10', 'KRT14', 'KRT5', 'KRT15',
              'KRT23', 'KRT77', 'KRT78', 'KRT80', 'CST6', 'CNFN', 'NCCRP1', 'SLURP1', # IFE Cornified - KRT23/KRT77/KRT78/KRT80/CST6/CNFN/NCCRP1/SLURP1
              'SLC35F1', 'SLC24A3', 'SOX5', 'TRAPPC5', # IFE Granular - SLC35F1/SLC24A3/SOX5/TRAPPC5
              'KRT35', 'KRT85', # KC - IRS/Cortex
              'JCHAIN', # MIX - KC/Plasma
              'GINS2', # KC - IFE Cycle - G2M
              'MKI67', # KC - IFE Cycle - S
              'MOXD1', 'TM4SF1' # KC - Infundibulum
              ], groupby='IntegratedAnnotation_3', layer='normalized', dendrogram=True, save='", this_celltype_name, "_IntegratedAnnotation_3_DEGs.pdf')"))


#######
# Write
#######
path_list[[this_celltype_name]][["IntegratedAnnotation_3"]] = paste0(py$this_output_dir, "/IntegratedAnnotation_3.h5ad")
py_run_string(paste0("adata_", this_celltype_name, "_imbalance.write('", path_list[[this_celltype_name]][["IntegratedAnnotation_3"]], "')"))
density_visualization(paste0("adata_", this_celltype_name, "_imbalance"), "umap_3_imbalance", "SS3_HumanSkin_20K_AnatomicalRegionLevel3", ".pdf", "density_color", "2", as.character(get_dot_size(cell_number)))

cell_annotation_SS3_HumanSkin_20K = cell_annotation[grep("SS3_HumanSkin_20K", names(cell_annotation))]
write.table(cell_annotation_SS3_HumanSkin_20K, file = paste0(py$this_output_dir, "/IntegratedAnnotation_3.tsv"), sep = "\t", col.names = F)

