###########
# Mast Cell
###########
this_celltype = "Mast Cell"
this_celltype_name = gsub("_", "", gsub("[ -]", "", this_celltype), fixed = T)
py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs_supported_using_ImmuneCell")
dir.create(py$this_output_dir, showWarnings = F, recursive = T)
py_run_string("this_figure_dir = this_output_dir + '/figures'")
py_run_string("sc.settings.figdir = this_figure_dir")
qc_state = "before_QC"
max_gene_number = 200
######
py_run_string("adata_ImmuneCell = adata_imbalance[adata_imbalance.obs['IntegratedAnnotation_1'].isin(['Lymphocyte', 'Myeloid Cell', 'Mast Cell']), :].copy()")
py_run_string("adata_ImmuneCell.X = adata_ImmuneCell.layers['raw_counts'].copy()")
py_run_string("del adata_ImmuneCell.uns['log1p']")
py_run_string("print(adata_ImmuneCell.X)")
###
py_run_string("adata_MastCell = adata_ImmuneCell[adata_ImmuneCell.obs['IntegratedAnnotation_1'] == 'Mast Cell'].copy()")
adata_str = "adata_MastCell"
clustering_embedding = "PCA_ImmuneCell_MastCell"
visualization_embedding = "UMAP_ImmuneCell_MastCell"
# HVG Selection specific for Mast Cell
input_counts = eval(parse(text = paste0("t(as(py$", adata_str, "$layers['raw_counts'], 'CsparseMatrix'))")))
dimnames(input_counts) = eval(parse(text = paste0("list(py$", adata_str, "$var_names$values, py$", adata_str, "$obs_names$values)")))
selected_metadata = eval(parse(text = paste0("as.data.frame(apply(py$", adata_str, "$obs[unlist(py$merge_by)], 2, as.character))")))
rownames(selected_metadata) = eval(parse(text = paste0("py$", adata_str, "$obs_names$values")))
#
HVG_subset = seurat_hvg_selection(input_counts = input_counts, metadata = selected_metadata, hvg_number = hvg_number, verbose = F)
rm(input_counts)
rm(selected_metadata)
#
py$tmp_adata_hvg = PCA_calculation_python(adata = py$adata_ImmuneCell, random_state = random_state, merge_by = unlist(py$merge_by), hvg = HVG_subset, pca_dim = pca_dim, return_all = FALSE)
py_run_string(paste0(adata_str, ".obsm['", clustering_embedding, "'] = tmp_adata_hvg[adata_ImmuneCell.obs['IntegratedAnnotation_1'] == 'Mast Cell'].obsm['PCA_use'].copy()"))
py_run_string("del tmp_adata_hvg")
py_run_string("gc.collect()")
# Clustering
clustering_keys = pipeline_post_PCA(adata_str, clustering_embedding, visualization_embedding,
                                    n_neighbors_clustering, n_neighbors_visualization, leiden_resolutions,
                                    UMAP_min_dist=0.5,
                                    clustering_prefix = "_leiden_",
                                    random_state = random_state)
#
h5ad_path = paste0(py$this_output_dir, "/clustering.h5ad")
py_run_string(paste0("", adata_str, ".write('", h5ad_path, "')"))
