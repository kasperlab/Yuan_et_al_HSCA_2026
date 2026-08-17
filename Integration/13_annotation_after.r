########
# Others
########
py_run_string("sc.settings.figdir = figure_dir")
cell_number = py$adata_imbalance$n_obs
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['OriginalAnnotation', 'SeparateAnnotation_1', 'SeparateAnnotation_2', 'IntegratedAnnotation_2_v35', 'IntegratedAnnotation_2_v65', 'IntegratedAnnotation_3_v65'], legend_loc='right margin', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=1, save='_annotation_2.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['KRT14', 'FBLN1', 'ACTA2', 'MLANA', 'MPZ', 'PECAM1', 'LYVE1', 'TPSB2', 'LYZ', 'CD3D'], legend_loc='right margin', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=5, save='_celltype_marker.pdf')"))
py_run_string(paste0("print(adata_imbalance.obs['StudyID'].cat.categories)"))
StudyID = py$adata_imbalance$obs[["StudyID"]]
py$adata_imbalance$obs[["StudyID_ordered"]] = factor(StudyID, levels = StudyID[!duplicated(StudyID)])
density_visualization("adata_imbalance", "umap_imbalance", "StudyID_ordered", "imbalance.pdf", "density_color", "6", as.character(get_dot_size(cell_number)))
density_visualization("adata_imbalance", "umap_imbalance", "StudyID_ordered", "imbalance_5x5.pdf", "density_color", "5", as.character(get_dot_size(cell_number)))
ii = "StudyID"
py_run_string(paste0("clustering_plot(adata_imbalance, '", ii, "', basis='X_umap_imbalance', size=16, colorbar_loc=None, ncols=7, save='_detail_", ii, "_big_dots.pdf')"))
py_run_string(paste0("clustering_plot(adata_imbalance, '", ii, "', basis='X_umap_imbalance', size=8, colorbar_loc=None, ncols=7, save='_detail_", ii, "_dot_size_8.pdf')"))
py_run_string(paste0("clustering_plot(adata_imbalance, '", ii, "', basis='X_umap_imbalance', size=4, colorbar_loc=None, ncols=7, save='_detail_", ii, "_dot_size_4.pdf')"))

py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['IntegratedAnnotation_1', 'CXCL14', 'KITLG', 'LAMA4'], legend_loc='right margin', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=1, save='_markers_20231106.pdf')"))
py_run_string("annotation_grant = adata_imbalance.obs['IntegratedAnnotation_1'].astype(str).copy()")
py_run_string("annotation_grant[annotation_grant.astype(str) == 'Keratinocyte'] = 'Other Keratinocytes'")
py_run_string("annotation_grant[adata_imbalance.obs['PreviousIntegratedAnnotation_2'].astype(str) == 'Endothelial-related Pericyte'] = 'Endothelial-related Pericyte'")
py_run_string("annotation_grant[adata_imbalance.obs['IntegratedAnnotation_2'].astype(str) == 'Hair Folicle Basal - 1'] = 'Hair Folicle Basal - 1'")
py_run_string("annotation_grant[adata_imbalance.obs['IntegratedAnnotation_2'].astype(str) == 'IFE Basal - 1'] = 'IFE Basal - 1'")
py_run_string("annotation_grant[adata_imbalance.obs['IntegratedAnnotation_2'].astype(str) == 'IFE Basal - 2'] = 'IFE Basal - 2'")
py_run_string("annotation_grant[adata_imbalance.obs['IntegratedAnnotation_2'].astype(str) == 'IFE Spinous - 1'] = 'IFE Spinous - 1'")
py_run_string("del adata_imbalance.obs['annotation_grant']")
py_run_string("adata_imbalance.obs['annotation_grant'] = annotation_grant.copy()")
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['annotation_grant', 'CXCL14'], legend_loc='right margin', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=1, save='_markers_20231108.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['KITLG'], vmax=3, legend_loc='right margin', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=1, save='_markers_20231108_1.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['LAMA4'], vmax=5, legend_loc='right margin', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=1, save='_markers_20231108_2.pdf')"))

py_run_string(paste0("sc.pl.stacked_violin(adata_imbalance, var_names=['CXCL14', 'KITLG', 'LAMA4'], groupby='annotation_grant', frameon=False, swap_axes=True, save='_markers_20231108_violin.pdf')"))



celltype_1 = as.character(unique(py$adata_imbalance$obs[["IntegratedAnnotation_1"]]))
celltype_1 = c("Mural Cell", "Lymphatic Endothelial Cell")
for(this_celltype in celltype_1){
  this_celltype_name = gsub("[ -]", "", this_celltype)
  ##############
  #Configuration
  ##############
  #py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name)
  py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs")
  dir.create(py$this_output_dir, showWarnings = F, recursive = T)
  py_run_string("this_figure_dir = this_output_dir + '/figures'")
  py_run_string("sc.settings.figdir = this_figure_dir")
  ###
  tryCatch({
    py_run_string(paste0("print(adata_", this_celltype_name, "_imbalance)"))
  }, error = function(x){
    py_run_string(paste0("adata_", this_celltype_name, "_imbalance = ad.read_h5ad('", py$output_dir_cache, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs/clustering.h5ad')"))
    print(paste0("Loaded ", this_celltype_name, " from cache."))
  })
  
  ###
  cell_number = eval(parse(text = paste0("py$adata_", this_celltype_name, "_imbalance$n_obs")))
  py_run_string(paste0("sc.pl.embedding(adata_", this_celltype_name, "_imbalance, 'umap_2_imbalance', color=['PreviousIntegratedAnnotation_2', 'CXCL14', 'KITLG', 'LAMA4'], legend_loc='right margin', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=1, save='_markers_20231106_", this_celltype_name, ".pdf')"))
}


##########
# Bodysite
##########
bodysite = as.character(eval(parse(text = paste0("py$adata_imbalance$obs[['AnatomicalRegionLevel3']]"))))
names(bodysite) = eval(parse(text = paste0("py$adata_imbalance$obs_names$values")))
study_id = eval(parse(text = paste0("py$adata_imbalance$obs[['StudyID']]")))
bodysite[study_id != "SS3_HumanSkin_20K"] = NA
bodysite_slot = "SS3_HumanSkin_20K_AnatomicalRegionLevel3"
eval(parse(text = paste0("py$adata_imbalance$obs[['", bodysite_slot, "']] = bodysite")))
#py_run_string(paste0("tmp_this_factor = pd.Series(adata_imbalance.obs['", bodysite_slot, "'].astype(str).values)"))
#for(jj in tolower(c("NA", "nan", "unknown", "True"))){
#  py_run_string(paste0("tmp_this_factor[tmp_this_factor.str.lower() == '", jj, "'] = pd.NA")) # Only pd.NA can make <NA> in Pandas (NA in R), np.nan will make NaN and None will make None.
#}
#py_run_string(paste0("del adata_imbalance.obs['", bodysite_slot, "']"))
#py_run_string(paste0("adata_imbalance.obs['", bodysite_slot, "'] = tmp_this_factor.values"))
#py_run_string("del tmp_this_factor")
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['", bodysite_slot, "'], size=8, ncols=1, legend_loc='right margin', legend_fontsize=6, frameon=False, save='_SS3_HumanSkin_20K_AnatomicalRegionLevel3.pdf')"))
py_run_string(paste0("clustering_plot(adata_imbalance, '", bodysite_slot, "', basis='umap_imbalance', size=8, ncols=2, save='_detail_", bodysite_slot, ".pdf')"))
density_visualization(paste0("adata_imbalance"), "umap_imbalance", bodysite_slot, ".pdf", "density_color", "2", "8")


py_run_string("adata_SS3_HumanSkin_20K = adata_imbalance[adata_imbalance.obs['StudyID'] == 'SS3_HumanSkin_20K', :].copy()")
py_run_string(paste0("del adata_SS3_HumanSkin_20K.uns['SeparateAnnotation_1_colors']"))
py$tmp_color = as.character(py$adata_imbalance$uns[["IntegratedAnnotation_1_colors"]])[-9]
py_run_string(paste0("adata_SS3_HumanSkin_20K.uns['SeparateAnnotation_1_colors'] = tmp_color"))
py_run_string(paste0("del tmp_color"))
#py_run_string(paste0("sc.tl.score_genes(adata_SS3_HumanSkin_20K, ['MPZ', 'SOX10', 'SOX2'], ctrl_size=50, gene_pool=None, n_bins=25, score_name='SchwannCellScore', random_state=0, copy=False, use_raw=None)"))
py_run_string(paste0("sc.tl.score_genes(adata_SS3_HumanSkin_20K, ['MPZ', 'SOX2'], ctrl_size=50, gene_pool=None, n_bins=25, score_name='SchwannCellScore', random_state=0, copy=False, use_raw=None)"))
py_run_string(paste0("sc.pl.embedding(adata_SS3_HumanSkin_20K, 'X_umap_imbalance', color=['SeparateAnnotation_1', 'IntegratedAnnotation_1', 'SchwannCellScore', 'SOX10', 'MPZ', 'SOX2'], size=8, ncols=1, legend_loc='right margin', legend_fontsize=6, frameon=False, save='_SS3_HumanSkin_20K_annotation_1.pdf')"))
py_run_string(paste0("adata_SS3_HumanSkin_20K.write('", paste0(py$output_dir, "/adata_SS3_HumanSkin_20K_1.h5ad"), "')"))


# 20240425
py_run_string(paste0("adata_imbalance = ad.read_h5ad('/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss/adata_imbalance_annotated_1.h5ad')"))

py_run_string("sc.settings.figdir = figure_dir")
cell_number = py$adata_imbalance$n_obs

py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['SOX10', 'NCMAP', 'GFRA3', 'BCHE'], legend_loc='right margin', legend_fontsize=6, size=", get_dot_size(cell_number), ", frameon=False, ncols=1, save='_markers_SC.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['SOX10', 'NCMAP', 'GFRA3', 'BCHE'], legend_loc='right margin', legend_fontsize=6, size=8, frameon=False, ncols=4, save='_markers_SC_big.pdf')"))

# 20240503
IntegratedAnnotation_1 = py$adata_imbalance$obs[["IntegratedAnnotation_1"]]
StudyID = py$adata_imbalance$obs[["StudyID"]]
cell_id = py$adata_imbalance$obs_names$values


table(StudyID[IntegratedAnnotation_1 == "Lymphocyte"])
table(StudyID[IntegratedAnnotation_1 == "Myeloid Cell"])
table(StudyID[IntegratedAnnotation_1 == "Mast Cell"])


reference_mask_1 = IntegratedAnnotation_1 %in% c("Lymphocyte", "Myeloid Cell", "Mast Cell")
reference_mask_2 = StudyID %in% c("Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020", "Rindler_Brunner_MolecularCancer_2021", "Tabib_Lafyatis_NatureCommunications_2021")

query_mask_1 = IntegratedAnnotation_1 %in% c("Lymphocyte", "Myeloid Cell", "Melanocyte")
query_mask_2 = StudyID %in% c("Alkon_Stingl_JournalofAllergyandClinicalImmunology_2022", "He_Guttman_JournalofAllergyandClinicalImmunology_2020", "Ji_Khavari_Cell_2020")

SCSEA_ref = rep("N", py$adata_imbalance$n_obs)
SCSEA_ref[reference_mask_1 & reference_mask_2] = "Y"
py$SCSEA_ref = SCSEA_ref
py_run_string("adata_SCSEA_ref = adata_imbalance[np.array(SCSEA_ref) == 'Y', :].copy()")
py_run_string("adata_SCSEA_ref.write(output_dir + '/adata_SCSEA_ref.h5ad')")

SCSEA_que = rep("N", py$adata_imbalance$n_obs)
SCSEA_que[query_mask_1 & reference_mask_2] = "Y"
py$SCSEA_que = SCSEA_que
py_run_string("adata_SCSEA_que = adata_imbalance[np.array(SCSEA_que) == 'Y', :].copy()")
py_run_string("adata_SCSEA_que.write(output_dir + '/adata_SCSEA_que.h5ad')")


#20240920
StudyID = py$adata_imbalance$obs[["StudyID"]]
py$adata_imbalance$obs[["StudyID_ordered"]] = factor(StudyID, levels = StudyID[!duplicated(StudyID)])
density_visualization("adata_imbalance", "umap_imbalance", "StudyID_ordered", "size_8.pdf", "density_color", "6", 8)
density_visualization("adata_imbalance", "umap_imbalance", "StudyID_ordered", "size_16.pdf", "density_color", "6", 16)
density_visualization("adata_imbalance", "umap_imbalance", "StudyID_ordered", "size_25.pdf", "density_color", "6", 25)

py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['KRT14', 'FBLN1', 'ACTA2', 'MLANA', 'MPZ', 'PECAM1', 'LYVE1', 'TPSB2', 'LYZ', 'CD3D'], cmap=gene_highlight_cmap, legend_loc='right margin', legend_fontsize=6, size=8, frameon=False, ncols=5, save='_celltype_marker_size_8.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['KRT14', 'FBLN1', 'ACTA2', 'MLANA', 'MPZ', 'PECAM1', 'LYVE1', 'TPSB2', 'LYZ', 'CD3D'], cmap=gene_highlight_cmap, legend_loc='right margin', legend_fontsize=6, size=8, frameon=False, ncols=5, save='_celltype_marker_size_16.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['KRT14', 'FBLN1', 'ACTA2', 'MLANA', 'MPZ', 'PECAM1', 'LYVE1', 'TPSB2', 'LYZ', 'CD3D'], cmap=gene_highlight_cmap, legend_loc='right margin', legend_fontsize=6, size=8, frameon=False, ncols=5, save='_celltype_marker_size_25.pdf')"))

py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['KRT10', 'COL1A1', 'TAGLN', 'DCT', 'SOX2', 'VWF', 'TFF3', 'TPSAB1', 'HLA-DRA', 'CD3D'], cmap=gene_highlight_cmap, legend_loc='right margin', legend_fontsize=6, size=8, frameon=False, ncols=5, save='_celltype_marker_C_size_8.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['KRT10', 'COL1A1', 'TAGLN', 'DCT', 'SOX2', 'VWF', 'TFF3', 'TPSAB1', 'HLA-DRA', 'CD3D'], cmap=gene_highlight_cmap, legend_loc='right margin', legend_fontsize=6, size=8, frameon=False, ncols=5, save='_celltype_marker_C_size_16.pdf')"))
py_run_string(paste0("sc.pl.embedding(adata_imbalance, 'X_umap_imbalance', color=['KRT10', 'COL1A1', 'TAGLN', 'DCT', 'SOX2', 'VWF', 'TFF3', 'TPSAB1', 'HLA-DRA', 'CD3D'], cmap=gene_highlight_cmap, legend_loc='right margin', legend_fontsize=6, size=8, frameon=False, ncols=5, save='_celltype_marker_C_size_25.pdf')"))


#
py_run_string("sc.settings.figdir = figure_dir")
ii = "IntegratedAnnotation_1"
py_run_string(paste0("clustering_plot(adata_imbalance, '", ii, "', basis='X_umap_imbalance', size=8, ncols=5, save='_detail_", ii, "_size_8.pdf')"))
py_run_string(paste0("clustering_plot(adata_imbalance, '", ii, "', basis='X_umap_imbalance', size=16, ncols=5, save='_detail_", ii, "_size_16.pdf')"))
py_run_string(paste0("clustering_plot(adata_imbalance, '", ii, "', basis='X_umap_imbalance', size=25, ncols=5, save='_detail_", ii, "_size_25.pdf')"))





