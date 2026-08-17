

adata_imbalance_annotated_path = paste0(py$output_dir_cache, "/IntegratedAnnotation_2.h5ad")
tryCatch({
  py_run_string(paste0("print(adata_imbalance)"))
}, error = function(x){
  py_run_string(paste0("adata_imbalance = ad.read_h5ad('", adata_imbalance_annotated_path, "')"))
  py_run_string("adata_imbalance.uns['log1p']['base'] = None")
  print("Loaded adata_imbalance from cache.")
})


py_run_string("print(adata_imbalance.obs['IntegratedAnnotation_1'].cat.categories)")
py_run_string("adata_imbalance.obs['Annotation_Level_1'] = adata_imbalance.obs['IntegratedAnnotation_1'].copy()")
py_run_string("new_order=['Keratinocyte', 'Fibroblast', 'Lymphocyte', 'Myeloid Cell', 'Mast Cell', 
                          'Melanocyte', 'Schwann Cell', 'Mural Cell', 'Vascular Endothelial Cell', 'Lymphatic Endothelial Cell' 
                          ]")
py_run_string("color_dict={'Keratinocyte': colorsys.hsv_to_rgb(0, 0.9, 0.5),
                           'Fibroblast': colorsys.hsv_to_rgb(0.3, 0.9, 0.5),
                           'Lymphocyte': colorsys.hsv_to_rgb(0.7, 0.9, 0.5),
                           'Myeloid Cell': colorsys.hsv_to_rgb(0.8, 0.7, 0.7),
                           'Mast Cell': colorsys.hsv_to_rgb(0.9, 0.5, 0.9),
                           'Melanocyte': colorsys.hsv_to_rgb(0.1, 0.7, 0.7),
                           'Schwann Cell': colorsys.hsv_to_rgb(0.2, 0.5, 0.9),
                           'Mural Cell': colorsys.hsv_to_rgb(0.4, 0.5, 0.9),
                           'Vascular Endothelial Cell': colorsys.hsv_to_rgb(0.5, 0.9, 0.5),
                           'Lymphatic Endothelial Cell': colorsys.hsv_to_rgb(0.6, 0.5, 0.9)
                          }")
py_run_string("adata_imbalance.obs['Annotation_Level_1'] = adata_imbalance.obs['Annotation_Level_1'].cat.reorder_categories(new_order)")
cell_number = py$adata_imbalance$n_obs
py$size = get_dot_size(cell_number)
py_run_string("sc.pl.embedding(adata_imbalance,
                               basis='X_umap_imbalance',
                               color=['Annotation_Level_1'],
                               title='Annotation_Level_1',
                               legend_loc='right margin',
                               legend_fontoutline=0.5,
                               legend_fontsize=6,
                               show=False,
                               ncols=5,
                               size=size,
                               frameon=False,
                               palette=color_dict,
                               add_outline=False,
                               save='_Annotation_Level_1.pdf'
)")

py_run_string("sc.pl.embedding(adata_imbalance,
                               basis='X_umap_imbalance',
                               color=['Annotation_Level_1'],
                               title='Annotation_Level_1',
                               legend_loc='on data',
                               legend_fontoutline=0.5,
                               legend_fontsize=6,
                               show=False,
                               ncols=5,
                               size=size,
                               frameon=False,
                               palette=color_dict,
                               add_outline=False,
                               save='_Annotation_Level_1_OD.pdf'
)")

py_run_string("sc.pl.embedding(adata_imbalance,
                               basis='X_umap_imbalance',
                               layer='normalized',
                               color=['SFN', 'COL1A1', 'CD3D', 'AIF1', 'CTSG', 'TYR', 'MPZ', 'ACTA2', 'PLVAP', 'TFF3'],
                               legend_loc='right margin',
                               legend_fontoutline=0.5,
                               legend_fontsize=6,
                               show=False,
                               ncols=5,
                               size=25,
                               frameon=False,
                               add_outline=False,
                               cmap=gene_highlight_cmap_gray_to_red,
                               save='_celltype_marker_genes_F1.pdf'
)")

DEG_key = "Annotation_Level_1"
py_run_string("marker_dict={'Keratinocyte': ['KRT14', 'KRT1', 'SFN', 'LY6D', 'DSC3'],
                            'Fibroblast': ['DCN', 'CFD', 'COL1A1', 'LUM', 'VCAN'],
                            'Lymphocyte': ['CD52', 'IL32', 'CXCR4', 'CD69', 'CD3D'],
                            'Myeloid Cell': ['AIF1', 'FCER1G', 'CD83', 'CD74', 'TYROBP'],
                            'Mast Cell': ['TPSAB1', 'HPGD', 'CTSG', 'TPSB2', 'SRGN'],
                            'Melanocyte': ['MLANA', 'TYRP1', 'DCT', 'PMEL', 'TYR'],
                            'Schwann Cell': ['MPZ', 'S100B', 'PLP1', 'CDH19', 'NRXN1'],
                            'Mural Cell': ['TAGLN', 'MYL9', 'ACTA2', 'TPM2', 'MYLK'],
                            'Vascular Endothelial Cell': ['PECAM1', 'AQP1', 'PLVAP', 'ACKR1', 'EMCN'],
                            'Lymphatic Endothelial Cell': ['TFF3', 'MMRN1', 'CLDN5', 'LYVE1', 'PROX1']
                          }")

py_run_string(paste0("fig = sc.pl.dotplot(adata_imbalance, var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=False, return_fig = True, show = False)"))
py_run_string("fig.add_totals()")
py_run_string("fig.savefig(f'{figure_dir}/dotplot_Annotation_Level_1_F1.pdf')")


py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(adata_imbalance, var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{figure_dir}/dotplot_adata_imbalance", "_", DEG_key, "_var.pdf', bbox_inches = 'tight')"))










py_run_string("sc.settings.figdir = figure_dir")



py_run_string("sc.pl.embedding(adata_imbalance,
                               basis='X_umap_imbalance',
                               color=['Annotation_Level_1'],
                               legend_loc='right margin',
                               legend_fontoutline=0.5,
                               legend_fontsize=6,
                               show=False,
                               ncols=1,
                               groups=['Lymphocyte', 'Mast Cell', 'Myeloid Cell'],
                               size=8,
                               frameon=False,
                               add_outline=False,
                               save='_IMMUNECELL_RM.pdf')")





py_run_string("sc.pl.embedding(adata_imbalance,
                               basis='X_umap_imbalance',
                               color=['Annotation_Level_1'],
                               legend_loc='right margin',
                               legend_fontoutline=0.5,
                               legend_fontsize=6,
                               show=False,
                               ncols=1,
                               groups=['Vascular Endothelial Cell', 'Lymphatic Endothelial Cell', 'Mural Cell'],
                               size=8,
                               frameon=False,
                               add_outline=False,
                               save='_Vessel_Associated_Cells_RM.pdf')")



# 20250409
# v67_gauss
py$rename_dict = list(
  "Keratinocyte" = "Keratinocytes",
  "Melanocyte" = "Melanocytes",
  "Schwann Cell" = "Schwann Cells",
  "Fibroblast" = "Fibroblasts",
  "Mural Cell" = "Mural Cells",
  "Vascular Endothelial Cell" = "Vascular Endothelial Cells",
  "Lymphatic Endothelial Cell" = "Lymphatic Endothelial Cells",
  "Erythrocyte" = "Erythrocytes",
  "Lymphocyte" = "Lymphocytes",
  "Myeloid Cell" = "Myeloid Cells",
  "Mast Cell" = "Mast Cells",
  "Plasma" = "Plasma Cells",
  "nan" = "Unknown"
)
py$index_name_dict = list(
  "Keratinocytes" = "KC",
  "Melanocytes" = "MEL",
  "Schwann Cells" = "SC",
  "Fibroblasts" = "FIB",
  "Mural Cells" = "MUR",
  "Vascular Endothelial Cells" = "VEC",
  "Lymphatic Endothelial Cells" = "LEC",
  "Erythrocytes" = "ERY",
  "Lymphocytes" = "LYM",
  "Myeloid Cells" = "MYE",
  "Mast Cells" = "MAST",
  "Plasma Cells" = "PC",
  "Unknown" = "NA"
)
py_run_string("color_dict_CT_FULL={'KC: Keratinocytes': colorsys.hsv_to_rgb(0, 0.9, 0.5),
                                   'MEL: Melanocytes': colorsys.hsv_to_rgb(0.1, 0.7, 0.7),
                                   'SC: Schwann Cells': colorsys.hsv_to_rgb(0.2, 0.5, 0.9),
                                   'FIB: Fibroblasts': colorsys.hsv_to_rgb(0.3, 0.9, 0.5),
                                   'MUR: Mural Cells': colorsys.hsv_to_rgb(0.4, 0.5, 0.9),
                                   'VEC: Vascular Endothelial Cells': colorsys.hsv_to_rgb(0.5, 0.9, 0.5),
                                   'LEC: Lymphatic Endothelial Cells': colorsys.hsv_to_rgb(0.6, 0.5, 0.9),
                                   'ERY: Erythrocytes': colorsys.hsv_to_rgb(0.85, 1, 1),
                                   'LYM: Lymphocytes': colorsys.hsv_to_rgb(0.7, 0.9, 0.5),
                                   'MYE: Myeloid Cells': colorsys.hsv_to_rgb(0.8, 0.7, 0.7),
                                   'MAST: Mast Cells': colorsys.hsv_to_rgb(0.9, 0.5, 0.9),
                                   'PC: Plasma Cells': colorsys.hsv_to_rgb(0.75, 1, 1),
                                   'NA: Unknown': colorsys.hsv_to_rgb(0, 0, 0.9)
}")

#
adata_str = "adata_imbalance"
#
py_run_string(paste0(adata_str, ".obs['AnnotationName'] = ", adata_str, ".obs['SeparateAnnotation_1'].astype(str).map({index: f'{name}' for index, name in rename_dict.items()})"))
py_run_string(paste0(adata_str, ".obs['AnnotationABBR'] = ", adata_str, ".obs['AnnotationName'].map({index: f'{name}' for index, name in index_name_dict.items()})"))
py_run_string(paste0(adata_str, ".obs['AnnotationFULL'] = ", adata_str, ".obs['AnnotationName'].map({index: f'{name}: {index}' for index, name in index_name_dict.items()})"))
py_run_string(paste0(adata_str, ".obs['AnnotationABBR_I'] = ", adata_str, ".obs['AnnotationABBR'].map({abbr: f'{index}:{abbr}' for index, abbr in enumerate(index_name_dict.values())})"))
py_run_string("order_AnnotationFULL = np.array([f'{name}: {index}' for index, name in index_name_dict.items()])")
py_run_string("order_AnnotationABBR_I = np.array([f'{index}:{abbr}' for index, abbr in enumerate(index_name_dict.values())])")
#
py_run_string(paste0(adata_str, ".obs['AnnotationFULL'] = ", adata_str, ".obs['AnnotationFULL'].cat.reorder_categories(order_AnnotationFULL[np.isin(order_AnnotationFULL, ", adata_str, ".obs['AnnotationFULL'])])"))
py_run_string(paste0(adata_str, ".obs['AnnotationABBR_I'] = ", adata_str, ".obs['AnnotationABBR_I'].cat.reorder_categories(order_AnnotationABBR_I[np.isin(order_AnnotationABBR_I, ", adata_str, ".obs['AnnotationABBR_I'])])"))
#
#
adata_str = "adata_imbalance"
x_axis = "AnnotationABBR_I"
y_axis = "StudyID_ordered"
#
py_run_string(paste0("count_matrix = pd.crosstab(", adata_str, ".obs['", y_axis, "'], ", adata_str, ".obs['", x_axis, "'])"))
py_run_string("count_matrix = count_matrix.loc[:, order_AnnotationABBR_I[np.isin(order_AnnotationABBR_I, count_matrix.columns)]]")
py_run_string("print(count_matrix)")
py_run_string("count_matrix_GEO = np.greater_equal(count_matrix.values, 1)")
py_run_string("normalized_matrix = pd.DataFrame(data=np.where(count_matrix_GEO, np.log10(np.where(count_matrix_GEO, count_matrix.values, 1)), np.nan), index=count_matrix.index, columns=count_matrix.columns)")
py_run_string("print(normalized_matrix)")
py_run_string("fig, ax = plt.subplots(figsize=(len(count_matrix.columns) * 0.5 + 3, len(count_matrix.index) * 0.5 + 3))")
py_run_string("sns.heatmap(normalized_matrix, square=True, annot=count_matrix, fmt='d', vmin=0, vmax=np.minimum(np.where(count_matrix_GEO, normalized_matrix, 0).max(), 6), cmap='YlOrRd', annot_kws={'fontsize':6}, ax=ax)")
py_run_string(paste0("fig.savefig(f'{figure_dir}/heatmap_", adata_str, "_", x_axis, "_", y_axis, "_ALL.pdf', bbox_inches = 'tight')"))
py_run_string("print(normalized_matrix)")
py_run_string("matplotlib.pyplot.close('all')")
py_run_string("gc.collect()")
#
py_run_string("print(adata_imbalance.obs['StudyID_ordered'].value_counts())")
py_run_string("print(adata_imbalance.obs['SeparateAnnotation_1'].astype(str).value_counts())")
py_run_string("print(adata_imbalance[adata_imbalance.obs['StudyID_ordered'] == 'Ning'].obs['SeparateAnnotation_1'].astype(str).value_counts())")
py_run_string("print(adata_imbalance[adata_imbalance.obs['StudyID_ordered'] == 'Zou_Liu_DevelopmentalCell_2020'].obs['SeparateAnnotation_1'].astype(str).value_counts())")
#
SeparateAnnotation_1 = py$adata_imbalance$obs[["SeparateAnnotation_1"]]
StudyID = py$adata_imbalance$obs[["StudyID"]]






py_run_string("color_list_tab20b = [matplotlib.colors.rgb2hex(matplotlib.colormaps['tab20b'](x)) for x in range(20)]")
py_run_string("color_list_tab20c = [matplotlib.colors.rgb2hex(matplotlib.colormaps['tab20c'](x)) for x in range(20)]")
py_run_string("print(color_list_tab20b)")
py_run_string("print(color_list_tab20c)")
py_run_string("color_list = color_list_tab20b + color_list_tab20c[:4]")
py_run_string("print(color_list)")
#
adata_str = "adata_imbalance"
x_axis = "AnnotationABBR_I"
y_axis = "StudyID_ordered"
#
py_run_string(paste0("count_matrix = pd.crosstab(", adata_str, ".obs['", y_axis, "'], ", adata_str, ".obs['", x_axis, "'])"))
py_run_string("count_matrix = count_matrix.loc[:, order_AnnotationABBR_I[np.isin(order_AnnotationABBR_I, count_matrix.columns)]]")
py_run_string("count_matrix.columns = count_matrix.columns.map({f'{index}:{abbr}': abbr for index, abbr in enumerate(index_name_dict.values())})")
py_run_string("print(count_matrix)")



py_run_string("from matplotlib import font_manager as fm")
py_run_string("fm.fontManager.addfont('/home/haoy/projects/Analysis/fonts/ARIAL.TTF')")
py_run_string("from matplotlib.gridspec import GridSpec")
py$figure_size = c(8.27 * 2, 11.69 * 2)
SW = 210
SH = 297
BT = 20
BB = 17
BL = 20
BR = 20
CW = SW - BL - BR
CH = SH - BT - BB
#
py_run_string("fig = plt.figure(figsize=figure_size)")
py_run_string(paste0("gs = GridSpec(", as.integer(CH), ", ", as.integer(CW), ", figure=fig, left=", BL/SW, ", bottom=", BB/SH, ", right=", 1 - BR/SW, ", top=", 1 - BT/SH, ", wspace=0, hspace=0, width_ratios=None, height_ratios=None)"))
######
gene_names = eval(parse(text = paste0("py$", adata_str, "$var_names$values")))
#
col_number = 6
WI = -2
HI = 5
start_LC = 0
start_TC = 10
WH = 30
#
for(ii in (seq(12) - 1)){
  this_row = ii %/% col_number
  this_col = ii %% col_number
  this_LC = start_LC + this_col * (WH + WI)
  this_TC = start_TC + this_row * (WH + HI)
  #
  py_run_string(paste0("ax = fig.add_subplot(gs[", as.integer(this_TC), ":", as.integer(this_TC + WH), ", ", as.integer(this_LC), ":", as.integer(this_LC + WH), "])"))
  #
  py_run_string(paste0("this_column = count_matrix.columns[", ii, "]"))
  py_run_string("sorted_series_raw = count_matrix[this_column].sort_values(ascending=False)")
  py_run_string("this_sum = sorted_series_raw.values.sum()")
  py_run_string("this_cutoff = this_sum * 0.05")
  py_run_string("sorted_series = sorted_series_raw[sorted_series_raw >= this_cutoff]")
  py_run_string("rest_cell_number = this_sum - sorted_series.sum()")
  if(py$rest_cell_number > 0){
    py_run_string("sorted_series['Others'] = this_sum - sorted_series.sum()")
  }
  py_run_string("this_label = [f'Other\\n' if np.isnan(y) else f'#{int(y)}\\n' for x, y in zip(sorted_series.values, sorted_series.index.map({name: index + 1 for index, name in enumerate(adata_imbalance.obs['StudyID_ordered'].cat.categories)}))]")
  py_run_string("this_color = ['#E8E8E8' if np.isnan(x) else color_list[int(x)] for x in sorted_series.index.map({name: index for index, name in enumerate(adata_imbalance.obs['StudyID_ordered'].cat.categories)})]")
  py_run_string("ax.pie(sorted_series.values, explode=[0.1] * len(sorted_series), labels=this_label, colors=this_color, autopct=None, labeldistance=0.75, counterclock=False, startangle=90, radius=1, textprops={'fontsize': 8, 'ha': 'center', 'va': 'center'}, rotatelabels=False)")
  py_run_string("ax.pie(sorted_series.values, explode=[0.05] * len(sorted_series), labels=None, colors='w', autopct=lambda x: '\\n' + str(round(x)) + '%', pctdistance=1.6, counterclock=False, startangle=90, radius=0.5, textprops={'fontsize': 6, 'ha': 'center', 'va': 'center'}, rotatelabels=False)")
  py_run_string(paste0("fig.text(", this_LC + BL + 0.5 * WH, " / ", SW, ", 1 - ", this_TC + BT + 0.5 * WH, "/",  SH,", this_column + '\\n', fontstyle='normal', fontfamily='Arial', fontsize=16, fontweight='normal', ha='center', va='center', color='black')"))
  py_run_string(paste0("fig.text(", this_LC + BL + 0.5 * WH, " / ", SW, ", 1 - ", this_TC + BT + 0.5 * WH, "/",  SH,", '\\n\\n' + str(this_sum) + '\\ncells', fontstyle='normal', fontfamily='Arial', fontsize=12, fontweight='normal', ha='center', va='center', color='black')"))
}
#
py_run_string(paste0("fig.savefig(f'{figure_dir}/AnnotationABBR_I_donut.pdf')"))
py_run_string("gc.collect()")
py_run_string("plt.close('all')")


