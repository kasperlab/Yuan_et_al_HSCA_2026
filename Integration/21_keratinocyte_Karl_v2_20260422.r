##############
# Keratinocyte
##############
this_celltype = "Keratinocyte"
this_celltype_name = gsub("_", "", gsub("[ -]", "", this_celltype), fixed = T)
py$this_output_dir = paste0(py$output_dir, "/", this_celltype_name, "_", hvg_number, "_hvgs_", pca_dim, "_pcs")
dir.create(py$this_output_dir, showWarnings = F, recursive = T)
qc_state = "before_QC"
max_gene_number = 200
######
adata_str = "adata_Keratinocyte_filtered"
clustering_embedding = "integrated_UMAP_L2"
visualization_embedding = "integrated_UMAP_L2"
######
py$this_output_dir = paste0(py$this_output_dir, "/KC_KARL")
dir.create(py$this_output_dir, showWarnings = F, recursive = T)
py_run_string("this_figure_dir = this_output_dir + '/figures'")
py_run_string("sc.settings.figdir = this_figure_dir")
######
tryCatch({
  stop()
  py_run_string(paste0("print(", adata_str, ")"))
}, error = function(x){
  py_run_string(paste0("adata_Keratinocyte_imbalance = ad.read_h5ad('", py$this_output_dir, "/clustering.h5ad')"))
  py_run_string("del adata_Keratinocyte_imbalance.uns")
  print("Loaded adata_Keratinocyte_imbalance from cache.")
})
tryCatch({
  stop()
  py_run_string("print(adata_KC_karl)")
}, error = function(x){
  py_run_string("adata_KC_karl = ad.read_h5ad('/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v70_gauss/adata_KC.h5ad')")
  print("Loaded adata_KC_karl from cache.")
})
######
py_run_string(paste0(adata_str, " = adata_Keratinocyte_imbalance[np.isin(adata_Keratinocyte_imbalance.obs_names, adata_KC_karl.obs_names)]"))
py_run_string(paste0("adata_KC_karl_filtered = adata_KC_karl[", adata_str, ".obs_names]"))
######
py_run_string(paste0(adata_str, ".obsm['integrated_UMAP_L2'] = adata_KC_karl_filtered.obsm['integrated_UMAP_L2'].copy()"))
py_run_string(paste0(adata_str, ".obs['Annotation_Level_2_FULL'] = adata_KC_karl_filtered.obs['Annotation_Level_2_FULL'].copy()"))
py_run_string(paste0(adata_str, ".uns = {'Annotation_Level_2_FULL_colors': adata_KC_karl_filtered.uns['Annotation_Level_2_FULL_colors'].copy()}"))
py_run_string("annotation_tmp = copy.deepcopy(annotation_1)")
cell_number = eval(parse(text = paste0("py$", adata_str, "$n_obs")))
py$dotsize = max(c(8, get_dot_size(cell_number)))
############
# Signatures
############
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[
              'LHX2', 'LGR5', 'KRT14', 'KRT15', 'KRT5', 'KRT1', 'KRT10'
              ], layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize,
                     ncols=4, frameon=False, save='_Signatures_IFE.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[
              'CXCL14', 'COL17A1', 'KRT5', 'KRT15', 'KRT16', 'TK1'
              ], layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize,
                     ncols=3, frameon=False, save='_Signatures_Basal.pdf')"))
genename = eval(parse(text = paste0("py$", adata_str)))$var_names$values
krt_genes = genename[grep("^KRT.+", genename)]
#py$rpsl_genes = genename[grep("^RP[LS].+", genename)]
for(ii in 1:ceiling(length(krt_genes) / 80)){
  this_index = seq(80) + 80 * (ii - 1)
  this_index = this_index[this_index <= length(krt_genes)]
  py$tmp_genes = krt_genes[this_index]
  py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=tmp_genes, layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=10, frameon=False, save='_krt_genes_", ii, ".pdf')"))
  py_run_string("del tmp_genes")
}
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[
              'KRT77', 'LOR', 'FLG2', 'ELOVL1', 'LCE1A'
              ], layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize,
                     ncols=3, frameon=False, save='_Signatures_differentiate.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[
              'PTTG1', 'CDC20', 'RRM2', 'HELLS', 'UHRF1', 'ASS1', 'COL17A1', 'POSTN', 'KRT19',
              'GJB2', 'KRT6A', 'KRT16', 'CCND1', 'DEFB1', 'FXYD3', 'CALML5', 'ZNF750', 'SPINK5'
              ], layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize,
                     ncols=3, frameon=False, save='_Wang_Atwood_NatureCommunications_2020_Fig_1D.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[
              'KRT14', 'KRT5', 'CDH3', 'CDH1', 'KRT1', 'KRT10', 'DSG1', 'DSC1', 'KRT2',
              'IVL', 'TGM3'
              ], layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize,
                     ncols=3, frameon=False, save='_Wang_Atwood_NatureCommunications_2020_Fig_2B.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[
              'MGST1', 'SCD', 'POSTN', 'CD34', 'KRT14', 'MT2A', 'KRT10', 'PTGS1',
              'KRT17', 'KRT79', 'KRT6A', 'KRT75'
              ], layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize,
                     ncols=2, frameon=False, save='_Joost_Kasper_CellSystems_2016_Fig_1C.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=[
              'KRT5', 'KRT15', 'KRT1', 'UBE2C', 'FASN', 'ATP1B1', 'GJB2', 'KRT6B'
              ], layer='normalized', legend_loc='on data', legend_fontsize=6, size=dotsize,
                     ncols=4, frameon=False, save='_He_Yassky_TheJournalofAllergyandClinicalImmunology_2020_Fig_2A.pdf')"))
######
py_run_string("annotation_tmp = annotation_tmp if np.isin('Annotation_Level_2_FULL', annotation_tmp) else annotation_tmp + ['Annotation_Level_2_FULL']")
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=annotation_tmp, legend_loc='on data', legend_fontsize=6, size=dotsize, frameon=False, ncols=4, save='_annotation.pdf')"))
py_run_string(paste0("sc.pl.embedding(", adata_str, ", '", visualization_embedding, "', color=annotation_tmp, legend_loc='right margin', legend_fontsize=6, size=dotsize, frameon=False, ncols=1, save='_annotation_1.pdf')"))
######
ii = "Annotation_Level_2_FULL"
py_run_string(paste0("clustering_plot(", adata_str, ", '", ii, "', basis='", visualization_embedding, "', size=dotsize, colorbar_loc=None, ncols=7, save='_detail_", ii, ".pdf')"))
#########
# Markers
#########
py_run_string(paste0("sc.pl.dotplot(", adata_str, ", var_names=[
              'KRT1', 'KRT10', 'LY6D', 'RHOV', 'SBSN', 'SPRR1B', 'SULT2B1', 'MT1G', # #6, #8, #5, #4, #0
              'CST6', 'FLG', 'SLURP1', # #19
              'SEMA3C', 'KRT23', 'GJB2', # #14
              'MMP7', 'S100P', 'KRT77', # #22
              'S100A7', 'S100A8', 'S100A9', 'CSTB', 'TMEM45A', # #11
              'KRT6B', 'ACTN1', 'SOSTDC1', # #20
              'APOD', # #18
              'SFRP1', 'FRZB', 'DIO2', # #13
              'C1QTNF12', 'MOXD1', # #12
              'KRT19', 'KRT18', 'KRT7', 'AQP5', 'DCD', 'MUCL1', # #16
              'MGST1', 'SAA1', 'APOC1', # #17
              'VIM', 'HBB', # #21, #18
              'PCNA', 'ASS1', 'LYPD3', # #3, #10
              'POSTN', 'CXCL14', 'DST', 'COL17A1', 'KRT15', 'KRT14', 'KRT5', # #1, #2, #15
              'GINS2', 'MCM5', # #7
              'TOP2A', 'MKI67', # #9
              ], groupby='Annotation_Level_2_FULL', layer='normalized', dendrogram=False, save='", this_celltype_name, "_Annotation_Level_2_FULL_DEGs.pdf')"))
# DEG
DEG_key = "Annotation_Level_2_FULL"
for(used_table in c("top_p_by_score", "top_p_by_FC", "top_p_by_P")){
  this_path = paste0(py$this_output_dir, "/DEG_", adata_str, "_" , DEG_key, ".rds")
  if(file.exists(this_path)){
    DEG_result = readRDS(this_path)
  }else{
    DEG_result = DEG_Analysis(eval(parse(text = paste0("py$", adata_str))), DEG_key, py$this_output_dir, prefix = paste0(adata_str, "_" , DEG_key, "_"), top_DE_range = top_DE_range)
    saveRDS(DEG_result, this_path)
  }
  this_table = DEG_result[["rest"]][[used_table]]
  this_clusters = unique(this_table$Cluster)
  top_DEGs = c()
  py_run_string("DEG_dict = {}")
  for(this_cluster in this_clusters){
    top_DEGs = c(top_DEGs, head(this_table[this_table$Cluster == this_cluster, ], 10)$Gene)
    py_run_string(paste0("DEG_dict['", this_cluster, "'] = ['", paste(head(this_table[this_table$Cluster == this_cluster, ], 10)$Gene, collapse = "', '"), "']"))
  }
  if(length(top_DEGs) <= max_gene_number){
    py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=DEG_dict, groupby='", DEG_key, "', layer='normalized', n_genes=10, swap_axes = False, highlight_genes = DEG_dict, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
    py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", used_table, ".pdf', bbox_inches = 'tight')"))
  }else{
    for(ii in 1:ceiling(length(top_DEGs) / max_gene_number)){
      this_index = seq(max_gene_number) + max_gene_number * (ii - 1)
      this_index = this_index[this_index <= length(top_DEGs)]
      py$tmp_genes = top_DEGs[this_index]
      py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=tmp_genes, groupby='", DEG_key, "', layer='normalized', n_genes=10, swap_axes = False, highlight_genes = DEG_dict, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
      py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_", ii, "_", used_table, ".pdf', bbox_inches = 'tight')"))
      py_run_string("del tmp_genes")
    }
  }
}
py_run_string("marker_genes_short = {
 '0: IFE_Basal_1': ['RETREG1','IGFBP3', 'KRT14', 'KRT5', 'SYT8'],
 '1: IFE_Basal_2': ['POSTN', 'LAMB4', 'KRT31', 'KRT15', 'DST'],
 '2: Infu_basal_1': ['C1QTNF12', 'MOXD1', 'EGR2', 'CDH13', 'GPX4'],
 '3: Infu_basal_2': ['S100A9', 'IL11RA', 'LAMC2', 'TAPBP', 'TNFAIP8'],
 '4: IFE_interferon': ['IFI44L', 'MX1', 'IFI6', 'IFI27', 'ISG15'],
 '5: IFE_Cycling_S': ['GINS2', 'HELLS', 'MCM6', 'MCM5', 'PCNA'],
 '6: IFE_Cycling_G2M_1': ['TK1', 'SPC25', 'RRM2', 'H4C3', 'CDK1'],
 '7: IFE_Cycling_G2M_2': ['DLGAP5', 'TROAP', 'ASPM', 'MKI67', 'CENPE'],
 '8: IFE_Spinous_1': ['MT1H'],
 '9: IFE_Spinous_2': ['MT1X', 'CCND1', 'DUSP23'],
 '10: IFE_Spinous_2.1': ['MT1X', 'MT1G', 'ADRB2', 'DSC3'],
 '11: IFE_Granular_1': ['RHOB', 'HES1', 'JUNB', 'ZFP36', 'MYC'],
 '12: IFE_Granular_1.1': ['KRT10', 'DMKN', 'SFN', 'RPLP1', 'MT-CO2'],
 '13: IFE_Granular_2': ['CHP2', 'KCNK7', 'KRTDAP', 'LAMP1', 'PYCARD'],
 '14: IFE_Granular_2.1': ['STX11', 'ZNF703', 'DCAF13', 'FOSL1', 'CHD1'],
 '15: IFE_Granular_3': ['PRSS3', 'KCNK7', 'CLEC2A'],
 '16: IFE_terminal_1': ['ACER1', 'SPRR1B', 'KRT2', 'THEM5', 'SCEL'],
 '17: IFE_terminal_2': ['BPIFC', 'PYDC1', 'CST6', 'FLG', 'LCE1C'],
 '18: Infu_diff': ['S100A7', 'S100A8', 'S100A9', 'SERPINB4', 'HBEGF'],
 '19: uHF': ['CARD18', 'PTN', 'ADGRL3', 'CYP27A1', 'DACH1'],
 '20: HF_1': ['KRT75', 'PDZRN3', 'KRT6C', 'KRT6B', 'KRT17'],
 '21: HF_2': ['DIO2', 'FRZB', 'SFRP1', 'LMCD1', 'DKK3'],
 '22: HF_3': ['LHX2', 'LGR5', 'BNC2', 'BCAM', 'BGN'],
 '23: HF_Anagen': ['KRT25', 'KRT35', 'KRT85', 'LEF1', 'BAMBI'],
 '24: SG_1': ['CYP4F8', 'SAA1', 'TINAGL1', 'SERPINF1', 'IL1R2'],
 '25: SG_2': ['THRSP', 'AWAT2', 'KRT79', 'ACSBG1', 'CIDEA'],
 '26: Gland_Channel_1': ['FST', 'TFCP2L1', 'TNC', 'FGF7', 'LBH'],
 '27: Gland_Channel_2': ['MGST1', 'CDA', 'KRT23', 'GJB6', 'SOSTDC1'],
 '28: Gland_Channel_3': ['MMP7', 'S100P', 'WFDC3', 'KRT77', 'KRT6B'],
 '29: Gland_Sweat_1': ['NCALD', 'CA6', 'SNORC', 'S100A1', 'AQP5'],
 '30: Gland_Sweat_2': ['CLDN3', 'SERHL2', 'TSPAN8', 'AZGP1', 'KRT19'],
 '31: Melanocytes': ['TYR', 'PLP1', 'DCT', 'TYRP1', 'MLANA'],
 '32: MIX_KC/Imm': ['PTPRC', 'IL32', 'SRGN', 'CD52', 'CD74'],
 '33: MIX_KC/FIB': ['VIM', 'COL3A1', 'COL1A2', 'LUM', 'CFD'],
 '34: MIX_KC/vSM': ['MYLK', 'TAGLN', 'MYL9', 'ACTA2', 'TPM2'],
 '35: MIX_KC/EC': ['PECAM1', 'EGFL7', 'RAMP2', 'CLDN5', 'GNG11'],
}")
py_run_string("marker_genes_long_1 = {
 '0: IFE_Basal_1': ['RETREG1', 'GDPD2', 'ASPN', 'WNT3', 'ASS1', 'SOX6', 'IGFBP3', 'KRT14', 'SYT8', 'KRT5', 'TGFBI', 'COL17A1', 'ITGB1', 'ARPC2', 'LSP1'],
 '1: IFE_Basal_2': ['POSTN', 'LAMB4', 'KRT31', 'KRT15', 'DST', 'IL18', 'NFKBIA', 'ERRFI1'],
 '2: Infu_basal_1': ['C1QTNF12', 'ABI3BP', 'MOXD1', 'EGR2', 'CDH13', 'CCN2', 'GPX4'],
 '3: Infu_basal_2': ['CCL20', 'S100A9', 'IQCG', 'HAS3', 'IL11RA', 'LAMC2', 'ARID5B', 'TAPBP', 'EDN1', 'TNFAIP8', 'GADD45A'],
 '4: IFE_interferon': ['IFI44L', 'MX1', 'IFI6', 'IFI44', 'IFI27', 'IFITM1', 'ISG15', 'IFIT1', 'OAS1', 'HERC6', 'PARP14'],
 '5: IFE_Cycling_S': ['E2F1', 'TK1', 'CLSPN', 'GINS2', 'HELLS', 'DHFR', 'CHEK1', 'MCM6', 'MCM5', 'TCF19', 'PCNA'],
 '6: IFE_Cycling_G2M_1': ['SPC25', 'RRM2', 'AURKB', 'PBK', 'CDCA3', 'KIF23', 'CDCA5', 'SGO1', 'CCNA2', 'MELK', 'TOP2A', 'MKI67', 'MND1', 'H4C3', 'PCLAF', 'STMN1', 'CENPW', 'HMGB2', 'CDK1', 'TUBB', 'MT2A', 'NUSAP1', 'CKS1B', 'PRC1'],
 '7: IFE_Cycling_G2M_2': ['DLGAP5', 'TROAP', 'NEK2', 'ASPM', 'HMMR', 'MKI67', 'CEP55', 'NUF2', 'CCNB2', 'CENPF', 'CDKN3', 'CENPE'],
 '8: IFE_Spinous_1': ['MT1H'],
 '9: IFE_Spinous_2': ['MT1X', 'CCND1', 'DUSP23'],
 '10: IFE_Spinous_2.1': ['CCND1', 'MT1X', 'MT1G', 'ADRB2', 'DSC3', 'GATM'],
 '11: IFE_Granular_1': ['NRARP', 'ZNF750', 'RHOB', 'HES1', 'JUNB', 'ZFP36', 'MYC', 'ARL4D', 'DNAJB1', 'DUSP1', 'HSPA1B', 'HSPA1A'],
 '12: IFE_Granular_1.1': ['KRT1', 'KRT10', 'DMKN', 'LY6D', 'DSP', 'MT1X', 'SFN', 'PERP', 'RPS26', 'LGALS3', 'RPLP1', 'MT-CO3', 'MT-CO2', 'NEAT1'],
 '13: IFE_Granular_2': ['CHP2', 'SBSN', 'KCNK7', 'KRTDAP', 'LAMP1', 'PYCARD'],
 '14: IFE_Granular_2.1': ['STX11', 'ZNF703', 'DCAF13', 'FOSL1', 'CHD1', 'ZPR1', 'IRF2BP2', 'PLIN3', 'NSG1', 'ELOC'],
 '15: IFE_Granular_3': ['PRSS3', 'KCNK7', 'CLEC2A'],
 '16: IFE_terminal_1': ['ACER1', 'SLURP1', 'FAM25A', 'SPRR1B', 'DEGS2', 'CLIC3', 'KRT2', 'THEM5', 'SPINK5', 'FABP5', 'CTNNBIP1', 'SCEL', 'CD24'],
 '17: IFE_terminal_2': ['BPIFC', 'SERPINB3', 'C1orf68', 'ATP6V1C2', 'ALOXE3', 'SDR9C7', 'NCCRP1', 'PYDC1', 'CST6', 'KRT78', 'CLIC3', 'SMIM5', 'FLG', 'CSTA', 'MIDN', 'CLDN4', 'LDLR', 'ORMDL1', 'ID4', 'LCE1C', 'LCE3D']
}")
py_run_string("marker_genes_long_2 = {
 '18: Infu_diff': ['S100A7', 'S100A8', 'S100A9', 'SERPINB4', 'HBEGF', 'FGFBP1', 'CSTB', 'CRABP2'],
 '19: uHF': ['SUSD2', 'CARD18', 'PTN', 'ADGRL3', 'CYP27A1', 'DACH1'],
 '20: HF_1': ['CTNND2', 'CBLN2', 'KRT75', 'VTCN1', 'PDZRN3', 'FXYD6', 'KRT6C', 'CHST2', 'KRT6B', 'CREB5', 'SOSTDC1'],
 '21: HF_2': ['DIO2', 'LHX2', 'WIF1', 'TCEAL2', 'FRZB', 'ISM1', 'SFRP1', 'KRT17', 'BCAM', 'LMCD1', 'DKK3', 'S100A6', 'LGALS1'],
 '22: HF_3': ['LHX2', 'SHISA2', 'LGR5', 'TENM2', 'BNC2', 'TNC', 'SOX4', 'BGN'],
 '23: HF_Anagen': ['KRT25', 'KRT35', 'KRT28', 'KRT85', 'DLX2', 'NRP2', 'LEF1', 'BAMBI', 'UNC5B', 'HOXC13', 'DLX3', 'KIF21A', 'KRTAP11-1', 'ODC1'],
 '24: SG_1': ['HSD11B1', 'CYP4F8', 'SAA1', 'NNAT', 'WFDC2', 'TINAGL1', 'SERPINF1', 'IL1R2', 'LY6E', 'MGST1'],
 '25: SG_2': ['SGK2', 'THRSP', 'AWAT2', 'KRT79', 'HAO2', 'FADS2', 'FAR2', 'ACSBG1', 'FBP1', 'TLCD4', 'CIDEA', 'ALOX15B', 'APOC1', 'CLSTN3', 'FDPS', 'ASCL1'],
 '26: Gland_Channel_1': ['IGFL2-AS1', 'PNLIPRP3', 'FST', 'SEMA3C', 'TFCP2L1', 'TNC', 'ECRG4', 'FGF7', 'GPX2', 'LBH', 'ATP1B1'],
 '27: Gland_Channel_2': ['MGST1', 'CDA', 'RSPO1', 'KRT77', 'KRT23', 'GJB6', 'SEMA3C', 'GJB2', 'SOSTDC1'],
 '28: Gland_Channel_3': ['CFTR', 'RHCG', 'MMP7', 'GPR12', 'S100P', 'WFDC3', 'MAL', 'KRT77', 'ATP6V1B1', 'KRT6B', 'CLIC3'],
 '29: Gland_Sweat_1': ['ELF5', 'NCALD', 'TESC', 'STAC2', 'CA6', 'SNORC', 'S100A1', 'MRAS', 'PPP1R1B', 'CLDN10', 'CAMK2N1', 'PDK3', 'COX7A1', 'KRT8', 'KRT7', 'AQP5'],
 '30: Gland_Sweat_2': ['SPDEF', 'OBP2B', 'GJC3', 'HMGCS2', 'CLDN3', 'LRRC26', 'ELAPOR1', 'SERHL2', 'TSPAN8', 'SCGB2A1', 'IQGAP2', 'AZGP1', 'KRT19'],
 '31: Melanocytes': ['TYR', 'TRPM1', 'PLP1', 'DCT', 'GPM6B', 'TYRP1', 'MLANA'],
 '32: MIX_KC/Imm': ['PTPRC', 'CORO1A', 'LCP1', 'CD69', 'CD37', 'CD53', 'GMFG', 'IL32', 'CXCR4', 'SRGN', 'CREM', 'CD52', 'CD74'],
 '33: MIX_KC/FIB': ['VIM', 'COL3A1', 'COL1A2', 'COL1A1', 'CFH', 'FBLN2', 'MFAP4', 'LUM', 'COL6A2', 'SPARC', 'CCDC80', 'SERPING1', 'CFD'],
 '34: MIX_KC/vSM': ['MYLK', 'TAGLN', 'MYL9', 'TIMP3', 'ACTA2', 'TPM2', 'SPARCL1', 'RGS16', 'GEM', 'TPM1'],
 '35: MIX_KC/EC': ['PECAM1', 'EGFL7', 'RAMP2', 'CLDN5', 'ACKR1', 'GNG11', 'SPARCL1', 'EMCN'],
}")
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_genes_short, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale='var', swap_axes=False, highlight_genes=None, highlight_params={'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line=False, cluster_line_highlight_params={'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_S_var.pdf', bbox_inches='tight')"))
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_genes_long_1, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale='var', swap_axes=False, highlight_genes=None, highlight_params={'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line=False, cluster_line_highlight_params={'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_L1_var.pdf', bbox_inches='tight')"))
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_genes_long_2, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale='var', swap_axes=False, highlight_genes=None, highlight_params={'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line=False, cluster_line_highlight_params={'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_L2_var.pdf', bbox_inches='tight')"))
#
py_run_string("marker_dict = {
  'IFE': ['KRT14', 'KRT5', 'KRT15', 'KRT1', 'KRT10'],
  'IFE_Basal': ['POSTN', 'IGFBP3', 'ASPN'],
  'IFE/Infu Basal': ['CCL2'],
  'Infu': ['S100A7', 'S100A8', 'S100A9'],
  'Infu_basal': ['C1QTNF12', 'MOXD1', 'CCN2', 'CCL19', 'LTF'],
  'Interferon': ['ISG15', 'MX1'],
  'IFE_Cycling': ['TK1', 'HELLS', 'MKI67'],
  'IFE_Spinous': ['MT1H', 'MT1G'],
  'IFE_Granular': ['NRARP', 'RHOB', 'SBSN', 'STX11', 'KCNK7', 'DSC1'],
  'IFE_terminal': ['ACER1', 'BPIFC', 'FLG', 'LCE1C', 'LOR'],
  'Infu_diff': ['SERPINB4'],
  'MIX(MEL)': ['TYR', 'MLANA'],
  'MIX(LYM)': ['PTPRC', 'CD3D'],
  'MIX(FIB)': ['COL1A1', 'LUM'],
  'MIX(MUR)': ['TAGLN', 'ACTA2'],
  'MIX(VEC)': ['PECAM1', 'ACKR1'],
}")
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(", adata_str, ", var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale='var', swap_axes=False, highlight_genes=None, highlight_params={'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line=False, cluster_line_highlight_params={'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_", adata_str, "_", DEG_key, "_var.pdf', bbox_inches='tight')"))
#######
# Write
#######
h5ad_path = paste0(py$this_output_dir, "/Annotation_Level_2_FULL.h5ad")
py_run_string(paste0("", adata_str, ".write('", h5ad_path, "')"))
#
# Further analysis
#
py_run_string(paste0("print(", adata_str, ".obs['Annotation_Level_2_FULL'].cat.categories)"))
py_run_string(paste0("name_index_dict = {x: x.split(': ')[0] for x in ", adata_str, ".obs['Annotation_Level_2_FULL'].cat.categories}"))
py_run_string("print(name_index_dict)")
py_run_string(paste0(adata_str, ".obs['Annotation_Level_2_renumber'] = ", adata_str, ".obs['Annotation_Level_2_FULL'].map(name_index_dict)"))
#
py_run_string(paste0("color_dict = {x.split(': ')[0]: y for x, y in zip(", adata_str, ".obs['Annotation_Level_2_FULL'].cat.categories, ", adata_str, ".uns['Annotation_Level_2_FULL_colors'])}"))
py_run_string("print(color_dict)")
py_run_string(paste0("color_dict_FULL = {x: y for x, y in zip(", adata_str, ".obs['Annotation_Level_2_FULL'].cat.categories, ", adata_str, ".uns['Annotation_Level_2_FULL_colors'])}"))
py_run_string("print(color_dict_FULL)")
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
py_run_string("adata_KC_no_HFG = adata_Keratinocyte_filtered[np.logical_or(adata_Keratinocyte_filtered.obs['Annotation_Level_2_renumber'].astype(np.int) < 19, adata_Keratinocyte_filtered.obs['Annotation_Level_2_renumber'].astype(np.int) > 30)]")
py_run_string("print(adata_KC_no_HFG.obs['Annotation_Level_2_FULL'].value_counts())")
py_run_string(paste0("dot_fig = plot_dotplot_with_annotations(adata_KC_no_HFG, var_names=marker_dict, groupby='", DEG_key, "', layer='normalized', dendrogram=False, n_genes=10, standard_scale = 'var', swap_axes = False, highlight_genes = None, highlight_params = {'facecolor': 'yellow', 'edgecolor': 'none', 'alpha': 0.5}, highlight_cluster_line = False, cluster_line_highlight_params = {'color': 'lightgreen', 'alpha': 0.2, 'zorder': -5},)"))
py_run_string(paste0("dot_fig.savefig(f'{this_figure_dir}/dotplot_adata_KC_no_HFG", "_", DEG_key, "_var.pdf', bbox_inches = 'tight')"))
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
py_run_string("matplotlib.pyplot.close('all')")
#
py$index_name_dict = list(
  "0: IFE_Basal_1" = "0: IFE_Basal_1",
  "1: IFE_Basal_2" = "1: IFE_Basal_2",
  "2: Infu_basal_1" = "2: Infu_basal_1",
  "3: Infu_basal_2" = "3: Infu_basal_2",
  "4: IFE_interferon" = "4: IFE_interferon",
  "5: IFE_Cycling_S" = "5: IFE_Cycling_S",
  "6: IFE_Cycling_G2M_1" = "6: IFE_Cycling_G2M_1",
  "7: IFE_Cycling_G2M_2" = "7: IFE_Cycling_G2M_2",
  "8: IFE_Spinous_1" = "8: IFE_Spinous_1",
  "9: IFE_Spinous_2" = "9: IFE_Spinous_2",
  "10: IFE_Spinous_2.1" = "10: IFE_Spinous_2.1",
  "11: IFE_Granular_1" = "11: IFE_Granular_1",
  "12: IFE_Granular_1.1" = "12: IFE_Granular_1.1",
  "13: IFE_Granular_2" = "13: IFE_Granular_2",
  "14: IFE_Granular_2.1" = "14: IFE_Granular_2.1",
  "15: IFE_Granular_3" = "15: IFE_Granular_3",
  "16: IFE_terminal_1" = "16: IFE_terminal_1",
  "17: IFE_terminal_2" = "17: IFE_terminal_2",
  "18: Infu_diff" = "18: Infu_diff",
  "19: uHF" = "HF&Gland",
  "20: HF_1" = "HF&Gland",
  "21: HF_2" = "HF&Gland",
  "22: HF_3" = "HF&Gland",
  "23: HF_Anagen" = "HF&Gland",
  "24: SG_1" = "HF&Gland",
  "25: SG_2" = "HF&Gland",
  "26: Gland_Channel_1" = "HF&Gland",
  "27: Gland_Channel_2" = "HF&Gland",
  "28: Gland_Channel_3" = "HF&Gland",
  "29: Gland_Sweat_1" = "HF&Gland",
  "30: Gland_Sweat_2" = "HF&Gland",
  "31: Melanocytes" = "19: MIX - MEL",
  "32: MIX_KC/Imm" = "20: MIX - LYM",
  "33: MIX_KC/FIB" = "21: MIX - FIB",
  "34: MIX_KC/vSM" = "22: MIX - MUR",
  "35: MIX_KC/EC" = "23: MIX - VEC"
)
py_run_string(paste0(adata_str, ".obs['Annotation_L2_FULL'] = ", adata_str, ".obs['Annotation_Level_2_FULL'].map({index: f'{name}' for index, name in index_name_dict.items()})"))
py_run_string(paste0("print(", adata_str, ".obs['Annotation_L2_FULL'].value_counts())"))
py_run_string(paste0("name_index_dict = {x: x.split(': ')[0] for x in ", adata_str, ".obs['Annotation_L2_FULL'].astype('category').cat.categories}"))
py_run_string("print(name_index_dict)")
py_run_string(paste0(adata_str, ".obs['Annotation_L2_renumber'] = ", adata_str, ".obs['Annotation_L2_FULL'].map(name_index_dict)"))
###
py_run_string(paste0("color_dict_FULL = {'0: IFE_Basal_1': '#0fc702ff', '1: IFE_Basal_2': '#1cfc0dff', '2: Infu_basal_1': '#01665eff', '3: Infu_basal_2': '#35978fff', '4: IFE_interferon': '#238b45ff', '5: IFE_Cycling_S': '#6bba18ff', '6: IFE_Cycling_G2M_1': '#42c908ff', '7: IFE_Cycling_G2M_2': '#009900ff', '8: IFE_Spinous_1': '#a2ff00ff', '9: IFE_Spinous_2': '#ffe600ff', '10: IFE_Spinous_2.1': '#eaff4aff', '11: IFE_Granular_1': '#fee08bff', '12: IFE_Granular_1.1': '#fc9a3aff', '13: IFE_Granular_2': '#fc4e2aff', '14: IFE_Granular_2.1': '#fd8d3cff', '15: IFE_Granular_3': '#cc3d10ff', '16: IFE_terminal_1': '#bd0026ff', '17: IFE_terminal_2': '#800026ff', '18: Infu_diff': '#80cdc1ff', '19: MIX - MEL': '#8b4513ff', '20: MIX - LYM': '#d3d3d3ff', '21: MIX - FIB': '#bfbfbfff', '22: MIX - MUR': '#a6a6a6ff', '23: MIX - VEC': '#8c8c8cff', 'HF&Gland': '#ffc8f0ff'}"))
py_run_string("print(color_dict_FULL)")
###
py_run_string(paste0("tmp_adata = ad.AnnData(obs=", adata_str, ".obs.copy(), obsm={'", visualization_embedding, "': ", adata_str, ".obsm['", visualization_embedding, "'].copy()})"))
py_run_string("np.random.seed(0)")
py_run_string("tmp_adata_random = tmp_adata[np.random.permutation(list(range(tmp_adata.n_obs))), :]")
py_run_string("del tmp_adata")
py_run_string("gc.collect()")
py_run_string(paste0("sc.pl.embedding(tmp_adata_random,
                                      basis='", visualization_embedding, "',
                                      color=['Annotation_L2_FULL'],
                                      legend_loc='none',
                                      legend_fontoutline=0.5,
                                      legend_fontsize=6,
                                      show=False,
                                      ncols=5,
                                      size=8,
                                      frameon=False,
                                      palette=color_dict_FULL,
                                      add_outline=False,
                                      save='_Annotation_L2_NL.png'
)"))
py_run_string(paste0("fig = sc.pl.embedding(tmp_adata_random,
                                            basis='", visualization_embedding, "',
                                            color=['Annotation_L2_renumber'],
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
                                      color=['Annotation_L2_FULL'],
                                      legend_loc='right margin',
                                      legend_fontsize=6,
                                      show=False,
                                      size=8,
                                      add_outline=False,
                                      palette=color_dict_FULL,
                                      ax=fig.axes[0],
                                      save='_Annotation_L2_RM.pdf'
)"))
py_run_string("del tmp_adata_random")
py_run_string("gc.collect()")








# 20251002
py_run_string(paste0("adata_Keratinocyte_filtered = ad.read_h5ad('", py$this_output_dir, "/IntegratedAnnotation_2_20250531.h5ad')"))
py_run_string("adata_Keratinocyte_HVG = ad.read_h5ad('/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss/Keratinocyte_2000_hvgs_50_pcs/clustering_hvg.h5ad')")
py_run_string("adata_Keratinocyte_filtered.obsm['PCA_use'] = adata_Keratinocyte_HVG[adata_Keratinocyte_filtered.obs_names].obsm['PCA_use'].copy()")


DEG_key = "Annotation_L2_FULL"
clustering_embedding = "PCA_use"
py_run_string(paste0(adata_str, ".obs['", DEG_key, "'] = ", adata_str, ".obs['", DEG_key, "'].astype('category')"))
py_run_string(paste0("print(", adata_str, ".obs['", DEG_key, "'].cat.categories.size)"))
py_run_string(paste0("sc.tl.dendrogram(", adata_str, ", groupby='", DEG_key, "', use_rep='", clustering_embedding, "')"))


py_run_string(paste0("fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(", adata_str, ".obs['", DEG_key, "'].cat.categories.size * 0.3, 20))"))
py_run_string(paste0("sc.pl.dendrogram(", adata_str, ", '", DEG_key, "', ax=axs)"))
py_run_string(paste0("fig.savefig(f'{this_figure_dir}/DendrogramPlot_", adata_str, "_", DEG_key, ".pdf')"))
py_run_string("plt.close('all')")
py_run_string("gc.collect()")





py_run_string("sc.pl.embedding(adata_KC_karl_filtered, 'integrated_UMAP_L2', color=['leiden_2_150', 'leiden_2_200', 'leiden_2_300', 'leiden_2_400', 'Annotation_Level_2_FULL'], legend_loc='on data', legend_fontsize=6, size=dotsize, frameon=False, ncols=3, save='_KC_clustering_annotation.pdf')")



# 20260422
py_run_string(paste0(adata_str, " = ad.read_h5ad('", py$this_output_dir, "/IntegratedAnnotation_2_20250531.h5ad')"))
py_run_string("annotation_tmp = copy.deepcopy(annotation_1)")
cell_number = eval(parse(text = paste0("py$", adata_str, "$n_obs")))
py$dotsize = max(c(8, get_dot_size(cell_number)))
#
py_run_string(paste0(adata_str, "_subset = ", adata_str, "[", adata_str, ".obs['Annotation_Level_2_renumber'].isin(np.array([2, 3, 18]).astype(str))].copy()"))
py_run_string(paste0("color_dict = {x.split(': ')[0]: y for x, y in zip(", adata_str, ".obs['Annotation_Level_2_FULL'].cat.categories, ", adata_str, ".uns['Annotation_Level_2_FULL_colors'])}"))
py_run_string("print(color_dict)")
visualization_embedding = "integrated_UMAP_L2"
py_run_string(paste0("sc.pl.embedding(", adata_str, "_subset, '", visualization_embedding, "', color='Annotation_Level_2_renumber', layer='normalized', cmap=gene_highlight_cmap, legend_loc='on data', legend_fontsize=6, size=dotsize, ncols=6, frameon=False, save='_subset_Annotation_Level_2_renumber.pdf')"))
py$adata_Keratinocyte_filtered_subset












