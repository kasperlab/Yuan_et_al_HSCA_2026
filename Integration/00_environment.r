r_utilities = "/home/haoy/projects/Analysis/code_library/utilities.r"
py_utilities = "/home/haoy/projects/Analysis/code_library/utilities.py"
#py_deseq2 = "/home/haoy/projects/Analysis/code_library/deseq2.py"
py_visualization = "/home/haoy/projects/Analysis/code_library/visualization.py"
#py_scsea = "/home/haoy/projects/Analysis/code_library/scsea.py"
py_markers = "/home/haoy/projects/Analysis/code_library/markers.py"
#
source(r_utilities)
library(reticulate)
py_run_file(py_utilities, convert = F)
#py_run_file(py_deseq2, convert = F)
py_run_file(py_visualization, convert = F)
#py_run_file(py_scsea, convert = F)
py_run_file(py_markers, convert = F)
######
######
######
py$output_dir = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss"
py$output_dir_cache = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Integration/v67_gauss"
dir.create(py$output_dir, showWarnings = F, recursive = T)
######
######
loom_path_2 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/SS3_HumanSkin_20K/ALL/v15_gauss/all.loom"
metadata_path_2 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/SS3_HumanSkin_20K/ALL/v15_gauss/metadata_all_subtype.tsv"
loom_path_3 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Ning/ALL/v10_gauss/all.loom"
metadata_path_3 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Ning/ALL/v10_gauss/metadata_all.tsv"
loom_path_4 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Gaydosik_Fuschiotti_ClinicalCancerResearch_2019/ALL/v10_gauss/all.loom"
metadata_path_4 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Gaydosik_Fuschiotti_ClinicalCancerResearch_2019/ALL/v10_gauss/metadata_all.tsv"
loom_path_5 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Wang_Atwood_NatureCommunications_2020/ALL/v10_gauss/all.loom"
metadata_path_5 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Wang_Atwood_NatureCommunications_2020/ALL/v10_gauss/metadata_all.tsv"
loom_path_6 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Takahashi_Lowry_JournalofInvestigativeDermatology_2020/ALL/v10_gauss/all.loom"
metadata_path_6 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Takahashi_Lowry_JournalofInvestigativeDermatology_2020/ALL/v10_gauss/metadata_all.tsv"
loom_path_7 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Belote_Torres_NatureCellBiology_2021/ALL/v10_gauss/all.loom"
metadata_path_7 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Belote_Torres_NatureCellBiology_2021/ALL/v10_gauss/metadata_all.tsv"
loom_path_8 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Boldo_Lyko_CommunicationsBiology_2020/ALL/v10_gauss/all.loom"
metadata_path_8 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Boldo_Lyko_CommunicationsBiology_2020/ALL/v10_gauss/metadata_all.tsv"
loom_path_9 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Cheng_Cho_CellReports_2018/ALL/v10_gauss/all.loom"
metadata_path_9 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Cheng_Cho_CellReports_2018/ALL/v10_gauss/metadata_all.tsv"
loom_path_10 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Ji_Khavari_Cell_2020/ALL/v10_gauss/all.loom"
metadata_path_10 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Ji_Khavari_Cell_2020/ALL/v10_gauss/metadata_all.tsv"
loom_path_11 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Tabib_Lafyatis_NatureCommunications_2021/ALL/v10_gauss/all.loom"
metadata_path_11 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Tabib_Lafyatis_NatureCommunications_2021/ALL/v10_gauss/metadata_all.tsv"
loom_path_12 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/He_Guttman_JournalofAllergyandClinicalImmunology_2020/ALL/v10_gauss/all.loom"
metadata_path_12 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/He_Guttman_JournalofAllergyandClinicalImmunology_2020/ALL/v10_gauss/metadata_all.tsv"
loom_path_13 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Zou_Liu_DevelopmentalCell_2020/ALL/v10_gauss/all.loom"
metadata_path_13 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Zou_Liu_DevelopmentalCell_2020/ALL/v10_gauss/metadata_all.tsv"
loom_path_14 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Wiedemann_Andersen_CellReports_2023/ALL/v10_gauss/all.loom"
metadata_path_14 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Wiedemann_Andersen_CellReports_2023/ALL/v10_gauss/metadata_all.tsv"
loom_path_15 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Alkon_Stingl_JournalofAllergyandClinicalImmunology_2022/ALL/v10_gauss/all.loom"
metadata_path_15 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Alkon_Stingl_JournalofAllergyandClinicalImmunology_2022/ALL/v10_gauss/metadata_all.tsv"
loom_path_16 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Dunlap_Rao_JCIInsight_2022/ALL/v10_gauss/all.loom"
metadata_path_16 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Dunlap_Rao_JCIInsight_2022/ALL/v10_gauss/metadata_all.tsv"
loom_path_17 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Boothby_Rosenblum_Nature_2021/ALL/v10_gauss/all.loom"
metadata_path_17 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Boothby_Rosenblum_Nature_2021/ALL/v10_gauss/metadata_all.tsv"
loom_path_18 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Gao_Hu_CellDeathDisease_2021/ALL/v10_gauss/all.loom"
metadata_path_18 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Gao_Hu_CellDeathDisease_2021/ALL/v10_gauss/metadata_all.tsv"
loom_path_19 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Hughes_Shalek_Immunity_2020/HC/v10_gauss/all.loom"
metadata_path_19 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Hughes_Shalek_Immunity_2020/HC/v10_gauss/metadata_all.tsv"
loom_path_20 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Kim_Nagao_NatureMedicine_2020/ALL/v10_gauss/all.loom"
metadata_path_20 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Kim_Nagao_NatureMedicine_2020/ALL/v10_gauss/metadata_all.tsv"
loom_path_21 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Rindler_Brunner_MolecularCancer_2021/ALL/v10_gauss/all.loom"
metadata_path_21 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Rindler_Brunner_MolecularCancer_2021/ALL/v10_gauss/metadata_all.tsv"
loom_path_22 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020/ALL/v10_gauss/all.loom"
metadata_path_22 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Rojahn_Brunner_JournalofAllergyandClinicalImmunology_2020/ALL/v10_gauss/metadata_all.tsv"
loom_path_23 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Theocharidis_Bhasin_NatureCommunications_2022/ALL/v10_gauss/all.loom"
metadata_path_23 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Theocharidis_Bhasin_NatureCommunications_2022/ALL/v10_gauss/metadata_all.tsv"
loom_path_24 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Vorstandlechner_Mildner_NatureCommunications_2021/ALL/v10_gauss/all.loom"
metadata_path_24 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Vorstandlechner_Mildner_NatureCommunications_2021/ALL/v10_gauss/metadata_all.tsv"
loom_path_25 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Xu_Chen_Nature_2022/ALL/v10_gauss/all.loom"
metadata_path_25 = "/home/haoy/projects/MetaStudiesAnalysis/Integration/Output/Xu_Chen_Nature_2022/ALL/v10_gauss/metadata_all.tsv"
######
DEG_leiden_path_list = list()
knn_result_list = list()
knn_result_path_list = list()
pca_result_list = list()
pca_result_path_list = list()
cluster_name_list = list()
split_adata_name_list = list()
StudyID_list = list()
hvgene_list = list()
pca_cal_list = list()
mean_std_list = list()
other_cells_list = list()
DEG_leiden_list = list()
path_list = list()
key_list = list("divide" = list())
######
knn_result_path_list[["All"]] = paste0(py$output_dir, "/knn_result.rds")
pca_result_path_list[["All"]] = paste0(py$output_dir, "/pca_result.rds")
# Clustering parameters
hvg_number = 2000
pca_dim = 50
random_state = 0
n_neighbors_clustering = 30
n_neighbors_visualization = 15
leiden_resolutions = c(0.01, 0.02, 0.03, 0.05, 0.08, 0.1, 0.12, 0.15, 0.2, 0.3, 0.5, 0.8, 1, 1.2, 1.5, 2, 3)
cell_cycle_regress_out = F
# DE parameters
top_DE_range = 50
# Other parameters
chunk_size = 500000
max_number = 4000
py_run_string("density_color = 'YlOrRd'")
# Other settings
do_DEG = F
use_altas_querying = F
use_cache = T
######
X_gene_path = "/home/haoy/projects/Analysis/useful_files/X_Gene_Homo_sapiens.GRCh38.103.rds"
Y_gene_path = "/home/haoy/projects/Analysis/useful_files/Y_Gene_Homo_sapiens.GRCh38.103.rds"
py$X_Genes = readRDS(X_gene_path)
py$Y_Genes = readRDS(Y_gene_path)
alias_gene_c_seurat_path = "/home/haoy/projects/Analysis/useful_files/hgnc_complete_set_c_seurat.rds"
CellCyele_gene_path = "/home/haoy/projects/MetaStudiesAnalysis/Integration/regev_lab_cell_cycle_genes.txt"
alias_gene_c_seurat = readRDS(alias_gene_c_seurat_path)
######
# Following https://nbviewer.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
CellCyele_gene = gsub(" ", "", toupper(as.matrix(read.table(CellCyele_gene_path))[, 1]))
CellCyele_gene[CellCyele_gene == "MLF1IP"] = "CENPU" # MLF1IP is coded by CENPU
py$s_genes = CellCyele_gene[1:43]
py$g2m_genes = CellCyele_gene[44:length(CellCyele_gene)]
# Settings
py_run_string("sc.settings.verbosity = 3")
py_run_string("sc.logging.print_header()")
py_run_string("sc.settings.set_figure_params(dpi=600, dpi_save=600, facecolor='white', fontsize=10)")
py_run_string("sc.settings.autoshow = False")
py_run_string("figure_dir = output_dir + '/figures'")
py_run_string("sc.settings.figdir = figure_dir")
dir.create(py$figure_dir, showWarnings = F, recursive = T)
######
summary_list = list("Selection" = list())
