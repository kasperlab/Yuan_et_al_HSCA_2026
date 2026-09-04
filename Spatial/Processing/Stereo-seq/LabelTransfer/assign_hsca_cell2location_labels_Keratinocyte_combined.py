import scanpy as sc
import numpy as np
import pandas as pd

data_directory = '/datadrive/data1/axel/'
# adata = sc.read_h5ad(f'{data_directory}/clavel25_ctl_occ_nonadi.h5ad')

ref_label = 'integrated_annotation_combined'

lineage = 'Keratinocyte'    
lineage_label = lineage.replace(' ', '-')

results_folder = f'./clavel25_ctl_occ_nonadi_cell2loc_{lineage_label}_{ref_label}/'

adata_lin = sc.read_h5ad(f'{data_directory}spatialdata/clavel25_ctl_occ_nonadi_{lineage_label}.h5ad')

# create paths and names to results folders for reference regression and cell2location models
run_name = f'{results_folder}cell2location_map'

celltype_abundances = []
for sample in adata_lin.obs['orig.ident'].unique():

    sample_file = f"{run_name}_{sample}/sp.h5ad"
    adata_cell2loc_sample = sc.read_h5ad(sample_file)

    celltype_abundances.append(adata_cell2loc_sample.obsm['q05_cell_abundance_w_sf'])

celltype_abundances_merged = pd.concat(celltype_abundances)
# celltype_abundances_merged = celltype_abundances_merged.loc[adata_lin.obs_names]
adata_lin = adata_lin[celltype_abundances_merged.index]

adata_lin.obsm[ref_label + '_cell_abundance_w_sf'] = celltype_abundances_merged

adata_lin.write(f'{data_directory}spatialdata/clavel25_ctl_occ_nonadi_{lineage_label}_{ref_label}.h5ad', compression='gzip')
