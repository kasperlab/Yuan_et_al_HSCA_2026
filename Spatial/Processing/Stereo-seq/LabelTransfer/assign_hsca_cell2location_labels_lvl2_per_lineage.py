import scanpy as sc
import numpy as np
import pandas as pd

data_directory = '/datadrive/data1/axel/spatialdata'
adata = sc.read_h5ad(f'{data_directory}/clavel25_ctl_occ_nonadi.h5ad')

ref_label = 'integrated_annotation_L2'
lineages_and_celltypes = {'Fibroblast': ['Fibroblast'],
                          'Immune': ['Myeloid Cell', 'Lymphocyte', 'Mast Cell'],
                          'Vascular Endothelial Cell': ['Vascular Endothelial Cell', 'Lymphatic Endothelial Cell'],
                          'Melanocyte': ['Melanocyte'],
                          'Schwann Cell': ['Schwann Cell'],
                           'Mural Cell': ['Mural Cell'],
                           'Keratinocyte': ['Keratinocyte']}

for lineage, celltypes in lineages_and_celltypes.items():

    adata_lin = adata[adata.obs['cell_lineage'] == lineage].copy()
    
    lineage_label = lineage.replace(' ', '-')

    results_folder = f'./clavel25_ctl_occ_nonadi_cell2loc_{ref_label}_{lineage_label}/'

    # create paths and names to results folders for reference regression and cell2location models
    run_name = f'{results_folder}cell2location_map'

    celltype_abundances = []
    for sample in adata.obs['orig.ident'].unique():

        sample_file = f"{run_name}_{lineage}_{sample}_nonadi/sp.h5ad"
        adata_cell2loc_sample = sc.read_h5ad(sample_file)

        celltype_abundances.append(adata_cell2loc_sample.obsm['q05_cell_abundance_w_sf'])

    celltype_abundances_merged = pd.concat(celltype_abundances)
    celltype_abundances_merged = celltype_abundances_merged.loc[adata_lin.obs_names] # Paranoia

    adata_lin.obsm[ref_label + '_cell_abundance_w_sf'] = celltype_abundances_merged

    adata_lin

    adata_lin.write(f'{data_directory}/clavel25_ctl_occ_nonadi_{lineage_label}.h5ad', compression='gzip')
