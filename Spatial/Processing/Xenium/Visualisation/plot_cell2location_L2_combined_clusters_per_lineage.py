
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
import os
from pathlib import Path

import matplotlib as mpl
from copy import copy
colour_map = copy(mpl.colormaps.get_cmap('hot_r'))
colour_map.set_under('lightgray')

# plt.rcParams['font.family'] = 'Arial'

root_directory = '/datadrive/data1/axel/'
spatial_directory = f'{root_directory}/spatialdata/Michigan_Xenium/'
figure_directory = f'{spatial_directory}/Figures/'
samples = ['Hip_5K', 'Multiple_480', 'Scalp_5K']

lineages = ['Keratinocyte', 'Fibroblast', 'Lymphocyte', 'Myeloid Cell', 'Lymphatic Endothelial Cell',\
             'Vascular Endothelial Cell', 'Melanocyte', 'Schwann Cell', 'Mural Cell', 'Mast Cell']

global_label = 'integrated_annotation_L1'
ref_label = 'integrated_annotation_L2'

sc.settings.figdir = figure_directory

for sample in samples:

    print(f"\nProcessing: {sample}")

    data_directory = os.path.join(spatial_directory, sample)
        
    for subfolder in Path(data_directory).iterdir():

        if not subfolder.is_dir():
            continue
            
        subfolder_path = subfolder.resolve()
        if subfolder.name.startswith('output'):
            print(f"\nProcessing: {subfolder.name}")

            adata = sc.read_h5ad(subfolder_path / "adata_proseg_filtered.h5ad")

            for lineage in lineages:

                lineage_label = lineage.replace(' ', '-')

                adata_lin = adata[adata.obs['cell_type'] == lineage].copy()
                adata_nonlin = adata[adata.obs['cell_type'] != lineage].copy()

                results_folder = str(subfolder_path / f'cell2loc_{ref_label}_{lineage_label}/')
                run_name = f'{results_folder}/cell2location_map'

                if f'cell2loc_{ref_label}_{lineage_label}' in os.listdir(subfolder_path):

                    print('Getting cell2location results for lineage: ' + lineage_label)
                    sample_file = f"{run_name}/sp.h5ad"
                    adata_cell2loc = sc.read_h5ad(sample_file)

                    celltype_abundances = adata_cell2loc.obsm['q05_cell_abundance_w_sf']

                    celltype_abundances = celltype_abundances.loc[adata_lin.obs_names] # Paranoia

                    adata_lin.obsm[ref_label + '_cell_abundance_w_sf'] = celltype_abundances

                    adata_lin.obs[ref_label + '_est'] = celltype_abundances.idxmax(axis=1).str.replace('q05cell_abundance_w_sf_', '').values   

                    adata_lin.obs[ref_label + '_est'].to_csv(subfolder_path / f'{sample}_{subfolder.name}_{ref_label}_{lineage_label}.csv', index=True)

                    # fig, ax = plt.subplots(figsize=(6, 6))

                    # adata_nonlin.obs['is_' + lineage] = 'False'
                    # adata_nonlin.uns['is_' + lineage + '_colors'] = ['#E4DDD4']

                    # sc.pl.spatial(adata_nonlin, 
                    #             color='is_' + lineage, spot_size=10, alpha=0.5, frameon=False,
                    #             show=False,
                    #             legend_loc = None,
                    #             ax = ax)

                    # sc.pl.spatial(adata_lin, 
                    #             color=ref_label + '_est', spot_size=20, frameon=False,
                    #             show=False,
                    #             ax = ax)
                    
                    # plt.savefig(figure_directory + f'bogle26_xenium_hsca_{ref_label}_{sample}_{subfolder.name}_{lineage_label}.pdf',
                    #             bbox_inches='tight',
                    #             transparent=True)


                    # for _label in sorted(adata_lin.obs[ref_label + '_est'].unique()):
                    #     adata.obs['is_' + _label] = 'False'
                    #     adata.obs.loc[adata_lin[adata_lin.obs[ref_label + '_est'] == _label].obs_names, 'is_' + _label] = 'True'

                    #     adata.obs['is_' + _label] = adata.obs['is_' + _label].astype('category')

                    #     adata.uns['is_' + _label + '_colors'] = ['#E4DDD4', '#ff3319']

                    #     fig, ax = plt.subplots(figsize=(6, 6))

                    #     sc.pl.spatial(adata[(adata.obs['is_' + _label] == 'False')], 
                    #                 color='is_' + _label, spot_size=20, frameon=False,
                    #                 ax=ax, show=False)
                    #     sc.pl.spatial(adata[(adata.obs['is_' + _label] == 'True')], 
                    #                 color='is_' + _label, spot_size=30, frameon=False,
                    #                 ax=ax, show=False)

                    #     plt.savefig(figure_directory + f'bogle26_xenium_hsca_{ref_label}_{sample}_{subfolder.name}_{lineage_label}_{_label.replace(" ", "-").replace(".", "-").replace("/", "-")}.pdf',
                    #                 bbox_inches='tight', transparent=True)
                
                    # adata.write(subfolder_path / f"adata_proseg_filtered_nonadi.h5ad", compression="gzip")
