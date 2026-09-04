
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

import matplotlib as mpl
from copy import copy
colour_map = copy(mpl.colormaps.get_cmap('hot_r'))
colour_map.set_under('lightgray')

plt.rcParams['font.family'] = 'Arial'

data_directory = '/Users/axelalmet/Documents/Data/Omics/ClavelLab_SkinStereoSeq/'
figure_directory = '/Users/axelalmet/Documents/PublicationSubmissions/Papers/InProgress/HSCAMeetsStereoseq_2025/Figures/'
sc.settings.figdir = figure_directory

lineages = ['Keratinocyte', 'Fibroblast', 'Immune', 'Vascular Endothelial Cell', 'Melanocyte', 'Schwann Cell', 'Mural Cell']

adata = sc.read_h5ad(f'{data_directory}clavel25_ctl_occ_nonadi.h5ad')
adata.obsm['spatial'] = adata.obs[['x', 'y']].values
adata.obsm['spatial'][adata.obs['orig.ident'] == 'HS00_CON_Occ'] = adata[adata.obs['orig.ident'] == 'HS00_CON_Occ'].obs[['y', 'x']].values
# adata.obsm['spatial'][adata.obs['orig.ident'] == 'HS00_CON_Occ', 1] *= -1 
# adata.obsm['spatial'][adata.obs['orig.ident'] == 'HS04_CON_Occ', 1] *= -1

for lin in lineages:

    lin_label = lin.replace(" ", "-")
    adata_lin = sc.read_h5ad(f'{data_directory}clavel25_ctl_occ_nonadi_{lin_label}.h5ad')


    # adata_lin.obsm['spatial'] = adata_lin.obs[['x', 'y']].values
    # adata_lin.obsm['spatial'][adata_lin.obs['orig.ident'] == 'HS00_CON_Occ'] = adata_lin[adata_lin.obs['orig.ident'] == 'HS00_CON_Occ'].obs[['y', 'x']].values
    # adata_lin.obsm['spatial'][adata_lin.obs['orig.ident'] == 'HS00_CON_Occ', 1] *= -1 
    # adata_lin.obsm['spatial'][adata_lin.obs['orig.ident'] == 'HS04_CON_Occ', 1] *= -1

    lvl2_results = adata_lin.obsm['integrated_annotation_L2_cell_abundance_w_sf']

    lvl2_results_norm = lvl2_results.div(lvl2_results.sum(1) + 1e-12, axis=0)

    lvl2_entropy =  -np.nansum(lvl2_results_norm * np.log(lvl2_results_norm + 1e-12), axis=1)
    entropy_max = lvl2_results.shape[1] * np.log(lvl2_results.shape[1])
    lvl2_evenness = lvl2_entropy / entropy_max

    adata_lin.obs['integrated_annotation_L2_est'] = lvl2_results.idxmax(axis=1).str.replace('q05cell_abundance_w_sf_', '').values   

    adata_lin.obs['integrated_annotation_L2_est'].to_csv(f'{data_directory}clavel25_stereoseq_ctl_occ_nonadi_' + lin_label + '_integrated_annotation_L2.csv')

    # for _label in sorted(adata_lin.obs['integrated_annotation_L2_est'].unique()):
    #     adata.obs['is_' + _label] = 'False'
    #     adata.obs.loc[adata_lin[adata_lin.obs['integrated_annotation_L2_est'] == _label].obs_names, 'is_' + _label] = 'True'

    #     adata.obs['is_' + _label] = adata.obs['is_' + _label].astype('category')

    #     fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    #     adata.uns['is_' + _label + '_colors'] = ['#E4DDD4', '#ff3319']

    #     for _sample, ax in zip(sorted(adata.obs['orig.ident'].unique()), axs.reshape(-1)):
    #         sc.pl.spatial(adata[(adata.obs['orig.ident'] == _sample ) & (adata.obs['is_' + _label] == 'False')], 
    #                       color='is_' + _label, spot_size=20, frameon=False,
    #                     ax=ax, show=False)
    #         sc.pl.spatial(adata[(adata.obs['orig.ident'] == _sample ) & (adata.obs['is_' + _label] == 'True')],
    #                        color='is_' + _label, spot_size=30, frameon=False,
    #                     ax=ax, show=False)
            
    #     plt.savefig(figure_directory + 'clavel25_stereoseq_hsca_integrated_annotation_L2_est_' + _label.replace(' ', '-').replace('.', '-').replace('/', '-')+ '.pdf')

    #     adata.write(subfolder_path / f"adata_proseg_filtered_nonadi.h5ad", compression="gzip")
