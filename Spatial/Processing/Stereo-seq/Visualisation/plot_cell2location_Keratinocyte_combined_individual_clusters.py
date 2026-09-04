
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

lineages = ['Keratinocyte']

adata = sc.read_h5ad(f'{data_directory}clavel25_ctl_occ_nonadi.h5ad')
adata.obsm['spatial'] = adata.obs[['x', 'y']].values
adata.obsm['spatial'][adata.obs['orig.ident'] == 'HS00_CON_Occ'] = adata[adata.obs['orig.ident'] == 'HS00_CON_Occ'].obs[['y', 'x']].values
# adata.obsm['spatial'][adata.obs['orig.ident'] == 'HS00_CON_Occ', 1] *= -1 
# adata.obsm['spatial'][adata.obs['orig.ident'] == 'HS04_CON_Occ', 1] *= -1

for lin in lineages:

    lin_label = lin.replace(" ", "-")
    adata_lin = sc.read_h5ad(f'{data_directory}clavel25_ctl_occ_nonadi_{lin_label}_integrated_annotation_combined.h5ad')

    # adata_lin.obsm['spatial'] = adata_lin.obs[['x', 'y']].values
    # adata_lin.obsm['spatial'][adata_lin.obs['orig.ident'] == 'HS00_CON_Occ'] = adata_lin[adata_lin.obs['orig.ident'] == 'HS00_CON_Occ'].obs[['y', 'x']].values
    # adata_lin.obsm['spatial'][adata_lin.obs['orig.ident'] == 'HS00_CON_Occ', 1] *= -1 
    # adata_lin.obsm['spatial'][adata_lin.obs['orig.ident'] == 'HS04_CON_Occ', 1] *= -1

    lvcombined_results = adata_lin.obsm['integrated_annotation_combined_cell_abundance_w_sf']

    lvcombined_results_norm = lvcombined_results.div(lvcombined_results.sum(1) + 1e-12, axis=0)

    lvcombined_entropy =  -np.nansum(lvcombined_results_norm * np.log(lvcombined_results_norm + 1e-12), axis=1)
    entropy_max = lvcombined_results.shape[1] * np.log(lvcombined_results.shape[1])
    lvcombined_evenness = lvcombined_entropy / entropy_max

    adata_lin.obs['integrated_annotation_combined_est'] = lvcombined_results.idxmax(axis=1).str.replace('q05cell_abundance_w_sf_', '').values   

    adata_lin.obs['integrated_annotation_combined_est'].to_csv(f'{data_directory}clavel25_stereoseq_ctl_occ_nonadi_' + lin_label + '_integrated_annotation_L2_L3_combined.csv')
    
    # for _label in sorted(adata_lin.obs['integrated_annotation_combined_est'].unique()):
    #     adata.obs[ _label] = 'False'
    #     adata.obs.loc[adata_lin[adata_lin.obs['integrated_annotation_combined_est'] == _label].obs_names, _label] = 'True'

    #     adata.obs[ _label] = adata.obs[_label].astype('category')

    #     fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    #     adata.uns[_label + '_colors'] = ['#E4DDD4', '#ff3319']

    #     for _sample, ax in zip(sorted(adata.obs['orig.ident'].unique()), axs.reshape(-1)):
    #         sc.pl.spatial(adata[(adata.obs['orig.ident'] == _sample ) & (adata.obs[_label] == 'False')], 
    #                       color=_label, spot_size=20, frameon=False,
    #                     ax=ax, show=False)
    #         sc.pl.spatial(adata[(adata.obs['orig.ident'] == _sample ) & (adata.obs[_label] == 'True')],
    #                        color=_label, spot_size=30, frameon=False,
    #                     ax=ax, show=False)
            
    #     plt.savefig(figure_directory + 'clavel25_stereoseq_hsca_integrated_annotation_combined_est_' + _label.replace(' ', '-').replace('.', '-').replace('/', '-')+ '.pdf')

   