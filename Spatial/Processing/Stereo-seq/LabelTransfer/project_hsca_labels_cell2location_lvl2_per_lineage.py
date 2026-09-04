import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scvi
import cell2location

import torch

import pyro


from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

torch.set_float32_matmul_precision('high')

data_directory = '/datadrive/data1/axel/'
adata = sc.read_h5ad(f'{data_directory}spatialdata/clavel25_ctl_occ_nonadi.h5ad')

ref_label = 'integrated_annotation_L2'
global_label = 'integrated_annotation_L1'
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
    ref_run_name = f'{results_folder}/reference_signatures'
    run_name = f'{results_folder}/cell2location_map'

    adata_file = f"{ref_run_name}/sc.h5ad"
    adata_ref = sc.read_h5ad(adata_file)
    mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

    # export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']
    inf_aver.iloc[0:5, 0:5]

    intersect = np.intersect1d(adata_lin.var_names, inf_aver.index)
    adata_lin = adata_lin[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    print(adata_lin)
    print(adata_lin.var_names)

    for sample in adata_lin.obs['orig.ident'].unique():
        adata_sample = adata_lin[adata_lin.obs['orig.ident'] == sample].copy()

        # prepare anndata for cell2location model
        cell2location.models.Cell2location.setup_anndata(adata=adata_sample, layer='counts')

        mod = cell2location.models.Cell2location(
            adata_sample, cell_state_df=inf_aver,
            # the expected average cell abundance: tissue-dependent
            # hyper-prior which can be estimated from paired histology:
            N_cells_per_location=1,
            # hyperparameter controlling normalisation of
            # within-experiment variation in RNA detection:
            detection_alpha=20
        )
        mod.view_anndata_setup()

        mod.train(max_epochs=20000,
                # train using full data (batch_size=None)
                batch_size=None,
                # use all data points in training because
                # we need to estimate cell abundance at all locations
                train_size=1,
                accelerator='gpu',
                )


        fig, ax = plt.subplots(figsize=(4, 3))
        mod.plot_history(iter_start = 1000, ax=ax)
        plt.savefig(f"{run_name}_{lineage}_{sample}_training_history.png",bbox_inches='tight')

        adata_sample = mod.export_posterior(
            adata_sample, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'accelerator': 'gpu'}
        )

        # Save model
        mod.save(f"{run_name}_{lineage}_{sample}_nonadi", overwrite=True)

        # mod = cell2location.models.Cell2location.load(f"{run_name}", adata)

        # Save anndata object with results
        adata_file = f"{run_name}_{lineage}_{sample}_nonadi/sp.h5ad"
        adata_sample.write(adata_file)
        adata_file