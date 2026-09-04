import os
from pathlib import Path
import scanpy as sc
import scvi
import numpy as np
from scipy.sparse import csr_matrix
import torch
import cell2location

from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

import matplotlib.pyplot as plt

torch.set_float32_matmul_precision('high')

root_directory = '/datadrive/data1/axel/'
spatial_directory = f'{root_directory}/spatialdata/Michigan_Xenium/'
sc_directory = f'{root_directory}/scdata/'

samples = ['Hip_5K', 'Multiple_480', 'Scalp_5K']
adata_sc = sc.read_h5ad(f'{sc_directory}/hsca_merged.h5ad')

ref_label = 'integrated_annotation_L3'
global_label = 'integrated_annotation_L2'

adata_sc = adata_sc[adata_sc.obs[ref_label] != 'nan'].copy()
global_clusters = ['KC - 19: uHF', \
                   'KC - 20: HF_1' 'KC - 21: HF_2', 'KC - 22: HF_3',\
                   'KC - 23: HF_Anagen', \
                    'KC - 24: SG_1', 'KC - 25: SG_2', \
                    'KC - 26: Gland_Channel_1', 'KC - 27: Gland_Channel_2', 'KC - 28: Gland_Channel_3', \
                    'KC - 29: Gland_Sweat_1', 'KC - 30: Gland_Sweat_2']

ref_folder = f'{sc_directory}/hsca_cell2location_{ref_label}'

for sample in samples:

    print(f"\nProcessing: {sample}")

    data_directory = os.path.join(spatial_directory, sample)
        
    for subfolder in Path(data_directory).iterdir():

        if not subfolder.is_dir():
            continue
            
        subfolder_path = subfolder.resolve()
        if subfolder.name.startswith('output'):
            print(f"\nProcessing: {subfolder.name}")

            adata =  sc.read_h5ad(subfolder_path / "adata_proseg_filtered.h5ad")

            if 'cell_type' in adata.obs.columns:

                adata_ker = adata[adata.obs['cell_type'] == 'Keratinocyte'].copy()

                results_folder = str(subfolder_path / f'cell2loc_{ref_label}/')

                if not os.path.exists(results_folder):

                    os.mkdir(results_folder)

                # create paths and names to results folders for reference regression and cell2location models
                ref_run_name = f'{ref_folder}/reference_signatures'
                run_name = f'{results_folder}/cell2location_map'

                if not os.path.exists(run_name):
                    os.mkdir(run_name)

                if os.path.exists(f"{run_name}/sp.h5ad"):
                    print(f"Results already exist for {ref_label} in {sample}, skipping...")
                    continue

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

                intersect = np.intersect1d(adata_ker.var_names, inf_aver.index)
                adata_ker = adata_ker[:, intersect].copy()
                inf_aver = inf_aver.loc[intersect, :].copy()

                print(adata_ker)
                print(adata_ker.var_names)

                # prepare anndata for cell2location model
                cell2location.models.Cell2location.setup_anndata(adata=adata_ker, layer='counts')

                mod = cell2location.models.Cell2location(
                    adata_ker, cell_state_df=inf_aver,
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
                
                adata_pos = mod.export_posterior(
                    adata_ker, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'accelerator': 'gpu'}
                )

                # Save anndata object with results
                adata_file = f"{run_name}/sp.h5ad"
                adata_pos.write(adata_file)
                adata_file
                
                # Save model
                mod.save(f"{run_name}", overwrite=True)

                fig, ax = plt.subplots(figsize=(4, 3))
                mod.plot_history(iter_start = 1000, ax=ax)
                plt.savefig(f"{run_name}/training_history.png",bbox_inches='tight')

                # mod = cell2location.models.Cell2location.load(f"{run_name}", adata)



