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
adata_ref = sc.read_h5ad(f'{data_directory}scdata/hsca_merged.h5ad')

ref_label = 'integrated_annotation_L1'

results_folder = f'{data_directory}/scdata/hsca_cell2location_{ref_label}/'

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'


selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()

# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='sample_ID',
                        # cell type, covariate used for constructing signatures
                        labels_key=ref_label,
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['assay', 'donor_id'],
                        layer = 'counts'
                       )

mod = RegressionModel(adata_ref)

mod.view_anndata_setup()

num_epochs = 300
mod.train(max_epochs=num_epochs, accelerator='gpu')


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'accelerator': 'gpu'}
)


# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file


fig, ax = plt.subplots(figsize=(4, 3))
mod.plot_history(iter_start = 20, ax=ax)
plt.savefig(f"{ref_run_name}_training_history_{num_epochs}.png",bbox_inches='tight')
