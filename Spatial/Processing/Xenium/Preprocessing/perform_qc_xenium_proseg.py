import os
from pathlib import Path
import spatialdata  as sd
import squidpy as sq
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import spotsweeper.local_outliers as lo

root_directory = '/datadrive/data1/axel/spatialdata/Michigan_Xenium/'

samples = ['Hip_5K', 'Multiple_480', 'Scalp_5K']

for sample in samples:

    print(f"\nProcessing: {sample}")

    data_directory = os.path.join(root_directory, sample)
        
    for subfolder in Path(data_directory).iterdir():

        if not subfolder.is_dir():
                continue
            
        subfolder_path = subfolder.resolve()
        print(f"\nProcessing: {subfolder.name}")

        sdata_proseg = None

        if (subfolder_path / "proseg-output.zarr").exists():
             sdata = sd.read_zarr(subfolder_path / "proseg-output.zarr")
             adata = sdata.tables['table']

        if sdata is None:
            print(f"Missing data for {subfolder.name}, skipping.")
            continue
        
        adata = adata[:, adata.X.sum(0) != 0]

        adata.var['mito'] = adata.var_names.str.startswith('MT-')

        n_vars = adata.n_vars
        percent_top = [n for n in (50, 100, 200, 500) if n < n_vars] # Accounts for the 480-gene panel data
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mito'], percent_top=percent_top, inplace=True)

        # Identify local spot outliers with excessively low total counts
        lo.local_outliers(
            adata,
            metric="log1p_total_counts",
            direction="lower",
            n_neighbors=36,
            sample_key="region",
            log=False,
            cutoff=3.0,
            coord_key="spatial",
        )

        # Identify local spot outliers with excessively low number of unique genes detected
        lo.local_outliers(
            adata,
            metric="log1p_n_genes_by_counts",
            direction="lower",
            n_neighbors=36,
            sample_key="region",
            log=False,
            cutoff=3.0,
            coord_key="spatial",
        )

        # Identify local spot outliers with excessively high percentage of mitochondrial gene counts
        lo.local_outliers(
            adata,
            metric="pct_counts_mito",
            direction="higher",
            n_neighbors=36,
            sample_key="region",
            log=False,
            cutoff=3.0,
            coord_key="spatial",
        )

        adata.obs['outliers'] = adata.obs['log1p_total_counts_outliers'] | adata.obs['log1p_n_genes_by_counts_outliers'] | adata.obs['pct_counts_mito_outliers']
        adata.write(subfolder_path / "adata_proseg.h5ad", compression="gzip")

        print(adata)

        adata_filtered = adata[~adata.obs['outliers']]
        adata_filtered = adata_filtered[:, adata_filtered.X.sum(0) != 0]

        print(adata_filtered)

        adata_filtered.write(subfolder_path / "adata_proseg_filtered.h5ad", compression="gzip")