import os
from pathlib import Path
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix

import matplotlib.pyplot as plt


root_directory = '/datadrive/data1/axel/'
spatial_directory = f'{root_directory}/spatialdata/Michigan_Xenium/'

samples = ['Hip_5K', 'Multiple_480', 'Scalp_5K']

ref_label = 'integrated_annotation_L1'

for sample in samples:

    print(f"\nProcessing: {sample}")

    data_directory = os.path.join(spatial_directory, sample)
        
    for subfolder in Path(data_directory).iterdir():

        if not subfolder.is_dir():
            continue
            
        subfolder_path = subfolder.resolve()
        print(f"\nProcessing: {subfolder.name}")

        if (subfolder_path / "adata_proseg_filtered_nonadi.h5ad").exists():

            adata = sc.read_h5ad(subfolder_path / "adata_proseg_filtered_nonadi.h5ad")
        
        else:
            adata =  sc.read_h5ad(subfolder_path / "adata_proseg_filtered.h5ad")

        results_folder = str(subfolder_path / f'cell2loc_{ref_label}/')

        # create paths and names to results folders for reference regression and cell2location models
        run_name = f'{results_folder}/cell2location_map'

        file_name = f"{run_name}/sp.h5ad"
        adata_cell2loc = sc.read_h5ad(file_name)

        celltype_abundances = adata_cell2loc.obsm['q05_cell_abundance_w_sf']
        celltype_abundances = celltype_abundances.loc[adata.obs_names] # Paranoia

        adata.obsm[ref_label + '_cell_abundance_w_sf'] = celltype_abundances

        adata.write(subfolder_path / f"adata_proseg_filtered_nonadi.h5ad", compression="gzip")
