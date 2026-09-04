import os
from pathlib import Path
import scanpy as sc
import scvi
from scipy.sparse import csr_matrix
import torch

torch.set_float32_matmul_precision('high')

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

        adata = sc.read_h5ad(subfolder_path / "adata_proseg_filtered.h5ad")

        # ResolVI requires us to filter out spots with fewer than 5 counts
        adata = adata[adata.X.sum(1) >= 5]
        adata = adata[:, adata.X.sum(0) != 0]

        adata.layers['counts'] = adata.X.copy() 
        adata.obsm['X_spatial'] = adata.obsm['spatial']

        scvi.external.RESOLVI.setup_anndata(adata, layer="counts")

        resolvi = scvi.external.RESOLVI(adata, semisupervised=False)

        resolvi.train(max_epochs=100)

        adata.obsm['X_resolVI'] = resolvi.get_latent_representation()

        sc.pp.neighbors(adata, use_rep="X_resolVI")
        sc.tl.leiden(adata, resolution=0.1, key_added="leiden_resolVI")

        adata.write(subfolder_path / "adata_proseg_filtered.h5ad", compression="gzip")

        resolvi.save(subfolder_path / "resolVI_model_proseg")

