
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

root_directory = '/Users/axelalmet/Documents/Data/Omics/'
spatial_directory = f'{root_directory}Bogle2026/'
samples = ['Hip_5K', 'Multiple_480', 'Scalp_5K']

ref_label = 'integrated_annotation_L1'


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
            
            if 'cell_type' in adata.obs:
                print('cell_type annotation present in adata.obs')
        