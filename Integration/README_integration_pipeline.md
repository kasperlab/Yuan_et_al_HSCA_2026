# Upload Pipeline Overview

This folder contains the upload-facing notebooks and helper files used to document, package, and plot the keratinocyte integration outputs.

HSCA core datasets were processed through the full pipeline, including subtype clustering and cell-type-specific annotation. HSCA extended datasets were processed through first-level clustering and main cell-type annotation only.

## Main iHSCA R Pipeline

- `00_environment.r` sets the output directory and core run parameters, including `hvg_number = 2000` and `pca_dim = 50`.
- `01_data_loading.r` loads and caches the individual study objects, harmonizes metadata, and records study-level cell selection summaries.
- `10_clustering.r` performs first-level atlas clustering.
- `12_annotation.r` assigns first-level integrated annotations.
- `20_clustering_subtype.r` performs second-level clustering for each first-level cell type. For keratinocytes, this creates the generic Hao outputs:
  - `Keratinocyte_2000_hvgs_50_pcs/clustering.h5ad`
  - `Keratinocyte_2000_hvgs_50_pcs/clustering_hvg.h5ad`
- The `21_*` scripts add cell-type-specific second-level annotations and write `IntegratedAnnotation_2_20250531.h5ad` objects for most cell types.

## Separate Keratinocyte Integration

Keratinocytes were additionally processed in the separate Python notebook `UPLOAD_KC_Integration.ipynb`.

That notebook starts from the integrated KC input, corrects sample identifiers, scales cells by sample, recombines them into `KC_full_samples_scaled.h5ad`, tests multiple integration and clustering settings, and saves the working integration objects:

- `KC_full_samples_scaled_Integration_tests.h5ad`
- `KC_full_samples_scaled_Integration_tests_RainbowDuck.h5ad`

The RainbowDuck object is the later version used when bringing Karl's KC integration results back into the manuscript pipeline.

## Merging Karl KC Results Back Into iHSCA

The R script `21_keratinocyte_Karl_v2_20260422.r` connects the separate KC integration back to the existing iHSCA structure.

It reads the cached Hao keratinocyte clustering object from the main pipeline, reads Karl's integrated KC object from `v70_gauss/adata_KC.h5ad`, keeps the shared cell set, transfers Karl's level-2 annotations and integrated UMAP coordinates, harmonizes annotation/color metadata, and writes:

- `Keratinocyte_2000_hvgs_50_pcs/KC_KARL/IntegratedAnnotation_2_20250531.h5ad`

For upload and plotting, these generic pipeline outputs are used under clearer exported names:

- `Keratinocyte_Hao_clustering_hvg.h5ad` corresponds to the Hao `clustering_hvg.h5ad` object.
- `Keratinocyte_Karl_IntegratedAnnotation_2_20250531.h5ad` corresponds to the Karl-merged `IntegratedAnnotation_2_20250531.h5ad` object.

## Plotting and Upload Notebooks

- `UPLOAD_KC_figure_panels.ipynb` reads `Keratinocyte_Hao_clustering_hvg.h5ad` and `Keratinocyte_Karl_IntegratedAnnotation_2_20250531.h5ad` to generate KC figure panels.
- `ED_Figure_1_KC_L1.ipynb` and `ED_Combined_FINAL_plot_generator.ipynb` collect final extended-data plotting outputs.
- `Cluster_annotation_colors_fixed.py` stores fixed annotation color dictionaries used by upload/plotting notebooks.

In short: the main R pipeline generates the standard Hao clustering and annotation structure; the separate KC notebook generates Karl's refined KC integration; `21_keratinocyte_Karl_v2_20260422.r` merges Karl's KC annotations/coordinates back into the iHSCA-style output objects used for upload and plotting.
