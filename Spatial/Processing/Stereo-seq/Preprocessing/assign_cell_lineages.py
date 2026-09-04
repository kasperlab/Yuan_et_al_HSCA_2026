import scanpy as sc

data_directory = '/Users/axelalmet/Documents/Data/Omics/ClavelLab_SkinStereoSeq/'

adata = sc.read_h5ad(data_directory + 'clavel25_ctl_occ_nonadi.h5ad')

lineages_and_celltypes = {'Fibroblast': ['Fibroblast I', 'Fibroblast II', 'Dermal Cup', 'Dermal Papilla', 'Dermal Sheath I', 'Dermal Sheath II'],
'Keratinocyte': ['Basal', 'Basal Progenitor', 'Suprabasal', 'Granular', 'Spinous', 'Corneum',
                    'Bulge', 'Cortex and Precortex', 'Hair Germ', 'Hair Keratins', 'Infundibulum', 'Isthmus', 'IRS', 'Lower Bulge', 
                     'Matrix I', 'Matrix II', 'Matrix III', 'Matrix IV', 'Matrix V', 'Medulla', 'ORS', 'Upper Follicle Barrier', 
                     'Basal Myoepithelial',
                    'Eccrine Ducts I', 'Eccrine Ducts II', 'Eccrine Secretory I', 'Eccrine Secretory II',
                    'Sebaceous I', 'Sebaceous II', 'Sebaceous III', 'Sebaceous Diff. I', 'Sebaceous Diff. II'], 
'Immune': ['B cell', 'T cell', 'Mast', 'Myeloid'],
'Schwann Cell': ['Schwann Cells'],
'Melanocyte': ['Me_Bulge', 'Me_Eccrine', 'Me_Epidermis', 'Me_Matrix', 'Me_ReteRidges'],
'Mural Cell': ['Smooth Muscle', 'Smooth Muscle Arrector Pili'],
'Vascular Endothelial Cell': ['Endothelial'],
'Schwann Cell': ['Schwann cell']}

for lineage, celltypes in lineages_and_celltypes.items():
    adata.obs.loc[adata.obs['labels7'].isin(celltypes), 'cell_lineage'] = lineage

print(adata.obs['cell_lineage'].value_counts())

adata.write(data_directory + 'clavel25_ctl_occ_nonadi.h5ad', compression='gzip')