

tryCatch({
  py_run_string("print(annotation_1)")
}, error = function(x){
  py$annotation_1 = readRDS(paste0(py$output_dir_cache, "/annotation_1.rds"))
})
tryCatch({
  py_run_string("print(adata_imbalance)")
}, error = function(x){
  py_run_string(paste0("adata_imbalance = ad.read_h5ad('", py$output_dir_cache, "/adata_imbalance_annotated_1.h5ad')"))
  py_run_string(paste0("adata_imbalance.uns['log1p']['base'] = None"))
  print(paste0("Loaded adata_imbalance from cache."))
})


