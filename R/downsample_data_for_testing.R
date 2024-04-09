nmf_test <- readRDS('~/Library/CloudStorage/Box-Box/neely.jdm/NMF/results_summary/nmf_ksweep_results_compiled.rds')

#Goal: downsample NMF results from JDM crosslong experiment to 1k cells each, same cells across all ranks
#do this later lol it might be faster to just save nmf output on 1k cell test seu
l

rownames(y) <- NULL # because as.matrix is faster without rownames
isRowSorted(as.matrix(y))
any(sapply(nmf_test$Bcells$consensus.results, function(X){
  cells = rownames(X[['H']])
}))
