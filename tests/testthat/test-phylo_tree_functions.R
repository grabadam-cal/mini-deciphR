test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
nmf_test <- readRDS('~/Library/CloudStorage/Box-Box/neely.jdm/NMF/results_summary/nmf_ksweep_results_compiled.rds')
test_phylo = readRDS('~/Library/CloudStorage/Box-Box/neely.jdm/NMF/DECIPHER_outputs/crosslong_phylo_trees.rds')
