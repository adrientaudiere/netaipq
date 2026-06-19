skip_on_cran()
library(MiscMetabar)
data(data_fungi_mini)

mini_pq <- function() {
  ps <- prune_samples(sample_names(data_fungi_mini)[1:5], data_fungi_mini)
  ps <- clean_pq(ps, silent = TRUE)
  top_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
  prune_taxa(top_taxa, ps)
}

test_that("netassoc_network_pq returns a network list with an igraph", {
  skip_if_not_installed("netassoc")
  res <- suppressWarnings(netassoc_network_pq(mini_pq(), numnulls = 20))
  expect_type(res, "list")
  expect_true("network_all" %in% names(res))
  expect_s3_class(res$network_all, "igraph")
})
