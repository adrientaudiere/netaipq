skip_on_cran()
library(MiscMetabar)
data(data_fungi_mini)

test_that("cooccurrence_network_pq returns a cooccur object", {
  skip_if_not_installed("cooccur")
  res <- cooccurrence_network_pq(data_fungi_mini, n_taxa = 20)
  expect_s3_class(res, "cooccur")
  expect_s3_class(res$results, "data.frame")
  expect_lte(res$species, 20)
})

test_that("cooccurrence_network_pq accepts abundance data", {
  skip_if_not_installed("cooccur")
  res <- cooccurrence_network_pq(data_fungi_mini, binary = FALSE, n_taxa = 15)
  expect_s3_class(res, "cooccur")
})
