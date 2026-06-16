skip_on_cran()
library(MiscMetabar)
data(data_fungi_mini)

test_that("ggclusternet_pq errors clearly when ggClusterNet is absent", {
  skip_if(requireNamespace("ggClusterNet", quietly = TRUE))
  expect_error(
    ggclusternet_pq(data_fungi_mini, group = "Height"),
    "ggClusterNet"
  )
})

test_that("ggclusternet_pq returns the network.2 list", {
  # ggClusterNet needs its full (GitHub + Bioconductor) dependency chain to run.
  skip_if_not_installed("ggClusterNet")
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("tidyfst")
  res <- ggclusternet_pq(data_fungi_mini, group = "Height", n = 100)
  expect_type(res, "list")
})
