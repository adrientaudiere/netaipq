skip_on_cran()
library(MiscMetabar)
data(data_fungi_mini)

test_that("ggclusternet_pq drops NA-in-group samples and returns the list", {
  # ggClusterNet needs its full (GitHub + Bioconductor) dependency chain to run.
  skip_if_not_installed("ggClusterNet")
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("tidyfst")
  # Regression: this direct call used to fail with "Component sample names do
  # not match" because `Height` contains NAs. It must now succeed and inform
  # the user that samples were dropped.
  expect_true(any(is.na(data_fungi_mini@sam_data$Height)))
  suppressWarnings(
    expect_message(
      res <- ggclusternet_pq(data_fungi_mini, group = "Height", n = 100),
      "Dropping"
    )
  )
  expect_type(res, "list")
})
