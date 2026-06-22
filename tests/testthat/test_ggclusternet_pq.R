skip_on_cran()
library(MiscMetabar)
data(data_fungi_mini)

test_that("ggclusternet_pq returns the network.2 list", {
  # ggClusterNet needs its full (GitHub + Bioconductor) dependency chain to run.
  skip_if_not_installed("ggClusterNet")
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("tidyfst")
  data_fungi_mini_woNA4height <- subset_samples(
    data_fungi_mini,
    !is.na(data_fungi_mini@sam_data$Height)
  )
  suppressWarnings(
    res <- ggclusternet_pq(
      data_fungi_mini_woNA4height,
      group = "Height",
      n = 100
    )
  )
  expect_type(res, "list")
})
