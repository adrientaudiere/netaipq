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


