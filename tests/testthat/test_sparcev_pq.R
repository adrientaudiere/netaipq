skip_on_cran()
library(MiscMetabar)
data(data_fungi_mini)

test_that("sparcev_pq errors with a clear message when CompoCor is absent", {
  skip_if(requireNamespace("CompoCor", quietly = TRUE))
  expect_error(sparcev_pq(data_fungi_mini, "Time"), "CompoCor")
})

test_that("sparcev_pq returns a per-taxon correlation table", {
  skip_if_not_installed("CompoCor")
  ps <- subset_samples(data_fungi_mini, !is.na(Time))
  res <- sparcev_pq(ps, variable = "Time")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("taxa", "cor", "correlated") %in% colnames(res)))
})
