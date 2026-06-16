skip_on_cran()
library(MiscMetabar)
data(data_fungi_mini)

test_that("dag_test_pq returns localTests results", {
  skip_if_not_installed("dagitty")
  dag <- "dag {
    Time -> nb_otu
    nb_seq -> nb_otu
    Time -> Height
    Height -> nb_otu
  }"
  res <- dag_test_pq(data_fungi_mini, dag, as_factor = c("Time", "Height"))
  expect_s3_class(res, "data.frame")
  expect_true("p.value" %in% colnames(res))
  expect_gt(nrow(res), 0)
})

test_that("dag_test_pq errors on a DAG variable absent from sample data", {
  skip_if_not_installed("dagitty")
  expect_error(
    dag_test_pq(data_fungi_mini, "dag { Time -> not_a_column }"),
    "not found in the sample data"
  )
})
