# Test a causal DAG against the sample data of a phyloseq object

################################################################################
#' Test a causal DAG against phyloseq sample data
#'
#' @description
#' \lifecycle{experimental}
#'
#' Confront a hypothesised causal directed acyclic graph (DAG) with the
#' sample data of a `phyloseq` object. The implied conditional independencies
#' of the DAG are tested against the data with [dagitty::localTests()]: each
#' independence that the DAG claims is checked, and a large estimate signals a
#' DAG that is inconsistent with the data.
#'
#' Diversity descriptors of the samples (`nb_seq`, the sequencing depth, and
#' `nb_otu`, the observed richness) are added to the sample data with
#' [MiscMetabar::add_info_to_sam_data()] so they can be used as DAG nodes
#' alongside the columns of `sam_data`.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   obtained using the `phyloseq` package.
#' @param dag (required) A causal DAG, either a `dagitty` object or a
#'   `dagitty` model string (e.g. `"dag { Time -> nb_otu nb_seq -> nb_otu }"`).
#'   Every node of the DAG must be a column of the (enriched) sample data.
#' @param as_factor (character vector, default NULL) Names of DAG variables to
#'   coerce to factor before the test (useful for categorical metadata stored
#'   as numbers or characters).
#' @param type (default "cis.chisq") The independence-test type passed to
#'   [dagitty::localTests()]. See its help for the available options.
#' @param na_remove (logical, default TRUE) Remove samples with a missing value
#'   in any DAG variable before the test.
#' @param ... Other parameters passed on to [dagitty::localTests()].
#'
#' @return The data frame returned by [dagitty::localTests()], one row per
#'   implied conditional independence, with the estimate and its confidence
#'   interval. Pass it to [dagitty::plotLocalTestResults()] to visualise it.
#' @export
#' @author Adrien Taudière
#' @references
#'   Ankan, A., Wortel, I. M. N. & Textor, J. (2021) Testing graphical causal
#'   models using the R package "dagitty". \doi{10.1002/cpz1.45}.
#' @examples
#' \donttest{
#' if (requireNamespace("dagitty", quietly = TRUE)) {
#'   dag <- "dag {
#'     Time -> nb_otu
#'     nb_seq -> nb_otu
#'     Time -> Height
#'     Height -> nb_otu
#'   }"
#'   dag_test_pq(data_fungi_mini, dag, as_factor = c("Time", "Height"))
#' }
#' }
#' \dontrun{
#' if (requireNamespace("dagitty", quietly = TRUE)) {
#'   # Add Hill number q = 2 (inverse Simpson) to sample data
#'   hill2 <- phyloseq::estimate_richness(data_fungi_mini, measures = "InvSimpson")
#'   phyloseq::sample_data(data_fungi_mini)$hill_2 <- hill2$InvSimpson
#'   dag2 <- "dag {
#'     Time -> hill_2
#'     nb_seq -> hill_2
#'     Diameter -> hill_2
#'     Time -> Diameter
#'   }"
#'   dag_test_pq(data_fungi_mini, dag2, as_factor = "Time")
#' }
#' }
dag_test_pq <- function(
  physeq,
  dag,
  as_factor = NULL,
  type = "cis.chisq",
  na_remove = TRUE,
  ...
) {
  if (!requireNamespace("dagitty", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg dagitty} is required for {.fn dag_test_pq}.",
      "i" = "Install it with {.code install.packages(\"dagitty\")}."
    ))
  }
  verify_pq(physeq)
  if (is.character(dag)) {
    dag <- dagitty::dagitty(dag)
  }

  physeq <- MiscMetabar::add_info_to_sam_data(physeq)
  df <- as(physeq@sam_data, "data.frame")

  dag_edges <- dagitty::edges(dag)
  dag_vars <- unique(as.character(c(dag_edges$v, dag_edges$w)))

  missing_vars <- setdiff(dag_vars, colnames(df))
  if (length(missing_vars) > 0) {
    cli::cli_abort(c(
      "DAG variable{?s} {.val {missing_vars}} not found in the sample data.",
      "i" = "Available variables: {.val {colnames(df)}}."
    ))
  }

  df <- df[, dag_vars, drop = FALSE]
  for (v in as_factor) {
    df[[v]] <- as.factor(df[[v]])
  }
  if (na_remove) {
    df <- df[stats::complete.cases(df), , drop = FALSE]
  }

  dagitty::localTests(x = dag, data = df, type = type, ...)
}
