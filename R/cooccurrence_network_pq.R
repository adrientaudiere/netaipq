# Probabilistic co-occurrence analysis for phyloseq objects

################################################################################
#' Probabilistic species co-occurrence from a phyloseq object
#'
#' @description
#' \lifecycle{experimental}
#'
#' Compute the probabilistic species co-occurrence model of Veech (2013) on
#' the taxa of a `phyloseq` object using [cooccur::cooccur()]. Each pair of
#' taxa is classified as positively associated, negatively associated or
#' randomly associated given the observed co-occurrence across samples.
#'
#' By default the OTU table is converted to presence/absence with
#' [MiscMetabar::as_binary_otu_table()] before the analysis: `cooccur` is a
#' binary (incidence) model and gives spurious results on raw abundance data.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   obtained using the `phyloseq` package.
#' @param binary (logical, default TRUE) Convert the OTU table to
#'   presence/absence before the analysis. Keep TRUE unless you know your data
#'   suit the (rarely appropriate) abundance behaviour of `cooccur`.
#' @param n_taxa (int, default NULL) If set, keep only the `n_taxa` most
#'   prevalent taxa before the analysis. `cooccur` scales quadratically with
#'   the number of taxa, so capping is recommended on rich datasets.
#' @param spp_names (logical, default TRUE) Passed to [cooccur::cooccur()]:
#'   keep taxa names in the output.
#' @param prob (default "comb") The probability method passed to
#'   [cooccur::cooccur()], either "comb" (combinatorics) or "hyper"
#'   (hypergeometric).
#' @param thresh (logical, default TRUE) Passed to [cooccur::cooccur()]:
#'   remove taxa pairs expected to co-occur less than one time.
#' @param ... Other parameters passed on to [cooccur::cooccur()].
#'
#' @return An object of class `cooccur` (see [cooccur::cooccur()]). Use
#'   [cooccur::effect.sizes()], [cooccur::prob.table()] or `summary()` to
#'   explore it, and `plot()` to visualise the association matrix.
#' @export
#' @author Adrien Taudière
#' @references
#'   Veech, J. A. (2013) A probabilistic model for analysing species
#'   co-occurrence. \doi{10.1111/j.1466-8238.2012.00789.x}.
#' @examples
#' if (requireNamespace("cooccur", quietly = TRUE)) {
#'   res <- cooccurrence_network_pq(data_fungi_mini, n_taxa = 30)
#'   summary(res)
#'   plot(res)
#' }
cooccurrence_network_pq <- function(
  physeq,
  binary = TRUE,
  n_taxa = NULL,
  spp_names = TRUE,
  prob = "comb",
  thresh = TRUE,
  ...
) {
  if (!requireNamespace("cooccur", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg cooccur} is required for {.fn cooccurrence_network_pq}.",
      "i" = "Install it with {.code install.packages(\"cooccur\")}."
    ))
  }
  verify_pq(physeq)
  physeq <- clean_pq(physeq, silent = TRUE)
  if (binary) {
    physeq <- MiscMetabar::as_binary_otu_table(physeq)
  }

  # cooccur expects a taxa (rows) x samples (columns) incidence matrix
  otu <- as(physeq@otu_table, "matrix")
  if (!physeq@otu_table@taxa_are_rows) {
    otu <- t(otu)
  }

  if (!is.null(n_taxa)) {
    keep <- order(rowSums(otu > 0), decreasing = TRUE)[
      seq_len(min(n_taxa, nrow(otu)))
    ]
    otu <- otu[keep, , drop = FALSE]
  }
  storage.mode(otu) <- "integer"

  cooccur::cooccur(
    mat = as.data.frame(otu),
    spp_names = spp_names,
    prob = prob,
    thresh = thresh,
    ...
  )
}
