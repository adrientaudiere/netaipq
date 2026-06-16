# Compositional correlation with an environmental variable (SparCEV)

################################################################################
#' Compositional taxa-environment correlation with SparCEV
#'
#' @description
#' \lifecycle{experimental}
#'
#' Compute the sparse correlation of each taxon with a continuous
#' environmental (sample) variable using `CompoCor::SparCEV()` (Sparse
#' Correlation with an Environmental Variable). SparCEV accounts for the
#' compositional nature of metabarcoding count data and flags the taxa whose
#' abundance is significantly correlated with the variable.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   obtained using the `phyloseq` package.
#' @param variable (required) Name of the continuous variable in the
#'   `sam_data` slot to correlate the taxa with.
#' @param add_tax_table (logical, default TRUE) Join the `tax_table` of
#'   `physeq` to the result so each taxon carries its taxonomy.
#' @param ... Other parameters passed on to `CompoCor::SparCEV()`.
#'
#' @return A data frame with one row per taxon: `taxa` (taxon name), `cor`
#'   (the SparCEV correlation with `variable`) and `correlated` (logical, TRUE
#'   when the correlation is deemed significant). When `add_tax_table = TRUE`
#'   the columns of the `tax_table` are appended. The rendering of this result
#'   (e.g. a labelled scatter plot) belongs in `ggplotpq`.
#' @export
#' @author Adrien Taudière
#' @seealso [cooccurrence_network_pq()]
#' @examples
#' \dontrun{
#' # `CompoCor` is a GitHub package:
#' # remotes::install_github("IbTJensen/CompoCor")
#' ps <- subset_samples(data_fungi_mini, !is.na(Time))
#' res <- sparcev_pq(ps, variable = "Time")
#'  res |>
#'    filter(correlated) |>
#'    ggplot(aes(x = cor, y = Family, shape = Guild, color = Order)) +
#'    geom_point() +
#'    geom_text(aes(label = taxa), nudge_y = 0.3, size = 2.5)
#' }
sparcev_pq <- function(
  physeq,
  variable,
  add_tax_table = TRUE,
  ...
) {
  if (!requireNamespace("CompoCor", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg CompoCor} is required for {.fn sparcev_pq}.",
      "i" = "Install it with {.code remotes::install_github(\"IbTJensen/CompoCor\")}."
    ))
  }
  verify_pq(physeq)
  physeq <- clean_pq(physeq, silent = TRUE)

  if (!variable %in% colnames(physeq@sam_data)) {
    cli::cli_abort(
      "{.val {variable}} is not a column of the {.field sam_data} slot."
    )
  }
  env <- as.vector(physeq@sam_data[[variable]])

  # SparCEV expects a samples (rows) x taxa (columns) matrix
  otu <- as(physeq@otu_table, "matrix")
  if (physeq@otu_table@taxa_are_rows) {
    otu <- t(otu)
  }

  res <- CompoCor::SparCEV(otu, env, ...)

  out <- data.frame(
    taxa = names(res$cor),
    cor = as.numeric(res$cor),
    correlated = res$Correlated,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  if (add_tax_table) {
    tax <- as.data.frame(unclass(physeq@tax_table))
    tax$taxa <- rownames(tax)
    out <- merge(out, tax, by = "taxa", all.x = TRUE)
  }
  out
}
