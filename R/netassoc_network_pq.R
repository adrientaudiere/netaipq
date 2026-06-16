# Species-association network via netassoc for phyloseq objects

################################################################################
#' Species-association network of a phyloseq object using netassoc
#'
#' @description
#' \lifecycle{experimental}
#'
#' Infer a species-association network from a `phyloseq` object with
#' [netassoc::make_netassoc_network()]. The method (Morueta-Holme et al. 2016)
#' compares the observed co-occurrence structure to a null model of random
#' assembly and keeps only the taxa pairs whose association departs
#' significantly from that expectation, returning a partial-correlation
#' network of the taxa.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   obtained using the `phyloseq` package.
#' @param numnulls (int, default 100) Number of null-model replicates used to
#'   build the expectation. Increase for real analyses; the default keeps the
#'   examples and tests fast.
#' @param method (default "partial_correlation") Association measure passed to
#'   [netassoc::make_netassoc_network()].
#' @param plot (logical, default FALSE) Draw the heatmaps produced by
#'   [netassoc::make_netassoc_network()]. FALSE by default to avoid a plotting
#'   side effect; the returned object still contains the networks.
#' @param verbose (logical, default FALSE) Passed to
#'   [netassoc::make_netassoc_network()].
#' @param ... Other parameters passed on to
#'   [netassoc::make_netassoc_network()].
#'
#' @return The list returned by [netassoc::make_netassoc_network()]. Its
#'   `network_all` element is an `igraph` graph of the significant
#'   associations.
#' @export
#' @author Adrien Taudière
#' @references
#'   Morueta-Holme, N. et al. (2016) A network approach for inferring species
#'   associations from co-occurrence data. \doi{10.1111/ecog.01892}.
#' @seealso [cooccurrence_network_pq()], [bipartite_network_pq()]
#' @examples
#' if (requireNamespace("netassoc", quietly = TRUE)) {
#'   ps <- prune_samples(sample_names(data_fungi_mini)[1:5], data_fungi_mini)
#'   ps <- clean_pq(ps, silent = TRUE)
#'   top_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
#'   ps <- prune_taxa(top_taxa, ps)
#'   res <- netassoc_network_pq(ps, numnulls = 20)
#'
#'   netassoc::plot_netassoc_network(res$network_all)
#' }
netassoc_network_pq <- function(
  physeq,
  numnulls = 100,
  method = "partial_correlation",
  plot = FALSE,
  verbose = FALSE,
  ...
) {
  if (!requireNamespace("netassoc", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg netassoc} is required for {.fn netassoc_network_pq}.",
      "i" = "Install it with {.code install.packages(\"netassoc\")}."
    ))
  }
  verify_pq(physeq)
  physeq <- clean_pq(physeq, force_taxa_as_rows = TRUE, silent = TRUE)

  # netassoc expects a species (rows) x sites (columns) matrix
  obs <- as(physeq@otu_table, "matrix")

  netassoc::make_netassoc_network(
    obs,
    method = method,
    numnulls = numnulls,
    plot = plot,
    verbose = verbose,
    ...
  )
}
