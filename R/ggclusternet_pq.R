# Correlation-network construction and layout via ggClusterNet

################################################################################
#' Microbial correlation network of a phyloseq object using ggClusterNet
#'
#' @description
#' \lifecycle{experimental}
#'
#' Build a per-group microbial correlation network from a `phyloseq` object
#' with `ggClusterNet::network.2()` (Wen et al. 2022). For each level of
#' `group`, taxa are correlated, thresholded on correlation strength and
#' significance, laid out with one of the `ggClusterNet` network layouts and
#' returned together with a `ggplot2` rendering.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   obtained using the `phyloseq` package.
#' @param group (required) Name of the `sam_data` column defining the groups
#'   for which a separate network is built.
#' @param n (int, default 500) Number of most-abundant taxa kept before the
#'   analysis (the `N` argument of `ggClusterNet::network.2()`).
#' @param r_threshold (float, default 0.6) Minimum absolute correlation for an
#'   edge to be kept.
#' @param p_threshold (float, default 0.05) Maximum p-value for an edge to be
#'   kept.
#' @param maxnode (int, default 5) Maximum node size in the layout.
#' @param layout_net (default "model_maptree2") The `ggClusterNet` layout
#'   algorithm.
#' @param big (logical, default TRUE) Use the large-network code path of
#'   `ggClusterNet::network.2()`.
#' @param select_layout (logical, default TRUE) Let `ggClusterNet` adjust the
#'   layout.
#' @param label (logical, default FALSE) Draw taxa labels on the network.
#' @param zipi (logical, default FALSE) Compute the within/among-module
#'   connectivity (Zi-Pi) classification.
#' @param ... Other parameters passed on to `ggClusterNet::network.2()`.
#'
#' @return The list returned by `ggClusterNet::network.2()`, whose first
#'   element is a `ggplot2` network plot and whose remaining elements hold the
#'   network and correlation data.
#' @export
#' @author Adrien Taudière
#' @references
#'   Wen, T. et al. (2022) ggClusterNet: an R package for microbiome network
#'   analysis and modularity-based multiple network layouts.
#'   \doi{10.1002/imt2.32}.
#' @examples
#' \dontrun{
#' # `ggClusterNet` is a GitHub package and needs WGCNA (Bioconductor):
#' # remotes::install_github("taowenmicro/ggClusterNet")
#' # pak::pkg_install("WGCNA")
#' res <- ggclusternet_pq(data_fungi_mini, group = "Height", n = 200)
#' res[[1]]
#' }
ggclusternet_pq <- function(
  physeq,
  group,
  n = 500,
  r_threshold = 0.6,
  p_threshold = 0.05,
  maxnode = 5,
  layout_net = "model_maptree2",
  big = TRUE,
  select_layout = TRUE,
  label = FALSE,
  zipi = FALSE,
  ...
) {
  if (!requireNamespace("ggClusterNet", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg ggClusterNet} is required for {.fn ggclusternet_pq}.",
      "i" = "Install it with {.code remotes::install_github(\"taowenmicro/ggClusterNet\")}."
    ))
  }
  verify_pq(physeq)

  if (!group %in% colnames(physeq@sam_data)) {
    cli::cli_abort(
      "{.val {group}} is not a column of the {.field sam_data} slot."
    )
  }
  physeq <- clean_pq(physeq, silent = TRUE)

  ggClusterNet::network.2(
    ps = physeq,
    N = n,
    group = group,
    big = big,
    maxnode = maxnode,
    select_layout = select_layout,
    layout_net = layout_net,
    r.threshold = r_threshold,
    p.threshold = p_threshold,
    label = label,
    zipi = zipi,
    ...
  )
}
