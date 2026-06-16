# Bipartite-network visualisation of a phyloseq object

# Resolve a label argument (a single TRUE/FALSE or a logical vector) into a
# logical vector aligned with `ids`.
resolve_node_label <- function(x, ids, arg_name) {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) {
    return(rep(x, length(ids)))
  }
  if (is.logical(x) && length(x) == length(ids)) {
    x[is.na(x)] <- FALSE
    return(x)
  }
  cli::cli_abort(c(
    "{.arg {arg_name}} must be a single {.code TRUE}/{.code FALSE} or a logical vector of length {length(ids)}.",
    "x" = "You supplied a {.cls {class(x)}} of length {length(x)}."
  ))
}

# Dissimilarity indices accepted by `projection` (from vegan::vegdist).
vegdist_methods <- c(
  "manhattan",
  "euclidean",
  "canberra",
  "clark",
  "bray",
  "kulczynski",
  "jaccard",
  "gower",
  "altGower",
  "morisita",
  "horn",
  "mountford",
  "raup",
  "binomial",
  "chao",
  "cao",
  "mahalanobis",
  "chisq",
  "chord",
  "hellinger",
  "aitchison",
  "robust.aitchison"
)

# Extract a samples x 2 coordinate matrix (row names = sample names) from a
# user-supplied projection. Accepts a numeric matrix with sample row names, or
# a data frame / tibble where the sample names are either the row names or one
# column, and the coordinates are a known axis-name pair (e.g. x_umap/y_umap,
# as returned by MiscMetabar::umap_pq()) or the first two numeric columns.
as_sample_coords <- function(projection, samples) {
  if (
    is.matrix(projection) &&
      is.numeric(projection) &&
      all(samples %in% rownames(projection))
  ) {
    return(projection[samples, 1:2, drop = FALSE])
  }

  df <- as.data.frame(projection, stringsAsFactors = FALSE)

  # Sample key: row names if they look like sample identifiers (overlap the
  # sample names), else the first column whose values cover every sample.
  if (any(samples %in% rownames(df))) {
    keys <- rownames(df)
  } else {
    key_col <- NULL
    for (cn in colnames(df)) {
      if (all(samples %in% as.character(df[[cn]]))) {
        key_col <- cn
        break
      }
    }
    if (is.null(key_col)) {
      cli::cli_abort(c(
        "Could not match {.arg projection} to the sample names.",
        "i" = "Supply a matrix or data frame whose row names (or one column) hold the sample names."
      ))
    }
    keys <- as.character(df[[key_col]])
    df[[key_col]] <- NULL
  }

  # Coordinate columns: a known axis-name pair if present, else the first two
  # numeric columns.
  coord_pairs <- list(
    c("x_umap", "y_umap"),
    c("x", "y"),
    c("NMDS1", "NMDS2"),
    c("MDS1", "MDS2"),
    c("Axis.1", "Axis.2"),
    c("Dim1", "Dim2"),
    c("PC1", "PC2"),
    c("V1", "V2")
  )
  xy_cols <- NULL
  for (pair in coord_pairs) {
    if (all(pair %in% colnames(df))) {
      xy_cols <- pair
      break
    }
  }
  if (is.null(xy_cols)) {
    num_cols <- colnames(df)[vapply(df, is.numeric, logical(1))]
    if (length(num_cols) < 2) {
      cli::cli_abort(
        "{.arg projection} must provide at least two numeric coordinate columns."
      )
    }
    xy_cols <- num_cols[1:2]
  }

  coords <- as.matrix(df[, xy_cols, drop = FALSE])
  rownames(coords) <- keys
  missing_samples <- setdiff(samples, rownames(coords))
  if (length(missing_samples) > 0) {
    cli::cli_abort(c(
      "{.arg projection} is missing coordinates for {length(missing_samples)} sample{?s}.",
      "x" = "First missing: {.val {utils::head(missing_samples, 3)}}."
    ))
  }
  coords[samples, , drop = FALSE]
}

# Compute 2D coordinates for every node from an ecological projection of the
# samples. `projection` is either a vegdist method name (the samples are
# ordinated to 2D with `ordination`) or a precomputed sample projection (a
# matrix or data frame, see `as_sample_coords`). Taxa are placed at the
# abundance-weighted average of the samples they occur in (vegan::wascores),
# co-embedding them in the sample ordination space.
resolve_projection_coords <- function(
  projection,
  ordination,
  otu,
  samples,
  taxa,
  seed
) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg vegan} is required to use {.arg projection}.",
      "i" = "Install it with {.code install.packages(\"vegan\")}."
    ))
  }

  if (is.character(projection) && length(projection) == 1L) {
    if (!projection %in% vegdist_methods) {
      cli::cli_abort(c(
        "{.arg projection} {.val {projection}} is not a {.fn vegan::vegdist} method.",
        "i" = "Valid methods: {.val {vegdist_methods}}."
      ))
    }
    d <- vegan::vegdist(otu, method = projection)
    if (!is.null(seed)) {
      set.seed(seed)
    }
    sample_coords <- switch(
      ordination,
      nmds = {
        invisible(utils::capture.output(
          ord <- suppressMessages(suppressWarnings(
            vegan::metaMDS(d, k = 2, trace = 0)
          ))
        ))
        ord$points
      },
      pcoa = stats::cmdscale(d, k = 2)
    )
  } else if (is.matrix(projection) || is.data.frame(projection)) {
    sample_coords <- as_sample_coords(projection, samples)
  } else {
    cli::cli_abort(
      "{.arg projection} must be NULL, a {.fn vegan::vegdist} method name, or a matrix / data frame of sample coordinates."
    )
  }

  sample_coords <- sample_coords[samples, , drop = FALSE]
  storage.mode(sample_coords) <- "double"
  # cmdscale can collapse to one axis on degenerate data; pad to 2D.
  if (ncol(sample_coords) < 2) {
    sample_coords <- cbind(sample_coords, 0)[, 1:2, drop = FALSE]
  }
  taxa_coords <- vegan::wascores(sample_coords, otu)

  data.frame(
    name = c(samples, taxa),
    x = c(
      sample_coords[, 1],
      taxa_coords[match(taxa, rownames(taxa_coords)), 1]
    ),
    y = c(
      sample_coords[, 2],
      taxa_coords[match(taxa, rownames(taxa_coords)), 2]
    ),
    stringsAsFactors = FALSE
  )
}

################################################################################
#' Bipartite network of samples and taxa from a phyloseq object
#'
#' @description
#' \lifecycle{experimental}
#'
#' Draw a bipartite network in which both samples and taxa are nodes and
#' edges link a sample to a taxon whenever that taxon occurs in that sample
#' (i.e. the corresponding cell of the OTU table is greater than
#' `min_abundance`). Node positions are computed with `igraph` and the result
#' is rendered as a static `ggplot2` object (no `htmlwidgets`), so it can be
#' further customised with the usual `+` syntax. Inspired by Figure 2 of
#' Taudière et al. (2018) and the vda-lab/snowflake project.
#'
#' Sample-based and taxa-based metrics can be mapped to the node colour and
#' node size aesthetics.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   obtained using the `phyloseq` package.
#' @param taxa_color (default NULL) Name of a `tax_table` column used to colour
#'   the taxa nodes. When NULL, taxa nodes share a single colour.
#' @param sample_color (default "grey20") Controls the colour of the sample
#'   nodes. Either the name of a `sam_data` column to colour samples by — a
#'   numeric column gives a continuous viridis gradient, a factor or character
#'   column a discrete viridis palette, drawn on its own scale so it stays
#'   independent of `taxa_color` — or a single colour (name or hex) applied to
#'   every sample node. A value matching a `sam_data` column name is treated as
#'   a column, otherwise as a colour.
#' @param taxa_size (default "abundance") Metric mapped to the size of taxa
#'   nodes. One of "abundance" (total number of sequences), "prevalence"
#'   (number of samples in which the taxon occurs) or "none" (constant size).
#' @param sample_size (default "abundance") Metric mapped to the size of sample
#'   nodes. One of "abundance" (sequencing depth), "richness" (number of taxa
#'   observed) or "none" (constant size).
#' @param min_abundance (int, default 0) An edge is drawn between a sample and a
#'   taxon only when the corresponding OTU-table value is strictly greater than
#'   this threshold.
#' @param layout (default "fr") The `igraph` layout algorithm used to position
#'   the nodes when `projection` is `NULL`. One of "fr" (Fruchterman-Reingold),
#'   "kk" (Kamada-Kawai), "bipartite" (two parallel rows of nodes) or "nicely"
#'   (igraph's automatic choice). Ignored when `projection` is set.
#' @param projection (default NULL) Position the nodes so that distances among
#'   the sample nodes reflect community dissimilarity, instead of the abstract
#'   `igraph` `layout`. Either a [vegan::vegdist()] dissimilarity method name
#'   (e.g. "bray", "jaccard"; the samples are then ordinated to two dimensions,
#'   see `ordination`) or a precomputed sample projection. The latter may be a
#'   numeric matrix with sample row names, or a data frame / tibble in which the
#'   sample names are the row names or one column, and the coordinates are an
#'   axis-name pair (`x_umap`/`y_umap` — so [MiscMetabar::umap_pq()] output
#'   works directly —, `NMDS1`/`NMDS2`, `PC1`/`PC2`, `Axis.1`/`Axis.2`, ...) or
#'   the first two numeric columns. In every case each taxon node is placed at
#'   the abundance-weighted average of the samples it occurs in
#'   ([vegan::wascores()]). Requires the `vegan` package.
#' @param ordination (default "nmds") Ordination used to project the
#'   dissimilarity matrix to two dimensions when `projection` is a distance
#'   method name. One of "nmds" (non-metric multidimensional scaling, via
#'   [vegan::metaMDS()]) or "pcoa" (principal coordinates analysis, via
#'   [stats::cmdscale()]). Ignored when `projection` is a coordinate matrix.
#' @param edge_color (default "grey70") Colour of the edges. When `edge_weight`
#'   encodes the weight as colour ("color" or "both"), `edge_color` defines the
#'   gradient: a single colour gives a gradient from "grey70" (fewest sequences)
#'   to that colour (most sequences), and a length-2 vector gives a gradient
#'   between the two colours (low to high). When edges are not colour-weighted,
#'   the last value is used as the uniform edge colour.
#' @param edge_alpha (float, `[0:1]`, default 0.3) Opacity of the edges.
#' @param edge_weight (default "none") Weight the edges by their number of
#'   sequences (the OTU-table value of the sample-taxon pair). One of "none"
#'   (uniform edges), "linewidth" (thicker edges for more sequences, with a
#'   legend), "color" (darker-grey edges for more sequences) or "both".
#' @param edge_trans (default "log10") Transformation applied to the per-edge
#'   sequence count before it is mapped to the edge aesthetic(s). One of
#'   "log10" or "identity".
#' @param point_size (default 4) Base point size, used directly when the
#'   corresponding `*_size` argument is "none" and as the upper bound of the
#'   size scale otherwise.
#' @param sample_label (logical, default FALSE) Print the sample names next to
#'   the sample nodes. Either a single `TRUE`/`FALSE` (label all or none) or a
#'   logical vector the length of the number of samples in `physeq` (in
#'   [phyloseq::sample_names()] order) to label only a subset.
#' @param taxa_label (logical, default FALSE) Print the taxa names next to the
#'   taxa nodes. Either a single `TRUE`/`FALSE` (label all or none) or a logical
#'   vector the length of the number of taxa in `physeq` (in
#'   [phyloseq::taxa_names()] order) to label only a subset.
#' @param label_size_sample (default 3) Text size of the sample-node labels.
#' @param label_size_taxa (default 3) Text size of the taxa-node labels.
#' @param na_remove (logical, default FALSE) If TRUE, samples with NA in
#'   `sample_color` are removed before plotting.
#' @param seed (int, default NULL) Optional seed passed to [set.seed()] to make
#'   the (stochastic) layout reproducible.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @export
#' @author Adrien Taudière
#' @seealso [MiscMetabar::plot_tax_pq()]
#' @examples
#' if (requireNamespace("igraph")) {
#'   # Keep the plot readable: a handful of samples and the commonest taxa
#'   ps <- prune_samples(sample_names(data_fungi_mini)[1:20], data_fungi_mini)
#'   ps <- clean_pq(ps, silent = TRUE)
#'   bipartite_network_pq(ps, taxa_color = "Class", seed = 1)
#' }
#' \donttest{
#'   ps <- prune_samples(sample_names(data_fungi_mini)[1:20], data_fungi_mini)
#'   ps <- clean_pq(ps, silent = TRUE)
#'   # Sample nodes positioned by Bray-Curtis dissimilarity; each taxon sits at
#'   # the weighted average of the samples it occurs in.
#'   bipartite_network_pq(ps, projection = "bray", seed = 1)
#'
#' # Position the samples from a precomputed UMAP (needs the `umap` package).
#'   ps <- prune_samples(sample_names(data_fungi_mini)[1:20], data_fungi_mini)
#'   ps <- clean_pq(ps, silent = TRUE)
#'   df_umap <- MiscMetabar::umap_pq(ps, seed = 1)
#'   bipartite_network_pq(ps, projection = df_umap, sample_color = "Height", seed = 1)
#'
#'   bipartite_network_pq(ps, projection = df_umap, sample_color = ("Height"), seed = 1,
#'    taxa_color="Order",
#'    sample_label= as.vector(sample_sums(ps)>10000),
#'    label_size_sample = 2,
#'    taxa_label = taxa_sums(ps)>10000
#'   )
#' }

bipartite_network_pq <- function(
  physeq,
  taxa_color = NULL,
  sample_color = "grey20",
  taxa_size = "abundance",
  sample_size = "abundance",
  min_abundance = 0,
  layout = "fr",
  projection = NULL,
  ordination = "nmds",
  edge_color = "grey70",
  edge_alpha = 0.3,
  edge_weight = "none",
  edge_trans = "log10",
  point_size = 4,
  sample_label = FALSE,
  taxa_label = FALSE,
  label_size_sample = 3,
  label_size_taxa = 3,
  na_remove = FALSE,
  seed = NULL
) {
  taxa_size <- match.arg(taxa_size, c("abundance", "prevalence", "none"))
  sample_size <- match.arg(sample_size, c("abundance", "richness", "none"))
  layout <- match.arg(layout, c("fr", "kk", "bipartite", "nicely"))
  ordination <- match.arg(ordination, c("nmds", "pcoa"))
  edge_weight <- match.arg(edge_weight, c("none", "linewidth", "color", "both"))
  edge_trans <- match.arg(edge_trans, c("log10", "identity"))
  size_is_constant <- taxa_size == "none" && sample_size == "none"

  verify_pq(physeq)
  physeq <- clean_pq(physeq, silent = TRUE)

  # `sample_color` is either a sam_data column to colour samples by, or a
  # literal colour applied to every sample node.
  sample_color_is_column <- !is.null(sample_color) &&
    length(sample_color) == 1L &&
    sample_color %in% colnames(physeq@sam_data)

  if (na_remove && sample_color_is_column) {
    keep <- !is.na(as.vector(physeq@sam_data[[sample_color]]))
    physeq <- prune_samples(sample_names(physeq)[keep], physeq)
    physeq <- clean_pq(physeq, silent = TRUE)
  }

  # OTU matrix oriented as samples (rows) x taxa (columns)
  otu <- as(physeq@otu_table, "matrix")
  if (physeq@otu_table@taxa_are_rows) {
    otu <- t(otu)
  }

  samples <- rownames(otu)
  taxa <- colnames(otu)

  # Colour values: a sam_data column for samples (numeric or categorical) and a
  # tax_table rank for taxa. Each is mapped to its own scale at render time.
  if (sample_color_is_column) {
    sample_values <- physeq@sam_data[[sample_color]][
      match(samples, sample_names(physeq))
    ]
    sample_is_numeric <- is.numeric(sample_values)
  }
  if (!is.null(taxa_color)) {
    tt <- unclass(physeq@tax_table)
    taxa_values <- as.character(tt[match(taxa, rownames(tt)), taxa_color])
  }

  # Resolve the label selectors into the node names to annotate
  sample_label <- resolve_node_label(sample_label, samples, "sample_label")
  taxa_label <- resolve_node_label(taxa_label, taxa, "taxa_label")
  labelled_names <- c(samples[sample_label], taxa[taxa_label])

  # Edge list from the non-zero cells of the OTU table
  idx <- which(otu > min_abundance, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    stop(
      "No edge to draw: every OTU-table value is <= `min_abundance` (",
      min_abundance,
      ").",
      call. = FALSE
    )
  }
  edges <- data.frame(
    from = samples[idx[, "row"]],
    to = taxa[idx[, "col"]],
    weight = otu[idx],
    stringsAsFactors = FALSE
  )

  # Node table: type, optional colour and size metrics ------------------------
  nodes <- data.frame(
    name = c(samples, taxa),
    type = c(
      rep("Sample", length(samples)),
      rep("Taxa", length(taxa))
    ),
    stringsAsFactors = FALSE
  )

  sample_metric <- switch(
    sample_size,
    abundance = rowSums(otu),
    richness = rowSums(otu > 0),
    none = rep(NA_real_, length(samples))
  )
  taxa_metric <- switch(
    taxa_size,
    abundance = colSums(otu),
    prevalence = colSums(otu > 0),
    none = rep(NA_real_, length(taxa))
  )
  nodes$node_size <- c(
    sample_metric[match(samples, names(sample_metric))],
    taxa_metric[match(taxa, names(taxa_metric))]
  )
  # When only one side requests a metric, the other side is NA; give those
  # nodes a constant (median) size so the size scale does not silently drop
  # them when the plot is rendered.
  if (!size_is_constant && anyNA(nodes$node_size)) {
    nodes$node_size[is.na(nodes$node_size)] <- median(
      nodes$node_size,
      na.rm = TRUE
    )
  }
  if (size_is_constant) {
    nodes$node_size <- point_size
  }

  # Node coordinates ----------------------------------------------------------
  # Either an ecological projection (sample distances reflect community
  # dissimilarity) or an abstract igraph layout of the graph structure.
  if (!is.null(projection)) {
    coords <- resolve_projection_coords(
      projection,
      ordination,
      otu,
      samples,
      taxa,
      seed
    )
  } else {
    graph <- igraph::graph_from_data_frame(
      edges,
      directed = FALSE,
      vertices = nodes
    )
    if (!is.null(seed)) {
      set.seed(seed)
    }
    coords <- switch(
      layout,
      fr = igraph::layout_with_fr(graph),
      kk = igraph::layout_with_kk(graph),
      bipartite = igraph::layout_as_bipartite(
        graph,
        types = igraph::V(graph)$type == "Taxa"
      ),
      nicely = igraph::layout_nicely(graph)
    )
    coords <- as.data.frame(coords)
    names(coords) <- c("x", "y")
    coords$name <- igraph::V(graph)$name
  }
  nodes <- merge(nodes, coords, by = "name")

  # Attach node coordinates to each edge end ----------------------------------
  edges$x <- nodes$x[match(edges$from, nodes$name)]
  edges$y <- nodes$y[match(edges$from, nodes$name)]
  edges$xend <- nodes$x[match(edges$to, nodes$name)]
  edges$yend <- nodes$y[match(edges$to, nodes$name)]

  # ggplot --------------------------------------------------------------------
  # Samples are drawn with the `fill` aesthetic and taxa with `colour`, so the
  # two node types carry independent colour scales (e.g. a numeric sample
  # gradient alongside a discrete taxa palette).
  nodes_s <- nodes[nodes$type == "Sample", , drop = FALSE]
  nodes_t <- nodes[nodes$type == "Taxa", , drop = FALSE]
  if (sample_color_is_column) {
    nodes_s$col_val <- sample_values[match(nodes_s$name, samples)]
  }
  if (!is.null(taxa_color)) {
    nodes_t$col_val <- taxa_values[match(nodes_t$name, taxa)]
  }

  # Edge layer: optionally weight edges by their (transformed) sequence count,
  # encoded as line width and/or a grey gradient. The gradient colour is applied
  # directly (not through a scale) so it never collides with the taxa colour
  # scale.
  edges$w <- switch(
    edge_trans,
    log10 = log10(edges$weight),
    identity = edges$weight
  )
  use_edge_lw <- edge_weight %in% c("linewidth", "both")
  use_edge_col <- edge_weight %in% c("color", "both")

  if (use_edge_col) {
    # Gradient endpoints from `edge_color`: a length-2 vector sets both ends
    # (low to high), a single colour goes from "grey70" (low) to that colour.
    edge_stops <- if (length(edge_color) == 1L) {
      c("grey70", edge_color)
    } else {
      edge_color
    }
    rng <- range(edges$w, finite = TRUE)
    wn <- if (diff(rng) > 0) {
      (edges$w - rng[1]) / diff(rng)
    } else {
      rep(0.5, nrow(edges))
    }
    edge_pal <- grDevices::colorRampPalette(edge_stops)(101)
    edge_col_vec <- edge_pal[round(wn * 100) + 1]
  }

  if (use_edge_lw) {
    seg_aes <- ggplot2::aes(
      x = .data$x,
      y = .data$y,
      xend = .data$xend,
      yend = .data$yend,
      linewidth = .data$w
    )
  } else {
    seg_aes <- ggplot2::aes(
      x = .data$x,
      y = .data$y,
      xend = .data$xend,
      yend = .data$yend
    )
  }
  p <- ggplot2::ggplot() +
    do.call(
      ggplot2::geom_segment,
      list(
        data = edges,
        mapping = seg_aes,
        alpha = edge_alpha,
        color = if (use_edge_col) {
          edge_col_vec
        } else {
          edge_color[[length(edge_color)]]
        }
      )
    )
  if (use_edge_lw) {
    p <- p +
      ggplot2::scale_linewidth_continuous(
        range = c(0.2, 1.5),
        name = if (edge_trans == "log10") "Sequences (log10)" else "Sequences"
      )
  }

  # Taxa nodes (colour scale)
  if (!is.null(taxa_color)) {
    p <- p +
      ggplot2::geom_point(
        data = nodes_t,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          size = .data$node_size,
          shape = .data$type,
          color = .data$col_val
        )
      )
  } else {
    p <- p +
      ggplot2::geom_point(
        data = nodes_t,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          size = .data$node_size,
          shape = .data$type
        ),
        color = "grey50"
      )
  }

  # Sample nodes (fill scale)
  if (sample_color_is_column) {
    p <- p +
      ggplot2::geom_point(
        data = nodes_s,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          size = .data$node_size,
          shape = .data$type,
          fill = .data$col_val
        ),
        color = "grey30"
      )
  } else {
    p <- p +
      ggplot2::geom_point(
        data = nodes_s,
        ggplot2::aes(
          x = .data$x,
          y = .data$y,
          size = .data$node_size,
          shape = .data$type
        ),
        color = "grey30",
        fill = if (is.null(sample_color)) "grey20" else sample_color
      )
  }

  label_df <- nodes[nodes$name %in% labelled_names, , drop = FALSE]
  label_s <- label_df[label_df$type == "Sample", , drop = FALSE]
  label_t <- label_df[label_df$type == "Taxa", , drop = FALSE]
  if (nrow(label_s) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = label_s,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$name),
        size = label_size_sample,
        vjust = -0.6,
        show.legend = FALSE
      )
  }
  if (nrow(label_t) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = label_t,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$name),
        size = label_size_taxa,
        vjust = -0.6,
        show.legend = FALSE
      )
  }

  p <- p +
    ggplot2::scale_shape_manual(
      values = c(Sample = 22, Taxa = 16),
      name = "Node type"
    )

  if (size_is_constant) {
    p <- p + ggplot2::scale_size_identity()
  } else {
    p <- p +
      ggplot2::scale_size_continuous(range = c(1, point_size), name = "Metric")
  }

  if (sample_color_is_column) {
    fill_scale <- if (sample_is_numeric) {
      ggplot2::scale_fill_viridis_c(name = sample_color, na.value = "grey80")
    } else {
      ggplot2::scale_fill_viridis_d(name = sample_color, na.value = "grey80")
    }
    p <- p + fill_scale
  }
  if (!is.null(taxa_color)) {
    p <- p + ggplot2::labs(color = taxa_color)
  }

  p +
    ggplot2::coord_equal() +
    ggplot2::theme_void()
}
