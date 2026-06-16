skip_on_cran()
library(MiscMetabar)
data(data_fungi_mini)

mini_net <- function() {
  ps <- prune_samples(sample_names(data_fungi_mini)[1:5], data_fungi_mini)
  ps <- clean_pq(ps, silent = TRUE)
  top_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
  ps <- prune_taxa(top_taxa, ps)
  clean_pq(ps, silent = TRUE)
}

test_that("bipartite_network_pq returns a ggplot that renders", {
  skip_if_not_installed("igraph")
  p <- bipartite_network_pq(mini_net(), seed = 1)
  expect_s3_class(p, "ggplot")
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("bipartite_network_pq keeps all nodes when only one side is sized", {
  skip_if_not_installed("igraph")
  p <- bipartite_network_pq(
    mini_net(),
    taxa_size = "none",
    sample_size = "abundance",
    seed = 1
  )
  # No node dropped by the size scale (NA sizes would warn here)
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("bipartite_network_pq maps colour and size metrics", {
  skip_if_not_installed("igraph")
  p <- bipartite_network_pq(
    mini_net(),
    taxa_color = "Phylum",
    taxa_size = "prevalence",
    sample_size = "richness",
    seed = 1
  )
  expect_s3_class(p, "ggplot")
})

test_that("bipartite_network_pq accepts constant node size", {
  skip_if_not_installed("igraph")
  p <- bipartite_network_pq(
    mini_net(),
    taxa_size = "none",
    sample_size = "none",
    layout = "bipartite",
    seed = 1
  )
  expect_s3_class(p, "ggplot")
})

test_that("bipartite_network_pq errors when no edge passes the threshold", {
  skip_if_not_installed("igraph")
  expect_error(
    bipartite_network_pq(mini_net(), min_abundance = 1e9),
    "No edge to draw"
  )
})

test_that("bipartite_network_pq labels nodes with TRUE and a logical vector", {
  skip_if_not_installed("igraph")
  ps <- mini_net()

  p_all <- bipartite_network_pq(
    ps,
    sample_label = TRUE,
    taxa_label = TRUE,
    seed = 1
  )
  # Sample and taxa labels live in two separate geom_text layers.
  labels_all <- unlist(lapply(
    seq_along(p_all$layers),
    function(i) ggplot2::layer_data(p_all, i)$label
  ))
  expect_setequal(labels_all, c(sample_names(ps), taxa_names(ps)))

  # A logical vector labels only the selected taxa, no samples
  sel <- c(TRUE, rep(FALSE, ntaxa(ps) - 1))
  p_sub <- bipartite_network_pq(
    ps,
    taxa_label = sel,
    sample_label = FALSE,
    seed = 1
  )
  labels_sub <- ggplot2::layer_data(p_sub, length(p_sub$layers))$label
  expect_equal(labels_sub, taxa_names(ps)[1])
})

test_that("bipartite_network_pq colours samples by a numeric sam_data column", {
  skip_if_not_installed("igraph")
  ps <- mini_net()
  sample_data(ps)$num_x <- seq_len(nsamples(ps))
  p <- bipartite_network_pq(ps, sample_color = "num_x", seed = 1)
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$scales$get_scales("fill"), "ScaleContinuous"))
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("bipartite_network_pq colours samples by a categorical column", {
  skip_if_not_installed("igraph")
  p <- bipartite_network_pq(mini_net(), sample_color = "Height", seed = 1)
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$scales$get_scales("fill"), "ScaleDiscrete"))
})

test_that("bipartite_network_pq fills NA-valued samples with a visible colour", {
  skip_if_not_installed("igraph")
  ps <- mini_net()
  skip_if(!anyNA(sample_data(ps)$Height))
  p <- bipartite_network_pq(ps, sample_color = "Height", seed = 1)
  layers <- ggplot2::ggplot_build(p)$data
  squares <- Filter(
    function(d) "fill" %in% names(d) && any(d$shape == 22),
    layers
  )[[1]]
  # NA-valued samples must be filled (grey80), not left as NA (invisible).
  expect_false(anyNA(squares$fill))
})

test_that("bipartite_network_pq treats a non-column value as a literal colour", {
  skip_if_not_installed("igraph")
  p <- bipartite_network_pq(mini_net(), sample_color = "red", seed = 1)
  expect_s3_class(p, "ggplot")
  expect_null(p$scales$get_scales("fill"))
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("bipartite_network_pq keeps sample and taxa colour scales separate", {
  skip_if_not_installed("igraph")
  ps <- mini_net()
  sample_data(ps)$num_x <- seq_len(nsamples(ps))
  p <- bipartite_network_pq(
    ps,
    sample_color = "num_x",
    taxa_color = "Phylum",
    seed = 1
  )
  # Samples carry an explicit continuous fill scale, independent of the taxa
  # colour aesthetic.
  expect_true(inherits(p$scales$get_scales("fill"), "ScaleContinuous"))
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("bipartite_network_pq na_remove works with a literal sample colour", {
  skip_if_not_installed("igraph")
  p <- bipartite_network_pq(
    mini_net(),
    sample_color = "grey20",
    na_remove = TRUE,
    seed = 1
  )
  expect_s3_class(p, "ggplot")
})

test_that("bipartite_network_pq weights edges by sequence count (linewidth)", {
  skip_if_not_installed("igraph")
  p <- bipartite_network_pq(mini_net(), edge_weight = "linewidth", seed = 1)
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$scales$get_scales("linewidth"), "ScaleContinuous"))
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("bipartite_network_pq weights edges by colour without scale clash", {
  skip_if_not_installed("igraph")
  p <- bipartite_network_pq(
    mini_net(),
    edge_weight = "color",
    edge_color = "darkred",
    taxa_color = "Phylum",
    seed = 1
  )
  expect_no_warning(ggplot2::ggplot_build(p))
  seg <- ggplot2::layer_data(p, 1)
  # Edges carry more than one colour (the sequence-count gradient).
  expect_gt(length(unique(seg$colour)), 1)
})

test_that("bipartite_network_pq edge colour gradient uses edge_color", {
  skip_if_not_installed("igraph")
  ps <- mini_net()

  # A single colour -> gradient from grey70 (low) to that colour (high).
  p1 <- bipartite_network_pq(
    ps,
    edge_weight = "color",
    edge_color = "red",
    seed = 1
  )
  seg1 <- ggplot2::layer_data(p1, 1)$colour
  pal1 <- grDevices::colorRampPalette(c("grey70", "red"))(101)
  expect_true(pal1[1] %in% seg1)
  expect_true(pal1[101] %in% seg1)

  # A length-2 vector -> gradient between the two colours.
  p2 <- bipartite_network_pq(
    ps,
    edge_weight = "color",
    edge_color = c("blue", "gold"),
    seed = 1
  )
  seg2 <- ggplot2::layer_data(p2, 1)$colour
  pal2 <- grDevices::colorRampPalette(c("blue", "gold"))(101)
  expect_true(pal2[1] %in% seg2 && pal2[101] %in% seg2)

  # Non-weighted edges with a 2-vector use the last colour, uniformly.
  p3 <- bipartite_network_pq(
    ps,
    edge_weight = "none",
    edge_color = c("blue", "gold"),
    seed = 1
  )
  expect_equal(unique(ggplot2::layer_data(p3, 1)$colour), "gold")
})

test_that("bipartite_network_pq edge_weight 'both' maps width and colour", {
  skip_if_not_installed("igraph")
  p <- bipartite_network_pq(
    mini_net(),
    edge_weight = "both",
    edge_color = "navy",
    edge_trans = "identity",
    seed = 1
  )
  expect_true(inherits(p$scales$get_scales("linewidth"), "ScaleContinuous"))
  seg <- ggplot2::layer_data(p, 1)
  expect_gt(length(unique(seg$colour)), 1)
  expect_gt(length(unique(seg$linewidth)), 1)
})

test_that("bipartite_network_pq rejects a mis-sized label vector", {
  skip_if_not_installed("igraph")
  expect_error(
    bipartite_network_pq(mini_net(), taxa_label = c(TRUE, FALSE)),
    "logical vector of length"
  )
})

test_that("bipartite_network_pq positions nodes from a precomputed projection", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("vegan")
  ps <- mini_net()
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  coords <- stats::cmdscale(vegan::vegdist(otu, "bray"), k = 2)
  p <- bipartite_network_pq(ps, projection = coords)
  expect_s3_class(p, "ggplot")
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("bipartite_network_pq projects from a vegdist method name", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("vegan")
  p <- bipartite_network_pq(
    mini_net(),
    projection = "bray",
    ordination = "pcoa",
    seed = 1
  )
  expect_s3_class(p, "ggplot")
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("bipartite_network_pq projects via NMDS", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("vegan")
  ps <- prune_samples(sample_names(data_fungi_mini)[1:15], data_fungi_mini)
  ps <- clean_pq(ps, silent = TRUE)
  ps <- prune_taxa(names(sort(taxa_sums(ps), decreasing = TRUE))[1:30], ps)
  ps <- clean_pq(ps, silent = TRUE)
  p <- suppressWarnings(
    bipartite_network_pq(ps, projection = "bray", ordination = "nmds", seed = 1)
  )
  expect_s3_class(p, "ggplot")
})

test_that("bipartite_network_pq accepts a tibble projection (umap_pq shape)", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("vegan")
  ps <- mini_net()
  # Mimic MiscMetabar::umap_pq(): coords in x_umap/y_umap, sample names in a
  # column, no row names.
  df <- data.frame(
    x_umap = stats::rnorm(nsamples(ps)),
    y_umap = stats::rnorm(nsamples(ps)),
    Sample = sample_names(ps),
    stringsAsFactors = FALSE
  )
  p <- bipartite_network_pq(ps, projection = df, seed = 1)
  expect_s3_class(p, "ggplot")
  expect_no_warning(ggplot2::ggplot_build(p))
})

test_that("bipartite_network_pq projects from a live umap_pq() result", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("umap")
  # UMAP needs n_neighbors (default 15) < n samples, so use a larger subset.
  ps <- prune_samples(sample_names(data_fungi_mini)[1:20], data_fungi_mini)
  ps <- clean_pq(ps, silent = TRUE)
  df_umap <- suppressWarnings(suppressMessages(
    MiscMetabar::umap_pq(ps, seed = 1)
  ))
  p <- bipartite_network_pq(ps, projection = df_umap, seed = 1)
  expect_s3_class(p, "ggplot")
})

test_that("bipartite_network_pq rejects an unknown projection distance", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("vegan")
  expect_error(
    bipartite_network_pq(mini_net(), projection = "not_a_distance"),
    "vegdist"
  )
})

test_that("bipartite_network_pq errors when projection lacks a sample", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("vegan")
  ps <- mini_net()
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  coords <- stats::cmdscale(vegan::vegdist(otu, "bray"), k = 2)
  expect_error(
    bipartite_network_pq(ps, projection = coords[-1, ]),
    "missing coordinates"
  )
})
