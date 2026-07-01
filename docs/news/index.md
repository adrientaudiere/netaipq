# Changelog

## netaipq 0.1.0 (Development version)

## netaipq 0.0.0

- Initial development version of the package, providing networks and
  machine-learning tools for ‘phyloseq’ objects within the ‘pqverse’
  ecosystem.
- Add a “Get started with netaipq” vignette and a pkgdown website
  skeleton.
- [`bipartite_network_pq()`](https://adrientaudiere.github.io/netaipq/reference/bipartite_network_pq.md)
  draws a static ‘ggplot2’ bipartite network in which samples and taxa
  are nodes and edges come from the OTU table, with sample- and
  taxa-based metrics mappable to node colour and size, and optional
  sample and taxa name labels via `sample_label` and `taxa_label` (a
  single `TRUE`/`FALSE` or a per-node logical vector). `sample_color`
  accepts either a `sam_data` column to colour samples by (numeric gives
  a continuous gradient, factor or character a discrete palette, on a
  scale independent of `taxa_color`) or a single fixed colour. Node
  positions can reflect community dissimilarity through `projection` (a
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  method ordinated by `ordination`, or a precomputed 2D sample
  projection such as the output of an NMDS or
  [`MiscMetabar::umap_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/umap_pq.html)),
  with each taxon placed at the abundance-weighted average of the
  samples it occurs in. Edges can be weighted by their number of
  sequences (log10-transformed by default) and encoded as line width, a
  grey colour gradient, or both via `edge_weight` and `edge_trans`.
- [`cooccurrence_network_pq()`](https://adrientaudiere.github.io/netaipq/reference/cooccurrence_network_pq.md)
  computes the probabilistic species co-occurrence model of Veech (2013)
  on the taxa of a ‘phyloseq’ object via
  [`cooccur::cooccur()`](https://rdrr.io/pkg/cooccur/man/cooccur.html),
  converting the OTU table to presence/absence by default.
- [`dag_test_pq()`](https://adrientaudiere.github.io/netaipq/reference/dag_test_pq.md)
  confronts a hypothesised causal DAG with the (diversity-enriched)
  sample data of a ‘phyloseq’ object using
  [`dagitty::localTests()`](https://rdrr.io/pkg/dagitty/man/localTests.html).
- [`ggclusternet_pq()`](https://adrientaudiere.github.io/netaipq/reference/ggclusternet_pq.md)
  builds a per-group microbial correlation network with
  [`ggClusterNet::network.2()`](https://rdrr.io/pkg/ggClusterNet/man/network.2.html).
- [`netassoc_network_pq()`](https://adrientaudiere.github.io/netaipq/reference/netassoc_network_pq.md)
  infers a species-association network from a ‘phyloseq’ object against
  a null model of random assembly via
  [`netassoc::make_netassoc_network()`](https://rdrr.io/pkg/netassoc/man/make_netassoc_network.html).
- [`sparcev_pq()`](https://adrientaudiere.github.io/netaipq/reference/sparcev_pq.md)
  correlates each taxon with a continuous environmental variable while
  accounting for compositionality via
  [`CompoCor::SparCEV()`](https://rdrr.io/pkg/CompoCor/man/SparCEV.html).
