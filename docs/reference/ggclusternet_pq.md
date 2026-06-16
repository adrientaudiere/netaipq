# Microbial correlation network of a phyloseq object using ggClusterNet

**\[experimental\]**

Build a per-group microbial correlation network from a `phyloseq` object
with
[`ggClusterNet::network.2()`](https://rdrr.io/pkg/ggClusterNet/man/network.2.html)
(Wen et al. 2022). For each level of `group`, taxa are correlated,
thresholded on correlation strength and significance, laid out with one
of the `ggClusterNet` network layouts and returned together with a
`ggplot2` rendering.

## Usage

``` r
ggclusternet_pq(
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
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- group:

  (required) Name of the `sam_data` column defining the groups for which
  a separate network is built.

- n:

  (int, default 500) Number of most-abundant taxa kept before the
  analysis (the `N` argument of
  [`ggClusterNet::network.2()`](https://rdrr.io/pkg/ggClusterNet/man/network.2.html)).

- r_threshold:

  (float, default 0.6) Minimum absolute correlation for an edge to be
  kept.

- p_threshold:

  (float, default 0.05) Maximum p-value for an edge to be kept.

- maxnode:

  (int, default 5) Maximum node size in the layout.

- layout_net:

  (default "model_maptree2") The `ggClusterNet` layout algorithm.

- big:

  (logical, default TRUE) Use the large-network code path of
  [`ggClusterNet::network.2()`](https://rdrr.io/pkg/ggClusterNet/man/network.2.html).

- select_layout:

  (logical, default TRUE) Let `ggClusterNet` adjust the layout.

- label:

  (logical, default FALSE) Draw taxa labels on the network.

- zipi:

  (logical, default FALSE) Compute the within/among-module connectivity
  (Zi-Pi) classification.

- ...:

  Other parameters passed on to
  [`ggClusterNet::network.2()`](https://rdrr.io/pkg/ggClusterNet/man/network.2.html).

## Value

The list returned by
[`ggClusterNet::network.2()`](https://rdrr.io/pkg/ggClusterNet/man/network.2.html),
whose first element is a `ggplot2` network plot and whose remaining
elements hold the network and correlation data.

## References

Wen, T. et al. (2022) ggClusterNet: an R package for microbiome network
analysis and modularity-based multiple network layouts.
[doi:10.1002/imt2.32](https://doi.org/10.1002/imt2.32) .

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
# `ggClusterNet` is a GitHub package and needs WGCNA (Bioconductor):
# remotes::install_github("taowenmicro/ggClusterNet")
# pak::pkg_install("WGCNA")
res <- ggclusternet_pq(data_fungi_mini, group = "Height", n = 200)
res[[1]]
} # }
```
