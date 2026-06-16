# Species-association network of a phyloseq object using netassoc

**\[experimental\]**

Infer a species-association network from a `phyloseq` object with
[`netassoc::make_netassoc_network()`](https://rdrr.io/pkg/netassoc/man/make_netassoc_network.html).
The method (Morueta-Holme et al. 2016) compares the observed
co-occurrence structure to a null model of random assembly and keeps
only the taxa pairs whose association departs significantly from that
expectation, returning a partial-correlation network of the taxa.

## Usage

``` r
netassoc_network_pq(
  physeq,
  numnulls = 100,
  method = "partial_correlation",
  plot = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- numnulls:

  (int, default 100) Number of null-model replicates used to build the
  expectation. Increase for real analyses; the default keeps the
  examples and tests fast.

- method:

  (default "partial_correlation") Association measure passed to
  [`netassoc::make_netassoc_network()`](https://rdrr.io/pkg/netassoc/man/make_netassoc_network.html).

- plot:

  (logical, default FALSE) Draw the heatmaps produced by
  [`netassoc::make_netassoc_network()`](https://rdrr.io/pkg/netassoc/man/make_netassoc_network.html).
  FALSE by default to avoid a plotting side effect; the returned object
  still contains the networks.

- verbose:

  (logical, default FALSE) Passed to
  [`netassoc::make_netassoc_network()`](https://rdrr.io/pkg/netassoc/man/make_netassoc_network.html).

- ...:

  Other parameters passed on to
  [`netassoc::make_netassoc_network()`](https://rdrr.io/pkg/netassoc/man/make_netassoc_network.html).

## Value

The list returned by
[`netassoc::make_netassoc_network()`](https://rdrr.io/pkg/netassoc/man/make_netassoc_network.html).
Its `network_all` element is an `igraph` graph of the significant
associations.

## References

Morueta-Holme, N. et al. (2016) A network approach for inferring species
associations from co-occurrence data.
[doi:10.1111/ecog.01892](https://doi.org/10.1111/ecog.01892) .

## See also

[`cooccurrence_network_pq()`](https://adrientaudiere.github.io/netaipq/reference/cooccurrence_network_pq.md),
[`bipartite_network_pq()`](https://adrientaudiere.github.io/netaipq/reference/bipartite_network_pq.md)

## Author

Adrien Taudière

## Examples

``` r
if (requireNamespace("netassoc", quietly = TRUE)) {
  ps <- prune_samples(sample_names(data_fungi_mini)[1:5], data_fungi_mini)
  ps <- clean_pq(ps, silent = TRUE)
  top_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
  ps <- prune_taxa(top_taxa, ps)
  res <- netassoc_network_pq(ps, numnulls = 20)

  netassoc::plot_netassoc_network(res$network_all)
}
#> .
#> Warning: 1 instances of variables with zero scale detected!
#> Warning: 1 instances of variables with zero scale detected!
#> ..
#> Warning: 1 instances of variables with zero scale detected!
#> Warning: 1 instances of variables with zero scale detected!
#> ....
#> Warning: 1 instances of variables with zero scale detected!
#> Warning: 1 instances of variables with zero scale detected!
#> ..
#> Warning: 1 instances of variables with zero scale detected!
#> Warning: 1 instances of variables with zero scale detected!
#> ....
#> Warning: 1 instances of variables with zero scale detected!
#> Warning: 1 instances of variables with zero scale detected!
#> ......
#> Warning: 1 instances of variables with zero scale detected!
#> Warning: 1 instances of variables with zero scale detected!
#> .
#> Warning: Non-positive edge weight found, ignoring all weights during graph layout.
```
