# netaipq

`netaipq` is the [pqverse](https://github.com/adrientaudiere) home for
**networks and machine-learning** analyses of
[`phyloseq`](https://joey711.github.io/phyloseq/) objects: co-occurrence
and association networks, bipartite-network visualisation, supervised
classification, indicator-species and differential-abundance methods,
and other *genuinely new analysis methods* that go beyond a simple
`ggplot2` wrapper.

By absorbing any feature that would add a brand-new machine-learning,
network, graph, or statistical-analysis dependency, `netaipq` keeps
[`MiscMetabar`](https://adrientaudiere.github.io/MiscMetabar/) lean.
Pure `ggplot2` visualisation belongs in `ggplotpq`; cross-study
comparators belong in `comparpq`.

## Installation

You can install the development version of `netaipq` from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("adrientaudiere/netaipq")
```

## Functions

| Function                                                                                               | Description                                            |
|--------------------------------------------------------------------------------------------------------|--------------------------------------------------------|
| [`bipartite_network_pq()`](https://adrientaudiere.github.io/netaipq/reference/bipartite_network_pq.md) | Static `ggplot2` bipartite network of samples and taxa |

## Example

``` r
library(netaipq)
data(data_fungi_mini)

# A readable subset: a few samples and the commonest taxa
ps <- prune_samples(sample_names(data_fungi_mini)[1:5], data_fungi_mini)
ps <- clean_pq(ps, silent = TRUE)
top_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
ps <- prune_taxa(top_taxa, ps)

bipartite_network_pq(ps, taxa_color = "Phylum", seed = 1)
```
