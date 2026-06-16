# Compositional taxa-environment correlation with SparCEV

**\[experimental\]**

Compute the sparse correlation of each taxon with a continuous
environmental (sample) variable using
[`CompoCor::SparCEV()`](https://rdrr.io/pkg/CompoCor/man/SparCEV.html)
(Sparse Correlation with an Environmental Variable). SparCEV accounts
for the compositional nature of metabarcoding count data and flags the
taxa whose abundance is significantly correlated with the variable.

## Usage

``` r
sparcev_pq(physeq, variable, add_tax_table = TRUE, ...)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- variable:

  (required) Name of the continuous variable in the `sam_data` slot to
  correlate the taxa with.

- add_tax_table:

  (logical, default TRUE) Join the `tax_table` of `physeq` to the result
  so each taxon carries its taxonomy.

- ...:

  Other parameters passed on to
  [`CompoCor::SparCEV()`](https://rdrr.io/pkg/CompoCor/man/SparCEV.html).

## Value

A data frame with one row per taxon: `taxa` (taxon name), `cor` (the
SparCEV correlation with `variable`) and `correlated` (logical, TRUE
when the correlation is deemed significant). When `add_tax_table = TRUE`
the columns of the `tax_table` are appended. The rendering of this
result (e.g. a labelled scatter plot) belongs in `ggplotpq`.

## See also

[`cooccurrence_network_pq()`](https://adrientaudiere.github.io/netaipq/reference/cooccurrence_network_pq.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
# `CompoCor` is a GitHub package:
# remotes::install_github("IbTJensen/CompoCor")
ps <- subset_samples(data_fungi_mini, !is.na(Time))
res <- sparcev_pq(ps, variable = "Time")
 res |>
   filter(correlated) |>
   ggplot(aes(x = cor, y = Family, shape = Guild, color = Order)) +
   geom_point() +
   geom_text(aes(label = taxa), nudge_y = 0.3, size = 2.5)
} # }
```
