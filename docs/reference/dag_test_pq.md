# Test a causal DAG against phyloseq sample data

**\[experimental\]**

Confront a hypothesised causal directed acyclic graph (DAG) with the
sample data of a `phyloseq` object. The implied conditional
independencies of the DAG are tested against the data with
[`dagitty::localTests()`](https://rdrr.io/pkg/dagitty/man/localTests.html):
each independence that the DAG claims is checked, and a large estimate
signals a DAG that is inconsistent with the data.

Diversity descriptors of the samples (`nb_seq`, the sequencing depth,
and `nb_otu`, the observed richness) are added to the sample data with
[`MiscMetabar::add_info_to_sam_data()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_info_to_sam_data.html)
so they can be used as DAG nodes alongside the columns of `sam_data`.

## Usage

``` r
dag_test_pq(
  physeq,
  dag,
  as_factor = NULL,
  type = "cis.chisq",
  na_remove = TRUE,
  ...
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- dag:

  (required) A causal DAG, either a `dagitty` object or a `dagitty`
  model string (e.g. `"dag { Time -> nb_otu nb_seq -> nb_otu }"`). Every
  node of the DAG must be a column of the (enriched) sample data.

- as_factor:

  (character vector, default NULL) Names of DAG variables to coerce to
  factor before the test (useful for categorical metadata stored as
  numbers or characters).

- type:

  (default "cis.chisq") The independence-test type passed to
  [`dagitty::localTests()`](https://rdrr.io/pkg/dagitty/man/localTests.html).
  See its help for the available options.

- na_remove:

  (logical, default TRUE) Remove samples with a missing value in any DAG
  variable before the test.

- ...:

  Other parameters passed on to
  [`dagitty::localTests()`](https://rdrr.io/pkg/dagitty/man/localTests.html).

## Value

The data frame returned by
[`dagitty::localTests()`](https://rdrr.io/pkg/dagitty/man/localTests.html),
one row per implied conditional independence, with the estimate and its
confidence interval. Pass it to
[`dagitty::plotLocalTestResults()`](https://rdrr.io/pkg/dagitty/man/plotLocalTestResults.html)
to visualise it.

## References

Ankan, A., Wortel, I. M. N. & Textor, J. (2021) Testing graphical causal
models using the R package "dagitty".
[doi:10.1002/cpz1.45](https://doi.org/10.1002/cpz1.45) .

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("dagitty", quietly = TRUE)) {
  dag <- "dag {
    Time -> nb_otu
    nb_seq -> nb_otu
    Time -> Height
    Height -> nb_otu
  }"
  dag_test_pq(data_fungi_mini, dag, as_factor = c("Time", "Height"))
}
#>                     rmsea       x2  df   p.value rmsea 2.5% rmsea 97.5%
#> Hght _||_ nb_s 0.00000000 123.3229 124 0.5002821 0.08891009   0.1397312
#> Time _||_ nb_s 0.02776546 196.6110 186 0.2828356 0.09810813   0.1389691
# }
if (FALSE) { # \dontrun{
if (requireNamespace("dagitty", quietly = TRUE)) {
  # Add Hill number q = 2 (inverse Simpson) to sample data
  hill2 <- phyloseq::estimate_richness(data_fungi_mini, measures = "InvSimpson")
  phyloseq::sample_data(data_fungi_mini)$hill_2 <- hill2$InvSimpson
  dag2 <- "dag {
    Time -> hill_2
    nb_seq -> hill_2
    Diameter -> hill_2
    Time -> Diameter
  }"
  dag_test_pq(data_fungi_mini, dag2, as_factor = "Time")
}
} # }
```
