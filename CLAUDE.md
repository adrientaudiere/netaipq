# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository — the **netaipq** sub-package of the pqverse.

## Package Overview

**netaipq** is the networks and machine-learning layer of the pqverse. It is
the home for *any* feature that adds a brand-new ML, network, graph,
causal-inference, or statistical-analysis dependency, keeping 'MiscMetabar'
lean. It also owns genuinely new analysis methods (alpha-diversity
estimation, dissimilarity modelling, spatial partitioning, indicator-taxa
analysis).

**Scope guard** (from `ROADMAP.md`):

- ✅ In scope: co-occurrence / association networks (cooccur, SpiecEasi,
  NetCoMi, rMAGMA), bipartite-network viz, DAGs (dagitty), supervised
  classification (PLS-DA via mixOmics), random-forest / Boruta / SHAP,
  indicator species & differential abundance (LefSe, ALDEx2 senAnalysis,
  adiv::dbMANOVAspecies, TITAN2), alpha-diversity estimation (DivNet +
  breakaway), robust ordination (rpca), dissimilarity & spatial methods
  (zetadiv, GDM, dbMEM), null-model assembly (NST, iCAMP).
- ❌ Out of scope: pure ggplot2 helpers (→ `ggplotpq`), multi-phyloseq
  comparators (→ `comparpq`), data-structure utilities (→ `tidypq`).

**Dependency rule.** netaipq absorbs heavy analysis/ML/network deps. Render
the resulting plots through `ggplotpq` where a pure-viz layer already exists;
keep the analysis wrapper here.

## Common Commands

```bash
# Run code with loaded package
Rscript -e "devtools::load_all(); code"

# Run all tests
Rscript -e "devtools::test()"

# Run tests for files starting with {name}
Rscript -e "devtools::test(filter = '^{name}')"

# Generate documentation
Rscript -e "devtools::document()"

# Full package check
Rscript -e "devtools::check()"
```

## Coding Conventions

- Use base pipe (`|>`) not magrittr (`%>%`)
- Use `function() {}` for anonymous functions (not `\()` for multi-statement)
- Tests for `R/{name}.R` go in `tests/testthat/test_{name}.R` (underscore)
- Every user-facing function must be exported with full roxygen2
  documentation (`@param`, `@return`, `@export`, `@examples`, `@author`)
- Wrap roxygen comments at 80 characters
- CRAN example constraints: primary example in `\donttest{}`, variants in
  `\dontrun{}`; cap per-sample work at 5 samples via
  `prune_samples(sample_names(data_fungi_mini)[1:5], data_fungi_mini)` and
  filter to the most abundant taxa for network/ordination examples
- Guard every Suggests-package call with `requireNamespace()` +
  `cli::cli_abort()`
- Air format the package: `air format .` (then scope the diff — revert
  incidental reformats to unrelated files)

## Cross-references

- Workspace CLAUDE.md: `pqverse/CLAUDE.md` (overall context)
- ROADMAP section: <https://github.com/adrientaudiere/pqverse/ROADMAP.md#netaipq>
- Sister packages: `pqverse_pkg/MiscMetabar/`, `pqverse_pkg/ggplotpq/`,
  `pqverse_pkg/comparpq/`, `pqverse_pkg/tidypq/`

## Agent skills

### Issue tracker

Issues and PRDs are tracked as GitHub issues via the `gh` CLI; external PRs are not a triage surface. See `docs/agents/issue-tracker.md`.

### Triage labels

Uses the five canonical triage labels (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`). See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: one `CONTEXT.md` + `docs/adr/` at the repo root. See `docs/agents/domain.md`.
