# Analysis runs

One folder per *version* of the fuzzy-RDD analysis. A version is fully
identified by five choices, which the folder name encodes:

```
instr-<instrument>_w<window>_trt-<treatment>_gap<score_gap_min>_illib<illiberal_cutoff>
```

| Part | Meaning |
|---|---|
| `instr-` | which V-Party score decides who the "illiberal" side of the top 2 is: `antiplural`, `popul`, `econleft` (negated `v2pariglef`, so economic LEFT is the high end), `anteli`, `galtan` |
| `w` | post-election window N, in years. Treatment *and* every outcome are measured over `(election_year - 1, election_year + N]` |
| `trt-` | `ert` (ERT autocratization episode), `ertOrDdcg` (that OR an Acemoglu et al. reversal), `polyarchy` (continuous decline in V-Dem polyarchy) |
| `gap` | minimum `score_gap_z`; `any` = no restriction |
| `illib` | minimum `illiberal_score`; `any` = no restriction |
| `_exclyr` | present only when `TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR = FALSE` (see below) |

## Sample-restriction thresholds

`SCORE_GAP_MIN` and `ILLIBERAL_CUTOFF` each accept three forms, and which one
you mean is decided by what you write — there is no separate "type" switch:

| Written as | Read as | Example |
|---|---|---|
| a number | absolute value on the variable's own scale | `0.6` |
| `"qNN"` | percentile of the variable's distribution | `"q50"`, `"q97.5"` |
| `-Inf` | no restriction | |

Anything else — `"Q50"`, `"p50"`, `"q 50"`, `NA`, `"q150"` — is a **hard error**,
not a fallback. An earlier version fell back to `as.numeric()`, so a typo became
`NA`, `is.finite(NA)` was `FALSE`, and the restriction was skipped silently: the
run reported itself as restricted and wasn't.

**Use the quantile form whenever you compare across instruments.** They are not
on a common scale:

| instrument | range of `illiberal_score` | median | what an absolute `0.6` keeps |
|---|---|---|---|
| `v2xpa_antiplural` | 0.02 – 1.00 | 0.57 | 48% |
| `v2xpa_popul` | 0.06 – 0.99 | 0.51 | 39% |
| `v2pariglef_neg` | −1.92 – 3.81 | 0.82 | 60% |
| `v2paanteli` | −2.37 – 4.41 | 0.33 | 42% |
| `ep_galtan` | 4.50 – 9.41 | 7.00 | **100% — a no-op** |

So one absolute number means five different things, and on `ep_galtan` it sits
below the whole range and restricts nothing. `"q50"` means the same *thing*
everywhere: the top half on whichever scale that instrument uses.

Every run prints a `Sample restrictions` block naming the spec, how it was read,
the absolute value it resolved to, and how many elections it cost:

```
Sample restrictions:
  SCORE_GAP_MIN    no restriction (1347 rows kept)
  ILLIBERAL_CUTOFF q50 = the 50th percentile, which is 0.33 on this instrument's scale
                     ->  keep illiberal_score > 0.33:  1347 -> 673 rows (674 dropped)
  FINAL SAMPLE     673 elections
```

Both thresholds are resolved against the **full** loaded build before either
filter is applied, so the two axes stay independent — otherwise changing
`SCORE_GAP_MIN` would silently move an `ILLIBERAL_CUTOFF` of `"q50"` as well.

Folder names use the **resolved absolute value** (`..._illib0p33`, not
`..._illibq50`), so a `"q50"` run and a hand-written run at the same resolved
number share a folder instead of duplicating. The spec as written is recorded
next to it in `run_config.csv` and `manifest.csv` as
`illiberal_cutoff_spec` / `illiberal_cutoff_form`, and appears in every table
subtitle — the folder name is never the only record.

## The treatment window convention

`TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR` (default `TRUE`) decides whether the
post-election window opens in the election year or the year after:

| | treatment window | treated (ERT, w=5) |
|---|---|---|
| `TRUE` (default) | `[election_year, election_year + N]` | 173 |
| `FALSE` (legacy) | `(election_year, election_year + N]` | 142 |

`TRUE` is the default because it is what aligns treatment with the outcomes:
`window_change()` measures every outcome from `election_year - 1` to
`election_year + N`, so the outcome window always spanned the election year,
while under `FALSE` the treatment window did not.

`FALSE` is kept so the earlier results reproduce, which they do exactly — the
unrestricted first stage is `-0.025 (0.058)`, bw 21.08, and the
`illiberal_score > 0.6` run is `0.181 (0.103)*`, bw 20.86.

Whichever way it is set, `prior_backsliding` is anchored to the same window
start, so the pre- and post-election windows stay a partition with no overlap
and no gap. Builds and run folders under the two conventions are kept apart by
the `_exclyr` suffix, and each build records its own convention as an attribute
that `12_rdd_analysis.R` checks before using it.

Note that `polyarchy_decline` is unaffected by the toggle: it is computed by
`window_change()` over `(election_year - 1, election_year + N]` and so was
always on the aligned window.

`manifest.csv` indexes every run with its configuration, sample size and
first-stage result. Each folder also carries its own `run_config.csv`.

## Folder contents

```
rdd_results.csv              one row per outcome (reduced form + fuzzy LATE)
rdd_first_stage_results.csv  the first stage
first_stage_table.html       the same, econ-paper formatted
outcomes_table.html
plots/first_stage.png
plots/running_var_density.png
plots/outcomes_<panel>.png   one figure per outcome panel
plots/outcomes_all.png       all panels on one sheet
```

## Runs vs sweeps

A **run** is one estimated specification: a single value for every one of the
five fields above. It lives in `output/runs/<slug>/` and records itself in
`run_config.csv`.

A **sweep** varies one or more of those fields and pools the results. It lives
in `output/runs/_sweeps/<name>/` and records itself in `sweep_config.csv`, which
lists what was held *fixed* and what was *swept*, with every level of every
swept axis.

The distinction matters because a `run_config.csv` asserts one value per field.
A sweep writing one would be asserting a single threshold for a grid that varies
it — so sweeps never call `run_dir()`. (`13_restriction_grid.R` used to, with
`-Inf` standing in for the swept thresholds, which both mislabelled its output
and parked it inside a folder `12_rdd_analysis.R` legitimately owns.)

Figures are produced only for a run's primary configuration; driver scripts set
`MAKE_PLOTS = FALSE` for the secondary cells they sweep over, so most folders
hold numbers only.

## Special folders

| Folder | What it holds |
|---|---|
| `_builds/<instrument>_w<N>/` | ERT-episode match accounting for one build (which episodes the election spine can and cannot reach, and why) |
| `_sweeps/restriction_grid_<instr>_w<N>_trt-<treatment>/` | to-do 1: marginal tables and colour-coded pair grids over the five sample-restriction axes |
| `_sweeps/window_sweep_<sample>/` | to-do 3: first stage / RD / fuzzy RD against window length, N = 1..10 |
| `_sweeps/alt_specs_<sample>/` | to-dos 4-6: instrument x treatment-definition grids, plus what the DDCG extension actually adds |
| `_sweeps/instrument_overlap/` | to-do 7: UpSet plots and Jaccard heatmaps measuring how much the instruments really differ |
| `_legacy/` | pre-run-folder output, kept for reference; see its own README for two caveats |

## Reproducing

```
Rscript --no-init-file scripts/01f_load_extra_outcomes.R   # extra outcomes
Rscript --no-init-file scripts/02a_build_panel.R           # country-year panel
Rscript --no-init-file scripts/11_build_rdd_data.R         # default build
Rscript --no-init-file scripts/12_rdd_analysis.R           # default run
Rscript --no-init-file scripts/13_restriction_grid.R       # to-do 1
Rscript --no-init-file scripts/14_window_sweep.R           # to-do 3
Rscript --no-init-file scripts/15_alt_specs.R              # to-dos 4-6
Rscript --no-init-file scripts/16_instrument_overlap.R     # to-do 7
```

`--no-init-file` is required: this machine's `~/.Rprofile` calls
`credentials::set_github_pat()`, which errors without a PAT. Do not use
`--vanilla` -- it also drops the user library, hiding most installed packages.

Every script guards its toggles with `if (!exists(...))`, so any of them can be
driven with overrides:

```r
Rscript --no-init-file -e 'ILLIBERAL_CUTOFF <- 0.6; source("scripts/14_window_sweep.R")'
```
