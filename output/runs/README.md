# Analysis runs

One folder per *version* of the fuzzy-RDD analysis. A version is fully
identified by the choices the folder name encodes:

```
instr-<instrument>_w<window>_trt-<treatment>_gap<score_gap_min>_illib<illiberal_cutoff>[_opp<other_cutoff_max>][_exclyr][_pre]
```

The last three parts appear only when they are not at their default, so every
folder written before they existed still has exactly the name it had.

| Part | Meaning |
|---|---|
| `instr-` | which V-Party score decides who the "illiberal" side of the top 2 is: `antiplural`, `popul`, `econleft` (negated `v2pariglef`, so economic LEFT is the high end), `anteli`, `galtan` |
| `w` | post-election window N, in years. Treatment *and* every outcome are measured over `(election_year - 1, election_year + N]` |
| `trt-` | `ert` (ERT autocratization episode), `ertOrDdcg` (that OR an Acemoglu et al. reversal), `polyarchy` (continuous decline in V-Dem polyarchy) |
| `gap` | minimum `score_gap_z`; `any` = no restriction |
| `illib` | minimum `illiberal_score` (the MORE illiberal of the top 2); `any` = no restriction |
| `_opp` | maximum `other_score` (the LESS illiberal of the top 2). Omitted entirely at its default of `Inf`. Paired with `illib` at the same number, this is the "one side illiberal, the other not" restriction |
| `_exclyr` | present only when `TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR = FALSE` (see below) |
| `_pre` | present only when `PLACEBO_PRE_WINDOW = TRUE`: the pre-election placebo, outcomes and treatment measured over the window BEFORE the election (see below) |

## Sample-restriction thresholds

`SCORE_GAP_MIN`, `ILLIBERAL_CUTOFF` and `OTHER_CUTOFF_MAX` each accept three
forms, and which one you mean is decided by what you write — there is no
separate "type" switch:

| Written as | Read as | Example |
|---|---|---|
| a number | absolute value on the variable's own scale | `0.6` |
| `"qNN"` | percentile of the variable's distribution | `"q50"`, `"q97.5"` |
| `-Inf` | no restriction, for a FLOOR (`SCORE_GAP_MIN`, `ILLIBERAL_CUTOFF`) | |
| `Inf` | no restriction, for a CEILING (`OTHER_CUTOFF_MAX`) | |
| a named external threshold | value read from a file another script wrote | `"popucut"`, `"popucut_pct"` |

The external forms are registered in `EXTERNAL_THRESHOLDS` in
`rdd_helpers.R`. Today there are two, both resolving out of
`data/populist_threshold.rds`, written by `01g_populist_threshold.R`:

| Spec | Resolves to | Meaning |
|---|---|---|
| `"popucut"` | 0.6535 | the PopuList-calibrated illiberality cut, accuracy criterion, carried across as a raw value |
| `"popucut_pct"` | 0.8220 | the same cut, transferred at its percentile (83.2nd) instead |

Setting `ILLIBERAL_CUTOFF` and `OTHER_CUTOFF_MAX` both to `"popucut"` is the
"one top-2 party illiberal, the other not" sample — 411 elections at w5. It is
a spec rather than a typed number so that re-running the calibration moves
every run that uses it; a hardcoded `0.6535` would not. Missing file is a hard
error naming the script to run.

Folder names still use the RESOLVED value (`..._illib0p6535_opp0p6535`), and
`run_config.csv` records the spec as written, so the provenance travels with
the run without the slug depending on a file.

The no-restriction sentinel differs by direction because nothing is below
`-Inf` and nothing is above `Inf`. `parse_threshold()` takes the right one via
its `none_value` argument, and `apply_threshold()` errors if a caller pairs a
ceiling operator with a floor sentinel — that mistake would drop every row
while the run still described itself as unrestricted.

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

## The pre-election placebo

`PLACEBO_PRE_WINDOW = TRUE` runs the whole design backwards, as a pre-trend
check: if growth was already falling *before* the narrow illiberal victory, the
headline reduced form is a trend rather than an effect.

|  | outcome window | treatment window (`_exclyr`) |
|---|---|---|
| `FALSE` (default) | years `ey` … `ey + N` | `ey + 1` … `ey + N` |
| `TRUE` | years `ey - N` … `ey - 1` | `ey - N + 1` … `ey - 1` |

It is a mirror, not a shift. The placebo outcome window *ends* where the real
one begins differencing (`ey - 1`), so the two abut with no overlap and no
year counted twice, and the treatment sits inside its own outcome window under
the same rule as the real design. Deriving the placebo by reflecting the window
*length* instead would put its treatment window a year later than its own
outcome window — testing a slightly different design from the one it is a
placebo for.

Nothing about the running variable, the cutoff or the sample changes; only
which years the outcomes and treatment are measured over. Because it changes
the data, it is part of the build identity: builds get a `_pre` suffix, carry a
`placebo_pre_window` attribute, and `12_rdd_analysis.R` refuses to estimate a
placebo build as a real one (a silent mix-up would report a pre-election
correlation as the headline effect, and nothing in the numbers would look
wrong). Every table subtitle says `PLACEBO (pre-election)` on its face.

At w = 5, `_exclyr`, full sample, the placebo is clean: GDP growth is
`+0.017 (0.025)` against the real `-0.053 (0.027)**`, and disposable Gini is
`-0.044 (0.314)` against the real `+0.663 (0.314)**`. Both pre-window estimates
are insignificant and of the opposite sign.

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

### Outcomes

37 outcomes in 15 panels, registered once in `rdd_helpers.R`
(`OUTCOME_PANELS`, `PANEL_TITLES`, `PANEL_YLABS`, `OUTCOME_FAMILIES`). A panel
is a set of series sharing one y axis, so a panel never mixes two units; a
family is a set a reader wants to look at together, and drives one window-sweep
sheet each.

Alongside the economic, inequality, fiscal and executive-power outcomes there
are 21 V-Dem democracy indices, listed in `scripts/vdem_indices.R` and grouped
on V-Dem's own taxonomy:

| Panel | Series |
|---|---|
| `vdem_high` | polyarchy, libdem, partipdem, delibdem, egaldem |
| `vdem_electoral` | elected officials, clean elections, association, suffrage, expression |
| `vdem_liberal` | rule of law (judicial and legislative constraints are in `institutions`) |
| `vdem_particip` | civil society, direct democracy, local and regional elections |
| `vdem_egal_delib` | deliberation, equal protection, equal access, equal distribution |

These run the OPPOSITE way from the economic outcomes: higher = more
democratic, so a **negative** RD estimate is the backsliding sign.

Adding an outcome means three edits, in order: get the source column into
`combined_panel.rds` (`01*` then `02a`), difference it in
`11_build_rdd_data.R`'s `build_outcomes()`, and register it here. Builds must
then be regenerated — `12_rdd_analysis.R` warns and drops outcomes a build
predates rather than failing, so a stale build silently under-reports.

One row to read with care: in a `trt-polyarchy` run, `Y_polyarchy` is the
negation of the treatment, so its fuzzy LATE is a tautology and its reduced
form *is* the first stage.

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
| `_sweeps/party_outcome_rdd_<instr>_w<N>/` | 1a: the same RD with the winner's OTHER party scores as the outcome — is a narrow anti-pluralist victory also a populist / left / minority-hostile one? |
| `_sweeps/vparty_jaccard_panels/` | 1c: the anti-pluralism x populism Jaccard heatmap over the full V-Party dataset, as a 2 x 5 OECD-by-decade grid |
| `_sweeps/vparty_ideology_quadrants/` | 1d: where the most common Wikipedia/Wikidata ideology tags sit on the illiberalism x populism plane, same 2 x 5 grid |
| `_sweeps/populist_threshold/` | 2b: the PopuList-calibrated cutoff for "illiberal", and the ROC it comes from |
| `_sweeps/cell_rdd_<instr>/` | reduced-form RD on the decade x OECD grid, every outcome, w1-10, full and PopuList-restricted samples |
| `_sweeps/econleft_split_rdd_<instr><restriction>/` | the main RD run separately on elections where the anti-pluralist party is the more RIGHT-wing of the top 2 and where it is the more LEFT-wing -- the sign-discordant test of whether the growth result is really about the economic right. `<restriction>` carries the sample the split is taken within, spelled as in a run slug and omitted axis by axis at its no-op value: no suffix = all 1,347 scored elections (822 right / 517 left), `_illib0p6535_opp0p6535` = the 411 where one top-2 party is illiberal and the other is not (241 / 170) |

The pre-run-folder output that used to sit in `_legacy/` has been deleted. It is
recoverable from the commit that preceded the cleanup, and the numbers in it
reproduce directly by setting `TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR = FALSE`.

## Reproducing

```
Rscript --no-init-file scripts/01f_load_extra_outcomes.R   # extra outcomes
Rscript --no-init-file scripts/01g_populist_threshold.R   # PopuList cut (downloads)
Rscript --no-init-file scripts/02a_build_panel.R           # country-year panel
Rscript --no-init-file scripts/11_build_rdd_data.R         # default build
Rscript --no-init-file scripts/12_rdd_analysis.R           # default run
Rscript --no-init-file scripts/13_restriction_grid.R       # to-do 1
Rscript --no-init-file scripts/14_window_sweep.R           # to-do 3
Rscript --no-init-file scripts/15_alt_specs.R              # to-dos 4-6
Rscript --no-init-file scripts/16_instrument_overlap.R     # to-do 7
Rscript --no-init-file scripts/17_party_outcomes_rdd.R     # 1a
Rscript --no-init-file scripts/18_vparty_jaccard_panels.R  # 1c
Rscript --no-init-file scripts/19_vparty_ideology_quadrants.R  # 1d
Rscript --no-init-file scripts/21_cell_rdd.R               # decade x OECD cells
Rscript --no-init-file scripts/22_econleft_split_rdd.R     # econ L-R sign split
SPLIT_SAMPLE=popucut \
  Rscript --no-init-file scripts/22_econleft_split_rdd.R   # ... within the PopuList 411

```

Scripts 18 and 19 read raw V-Party and never touch a build, so they source
`scripts/vparty_helpers.R` rather than `scripts/rdd_helpers.R`. Each plots a
LIST of score pairs and suffixes its outputs with a pair slug
(`antiplural_vs_popul`, `econlr_vs_antiplural`); the axes, the Jaccard bin
mode and whether `coord_fixed()` applies are all derived from `VPARTY_SCORES`
in that helper. The left-right figures use the RAW right-positive
`v2pariglef`, not the negated `v2pariglef_neg` the RDD carries, so the axis
reads left-to-right conventionally -- the axis title says so.
`01g_populist_threshold.R` is a calibration step and sources neither; it sits
in the loader tier because `12_rdd_analysis.R` now consumes its output through
the `"popucut"` spec, so it has to run before the build.

`--no-init-file` is required: this machine's `~/.Rprofile` calls
`credentials::set_github_pat()`, which errors without a PAT. Do not use
`--vanilla` -- it also drops the user library, hiding most installed packages.

Every script guards its toggles with `if (!exists(...))`, so any of them can be
driven with overrides:

```r
Rscript --no-init-file -e 'ILLIBERAL_CUTOFF <- 0.6; source("scripts/14_window_sweep.R")'
```
