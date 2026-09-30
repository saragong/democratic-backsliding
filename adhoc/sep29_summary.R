# ==============================================================================
# One table answering the Sep 29 asks: the headline RD (narrow anti-pluralist
# victory -> PWT GDP growth since the pre-election year) under every new
# sample, restriction and adjustment, next to the baseline.
#
# Reads the pooled window-sweep output of each spec (14_window_sweep.R) and the
# left/right split (22_econleft_split_rdd.R); estimates nothing. Each spec's
# folder is computed with spec_slug() from the resolved cuts, so this reads
# exactly the folders the sweeps wrote.
#
# Output: output/adhoc/sep29_summary.csv / .html
# ==============================================================================

library(tidyverse)
library(here)
library(gt)

source(here::here("scripts", "rdd_helpers.R"))

if (!exists("SUM_WINDOW")) SUM_WINDOW <- 5
if (!exists("SUM_OUTCOME")) SUM_OUTCOME <- "Y_gdp_growth"
if (!exists("SUM_TREATMENT")) SUM_TREATMENT <- DEFAULT_TREATMENT

thr <- readRDS(here::here("data", "populist_threshold.rds"))
spec <- function(...) {
  base <- list(
    instrument = DEFAULT_INSTRUMENT, score_gap_min = -Inf,
    illiberal_cutoff = -Inf, other_cutoff_max = Inf,
    incl_election_year = DEFAULT_INCL_ELECTION_YEAR, placebo = FALSE
  )
  modifyList(base, list(...))
}
pair <- function(cut, scale = "raw") {
  list(illiberal_cutoff = cut, other_cutoff_max = cut, threshold_scale = scale)
}

SPECS <- list(
  list(group = "Baseline", label = "All scored elections", cfg = spec()),
  list(group = "Baseline", label = "PopuList cut 0.65, one side above, one below", cfg = do.call(spec, pair(thr$threshold_abs))),
  list(group = "Thresholds", label = sprintf("Within-country percentile cut, accuracy (%.1fth pct)", 100 * thr$threshold_pct_ctry),
       cfg = do.call(spec, pair(thr$threshold_pct_ctry, "ctry_pct"))),
  list(group = "Thresholds", label = sprintf("Within-country percentile cut, Youden (%.1fth pct)", 100 * thr$threshold_pct_ctry_youden),
       cfg = do.call(spec, pair(thr$threshold_pct_ctry_youden, "ctry_pct"))),
  list(group = "Pooled comparisons for the split", label = "Left-right gap >= 1", cfg = spec(lr_gap_min = 1)),
  list(group = "Pooled comparisons for the split", label = "Parties on opposite sides of the left-right centre", cfg = spec(lr_straddle = TRUE)),
  list(group = "Placebo: left-right gap <= 1", label = "All scored elections", cfg = spec(lr_gap_max = 1)),
  list(group = "Placebo: left-right gap <= 1", label = "PopuList cut 0.65", cfg = do.call(spec, c(pair(thr$threshold_abs), lr_gap_max = 1))),
  list(group = "Placebo: left-right gap <= 1", label = "Within-country cut, accuracy",
       cfg = do.call(spec, c(pair(thr$threshold_pct_ctry, "ctry_pct"), lr_gap_max = 1))),
  list(group = "Placebo: left-right gap <= 1", label = "Within-country cut, Youden",
       cfg = do.call(spec, c(pair(thr$threshold_pct_ctry_youden, "ctry_pct"), lr_gap_max = 1))),
  list(group = "Funke populist spells", label = "Excluding windows that overlap a spell", cfg = spec(exclude_funke = TRUE)),
  list(group = "Covariate-adjusted (local projection)", label = "All scored elections", cfg = spec(rd_covariates = "lp")),
  list(group = "Covariate-adjusted (local projection)", label = "PopuList cut 0.65", cfg = do.call(spec, c(pair(thr$threshold_abs), rd_covariates = "lp")))
)

read_pooled <- function(s) {
  dir <- file.path(spec_path(s$cfg), "pooled")
  f <- file.path(dir, "window_sweep_results.csv")
  if (!file.exists(f)) {
    return(tibble(group = s$group, label = s$label, spec = spec_slug(s$cfg), note = "not run"))
  }
  r <- read_csv(f, show_col_types = FALSE) |>
    filter(window == SUM_WINDOW, treatment_def == SUM_TREATMENT, outcome == SUM_OUTCOME)
  fsr <- read_csv(file.path(dir, "window_sweep_first_stage.csv"), show_col_types = FALSE) |>
    filter(window == SUM_WINDOW, treatment_def == SUM_TREATMENT)
  tibble(
    group = s$group, label = s$label, spec = spec_slug(s$cfg),
    n = r$n_outcome, n_treated = r$n_outcome_treated,
    fs = fsr$first_stage_coef, fs_se = fsr$first_stage_se, fs_p = fsr$first_stage_pval,
    rf = r$rd_estimate, rf_se = r$rd_se, rf_p = r$rd_pval,
    late = r$late_estimate, late_se = r$late_se, late_p = r$late_pval,
    note = NA_character_
  )
}
tab <- map_dfr(SPECS, read_pooled)

# The left/right split: reduced form within each subset.
split_rows <- function(cfg, label) {
  f <- file.path(spec_path(cfg), "econleft_split", "econleft_split_results.csv")
  if (!file.exists(f)) {
    return(NULL)
  }
  read_csv(f, show_col_types = FALSE) |>
    filter(window == SUM_WINDOW, outcome == SUM_OUTCOME, treatment == SUM_TREATMENT) |>
    transmute(
      group = "Left/right split (anti-pluralist is the more ... of the two)",
      label = sprintf("%s: %s", label, if_else(subset == "right", "more RIGHT-wing", "more LEFT-wing")),
      spec = spec_slug(cfg), n = n_outcome, n_treated = n_outcome_treated,
      fs = first_stage_coef, fs_se = first_stage_se, fs_p = first_stage_pval,
      rf = rd_estimate, rf_se = rd_se, rf_p = rd_pval,
      late = late_estimate, late_se = late_se, late_p = late_pval,
      note = NA_character_
    )
}
tab <- bind_rows(
  tab,
  split_rows(spec(), "All, no gap floor"),
  split_rows(spec(lr_gap_min = 1), "All, left-right gap >= 1"),
  split_rows(spec(lr_straddle = TRUE), "Opposite sides of the left-right centre")
)

out_csv <- adhoc_path("sep29_summary.csv")
write_csv(tab, out_csv)

g <- tab |>
  mutate(
    `First stage` = fmt_est(fs, fs_se, fs_p),
    `Reduced form` = fmt_est(rf, rf_se, rf_p),
    `Fuzzy LATE` = fmt_est(late, late_se, late_p)
  ) |>
  select(group, Sample = label, N = n, Treated = n_treated, `First stage`, `Reduced form`, `Fuzzy LATE`) |>
  gt(groupname_col = "group") |>
  sub_missing(missing_text = "--") |>
  tab_header(
    title = "The Sep 29 asks: the headline RD under each new sample and adjustment",
    subtitle = sprintf(
      "Outcome: %s. Treatment: %s. w%d. Robust bias-corrected rdrobust estimates, MSE-optimal bandwidth.",
      outcome_full_label(SUM_OUTCOME), TREATMENT_DISPLAY[[SUM_TREATMENT]], SUM_WINDOW
    )
  ) |>
  tab_source_note(paste(
    SIG_FOOTNOTE,
    "Left-right gap: |difference in V-Party v2pariglef| between the top-2 parties.",
    "Funke exclusion conditions on the post-election window (see add_funke_overlap()).",
    "Covariate adjustment: each outcome's leave-country-out local projection on pre-election history, entered linearly (Calonico et al. 2019).",
    "The split's first stage is estimated within each subset."
  )) |>
  apply_table_style(row_group_borders = TRUE)
gtsave(g, adhoc_path("sep29_summary.html"))
print(tab |> select(group, label, n, rf, rf_se, rf_p, late, note), n = Inf, width = Inf)
message("Wrote ", out_csv)
