# ==============================================================================
# Fuzzy RDD estimation: illiberal-party elections -> democratic backsliding
# -> economic, fiscal and institutional outcomes
#
# Loads a build produced by 11_build_rdd_data.R and runs, for every outcome, on
# the sample selected by the toggles below:
#   - a first stage: rdrobust(y = <treatment>, x = running_var)
#   - a fuzzy RD:    rdrobust(y = Y_<outcome>, x = running_var, fuzzy = <treatment>)
#   - a plain reduced-form RD on Y (fuzzy = NULL) for comparison
#
# Everything this script writes goes into ONE run folder,
# output/runs/<slug>/, where the slug encodes the instrument, window,
# treatment definition and both sample-restriction thresholds. Different
# versions of the analysis therefore never overwrite each other, and each
# folder is self-describing via its run_config.csv.
#
# Data:   data/rdd_build/rdd_<ILLIBERALISM_VAR>_w<N>.rds
# Output: output/runs/<slug>/
#           run_config.csv
#           rdd_results.csv, rdd_first_stage_results.csv
#           first_stage_table.html, outcomes_table.html
#           plots/first_stage.png, plots/running_var_density.png,
#           plots/outcomes_<panel>.png, plots/outcomes_all.png
# ==============================================================================

library(tidyverse)
library(rdrobust)
library(here)
library(patchwork)
library(gt)

source(here::here("scripts", "rdd_helpers.R"))

data_dir <- here::here("data")

# ------------------------------------------------------------------------------
# Toggles
#
# All set with `if (!exists(...))` so a driver script (14_window_sweep.R,
# 15_alt_specs.R) can source() this file into an environment that already
# defines some of them. Running standalone uses the defaults.
# ------------------------------------------------------------------------------

# Which build to read: these two must match a build produced by
# 11_build_rdd_data.R (they select the file, they don't re-derive anything).
if (!exists("ILLIBERALISM_VAR")) ILLIBERALISM_VAR <- "v2xpa_antiplural"
if (!exists("BACKSLIDING_WINDOW_YEARS")) BACKSLIDING_WINDOW_YEARS <- 5
# Selects which build to read, and is echoed into the run slug. Set in
# 11_build_rdd_data.R -- see the long comment there for what it does and why
# TRUE is the default.
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- TRUE
}
# Selects the PLACEBO build -- outcomes and treatment measured over the window
# BEFORE the election rather than after it, as a pre-trend check. Set in
# 11_build_rdd_data.R; see the long comment there. Like the convention above it
# selects a build, it does not re-derive anything, and a mismatch between what
# is asked for here and what the build actually is hard-errors below.
if (!exists("PLACEBO_PRE_WINDOW")) PLACEBO_PRE_WINDOW <- FALSE
stopifnot(is.logical(PLACEBO_PRE_WINDOW))

# Which treatment the fuzzy RD instruments for:
#   backsliding_Nyr        an ERT autocratization episode starts in the window
#   backsliding_union_Nyr  that OR an Acemoglu et al. (DDCG) democratic reversal.
#                          DDCG only covers 1960-2010, so this can only add
#                          events for elections whose window overlaps that span.
#   polyarchy_decline      continuous: the FALL in V-Dem polyarchy over the
#                          window, negated so higher = more backsliding, exactly
#                          like the two binary definitions. rdrobust's fuzzy=
#                          argument is a Wald ratio either way, so no special
#                          handling is needed -- but note the LATE is then "per
#                          one unit of polyarchy decline", i.e. per a full 0-1
#                          swing of the index, not per episode.
if (!exists("TREATMENT_VAR")) TREATMENT_VAR <- "backsliding_Nyr"
stopifnot(TREATMENT_VAR %in% names(TREATMENT_LABELS))

# ---- Sample restrictions -----------------------------------------------------
#
# BOTH thresholds below accept THREE forms, and which one you mean is decided by
# what you write -- there is no separate "type" switch to keep in sync:
#
#     0.6      ABSOLUTE  a value on the variable's own scale
#     "q50"    QUANTILE  a percentile of the variable's distribution
#     -Inf     NONE      no restriction
#
# Anything else is an error, not a fallback (see parse_threshold() in
# rdd_helpers.R for why). Whichever form is used, the run prints a "Sample
# restrictions" block naming the spec, how it was read, the absolute value it
# resolved to, and how many elections it cost; the same information lands in the
# run's run_config.csv and in every table subtitle.
#
# Use the QUANTILE form whenever you are comparing across instruments. The five
# candidate instruments are not on a common scale (v2xpa_* are [0,1] indices,
# v2pariglef_neg and v2paanteli are expert scales running about -2 to +4,
# ep_galtan runs 4.5 to 9.4), so one absolute number means a different thing for
# each -- and for ep_galtan an ILLIBERAL_CUTOFF of 0.6 sits below the whole
# range and restricts nothing at all.

# Minimum standardized top-2 score gap. score_gap_z = score_gap /
# sd(illiberalism_score across every top-2 party-year in that country) -- the
# raw illiberal-vs-other gap scaled by how much that country's own parties
# typically differ in ideology. Since illiberal_score is defined as the HIGHER
# of the top-2's scores, score_gap_z is always >= 0, so an absolute 0 is a
# near-no-op; meaningful absolute thresholds start above 0. -Inf differs from 0
# in one way worth knowing: it also keeps the countries with too little data to
# compute their own score SD, where score_gap_z is NA.
#
# Restricting this way tests whether the first stage is diluted by low-contrast
# top-2 pairs: in the full sample 29% of elections have a top-2 score_gap under
# 0.05, i.e. the two parties are barely distinguishable in illiberalism, so
# crossing the vote-share cutoff there isn't a real treatment contrast.
# Compared with `>=` (a gap of exactly the threshold is kept).
if (!exists("SCORE_GAP_MIN")) SCORE_GAP_MIN <- 0

# Minimum illiberal_score -- the more-illiberal top-2 member's raw LEVEL, not
# the gap between the two. Restricting to elections where that party clears an
# ideological bar in absolute terms, rather than merely edging out the other
# top-2 member, is the restriction that most moves the first stage.
# Compared with `>` (strictly above the threshold).
if (!exists("ILLIBERAL_CUTOFF")) ILLIBERAL_CUTOFF <- 0.6

# Maximum other_score -- a CEILING on the LESS-illiberal top-2 member, the
# mirror of ILLIBERAL_CUTOFF's floor on the more-illiberal one. Same three
# forms, but the no-restriction sentinel is +Inf rather than -Inf, because
# nothing is above +Inf.
#
# On its own it is "the opponent is not itself illiberal". Paired with
# ILLIBERAL_CUTOFF at the same number it is the restriction the Sep 7 memo
# actually wanted: ONE side classified illiberal and the OTHER not. That
# addresses the standing objection that "narrow" constrains the vote margin
# while the illiberality GAP between the two parties can still be near zero --
# 29% of elections have a top-2 gap under 0.05, where crossing the cutoff is
# not a real treatment contrast at all.
#
# Where the number comes from is deliberately not this script's business.
# 20_populist_threshold.R derives one by calibrating against The PopuList, but
# a quantile ("q50") or a hand-picked value are equally valid ways to set it,
# and 13_restriction_grid.R sweeps it as axis R6.
# Compared with `<=` (a score exactly at the threshold is kept).
if (!exists("OTHER_CUTOFF_MAX")) OTHER_CUTOFF_MAX <- Inf

# Driver scripts that only need the numbers (14_window_sweep.R's secondary
# treatment definitions, 15_alt_specs.R's grid) can set this FALSE to skip the
# figures, which are the slow part of a run.
if (!exists("MAKE_PLOTS")) MAKE_PLOTS <- TRUE

# ------------------------------------------------------------------------------
# Load + restrict
# ------------------------------------------------------------------------------

build_suffix <- paste0(
  if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "" else "_exclyr",
  if (PLACEBO_PRE_WINDOW) "_pre" else ""
)
build_path <- file.path(
  data_dir,
  "rdd_build",
  sprintf(
    "rdd_%s_w%d%s.rds",
    ILLIBERALISM_VAR, BACKSLIDING_WINDOW_YEARS, build_suffix
  )
)
if (!file.exists(build_path)) {
  stop(
    "No build at ", build_path, ".\n",
    "Run 11_build_rdd_data.R with ILLIBERALISM_VAR = '", ILLIBERALISM_VAR,
    "', BACKSLIDING_WINDOW_YEARS = ", BACKSLIDING_WINDOW_YEARS,
    ", TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR = ",
    TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR,
    " and PLACEBO_PRE_WINDOW = ", PLACEBO_PRE_WINDOW, " first."
  )
}
d <- readRDS(build_path)
# Builds written before this toggle existed carry no attribute; they all used
# the old exclude-the-election-year convention, so treat a missing attribute as
# FALSE rather than assuming it matches the current default.
build_incl <- attr(d, "includes_election_year") %||% FALSE
if (!identical(build_incl, TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR)) {
  stop(
    "Build at ", build_path, " was made with ",
    "TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR = ", build_incl,
    " but this run asked for ", TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR,
    ". Rebuild it with 11_build_rdd_data.R."
  )
}
# Same guard for the placebo, and it matters more: a placebo build estimated as
# though it were the real one would report a pre-election correlation as the
# headline post-election effect, and nothing in the numbers would look wrong.
# Builds predating the attribute are all real (non-placebo) ones.
build_placebo <- attr(d, "placebo_pre_window") %||% FALSE
if (!identical(build_placebo, PLACEBO_PRE_WINDOW)) {
  stop(
    "Build at ", build_path, " was made with PLACEBO_PRE_WINDOW = ",
    build_placebo, " but this run asked for ", PLACEBO_PRE_WINDOW,
    ". Rebuild it with 11_build_rdd_data.R."
  )
}
if (!TREATMENT_VAR %in% names(d)) {
  stop(
    "Build at ", build_path, " has no column '", TREATMENT_VAR,
    "'. It predates that treatment definition -- rebuild it with ",
    "11_build_rdd_data.R."
  )
}
cat(sprintf("Loaded %s (%d elections)\n", basename(build_path), nrow(d)))

# Both specs are parsed and resolved against the FULL loaded build, before
# either filter is applied. That ordering is load-bearing: resolving a quantile
# after the other restriction had bitten would make the two axes interact, so
# changing SCORE_GAP_MIN would silently move an ILLIBERAL_CUTOFF of "q50" too.
score_gap_thr <- resolve_threshold(
  parse_threshold(SCORE_GAP_MIN, "SCORE_GAP_MIN"),
  d$score_gap_z, "SCORE_GAP_MIN"
)
illiberal_thr <- resolve_threshold(
  parse_threshold(ILLIBERAL_CUTOFF, "ILLIBERAL_CUTOFF"),
  d$illiberal_score, "ILLIBERAL_CUTOFF"
)
# none_value = Inf: this one is a ceiling, so "no restriction" is +Inf. It is
# resolved against other_score, the LESS-illiberal member -- so a quantile spec
# means a percentile of the opponents' distribution, not of the illiberal
# parties'.
other_thr <- resolve_threshold(
  parse_threshold(OTHER_CUTOFF_MAX, "OTHER_CUTOFF_MAX", none_value = Inf),
  d$other_score, "OTHER_CUTOFF_MAX"
)

# The run folder is named by the RESOLVED absolute values, not the specs, so a
# "q50" run and a hand-written run at the same resolved number correctly share a
# folder instead of duplicating. The specs as written are recorded alongside
# them in run_config.csv, so the folder name is never the only record of how the
# threshold was expressed.
cfg <- list(
  instrument = ILLIBERALISM_VAR,
  window = BACKSLIDING_WINDOW_YEARS,
  treatment = TREATMENT_VAR,
  score_gap_min = score_gap_thr$absolute,
  illiberal_cutoff = illiberal_thr$absolute,
  other_cutoff_max = other_thr$absolute,
  incl_election_year = TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR,
  placebo = PLACEBO_PRE_WINDOW,
  score_gap_spec = score_gap_thr$spec,
  score_gap_form = score_gap_thr$kind,
  illiberal_cutoff_spec = illiberal_thr$spec,
  illiberal_cutoff_form = illiberal_thr$kind,
  other_cutoff_spec = other_thr$spec,
  other_cutoff_form = other_thr$kind
)
slug <- run_slug(cfg)
out_run <- run_dir(cfg)
plots_dir <- file.path(out_run, "plots")

cat("\nSample restrictions:\n")
d <- apply_threshold(d, "score_gap_z", score_gap_thr, "SCORE_GAP_MIN", op = ">=")
d <- apply_threshold(d, "illiberal_score", illiberal_thr, "ILLIBERAL_CUTOFF", op = ">")
d <- apply_threshold(d, "other_score", other_thr, "OTHER_CUTOFF_MAX", op = "<=")
cat(sprintf("  %-16s %d elections\n\n", "FINAL SAMPLE", nrow(d)))

if (nrow(d) < RD_MIN_OBS) {
  stop(
    "Only ", nrow(d), " elections survive the sample restrictions -- too few to ",
    "estimate anything. Check that the thresholds are on the right scale for ",
    "instrument '", ILLIBERALISM_VAR, "'."
  )
}

# One human-readable sentence describing the sample, reused in every table
# subtitle so the restriction travels with the output rather than living only in
# the folder name.
restriction_label <- {
  parts <- c(
    if (score_gap_thr$kind != "none") {
      sprintf("score_gap_z >= %.4g [%s]", score_gap_thr$absolute, score_gap_thr$spec)
    },
    if (illiberal_thr$kind != "none") {
      sprintf("illiberal_score > %.4g [%s]", illiberal_thr$absolute, illiberal_thr$spec)
    },
    if (other_thr$kind != "none") {
      sprintf("other_score <= %.4g [%s]", other_thr$absolute, other_thr$spec)
    }
  )
  if (length(parts) == 0) "all scored elections" else paste(parts, collapse = ", ")
}

sample_label <- sprintf("%s, N = %d", slug, nrow(d))

# ------------------------------------------------------------------------------
# Outcomes, grouped into plot panels
#
# The results TABLE stays one row per outcome (flat and scannable); the panel
# grouping only controls the figures, where several same-unit series belong
# together on one set of axes:
#   growth       three log GDP-per-capita series (PWT / World Bank / IMF WEO),
#                all log-differences over the same window, so directly comparable
#   gini         the two SWIID Ginis, both 0-100 points
#   fiscal       debt and deficit, both % of GDP
#   institutions the V-Dem constraint/power indices, all on [0,1]
# Every panel is single-unit by construction -- no panel ever mixes two scales.
# ------------------------------------------------------------------------------

outcome_vars <- unlist(lapply(OUTCOME_PANELS, names), use.names = FALSE)
outcome_labels <- unlist(OUTCOME_PANELS, use.names = FALSE)
names(outcome_labels) <- outcome_vars

# Drop any outcome missing from this build (older builds, or a variable whose
# source data didn't cover this sample at all) rather than erroring later.
missing_outcomes <- setdiff(outcome_vars, names(d))
if (length(missing_outcomes) > 0) {
  warning(
    "Build is missing outcomes, dropping: ",
    paste(missing_outcomes, collapse = ", ")
  )
  OUTCOME_PANELS <- lapply(OUTCOME_PANELS, function(p) p[names(p) %in% names(d)])
  OUTCOME_PANELS <- OUTCOME_PANELS[lengths(OUTCOME_PANELS) > 0]
  outcome_vars <- setdiff(outcome_vars, missing_outcomes)
}

# ------------------------------------------------------------------------------
# Estimation
# ------------------------------------------------------------------------------

# Estimates the first stage once per subsample, then a reduced-form and a fuzzy
# RD for every outcome reusing that same treatment column. Returns BOTH a
# standalone first-stage row (the way a paper's first-stage table reports it,
# rather than folding the same numbers redundantly into every outcome row) and
# the per-outcome table.
run_spec <- function(data, spec_label, treatment_col = TREATMENT_VAR) {
  fs_fit <- safe_rdrobust(data[[treatment_col]], data$running_var)
  fs <- extract_rd(fs_fit)
  # The first stage regresses the treatment ON the running variable, so its own
  # y IS the treatment -- the treated count is over that same mask.
  fs_treated <- rd_n_treated(
    data[[treatment_col]], data$running_var, data[[treatment_col]]
  )

  first_stage <- tibble(
    spec = spec_label,
    treatment = treatment_col,
    n_first_stage = fs$N,
    n_first_stage_treated = fs_treated,
    first_stage_coef = fs$coef,
    first_stage_se = fs$se,
    first_stage_pval = fs$pval,
    first_stage_bandwidth = fs$bw
  )

  outcomes <- map_dfr(outcome_vars, function(oc) {
    rf <- extract_rd(safe_rdrobust(data[[oc]], data$running_var))
    late <- extract_rd(safe_rdrobust(
      data[[oc]],
      data$running_var,
      fuzzy = data[[treatment_col]]
    ))

    tibble(
      outcome = oc,
      outcome_label = unname(outcome_labels[oc]),
      panel = names(OUTCOME_PANELS)[
        vapply(OUTCOME_PANELS, function(p) oc %in% names(p), logical(1))
      ][1],
      spec = spec_label,
      treatment = treatment_col,
      n_first_stage = fs$N,
      n_first_stage_treated = fs_treated,
      first_stage_coef = fs$coef,
      first_stage_pval = fs$pval,
      n_outcome = rf$N,
      # Per outcome, not per spec: outcomes differ in missingness, so the
      # treated count in each estimation sample differs too.
      n_outcome_treated = rd_n_treated(
        data[[oc]], data$running_var, data[[treatment_col]]
      ),
      rd_estimate = rf$coef,
      rd_se = rf$se,
      rd_pval = rf$pval,
      late_estimate = late$coef,
      late_se = late$se,
      late_pval = late$pval,
      bandwidth = late$bw
    )
  })

  list(first_stage = first_stage, outcomes = outcomes)
}

# ------------------------------------------------------------------------------
# Plotting
#
# side_fit(), binned_means() and build_panel_plot() live in rdd_helpers.R: they
# depend only on SERIES_COLORS/SERIES_SHAPES/RD_MIN_OBS/safe_rdrobust, and
# 17_party_outcomes_rdd.R draws the same RD panel figure. make_rd_plots() stays
# here because it is the one part that reads this run's own state (OUTCOME_PANELS
# as narrowed above, plots_dir, sample_label, TREATMENT_DISPLAY).
# ------------------------------------------------------------------------------

make_rd_plots <- function(data, treatment_col = TREATMENT_VAR) {
  # First stage. y.lim is fixed to [0,1] for the binary treatments (it's a
  # probability, and a wide-CI bin at small N would otherwise blow the scale
  # out); the continuous treatment gets an automatic scale.
  fs_ylim <- if (treatment_col == "polyarchy_decline") NULL else c(0, 1)
  fs_label <- if (treatment_col == "polyarchy_decline") {
    "Polyarchy decline over window"
  } else {
    "P(backsliding within window)"
  }
  p_fs <- build_panel_plot(
    data,
    setNames("First stage", treatment_col),
    sprintf("First stage: %s\n%s", TREATMENT_DISPLAY[[treatment_col]], sample_label),
    fs_label,
    y_lim = fs_ylim
  )
  if (!is.null(p_fs)) {
    ggsave(file.path(plots_dir, "first_stage.png"), p_fs, width = 7, height = 3.6, dpi = 150)
    cat("Saved plots/first_stage.png\n")
  }

  # One figure per outcome panel, plus a combined sheet. Each panel keeps its
  # own independently-selected bandwidth and bins (outcomes differ in
  # missingness and variance, so a shared binning across a facet_wrap would be
  # wrong), assembled with patchwork rather than faceting.
  panels <- list()
  for (nm in names(OUTCOME_PANELS)) {
    p <- build_panel_plot(
      data,
      OUTCOME_PANELS[[nm]],
      unname(PANEL_TITLES[nm]),
      unname(PANEL_YLABS[nm])
    )
    if (is.null(p)) {
      next
    }
    panels[[nm]] <- p
    ggsave(
      file.path(plots_dir, sprintf("outcomes_%s.png", nm)),
      p, width = 6, height = 3.6, dpi = 150
    )
  }
  if (length(panels) > 0) {
    combined <- wrap_plots(panels, ncol = 3) +
      plot_annotation(
        title = sprintf("Reduced-form RD by outcome (%s)", sample_label),
        subtitle = "Shaded band = 95% CI of the local-linear fit (conventional, not bias-corrected)"
      )
    n_rows <- ceiling(length(panels) / 3)
    ggsave(
      file.path(plots_dir, "outcomes_all.png"),
      combined, width = 16, height = 3.6 * n_rows, dpi = 150, limitsize = FALSE
    )
    cat(sprintf("Saved plots/outcomes_all.png (%d panels)\n", length(panels)))
  }
}

# ------------------------------------------------------------------------------
# Diagnostic: density of the running variable
# Standard RD sanity check -- if elections near the cutoff were selected or
# sorted (the illiberal side systematically squeaking out narrow wins), the
# density would bunch or jump at 0. Visual only, not a formal manipulation test.
# ------------------------------------------------------------------------------

density_plot <- ggplot(d, aes(x = running_var)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "grey85", colour = "white") +
  geom_density(colour = SERIES_COLORS[1], linewidth = 0.8) +
  geom_vline(xintercept = 0, colour = "grey35", linetype = "dashed", linewidth = 0.4) +
  labs(
    title = sprintf("Density of the running variable (%s)", sample_label),
    subtitle = "Dashed line = RD cutoff",
    x = "Running variable (illiberal - other vote/seat share, pp)",
    y = "Density"
  ) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank())

if (MAKE_PLOTS) {
  ggsave(file.path(plots_dir, "running_var_density.png"), density_plot, width = 8, height = 3, dpi = 150)
  cat("Saved plots/running_var_density.png\n")
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

cat(sprintf("\n=== %s ===\n", slug))
results <- run_spec(d, slug)

# The placebo has to announce itself on the face of every table. A reader who
# picks up an outcomes_table.html without the folder name in front of them must
# not be able to mistake a pre-election correlation for the headline effect.
window_label <- if (PLACEBO_PRE_WINDOW) {
  sprintf(
    "PLACEBO (pre-election): [election_year - %d, election_year - 1]",
    BACKSLIDING_WINDOW_YEARS - (!TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR)
  )
} else {
  sprintf(
    "[%s, election_year + %d]",
    if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) {
      "election_year"
    } else {
      "election_year + 1"
    },
    BACKSLIDING_WINDOW_YEARS
  )
}

subtitle <- sprintf(
  "Instrument: %s | Treatment: %s | Window: %s | Sample: %s (N = %d)",
  INSTRUMENT_DISPLAY[[ILLIBERALISM_VAR]],
  TREATMENT_DISPLAY[[TREATMENT_VAR]],
  window_label,
  restriction_label,
  nrow(d)
)

save_table_html(
  results$first_stage |>
    transmute(
      Spec = spec,
      N = n_first_stage,
      `N treated` = n_first_stage_treated,
      `First stage` = fmt_est(first_stage_coef, first_stage_se, first_stage_pval),
      Bandwidth = round(first_stage_bandwidth, 2)
    ),
  file.path(out_run, "first_stage_table.html"),
  "First stage",
  subtitle,
  note = paste(SIG_FOOTNOTE, TREATED_FOOTNOTE)
)

save_table_html(
  results$outcomes |>
    transmute(
      Panel = panel,
      Outcome = outcome_label,
      Variable = outcome,
      N = n_outcome,
      `N treated` = n_outcome_treated,
      `Reduced form` = fmt_est(rd_estimate, rd_se, rd_pval),
      `Fuzzy RD (LATE)` = fmt_est(late_estimate, late_se, late_pval)
    ),
  file.path(out_run, "outcomes_table.html"),
  "Reduced form / fuzzy RD by outcome",
  subtitle,
  note = paste(SIG_FOOTNOTE, TREATED_FOOTNOTE)
)

if (MAKE_PLOTS) {
  make_rd_plots(d)
}

write_csv(results$outcomes, file.path(out_run, "rdd_results.csv"))
write_csv(results$first_stage, file.path(out_run, "rdd_first_stage_results.csv"))

append_run_manifest(cfg, list(
  n_elections = nrow(d),
  n_treated = if (TREATMENT_VAR == "polyarchy_decline") NA else sum(d[[TREATMENT_VAR]], na.rm = TRUE),
  first_stage_coef = round(results$first_stage$first_stage_coef, 4),
  first_stage_pval = round(results$first_stage$first_stage_pval, 4)
))

n_failed <- sum(is.na(results$outcomes$late_estimate))
cat(sprintf(
  "\n%d / %d outcomes produced a fuzzy RD estimate (%d skipped -- too few observations near the cutoff, or rdrobust's internal checks failed)\n",
  nrow(results$outcomes) - n_failed, nrow(results$outcomes), n_failed
))
cat(sprintf(
  "First stage: %s  (N = %d, bw = %.2f)\n",
  fmt_est(
    results$first_stage$first_stage_coef,
    results$first_stage$first_stage_se,
    results$first_stage$first_stage_pval
  ),
  results$first_stage$n_first_stage,
  results$first_stage$first_stage_bandwidth
))
message("All output written to ", out_run)
