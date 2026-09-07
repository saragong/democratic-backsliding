# ==============================================================================
# Window-length sweep: how do the first stage, the reduced-form RD and the
# fuzzy RD move as the post-election window N runs from 1 to 10 years?
#
# Both the treatment and every outcome are anchored to the SAME
# (election_year - 1) -> (election_year + N) span, so changing N moves the whole
# design consistently rather than re-timing the treatment against a fixed
# outcome window. That alignment is enforced in 11_build_rdd_data.R; this script
# just drives it.
#
# For each N it (a) builds data/rdd_build/rdd_<instr>_w<N>.rds if absent, then
# (b) runs 12_rdd_analysis.R once per treatment definition, each into its own
# run folder. Finally it pools every run's CSVs into one comparison set.
#
# Output: output/runs/_sweeps/window_sweep/
#           window_sweep_results.csv        one row per outcome x window x treatment
#           window_sweep_first_stage.csv    one row per window x treatment
#           window_first_stage.png/.html
#           window_rdd_by_outcome.png, window_fuzzy_by_outcome.png
# ==============================================================================

library(tidyverse)
library(here)
library(gt)

source(here::here("scripts", "rdd_helpers.R"))

data_dir <- here::here("data")
scripts_dir <- here::here("scripts")

if (!exists("WINDOWS")) WINDOWS <- 1:10
if (!exists("SWEEP_INSTRUMENT")) SWEEP_INSTRUMENT <- "v2xpa_antiplural"
if (!exists("SWEEP_TREATMENTS")) {
  SWEEP_TREATMENTS <- c(
    "backsliding_Nyr", "backsliding_union_Nyr", "polyarchy_decline"
  )
}
# The sample restrictions held fixed across the sweep. -Inf on both keeps the
# full sample, so the only thing varying is the window.
if (!exists("SWEEP_SCORE_GAP_MIN")) SWEEP_SCORE_GAP_MIN <- -Inf
if (!exists("SWEEP_ILLIBERAL_CUTOFF")) SWEEP_ILLIBERAL_CUTOFF <- -Inf

# Run one of the pipeline scripts in a fresh environment with the given toggle
# overrides pre-defined. The scripts all guard their toggles with
# `if (!exists(...))`, so anything set here wins and everything else falls back
# to that script's own default. A fresh env per call also means no state leaks
# between windows.
run_script_with <- function(script, overrides) {
  env <- new.env(parent = globalenv())
  for (nm in names(overrides)) assign(nm, overrides[[nm]], envir = env)
  sys.source(file.path(scripts_dir, script), envir = env)
  invisible(env)
}

# ------------------------------------------------------------------------------
# Build + estimate
# ------------------------------------------------------------------------------

for (n in WINDOWS) {
  build_path <- file.path(
    data_dir, "rdd_build", sprintf("rdd_%s_w%d.rds", SWEEP_INSTRUMENT, n)
  )
  if (file.exists(build_path)) {
    cat(sprintf("[w=%d] build cached, skipping\n", n))
  } else {
    cat(sprintf("[w=%d] building...\n", n))
    run_script_with("11_build_rdd_data.R", list(
      ILLIBERALISM_VAR = SWEEP_INSTRUMENT,
      BACKSLIDING_WINDOW_YEARS = n
    ))
  }

  for (trt in SWEEP_TREATMENTS) {
    cat(sprintf("[w=%d] estimating %s...\n", n, trt))
    run_script_with("12_rdd_analysis.R", list(
      ILLIBERALISM_VAR = SWEEP_INSTRUMENT,
      BACKSLIDING_WINDOW_YEARS = n,
      TREATMENT_VAR = trt,
      SCORE_GAP_MIN = SWEEP_SCORE_GAP_MIN,
      ILLIBERAL_CUTOFF = SWEEP_ILLIBERAL_CUTOFF,
      # Only the main treatment definition gets the full figure set; the other
      # two contribute numbers to the sweep plots below.
      MAKE_PLOTS = identical(trt, SWEEP_TREATMENTS[1])
    ))
  }
}

# ------------------------------------------------------------------------------
# Pool every run's results
# ------------------------------------------------------------------------------

# Name the sweep folder after the sample restriction it holds fixed, so the
# unrestricted sweep and a restricted one (where the first stage actually has
# any signal) don't overwrite each other.
sweep_name <- paste0(
  "window_sweep_gap", fmt_slug_num(SWEEP_SCORE_GAP_MIN),
  "_illib", fmt_slug_num(SWEEP_ILLIBERAL_CUTOFF)
)
out_dir <- sweep_dir(sweep_name)

collect <- function(file_name) {
  map_dfr(WINDOWS, function(n) {
    map_dfr(SWEEP_TREATMENTS, function(trt) {
      cfg <- list(
        instrument = SWEEP_INSTRUMENT, window = n, treatment = trt,
        score_gap_min = SWEEP_SCORE_GAP_MIN,
        illiberal_cutoff = SWEEP_ILLIBERAL_CUTOFF
      )
      path <- file.path(RUNS_ROOT, run_slug(cfg), file_name)
      if (!file.exists(path)) {
        return(NULL)
      }
      read_csv(path, show_col_types = FALSE) |>
        mutate(window = n, treatment_def = trt, .before = 1)
    })
  })
}

results <- collect("rdd_results.csv")
first_stage <- collect("rdd_first_stage_results.csv")

write_csv(results, file.path(out_dir, "window_sweep_results.csv"))
write_csv(first_stage, file.path(out_dir, "window_sweep_first_stage.csv"))

# c() drops NULLs; paste(NULL, "x", sep = ", ") would leave a leading comma.
restriction_parts <- c(
  if (is.finite(SWEEP_SCORE_GAP_MIN)) {
    sprintf("score_gap_z >= %s", SWEEP_SCORE_GAP_MIN)
  },
  if (is.finite(SWEEP_ILLIBERAL_CUTOFF)) {
    sprintf("illiberal_score > %s", SWEEP_ILLIBERAL_CUTOFF)
  }
)
restriction_label <- if (length(restriction_parts) == 0) {
  "all elections"
} else {
  paste(restriction_parts, collapse = ", ")
}

trt_display <- unname(TREATMENT_DISPLAY[SWEEP_TREATMENTS])
label_treatment <- function(x) factor(TREATMENT_DISPLAY[x], levels = trt_display)

# Okabe-Ito, same fixed order as 12_rdd_analysis.R so a treatment definition
# keeps its colour across every figure in the project.
TRT_COLORS <- setNames(c("#0072B2", "#D55E00", "#009E73")[seq_along(SWEEP_TREATMENTS)], trt_display)
TRT_SHAPES <- setNames(c(16, 17, 15)[seq_along(SWEEP_TREATMENTS)], trt_display)

# ---- first stage by window ---------------------------------------------------
# polyarchy_decline is a CONTINUOUS treatment, so its first-stage coefficient is
# in index units (a fall in polyarchy) while the two binary definitions' are in
# probability units. Plotting them on one axis would be a dual-scale chart in
# disguise, so they get separate facets with free y scales.
fs_plot_data <- first_stage |>
  mutate(
    treatment_label = label_treatment(treatment_def),
    lo = first_stage_coef - 1.96 * first_stage_se,
    hi = first_stage_coef + 1.96 * first_stage_se
  )

p_fs <- ggplot(fs_plot_data, aes(x = window, y = first_stage_coef, colour = treatment_label)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_ribbon(
    aes(ymin = lo, ymax = hi, fill = treatment_label),
    alpha = 0.13, colour = NA
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(shape = treatment_label), size = 1.9) +
  facet_wrap(~treatment_label, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = TRT_COLORS) +
  scale_fill_manual(values = TRT_COLORS) +
  scale_shape_manual(values = TRT_SHAPES) +
  scale_x_continuous(breaks = WINDOWS) +
  labs(
    title = "First stage by post-election window length",
    subtitle = sprintf(
      "Instrument: %s | Sample: %s. Band = 95%% CI. Treatment and every outcome share the same (t-1, t+N] window.",
      INSTRUMENT_DISPLAY[[SWEEP_INSTRUMENT]], restriction_label
    ),
    x = "Window length N (years after the election)",
    y = "RD estimate at the cutoff",
    colour = NULL, fill = NULL, shape = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 7.5)
  )

ggsave(file.path(out_dir, "window_first_stage.png"), p_fs, width = 11, height = 3.4, dpi = 150)
cat("Saved window_first_stage.png\n")

save_table_html(
  first_stage |>
    transmute(
      Window = window,
      Treatment = unname(TREATMENT_DISPLAY[treatment_def]),
      N = n_first_stage,
      `First stage` = fmt_est(first_stage_coef, first_stage_se, first_stage_pval),
      Bandwidth = round(first_stage_bandwidth, 2)
    ) |>
    arrange(Treatment, Window),
  file.path(out_dir, "window_first_stage.html"),
  "First stage by window length",
  sprintf(
    "Instrument: %s | Sample: %s",
    INSTRUMENT_DISPLAY[[SWEEP_INSTRUMENT]], restriction_label
  )
)

# ---- outcomes by window ------------------------------------------------------
# A fuzzy LATE is a Wald ratio (reduced form / first stage), so wherever the
# first stage is near zero the estimate and its CI blow up by several orders of
# magnitude -- a handful of such points otherwise squash every readable facet
# into a flat line at 0. Clip each facet's y range to the span of its POINT
# estimates (padded), and say so on the figure; the unclipped numbers are all in
# window_sweep_results.csv. Point estimates are never moved, only the ribbon is
# cut off at the panel edge, so nothing is silently rescaled.
clip_to_point_range <- function(dat, pad = 0.6) {
  dat |>
    group_by(outcome) |>
    mutate(
      .rng_lo = min(est, na.rm = TRUE),
      .rng_hi = max(est, na.rm = TRUE),
      .pad = pmax((.rng_hi - .rng_lo) * pad, abs(.rng_hi) * 0.1, 1e-9),
      lo = pmax(lo, .rng_lo - .pad),
      hi = pmin(hi, .rng_hi + .pad)
    ) |>
    ungroup() |>
    select(-.rng_lo, -.rng_hi, -.pad)
}

# by_treatment = FALSE for the reduced form: Y on running_var with no fuzzy
# argument doesn't reference the treatment at all, so all three treatment
# definitions return the IDENTICAL estimate and plotting them as three series
# just overplots one line three times, implying a comparison that isn't there.
outcome_sweep_plot <- function(
  est_col, se_col, pval_col, title, file_name,
  clip = FALSE, by_treatment = TRUE
) {
  results_used <- if (by_treatment) {
    results
  } else {
    results |> filter(treatment_def == SWEEP_TREATMENTS[1])
  }
  dat <- results_used |>
    mutate(
      treatment_label = label_treatment(treatment_def),
      # Facets are ordered by the canonical panel order, not alphabetically, and
      # labelled with the panel title so a series name like "PWT" or
      # "Disposable" is intelligible outside its own panel.
      outcome_facet = factor(
        outcome_full_label(outcome),
        levels = unique(outcome_full_label(ALL_OUTCOME_VARS))
      ),
      est = .data[[est_col]],
      lo = .data[[est_col]] - 1.96 * .data[[se_col]],
      hi = .data[[est_col]] + 1.96 * .data[[se_col]]
    ) |>
    filter(!is.na(est))
  if (nrow(dat) == 0) {
    return(invisible(NULL))
  }
  if (clip) dat <- clip_to_point_range(dat)
  if (!by_treatment) {
    dat$treatment_label <- factor("Reduced form", levels = "Reduced form")
  }
  series_colors <- if (by_treatment) TRT_COLORS else c("Reduced form" = SERIES_COLORS[1])
  series_shapes <- if (by_treatment) TRT_SHAPES else c("Reduced form" = SERIES_SHAPES[1])

  p <- ggplot(dat, aes(x = window, y = est, colour = treatment_label)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.35) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = treatment_label), alpha = 0.10, colour = NA) +
    geom_line(linewidth = 0.6) +
    geom_point(aes(shape = treatment_label), size = 1.2) +
    facet_wrap(~outcome_facet, scales = "free_y", ncol = 4, labeller = label_wrap_gen(38)) +
    scale_colour_manual(values = series_colors) +
    scale_fill_manual(values = series_colors) +
    scale_shape_manual(values = series_shapes) +
    scale_x_continuous(breaks = WINDOWS) +
    labs(
      title = title,
      subtitle = sprintf(
        "Instrument: %s | Sample: %s. Band = 95%% CI, free y scale per outcome.%s%s",
        INSTRUMENT_DISPLAY[[SWEEP_INSTRUMENT]], restriction_label,
        if (clip) " CIs clipped to the span of the point estimates; full values in window_sweep_results.csv." else "",
        if (!by_treatment) " The reduced form does not depend on the treatment definition, so there is one line." else ""
      ),
      x = "Window length N (years after the election)",
      y = "Estimate at the cutoff",
      colour = NULL, fill = NULL, shape = NULL
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = if (by_treatment) "bottom" else "none",
      strip.text = element_text(size = 7)
    )
  n_rows <- ceiling(n_distinct(dat$outcome_facet) / 4)
  ggsave(file.path(out_dir, file_name), p, width = 12, height = 2.1 * n_rows + 1.2, dpi = 150, limitsize = FALSE)
  cat(sprintf("Saved %s\n", file_name))
}

outcome_sweep_plot(
  "rd_estimate", "rd_se", "rd_pval",
  "Reduced-form RD by outcome and window length",
  "window_rdd_by_outcome.png",
  by_treatment = FALSE
)
outcome_sweep_plot(
  "late_estimate", "late_se", "late_pval",
  "Fuzzy RD (LATE) by outcome and window length",
  "window_fuzzy_by_outcome.png",
  clip = TRUE
)

message("Window sweep written to ", out_dir)
