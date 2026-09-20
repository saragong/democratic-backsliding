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

if (!exists("WINDOWS")) {
  WINDOWS <- 1:10
}
if (!exists("SWEEP_INSTRUMENT")) {
  SWEEP_INSTRUMENT <- "v2xpa_antiplural"
}
if (!exists("SWEEP_TREATMENTS")) {
  SWEEP_TREATMENTS <- c(
    "backsliding_Nyr",
    "backsliding_union_Nyr",
    "polyarchy_decline"
  )
}
# The sample restrictions held fixed across the sweep. -Inf on both keeps the
# full sample, so the only thing varying is the window.
if (!exists("SWEEP_SCORE_GAP_MIN")) {
  SWEEP_SCORE_GAP_MIN <- -Inf
}
if (!exists("SWEEP_ILLIBERAL_CUTOFF")) {
  SWEEP_ILLIBERAL_CUTOFF <- -Inf
}
# Passed explicitly into every child script below. run_script_with() builds a
# FRESH environment per call, so a toggle merely set in this script's own scope
# would not reach 11 or 12 -- they would silently fall back to their own default
# and the sweep would mix conventions.
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- FALSE
}
build_suffix <- if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "" else "_exclyr"

# Run one of the pipeline scripts in a fresh environment with the given toggle
# overrides pre-defined. The scripts all guard their toggles with
# `if (!exists(...))`, so anything set here wins and everything else falls back
# to that script's own default. A fresh env per call also means no state leaks
# between windows.
run_script_with <- function(script, overrides) {
  env <- new.env(parent = globalenv())
  for (nm in names(overrides)) {
    assign(nm, overrides[[nm]], envir = env)
  }
  sys.source(file.path(scripts_dir, script), envir = env)
  invisible(env)
}

# Re-pool an existing sweep without re-running anything. The per-run output is
# deterministic (verified: rebuilding a build gives a byte-identical file), so
# when the run folders are already on disk and only the POOLING needs redoing --
# after a change to a figure, a table column, or the sweep's own folder name --
# re-estimating 30 identical specs is pure waste. Set FALSE for that case; the
# collection step below then checks every expected run is present and errors
# rather than silently pooling a partial set.
if (!exists("SWEEP_REESTIMATE")) {
  SWEEP_REESTIMATE <- FALSE
}

# ------------------------------------------------------------------------------
# Build + estimate
# ------------------------------------------------------------------------------

for (n in if (SWEEP_REESTIMATE) WINDOWS else integer(0)) {
  build_path <- file.path(
    data_dir,
    "rdd_build",
    sprintf("rdd_%s_w%d%s.rds", SWEEP_INSTRUMENT, n, build_suffix)
  )
  if (file.exists(build_path)) {
    cat(sprintf("[w=%d] build cached, skipping\n", n))
  } else {
    cat(sprintf("[w=%d] building...\n", n))
    run_script_with(
      "11_build_rdd_data.R",
      list(
        ILLIBERALISM_VAR = SWEEP_INSTRUMENT,
        BACKSLIDING_WINDOW_YEARS = n,
        TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR = TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR
      )
    )
  }

  for (trt in SWEEP_TREATMENTS) {
    cat(sprintf("[w=%d] estimating %s...\n", n, trt))
    run_script_with(
      "12_rdd_analysis.R",
      list(
        ILLIBERALISM_VAR = SWEEP_INSTRUMENT,
        BACKSLIDING_WINDOW_YEARS = n,
        TREATMENT_VAR = trt,
        TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR = TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR,
        SCORE_GAP_MIN = SWEEP_SCORE_GAP_MIN,
        ILLIBERAL_CUTOFF = SWEEP_ILLIBERAL_CUTOFF,
        # Only the main treatment definition gets the full figure set; the other
        # two contribute numbers to the sweep plots below.
        MAKE_PLOTS = identical(trt, SWEEP_TREATMENTS[1])
      )
    )
  }
}

# ------------------------------------------------------------------------------
# Pool every run's results
# ------------------------------------------------------------------------------

# The folder name must carry every axis this sweep holds fixed: the instrument,
# both sample-restriction thresholds, and the window convention. Omitting the
# instrument (as an earlier version did, when there was only ever one) makes
# consecutive sweeps over different instruments silently overwrite one another
# -- three runs leave one folder holding only the last instrument's numbers.
sweep_name <- paste0(
  "window_sweep_",
  INSTRUMENT_LABELS[[SWEEP_INSTRUMENT]],
  "_gap",
  fmt_slug_num(SWEEP_SCORE_GAP_MIN),
  "_illib",
  fmt_slug_num(SWEEP_ILLIBERAL_CUTOFF),
  build_suffix
)
out_dir <- sweep_dir(sweep_name)

collect <- function(file_name) {
  missing <- character(0)
  out <- map_dfr(WINDOWS, function(n) {
    map_dfr(SWEEP_TREATMENTS, function(trt) {
      cfg <- list(
        instrument = SWEEP_INSTRUMENT,
        window = n,
        treatment = trt,
        score_gap_min = SWEEP_SCORE_GAP_MIN,
        illiberal_cutoff = SWEEP_ILLIBERAL_CUTOFF,
        incl_election_year = TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR
      )
      path <- file.path(RUNS_ROOT, run_slug(cfg), file_name)
      if (!file.exists(path)) {
        missing <<- c(missing, run_slug(cfg))
        return(NULL)
      }
      read_csv(path, show_col_types = FALSE) |>
        mutate(window = n, treatment_def = trt, .before = 1)
    })
  })
  if (length(missing) > 0) {
    stop(
      "Cannot pool: ",
      length(missing),
      " expected run(s) are missing ",
      "their ",
      file_name,
      ", e.g. ",
      missing[1],
      if (SWEEP_REESTIMATE) {
        ". A run must have failed above."
      } else {
        ". SWEEP_REESTIMATE is FALSE, so nothing was estimated -- set it TRUE."
      },
      call. = FALSE
    )
  }
  out
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
window_label <- sprintf(
  "[%s, election_year + N]",
  if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) {
    "election_year"
  } else {
    "election_year + 1"
  }
)

trt_display <- unname(TREATMENT_DISPLAY[SWEEP_TREATMENTS])
label_treatment <- function(x) {
  factor(TREATMENT_DISPLAY[x], levels = trt_display)
}

# Okabe-Ito, same fixed order as 12_rdd_analysis.R so a treatment definition
# keeps its colour across every figure in the project.
TRT_COLORS <- setNames(
  c("#0072B2", "#D55E00", "#009E73")[seq_along(SWEEP_TREATMENTS)],
  trt_display
)
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

p_fs <- ggplot(
  fs_plot_data,
  aes(x = window, y = first_stage_coef, colour = treatment_label)
) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_ribbon(
    aes(ymin = lo, ymax = hi, fill = treatment_label),
    alpha = 0.13,
    colour = NA
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
      "Instrument: %s | Sample: %s | Treatment window: %s. Band = 95%% CI.",
      INSTRUMENT_DISPLAY[[SWEEP_INSTRUMENT]],
      restriction_label,
      window_label
    ),
    x = "Window length N (years after the election)",
    y = "RD estimate at the cutoff",
    colour = NULL,
    fill = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 7.5)
  )

ggsave(
  file.path(out_dir, "window_first_stage.png"),
  p_fs,
  width = 11,
  height = 3.4,
  dpi = 150
)
cat("Saved window_first_stage.png\n")

save_table_html(
  first_stage |>
    transmute(
      Window = window,
      Treatment = unname(TREATMENT_DISPLAY[treatment_def]),
      N = n_first_stage,
      `N treated` = n_first_stage_treated,
      `First stage` = fmt_est(
        first_stage_coef,
        first_stage_se,
        first_stage_pval
      ),
      Bandwidth = round(first_stage_bandwidth, 2)
    ) |>
    arrange(Treatment, Window),
  file.path(out_dir, "window_first_stage.html"),
  "First stage by window length",
  sprintf(
    "Instrument: %s | Sample: %s | Treatment window: %s",
    INSTRUMENT_DISPLAY[[SWEEP_INSTRUMENT]],
    restriction_label,
    window_label
  ),
  note = paste(SIG_FOOTNOTE, TREATED_FOOTNOTE)
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
#
# `vars` restricts the figure to one OUTCOME FAMILY. A single sheet carrying all
# 18 outcomes is too dense to read, and the families (economic, inequality,
# fiscal, institutions) are the groupings a reader actually compares within.
outcome_sweep_plot <- function(
  est_col,
  se_col,
  pval_col,
  title,
  file_name,
  clip = FALSE,
  by_treatment = TRUE,
  vars = ALL_OUTCOME_VARS
) {
  results_used <- if (by_treatment) {
    results
  } else {
    results |> filter(treatment_def == SWEEP_TREATMENTS[1])
  }
  dat <- results_used |>
    filter(outcome %in% vars) |>
    mutate(
      treatment_label = label_treatment(treatment_def),
      # Facets are ordered by the family's own order, not alphabetically, and
      # labelled with the panel title so a series name like "PWT" or
      # "Disposable" is intelligible outside its own panel.
      outcome_facet = factor(
        outcome_full_label(outcome),
        levels = unique(outcome_full_label(vars))
      ),
      est = .data[[est_col]],
      lo = .data[[est_col]] - 1.96 * .data[[se_col]],
      hi = .data[[est_col]] + 1.96 * .data[[se_col]]
    ) |>
    filter(!is.na(est))
  if (nrow(dat) == 0) {
    return(invisible(NULL))
  }
  if (clip) {
    dat <- clip_to_point_range(dat)
  }
  if (!by_treatment) {
    dat$treatment_label <- factor("Reduced form", levels = "Reduced form")
  }
  series_colors <- if (by_treatment) {
    TRT_COLORS
  } else {
    c("Reduced form" = SERIES_COLORS[1])
  }
  series_shapes <- if (by_treatment) {
    TRT_SHAPES
  } else {
    c("Reduced form" = SERIES_SHAPES[1])
  }

  p <- ggplot(dat, aes(x = window, y = est, colour = treatment_label)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.35) +
    geom_ribbon(
      aes(ymin = lo, ymax = hi, fill = treatment_label),
      alpha = 0.10,
      colour = NA
    ) +
    geom_line(linewidth = 0.6) +
    geom_point(aes(shape = treatment_label), size = 1.2) +
    facet_wrap(
      ~outcome_facet,
      scales = "free_y",
      ncol = min(3, n_distinct(dat$outcome_facet)),
      labeller = label_wrap_gen(30)
    ) +
    scale_colour_manual(values = series_colors) +
    scale_fill_manual(values = series_colors) +
    scale_shape_manual(values = series_shapes) +
    scale_x_continuous(breaks = WINDOWS) +
    labs(
      title = title,
      # Hard-wrapped: the per-family figures are much narrower than the old
      # all-outcomes sheet, and ggplot silently truncates an over-long subtitle
      # at the panel edge rather than wrapping it.
      subtitle = paste(
        strwrap(
          sprintf(
            "Instrument: %s | Sample: %s | Treatment window: %s. Band = 95%% CI, free y scale per outcome.%s%s",
            INSTRUMENT_DISPLAY[[SWEEP_INSTRUMENT]],
            restriction_label,
            window_label,
            if (clip) {
              " CIs clipped to the span of the point estimates; full values in window_sweep_results.csv."
            } else {
              ""
            },
            if (!by_treatment) {
              " The reduced form does not depend on the treatment definition, so there is one line."
            } else {
              ""
            }
          ),
          width = max(60, 22 * min(3, length(vars)))
        ),
        collapse = "\n"
      ),
      x = "Window length N (years after the election)",
      y = "Estimate at the cutoff",
      colour = NULL,
      fill = NULL,
      shape = NULL
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = if (by_treatment) "bottom" else "none",
      strip.text = element_text(size = 7)
    )
  # Size to the family: a 2-outcome family shouldn't be stretched across the
  # same 12 inches a 7-outcome one needs.
  n_facets <- n_distinct(dat$outcome_facet)
  n_cols <- min(3, n_facets)
  n_rows <- ceiling(n_facets / n_cols)
  ggsave(
    file.path(out_dir, file_name),
    p,
    width = 3.1 * n_cols + 0.8,
    height = 2.3 * n_rows + 1.3,
    dpi = 150,
    limitsize = FALSE
  )
  cat(sprintf("Saved %s\n", file_name))
}

for (fam in names(OUTCOME_FAMILIES)) {
  fam_vars <- OUTCOME_FAMILIES[[fam]]
  fam_title <- unname(OUTCOME_FAMILY_TITLES[fam])

  # The reduced form regresses the outcome on the running variable with no
  # treatment term, so it is numerically IDENTICAL across all three treatment
  # definitions -- verified: 0 of 180 window x outcome cells differ. Drawn as
  # one series rather than three copies of the same line.
  outcome_sweep_plot(
    "rd_estimate",
    "rd_se",
    "rd_pval",
    sprintf("Reduced-form RD by window length: %s", fam_title),
    sprintf("window_rdd_%s.png", fam),
    by_treatment = FALSE,
    vars = fam_vars
  )

  # The fuzzy LATE rescales the reduced form by the first stage, so it DOES
  # depend on the treatment definition -- all 180 cells differ, so all three
  # are drawn.
  outcome_sweep_plot(
    "late_estimate",
    "late_se",
    "late_pval",
    sprintf("Fuzzy RD (LATE) by window length: %s", fam_title),
    sprintf("window_fuzzy_%s.png", fam),
    clip = TRUE,
    vars = fam_vars
  )
}

message("Window sweep written to ", out_dir)
