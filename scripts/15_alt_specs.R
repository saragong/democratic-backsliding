# ==============================================================================
# Alternative instrument and treatment definitions.
#
# Crosses every candidate INSTRUMENT (which party-level score decides who the
# "illiberal" side of the top 2 is, and therefore the sign of the running
# variable) against every candidate TREATMENT definition, at a fixed window, and
# lays the first stage out as one grid. The question this answers is narrow and
# specific: does redefining the instrument, or extending/continuous-ifying the
# treatment, do anything for the first-stage power problem?
#
#   Instruments (to-do 6): anti-pluralism (the original), populism, economic
#     left (negated v2pariglef, so LEFT is the high end -- flagged everywhere it
#     appears), anti-elitism, GAL-TAN.
#   Treatments (to-dos 4 and 5): ERT episode; ERT OR DDCG reversal; continuous
#     decline in V-Dem polyarchy.
#
# Run twice by default, at two sample definitions:
#   "all"  every scored election
#   "q50"  only elections where the more-illiberal party is in the top half of
#          THAT INSTRUMENT's own score distribution. A fixed absolute cutoff
#          (the 0.6 used for v2xpa_antiplural) is meaningless across
#          instruments on different scales, so the comparable restriction is a
#          quantile of each instrument's own distribution.
#
# Output: output/runs/_sweeps/alt_specs_<sample>/
#           alt_specs_results.csv, alt_specs_first_stage.csv
#           grid_first_stage.html, grid_gdp_late.html, grid_gdp_rf.html
#           ddcg_contribution.{csv,html}
# ==============================================================================

library(tidyverse)
library(here)
library(gt)

source(here::here("scripts", "rdd_helpers.R"))

data_dir <- here::here("data")
scripts_dir <- here::here("scripts")

if (!exists("ALT_INSTRUMENTS")) {
  ALT_INSTRUMENTS <- c(
    "v2xpa_antiplural", "v2xpa_popul", "v2pariglef_neg", "v2paanteli", "ep_galtan"
  )
}
if (!exists("ALT_TREATMENTS")) {
  ALT_TREATMENTS <- c(
    "backsliding_Nyr", "backsliding_union_Nyr", "polyarchy_decline"
  )
}
if (!exists("ALT_WINDOW")) ALT_WINDOW <- 5
# Sample definitions to run the whole grid at. Each entry is either the literal
# "all" (no restriction) or a threshold spec that parse_threshold() understands
# -- in practice a quantile such as "q50", which is the only form that means the
# same thing across instruments on different scales. An ABSOLUTE number would
# also parse here, but it makes the grid a comparison of arbitrary numeric
# thresholds rather than of instruments, so it is not the default.
if (!exists("ALT_SAMPLES")) ALT_SAMPLES <- c("all", "q50")
if (!exists("ALT_OUTCOME")) ALT_OUTCOME <- "Y_gdp_growth"
# Passed explicitly into every child script below: run_script_with() builds a
# FRESH environment per call, so a toggle merely set in this script's scope
# would not reach 11 or 12 and the grid would silently mix conventions.
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- TRUE
}
build_suffix <- if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "" else "_exclyr"

# Every read of a build goes through this, so the suffix can never be forgotten
# at one of the three call sites below.
alt_build_path <- function(instr) {
  file.path(
    data_dir, "rdd_build",
    sprintf("rdd_%s_w%d%s.rds", instr, ALT_WINDOW, build_suffix)
  )
}

run_script_with <- function(script, overrides) {
  env <- new.env(parent = globalenv())
  for (nm in names(overrides)) assign(nm, overrides[[nm]], envir = env)
  sys.source(file.path(scripts_dir, script), envir = env)
  invisible(env)
}

# ------------------------------------------------------------------------------
# Builds: one per instrument at ALT_WINDOW (all three treatment definitions and
# every outcome are computed inside a single build, so treatment is free).
# ------------------------------------------------------------------------------

for (instr in ALT_INSTRUMENTS) {
  build_path <- alt_build_path(instr)
  if (file.exists(build_path)) {
    cat(sprintf("[%s] build cached\n", instr))
    next
  }
  cat(sprintf("[%s] building...\n", instr))
  run_script_with("11_build_rdd_data.R", list(
    ILLIBERALISM_VAR = instr,
    BACKSLIDING_WINDOW_YEARS = ALT_WINDOW,
    TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR =
      TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR
  ))
}

# ------------------------------------------------------------------------------
# Estimate every (instrument x treatment) cell, once per sample definition
# ------------------------------------------------------------------------------

for (samp in ALT_SAMPLES) {
  # "all" is a label, not a spec; everything else is passed through verbatim so
  # parse_threshold() is the single arbiter of what is a valid threshold.
  cutoff_spec <- if (identical(samp, "all")) -Inf else samp

  for (instr in ALT_INSTRUMENTS) {
    for (trt in ALT_TREATMENTS) {
      cat(sprintf("[%s | %s | %s] estimating...\n", samp, instr, trt))
      run_script_with("12_rdd_analysis.R", list(
        ILLIBERALISM_VAR = instr,
        BACKSLIDING_WINDOW_YEARS = ALT_WINDOW,
        TREATMENT_VAR = trt,
        TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR =
          TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR,
        SCORE_GAP_MIN = -Inf,
        ILLIBERAL_CUTOFF = cutoff_spec,
        # Only the reference cell gets figures; the rest contribute numbers.
        MAKE_PLOTS = identical(instr, ALT_INSTRUMENTS[1]) &&
          identical(trt, ALT_TREATMENTS[1])
      ))
    }
  }

  # --- pool -------------------------------------------------------------------
  # 12_rdd_analysis.R names its folder by the RESOLVED absolute threshold, not
  # the spec, so re-resolve each instrument's spec here to find the run again.
  # This goes through the same parse_threshold()/resolve_threshold() pair that
  # script 12 uses -- a second, local copy of the quantile arithmetic would be
  # free to drift out of step with it and send us looking in the wrong folder.
  resolved_cutoff <- function(instr) {
    d <- readRDS(alt_build_path(instr))
    resolve_threshold(
      parse_threshold(cutoff_spec, "ILLIBERAL_CUTOFF"),
      d$illiberal_score, "ILLIBERAL_CUTOFF"
    )$absolute
  }
  cutoffs <- setNames(vapply(ALT_INSTRUMENTS, resolved_cutoff, numeric(1)), ALT_INSTRUMENTS)

  collect <- function(file_name) {
    map_dfr(ALT_INSTRUMENTS, function(instr) {
      map_dfr(ALT_TREATMENTS, function(trt) {
        cfg <- list(
          instrument = instr, window = ALT_WINDOW, treatment = trt,
          score_gap_min = -Inf, illiberal_cutoff = cutoffs[[instr]],
          incl_election_year = TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR
        )
        path <- file.path(RUNS_ROOT, run_slug(cfg), file_name)
        if (!file.exists(path)) {
          warning("Missing run output: ", path)
          return(NULL)
        }
        read_csv(path, show_col_types = FALSE) |>
          mutate(instrument = instr, treatment_def = trt, .before = 1)
      })
    })
  }

  out_dir <- sweep_dir(paste0("alt_specs_", samp, build_suffix))
  results <- collect("rdd_results.csv")
  first_stage <- collect("rdd_first_stage_results.csv")
  write_csv(results, file.path(out_dir, "alt_specs_results.csv"))
  write_csv(first_stage, file.path(out_dir, "alt_specs_first_stage.csv"))

  sample_note <- if (samp == "all") {
    "All scored elections."
  } else {
    sprintf(
      "Elections where the more-illiberal party is at or above the %s percentile of that instrument's own score distribution.",
      sub("^q", "", samp)
    )
  }

  # --- grids ------------------------------------------------------------------
  instr_labels <- unname(INSTRUMENT_DISPLAY[ALT_INSTRUMENTS])
  trt_labels <- unname(TREATMENT_DISPLAY[ALT_TREATMENTS])

  fs_mat <- function(col) {
    m <- matrix(
      NA_real_, nrow = length(ALT_INSTRUMENTS), ncol = length(ALT_TREATMENTS),
      dimnames = list(instr_labels, trt_labels)
    )
    for (i in seq_along(ALT_INSTRUMENTS)) {
      for (j in seq_along(ALT_TREATMENTS)) {
        v <- first_stage |>
          filter(instrument == ALT_INSTRUMENTS[i], treatment_def == ALT_TREATMENTS[j])
        if (nrow(v) == 1) m[i, j] <- v[[col]]
      }
    }
    m
  }

  color_grid_table(
    est = fs_mat("first_stage_coef"),
    se = fs_mat("first_stage_se"),
    pval = fs_mat("first_stage_pval"),
    n = fs_mat("n_first_stage"),
    path = file.path(out_dir, "grid_first_stage.html"),
    title = "First stage: instrument x treatment definition",
    subtitle = sprintf(
      "Treatment window: [%s, election_year + %d]. %s",
      if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "election_year" else "election_year + 1",
      ALT_WINDOW, sample_note
    ),
    row_label = "Instrument",
    note = paste(
      "The continuous treatment's coefficient is in polyarchy-index units, not",
      "probability units, so its column is NOT comparable in magnitude to the two",
      "binary columns -- compare significance and sign across that column, magnitude",
      "only down it.",
      "GAL-TAN comes from the CHES bridge (European parties, recent decades) and has",
      "both top-2 scores for only ~4% of the election spine, so its row is typically",
      "all '--': rdrobust cannot fit. See instrument_coverage.html in the",
      "instrument_overlap sweep."
    )
  )

  outcome_mat <- function(col) {
    m <- matrix(
      NA_real_, nrow = length(ALT_INSTRUMENTS), ncol = length(ALT_TREATMENTS),
      dimnames = list(instr_labels, trt_labels)
    )
    for (i in seq_along(ALT_INSTRUMENTS)) {
      for (j in seq_along(ALT_TREATMENTS)) {
        v <- results |>
          filter(
            instrument == ALT_INSTRUMENTS[i],
            treatment_def == ALT_TREATMENTS[j],
            outcome == ALT_OUTCOME
          )
        if (nrow(v) == 1) m[i, j] <- v[[col]]
      }
    }
    m
  }

  for (mt in list(
    list(prefix = "late", coef = "late_estimate", se = "late_se", pval = "late_pval",
         title = "Fuzzy RD (LATE)"),
    list(prefix = "rf", coef = "rd_estimate", se = "rd_se", pval = "rd_pval",
         title = "Reduced-form RD")
  )) {
    color_grid_table(
      est = outcome_mat(mt$coef),
      se = outcome_mat(mt$se),
      pval = outcome_mat(mt$pval),
      n = outcome_mat("n_outcome"),
      path = file.path(out_dir, sprintf("grid_gdp_%s.html", mt$prefix)),
      title = sprintf("%s on %s: instrument x treatment definition", mt$title, ALT_OUTCOME),
      subtitle = sprintf(
        "Treatment window: [%s, election_year + %d]. %s",
        if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "election_year" else "election_year + 1",
        ALT_WINDOW, sample_note
      ),
      row_label = "Instrument",
      note = if (mt$prefix == "rf") {
        "The reduced form does not use the treatment, so every column of this grid is identical by construction; it is shown as a check, not a comparison."
      } else {
        NULL
      }
    )
  }

  message("Alt-spec grids written to ", out_dir)
}

# ------------------------------------------------------------------------------
# How much does the DDCG extension actually add?
#
# The honest denominator for "did extending the treatment fix the power
# problem": DDCG's panel only runs 1960-2010, so for any election whose
# post-election window falls entirely outside that span the union treatment is
# identical to ERT alone by construction.
# ------------------------------------------------------------------------------

ddcg_rows <- map_dfr(ALT_INSTRUMENTS, function(instr) {
  d <- readRDS(alt_build_path(instr))
  tibble(
    instrument = unname(INSTRUMENT_DISPLAY[instr]),
    n_elections = nrow(d),
    n_in_ddcg_coverage = sum(d$ddcg_covered == 1),
    n_treated_ert = sum(d$backsliding_Nyr),
    n_treated_union = sum(d$backsliding_union_Nyr),
    n_added_by_ddcg = sum(d$backsliding_union_Nyr == 1 & d$backsliding_Nyr == 0),
    pct_added = round(
      100 * (sum(d$backsliding_union_Nyr) - sum(d$backsliding_Nyr)) /
        sum(d$backsliding_Nyr), 1
    )
  )
})

ddcg_dir <- sweep_dir(paste0("alt_specs_", ALT_SAMPLES[1], build_suffix))
write_csv(ddcg_rows, file.path(ddcg_dir, "ddcg_contribution.csv"))
save_table_html(
  ddcg_rows |>
    rename(
      Instrument = instrument,
      `N elections` = n_elections,
      `N in DDCG window` = n_in_ddcg_coverage,
      `Treated (ERT)` = n_treated_ert,
      `Treated (ERT or DDCG)` = n_treated_union,
      `Added by DDCG` = n_added_by_ddcg,
      `% added` = pct_added
    ),
  file.path(ddcg_dir, "ddcg_contribution.html"),
  "What the DDCG extension adds to the treatment",
  sprintf("Window: %d yr", ALT_WINDOW),
  note = paste(
    "DDCG's reversal events only span 1960-2010. For elections whose",
    "post-election window falls outside that span the union treatment equals the",
    "ERT treatment by construction, so 'N in DDCG window' is the real denominator",
    "for any power gain."
  )
)
