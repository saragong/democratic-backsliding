# ==============================================================================
# Sample-restriction sweep: how does the first stage move as we restrict the
# election sample along different dimensions?
#
# This is an EXPLICIT SPECIFICATION SEARCH, laid out so it can be read as one.
# The point is to see the whole surface of first-stage estimates a researcher
# could have reached by choosing restrictions, rather than reporting whichever
# cell looks best. Every table says so on its face.
#
# Restriction axes (all orthogonal to which side of the cutoff an election
# landed on, so none of them breaks the RD by selecting on treatment):
#   R1  score_gap        the raw illiberality gap between the top 2
#   R2  illiberal_score  the more-illiberal party's absolute level
#   R3  score_gap_z      that gap standardized by the country's own party
#                        ideology spread
#   R4  prior_backsliding whether to drop elections already following a recent
#                        backsliding episode
#   R5  election_type    presidential / parliamentary / both
#
# Note on R1 and R2: "the gap between the winner and the loser" and "the
# winner's score" are deliberately NOT keyed to who actually won. illiberal_
# score is the higher of the top 2 and score_gap is always >= 0, so both
# restrictions keep elections on BOTH sides of the cutoff. Keying them to the
# actual winner would condition on treatment status and leave one side of the
# cutoff empty.
#
# Output: output/runs/<slug>/restriction_grid/
#           marginal_<axis>.{csv,html}   one row per level of a single axis
#           grid_<A>x<B>_<metric>.{csv,html}  colour-coded pair grids
# ==============================================================================

library(tidyverse)
library(rdrobust)
library(here)
library(gt)

source(here::here("scripts", "rdd_helpers.R"))

data_dir <- here::here("data")

# What this script holds FIXED. The sample restrictions are deliberately absent:
# they are the axes being swept (see AXES below), which is exactly why this
# script must not build a `cfg` and call run_dir(). An earlier version did, with
# score_gap_min = illiberal_cutoff = -Inf standing in for "unrestricted" -- but
# that wrote a run_config.csv asserting a single threshold for a grid that
# varies it, and parked the output inside a folder 12_rdd_analysis.R
# legitimately owns for that configuration. This is a sweep, so it writes to
# _sweeps/ and records its fixed/swept axes in sweep_config.csv.
if (!exists("ILLIBERALISM_VAR")) ILLIBERALISM_VAR <- "v2xpa_antiplural"
if (!exists("BACKSLIDING_WINDOW_YEARS")) BACKSLIDING_WINDOW_YEARS <- 5
if (!exists("TREATMENT_VAR")) TREATMENT_VAR <- "backsliding_Nyr"
# Selects which build to read and is echoed into the sweep folder name -- see
# 11_build_rdd_data.R for what it does.
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- TRUE
}
# The headline outcome the grids report alongside the first stage.
if (!exists("GRID_OUTCOME")) GRID_OUTCOME <- "Y_gdp_growth"

build_suffix <- if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "" else "_exclyr"
build_path <- file.path(
  data_dir, "rdd_build",
  sprintf(
    "rdd_%s_w%d%s.rds",
    ILLIBERALISM_VAR, BACKSLIDING_WINDOW_YEARS, build_suffix
  )
)
if (!file.exists(build_path)) {
  stop("No build at ", build_path, " -- run 11_build_rdd_data.R first.")
}
d_full <- readRDS(build_path)
# Same check 12_rdd_analysis.R makes: a build written before the window toggle
# existed carries no attribute and used the old convention, so a missing
# attribute means FALSE rather than the current default.
build_incl <- attr(d_full, "includes_election_year") %||% FALSE
if (!identical(build_incl, TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR)) {
  stop(
    "Build at ", build_path, " was made with ",
    "TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR = ", build_incl,
    " but this sweep asked for ", TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR,
    ". Rebuild it with 11_build_rdd_data.R."
  )
}
cat(sprintf("Loaded %s (%d elections)\n", basename(build_path), nrow(d_full)))

grid_dir <- sweep_dir(sprintf(
  "restriction_grid_%s_w%d_trt-%s%s",
  INSTRUMENT_LABELS[[ILLIBERALISM_VAR]],
  BACKSLIDING_WINDOW_YEARS,
  TREATMENT_LABELS[[TREATMENT_VAR]],
  build_suffix
))

# ------------------------------------------------------------------------------
# Restriction axes
#
# An axis is a name, a display label, and a list of levels; each level is a
# short label plus a predicate over the data frame. Written this way so adding
# a sixth axis is one list entry, not a new code path.
# ------------------------------------------------------------------------------

# Quintile thresholds of a continuous variable, with "no restriction" first.
# Quantiles come from the FULL build, so they don't drift as other axes bite.
threshold_axis <- function(var, label, data = d_full, probs = c(0.2, 0.4, 0.6, 0.8)) {
  cuts <- unname(quantile(data[[var]], probs = probs, na.rm = TRUE))
  levels <- list(list(label = "all", cutoff = -Inf, fn = function(df) df))
  for (i in seq_along(cuts)) {
    local({
      cut_i <- cuts[i]
      p_i <- probs[i]
      levels[[length(levels) + 1]] <<- list(
        label = sprintf("q%02d (\u2265 %.3g)", round(p_i * 100), cut_i),
        cutoff = cut_i,
        fn = function(df) df[!is.na(df[[var]]) & df[[var]] >= cut_i, , drop = FALSE]
      )
    })
  }
  list(name = var, label = label, levels = levels)
}

AXES <- list(
  R1 = threshold_axis("score_gap", "R1: top-2 illiberality gap"),
  R2 = threshold_axis("illiberal_score", "R2: more-illiberal party's score"),
  R3 = threshold_axis("score_gap_z", "R3: gap, within-country standardized"),
  R4 = list(
    name = "prior_backsliding",
    label = "R4: prior backsliding",
    levels = list(
      list(label = "all", cutoff = NA, fn = function(df) df),
      list(
        label = "no prior episode",
        cutoff = NA,
        fn = function(df) df[df$prior_backsliding == 0, , drop = FALSE]
      )
    )
  ),
  R5 = list(
    name = "election_type",
    label = "R5: election type",
    levels = list(
      list(label = "both", cutoff = NA, fn = function(df) df),
      list(
        label = "presidential",
        cutoff = NA,
        fn = function(df) df[df$election_type == "presidential", , drop = FALSE]
      ),
      list(
        label = "parliamentary",
        cutoff = NA,
        fn = function(df) df[df$election_type == "parliamentary", , drop = FALSE]
      )
    )
  )
)

# ------------------------------------------------------------------------------
# Estimation for one cell
# ------------------------------------------------------------------------------

# First stage, plus the reduced form and fuzzy LATE for GRID_OUTCOME. Returns
# NA-filled fields when rdrobust can't fit -- which happens routinely in the
# thin corners of the grid and is itself informative, so it is rendered as
# "--" rather than dropped.
estimate_cell <- function(df) {
  n <- nrow(df)
  if (n < RD_MIN_OBS) {
    return(tibble(
      n = n, n_treated = NA_integer_,
      fs_coef = NA_real_, fs_se = NA_real_, fs_pval = NA_real_, fs_n = NA_integer_,
      rf_coef = NA_real_, rf_se = NA_real_, rf_pval = NA_real_,
      late_coef = NA_real_, late_se = NA_real_, late_pval = NA_real_
    ))
  }
  fs <- extract_rd(safe_rdrobust(df[[TREATMENT_VAR]], df$running_var))
  rf <- extract_rd(safe_rdrobust(df[[GRID_OUTCOME]], df$running_var))
  late <- extract_rd(safe_rdrobust(
    df[[GRID_OUTCOME]], df$running_var, fuzzy = df[[TREATMENT_VAR]]
  ))
  tibble(
    n = n,
    n_treated = if (TREATMENT_VAR == "polyarchy_decline") {
      NA_integer_
    } else {
      as.integer(sum(df[[TREATMENT_VAR]], na.rm = TRUE))
    },
    fs_coef = fs$coef, fs_se = fs$se, fs_pval = fs$pval, fs_n = fs$N,
    rf_coef = rf$coef, rf_se = rf$se, rf_pval = rf$pval,
    late_coef = late$coef, late_se = late$se, late_pval = late$pval
  )
}

write_sweep_config(
  grid_dir,
  fixed = list(
    instrument = ILLIBERALISM_VAR,
    window = BACKSLIDING_WINDOW_YEARS,
    treatment = TREATMENT_VAR,
    treatment_window_includes_election_year = TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR,
    grid_outcome = GRID_OUTCOME,
    n_elections_unrestricted = nrow(d_full)
  ),
  swept = lapply(AXES, function(ax) {
    sprintf(
      "%s: %s",
      ax$name,
      paste(vapply(ax$levels, function(l) l$label, character(1)), collapse = " | ")
    )
  })
)

SPEC_SEARCH_NOTE <- paste(
  "This table is an explicit specification search: every cell is a different",
  "sample restriction a researcher could have chosen. Read the surface, not the",
  "best cell. No multiple-testing correction is applied."
)

# ------------------------------------------------------------------------------
# Marginal tables, one per axis
# ------------------------------------------------------------------------------

cat("\n=== Marginal tables ===\n")
for (ax_name in names(AXES)) {
  ax <- AXES[[ax_name]]
  rows <- map_dfr(ax$levels, function(lv) {
    bind_cols(tibble(level = lv$label), estimate_cell(lv$fn(d_full)))
  })
  write_csv(rows, file.path(grid_dir, sprintf("marginal_%s.csv", ax_name)))
  save_table_html(
    rows |>
      transmute(
        Level = level,
        N = n,
        `N treated` = n_treated,
        `First stage` = fmt_est(fs_coef, fs_se, fs_pval),
        `RD (reduced form)` = fmt_est(rf_coef, rf_se, rf_pval),
        `Fuzzy LATE` = fmt_est(late_coef, late_se, late_pval)
      ),
    file.path(grid_dir, sprintf("marginal_%s.html", ax_name)),
    ax$label,
    sprintf(
      "Outcome: %s | Treatment: %s | Instrument: %s | Window: %d yr",
      GRID_OUTCOME, TREATMENT_DISPLAY[[TREATMENT_VAR]],
      INSTRUMENT_DISPLAY[[ILLIBERALISM_VAR]], BACKSLIDING_WINDOW_YEARS
    ),
    note = paste(SIG_FOOTNOTE, SPEC_SEARCH_NOTE)
  )
}

# ------------------------------------------------------------------------------
# Pair grids
# ------------------------------------------------------------------------------

# R1 x R3 pairs the raw gap against its within-country standardization -- they
# measure the same idea on different scales, so that grid shows whether the
# standardization is doing any work. The rest pair a gap/level restriction
# against a structural one.
PAIRS <- list(
  c("R1", "R2"),
  c("R1", "R3"),
  c("R2", "R4"),
  c("R2", "R5")
)

METRICS <- list(
  first_stage = list(
    coef = "fs_coef", se = "fs_se", pval = "fs_pval",
    title = "First-stage coefficient"
  ),
  rf = list(
    coef = "rf_coef", se = "rf_se", pval = "rf_pval",
    title = "Reduced-form RD"
  ),
  late = list(
    coef = "late_coef", se = "late_se", pval = "late_pval",
    title = "Fuzzy RD (LATE)"
  )
)

cat("\n=== Pair grids ===\n")
for (pr in PAIRS) {
  ax_r <- AXES[[pr[1]]]
  ax_c <- AXES[[pr[2]]]
  row_labels <- vapply(ax_r$levels, function(l) l$label, character(1))
  col_labels <- vapply(ax_c$levels, function(l) l$label, character(1))

  cells <- map_dfr(seq_along(ax_r$levels), function(i) {
    map_dfr(seq_along(ax_c$levels), function(j) {
      df <- ax_c$levels[[j]]$fn(ax_r$levels[[i]]$fn(d_full))
      bind_cols(
        tibble(row = row_labels[i], col = col_labels[j], i = i, j = j),
        estimate_cell(df)
      )
    })
  })
  write_csv(
    cells |> select(-i, -j),
    file.path(grid_dir, sprintf("grid_%sx%s.csv", pr[1], pr[2]))
  )

  to_matrix <- function(col) {
    m <- matrix(
      NA_real_,
      nrow = length(row_labels), ncol = length(col_labels),
      dimnames = list(row_labels, col_labels)
    )
    m[cbind(cells$i, cells$j)] <- cells[[col]]
    m
  }
  n_mat <- to_matrix("n")

  for (mt_name in names(METRICS)) {
    mt <- METRICS[[mt_name]]
    color_grid_table(
      est = to_matrix(mt$coef),
      se = to_matrix(mt$se),
      pval = to_matrix(mt$pval),
      n = n_mat,
      path = file.path(
        grid_dir, sprintf("grid_%sx%s_%s.html", pr[1], pr[2], mt_name)
      ),
      title = sprintf("%s: %s x %s", mt$title, ax_r$label, ax_c$label),
      subtitle = sprintf(
        "%s | Outcome: %s | Treatment: %s | Window: %d yr",
        INSTRUMENT_DISPLAY[[ILLIBERALISM_VAR]], GRID_OUTCOME,
        TREATMENT_DISPLAY[[TREATMENT_VAR]], BACKSLIDING_WINDOW_YEARS
      ),
      row_label = ax_r$label,
      note = SPEC_SEARCH_NOTE
    )
  }
}

message("Restriction grid written to ", grid_dir)
