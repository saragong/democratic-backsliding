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
# Output: output/runs/_sweeps/restriction_grid_<instr>_w<N>_trt-<treatment>/
#           sweep_config.csv     what was held fixed, what was swept
#           marginal_all.{csv,html}
#             ONE table, one section per axis: N, N treated, first stage,
#             reduced-form RD and fuzzy LATE at each level of that axis
#           grid_first_stage.html, grid_rf.html, grid_late.html
#             ONE table per analysis type, one section per pair of axes, each
#             cell colour-coded by coefficient magnitude
#           grid_all_pairs.csv   every pair-grid cell in long form
#
#         output/runs/_sweeps/sample_composition_<instr>_w<N>/
#           composition_all.{csv,html}, composition_narrow<PP>.{csv,html}
#             For each restriction level: N surviving, and how it splits between
#             the illiberal side winning and losing -- i.e. the two sides of the
#             RD cutoff. Separate folder because these depend only on the
#             instrument, the window and the window convention, NOT on the
#             treatment definition or the outcome, so keying them by treatment
#             would write three identical copies.
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

# The CEILING counterpart: tightens from above, so the levels run from "no
# restriction" down to a low quantile. Used for the less-illiberal party's
# score, where the restriction of interest is "the opponent is NOT itself
# illiberal" -- the mirror of R2.
#
# Levels are the LOW quantiles (q80 down to q20) so that, as in every other
# axis here, moving down the rows means a tighter sample.
ceiling_axis <- function(var, label, data = d_full, probs = c(0.8, 0.6, 0.4, 0.2)) {
  cuts <- unname(quantile(data[[var]], probs = probs, na.rm = TRUE))
  levels <- list(list(label = "all", cutoff = Inf, fn = function(df) df))
  for (i in seq_along(cuts)) {
    local({
      cut_i <- cuts[i]
      p_i <- probs[i]
      levels[[length(levels) + 1]] <<- list(
        label = sprintf("q%02d (\u2264 %.3g)", round(p_i * 100), cut_i),
        cutoff = cut_i,
        fn = function(df) df[!is.na(df[[var]]) & df[[var]] <= cut_i, , drop = FALSE]
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
  # The opponent-side ceiling. On its own it asks "does the result survive
  # dropping elections where BOTH top-2 parties are illiberal?"; crossed with
  # R2 (a floor on the more-illiberal party) it is the "one illiberal, one
  # not" restriction that 01g_populist_threshold.R calibrates a number for.
  # Here it is swept over quantiles instead, so the answer does not depend on
  # PopuList, which covers only 31 European countries.
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
  ),
  # The opponent-side ceiling. On its own it asks "does the result survive
  # dropping elections where BOTH top-2 parties are illiberal?"; crossed with
  # R2 (a floor on the more-illiberal party) it is the "one illiberal, one
  # not" restriction that 01g_populist_threshold.R calibrates a number for.
  # Swept over quantiles here, so this axis does not depend on PopuList, which
  # covers only 31 European countries.
  R6 = ceiling_axis("other_score", "R6: less-illiberal party's score (ceiling)")
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
    # Via the shared helper, not a local sum: it handles the continuous
    # treatment (counting elections where polyarchy strictly fell) instead of
    # returning NA, and it applies the same non-missing mask the first-stage
    # fit above uses, so N treated is always a subset of the reported N.
    n_treated = rd_n_treated(
      df[[TREATMENT_VAR]],
      df$running_var,
      df[[TREATMENT_VAR]]
    ),
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

# ------------------------------------------------------------------------------
# One marginal table across every axis
#
# Previously one file per axis. They all share the same columns, so they belong
# in a single table with the restriction as a left-hand column and each axis as
# its own section.
# ------------------------------------------------------------------------------

cat("\n=== Marginal table ===\n")

marginal <- map_dfr(names(AXES), function(ax_name) {
  ax <- AXES[[ax_name]]
  map_dfr(ax$levels, function(lv) {
    bind_cols(
      tibble(axis = ax_name, axis_label = ax$label, level = lv$label),
      estimate_cell(lv$fn(d_full))
    )
  })
})
write_csv(marginal, file.path(grid_dir, "marginal_all.csv"))

marginal_tbl <- marginal |>
  transmute(
    Restriction = axis_label,
    Level = level,
    N = n,
    `N treated` = n_treated,
    `First stage` = fmt_est(fs_coef, fs_se, fs_pval),
    `RD (reduced form)` = fmt_est(rf_coef, rf_se, rf_pval),
    `Fuzzy RD (LATE)` = fmt_est(late_coef, late_se, late_pval)
  )

gt_marginal <- marginal_tbl |>
  gt(groupname_col = "Restriction") |>
  tab_header(
    title = "Sample-restriction sweep: how the first stage and the GDP-growth RD move as each restriction tightens",
    subtitle = sprintf(
      "One section per restriction axis, each tightening from no restriction downwards. Instrument: %s | Treatment: %s | Window: [%s, election_year + %d] | Outcome: %s",
      INSTRUMENT_DISPLAY[[ILLIBERALISM_VAR]],
      TREATMENT_DISPLAY[[TREATMENT_VAR]],
      if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "election_year" else "election_year + 1",
      BACKSLIDING_WINDOW_YEARS,
      GRID_OUTCOME
    )
  ) |>
  cols_align(align = "center", columns = c(N, `N treated`)) |>
  cols_align(align = "left", columns = Level) |>
  # Section headers spanning the full width, rather than a left-hand column --
  # the header names the restriction clearly enough on its own.
  tab_options(
    row_group.as_column = FALSE,
    row_group.font.weight = "bold",
    row_group.background.color = "#f0f0f0",
    row_group.border.top.width = px(2),
    row_group.border.top.color = "black",
    row_group.border.bottom.width = px(1),
    row_group.border.bottom.color = "#888888"
  ) |>
  opt_row_striping() |>
  apply_table_style() |>
  tab_source_note(source_note = paste(SIG_FOOTNOTE, SPEC_SEARCH_NOTE))

gtsave(gt_marginal, file.path(grid_dir, "marginal_all.html"))
cat(sprintf("Saved %s\n", file.path(grid_dir, "marginal_all.html")))

# ------------------------------------------------------------------------------
# Sample composition across the cutoff
#
# For every restriction level: how many elections survive, and how they split
# between the illiberal side WINNING and LOSING -- i.e. the two sides of the RD
# cutoff. running_var is the illiberal side's share minus the other's, so
# running_var > 0 is exactly "the more-illiberal of the top 2 won".
#
# This is the diagnostic that says whether a restriction is quietly starving one
# side of the comparison. None of these restrictions conditions on who won (they
# are all functions of the top-2 SCORES, the country's score spread, prior
# episodes, or the office type), so identification survives any of them. But a
# more-illiberal candidate with a very high absolute score tends to face a
# weaker opponent, so those races are less often close -- tightening R2 thins
# the losing side faster than the winning side, and the standard errors in the
# grid corners blow up faster than the point estimates as a result.
# ------------------------------------------------------------------------------

cat("\n=== Sample composition ===\n")

# These tables go in their OWN sweep folder, not the treatment-keyed
# restriction_grid one, because they don't depend on the treatment definition or
# on the outcome -- nothing below reads TREATMENT_VAR or GRID_OUTCOME. Writing
# them alongside the grids would produce three identical copies, one per
# treatment.
#
# What they do depend on, verified against the builds:
#   instrument   R1 (score_gap), R2 (illiberal_score), R3 (score_gap_z) and
#                running_var itself are all functions of the instrument.
#   window +     R4 (prior_backsliding) is counted over the N years before the
#   convention   treatment window opens, so it moves with both (105 flagged at
#                w5/exclyr, 70 at w3/exclyr, 107 at w5/inclyr).
#   neither      R5 (election_type).
# So the folder is keyed on instrument x window x convention.
composition_dir <- sweep_dir(sprintf(
  "sample_composition_%s_w%d%s",
  INSTRUMENT_LABELS[[ILLIBERALISM_VAR]],
  BACKSLIDING_WINDOW_YEARS,
  build_suffix
))
write_sweep_config(
  composition_dir,
  fixed = list(
    instrument = ILLIBERALISM_VAR,
    window = BACKSLIDING_WINDOW_YEARS,
    treatment_window_includes_election_year = TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR,
    depends_on_treatment = FALSE,
    depends_on_outcome = FALSE,
    n_elections_unrestricted = nrow(d_full)
  ),
  swept = list(
    restrictions = paste(names(AXES), collapse = ", "),
    samples = "all elections; narrow window"
  )
)

# Half-width of the narrow-election window, in percentage points of vote or seat
# share. The second composition table restricts to |running_var| <= this, i.e.
# the elections near the cutoff that a local RD estimate is actually built from.
if (!exists("COMPOSITION_NARROW_PP")) COMPOSITION_NARROW_PP <- 10

composition_for <- function(base) {
  map_dfr(names(AXES), function(ax_name) {
    ax <- AXES[[ax_name]]
    map_dfr(ax$levels, function(lv) {
      # The restriction is applied to the already-narrowed frame, so each row is
      # "this restriction AND within the margin", not one or the other.
      df <- lv$fn(base)
      # Exact ties are impossible in this spine (checked: 0 elections have
      # running_var == 0), but split strictly and carry a tie column so a future
      # spine with ties would show them rather than silently folding them in.
      tibble(
        axis = ax_name,
        axis_label = ax$label,
        level = lv$label,
        n = nrow(df),
        n_illiberal_won = sum(df$running_var > 0),
        n_illiberal_lost = sum(df$running_var < 0),
        n_tied = sum(df$running_var == 0),
        pct_illiberal_won = if (nrow(df) > 0) {
          round(100 * mean(df$running_var > 0), 1)
        } else {
          NA_real_
        }
      )
    })
  })
}

save_composition <- function(composition, file_stem, sample_note) {
  write_csv(composition, file.path(composition_dir, paste0(file_stem, ".csv")))

  gt_tbl <- composition |>
    transmute(
      Restriction = axis_label,
      Level = level,
      N = n,
      `Illiberal won` = n_illiberal_won,
      `Illiberal lost` = n_illiberal_lost,
      `% won` = pct_illiberal_won
    ) |>
    gt(groupname_col = "Restriction") |>
    tab_header(
      title = "Sample composition either side of the RD cutoff, by restriction level",
      subtitle = sprintf(
        "%s The illiberal side 'won' when running_var > 0, i.e. the more-illiberal of the top 2 took the larger vote or seat share. Instrument: %s | Window: [%s, election_year + %d]",
        sample_note,
        INSTRUMENT_DISPLAY[[ILLIBERALISM_VAR]],
        if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) {
          "election_year"
        } else {
          "election_year + 1"
        },
        BACKSLIDING_WINDOW_YEARS
      )
    ) |>
  cols_align(align = "center", columns = -c(1, 2)) |>
  cols_align(align = "left", columns = Level) |>
  tab_options(
    row_group.as_column = FALSE,
    row_group.font.weight = "bold",
    row_group.background.color = "#f0f0f0",
    row_group.border.top.width = px(2),
    row_group.border.top.color = "black",
    row_group.border.bottom.width = px(1),
    row_group.border.bottom.color = "#888888"
  ) |>
    opt_row_striping() |>
    apply_table_style() |>
    tab_source_note(source_note = paste(
      "None of these restrictions conditions on who won, so all of them leave the",
      "RD identified. The columns show whether a restriction nevertheless thins one",
      "side of the cutoff faster than the other, which drives how quickly the",
      "standard errors deteriorate as it tightens.",
      if (sum(composition$n_tied) == 0) {
        "No election in this spine has an exact tie in vote or seat share."
      } else {
        sprintf(
          "%d election-level exact ties are excluded from both columns.",
          sum(composition$n_tied)
        )
      }
    ))

  gtsave(gt_tbl, file.path(composition_dir, paste0(file_stem, ".html")))
  cat(sprintf("Saved %s\n", file.path(composition_dir, paste0(file_stem, ".html"))))
  invisible(composition)
}

# Full sample.
save_composition(
  composition_for(d_full),
  "composition_all",
  "All scored elections."
)

# Narrow elections only -- the local neighbourhood an RD estimate is actually
# built from, so this is the split that matters for whether a restriction leaves
# enough observations on BOTH sides of the cutoff to compare.
d_narrow <- d_full |> filter(abs(running_var) <= COMPOSITION_NARROW_PP)
cat(sprintf(
  "Narrow window |running_var| <= %g pp: %d / %d elections\n",
  COMPOSITION_NARROW_PP,
  nrow(d_narrow),
  nrow(d_full)
))
save_composition(
  composition_for(d_narrow),
  sprintf("composition_narrow%g", COMPOSITION_NARROW_PP),
  sprintf(
    "Narrow elections only: |vote or seat share margin| <= %g pp (%d of %d elections).",
    COMPOSITION_NARROW_PP,
    nrow(d_narrow),
    nrow(d_full)
  )
)

# ------------------------------------------------------------------------------
# Pair grids
# ------------------------------------------------------------------------------

# The three pairs that actually carry information about the first stage.
#
# R2 (the more-illiberal party's absolute score) is the restriction that moves
# the first stage, so it is crossed against both of the others. R3 (the
# within-country standardized gap) is preferred over R1 (the raw gap) here:
# the two measure the same idea on different scales, so R1 x R3 was largely
# degenerate -- tightening one stopped biting once the other had, leaving
# identical cells repeated across a row. R1 and R5 are still swept in the
# marginal table, which covers every axis; they are just not worth a grid.
PAIRS <- list(
  c("R2", "R3"),
  c("R2", "R4"),
  c("R3", "R4"),
  # The floor on the illiberal side crossed with the ceiling on the other: the
  # cells below the diagonal of this grid are the "one illiberal, one not"
  # samples, at every threshold rather than at the one PopuList happens to
  # imply.
  c("R2", "R6")
)

cat("\n=== Pair grids ===\n")

# Estimate every cell of every pair once, then render the same numbers three
# times (one table per analysis type) rather than re-fitting per metric.
pair_blocks <- map(PAIRS, function(pr) {
  ax_r <- AXES[[pr[1]]]
  ax_c <- AXES[[pr[2]]]
  row_labels <- vapply(ax_r$levels, function(l) l$label, character(1))
  col_labels <- vapply(ax_c$levels, function(l) l$label, character(1))

  cells <- map_dfr(seq_along(ax_r$levels), function(i) {
    map_dfr(seq_along(ax_c$levels), function(j) {
      df <- ax_c$levels[[j]]$fn(ax_r$levels[[i]]$fn(d_full))
      bind_cols(
        tibble(
          pair = paste0(pr[1], " x ", pr[2]),
          row = row_labels[i], col = col_labels[j], i = i, j = j
        ),
        estimate_cell(df)
      )
    })
  })

  to_matrix <- function(col) {
    m <- matrix(
      NA_real_,
      nrow = length(row_labels), ncol = length(col_labels),
      dimnames = list(row_labels, col_labels)
    )
    m[cbind(cells$i, cells$j)] <- cells[[col]]
    m
  }

  list(
    pair = paste0(pr[1], " x ", pr[2]),
    row_axis = ax_r$label,
    col_axis = ax_c$label,
    row_labels = row_labels,
    col_labels = col_labels,
    cells = cells |> select(-i, -j),
    mat = to_matrix
  )
})

# Long-format record of every cell of every pair, for anything that wants the
# numbers rather than the rendered table.
write_csv(
  map_dfr(pair_blocks, function(b) b$cells),
  file.path(grid_dir, "grid_all_pairs.csv")
)

METRICS <- list(
  first_stage = list(
    coef = "fs_coef", se = "fs_se", pval = "fs_pval",
    title = "First-stage coefficient across pairs of sample restrictions"
  ),
  rf = list(
    coef = "rf_coef", se = "rf_se", pval = "rf_pval",
    title = "Reduced-form RD across pairs of sample restrictions"
  ),
  late = list(
    coef = "late_coef", se = "late_se", pval = "late_pval",
    title = "Fuzzy RD (LATE) across pairs of sample restrictions"
  )
)

for (mt_name in names(METRICS)) {
  mt <- METRICS[[mt_name]]
  blocks <- map(pair_blocks, function(b) {
    list(
      label = b$pair,
      row_axis = b$row_axis,
      col_axis = b$col_axis,
      row_labels = b$row_labels,
      col_labels = b$col_labels,
      est = b$mat(mt$coef),
      se = b$mat(mt$se),
      pval = b$mat(mt$pval),
      n = b$mat("n")
    )
  })
  color_grid_blocks(
    blocks,
    path = file.path(grid_dir, sprintf("grid_%s.html", mt_name)),
    title = mt$title,
    subtitle = sprintf(
      "%s | Outcome: %s | Treatment: %s | Window: [%s, election_year + %d]",
      INSTRUMENT_DISPLAY[[ILLIBERALISM_VAR]], GRID_OUTCOME,
      TREATMENT_DISPLAY[[TREATMENT_VAR]],
      if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "election_year" else "election_year + 1",
      BACKSLIDING_WINDOW_YEARS
    ),
    note = SPEC_SEARCH_NOTE
  )
}

message("Restriction grid written to ", grid_dir)
