# ==============================================================================
# Reduced-form RDD on the decade x OECD grid, for every outcome.
#
# The headline RDD pools 1970-2019 and both country groups. Scripts 18 and 19
# showed those two regimes are not alike -- anti-pluralism and populism
# correlate +0.53 in the OECD by the 2010s and about zero outside it, and the
# OECD sits an entire 0.38 lower on anti-pluralism -- so a pooled estimate may
# be averaging over genuinely different worlds. This cuts the main outcomes
# the same 2 x 5 way as 18 and 19, so the three are read together.
#
# REDUCED FORM ONLY, and the reason is in the data rather than a preference.
# Treated elections per cell are 0, 1, 1, 3, 13 across the OECD decades and
# 8, 6, 19, 38, 46 outside it. A first stage needs variation in treatment
# near the cutoff; in half these cells there is essentially no treatment at
# all, so neither the first stage nor the fuzzy LATE exists. Treated counts
# are written out anyway, so the omission is evidenced rather than asserted.
#
# Two samples, ten windows, ten cells, every outcome. The restriction sample
# goes through the same parse_threshold/apply_threshold path 12_rdd_analysis.R
# uses, so "one side illiberal, the other not" is defined in exactly one place
# and cannot drift between the pooled runs and these.
#
# Plotting is deliberately not here. The CSV is the deliverable and is shaped
# for a plotting pass to read directly.
#
#   Rscript --no-init-file scripts/21_cell_rdd.R
#
# Output: output/runs/_sweeps/cell_rdd_<instr><suffix>/
#           cell_rdd_results.csv   sample x window x group x decade x outcome
#           cell_counts.csv        cell sizes and treated counts per window
#           grid_cells_w5_<family>.html   a first look at w5
# ==============================================================================

library(tidyverse)
library(rdrobust)
library(here)
library(gt)

source(here::here("scripts", "rdd_helpers.R"))
# oecd_group() and decade_label(): general splits on a country code and a
# year, living in vparty_helpers.R because the raw-V-Party scripts needed
# them first.
source(here::here("scripts", "vparty_helpers.R"))

data_dir <- here::here("data")

# ---- toggles -----------------------------------------------------------------

if (!exists("CELL_INSTRUMENT")) CELL_INSTRUMENT <- "v2xpa_antiplural"
if (!exists("CELL_WINDOWS")) CELL_WINDOWS <- 1:10
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- FALSE
}
if (!exists("PLACEBO_PRE_WINDOW")) PLACEBO_PRE_WINDOW <- FALSE

# The samples to cut each cell to. Each is a (score_gap, illiberal, other)
# triple in the same spec vocabulary 12_rdd_analysis.R accepts, so "popucut"
# resolves out of data/populist_threshold.rds exactly as it does there.
if (!exists("CELL_SAMPLES")) {
  CELL_SAMPLES <- list(
    full = list(
      label = "All scored elections",
      score_gap_min = -Inf, illiberal_cutoff = -Inf, other_cutoff_max = Inf
    ),
    popucut = list(
      label = "One side illiberal, the other not (PopuList-calibrated)",
      score_gap_min = -Inf,
      illiberal_cutoff = "popucut", other_cutoff_max = "popucut"
    )
  )
}

# Bounded to the five decades 18 and 19 use, for the same reason: it is
# exactly five decades and it is where V-Party coverage is not thin.
if (!exists("CELL_YEAR_MIN")) CELL_YEAR_MIN <- VPARTY_YEAR_MIN
if (!exists("CELL_YEAR_MAX")) CELL_YEAR_MAX <- VPARTY_YEAR_MAX

# Minimum USABLE observations for a given outcome in a given cell -- not cell
# size. Coverage differs sharply by outcome (the Ginis are missing for whole
# country-decades), so gating on cell size would hand rdrobust a handful of
# points and draw a number around whatever came back. This is the mistake
# caught in 17, where ep_galtan had 0-21 usable inside cells of 66-236.
if (!exists("CELL_MIN_N")) CELL_MIN_N <- 40

build_suffix <- paste0(
  if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "" else "_exclyr",
  if (PLACEBO_PRE_WINDOW) "_pre" else ""
)

out_dir <- sweep_dir(sprintf(
  "cell_rdd_%s%s", INSTRUMENT_LABELS[[CELL_INSTRUMENT]], build_suffix
))

build_path <- function(n) {
  file.path(
    data_dir, "rdd_build",
    sprintf("rdd_%s_w%d%s.rds", CELL_INSTRUMENT, n, build_suffix)
  )
}

missing_builds <- CELL_WINDOWS[!file.exists(vapply(CELL_WINDOWS, build_path, character(1)))]
if (length(missing_builds) > 0) {
  stop(
    "No build for window(s) ", paste(missing_builds, collapse = ", "),
    ". Run 11_build_rdd_data.R at those windows first (14_window_sweep.R ",
    "does it as a side effect).",
    call. = FALSE
  )
}

# ---- one window x sample -----------------------------------------------------

# Cut a build to one sample, using the same machinery as 12_rdd_analysis.R.
# Thresholds are resolved against the FULL build before either filter bites,
# so the two axes stay independent and a quantile spec means the same thing
# here as it does there.
apply_sample <- function(d, spec) {
  gap <- resolve_threshold(
    parse_threshold(spec$score_gap_min, "score_gap_min"),
    d$score_gap_z, "score_gap_min"
  )
  ill <- resolve_threshold(
    parse_threshold(spec$illiberal_cutoff, "illiberal_cutoff"),
    d$illiberal_score, "illiberal_cutoff"
  )
  oth <- resolve_threshold(
    parse_threshold(spec$other_cutoff_max, "other_cutoff_max", none_value = Inf),
    d$other_score, "other_cutoff_max"
  )
  d <- apply_threshold(d, "score_gap_z", gap, "score_gap_min", op = ">=")
  d <- apply_threshold(d, "illiberal_score", ill, "illiberal_cutoff", op = ">")
  d <- apply_threshold(d, "other_score", oth, "other_cutoff_max", op = "<=")
  attr(d, "resolved") <- c(
    score_gap_min = gap$absolute,
    illiberal_cutoff = ill$absolute,
    other_cutoff_max = oth$absolute
  )
  d
}

cell_keys <- expand_grid(
  group = c("OECD", "Non-OECD"),
  decade = paste0(seq((CELL_YEAR_MIN %/% 10) * 10, (CELL_YEAR_MAX %/% 10) * 10, 10), "s")
)

outcome_vars <- ALL_OUTCOME_VARS

# ---- run ---------------------------------------------------------------------

results <- list()
counts <- list()

for (samp in names(CELL_SAMPLES)) {
  spec <- CELL_SAMPLES[[samp]]
  for (n in CELL_WINDOWS) {
    d <- readRDS(build_path(n))
    cat(sprintf("\n[%s | w=%d]\n", samp, n))
    d <- apply_sample(d, spec)
    res <- attr(d, "resolved")

    d <- d |>
      filter(election_year >= CELL_YEAR_MIN, election_year <= CELL_YEAR_MAX) |>
      mutate(
        group = oecd_group(country_text_id),
        decade = decade_label(election_year)
      )

    for (k in seq_len(nrow(cell_keys))) {
      g <- cell_keys$group[k]
      dd0 <- cell_keys$decade[k]
      cell <- d[d$group == g & d$decade == dd0, , drop = FALSE]

      counts[[length(counts) + 1]] <- tibble(
        sample = samp, window = n, group = g, decade = dd0,
        n_cell = nrow(cell),
        # Recorded so "no fuzzy arm" is a fact on the page rather than a
        # claim in a comment.
        n_treated_ert = sum(cell$backsliding_Nyr, na.rm = TRUE),
        n_treated_union = sum(cell$backsliding_union_Nyr, na.rm = TRUE)
      )

      for (v in outcome_vars) {
        if (!v %in% names(cell)) next
        nu <- sum(!is.na(cell[[v]]))
        base <- tibble(
          sample = samp, sample_label = spec$label, window = n,
          group = g, decade = dd0,
          outcome = v, outcome_label = outcome_full_label(v),
          panel = names(OUTCOME_PANELS)[
            vapply(OUTCOME_PANELS, function(p) v %in% names(p), logical(1))
          ][1],
          n_cell = nrow(cell), n_usable = nu, fitted = nu >= CELL_MIN_N,
          score_gap_min = res[["score_gap_min"]],
          illiberal_cutoff = res[["illiberal_cutoff"]],
          other_cutoff_max = res[["other_cutoff_max"]]
        )
        fit <- if (nu >= CELL_MIN_N) {
          extract_rd(safe_rdrobust(cell[[v]], cell$running_var))
        } else {
          extract_rd(NULL)
        }
        results[[length(results) + 1]] <- bind_cols(
          base,
          fit |> rename(n = N, est = coef, ci_lo = ci_lo, ci_hi = ci_hi)
        )
      }
    }
  }
}

results <- bind_rows(results)
counts <- bind_rows(counts) |> distinct()

write_csv(results, file.path(out_dir, "cell_rdd_results.csv"))
write_csv(counts, file.path(out_dir, "cell_counts.csv"))

write_sweep_config(
  out_dir,
  fixed = list(
    instrument = CELL_INSTRUMENT,
    incl_election_year = TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR,
    placebo = PLACEBO_PRE_WINDOW,
    estimand = "reduced form only (no first stage / fuzzy: see header)",
    years = sprintf("%d-%d", CELL_YEAR_MIN, CELL_YEAR_MAX),
    cell_min_n = CELL_MIN_N
  ),
  swept = list(
    sample = names(CELL_SAMPLES),
    window = CELL_WINDOWS,
    cell = paste(cell_keys$group, cell_keys$decade),
    outcome = outcome_vars
  )
)

cat(sprintf(
  "\n%d estimates (%d fitted, %d below the %d-observation floor)\n",
  nrow(results), sum(results$fitted), sum(!results$fitted), CELL_MIN_N
))

# ---- a first look at w5 ------------------------------------------------------
#
# One colour grid per outcome family, rows = cells, columns = outcomes, one
# block per sample. Not the eventual figure -- just enough to see the shape
# before deciding how to plot this properly.

for (fam in names(OUTCOME_FAMILIES)) {
  fam_vars <- intersect(OUTCOME_FAMILIES[[fam]], outcome_vars)
  if (length(fam_vars) == 0) next

  blocks <- list()
  for (samp in names(CELL_SAMPLES)) {
    sub <- results |> filter(sample == samp, window == 5, outcome %in% fam_vars)
    if (nrow(sub) == 0 || all(!sub$fitted)) next
    row_labels <- setNames(
      vapply(seq_len(nrow(cell_keys)), function(k) {
        cc <- counts |>
          filter(sample == samp, window == 5,
                 group == cell_keys$group[k], decade == cell_keys$decade[k])
        sprintf("%s %s (N = %d)", cell_keys$group[k], cell_keys$decade[k],
                if (nrow(cc)) cc$n_cell[1] else 0L)
      }, character(1)),
      paste(cell_keys$group, cell_keys$decade)
    )
    sub <- sub |> mutate(cell = paste(group, decade))

    mat <- function(col) {
      m <- matrix(
        NA_real_, nrow = length(row_labels), ncol = length(fam_vars),
        dimnames = list(unname(row_labels), unname(outcome_full_label(fam_vars)))
      )
      for (i in seq_along(row_labels)) {
        for (j in seq_along(fam_vars)) {
          x <- sub[sub$cell == names(row_labels)[i] & sub$outcome == fam_vars[j], ]
          if (nrow(x) == 1 && isTRUE(x$fitted)) m[i, j] <- x[[col]]
        }
      }
      m
    }
    blocks[[length(blocks) + 1]] <- list(
      label = CELL_SAMPLES[[samp]]$label,
      row_axis = CELL_SAMPLES[[samp]]$label,
      col_axis = "Outcome",
      row_labels = unname(row_labels),
      col_labels = unname(outcome_full_label(fam_vars)),
      est = mat("est"), se = mat("se"), pval = mat("pval"), n = mat("n")
    )
  }
  if (length(blocks) == 0) next

  color_grid_blocks(
    blocks,
    path = file.path(out_dir, sprintf("grid_cells_w5_%s.html", fam)),
    title = sprintf(
      "Reduced-form RD by decade and OECD, w = 5: %s",
      unname(OUTCOME_FAMILY_TITLES[fam])
    ),
    subtitle = sprintf(
      "Instrument: %s | Window: [%s, election_year + 5] | %d-%d | Reduced form only",
      INSTRUMENT_DISPLAY[[CELL_INSTRUMENT]],
      if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "election_year" else "election_year + 1",
      CELL_YEAR_MIN, CELL_YEAR_MAX
    ),
    note = paste(
      SIG_FOOTNOTE,
      sprintf(
        "Blank cells have fewer than %d usable observations for that outcome.",
        CELL_MIN_N
      ),
      "No first stage or fuzzy LATE: treated elections per cell run 0-13 in",
      "the OECD and 6-46 outside it, so the fuzzy arm is not estimable. See",
      "cell_counts.csv."
    )
  )
  cat(sprintf("Saved grid_cells_w5_%s.html\n", fam))
}

message("\nCell RDD written to ", out_dir)
