# ==============================================================================
# Shared helpers for the fuzzy-RDD pipeline (scripts 11-17).
#
# Six groups of things live here:
#   1. Run-folder machinery -- every analysis run (a combination of instrument,
#      window, treatment definition and sample restrictions) gets its own
#      output/runs/<slug>/ directory, so the many "versions" of the analysis
#      stop colliding in a flat output/ distinguished only by filename suffix.
#      Also the instrument / party-score / treatment registries.
#   2. gt() table styling + econ-paper coefficient formatting, previously
#      copy-pasted between 11_build_rdd_data.R and 12_rdd_analysis.R.
#   3. Sample-restriction thresholds: one strict parser for all three accepted
#      forms, in both the floor and the ceiling direction.
#   4. The outcome registry -- which Y_ columns exist, how they are labelled,
#      and which share a plot panel.
#   5. rdrobust wrappers (safe_rdrobust / extract_rd), needed identically by
#      scripts 12, 13, 14 and 15.
#   6. RD panel figures, shared by 12_rdd_analysis.R and
#      17_party_outcomes_rdd.R.
#
# Analyses over RAW V-Party -- no election spine, no top-2 filter -- source
# scripts/vparty_helpers.R instead; it holds the file, the OECD and decade
# splits, and the Jaccard binning.
#
# Sourced, not run. Assumes tidyverse + gt + rdrobust are attached by the
# caller (same convention as scripts/plot_helpers.R).
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Run folders
# ------------------------------------------------------------------------------

source(here::here("scripts", "vdem_indices.R"))

RUNS_ROOT <- here::here("output", "runs")

# Where the one-off scripts in adhoc/ write. Kept out of output/runs/, which
# belongs to the pipeline: a run folder is named by run_slug() and describes a
# configuration of the analysis, and a _sweeps folder aggregates a set of
# those. An adhoc script is neither -- it answers a single question once, and
# filing its output under runs/ would imply a provenance it does not have.
ADHOC_ROOT <- here::here("output", "adhoc")

# Short, filesystem-safe label for an instrument variable. Keeps slugs
# readable ("instr-antiplural" rather than "instr-v2xpa_antiplural").
INSTRUMENT_LABELS <- c(
  v2xpa_antiplural = "antiplural",
  v2xpa_popul = "popul",
  ep_galtan = "galtan",
  v2pariglef_neg = "econleft",
  v2paanteli = "anteli"
)

# Human-readable names for plot titles / table headers. v2pariglef is NEGATED
# (see 11_build_rdd_data.R Step 2) so that economic LEFT scores high, keeping
# "higher = the side we expect to erode democracy" consistent across every
# instrument; say so explicitly everywhere it is displayed rather than leaving
# the orientation implicit.
INSTRUMENT_DISPLAY <- c(
  v2xpa_antiplural = "Anti-pluralism (v2xpa_antiplural)",
  v2xpa_popul = "Populism (v2xpa_popul)",
  ep_galtan = "GAL-TAN (ep_galtan)",
  v2pariglef_neg = "Economic left (negated v2pariglef)",
  v2paanteli = "Anti-elitism (v2paanteli)"
)

# The candidate instruments, in canonical order. Names must match
# INSTRUMENT_LABELS / INSTRUMENT_DISPLAY above.
ILLIBERALISM_VARS <- names(INSTRUMENT_LABELS)

# Party-level scores CARRIED IN THE BUILD, which is a strictly larger set than
# the scores usable AS THE INSTRUMENT.
#
# The distinction matters. 17_party_outcomes_rdd.R asks what else moves when the
# more-anti-pluralist party narrowly wins, and minority rights is one of the
# things we want to see move -- but v2paminor is not a candidate instrument, and
# adding it to ILLIBERALISM_VARS would silently enrol it in 15_alt_specs.R's
# instrument grid and in 11's `stopifnot(ILLIBERALISM_VAR %in% ILLIBERALISM_VARS)`
# as though we had proposed narrow-minority-rights-victory as a design.
#
# So: 11_build_rdd_data.R carries every PARTY_SCORE_VAR into the
# <score>__winner / <score>__loser columns and the _parties.rds companion;
# only ILLIBERALISM_VARS may be chosen as ILLIBERALISM_VAR.
#
# v2paminor is V-Party's "minority rights" item: HIGHER = more supportive of
# minority rights, i.e. it runs the OPPOSITE way from every other score here,
# where higher = the side hypothesized to erode democracy. It is deliberately
# NOT negated -- it is never an instrument, so no sign convention depends on it,
# and flipping it would make the column disagree with the V-Party codebook.
# Every display of it says which way it runs.
PARTY_SCORE_EXTRA <- c(v2paminor = "minor")
PARTY_SCORE_EXTRA_DISPLAY <- c(
  v2paminor = "Minority rights (v2paminor; higher = MORE supportive)"
)

PARTY_SCORE_VARS <- c(ILLIBERALISM_VARS, names(PARTY_SCORE_EXTRA))
PARTY_SCORE_LABELS <- c(INSTRUMENT_LABELS, PARTY_SCORE_EXTRA)
PARTY_SCORE_DISPLAY <- c(INSTRUMENT_DISPLAY, PARTY_SCORE_EXTRA_DISPLAY)

stopifnot(
  all(ILLIBERALISM_VARS %in% PARTY_SCORE_VARS),
  identical(names(PARTY_SCORE_LABELS), PARTY_SCORE_VARS),
  identical(names(PARTY_SCORE_DISPLAY), PARTY_SCORE_VARS)
)

TREATMENT_LABELS <- c(
  backsliding_Nyr = "ert",
  backsliding_union_Nyr = "ertOrDdcg",
  polyarchy_decline = "polyarchy",
  polyarchy_declined = "polyarchyBin"
)

# The last two are the same underlying quantity at two levels of measurement:
# how far polyarchy fell over the window, and simply whether it fell at all.
TREATMENT_DISPLAY <- c(
  backsliding_Nyr = "ERT autocratization episode starts within window",
  backsliding_union_Nyr = "ERT episode OR DDCG reversal within window",
  polyarchy_decline = "Size of V-Dem polyarchy decline over window (continuous)",
  polyarchy_declined = "Any V-Dem polyarchy decline over window (binary)"
)

# Numbers inside a slug: "0.6" -> "0p6", "-1.5" -> "neg1p5", non-finite -> "any".
fmt_slug_num <- function(x) {
  if (!is.finite(x)) {
    return("any")
  }
  gsub("-", "neg", gsub("\\.", "p", format(x, trim = TRUE)))
}

# A run is fully identified by these five fields; the slug is a deterministic
# function of them, so re-running the same configuration overwrites its own
# folder rather than creating a near-duplicate.
run_slug <- function(cfg) {
  # cfg$incl_election_year is optional: absent means the default (TRUE), so
  # callers that don't care about the window convention are unaffected. Only the
  # non-default convention gets a suffix, keeping default slugs clean -- but it
  # MUST get one, or a run under each convention would share a folder.
  incl <- cfg$incl_election_year %||% TRUE
  # Same pattern for the two axes added later. Both are absent from most cfgs
  # and both are omitted at their no-op value, so every slug written before they
  # existed is reproduced byte-identically and no run folder moves.
  #   other_cutoff_max  a ceiling on the LESS-illiberal party's score. No-op is
  #                     +Inf (nothing excluded), which fmt_slug_num renders
  #                     "any" -- so the suffix is dropped rather than written as
  #                     "_oppany", which would be noise on every existing run.
  #   placebo           TRUE = the pre-election placebo window.
  opp <- cfg$other_cutoff_max %||% Inf
  placebo <- cfg$placebo %||% FALSE
  paste0(
    "instr-",
    INSTRUMENT_LABELS[[cfg$instrument]],
    "_w",
    cfg$window,
    "_trt-",
    TREATMENT_LABELS[[cfg$treatment]],
    "_gap",
    fmt_slug_num(cfg$score_gap_min),
    "_illib",
    fmt_slug_num(cfg$illiberal_cutoff),
    if (is.finite(opp)) paste0("_opp", fmt_slug_num(opp)) else "",
    if (isTRUE(incl)) "" else "_exclyr",
    if (isTRUE(placebo)) "_pre" else ""
  )
}

# Create output/runs/<slug>/ (plus plots/) and drop a run_config.csv inside so
# a folder is self-describing without having to decode its own name.
run_dir <- function(cfg, subdirs = c("plots")) {
  slug <- run_slug(cfg)
  dir <- file.path(RUNS_ROOT, slug)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (sd in subdirs) {
    dir.create(file.path(dir, sd), recursive = TRUE, showWarnings = FALSE)
  }
  readr::write_csv(
    tibble::tibble(
      slug = slug,
      key = names(cfg),
      value = vapply(cfg, function(v) format(v, trim = TRUE), character(1))
    ),
    file.path(dir, "run_config.csv")
  )
  dir
}

# One row per run in output/runs/manifest.csv, so every version of the analysis
# is discoverable from a single index. Re-running a slug replaces its row
# rather than appending a duplicate.
append_run_manifest <- function(cfg, extra = list()) {
  dir.create(RUNS_ROOT, recursive = TRUE, showWarnings = FALSE)
  slug <- run_slug(cfg)
  row <- tibble::as_tibble(c(
    list(slug = slug),
    lapply(cfg, function(v) format(v, trim = TRUE)),
    lapply(extra, function(v) format(v, trim = TRUE)),
    list(run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  ))
  path <- file.path(RUNS_ROOT, "manifest.csv")
  if (file.exists(path)) {
    prev <- readr::read_csv(path, show_col_types = FALSE) |>
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
      dplyr::filter(slug != !!slug)
    row <- dplyr::bind_rows(prev, dplyr::mutate(row, dplyr::across(
      dplyr::everything(),
      as.character
    )))
  }
  readr::write_csv(row, path)
  invisible(path)
}

# Where a sweep across many runs (window sweep, alt-spec grid, restriction grid,
# instrument overlap) writes its aggregated comparison output.
#
# A sweep is NOT a run and must not borrow run_dir(): run_dir() writes a
# run_config.csv asserting a single value for every field, which is a lie about
# the axis a sweep is varying, and it would also collide with the folder that
# 12_rdd_analysis.R legitimately owns for that configuration.
sweep_dir <- function(name) {
  dir <- file.path(RUNS_ROOT, "_sweeps", name)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

# Output path for an adhoc script. Flat rather than a folder per script: each
# writes only a handful of files and they are prefixed by subject, so a
# directory per script would be one folder per file.
adhoc_path <- function(...) {
  dir.create(ADHOC_ROOT, recursive = TRUE, showWarnings = FALSE)
  file.path(ADHOC_ROOT, ...)
}

# Record what a sweep held FIXED and what it VARIED, so a sweep folder is
# self-describing in the way run_config.csv makes a run folder self-describing
# -- without pretending the swept axes have a single value.
write_sweep_config <- function(dir, fixed = list(), swept = list()) {
  rows <- dplyr::bind_rows(
    tibble::tibble(
      role = "fixed",
      key = names(fixed),
      value = vapply(fixed, function(v) paste(format(v, trim = TRUE), collapse = ", "), character(1))
    ),
    tibble::tibble(
      role = "swept",
      key = names(swept),
      value = vapply(swept, function(v) paste(format(v, trim = TRUE), collapse = ", "), character(1))
    )
  )
  readr::write_csv(rows, file.path(dir, "sweep_config.csv"))
  invisible(rows)
}

# ------------------------------------------------------------------------------
# 2. Tables
# ------------------------------------------------------------------------------

# font_size exists for one caller: 11_build_rdd_data.R's ERT miss table is 262
# rows x 8 columns and is set a point smaller so it fits. That used to be a
# second, shadowing copy of this whole function in script 11, which meant any
# change here silently failed to reach that table.
apply_table_style <- function(gt_tbl, font_size = 12) {
  gt_tbl |>
    gt::tab_options(
      table.font.size = gt::px(font_size),
      table.border.top.style = "solid",
      table.border.top.width = gt::px(2),
      table.border.top.color = "black",
      table.border.bottom.style = "solid",
      table.border.bottom.width = gt::px(2),
      table.border.bottom.color = "black",
      column_labels.border.top.style = "solid",
      column_labels.border.top.width = gt::px(2),
      column_labels.border.top.color = "black",
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = gt::px(1.5),
      column_labels.border.bottom.color = "black",
      table_body.hlines.style = "solid",
      table_body.hlines.width = gt::px(0.5),
      table_body.hlines.color = "#cccccc"
    )
}

# "-0.679 [-0.802, -0.555]***" -- estimate, robust 95% CI, significance. For
# plot titles, where a CI says more than a standard error: a reader can see at
# a glance whether the interval clears zero and how wide it is, which on a thin
# subsample is the thing worth knowing.
fmt_est_ci <- function(coef, lo, hi, pval, digits = 3) {
  ifelse(
    is.na(coef),
    "--",
    sprintf(
      "%.*f [%.*f, %.*f]%s",
      digits, coef, digits, lo, digits, hi, sig_stars(pval)
    )
  )
}

SIG_FOOTNOTE <- "Significance: * p<0.10, ** p<0.05, *** p<0.01. SE in parentheses."

CI_FOOTNOTE <- paste(
  "Brackets are rdrobust's robust (bias-corrected) 95% confidence interval,",
  "the same inference the tables report."
)

# The "N treated" column means something slightly different for the continuous
# treatment, so any table carrying that column says so on its face.
TREATED_FOOTNOTE <- paste(
  "N treated counts, within each estimation sample, the elections with a",
  "backsliding episode in the window; for the continuous treatment (decline in",
  "V-Dem polyarchy) it counts elections where polyarchy strictly fell."
)

sig_stars <- function(pval) {
  dplyr::case_when(
    is.na(pval) ~ "",
    pval < 0.01 ~ "***",
    pval < 0.05 ~ "**",
    pval < 0.10 ~ "*",
    TRUE ~ ""
  )
}

# Standard economics-paper convention: coefficient, SE in parentheses,
# significance stars -- e.g. "-0.045 (0.023)*".
fmt_est <- function(coef, se, pval, digits = 3) {
  dplyr::if_else(
    is.na(coef),
    "--",
    sprintf(
      paste0("%.", digits, "f (%.", digits, "f)%s"),
      coef,
      se,
      sig_stars(pval)
    )
  )
}

# Saved as HTML rather than PNG: gt's image export needs a headless-browser
# backend (webshot2/chromote) that isn't installed here. HTML also sizes to
# its content, so there's no width/clipping bookkeeping to get wrong.
save_table_html <- function(tbl, path, title, subtitle = NULL, note = SIG_FOOTNOTE) {
  gt_tbl <- tbl |>
    gt::gt() |>
    gt::tab_header(title = title, subtitle = subtitle) |>
    gt::opt_row_striping() |>
    apply_table_style()
  if (!is.null(note)) {
    gt_tbl <- gt_tbl |> gt::tab_source_note(source_note = note)
  }
  gt::gtsave(gt_tbl, path)
  cat(sprintf("Saved %s\n", path))
  invisible(gt_tbl)
}

# Diverging red -> white -> green colour ramp over a SYMMETRIC domain centred
# at zero, so a cell's hue encodes the sign and its intensity the magnitude,
# comparable across the whole grid. Returns a hex colour per value; NA values
# get white (they render as "--" text anyway).
diverging_fill <- function(values, max_abs = NULL) {
  finite <- values[is.finite(values)]
  if (is.null(max_abs)) {
    max_abs <- if (length(finite)) max(abs(finite)) else 1
  }
  if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  # green = positive, red = negative
  ramp <- grDevices::colorRamp(
    c("#b2182b", "#f4a582", "#ffffff", "#a6d96a", "#1a9850"),
    space = "Lab"
  )
  scaled <- (pmax(pmin(values, max_abs), -max_abs) + max_abs) / (2 * max_abs)
  out <- rep("#ffffff", length(values))
  ok <- is.finite(scaled)
  if (any(ok)) {
    rgbmat <- ramp(scaled[ok])
    out[ok] <- grDevices::rgb(rgbmat, maxColorValue = 255)
  }
  out
}

# One grid cell: estimate with significance stars, then the standard error in
# parentheses, then the subsample size -- the same ordering a regression table
# uses, so the cell reads "coef*** (se)" before the bookkeeping.
#
# Single line rather than stacked: a literal "\n" renders as a space in HTML,
# and wrapping every cell in gt::html() just to get a <br> is not worth it.
grid_cell_text <- function(est, se, pval, n, digits = 3) {
  dplyr::if_else(
    is.na(est),
    "--",
    sprintf(
      paste0("%.", digits, "f%s (%.", digits, "f)  (N=%s)"),
      est,
      sig_stars(pval),
      se,
      trimws(format(n, trim = TRUE, big.mark = ","))
    )
  )
}

# A coefficient grid: `est`, `se`, `pval` and `n` are matrices with matching
# dimnames (rows = one restriction axis, cols = the other). Cell text is
# "coef*** (N=123)"; cell fill is the diverging ramp on the coefficient.
color_grid_table <- function(
  est,
  se,
  pval,
  n,
  path,
  title,
  subtitle = NULL,
  row_label = "R1",
  digits = 3,
  note = NULL
) {
  stopifnot(identical(dim(est), dim(pval)))
  max_abs <- suppressWarnings(max(abs(est[is.finite(est)])))
  cell_text <- matrix(
    grid_cell_text(as.vector(est), as.vector(se), as.vector(pval), as.vector(n), digits),
    nrow = nrow(est),
    dimnames = dimnames(est)
  )
  fills <- matrix(
    diverging_fill(as.vector(est), max_abs),
    nrow = nrow(est),
    dimnames = dimnames(est)
  )

  tbl <- tibble::as_tibble(cell_text) |>
    dplyr::mutate(!!row_label := rownames(est), .before = 1)

  gt_tbl <- tbl |>
    gt::gt() |>
    gt::tab_header(title = title, subtitle = subtitle) |>
    gt::cols_align(align = "center", columns = -1) |>
    apply_table_style() |>
    gt::tab_source_note(source_note = paste(
      SIG_FOOTNOTE,
      "Cells show the coefficient with significance stars, the standard error in parentheses, then the subsample size.",
      "Cell colour = coefficient magnitude (green positive, red negative), on a symmetric scale.",
      "'--' = rdrobust could not fit (too few observations near the cutoff)."
    ))
  if (!is.null(note)) {
    gt_tbl <- gt_tbl |> gt::tab_source_note(source_note = note)
  }

  for (j in seq_len(ncol(est))) {
    col <- colnames(est)[j]
    for (i in seq_len(nrow(est))) {
      gt_tbl <- gt_tbl |>
        gt::tab_style(
          style = gt::cell_fill(color = fills[i, j]),
          locations = gt::cells_body(columns = col, rows = i)
        )
    }
  }

  gt::gtsave(gt_tbl, path)
  cat(sprintf("Saved %s\n", path))
  invisible(gt_tbl)
}

# ------------------------------------------------------------------------------
# 3. Sample-restriction thresholds
#
# A restriction threshold (SCORE_GAP_MIN, ILLIBERAL_CUTOFF) can be written two
# ways, and which one is meant must be unambiguous from what was written:
#
#   0.6      ABSOLUTE  -- a value on the variable's own scale
#   "q50"    QUANTILE  -- a percentile of the variable's distribution
#   -Inf     NONE      -- no restriction
#
# Both forms are needed because the five candidate instruments are NOT on a
# common scale. illiberal_score ranges observed in this project:
#
#   v2xpa_antiplural   0.02 .. 1.00   (a [0,1] index)      median 0.57
#   v2xpa_popul        0.06 .. 0.99   (a [0,1] index)      median 0.51
#   v2pariglef_neg    -1.92 .. 3.81   (expert scale)       median 0.82
#   v2paanteli        -2.37 .. 4.41   (expert scale)       median 0.33
#   ep_galtan          4.50 .. 9.41   (CHES scale)         median 7.00
#
# So a single absolute number means five different things: 0.6 is just above the
# median for anti-pluralism, BELOW the median for economic left, and below the
# entire range of ep_galtan -- where it is a silent no-op that keeps 100% of
# rows while the output still says "restricted". "q50" instead means the same
# THING everywhere ("the top half on whichever scale this instrument uses"),
# which is what makes 15_alt_specs.R a comparison of instruments rather than of
# arbitrary numbers.
#
# Anything not matching one of the three forms is a hard ERROR, deliberately.
# An earlier version fell back to as.numeric(), so a plausible typo -- "Q50",
# "p50", "q 50" -- became NA, is.finite(NA) was FALSE, and the restriction was
# skipped without a word: the run reported itself as restricted and was not.
# ------------------------------------------------------------------------------

# Threshold specs whose VALUE lives in a file another script computed, rather
# than being typed or derived from the data at hand.
#
# There is one today: the PopuList-calibrated cut that separates "illiberal"
# from "not", written by 01g_populist_threshold.R. It is a spec rather than a
# number typed into a run because the number is derived -- re-running the
# calibration on a new PopuList release should move every run that uses it,
# and a hardcoded 0.6535 would not.
#
# Each entry names the file, the field to read, and a sentence for the run's
# restriction label. Resolution happens in resolve_threshold(), not here:
# parse_threshold()'s contract is that it classifies a spec WITHOUT touching
# data, and reading an RDS would break that.
EXTERNAL_THRESHOLDS <- list(
  popucut = list(
    file = file.path("data", "populist_threshold.rds"),
    field = "threshold_abs",
    source_script = "01g_populist_threshold.R",
    describe = function(x) sprintf(
      "PopuList-calibrated illiberality cut (%s criterion, absolute transfer; fitted on %d party-years, AUC %.3f)",
      x$method, x$n, x$auc
    )
  ),
  popucut_pct = list(
    file = file.path("data", "populist_threshold.rds"),
    field = "threshold_pct",
    source_script = "01g_populist_threshold.R",
    describe = function(x) sprintf(
      "PopuList-calibrated illiberality cut (%s criterion, %.1fth-percentile transfer; fitted on %d party-years, AUC %.3f)",
      x$method, x$percentile, x$n, x$auc
    )
  )
)

# Classify a threshold spec WITHOUT looking at data. Returns a list with
# $kind ("none" | "absolute" | "quantile"), $quantile, $absolute and $spec (the
# spec as written, for the record).
#
# `none_value` is the absolute value that "no restriction" resolves to, and it
# depends on which DIRECTION the threshold cuts:
#   -Inf  for a FLOOR   (op ">=" / ">"), the original and still the default --
#         nothing is below -Inf, so nothing is excluded
#   +Inf  for a CEILING (op "<=" / "<"), used by OTHER_CUTOFF_MAX -- nothing is
#         above +Inf, so nothing is excluded
# Getting this wrong is silent and total: a ceiling whose no-op resolved to -Inf
# would keep zero rows while reporting itself as unrestricted. It is a parameter
# rather than a second parser so the three accepted forms, the error messages
# and the strictness stay defined exactly once.
parse_threshold <- function(spec, name = "threshold", none_value = -Inf) {
  stopifnot(length(none_value) == 1, is.numeric(none_value), !is.na(none_value))
  none <- function(label) {
    list(
      kind = "none", quantile = NA_real_, absolute = none_value,
      none_value = none_value, spec = label
    )
  }
  if (is.null(spec)) {
    return(none("none"))
  }
  if (length(spec) != 1) {
    stop(
      name, " must be a single value, got length ", length(spec), ".",
      call. = FALSE
    )
  }
  if (is.numeric(spec)) {
    if (is.na(spec)) {
      stop(
        name, " is NA. Use ", format(none_value, trim = TRUE),
        " for no restriction; NA is too easy to produce by accident to ",
        "accept as 'no restriction'.",
        call. = FALSE
      )
    }
    if (!is.finite(spec)) {
      return(none(format(spec, trim = TRUE)))
    }
    return(list(
      kind = "absolute", quantile = NA_real_,
      absolute = as.numeric(spec), none_value = none_value,
      spec = format(spec, trim = TRUE)
    ))
  }
  if (is.character(spec)) {
    # A named external threshold, resolved from file in resolve_threshold().
    if (spec %in% names(EXTERNAL_THRESHOLDS)) {
      return(list(
        kind = "external", quantile = NA_real_, absolute = NA_real_,
        none_value = none_value, external = spec, spec = spec
      ))
    }
    # Lowercase q, then a number: "q50", "q97.5". Nothing else.
    if (grepl("^q[0-9]+(\\.[0-9]+)?$", spec)) {
      p <- as.numeric(sub("^q", "", spec))
      if (p < 0 || p > 100) {
        stop(
          name, ' = "', spec, '" is out of range: a percentile must be 0-100.',
          call. = FALSE
        )
      }
      return(list(
        kind = "quantile", quantile = p, absolute = NA_real_,
        none_value = none_value, spec = spec
      ))
    }
    stop(
      name, ' = "', spec, '" is not a recognized threshold.\n',
      "  Accepted forms:\n",
      "    a number  -- absolute cutoff on the variable's own scale, e.g. 0.6\n",
      '    "qNN"     -- percentile of its distribution, e.g. "q50", "q97.5"\n',
      "                 (lowercase q, no space, NN between 0 and 100)\n",
      "    a named external threshold, resolved from a file:\n",
      "                 ", paste(names(EXTERNAL_THRESHOLDS), collapse = ", "), "\n",
      "    ", format(none_value, trim = TRUE),
      "      -- no restriction (this threshold's no-op sentinel)\n",
      "  This is an error rather than a fallback on purpose: silently ignoring an\n",
      "  unparseable threshold would produce an unrestricted run labelled as a\n",
      "  restricted one.",
      call. = FALSE
    )
  }
  stop(
    name, ' must be a number or a "qNN" string, got ', class(spec)[1], ".",
    call. = FALSE
  )
}

# Turn a parsed spec into a concrete absolute threshold, using `values`.
#
# Call this on the FULL sample, before any other restriction is applied.
# Resolving a quantile after another filter has bitten would make the two
# restrictions interact -- changing SCORE_GAP_MIN would silently move an
# ILLIBERAL_CUTOFF of "q50" as well -- so the axes must be resolved up front to
# stay independent.
resolve_threshold <- function(parsed, values, name = "threshold") {
  if (parsed$kind == "none") {
    parsed$absolute <- parsed$none_value %||% -Inf
    parsed$label <- "no restriction"
    return(parsed)
  }
  if (parsed$kind == "external") {
    ext <- EXTERNAL_THRESHOLDS[[parsed$external]]
    path <- here::here(ext$file)
    if (!file.exists(path)) {
      stop(
        name, ' = "', parsed$spec, '" needs ', ext$file,
        ", which does not exist. Run scripts/", ext$source_script,
        " first -- it downloads The PopuList and calibrates the cut.",
        call. = FALSE
      )
    }
    x <- readRDS(path)
    if (is.null(x[[ext$field]]) || !is.finite(x[[ext$field]])) {
      stop(
        name, ' = "', parsed$spec, '": ', ext$file, " has no usable '",
        ext$field, "'. Re-run scripts/", ext$source_script, ".",
        call. = FALSE
      )
    }
    parsed$absolute <- as.numeric(x[[ext$field]])
    parsed$label <- sprintf("%.4g -- %s", parsed$absolute, ext$describe(x))
    return(parsed)
  }
  if (parsed$kind == "quantile") {
    ok <- values[!is.na(values)]
    if (length(ok) == 0) {
      stop(
        name, ': cannot resolve "', parsed$spec,
        '" -- the variable it restricts is missing for every row.',
        call. = FALSE
      )
    }
    parsed$absolute <- unname(quantile(ok, probs = parsed$quantile / 100))
    parsed$label <- sprintf(
      "%s = the %sth percentile, which is %.4g on this instrument's scale",
      parsed$spec, format(parsed$quantile, trim = TRUE), parsed$absolute
    )
    parsed$n_resolved_on <- length(ok)
  } else {
    parsed$label <- sprintf("%.4g (absolute)", parsed$absolute)
  }
  parsed
}

# Apply a resolved threshold and say plainly what happened: the spec as written,
# how it was interpreted, the absolute value used, and how many rows it cost.
#
# `op` is a FLOOR (">=" or ">") or a CEILING ("<=" or "<"). The ceiling forms
# exist for OTHER_CUTOFF_MAX, which bounds the LESS-illiberal party's score from
# above; paired with a floor on the more-illiberal party's score they express
# "one side is illiberal and the other is not". The `op` a caller passes must
# agree in direction with the `none_value` it parsed with -- see parse_threshold.
apply_threshold <- function(data, var, parsed, name,
                            op = c(">=", ">", "<=", "<")) {
  op <- match.arg(op)
  if (parsed$kind == "none") {
    cat(sprintf("  %-16s no restriction (%d rows kept)\n", name, nrow(data)))
    return(data)
  }
  # A floor whose no-op is +Inf (or a ceiling whose no-op is -Inf) would drop
  # every row the moment the spec went back to "no restriction", and the run
  # would still describe itself as unrestricted. Catch the mismatch here rather
  # than letting it surface as an inexplicably empty sample.
  expected_none <- if (op %in% c(">=", ">")) -Inf else Inf
  none_value <- parsed$none_value %||% -Inf
  if (!identical(none_value, expected_none)) {
    stop(
      name, ": direction mismatch. op = \"", op, "\" is a ",
      if (op %in% c(">=", ">")) "floor" else "ceiling",
      ", whose no-restriction sentinel is ", format(expected_none, trim = TRUE),
      ", but the spec was parsed with none_value = ",
      format(none_value, trim = TRUE),
      ". Pass none_value = ", format(expected_none, trim = TRUE),
      " to parse_threshold().",
      call. = FALSE
    )
  }
  n_before <- nrow(data)
  keep <- !is.na(data[[var]]) & switch(op,
    ">=" = data[[var]] >= parsed$absolute,
    ">" = data[[var]] > parsed$absolute,
    "<=" = data[[var]] <= parsed$absolute,
    "<" = data[[var]] < parsed$absolute
  )
  out <- data[keep, , drop = FALSE]
  cat(sprintf(
    "  %-16s %s  ->  keep %s %s %.4g:  %d -> %d rows (%d dropped)\n",
    name, parsed$label, var, op, parsed$absolute,
    n_before, nrow(out), n_before - nrow(out)
  ))
  if (nrow(out) == 0) {
    warning(
      name, " = ", parsed$spec, " (", parsed$absolute,
      ") removed every row. Check whether the threshold is on the right scale ",
      "for this instrument -- a quantile spec such as \"q50\" adapts to the ",
      "scale, an absolute number does not.",
      call. = FALSE
    )
  }
  out
}

# parse + resolve in one step, returning just the NUMBER.
#
# Callers that need a threshold as a number rather than a spec object -- a
# folder name via fmt_slug_num(), a cfg that run_slug() will turn back into a
# folder name, a sentence in a subtitle -- all want the same two calls in the
# same order, and getting the order wrong (resolving after a filter has bitten)
# is the mistake resolve_threshold() warns about. One entry point so there is
# one place to get it right.
resolve_threshold_abs <- function(spec, name, values, none_value = -Inf) {
  resolve_threshold(
    parse_threshold(spec, name, none_value = none_value), values, name
  )$absolute
}

# One English sentence naming every ACTIVE restriction, from the three resolved
# absolute thresholds. Every table subtitle and figure subtitle that describes a
# sample goes through here, so the sample restriction travels with the output
# instead of living only in a folder name.
#
# It must name every active axis. This list was written when there were two axes
# and was not extended when other_cutoff_max was added, so the PopuList sweeps
# described themselves as the one-sided "illiberal_score > 0.6535" (603
# elections) while actually estimating the two-sided pair condition (411). For
# that specification the pair condition IS the design, so the label was not a
# cosmetic slip -- it named a different design from the one being run.
#
# When the floor and the ceiling sit at the same value the pair reads as one
# statement rather than two, and saying so is clearer than making the reader
# notice that the two numbers coincide.
restriction_sentence <- function(score_gap_min = -Inf,
                                 illiberal_cutoff = -Inf,
                                 other_cutoff_max = Inf,
                                 none = "all elections") {
  # c() drops NULLs, so an inactive axis contributes nothing rather than an
  # empty string that paste() would render as a stray comma.
  parts <- c(
    if (is.finite(score_gap_min)) {
      sprintf("score_gap_z >= %s", format(score_gap_min, trim = TRUE))
    },
    if (is.finite(illiberal_cutoff) && is.finite(other_cutoff_max) &&
          isTRUE(all.equal(illiberal_cutoff, other_cutoff_max))) {
      sprintf(
        "one top-2 party illiberal and the other not (illiberal_score > %s >= other_score)",
        format(illiberal_cutoff, trim = TRUE)
      )
    } else {
      c(
        if (is.finite(illiberal_cutoff)) {
          sprintf("illiberal_score > %s", format(illiberal_cutoff, trim = TRUE))
        },
        if (is.finite(other_cutoff_max)) {
          sprintf("other_score <= %s", format(other_cutoff_max, trim = TRUE))
        }
      )
    }
  )
  if (length(parts) == 0) none else paste(parts, collapse = ", ")
}


# Several coefficient grids stacked into ONE table, one section per grid.
#
# The awkward part is that each grid has a different column axis -- 5 quantile
# levels for a continuous restriction, 2 for a yes/no one, 3 for election type
# -- so there is no single set of column headers that means the same thing in
# every section. Rather than pad to a 15-column union (mostly blank) or throw
# away the grid shape by going long, the columns are POSITIONAL (Level 1..k) and
# each section opens with its own label row naming what its columns are. The
# section header states the row axis and the column axis.
#
# `blocks` is a list of lists, each with: label, row_axis, col_axis, row_labels,
# col_labels, and the est/se/pval/n matrices.
#
# The colour scale is shared across every section, on one symmetric domain, so a
# hue means the same magnitude everywhere in the table.
color_grid_blocks <- function(
  blocks, path, title, subtitle = NULL, digits = 3, note = NULL
) {
  n_col_max <- max(vapply(blocks, function(b) length(b$col_labels), integer(1)))
  col_ids <- paste0("C", seq_len(n_col_max))

  all_est <- unlist(lapply(blocks, function(b) as.vector(b$est)))
  max_abs <- suppressWarnings(max(abs(all_est[is.finite(all_est)])))

  rows <- list()
  fills <- list() # (row index, column id) -> colour, data cells only
  label_row_ids <- integer(0)
  r <- 0L

  pad <- function(x) c(x, rep("", n_col_max - length(x)))

  for (b in blocks) {
    section <- sprintf(
      "%s  (rows)   \u00d7   %s  (columns)", b$row_axis, b$col_axis
    )

    # The section's own column-label row.
    r <- r + 1L
    label_row_ids <- c(label_row_ids, r)
    rows[[length(rows) + 1L]] <- c(
      list(section = section, level = "columns:"),
      setNames(as.list(pad(b$col_labels)), col_ids)
    )

    for (i in seq_along(b$row_labels)) {
      r <- r + 1L
      txt <- grid_cell_text(
        b$est[i, ], b$se[i, ], b$pval[i, ], b$n[i, ], digits
      )
      rows[[length(rows) + 1L]] <- c(
        list(section = section, level = b$row_labels[i]),
        setNames(as.list(pad(txt)), col_ids)
      )
      cols_here <- diverging_fill(b$est[i, ], max_abs)
      for (j in seq_along(cols_here)) {
        fills[[length(fills) + 1L]] <- list(
          row = r, col = col_ids[j], fill = cols_here[j],
          na = is.na(b$est[i, j])
        )
      }
    }
  }

  tbl <- dplyr::bind_rows(lapply(rows, tibble::as_tibble))

  gt_tbl <- tbl |>
    gt::gt(groupname_col = "section") |>
    gt::tab_header(title = title, subtitle = subtitle) |>
    gt::cols_label(.list = setNames(as.list(rep("", n_col_max)), col_ids)) |>
    gt::cols_label(level = "") |>
    gt::cols_align(align = "center", columns = dplyr::all_of(col_ids)) |>
    gt::cols_align(align = "left", columns = "level") |>
    apply_table_style() |>
    gt::tab_options(
      row_group.font.weight = "bold",
      row_group.background.color = "#f0f0f0",
      row_group.border.top.width = gt::px(2),
      row_group.border.top.color = "black",
      row_group.border.bottom.width = gt::px(1),
      row_group.border.bottom.color = "#888888"
    ) |>
    gt::tab_source_note(source_note = paste(
      SIG_FOOTNOTE,
      "Cells show the coefficient with significance stars, the standard error in parentheses, then the subsample size.",
      "Cell colour = coefficient magnitude (green positive, red negative), on one symmetric scale shared by every section.",
      "'--' = rdrobust could not fit (too few observations near the cutoff)."
    ))
  if (!is.null(note)) {
    gt_tbl <- gt_tbl |> gt::tab_source_note(source_note = note)
  }

  # Column-label rows read as sub-headers, not data.
  gt_tbl <- gt_tbl |>
    gt::tab_style(
      style = list(
        gt::cell_text(style = "italic", color = "#333333", size = gt::px(11)),
        gt::cell_fill(color = "#fafafa")
      ),
      locations = gt::cells_body(rows = label_row_ids)
    )

  for (f in fills) {
    if (f$na) {
      next
    }
    gt_tbl <- gt_tbl |>
      gt::tab_style(
        style = gt::cell_fill(color = f$fill),
        locations = gt::cells_body(columns = f$col, rows = f$row)
      )
  }

  gt::gtsave(gt_tbl, path)
  cat(sprintf("Saved %s\n", path))
  invisible(gt_tbl)
}

# ------------------------------------------------------------------------------
# 4. Outcome panels and series palette
#
# Shared by 12_rdd_analysis.R (which plots them) and 14/15 (which label pooled
# results with them), so the grouping and the labels are defined once.
# ------------------------------------------------------------------------------

OUTCOME_PANELS <- list(
  growth = c(
    Y_gdp_growth = "PWT",
    Y_gdp_growth_wb = "World Bank",
    Y_gdp_growth_imf = "IMF WEO"
  ),
  inflation = c(Y_inflation = "CPI, cumulative log"),
  unemployment = c(Y_unemployment = "Unemployment rate"),
  trade = c(Y_trade_pct_gdp = "Trade / GDP"),
  top10 = c(Y_top10_share = "Top-10% income share"),
  gini = c(Y_gini_disp = "Disposable", Y_gini_mkt = "Market"),
  fiscal = c(Y_debt = "Debt / GDP", Y_deficit = "Deficit / GDP"),
  institutions = c(
    Y_checks_balances = "Checks & balances",
    Y_jucon = "Judicial constraints",
    Y_legcon = "Legislative constraints"
  ),
  exec_power = c(
    Y_hos_power = "HOS",
    Y_hog_power = "HOG"
  ),
  exec_power_alt = c(
    Y_hos_power_vdem = "HOS (v2ex_hosw)",
    Y_hog_power_vdem = "HOG (v2ex_hogw)"
  ),
  # ---- V-Dem democracy indices --------------------------------------------
  #
  # Grouped on V-DEM'S OWN taxonomy: the five high-level indices, then the
  # mid-level indices that aggregate into each of them. A principled grouping
  # beats a convenience one here, and it keeps every panel at or under five
  # series, which is all SERIES_COLORS carries.
  #
  # Y_jucon and Y_legcon are not new -- they were already outcomes in the
  # `institutions` panel above, where they sit alongside checks_balances (their
  # own mean). They appear ONCE in OUTCOME_PANELS, under `institutions`, so the
  # results table keeps one row per outcome; vdem_liberal therefore shows only
  # the rule-of-law index that completes V-Dem's liberal component, and the
  # panel title says where the other two are.
  #
  # Direction: every index runs HIGHER = MORE DEMOCRATIC, so a NEGATIVE RD
  # estimate is the backsliding sign -- the opposite of the economic outcomes.
  vdem_high = c(
    Y_polyarchy = "Electoral (polyarchy)",
    Y_libdem = "Liberal",
    Y_partipdem = "Participatory",
    Y_delibdem = "Deliberative",
    Y_egaldem = "Egalitarian"
  ),
  vdem_electoral = c(
    Y_elecoff = "Elected officials",
    Y_frefair = "Clean elections",
    Y_frassoc = "Association",
    Y_suffr = "Suffrage",
    Y_freexp = "Expression & alt. info"
  ),
  vdem_liberal = c(
    Y_cl_rol = "Rule of law"
  ),
  vdem_particip = c(
    Y_cspart = "Civil society",
    Y_dd = "Direct democracy",
    Y_locelec = "Local elections",
    Y_regelec = "Regional elections"
  ),
  vdem_egal_delib = c(
    Y_delib = "Deliberation",
    Y_eqprotec = "Equal protection",
    Y_eqaccess = "Equal access",
    Y_eqdr = "Equal distribution"
  )
)

PANEL_TITLES <- c(
  growth = "GDP per capita growth",
  inflation = "Inflation",
  unemployment = "Unemployment",
  trade = "Trade openness",
  top10 = "Top-10% income share",
  gini = "Income inequality (Gini)",
  fiscal = "Public debt and deficit",
  institutions = "Executive constraints",
  exec_power = "Executive power (ET)",
  exec_power_alt = "Executive power (V-Dem)",
  vdem_high = "V-Dem high-level democracy indices",
  vdem_electoral = "V-Dem mid-level: electoral component",
  vdem_liberal = "V-Dem mid-level: liberal component (see also Executive constraints)",
  vdem_particip = "V-Dem mid-level: participatory component",
  vdem_egal_delib = "V-Dem mid-level: deliberative and egalitarian components"
)

# The y axis carries the UNIT; the panel title carries the subject. Repeating
# the title down the y axis (an earlier version did) wastes the axis and tells
# the reader nothing about what a value of 0.1 means.
PANEL_YLABS <- c(
  growth = "Cumulative log change",
  inflation = "Cumulative log change",
  unemployment = "Change, pp of labour force",
  trade = "Change, pp of GDP",
  top10 = "Change, pp of national income",
  gini = "Change, Gini points",
  fiscal = "Change, pp of GDP",
  institutions = "Change, index (0-1)",
  exec_power = "Change, index (0-1)",
  exec_power_alt = "Change, index (0-1)",
  vdem_high = "Change, index (0-1)",
  vdem_electoral = "Change, index (0-1)",
  vdem_liberal = "Change, index (0-1)",
  vdem_particip = "Change, index (0-1)",
  vdem_egal_delib = "Change, index (0-1)"
)


# Okabe-Ito blue / vermillion / green / reddish-purple / black. Colourblind-safe:
# worst pairwise separation across all five is dE 19.0 under deuteranopia and
# 23.7 under protanopia, and every one clears 3:1 contrast against a white panel
# (5.19, 3.87, 3.42, 3.06, 21.0). Series are ALSO distinguished by point shape,
# so identity never rests on colour alone.
#
# Assigned in FIXED order and never cycled, so extending the palette from three
# to five leaves every existing one-, two- and three-series panel byte-identical.
#
# Black is the fifth rather than Okabe-Ito's orange (#E69F00) on purpose: orange
# only reaches 2.25:1 against white, and every darkened variant that fixes the
# contrast collapses onto vermillion under deuteranopia (dE 1.1-13.5). Black is
# itself an Okabe-Ito colour, and the only cutoff/reference marks it could be
# confused with are grey35 dashed, never solid.
#
# Shapes: base R has only four solid glyphs, so the fifth is the asterisk. It is
# the weakest of the five at size 1.6, which is why it sits alongside the most
# distinctive colour.
SERIES_COLORS <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#000000")
SERIES_SHAPES <- c(16, 17, 15, 18, 8)


# Coarser groupings than OUTCOME_PANELS, for figures that need one sheet per
# family of outcomes rather than one per plot panel. OUTCOME_PANELS is about
# which series share an axis (same units); this is about which outcomes a reader
# wants to look at together.
OUTCOME_FAMILIES <- list(
  economic = c(
    "Y_gdp_growth", "Y_gdp_growth_wb", "Y_gdp_growth_imf",
    "Y_inflation", "Y_unemployment", "Y_trade_pct_gdp"
  ),
  inequality = c("Y_top10_share", "Y_gini_disp", "Y_gini_mkt"),
  fiscal = c("Y_debt", "Y_deficit"),
  institutions = c(
    "Y_checks_balances", "Y_jucon", "Y_legcon",
    "Y_hos_power", "Y_hog_power", "Y_hos_power_vdem", "Y_hog_power_vdem"
  ),
  # Their own families rather than folded into `institutions`: these are the
  # mechanism outcomes, a different question from executive constraints.
  #
  # Split high from mid on the same line the ask draws. One combined family
  # would be 19 facets on a single window-sweep sheet -- seven rows deep, and
  # the five headline indices, which are the ones actually being asked about,
  # would be buried among fourteen components.
  #
  # Y_jucon and Y_legcon are excluded from vdem_mid because they already appear
  # in the `institutions` family; OUTCOME_FAMILIES may repeat an outcome across
  # families, but showing the same two panels twice on adjacent sheets is
  # clutter, not emphasis.
  vdem_high = names(OUTCOME_PANELS$vdem_high),
  vdem_mid = setdiff(
    names(VDEM_INDEX_VARS),
    c(names(OUTCOME_PANELS$vdem_high), "Y_jucon", "Y_legcon")
  )
)

OUTCOME_FAMILY_TITLES <- c(
  economic = "Economic outcomes",
  inequality = "Inequality",
  fiscal = "Public finances",
  institutions = "Institutions and executive power",
  vdem_high = "V-Dem high-level democracy indices",
  vdem_mid = "V-Dem mid-level democracy indices"
)

# Flat outcome list, in panel order (NOT alphabetical), so facets and table rows
# read in the same sequence as the figures.
ALL_OUTCOME_VARS <- unlist(lapply(OUTCOME_PANELS, names), use.names = FALSE)

# A self-describing label for one outcome. Inside a multi-series panel the
# series label alone ("PWT", "Disposable") means nothing out of context, so
# qualify it with the panel title; a lone series just takes the panel title.
outcome_full_label <- function(var) {
  vapply(var, function(v) {
    pn <- names(OUTCOME_PANELS)[
      vapply(OUTCOME_PANELS, function(p) v %in% names(p), logical(1))
    ]
    if (length(pn) == 0) {
      return(v)
    }
    pn <- pn[1]
    panel <- OUTCOME_PANELS[[pn]]
    if (length(panel) > 1) {
      sprintf("%s: %s", unname(PANEL_TITLES[pn]), unname(panel[v]))
    } else {
      unname(PANEL_TITLES[pn])
    }
  }, character(1), USE.NAMES = FALSE)
}

# ------------------------------------------------------------------------------
# 5. rdrobust wrappers
# ------------------------------------------------------------------------------

# bwselect = "mserd" (Calonico, Cattaneo & Titiunik 2014's MSE-optimal,
# common-bandwidth selector) is rdrobust()'s own default -- named explicitly so
# the bandwidth-selection method is documented in code, and so plotting code can
# request the identical procedure when it re-fits to get a bandwidth.
RD_BWSELECT <- "mserd"

# Minimum non-missing observations before we even attempt a fit. Below this
# rdrobust either errors or returns something not worth reporting.
RD_MIN_OBS <- 20

safe_rdrobust <- function(y, x, fuzzy = NULL) {
  tryCatch(
    {
      # is.null(fuzzy) | !is.na(fuzzy) would silently collapse to logical(0)
      # when fuzzy is NULL (!is.na(NULL) is zero-length), zeroing out `keep`
      # via vectorized `&` -- branch explicitly instead.
      fuzzy_ok <- if (is.null(fuzzy)) rep(TRUE, length(y)) else !is.na(fuzzy)
      keep <- !is.na(y) & !is.na(x) & fuzzy_ok
      if (sum(keep) < RD_MIN_OBS) {
        return(NULL)
      }
      if (is.null(fuzzy)) {
        rdrobust::rdrobust(y = y[keep], x = x[keep], bwselect = RD_BWSELECT)
      } else {
        # A degenerate treatment (no variation, or none within the bandwidth)
        # makes the Wald denominator zero; rdrobust's own error is opaque, so
        # bail out early instead.
        if (length(unique(fuzzy[keep])) < 2) {
          return(NULL)
        }
        rdrobust::rdrobust(
          y = y[keep],
          x = x[keep],
          fuzzy = fuzzy[keep],
          bwselect = RD_BWSELECT
        )
      }
    },
    # NOTE: errors only. rdrobust warns routinely (mass points in the running
    # variable, small effective samples); treating a warning as a failure here
    # would silently discard perfectly good fits.
    error = function(e) NULL
  )
}

# How many of the observations a given rdrobust call actually used are treated.
#
# Must apply the SAME non-missing mask safe_rdrobust() applies, otherwise the
# reported "N treated" would not be a subset of the reported "N" -- outcomes
# differ in missingness, so a treated count taken over the whole sample would
# exceed N for the sparser ones.
#
# For a 0/1 treatment this is just the number of 1s. For the CONTINUOUS
# treatment (polyarchy_decline, the negated change in V-Dem polyarchy over the
# window) there is no 0/1 split, but there is still a meaningful treated set:
# the elections where polyarchy actually FELL, i.e. where the decline is
# strictly positive. Exact zeros are untreated -- no change is no backsliding.
# Reporting NA here instead, as an earlier version did, threw away the one
# number that says how much of the sample carries any backsliding signal at all.
rd_n_treated <- function(y, x, treatment) {
  if (is.null(treatment)) {
    return(NA_integer_)
  }
  keep <- !is.na(y) & !is.na(x) & !is.na(treatment)
  vals <- treatment[keep]
  if (length(vals) == 0) {
    return(NA_integer_)
  }
  if (all(vals %in% c(0, 1))) {
    return(as.integer(sum(vals)))
  }
  as.integer(sum(vals > 0))
}

# Pull the robust bias-corrected coefficient/SE/p-value/bandwidth out of an
# rdrobust fit (or NA placeholders if the fit is NULL/failed).
extract_rd <- function(fit) {
  if (is.null(fit)) {
    return(tibble::tibble(
      N = NA_integer_,
      coef = NA_real_,
      se = NA_real_,
      pval = NA_real_,
      ci_lo = NA_real_,
      ci_hi = NA_real_,
      bw = NA_real_
    ))
  }
  tibble::tibble(
    N = sum(fit$N),
    coef = unname(fit$coef["Robust", 1]),
    se = unname(fit$se["Robust", 1]),
    pval = unname(fit$pv["Robust", 1]),
    # rdrobust's own robust interval rather than coef +/- 1.96*se. The two
    # agree exactly today (checked), but taking the package's number means the
    # interval stays correct if it ever changes how it forms one -- and the
    # robust interval is NOT centred on the conventional estimate, so
    # reconstructing it by hand from the wrong row is an easy mistake.
    ci_lo = unname(fit$ci["Robust", 1]),
    ci_hi = unname(fit$ci["Robust", 2]),
    bw = unname(fit$bws[1, 1])
  )
}

# ------------------------------------------------------------------------------
# 6. RD panel figures
#
# Moved here from 12_rdd_analysis.R so 17_party_outcomes_rdd.R draws the same
# figure rather than a lookalike. Nothing below reads any script-local state --
# only SERIES_COLORS, SERIES_SHAPES, RD_MIN_OBS and safe_rdrobust(), all defined
# above. The caller supplies the data, the series and the labels.
# ------------------------------------------------------------------------------

# The weighted least-squares regression rdrobust()/rdplot() use for the point
# estimate: OLS with triangular-kernel weights (1 - |x|/h, clipped at 0)
# restricted to |x| <= h, fitted separately on each side of the cutoff. This is
# not an approximation -- the weighted-OLS intercept at 0 reproduces
# rdrobust()'s own "Conventional" jump estimate exactly. The SEs are ordinary
# WLS prediction SEs, NOT CCT's bias-corrected "Robust" SEs reported in the
# table, so the ribbon visualizes the conventional fit's uncertainty rather
# than substituting for the table's formal inference.
side_fit <- function(y, x, side, bandwidth, grid) {
  w <- pmax(1 - abs(x) / bandwidth, 0)
  idx <- which(side & w > 0 & !is.na(y))
  if (length(idx) < 3) {
    return(NULL)
  }
  df <- data.frame(xx = x[idx], yy = y[idx], ww = w[idx])
  fit <- tryCatch(lm(yy ~ xx, data = df, weights = ww), error = function(e) NULL)
  if (is.null(fit)) {
    return(NULL)
  }
  pred <- predict(fit, newdata = data.frame(xx = grid), se.fit = TRUE)
  tcrit <- qt(0.975, df = fit$df.residual)
  tibble(
    xx = grid,
    yhat = pred$fit,
    ymin = pred$fit - tcrit * pred$se.fit,
    ymax = pred$fit + tcrit * pred$se.fit
  )
}

# Binned local means of y against the running variable, as rdplot() computes
# them. hide = TRUE computes the bins without rendering (rdplot() otherwise
# prints to whatever device is open, which errors in a non-interactive Rscript
# session). binselect defaults to "esmv" (mimicking-variance, evenly spaced);
# RD_BIN_SCALE multiplies its automatically-chosen bin count rather than
# hardcoding an nbins that wouldn't adapt across very different subsample sizes.
RD_BIN_SCALE <- 2

binned_means <- function(y, x, bandwidth) {
  probe <- tryCatch(
    rdplot(
      y = y, x = x, p = 1, h = bandwidth,
      kernel = "triangular", scale = RD_BIN_SCALE, hide = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(probe)) {
    return(NULL)
  }
  as_tibble(probe$vars_bins) |>
    select(xx = rdplot_mean_bin, yy = rdplot_mean_y) |>
    filter(is.finite(xx), is.finite(yy))
}

# One RD panel, one or more series sharing a single y-axis.
#
# All series in a panel share ONE bandwidth -- the MSE-optimal bandwidth of the
# first (primary) series -- so the binning and the fitted window are identical
# across series and the lines are visually comparable. Letting each series pick
# its own bandwidth would make three curves drawn over different x-ranges look
# like a substantive difference when it is a bandwidth difference.
#
# p = 1 (local linear) matches what rdrobust() actually estimates;
# rdplot()'s own default p = 4 snakes through nearly every bin with only ~15-20
# binned points per side. kernel = "triangular" matches rdrobust()'s weighting
# (rdplot() otherwise defaults to uniform).
build_panel_plot <- function(
  data, series, title, y_label,
  y_lim = NULL, show_legend = TRUE
) {
  x <- data$running_var
  vars <- names(series)

  primary_fit <- safe_rdrobust(data[[vars[1]]], x)
  h <- if (!is.null(primary_fit)) max(unname(primary_fit$bws["h", ])) else NULL
  if (is.null(h) || !is.finite(h)) {
    h <- max(abs(x), na.rm = TRUE)
  }

  grid_left <- seq(-h, 0, length.out = 100)
  grid_right <- seq(0, h, length.out = 100)

  bins <- list()
  fits <- list()
  for (i in seq_along(vars)) {
    v <- vars[i]
    y <- data[[v]]
    keep <- !is.na(y) & !is.na(x)
    if (sum(keep) < RD_MIN_OBS) {
      next
    }
    yk <- y[keep]
    xk <- x[keep]
    b <- binned_means(yk, xk, h)
    if (!is.null(b)) {
      bins[[v]] <- b |>
        filter(abs(xx) <= h) |>
        mutate(series = unname(series[v]))
    }
    # side_fit() returns NULL when a side has fewer than three usable points
    # inside the bandwidth, and mutate() on NULL errors. Piping straight into
    # mutate() therefore crashed the whole figure whenever ONE side was thin --
    # which never happened on the pooled samples this was written for, but does
    # on a sparse subsample (an ep_galtan decade x OECD cell, say, where both
    # top-2 members are scored in about 4% of elections). Tag each side only if
    # it fitted; a panel with one side is still worth drawing.
    tag_side <- function(fit, side) {
      if (is.null(fit)) NULL else mutate(fit, side = side)
    }
    f <- bind_rows(
      tag_side(side_fit(yk, xk, xk < 0, h, grid_left), "left"),
      tag_side(side_fit(yk, xk, xk >= 0, h, grid_right), "right")
    )
    if (nrow(f) > 0) {
      fits[[v]] <- f |> mutate(series = unname(series[v]))
    }
  }

  if (length(bins) == 0 && length(fits) == 0) {
    message("Skipping panel (too few observations): ", title)
    return(NULL)
  }

  bins_df <- bind_rows(bins)
  fits_df <- bind_rows(fits)
  # Fixed factor order so a series keeps its colour regardless of which
  # series happen to survive in a given subsample.
  lvls <- unname(series)
  if (nrow(bins_df)) bins_df$series <- factor(bins_df$series, levels = lvls)
  if (nrow(fits_df)) fits_df$series <- factor(fits_df$series, levels = lvls)

  n_series <- length(lvls)
  ribbon_alpha <- if (n_series >= 3) 0.09 else 0.14

  if (is.null(y_lim)) {
    bounds <- c(bins_df$yy, fits_df$ymin, fits_df$ymax)
    bounds <- bounds[is.finite(bounds)]
    if (length(bounds) > 0) {
      rng <- range(bounds)
      pad <- max(diff(rng) * 0.1, 1e-9)
      y_lim <- c(rng[1] - pad, rng[2] + pad)
    }
  }

  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35", linewidth = 0.4)

  if (nrow(fits_df)) {
    p <- p +
      geom_ribbon(
        data = fits_df,
        aes(x = xx, ymin = ymin, ymax = ymax, group = interaction(series, side), fill = series),
        alpha = ribbon_alpha, colour = NA
      ) +
      geom_line(
        data = fits_df,
        aes(x = xx, y = yhat, group = interaction(series, side), colour = series),
        linewidth = 0.7
      )
  }
  if (nrow(bins_df)) {
    p <- p +
      geom_point(
        data = bins_df,
        aes(x = xx, y = yy, colour = series, shape = series),
        size = 1.6, alpha = 0.85
      )
  }

  p +
    scale_colour_manual(values = setNames(SERIES_COLORS[seq_len(n_series)], lvls), drop = FALSE) +
    scale_fill_manual(values = setNames(SERIES_COLORS[seq_len(n_series)], lvls), drop = FALSE) +
    scale_shape_manual(values = setNames(SERIES_SHAPES[seq_len(n_series)], lvls), drop = FALSE) +
    coord_cartesian(xlim = c(-h, h), ylim = y_lim) +
    labs(
      title = title,
      x = "Running variable (illiberal - other vote/seat share, pp)",
      y = y_label,
      colour = NULL, fill = NULL, shape = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      plot.title = element_text(size = 9, face = "bold"),
      # A legend is always present for >= 2 series so identity is never carried
      # by colour alone; a single-series panel is named by its own title.
      legend.position = if (show_legend && n_series > 1) "bottom" else "none",
      legend.key.size = unit(0.35, "cm"),
      legend.margin = margin(t = -4)
    )
}
