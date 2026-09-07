# ==============================================================================
# Shared helpers for the fuzzy-RDD pipeline (scripts 11-16).
#
# Three groups of things live here:
#   1. Run-folder machinery -- every analysis run (a combination of instrument,
#      window, treatment definition and sample restrictions) gets its own
#      output/runs/<slug>/ directory, so the many "versions" of the analysis
#      stop colliding in a flat output/ distinguished only by filename suffix.
#   2. gt() table styling + econ-paper coefficient formatting, previously
#      copy-pasted between 11_build_rdd_data.R and 12_rdd_analysis.R.
#   3. rdrobust wrappers (safe_rdrobust / extract_rd), needed identically by
#      scripts 12, 13, 14 and 15.
#
# Sourced, not run. Assumes tidyverse + gt + rdrobust are attached by the
# caller (same convention as scripts/plot_helpers.R).
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Run folders
# ------------------------------------------------------------------------------

RUNS_ROOT <- here::here("output", "runs")

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

TREATMENT_LABELS <- c(
  backsliding_Nyr = "ert",
  backsliding_union_Nyr = "ertOrDdcg",
  polyarchy_decline = "polyarchy"
)

TREATMENT_DISPLAY <- c(
  backsliding_Nyr = "ERT autocratization episode starts within window",
  backsliding_union_Nyr = "ERT episode OR DDCG reversal within window",
  polyarchy_decline = "Decline in V-Dem polyarchy over window (continuous)"
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
    if (isTRUE(incl)) "" else "_exclyr"
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

apply_table_style <- function(gt_tbl) {
  gt_tbl |>
    gt::tab_options(
      table.font.size = gt::px(12),
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

SIG_FOOTNOTE <- "Significance: * p<0.10, ** p<0.05, *** p<0.01. SE in parentheses."

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
  # Single line: a literal "\n" would render as a space in HTML, and pulling in
  # gt::html() per cell just to get a <br> is not worth it for a two-part label.
  cell_text <- matrix(
    dplyr::if_else(
      is.na(as.vector(est)),
      "--",
      sprintf(
        paste0("%.", digits, "f%s  (N=%s)"),
        as.vector(est),
        sig_stars(as.vector(pval)),
        trimws(format(as.vector(n), trim = TRUE, big.mark = ","))
      )
    ),
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

# Classify a threshold spec WITHOUT looking at data. Returns a list with
# $kind ("none" | "absolute" | "quantile"), $quantile, $absolute and $spec (the
# spec as written, for the record).
parse_threshold <- function(spec, name = "threshold") {
  none <- function(label) {
    list(kind = "none", quantile = NA_real_, absolute = -Inf, spec = label)
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
        name, " is NA. Use -Inf for no restriction; NA is too easy to produce ",
        "by accident to accept as 'no restriction'.",
        call. = FALSE
      )
    }
    if (!is.finite(spec)) {
      return(none(format(spec, trim = TRUE)))
    }
    return(list(
      kind = "absolute", quantile = NA_real_,
      absolute = as.numeric(spec), spec = format(spec, trim = TRUE)
    ))
  }
  if (is.character(spec)) {
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
        kind = "quantile", quantile = p, absolute = NA_real_, spec = spec
      ))
    }
    stop(
      name, ' = "', spec, '" is not a recognized threshold.\n',
      "  Accepted forms:\n",
      "    a number  -- absolute cutoff on the variable's own scale, e.g. 0.6\n",
      '    "qNN"     -- percentile of its distribution, e.g. "q50", "q97.5"\n',
      "                 (lowercase q, no space, NN between 0 and 100)\n",
      "    -Inf      -- no restriction\n",
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
    parsed$absolute <- -Inf
    parsed$label <- "no restriction"
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
# `op` is ">=" or ">", matching whichever convention the caller wants.
apply_threshold <- function(data, var, parsed, name, op = c(">=", ">")) {
  op <- match.arg(op)
  if (parsed$kind == "none") {
    cat(sprintf("  %-16s no restriction (%d rows kept)\n", name, nrow(data)))
    return(data)
  }
  n_before <- nrow(data)
  keep <- !is.na(data[[var]]) &
    if (op == ">=") data[[var]] >= parsed$absolute else data[[var]] > parsed$absolute
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
  exec_power_alt = "Executive power (V-Dem)"
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
  exec_power_alt = "Change, index (0-1)"
)


# Okabe-Ito blue / vermillion / green. Colourblind-safe: worst adjacent-pair
# separation is dE 11.0 under deuteranopia, 25.8 under normal vision, and all
# three clear 3:1 contrast against a white panel. Series are ALSO distinguished
# by point shape, so identity never rests on colour alone. Assigned in fixed
# order, never cycled -- no panel here has more than three series.
SERIES_COLORS <- c("#0072B2", "#D55E00", "#009E73")
SERIES_SHAPES <- c(16, 17, 15)


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

# Pull the robust bias-corrected coefficient/SE/p-value/bandwidth out of an
# rdrobust fit (or NA placeholders if the fit is NULL/failed).
extract_rd <- function(fit) {
  if (is.null(fit)) {
    return(tibble::tibble(
      N = NA_integer_,
      coef = NA_real_,
      se = NA_real_,
      pval = NA_real_,
      bw = NA_real_
    ))
  }
  tibble::tibble(
    N = sum(fit$N),
    coef = unname(fit$coef["Robust", 1]),
    se = unname(fit$se["Robust", 1]),
    pval = unname(fit$pv["Robust", 1]),
    bw = unname(fit$bws[1, 1])
  )
}
