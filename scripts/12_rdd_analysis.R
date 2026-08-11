# ==============================================================================
# Fuzzy RDD estimation: illiberal-party elections -> democratic backsliding
# -> economic outcomes
#
# Loads data/rdd_analysis_data.rds (built by 11_build_rdd_data.R) and runs,
# for every outcome, on the sample selected by the toggles below:
#   - a first stage: rdrobust(y = backsliding_Nyr, x = running_var)
#   - a fuzzy RD:    rdrobust(y = Y_<outcome>, x = running_var, fuzzy = backsliding_Nyr)
#   - a plain reduced-form RD on Y (fuzzy = NULL) for comparison
#
# Used to be baseline + three heterogeneity splits (illiberal_score above
# the sample median, excluding indirectly-elected offices, and a Bermeo-
# taxonomy category split) -- removed after finding (a) the illiberal-
# score restriction meaningfully powers up the first stage (promoted to a
# first-class sample-restriction toggle, ILLIBERAL_CUTOFF, rather than a
# one-off comparison) and (b)/(c) didn't move the estimates much, so
# weren't worth the added output/runtime.
#
# Data: data/rdd_analysis_data.rds
# Output (filenames get suffixes reflecting whichever of SCORE_GAP_MIN /
# ILLIBERAL_CUTOFF are finite, so different threshold runs never overwrite
# each other):
#         output/rdd_results.csv (one row per outcome x spec)
#         output/rdd_first_stage_results.csv (one row per spec -- the
#         standalone first-stage table: N, coefficient, SE, p-value,
#         bandwidth, the way a paper would report it, rather than folding
#         the same first-stage numbers redundantly into every outcome row)
#         output/rdd_plots/*.png -- baseline rdplot() for the first stage
#         and every outcome (binned running-variable trend on each side of
#         the cutoff, for visually checking the discontinuity)
#         output/rdd_plots/*_table.html -- the first-stage and outcomes
#         tables rendered as gt() tables (same convention as
#         11_build_rdd_data.R's tables), econ-paper style (coefficient,
#         SE in parentheses, significance stars)
# ==============================================================================

library(tidyverse)
library(rdrobust)
library(here)
library(patchwork)
library(gt)

data_dir <- here::here("data")
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE)
plots_dir <- file.path(out_dir, "rdd_plots")
dir.create(plots_dir, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Toggles
# ------------------------------------------------------------------------------

# Minimum standardized top-2 score gap (score_gap_z, built in
# 11_build_rdd_data.R) an election must clear to be included. score_gap_z
# = score_gap / sd(illiberalism_score across every top-2 party-year in
# that country) -- i.e. the raw illiberal-vs-other gap scaled by how much
# that country's own parties typically differ in ideology, pooling BOTH
# top-2 members across every scored election in the country (not the
# distribution of election-level gaps -- an earlier version did that and
# had a mechanical ceiling: a 2-election country could never exceed |z| ~
# 0.7 no matter how large its actual gaps were, since the gap was being
# standardized against its own tiny sample). Since illiberal_score is
# defined as the HIGHER of the top-2's scores, score_gap -- and therefore
# score_gap_z -- is always >= 0, so SCORE_GAP_MIN <- 0 is a NO-OP (keeps
# ~every election); meaningful thresholds start above 0, e.g. 1 means
# "this election's gap is at least one country-scaled party-ideology-SD."
# -Inf = no filter (every election, including the few countries with too
# little data to compute their own score SD, where score_gap_z is NA).
# Restricting the ENTIRE analysis (baseline, every heterogeneity split,
# and all plots) this way tests whether the first stage is diluted by
# low-contrast top-2 pairs -- in the full sample, 29% of elections have a
# top-2 score_gap under 0.05, i.e. the two "top-2" parties are barely
# distinguishable in illiberalism, so crossing the vote-share cutoff there
# doesn't correspond to a real liberal/illiberal treatment contrast, just
# noise. Unlike a one-off heterogeneity subsample, this runs the SAME full
# plot pipeline as the main analysis (not just the results table) on the
# restricted sample, so the discontinuity (or lack of one) can be
# inspected visually.
SCORE_GAP_MIN <- 0
stopifnot(is.numeric(SCORE_GAP_MIN), length(SCORE_GAP_MIN) == 1)

# Minimum illiberal_score (the WINNING/more-illiberal top-2 member's raw
# level on ILLIBERALISM_VAR, e.g. v2xpa_antiplural in [0,1] -- not the
# gap between the top-2, which SCORE_GAP_MIN above already restricts) an
# election must clear to be included. Originally a one-off heterogeneity
# comparison (elections above the SAMPLE MEDIAN illiberal_score); promoted
# to a first-class sample restriction after finding it meaningfully powers
# up the first stage -- restricting to elections where the illiberal
# party actually clears an absolute ideological bar (not just edges out
# the other top-2 member on a relative basis) seems to matter for whether
# crossing the vote-share cutoff corresponds to a real change in governing
# ideology. -Inf = no filter.
ILLIBERAL_CUTOFF <- 0.6
stopifnot(is.numeric(ILLIBERAL_CUTOFF), length(ILLIBERAL_CUTOFF) == 1)

d <- readRDS(file.path(data_dir, "rdd_analysis_data.rds"))

# Label every output (specs, plot/CSV filenames) with whichever thresholds
# are active so different runs never silently overwrite each other.
fmt_threshold_suffix <- function(prefix, value) {
  if (is.finite(value)) {
    paste0(
      prefix,
      gsub("-", "neg", gsub("\\.", "p", format(value, trim = TRUE)))
    )
  } else {
    ""
  }
}
score_gap_suffix <- fmt_threshold_suffix("_scoregap_ge_", SCORE_GAP_MIN)
illiberal_suffix <- fmt_threshold_suffix("_illib_gt_", ILLIBERAL_CUTOFF)
baseline_label <- paste0("baseline", score_gap_suffix, illiberal_suffix)

if (is.finite(SCORE_GAP_MIN)) {
  n_before <- nrow(d)
  d <- d |> filter(!is.na(score_gap_z), score_gap_z >= SCORE_GAP_MIN)
  cat(sprintf(
    "SCORE_GAP_MIN = %s -- restricted from %d to %d elections (score_gap_z >= %s within-country)\n",
    SCORE_GAP_MIN,
    n_before,
    nrow(d),
    SCORE_GAP_MIN
  ))
}

if (is.finite(ILLIBERAL_CUTOFF)) {
  n_before <- nrow(d)
  d <- d |> filter(illiberal_score > ILLIBERAL_CUTOFF)
  cat(sprintf(
    "ILLIBERAL_CUTOFF = %s -- restricted from %d to %d elections (illiberal_score > %s)\n",
    ILLIBERAL_CUTOFF,
    n_before,
    nrow(d),
    ILLIBERAL_CUTOFF
  ))
}

# ------------------------------------------------------------------------------
# Diagnostic: density of the running variable (on whatever sample
# SCORE_GAP_MIN above selects). Standard RD sanity check: if elections
# near the cutoff were somehow selected or sorted (e.g. the illiberal side
# systematically squeaking out narrow wins/losses via manipulation), the
# running variable's density would show bunching or a jump right at 0. A
# smooth, continuous-looking density on both sides is consistent with the
# "no precise sorting around the cutoff" identifying assumption a fuzzy RD
# relies on -- this is a visual check only, not a formal manipulation test
# (e.g. McCrary 2008 / Cattaneo-Jansson-Ma's rddensity).
# ------------------------------------------------------------------------------

density_plot_file <- paste0(
  "running_var_density",
  score_gap_suffix,
  illiberal_suffix,
  ".png"
)

running_var_density_plot <- ggplot(d, aes(x = running_var)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 60,
    fill = "grey80",
    color = "white"
  ) +
  geom_density(color = "darkblue", linewidth = 1) +
  geom_vline(
    xintercept = 0,
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    title = sprintf(
      "Density of the running variable (%s, N = %d)",
      baseline_label,
      nrow(d)
    ),
    subtitle = "Dashed red line = RD cutoff.",
    x = "Running variable (illiberal - other vote/seat share)",
    y = "Density"
  ) +
  theme_minimal()

ggsave(
  file.path(plots_dir, density_plot_file),
  running_var_density_plot,
  width = 8,
  height = 3
)
cat(sprintf("Saved output/rdd_plots/%s\n", density_plot_file))

outcome_vars <- c(
  "Y_gdp_growth",
  "Y_inflation",
  "Y_unemployment",
  "Y_trade_pct_gdp",
  "Y_top10_share",
  "Y_gini_disp",
  "Y_gini_mkt"
)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# Pull the robust bias-corrected coefficient/p-value/bandwidth out of an
# rdrobust fit (or NA placeholders if the fit is NULL/failed).
extract_rd <- function(fit) {
  if (is.null(fit)) {
    return(tibble(
      N = NA_integer_,
      coef = NA_real_,
      se = NA_real_,
      pval = NA_real_,
      bw = NA_real_
    ))
  }
  tibble(
    N = sum(fit$N),
    coef = unname(fit$coef["Robust", 1]),
    se = unname(fit$se["Robust", 1]),
    pval = unname(fit$pv["Robust", 1]),
    bw = unname(fit$bws[1, 1])
  )
}

# bwselect = "mserd" (Calonico, Cattaneo & Titiunik 2014's MSE-optimal,
# common-bandwidth selector) is rdrobust()'s own default -- named
# explicitly here so the bandwidth-selection method is documented in code,
# not left as an implicit default, and so save_rdplot() below can request
# the identical procedure when it separately re-fits to get a bandwidth
# for plotting.
RD_BWSELECT <- "mserd"

# rdplot()'s default binselect = "esmv" (mimicking-variance evenly-spaced
# bins) picks a relatively coarse number of bins per side on its own --
# e.g. 24 left / 13 right for the full Y_gdp_growth sample. `scale`
# multiplies that automatically-selected count (2 = twice as many bins per
# side) rather than hand-picking a fixed nbins that wouldn't adapt across
# very different subsample sizes (baseline N=1347 vs. a SCORE_GAP_MIN-
# restricted N=221).
RD_BIN_SCALE <- 2

safe_rdrobust <- function(y, x, fuzzy = NULL) {
  tryCatch(
    {
      # is.null(fuzzy) | !is.na(fuzzy) would silently collapse to logical(0)
      # when fuzzy is NULL (!is.na(NULL) is zero-length), zeroing out `keep`
      # via vectorized `&` -- branch explicitly instead.
      fuzzy_ok <- if (is.null(fuzzy)) rep(TRUE, length(y)) else !is.na(fuzzy)
      keep <- !is.na(y) & !is.na(x) & fuzzy_ok
      if (sum(keep) < 20) {
        return(NULL)
      }
      if (is.null(fuzzy)) {
        rdrobust(y = y[keep], x = x[keep], bwselect = RD_BWSELECT)
      } else {
        rdrobust(
          y = y[keep],
          x = x[keep],
          fuzzy = fuzzy[keep],
          bwselect = RD_BWSELECT
        )
      }
    },
    error = function(e) NULL
  )
}

# Builds a 95% confidence ribbon around the local-linear fit, one side of
# the cutoff at a time, by directly fitting the SAME weighted least-squares
# regression rdrobust()/rdplot() use for the point estimate: OLS with
# triangular-kernel weights (1 - |x|/h, clipped at 0) restricted to |x| <=
# h. This is not an approximation -- verified against a live fit, the
# weighted-OLS intercept-at-0 reproduces rdrobust()'s own "Conventional"
# jump estimate exactly (both gave -0.04603014 for a test outcome). The
# resulting ribbon uses the standard weighted-least-squares prediction SE,
# which is NOT the same as CCT's more sophisticated bias-corrected
# "Robust" SE reported in the results table -- this is a faithful
# visualization of the conventional local-linear fit's uncertainty, not a
# substitute for the table's formal inference.
build_ci_ribbon <- function(y, x, h_left, h_right, x_zoom, n_grid = 100) {
  fit_side <- function(side, bandwidth, grid) {
    w <- pmax(1 - abs(x) / bandwidth, 0)
    idx <- which(side & w > 0)
    if (length(idx) < 3) {
      return(NULL)
    }
    df <- data.frame(xx = x[idx], yy = y[idx], ww = w[idx])
    fit <- tryCatch(lm(yy ~ xx, data = df, weights = ww), error = function(e) {
      NULL
    })
    if (is.null(fit)) {
      return(NULL)
    }
    pred <- predict(fit, newdata = data.frame(xx = grid), se.fit = TRUE)
    tcrit <- qt(0.975, df = fit$df.residual)
    tibble(
      xx = grid,
      ymin = pred$fit - tcrit * pred$se.fit,
      ymax = pred$fit + tcrit * pred$se.fit
    )
  }
  dplyr::bind_rows(
    fit_side(x < 0, h_left, seq(-x_zoom, 0, length.out = n_grid)),
    fit_side(x >= 0, h_right, seq(0, x_zoom, length.out = n_grid))
  )
}

# rdrobust::rdplot() is the standard RD visualization: binned local means of
# y against the running variable, with a separate local-polynomial fit on
# each side of the cutoff -- lets you see the discontinuity (or lack of
# one) directly, rather than just reading a coefficient/p-value off the
# results table. pdf(NULL) swallows its default auto-print (rdplot() prints
# to whatever device is open, which errors/hangs in a non-interactive
# Rscript session with no device). Returns the finished ggplot object (or
# NULL if there's too little data / the fit fails) -- callers decide
# whether to ggsave() it standalone (save_rdplot(), for the first stage)
# or combine several into one multi-panel figure (build_combined_outcomes_
# plot(), for the outcomes).
build_rdplot <- function(y, x, title, y_label, y.lim = NULL) {
  keep <- !is.na(y) & !is.na(x)
  if (sum(keep) < 20) {
    message("Skipping plot (too few observations): ", title)
    return(NULL)
  }
  y <- y[keep]
  x <- x[keep]

  # rdplot()'s own h defaults to spanning the FULL support of x on each
  # side -- a much wider window than what rdrobust() actually uses for the
  # point estimate/inference. Re-fit rdrobust() here (same RD_BWSELECT =
  # "mserd" CCT procedure used in run_spec()) purely to pull its
  # MSE-optimal bandwidth, then pass that same h into rdplot() so the
  # plotted local-linear fit reflects the identical window the results
  # table's estimate is actually computed from, rather than a
  # wider/unrelated one. Falls back to rdplot()'s own default (full
  # support) if that fit fails.
  rd_fit <- safe_rdrobust(y, x)
  h <- if (!is.null(rd_fit)) unname(rd_fit$bws["h", ]) else NULL
  # Zoom the x-axis to exactly the optimal bandwidth -- the fitted line
  # never extends past +/-h anyway, so padding the zoom out further (an
  # earlier version used 3x) only wastes plot area on a region the
  # estimate doesn't use.
  x_zoom <- if (!is.null(h)) max(h) else max(abs(x), na.rm = TRUE)

  h_left <- if (!is.null(h)) h[1] else NULL
  h_right <- if (!is.null(h)) h[length(h)] else NULL
  ribbon_data <- if (!is.null(h)) {
    build_ci_ribbon(y, x, h_left, h_right, x_zoom)
  } else {
    NULL
  }

  # rdplot()'s x.lim only crops the DISPLAY viewport (coord_cartesian-style)
  # -- bins outside that window still feed the automatic y-axis scaling.
  # A single sparse, high-variance bin far from the cutoff (common with
  # this sample's uneven coverage) can carry a mean wildly off from the
  # visible region, which would silently squash the whole well-behaved
  # area into a sliver. Probe with hide=TRUE (computes bins, skips
  # rendering) to see only the bins that will actually be visible, and
  # derive y.lim from those PLUS the ribbon's own extent -- unless the
  # caller already passed an explicit y.lim (the first-stage plot
  # hard-codes [0,1] since it's a probability, and shouldn't be overridden
  # here).
  if (is.null(y.lim)) {
    bins_probe <- tryCatch(
      rdplot(
        y = y,
        x = x,
        p = 1,
        h = h,
        kernel = "triangular",
        scale = RD_BIN_SCALE,
        hide = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(bins_probe)) {
      visible <- bins_probe$vars_bins |>
        dplyr::filter(abs(rdplot_mean_bin) <= x_zoom)
      bounds <- c(visible$rdplot_mean_y, ribbon_data$ymin, ribbon_data$ymax)
      bounds <- bounds[is.finite(bounds)]
      if (length(bounds) > 0) {
        rng <- range(bounds)
        pad <- diff(rng) * 0.1
        y.lim <- c(rng[1] - pad, rng[2] + pad)
      }
    }
  }

  tryCatch(
    {
      pdf(NULL)
      fit <- rdplot(
        y = y,
        x = x,
        title = title,
        x.label = "Running variable (illiberal - other vote/seat share)",
        y.label = y_label,
        x.lim = c(-1, 1) * x_zoom,
        # rdplot()'s default p=4 (quartic) fits a separate 4th-degree
        # polynomial on each side of the cutoff -- with only ~15-20 binned
        # points per side here, that's enough flexibility to snake through
        # nearly every bin, producing the wildly oscillating "spine" seen
        # at p=4. p=1 (local linear) matches what rdrobust() actually
        # estimates for the point estimate/inference, so the plotted curve
        # reflects the same model the results table reports, instead of a
        # more flexible and misleading one. kernel="triangular" matches
        # rdrobust()'s own default weighting (rdplot() otherwise defaults
        # to "uniform") -- without this, the plotted fit used the right
        # bandwidth and polynomial order but a different weighting scheme
        # than the actual point estimate it's meant to visualize.
        p = 1,
        h = h,
        kernel = "triangular",
        scale = RD_BIN_SCALE,
        y.lim = y.lim
      )
      dev.off()
      # Per-bin CI whiskers were tried and dropped in favor of a single
      # shaded confidence ribbon around the fit line itself (see
      # build_ci_ribbon() above) -- less visual clutter, and it directly
      # answers "how uncertain is the estimated jump" rather than "how
      # uncertain is each individual bin's mean." Prepend (not append) the
      # ribbon layer so it renders BEHIND the bin dots and fit line rather
      # than covering them.
      if (!is.null(ribbon_data) && nrow(ribbon_data) > 0) {
        ribbon_layer <- geom_ribbon(
          data = ribbon_data,
          aes(x = xx, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill = "red",
          alpha = 0.15
        )
        fit$rdplot$layers <- append(
          fit$rdplot$layers,
          list(ribbon_layer),
          after = 0
        )
      }
      fit$rdplot
    },
    error = function(e) {
      if (!is.null(dev.list())) {
        dev.off()
      }
      message(
        "Skipping plot (rdplot failed): ",
        title,
        " -- ",
        conditionMessage(e)
      )
      NULL
    }
  )
}

# Thin wrapper around build_rdplot() for plots saved standalone (just the
# first stage now -- the outcomes are combined into one multi-panel figure
# by build_combined_outcomes_plot() below).
save_rdplot <- function(y, x, file_name, title, y_label, y.lim = NULL) {
  p <- build_rdplot(y, x, title, y_label, y.lim)
  if (!is.null(p)) {
    ggsave(file.path(plots_dir, file_name), p, width = 7, height = 3)
  }
  invisible(p)
}

# One combined figure for all outcomes in outcome_vars, arranged in a grid
# via patchwork rather than one PNG per outcome -- each panel keeps its
# own independently-selected bandwidth/bins (outcomes have different
# missingness and variance, so a shared binning scheme across a facet_wrap
# would be wrong here), but they're assembled into a single image with one
# overall title.
build_combined_outcomes_plot <- function(data, sample_label, ncol = 3) {
  panels <- purrr::map(outcome_vars, function(oc) {
    build_rdplot(data[[oc]], data$running_var, oc, oc)
  })
  panels <- purrr::compact(panels) # drop any outcome that failed/had too few obs
  if (length(panels) == 0) {
    return(NULL)
  }
  patchwork::wrap_plots(panels, ncol = ncol) +
    patchwork::plot_annotation(
      title = sprintf("Reduced-form RD by outcome (%s)", sample_label),
      subtitle = "Shaded band = 95% CI of the local-linear fit"
    )
}

# One RD plot per outcome (plus the first stage) for a given subsample --
# called once for baseline; call again for any other spec worth visualizing.
make_rd_plots <- function(data, spec_label, treatment_col = "backsliding_Nyr") {
  # Label every plot with both the spec and the sample size it was
  # estimated on, so a plot is self-describing without needing to cross-
  # reference its filename or the console output.
  sample_label <- sprintf("%s, N = %d", spec_label, nrow(data))
  save_rdplot(
    data[[treatment_col]],
    data$running_var,
    sprintf("%s_first_stage.png", spec_label),
    sprintf("First stage (%s)", sample_label),
    "P(backsliding within window)",
    # treatment_col is a 0/1 probability -- fix the y-axis to its actual
    # range instead of letting a wide-CI bin (common with small N) blow
    # the scale out to +/-5, which was making every first-stage plot look
    # far more zoomed-out than the data warrants.
    y.lim = c(0, 1)
  )
  combined <- build_combined_outcomes_plot(data, sample_label, ncol = 3)
  if (!is.null(combined)) {
    n_rows <- ceiling(length(outcome_vars) / 2)
    ggsave(
      file.path(plots_dir, sprintf("%s_outcomes.png", spec_label)),
      combined,
      width = 16,
      height = 2.5 * n_rows
    )
  }
}

# Standard economics-paper table conventions: significance stars on the
# coefficient, standard error in parentheses immediately after -- e.g.
# "-0.045 (0.023)*".
sig_stars <- function(pval) {
  dplyr::case_when(
    is.na(pval) ~ "",
    pval < 0.01 ~ "***",
    pval < 0.05 ~ "**",
    pval < 0.10 ~ "*",
    TRUE ~ ""
  )
}

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

# Same gt() + gtsave() convention used in 11_build_rdd_data.R's tables --
# saved as HTML rather than PNG, since PNG export needs a headless-browser
# backend (webshot2/chromote, or webshot+PhantomJS) that isn't set up in
# this environment. HTML also sizes to its content automatically, so there
# is no manual width/clipping bookkeeping to get wrong.
apply_table_style <- function(gt_tbl) {
  gt_tbl |>
    tab_options(
      table.font.size = px(12),
      table.border.top.style = "solid",
      table.border.top.width = px(2),
      table.border.top.color = "black",
      table.border.bottom.style = "solid",
      table.border.bottom.width = px(2),
      table.border.bottom.color = "black",
      column_labels.border.top.style = "solid",
      column_labels.border.top.width = px(2),
      column_labels.border.top.color = "black",
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = px(1.5),
      column_labels.border.bottom.color = "black",
      table_body.hlines.style = "solid",
      table_body.hlines.width = px(0.5),
      table_body.hlines.color = "#cccccc"
    )
}

SIG_FOOTNOTE <- "Significance: * p<0.10, ** p<0.05, *** p<0.01. SE in parentheses."

save_table_html <- function(tbl, file_name, title) {
  gt_tbl <- tbl |>
    gt::gt() |>
    gt::tab_header(title = title) |>
    gt::tab_source_note(source_note = SIG_FOOTNOTE) |>
    gt::opt_row_striping() |>
    apply_table_style()
  gt::gtsave(gt_tbl, file.path(plots_dir, file_name))
  cat(sprintf("Saved output/rdd_plots/%s\n", file_name))
}

save_first_stage_table <- function(first_stage, file_name, title) {
  tbl <- first_stage |>
    transmute(
      Spec = spec,
      N = n_first_stage,
      `First stage` = fmt_est(
        first_stage_coef,
        first_stage_se,
        first_stage_pval
      ),
      Bandwidth = round(first_stage_bandwidth, 2)
    )
  save_table_html(tbl, file_name, title)
}

save_outcomes_table <- function(outcomes, file_name, title) {
  tbl <- outcomes |>
    transmute(
      Outcome = outcome,
      N = n_outcome,
      `Reduced form` = fmt_est(rd_estimate, rd_se, rd_pval),
      `Fuzzy RD (LATE)` = fmt_est(late_estimate, late_se, late_pval)
    )
  save_table_html(tbl, file_name, title)
}

# Estimates a full rdrobust() fit for the first stage (once per subsample)
# plus a fuzzy RD for every outcome (reusing that same treatment_col), and
# returns BOTH a standalone first-stage row (full N/coef/se/pval/bandwidth,
# the way a paper's first-stage table would report it -- not just a
# coefficient and p-value folded redundantly into every outcome row) and
# the per-outcome table.
run_spec <- function(data, spec_label, treatment_col = "backsliding_Nyr") {
  fs_fit <- safe_rdrobust(data[[treatment_col]], data$running_var)
  fs <- extract_rd(fs_fit)

  first_stage <- tibble(
    spec = spec_label,
    n_first_stage = fs$N,
    first_stage_coef = fs$coef,
    first_stage_se = fs$se,
    first_stage_pval = fs$pval,
    first_stage_bandwidth = fs$bw
  )

  outcomes <- map_dfr(outcome_vars, function(oc) {
    rf_fit <- safe_rdrobust(data[[oc]], data$running_var)
    fuzzy_fit <- safe_rdrobust(
      data[[oc]],
      data$running_var,
      fuzzy = data[[treatment_col]]
    )
    rf <- extract_rd(rf_fit)
    late <- extract_rd(fuzzy_fit)

    tibble(
      outcome = oc,
      spec = spec_label,
      n_first_stage = fs$N,
      first_stage_coef = fs$coef,
      first_stage_pval = fs$pval,
      n_outcome = rf$N,
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
# Baseline
# ------------------------------------------------------------------------------

cat(sprintf("=== Baseline (%s) ===\n", baseline_label))
results_baseline <- run_spec(d, baseline_label)
save_first_stage_table(
  results_baseline$first_stage,
  sprintf("%s_first_stage_table.html", baseline_label),
  sprintf("First stage (%s)", baseline_label)
)
save_outcomes_table(
  results_baseline$outcomes,
  sprintf("%s_outcomes_table.html", baseline_label),
  sprintf("Reduced form / fuzzy RD by outcome (%s)", baseline_label)
)
make_rd_plots(d, baseline_label)
cat(sprintf("Saved %s RD plots to %s/\n", baseline_label, plots_dir))

# ------------------------------------------------------------------------------
# Combine + save
# ------------------------------------------------------------------------------

all_results <- results_baseline$outcomes
all_first_stage <- results_baseline$first_stage

n_failed <- sum(is.na(all_results$late_estimate))
cat(sprintf(
  "\n=== Overall: %d / %d outcome x spec combinations produced a fuzzy RD estimate (%d skipped/failed -- too few observations near the cutoff or rdrobust's internal checks failed) ===\n",
  nrow(all_results) - n_failed,
  nrow(all_results),
  n_failed
))

results_file <- paste0(
  "rdd_results",
  score_gap_suffix,
  illiberal_suffix,
  ".csv"
)
first_stage_file <- paste0(
  "rdd_first_stage_results",
  score_gap_suffix,
  illiberal_suffix,
  ".csv"
)

write_csv(all_results, file.path(out_dir, results_file))
write_csv(all_first_stage, file.path(out_dir, first_stage_file))
message(sprintf(
  "Saved output/%s and output/%s",
  results_file,
  first_stage_file
))
