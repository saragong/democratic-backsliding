# ==============================================================================
# What ELSE does a narrow anti-pluralist victory deliver?
#
# The reduced-form result is that when the more anti-pluralist of the top 2
# narrowly wins, GDP growth over the next five years is lower and inequality
# higher. The obvious objection is that anti-pluralism travels with something
# else -- populism, economic left, social conservatism -- and that the
# "something else" is doing the work.
#
# This script measures the travelling directly. Same running variable, same
# cutoff, same estimator as the headline RD; only the OUTCOME changes, to a
# property of the winner rather than a property of the country:
#
#   binary      is the winner also the MORE populist / more left / more TAN /
#               more anti-elite / less minority-friendly of the two?
#   continuous  the winner's own score on that index.
#
# Reading it: a coefficient near zero says the narrow anti-pluralist victory is
# NOT also a narrow victory on that other dimension, so the growth result
# cannot be that dimension in disguise. A large coefficient says the two are
# bundled at the cutoff and the reduced form cannot separate them.
#
# The binary outcome has a useful property: at the cutoff the anti-pluralist
# wins by construction, so E[binary] is "given the anti-pluralist just won, how
# often is it also the more populist party". A jump of zero means crossing the
# cutoff does not change that -- the two indices are orthogonal AT the margin,
# which is a sharper statement than the 0.115 unconditional correlation.
#
# There is no treatment here and no fuzzy RD: these are properties of the
# election, not outcomes that backsliding could mediate. Reduced form only.
#
# The window is irrelevant to every outcome above -- they are all measured at
# the election -- so the script reads a single build and says so.
#
#   Rscript --no-init-file scripts/17_party_outcomes_rdd.R
#
# Output: output/runs/_sweeps/party_outcome_rdd_<instr>_w<N><suffix>/
#           party_outcome_results.csv    one row per score x form x restriction
#           party_outcomes_table.html    the unrestricted sample, flat
#           grid_party_outcomes.html     restriction x score colour grid
#           plots/party_outcomes_<form>.png
# ==============================================================================

library(tidyverse)
library(rdrobust)
library(here)
library(patchwork)
library(gt)

source(here::here("scripts", "rdd_helpers.R"))

data_dir <- here::here("data")

# ---- toggles -----------------------------------------------------------------

if (!exists("ILLIBERALISM_VAR")) ILLIBERALISM_VAR <- "v2xpa_antiplural"
# Any window gives the same answer -- see the header -- so this only selects
# which build file to open.
if (!exists("BACKSLIDING_WINDOW_YEARS")) BACKSLIDING_WINDOW_YEARS <- 5
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- FALSE
}

# Every party score except the one currently defining the running variable.
# Including it would produce a mechanical result: illiberal_score is the max of
# the two by construction, so "is the winner the more anti-pluralist one" is
# exactly the treatment indicator and jumps from 0 to 1 at the cutoff.
if (!exists("COMPARISON_SCORES")) {
  COMPARISON_SCORES <- setdiff(PARTY_SCORE_VARS, ILLIBERALISM_VAR)
}

# Restriction levels swept on each axis. Quantiles of the FULL build, matching
# the vocabulary 13_restriction_grid.R uses, so "q60" means the same thing in
# both places.
if (!exists("RESTRICTION_SPECS")) {
  RESTRICTION_SPECS <- list(-Inf, "q20", "q40", "q60", "q80")
}

# The two ways to say "the top 2 are far enough apart in illiberality to make
# crossing the cutoff a real contrast" -- the ask's "minimum illiberality score
# gap". Raw and within-country-standardized, because the raw gap is not
# comparable between a country whose parties all cluster and one whose parties
# are spread out.
if (!exists("GAP_AXES")) {
  GAP_AXES <- c(
    score_gap = "Top-2 illiberality gap (raw)",
    score_gap_z = "Top-2 illiberality gap (within-country SDs)"
  )
}

build_suffix <- if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "" else "_exclyr"
build_path <- file.path(
  data_dir, "rdd_build",
  sprintf(
    "rdd_%s_w%d%s.rds",
    ILLIBERALISM_VAR, BACKSLIDING_WINDOW_YEARS, build_suffix
  )
)
if (!file.exists(build_path)) {
  stop("No build at ", build_path, ". Run 11_build_rdd_data.R first.")
}
d_full <- readRDS(build_path)

out_dir <- sweep_dir(sprintf(
  "party_outcome_rdd_%s_w%d%s",
  INSTRUMENT_LABELS[[ILLIBERALISM_VAR]], BACKSLIDING_WINDOW_YEARS, build_suffix
))
plots_dir <- file.path(out_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Loaded %s (%d elections)\n", basename(build_path), nrow(d_full)))

# ---- outcomes ----------------------------------------------------------------

# Built from the __winner / __loser columns the build already carries, so
# nothing here re-derives a party match.
missing_cols <- unlist(lapply(
  COMPARISON_SCORES, \(s) paste0(s, c("__winner", "__loser"))
))
missing_cols <- setdiff(missing_cols, names(d_full))
if (length(missing_cols) > 0) {
  stop(
    "Build predates PARTY_SCORE_VARS and is missing: ",
    paste(missing_cols, collapse = ", "),
    ". Rebuild it with 11_build_rdd_data.R."
  )
}

for (s in COMPARISON_SCORES) {
  w <- d_full[[paste0(s, "__winner")]]
  l <- d_full[[paste0(s, "__loser")]]
  # NA if either side is unscored, rather than FALSE: "we don't know which was
  # higher" is not "the winner was not higher". ep_galtan is scored for both
  # top-2 members in only ~4% of elections, so this matters.
  d_full[[paste0("Zbin_", s)]] <- ifelse(is.na(w) | is.na(l), NA_real_,
                                         as.numeric(w > l))
  d_full[[paste0("Zcont_", s)]] <- w
}

FORMS <- list(
  binary = list(
    prefix = "Zbin_",
    title = "Is the winner also the higher-scoring of the top 2?",
    ylab = "P(winner scores higher)",
    ylim = c(0, 1)
  ),
  continuous = list(
    prefix = "Zcont_",
    title = "The winner's own score",
    ylab = "Score (index units)",
    ylim = NULL
  )
)

# ---- estimation --------------------------------------------------------------

# Reduced form only: no treatment, so no first stage and no fuzzy RD.
estimate_cell <- function(data, var) {
  fit <- extract_rd(safe_rdrobust(data[[var]], data$running_var))
  tibble(
    n = fit$N, est = fit$coef, se = fit$se, pval = fit$pval, bw = fit$bw
  )
}

# One restriction level: resolve against the FULL build (never the already-cut
# frame -- otherwise the levels would not be the quantiles they claim to be),
# apply, estimate every score x form.
run_level <- function(axis_var, spec) {
  thr <- resolve_threshold(
    parse_threshold(spec, axis_var), d_full[[axis_var]], axis_var
  )
  dd <- d_full[!is.na(d_full[[axis_var]]) & d_full[[axis_var]] >= thr$absolute, ,
               drop = FALSE]
  if (thr$kind == "none") dd <- d_full

  if (nrow(dd) < RD_MIN_OBS) {
    return(tibble())
  }
  expand_grid(form = names(FORMS), score = COMPARISON_SCORES) |>
    pmap_dfr(function(form, score) {
      var <- paste0(FORMS[[form]]$prefix, score)
      bind_cols(
        tibble(
          axis = axis_var,
          axis_label = unname(GAP_AXES[axis_var]),
          level = thr$spec,
          cutoff = thr$absolute,
          n_sample = nrow(dd),
          form = form,
          score = score,
          score_label = unname(PARTY_SCORE_DISPLAY[score])
        ),
        estimate_cell(dd, var)
      )
    })
}

results <- expand_grid(axis = names(GAP_AXES), spec = RESTRICTION_SPECS) |>
  pmap_dfr(function(axis, spec) run_level(axis, spec))

write_csv(results, file.path(out_dir, "party_outcome_results.csv"))

# ---- flat table, unrestricted sample ----------------------------------------

subtitle_base <- sprintf(
  "Instrument: %s | Running variable and cutoff exactly as in the headline RD | Build: %s",
  INSTRUMENT_DISPLAY[[ILLIBERALISM_VAR]], basename(build_path)
)

# Anti-elitism correlates 0.966 with populism in this sample, so its two rows
# are very close to the populism rows by construction. Saying so on the table
# stops that being read as two independent confirmations.
REDUNDANCY_NOTE <- paste(
  "v2paanteli (anti-elitism) correlates 0.966 with v2xpa_popul (populism)",
  "across the top-2 sample, so its rows are near-duplicates of the populism",
  "rows rather than independent evidence.",
  "v2paminor runs the opposite way from the other scores: HIGHER = more",
  "supportive of minority rights.",
  "ep_galtan is scored for both top-2 members in only about 4% of elections,",
  "so its N is far smaller than the rest."
)

flat <- results |>
  filter(level == "-Inf", axis == names(GAP_AXES)[1]) |>
  transmute(
    Form = if_else(form == "binary", "Winner scores higher (0/1)", "Winner's score"),
    Score = score_label,
    N = n,
    `Reduced-form RD` = fmt_est(est, se, pval),
    Bandwidth = round(bw, 2)
  ) |>
  arrange(Form, Score)

save_table_html(
  flat,
  file.path(out_dir, "party_outcomes_table.html"),
  "Does a narrow anti-pluralist victory also shift the winner's other scores?",
  sprintf("%s | Sample: all scored elections (N = %d)", subtitle_base, nrow(d_full)),
  note = paste(SIG_FOOTNOTE, REDUNDANCY_NOTE)
)

# ---- restriction x score grid -----------------------------------------------

# One block per (gap axis x outcome form): rows are restriction levels, columns
# are the comparison scores. Same idiom as 13_restriction_grid.R's pair grids.
blocks <- list()
for (ax in names(GAP_AXES)) {
  for (fm in names(FORMS)) {
    sub <- results |> filter(axis == ax, form == fm)
    if (nrow(sub) == 0) next
    lvls <- unique(sub$level)
    row_labels <- vapply(lvls, function(l) {
      r <- sub |> filter(level == l) |> slice(1)
      if (l == "-Inf") {
        sprintf("all (N = %d)", r$n_sample)
      } else {
        sprintf("%s: >= %.3g (N = %d)", l, r$cutoff, r$n_sample)
      }
    }, character(1))

    mat <- function(col) {
      m <- matrix(
        NA_real_, nrow = length(lvls), ncol = length(COMPARISON_SCORES),
        dimnames = list(row_labels, unname(PARTY_SCORE_LABELS[COMPARISON_SCORES]))
      )
      for (i in seq_along(lvls)) {
        for (j in seq_along(COMPARISON_SCORES)) {
          v <- sub |> filter(level == lvls[i], score == COMPARISON_SCORES[j])
          if (nrow(v) == 1) m[i, j] <- v[[col]]
        }
      }
      m
    }
    blocks[[length(blocks) + 1]] <- list(
      label = sprintf("%s -- %s", unname(GAP_AXES[ax]), FORMS[[fm]]$title),
      row_axis = "Restriction level",
      col_axis = "Comparison score",
      row_labels = row_labels,
      col_labels = unname(PARTY_SCORE_LABELS[COMPARISON_SCORES]),
      est = mat("est"), se = mat("se"), pval = mat("pval"), n = mat("n")
    )
  }
}

if (length(blocks) > 0) {
  color_grid_blocks(
    blocks,
    path = file.path(out_dir, "grid_party_outcomes.html"),
    title = "Reduced-form RD on the winner's other party scores, by illiberality-gap restriction",
    subtitle = subtitle_base,
    note = paste(
      SIG_FOOTNOTE, REDUNDANCY_NOTE,
      "Every cell is a different sample restriction a researcher could have",
      "chosen. Read the surface, not the best cell. No multiple-testing",
      "correction is applied."
    )
  )
}

# ---- figures -----------------------------------------------------------------

# build_panel_plot() is the same function 12_rdd_analysis.R uses, so these RD
# panels are drawn identically to the headline outcome panels.
for (fm in names(FORMS)) {
  spec <- FORMS[[fm]]
  series <- setNames(
    unname(PARTY_SCORE_LABELS[COMPARISON_SCORES]),
    paste0(spec$prefix, COMPARISON_SCORES)
  )
  # One panel per score rather than all five on shared axes: the binary and
  # continuous scores are on genuinely different scales (v2xpa_* on [0, 1],
  # v2pariglef_neg and v2paanteli on expert scales running past 4), and
  # build_panel_plot() is explicit that a panel must never mix two units.
  panels <- list()
  for (s in COMPARISON_SCORES) {
    v <- paste0(spec$prefix, s)
    p <- build_panel_plot(
      d_full,
      setNames(unname(PARTY_SCORE_LABELS[s]), v),
      unname(PARTY_SCORE_DISPLAY[s]),
      spec$ylab,
      y_lim = spec$ylim
    )
    if (!is.null(p)) panels[[s]] <- p
  }
  if (length(panels) == 0) next
  combined <- wrap_plots(panels, ncol = 3) +
    plot_annotation(
      title = sprintf("%s (%s)", spec$title, INSTRUMENT_DISPLAY[[ILLIBERALISM_VAR]]),
      subtitle = paste(
        "All scored elections. Shaded band = 95% CI of the local-linear fit",
        "(conventional, not bias-corrected)."
      )
    )
  n_rows <- ceiling(length(panels) / 3)
  ggsave(
    file.path(plots_dir, sprintf("party_outcomes_%s.png", fm)),
    combined, width = 15, height = 3.6 * n_rows + 0.6, dpi = 150,
    limitsize = FALSE
  )
  cat(sprintf("Saved plots/party_outcomes_%s.png\n", fm))
}

cat("\n=== Unrestricted sample ===\n")
print(as.data.frame(flat), row.names = FALSE)
message("\nParty-outcome RDs written to ", out_dir)
