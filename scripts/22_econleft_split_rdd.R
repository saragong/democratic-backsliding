# ==============================================================================
# The main RDD, split by whether the anti-pluralist is the RIGHT or the LEFT
# party of the top 2.
#
# Script 17 established that a narrow anti-pluralist victory is strongly
# bundled with the economic right: at the margin the more anti-pluralist party
# is also the more right-wing one about two thirds of the time (binary RD
# -0.439, -0.75 SD continuous). That says a confound EXISTS. It cannot say
# whether the confound is what lowers growth.
#
# This can. Split the sample on the sign of the left-right gap between the two
# top-2 parties:
#
#   right  the more anti-pluralist party is also the more RIGHT-wing one.
#          A narrow anti-pluralist victory is also a narrow right victory, so
#          "anti-pluralism hurts growth" and "the economic right hurts growth"
#          predict the SAME sign here and cannot be told apart.
#   left   the more anti-pluralist party is the more LEFT-wing one. Now a
#          narrow anti-pluralist victory is simultaneously a narrow economic
#          LEFT victory, so the two explanations predict OPPOSITE signs.
#
# The left subset is the informative one. If growth still falls there, the
# effect is anti-pluralism rather than the economic right; if it reverses, the
# headline is about economic policy and the anti-pluralism framing is wrong.
#
# Everything else is the main spec untouched: the ~1,347-election sample, the
# anti-pluralism running variable, the election year excluded from the window.
# No new build is needed -- the split is computed from columns the build
# already carries.
#
#   Rscript --no-init-file scripts/22_econleft_split_rdd.R
#
# Output: output/runs/_sweeps/econleft_split_rdd_<instr><suffix>/
#           econleft_split_results.csv   subset x window x outcome
#           subset_counts.csv            sample accounting per subset x window
#           comparison_w5.html           the headline numbers side by side
#           econleft_split_growth.png    two panels, window sweep, 3 GDP series
# ==============================================================================

library(tidyverse)
library(rdrobust)
library(here)
library(gt)

source(here::here("scripts", "rdd_helpers.R"))

data_dir <- here::here("data")

# ---- toggles -----------------------------------------------------------------

if (!exists("SPLIT_INSTRUMENT")) SPLIT_INSTRUMENT <- "v2xpa_antiplural"
if (!exists("SPLIT_WINDOWS")) SPLIT_WINDOWS <- 1:10
if (!exists("SPLIT_TREATMENT")) SPLIT_TREATMENT <- "backsliding_Nyr"
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- FALSE
}
if (!exists("PLACEBO_PRE_WINDOW")) PLACEBO_PRE_WINDOW <- FALSE

# The dimension the sample is split on. Carried in the build as the NEGATED
# variable (see 11_build_rdd_data.R Step 2), so higher = more LEFT.
SPLIT_SCORE <- "v2pariglef_neg"

build_suffix <- paste0(
  if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "" else "_exclyr",
  if (PLACEBO_PRE_WINDOW) "_pre" else ""
)

out_dir <- sweep_dir(sprintf(
  "econleft_split_rdd_%s%s", INSTRUMENT_LABELS[[SPLIT_INSTRUMENT]], build_suffix
))

build_path <- function(n) {
  file.path(
    data_dir, "rdd_build",
    sprintf("rdd_%s_w%d%s.rds", SPLIT_INSTRUMENT, n, build_suffix)
  )
}
missing_builds <- SPLIT_WINDOWS[
  !file.exists(vapply(SPLIT_WINDOWS, build_path, character(1)))
]
if (length(missing_builds) > 0) {
  stop(
    "No build for window(s) ", paste(missing_builds, collapse = ", "),
    ". Run 11_build_rdd_data.R at those windows first.",
    call. = FALSE
  )
}

# ---- the split ---------------------------------------------------------------

# The build's <score>__winner / __loser columns are ordered by VOTE SHARE, not
# by anti-pluralism, so the more-anti-pluralist party's left-right score has to
# be recovered from the sign of the running variable: that party is the winner
# exactly when running_var > 0. (Same reconstruction as script 17's balance
# test. Reading __winner as "the anti-pluralist" would silently mix the two
# subsets together.)
#
# SIGN DISCIPLINE. v2pariglef_neg is the NEGATION of V-Party's right-positive
# v2pariglef, so higher = more LEFT. The more anti-pluralist party is therefore
# more RIGHT-wing when its negated score is LOWER. Getting this backwards
# inverts the entire result and would look perfectly plausible, so it is
# asserted below against the raw right-positive direction rather than trusted
# to this comment.
add_split <- function(d) {
  a <- ifelse(d$running_var > 0,
              d[[paste0(SPLIT_SCORE, "__winner")]],
              d[[paste0(SPLIT_SCORE, "__loser")]])
  b <- ifelse(d$running_var > 0,
              d[[paste0(SPLIT_SCORE, "__loser")]],
              d[[paste0(SPLIT_SCORE, "__winner")]])
  d$ap_leftscore <- a       # more anti-pluralist party, higher = more left
  d$other_leftscore <- b
  d$split <- dplyr::case_when(
    is.na(a) | is.na(b) ~ NA_character_,
    a < b ~ "right",        # anti-pluralist is the more RIGHT-wing of the two
    a > b ~ "left",         # anti-pluralist is the more LEFT-wing of the two
    TRUE ~ NA_character_    # exact tie: neither, and there are only a handful
  )
  d
}

SPLIT_LABELS <- c(
  right = "Anti-pluralist party is the more RIGHT-wing of the top 2",
  left = "Anti-pluralist party is the more LEFT-wing of the top 2"
)

outcome_vars <- ALL_OUTCOME_VARS

# ---- run ---------------------------------------------------------------------

results <- list()
counts <- list()

for (n in SPLIT_WINDOWS) {
  d <- add_split(readRDS(build_path(n)))

  n_tie <- sum(
    !is.na(d$ap_leftscore) & !is.na(d$other_leftscore) &
      d$ap_leftscore == d$other_leftscore
  )
  n_miss <- sum(is.na(d$ap_leftscore) | is.na(d$other_leftscore))
  n_right <- sum(d$split == "right", na.rm = TRUE)
  n_left <- sum(d$split == "left", na.rm = TRUE)

  # The sample must account for itself exactly.
  stopifnot(n_right + n_left + n_tie + n_miss == nrow(d))

  cat(sprintf(
    "\n[w=%d] %d elections = %d right + %d left + %d ties + %d missing\n",
    n, nrow(d), n_right, n_left, n_tie, n_miss
  ))

  for (sp in names(SPLIT_LABELS)) {
    dd <- d[which(d$split == sp), , drop = FALSE]

    # THE SIGN CHECK. Re-derived from the RAW right-positive scale by undoing
    # the negation, so it is independent of the case_when above: in the
    # "right" subset the anti-pluralist party's raw v2pariglef must be HIGHER
    # than its opponent's, and lower in "left". A sign error anywhere in the
    # reconstruction fails here rather than silently swapping the two panels.
    raw_ap <- -dd$ap_leftscore
    raw_other <- -dd$other_leftscore
    if (nrow(dd) > 0) {
      if (sp == "right") {
        stopifnot(all(raw_ap > raw_other))
      } else {
        stopifnot(all(raw_ap < raw_other))
      }
    }

    counts[[length(counts) + 1]] <- tibble(
      subset = sp, subset_label = unname(SPLIT_LABELS[sp]), window = n,
      n_elections = nrow(dd),
      n_treated = sum(dd[[SPLIT_TREATMENT]], na.rm = TRUE),
      mean_raw_leftright_antipluralist = mean(raw_ap, na.rm = TRUE),
      mean_raw_leftright_other = mean(raw_other, na.rm = TRUE),
      n_tie_excluded = n_tie, n_missing_excluded = n_miss
    )

    if (nrow(dd) < RD_MIN_OBS) next

    # First stage once per subsample, reused across outcomes -- the same
    # structure as run_spec() in 12_rdd_analysis.R. Unlike the decade x OECD
    # cells in script 21, these subsets are large enough (822 and 517) for the
    # fuzzy arm to be worth estimating.
    fs <- extract_rd(safe_rdrobust(dd[[SPLIT_TREATMENT]], dd$running_var))

    for (v in outcome_vars) {
      if (!v %in% names(dd)) next
      rf <- extract_rd(safe_rdrobust(dd[[v]], dd$running_var))
      late <- extract_rd(safe_rdrobust(
        dd[[v]], dd$running_var, fuzzy = dd[[SPLIT_TREATMENT]]
      ))
      results[[length(results) + 1]] <- tibble(
        subset = sp, subset_label = unname(SPLIT_LABELS[sp]), window = n,
        treatment = SPLIT_TREATMENT,
        outcome = v, outcome_label = outcome_full_label(v),
        panel = names(OUTCOME_PANELS)[
          vapply(OUTCOME_PANELS, function(p) v %in% names(p), logical(1))
        ][1],
        n_elections = nrow(dd),
        first_stage_coef = fs$coef, first_stage_se = fs$se,
        first_stage_pval = fs$pval, n_first_stage = fs$N,
        n_outcome = rf$N,
        rd_estimate = rf$coef, rd_se = rf$se, rd_pval = rf$pval,
        rd_ci_lo = rf$ci_lo, rd_ci_hi = rf$ci_hi, bandwidth = rf$bw,
        late_estimate = late$coef, late_se = late$se, late_pval = late$pval,
        late_ci_lo = late$ci_lo, late_ci_hi = late$ci_hi
      )
    }
  }
}

results <- bind_rows(results)
counts <- bind_rows(counts)

write_csv(results, file.path(out_dir, "econleft_split_results.csv"))
write_csv(counts, file.path(out_dir, "subset_counts.csv"))

write_sweep_config(
  out_dir,
  fixed = list(
    instrument = SPLIT_INSTRUMENT,
    treatment = SPLIT_TREATMENT,
    split_on = sprintf("%s (negated; higher = more LEFT)", SPLIT_SCORE),
    incl_election_year = TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR,
    placebo = PLACEBO_PRE_WINDOW,
    sample = "no restriction beyond the split"
  ),
  swept = list(
    subset = names(SPLIT_LABELS), window = SPLIT_WINDOWS, outcome = outcome_vars
  )
)

cat("\nSample accounting:\n")
print(
  counts |>
    filter(window == 5) |>
    transmute(
      subset, n_elections, n_treated,
      `mean raw L-R, anti-pluralist` = round(mean_raw_leftright_antipluralist, 3),
      `mean raw L-R, other` = round(mean_raw_leftright_other, 3)
    ) |>
    as.data.frame(),
  row.names = FALSE
)

# ---- headline table ----------------------------------------------------------

POOLED_RUN <- file.path(
  RUNS_ROOT,
  sprintf(
    "instr-%s_w5_trt-%s_gapany_illibany%s",
    INSTRUMENT_LABELS[[SPLIT_INSTRUMENT]], TREATMENT_LABELS[[SPLIT_TREATMENT]],
    build_suffix
  ),
  "rdd_results.csv"
)

w5 <- results |> filter(window == 5)
tbl <- w5 |>
  transmute(
    Outcome = outcome_label, Variable = outcome,
    subset, cell = fmt_est(rd_estimate, rd_se, rd_pval), n = n_outcome
  ) |>
  pivot_wider(names_from = subset, values_from = c(cell, n))

if (file.exists(POOLED_RUN)) {
  pooled <- read_csv(POOLED_RUN, show_col_types = FALSE) |>
    transmute(
      Variable = outcome,
      `Pooled` = fmt_est(rd_estimate, rd_se, rd_pval), n_pooled = n_outcome
    )
  tbl <- tbl |> left_join(pooled, by = "Variable")
}

save_table_html(
  tbl |>
    transmute(
      Outcome, Variable,
      `Pooled` = if ("Pooled" %in% names(tbl)) Pooled else NA_character_,
      `N pooled` = if ("n_pooled" %in% names(tbl)) n_pooled else NA_integer_,
      `Anti-pluralist more RIGHT` = cell_right, `N right` = n_right,
      `Anti-pluralist more LEFT` = cell_left, `N left` = n_left
    ),
  file.path(out_dir, "comparison_w5.html"),
  "Reduced-form RD split by the top-2 left-right ordering, w = 5",
  sprintf(
    "Instrument: %s | Treatment: %s | Window: [%s, election_year + 5] | Split on %s",
    INSTRUMENT_DISPLAY[[SPLIT_INSTRUMENT]], TREATMENT_DISPLAY[[SPLIT_TREATMENT]],
    if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "election_year" else "election_year + 1",
    "V-Party economic left-right (v2pariglef)"
  ),
  note = paste(
    SIG_FOOTNOTE,
    "In the LEFT subset a narrow anti-pluralist victory is simultaneously a",
    "narrow economic-LEFT victory, so 'anti-pluralism lowers growth' and 'the",
    "economic right lowers growth' predict OPPOSITE signs there; in the RIGHT",
    "subset they predict the same sign and cannot be separated.",
    sprintf(
      "%d elections are excluded from both subsets at w5 (%d exact ties, %d missing a score).",
      counts$n_tie_excluded[1] + counts$n_missing_excluded[1],
      counts$n_tie_excluded[1], counts$n_missing_excluded[1]
    )
  )
)

# ---- figure ------------------------------------------------------------------

GROWTH_VARS <- names(OUTCOME_PANELS$growth)

growth <- results |>
  filter(outcome %in% GROWTH_VARS, !is.na(rd_estimate)) |>
  mutate(
    series = factor(unname(OUTCOME_PANELS$growth[outcome]),
                    levels = unname(OUTCOME_PANELS$growth)),
    subset = factor(subset, levels = c("right", "left"))
  )

# Same ribbon clipping as 21_cell_rdd.R: a few wide CIs on a shared axis
# flatten every line. Estimates are never clipped and the unclipped bounds
# stay in the CSV.
rng <- range(growth$rd_estimate, na.rm = TRUE)
pad <- max(diff(rng) * 0.15, abs(rng[2]) * 0.1, 1e-9)
y_lim <- c(rng[1] - pad, rng[2] + pad)
n_clipped <- sum(
  growth$rd_ci_lo < y_lim[1] | growth$rd_ci_hi > y_lim[2], na.rm = TRUE
)
growth <- growth |>
  mutate(rd_ci_lo = pmax(rd_ci_lo, y_lim[1]), rd_ci_hi = pmin(rd_ci_hi, y_lim[2]))

strip_n <- counts |>
  filter(window == 5) |>
  mutate(strip = sprintf("%s\n(N = %d at w5)", unname(SPLIT_LABELS[subset]), n_elections))
growth <- growth |>
  left_join(strip_n |> select(subset, strip), by = "subset") |>
  mutate(strip = factor(strip, levels = strip_n$strip[match(c("right", "left"), strip_n$subset)]))

p <- ggplot(growth, aes(x = window, y = rd_estimate, colour = series)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.35) +
  geom_ribbon(
    aes(ymin = rd_ci_lo, ymax = rd_ci_hi, fill = series),
    alpha = 0.10, colour = NA
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(shape = series), size = 1.6) +
  facet_wrap(~strip, ncol = 2) +
  scale_colour_manual(values = SERIES_COLORS[1:3], name = NULL) +
  scale_fill_manual(values = SERIES_COLORS[1:3], name = NULL) +
  scale_shape_manual(values = SERIES_SHAPES[1:3], name = NULL) +
  scale_x_continuous(breaks = SPLIT_WINDOWS) +
  coord_cartesian(ylim = y_lim) +
  labs(
    title = "GDP per capita growth, split by which top-2 party is further right",
    subtitle = paste(
      strwrap(sprintf(paste(
        "Instrument: %s. Reduced-form RD at each window length. Decades and",
        "OECD status pooled. In the LEFT panel a narrow anti-pluralist victory",
        "is also a narrow economic-LEFT victory, so 'anti-pluralism lowers",
        "growth' and 'the economic right lowers growth' predict opposite",
        "signs; in the RIGHT panel they predict the same sign and cannot be",
        "told apart. Band = robust 95%% CI, shared y across both panels;",
        "%d ribbon bound(s) clipped to it, unclipped values in",
        "econleft_split_results.csv. %d elections excluded from both subsets",
        "at w5 (%d exact ties, %d missing a left-right score)."
      ), INSTRUMENT_DISPLAY[[SPLIT_INSTRUMENT]], n_clipped,
      counts$n_tie_excluded[1] + counts$n_missing_excluded[1],
      counts$n_tie_excluded[1], counts$n_missing_excluded[1]), 118),
      collapse = "\n"
    ),
    x = "Window length N (years after the election)",
    y = "Cumulative log change in GDP per capita at the cutoff"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(size = 8.5, face = "bold"),
    legend.position = "bottom",
    legend.margin = margin(t = -2),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 7.5, colour = "grey25")
  )

ggsave(
  file.path(out_dir, "econleft_split_growth.png"), p,
  width = 10, height = 5.4, dpi = 150
)
cat("Saved econleft_split_growth.png\n")

message("\nEcon-L-R split written to ", out_dir)
