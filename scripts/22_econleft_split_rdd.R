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
# Everything else is the main spec untouched: the anti-pluralism running
# variable, the election year excluded from the window. No new build is needed
# -- the split is computed from columns the build already carries.
#
# The split composes with the sample restrictions, which take the same spec
# vocabulary 12_rdd_analysis.R accepts, so the same two panels can be drawn on
# the full ~1,347-election sample or on the ~411 elections where one top-2
# party is illiberal and the other is not:
#
#   Rscript --no-init-file scripts/22_econleft_split_rdd.R
#   SPLIT_SAMPLE=popucut Rscript --no-init-file scripts/22_econleft_split_rdd.R
#
# Output: output/runs/_sweeps/econleft_split_rdd_<instr><restriction><suffix>/
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

# The sample the split is taken WITHIN. A (score_gap, illiberal, other) triple
# in the spec vocabulary of 12_rdd_analysis.R, so "popucut" resolves out of
# data/populist_threshold.rds exactly as it does there. The named presets are
# the two samples 21_cell_rdd.R uses, under the same two names; the three
# SPLIT_* thresholds can also be set directly for anything else.
SPLIT_SAMPLES <- list(
  full = list(score_gap_min = -Inf, illiberal_cutoff = -Inf, other_cutoff_max = Inf),
  popucut = list(
    score_gap_min = -Inf,
    illiberal_cutoff = "popucut", other_cutoff_max = "popucut"
  )
)
if (!exists("SPLIT_SAMPLE")) SPLIT_SAMPLE <- Sys.getenv("SPLIT_SAMPLE", "full")
if (!SPLIT_SAMPLE %in% names(SPLIT_SAMPLES)) {
  stop(
    'SPLIT_SAMPLE = "', SPLIT_SAMPLE, '" is not one of: ',
    paste(names(SPLIT_SAMPLES), collapse = ", "), ".",
    call. = FALSE
  )
}
if (!exists("SPLIT_SCORE_GAP_MIN")) {
  SPLIT_SCORE_GAP_MIN <- SPLIT_SAMPLES[[SPLIT_SAMPLE]]$score_gap_min
}
if (!exists("SPLIT_ILLIBERAL_CUTOFF")) {
  SPLIT_ILLIBERAL_CUTOFF <- SPLIT_SAMPLES[[SPLIT_SAMPLE]]$illiberal_cutoff
}
if (!exists("SPLIT_OTHER_CUTOFF_MAX")) {
  SPLIT_OTHER_CUTOFF_MAX <- SPLIT_SAMPLES[[SPLIT_SAMPLE]]$other_cutoff_max
}

build_suffix <- paste0(
  if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "" else "_exclyr",
  if (PLACEBO_PRE_WINDOW) "_pre" else ""
)

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

# ---- the sample restriction --------------------------------------------------

# Resolve the three specs to NUMBERS once, against the w5 build, BEFORE any
# filter bites -- the same discipline, and the same arbitrary-but-fixed
# reference build, as 14_window_sweep.R. Two things downstream need a number
# rather than a spec: this sweep's folder name via fmt_slug_num(), and the cfg
# whose run_slug() names the pooled run the comparison table reads against.
# Resolving per window instead would let a quantile spec mean a different
# sample in each panel of one figure.
local({
  ref_path <- build_path(5L)
  if (!file.exists(ref_path)) {
    stop(
      "Cannot resolve the sample restrictions: no reference build at ",
      ref_path, ". Build w5 first.",
      call. = FALSE
    )
  }
  ref <- readRDS(ref_path)
  SPLIT_SCORE_GAP_MIN <<- resolve_threshold_abs(
    SPLIT_SCORE_GAP_MIN, "SPLIT_SCORE_GAP_MIN", ref$score_gap_z
  )
  SPLIT_ILLIBERAL_CUTOFF <<- resolve_threshold_abs(
    SPLIT_ILLIBERAL_CUTOFF, "SPLIT_ILLIBERAL_CUTOFF", ref$illiberal_score
  )
  SPLIT_OTHER_CUTOFF_MAX <<- resolve_threshold_abs(
    SPLIT_OTHER_CUTOFF_MAX, "SPLIT_OTHER_CUTOFF_MAX", ref$other_score,
    none_value = Inf
  )
})
stopifnot(
  is.numeric(SPLIT_SCORE_GAP_MIN),
  is.numeric(SPLIT_ILLIBERAL_CUTOFF),
  is.numeric(SPLIT_OTHER_CUTOFF_MAX)
)

# Now that the three are plain numbers, parsing them again yields only "none"
# (non-finite) or "absolute", never "quantile" or "external" -- which is why
# resolve_threshold() can be handed NULL values here. It is called purely for
# the human-readable label apply_threshold() prints.
RESTRICTIONS <- list(
  list(var = "score_gap_z", name = "score_gap_min", op = ">=", thr = resolve_threshold(
    parse_threshold(SPLIT_SCORE_GAP_MIN, "score_gap_min"), NULL
  )),
  list(var = "illiberal_score", name = "illiberal_cutoff", op = ">", thr = resolve_threshold(
    parse_threshold(SPLIT_ILLIBERAL_CUTOFF, "illiberal_cutoff"), NULL
  )),
  list(var = "other_score", name = "other_cutoff_max", op = "<=", thr = resolve_threshold(
    parse_threshold(SPLIT_OTHER_CUTOFF_MAX, "other_cutoff_max", none_value = Inf), NULL
  ))
)

# The restriction, applied with the same three calls 12_rdd_analysis.R makes,
# in the same order and with the same operators. A hand-rolled filter here
# would be one more place for the pair condition to drift out of agreement
# with the main spec -- the mistake already caught once in 17's run_level().
apply_restrictions <- function(d) {
  for (r in RESTRICTIONS) d <- apply_threshold(d, r$var, r$thr, r$name, op = r$op)
  d
}

restriction_label <- restriction_sentence(
  SPLIT_SCORE_GAP_MIN, SPLIT_ILLIBERAL_CUTOFF, SPLIT_OTHER_CUTOFF_MAX,
  none = "no restriction beyond the split"
)

# Each active axis contributes a slug part, and only an active one: at the
# no-op values this reproduces the unrestricted run's folder name exactly, so
# the run already on disk is not orphaned by adding the restriction axis. Same
# rule, and the same spelling, as run_slug().
# `else ""` is load-bearing: a bare `if` with no else yields NULL, and
# paste0() of three NULLs is character(0), not "" -- sprintf() would then
# return character(0) and the folder name would vanish on the full sample.
restriction_slug <- paste0(
  if (is.finite(SPLIT_SCORE_GAP_MIN)) paste0("_gap", fmt_slug_num(SPLIT_SCORE_GAP_MIN)) else "",
  if (is.finite(SPLIT_ILLIBERAL_CUTOFF)) paste0("_illib", fmt_slug_num(SPLIT_ILLIBERAL_CUTOFF)) else "",
  if (is.finite(SPLIT_OTHER_CUTOFF_MAX)) paste0("_opp", fmt_slug_num(SPLIT_OTHER_CUTOFF_MAX)) else ""
)

out_dir <- sweep_dir(sprintf(
  "econleft_split_rdd_%s%s%s",
  INSTRUMENT_LABELS[[SPLIT_INSTRUMENT]], restriction_slug, build_suffix
))
cat("Sample within which the split is taken: ", restriction_label, "\n", sep = "")

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
  cat(sprintf("\n[w=%d] sample restrictions:\n", n))
  # Restrict first, split second. The other order would compute the tie and
  # missing counts on elections the sample does not contain, so the accounting
  # printed below would not add up to the sample actually estimated on.
  d <- add_split(apply_restrictions(readRDS(build_path(n))))

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
    "[w=%d] %d elections = %d right + %d left + %d ties + %d missing\n",
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
    score_gap_min = SPLIT_SCORE_GAP_MIN,
    illiberal_cutoff = SPLIT_ILLIBERAL_CUTOFF,
    other_cutoff_max = SPLIT_OTHER_CUTOFF_MAX,
    sample = restriction_label
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

# The pooled column must be the SAME sample, unsplit -- not the unrestricted
# main run. Naming it with run_slug() off the resolved numbers rather than
# writing the slug out by hand is what guarantees that: a restriction that
# changes the folder here changes the folder read there too, so the comparison
# cannot silently become apples-to-oranges. Missing folder = no pooled column,
# which is why this is a file.exists() check and not a stop().
POOLED_RUN <- file.path(
  RUNS_ROOT,
  run_slug(list(
    instrument = SPLIT_INSTRUMENT, window = 5L, treatment = SPLIT_TREATMENT,
    score_gap_min = SPLIT_SCORE_GAP_MIN,
    illiberal_cutoff = SPLIT_ILLIBERAL_CUTOFF,
    other_cutoff_max = SPLIT_OTHER_CUTOFF_MAX,
    incl_election_year = TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR,
    placebo = PLACEBO_PRE_WINDOW
  )),
  "rdd_results.csv"
)
if (!file.exists(POOLED_RUN)) {
  message(
    "No pooled run at ", POOLED_RUN,
    " -- comparison_w5.html will omit the pooled column. Run ",
    "14_window_sweep.R under the same restriction to produce it."
  )
}

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
    "Instrument: %s | Treatment: %s | Window: [%s, election_year + 5] | Split on %s | Sample: %s",
    INSTRUMENT_DISPLAY[[SPLIT_INSTRUMENT]], TREATMENT_DISPLAY[[SPLIT_TREATMENT]],
    if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "election_year" else "election_year + 1",
    "V-Party economic left-right (v2pariglef)", restriction_label
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
        "Instrument: %s. Sample: %s. Reduced-form RD at each window length.",
        "Decades and OECD status pooled. In the LEFT panel a narrow anti-pluralist victory",
        "is also a narrow economic-LEFT victory, so 'anti-pluralism lowers",
        "growth' and 'the economic right lowers growth' predict opposite",
        "signs; in the RIGHT panel they predict the same sign and cannot be",
        "told apart. Band = robust 95%% CI, shared y across both panels;",
        "%d ribbon bound(s) clipped to it, unclipped values in",
        "econleft_split_results.csv. %d elections excluded from both subsets",
        "at w5 (%d exact ties, %d missing a left-right score)."
      ), INSTRUMENT_DISPLAY[[SPLIT_INSTRUMENT]], restriction_label, n_clipped,
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
