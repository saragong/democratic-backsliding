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
#   binary      is the winner also the MORE populist / more left / more
#               anti-elite / less minority-friendly of the two?
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
# WHAT THESE ESTIMATES ARE, AND ARE NOT
#
# They are not causal, and the arithmetic makes that concrete. Which party is
# "the winner" flips AT the cutoff by construction: above it the winner is the
# more anti-pluralist party, below it the other one. So
#
#   binary jump      ~  2 * P(the more anti-pluralist party is also the
#                            higher-scoring one on this index) - 1
#   continuous jump  ~  E[score of the more anti-pluralist party
#                         - score of the other one]
#
# both evaluated at the cutoff. Verified numerically: the binary populism
# estimate is -0.176, and 2P-1 computed directly in shrinking windows around
# the cutoff runs -0.12 to -0.22.
#
# So each number is a COMPOSITION statistic about the elections at the margin,
# not an effect of anything. That is still the thing worth knowing here: it
# answers whether the treatment is separable from a correlated trait at the
# margin, which is the precondition for the headline growth result meaning
# what it says. If populism were perfectly bundled (P = 1, jump = +1) the
# design could not tell the two apart at all.
#
# What it does NOT do is apportion the growth effect. Bundling says a confound
# exists; it does not say the confound drives growth. The test that would is
# the sign-discordant subsample -- elections where the more anti-pluralist
# party is the LESS right-wing one, so the two explanations predict
# opposite-signed growth effects -- which belongs in the growth RD, not here.
#
# The window is irrelevant to every outcome above -- they are all measured at
# the election -- so the script reads a single build and says so.
#
#   Rscript --no-init-file scripts/17_party_outcomes_rdd.R
#
# Output: output/runs/_sweeps/party_outcome_rdd_<instr>_w<N><suffix>/
#           party_outcome_results.csv    one row per score x form x restriction
#           party_outcome_cells.csv      one row per score x form x decade x OECD
#           party_outcomes_table.html    the unrestricted sample, flat
#           grid_party_outcomes_binary.html      restriction x score colour
#           grid_party_outcomes_continuous.html  grids, one per outcome form
#                                        (continuous cells are in SD units --
#                                        see that section for why they are two
#                                        files rather than one)
#           grid_cells_<form>.html       decade x OECD x score colour grids
#           plots/party_outcomes_<form>.png       pooled, one panel per score
#           plots/cells_<form>_<score>.png        the 2 x 5 decade x OECD grid,
#                                        one figure per score x form
# ==============================================================================

library(tidyverse)
library(rdrobust)
library(here)
library(patchwork)
library(gt)

source(here::here("scripts", "rdd_helpers.R"))
# For oecd_group() and decade_label(). They live in vparty_helpers.R because
# the raw-V-Party scripts needed them first, but they are general splits on a
# country code and a year, not anything V-Party-specific.
source(here::here("scripts", "vparty_helpers.R"))

data_dir <- here::here("data")

# ---- toggles -----------------------------------------------------------------

if (!exists("ILLIBERALISM_VAR")) ILLIBERALISM_VAR <- "v2xpa_antiplural"
# Any window gives the same answer -- see the header -- so this only selects
# which build file to open.
if (!exists("BACKSLIDING_WINDOW_YEARS")) BACKSLIDING_WINDOW_YEARS <- 5
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- FALSE
}

# Every party score except two.
#
# ILLIBERALISM_VAR itself is excluded because including it would produce a
# mechanical result: illiberal_score is the max of the two by construction, so
# "is the winner the more anti-pluralist one" is exactly the treatment
# indicator and jumps from 0 to 1 at the cutoff.
#
# ep_galtan is excluded because it is too thinly coded to support this
# analysis. It comes from the CHES merge and both top-2 members carry a score
# in only about 4% of elections -- 53 of 1,347 pooled, and 0 to 21 usable
# observations per decade x OECD cell against a CELL_MIN_N of 40, so every one
# of the ten cells came back blank. The pooled estimate it did produce was
# insignificant on n = 53 with a bandwidth less than half the others'.
# Carrying it through only added an empty column to every grid and two blank
# figures. Put it back by setting COMPARISON_SCORES explicitly.
if (!exists("COMPARISON_SCORES")) {
  COMPARISON_SCORES <- setdiff(PARTY_SCORE_VARS, c(ILLIBERALISM_VAR, "ep_galtan"))
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

# ---- decade x OECD cells -----------------------------------------------------
#
# The same 2 x 5 grid scripts 18 and 19 use, so the three figures are read
# together. Bounded to VPARTY_YEAR_MIN..MAX (1970-2019) for the same reason:
# it is exactly five decades, and it is the window where V-Party coverage is
# not thin and lopsided. 87 of the 1,347 elections fall outside it and are in
# the pooled sample but not in any cell -- the figures say so.
if (!exists("CELL_YEAR_MIN")) CELL_YEAR_MIN <- VPARTY_YEAR_MIN
if (!exists("CELL_YEAR_MAX")) CELL_YEAR_MAX <- VPARTY_YEAR_MAX

# A cell below this many usable observations is reported in the CSV but not
# fitted or drawn. The thinnest cell (non-OECD 1970s) has 66 elections and only
# ~21 within a typical bandwidth, which is at the edge of what rdrobust will
# attempt and well past what it will do reliably.
#
# Counted PER OUTCOME, not per cell. Coverage differs by score, so a large
# cell can still hold almost no usable observations for one of them -- which
# is exactly what ruled ep_galtan out of COMPARISON_SCORES above. Gating on
# cell size alone would hand rdrobust a 3-observation sample and draw a panel
# around whatever came back. The remaining scores all clear the bar in every
# cell, so this is now a guard rather than an active filter, but it is what
# makes adding a thinly-coded score back safe.
if (!exists("CELL_MIN_N")) CELL_MIN_N <- 40

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
  # higher" is not "the winner was not higher". Coding it FALSE would silently
  # push every unscored pair into the "winner does not score higher" bin and
  # bias the level of the binary outcome downwards.
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
    n = fit$N, est = fit$coef, se = fit$se, pval = fit$pval,
    ci_lo = fit$ci_lo, ci_hi = fit$ci_hi, bw = fit$bw
  )
}

# One restriction level: resolve against the FULL build (never the already-cut
# frame -- otherwise the levels would not be the quantiles they claim to be),
# apply, estimate every score x form.
#
# apply_threshold() rather than a hand-rolled filter, so these levels behave
# exactly like every other restricted run in the pipeline and print what each
# one cost. It also gets the no-restriction case right by construction: "none"
# returns the data UNTOUCHED, whereas filtering on `>= -Inf` would quietly drop
# the rows where the axis variable is NA -- which for score_gap_z is every
# country with too few party-years to compute its own SD.
run_level <- function(axis_var, spec) {
  thr <- resolve_threshold(
    parse_threshold(spec, axis_var), d_full[[axis_var]], axis_var
  )
  cat(sprintf("\n[%s = %s]\n", axis_var, thr$spec))
  dd <- apply_threshold(d_full, axis_var, thr, axis_var, op = ">=")

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
  "ep_galtan is excluded: both top-2 members carry a CHES GAL-TAN score in",
  "only about 4% of elections, too few to fit per decade x OECD cell."
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

# ---- restriction x score grids ----------------------------------------------
#
# ONE FILE PER OUTCOME FORM, not one file with both. Two reasons, both about
# color_grid_blocks():
#
#   1. It groups rows with gt(groupname_col = "section") and builds `section`
#      from row_axis and col_axis alone -- it ignores each block's `label`.
#      Four blocks that share a row/column axis therefore collapse into ONE
#      row group with repeated row labels and nothing saying which is which.
#      So the gap axis goes INTO row_axis, making the two sections distinct.
#   2. It normalizes the colour scale over every block in the file ("one
#      symmetric scale shared by every section", per its own footnote). The
#      binary estimates are probabilities (|est| <= 0.68); the continuous ones
#      include expert-scale indices reaching 2.14. Sharing a scale washes the
#      binary cells out to near-white, so a 0.68 jump in probability -- which
#      is enormous -- reads as weaker than a trivial expert-scale cell.
#
# Within the continuous form the columns are STILL on different scales
# (v2xpa_popul on [0, 1], v2pariglef_neg and v2paanteli on expert scales
# running about -2 to +4, v2paminor likewise), so its cells are reported in SD
# units of each score's own distribution across top-2 parties. That makes both
# the number and the colour comparable across columns; the raw estimates stay
# in the CSV and in party_outcomes_table.html.
# ------------------------------------------------------------------------------

# Pooled SD of each score over both top-2 members of every scored election --
# the same grain the scores are measured at.
score_sd <- vapply(COMPARISON_SCORES, function(s) {
  sd(
    c(d_full[[paste0(s, "__winner")]], d_full[[paste0(s, "__loser")]]),
    na.rm = TRUE
  )
}, numeric(1))

# One color_grid_blocks() block: rows are whatever `row_labels` names, columns
# are the comparison scores. Used by both the restriction grids and the
# decade x OECD grids, so the SD rescaling and the column ordering are defined
# once.
#
# `row_labels` is a NAMED vector: names are the keys to look up in `df`, values
# are what the table shows.
make_block <- function(df, key_col, row_labels, row_axis, fm) {
  keys <- names(row_labels)
  rescale <- identical(fm, "continuous")
  mat <- function(col, do_rescale) {
    m <- matrix(
      NA_real_, nrow = length(keys), ncol = length(COMPARISON_SCORES),
      dimnames = list(
        unname(row_labels), unname(PARTY_SCORE_LABELS[COMPARISON_SCORES])
      )
    )
    for (i in seq_along(keys)) {
      for (j in seq_along(COMPARISON_SCORES)) {
        v <- df[df[[key_col]] == keys[i] & df$score == COMPARISON_SCORES[j], ]
        if (nrow(v) == 1 && !is.na(v[[col]])) {
          m[i, j] <- if (do_rescale) {
            v[[col]] / score_sd[[COMPARISON_SCORES[j]]]
          } else {
            v[[col]]
          }
        }
      }
    }
    m
  }
  list(
    # `label` is never rendered by color_grid_blocks() -- it names its gt row
    # group from row_axis and col_axis -- so anything that must distinguish
    # two blocks has to go in row_axis.
    label = row_axis,
    row_axis = row_axis,
    col_axis = "Comparison score",
    row_labels = unname(row_labels),
    col_labels = unname(PARTY_SCORE_LABELS[COMPARISON_SCORES]),
    est = mat("est", rescale),
    se = mat("se", rescale),
    pval = mat("pval", FALSE),
    n = mat("n", FALSE)
  )
}

# The unit note that must travel with any grid of these cells.
unit_note_for <- function(fm) {
  if (identical(fm, "continuous")) {
    paste(
      "Cells are in SD units of each score's own distribution across top-2",
      "parties, so columns on different native scales are comparable.",
      "Raw estimates are in the CSVs and in party_outcomes_table.html."
    )
  } else {
    "Cells are changes in probability, so all columns share a scale."
  }
}

for (fm in names(FORMS)) {
  blocks <- list()
  for (ax in names(GAP_AXES)) {
    sub <- results |> filter(axis == ax, form == fm)
    if (nrow(sub) == 0) next
    lvls <- unique(sub$level)
    row_labels <- setNames(
      vapply(lvls, function(l) {
        r <- sub |> filter(level == l) |> slice(1)
        if (l == "-Inf") {
          sprintf("all (N = %d)", r$n_sample)
        } else {
          sprintf("%s: >= %.3g (N = %d)", l, r$cutoff, r$n_sample)
        }
      }, character(1)),
      lvls
    )
    blocks[[length(blocks) + 1]] <- make_block(
      sub, "level", row_labels, unname(GAP_AXES[ax]), fm
    )
  }
  if (length(blocks) == 0) next

  color_grid_blocks(
    blocks,
    path = file.path(out_dir, sprintf("grid_party_outcomes_%s.html", fm)),
    title = sprintf(
      "Reduced-form RD on the winner's other party scores -- %s",
      FORMS[[fm]]$title
    ),
    subtitle = subtitle_base,
    note = paste(
      SIG_FOOTNOTE, unit_note_for(fm), REDUNDANCY_NOTE,
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
    # Same convention as the cell figures: the estimate and its robust CI go
    # in the panel title, so the figure can be read without the table.
    e <- results |>
      filter(level == "-Inf", axis == names(GAP_AXES)[1], form == fm, score == s)
    ttl <- sprintf(
      "%s%s",
      unname(PARTY_SCORE_DISPLAY[s]),
      if (nrow(e) == 1 && !is.na(e$est)) {
        sprintf("\nRD %s", fmt_est_ci(e$est, e$ci_lo, e$ci_hi, e$pval))
      } else {
        ""
      }
    )
    p <- build_panel_plot(
      d_full,
      setNames(unname(PARTY_SCORE_LABELS[s]), v),
      ttl,
      spec$ylab,
      y_lim = spec$ylim
    )
    if (!is.null(p)) panels[[s]] <- p
  }
  if (length(panels) == 0) next
  # 3 across reads well up to 6 panels, but leaves a lone panel stranded on a
  # second row at exactly 4 -- which is the current count. Square it instead.
  n_col <- if (length(panels) == 4) 2 else min(3, length(panels))
  combined <- wrap_plots(panels, ncol = n_col) +
    plot_annotation(
      title = sprintf("%s (%s)", spec$title, INSTRUMENT_DISPLAY[[ILLIBERALISM_VAR]]),
      subtitle = paste(
        "All scored elections. Panel titles give the RD estimate and its",
        "robust 95% CI. Shaded band = 95% CI of the LOCAL-LINEAR FIT",
        "(conventional, not bias-corrected), which is a different thing."
      )
    )
  n_rows <- ceiling(length(panels) / n_col)
  ggsave(
    file.path(plots_dir, sprintf("party_outcomes_%s.png", fm)),
    combined, width = 5 * n_col, height = 3.8 * n_rows + 0.6, dpi = 150,
    limitsize = FALSE
  )
  cat(sprintf("Saved plots/party_outcomes_%s.png\n", fm))
}

# ==============================================================================
# Decade x OECD cells
#
# The same 2 x 5 grid as 18_vparty_jaccard_panels.R and
# 19_vparty_ideology_quadrants.R, so all three are read together: rows are
# OECD / non-OECD, columns are the 1970s to the 2010s.
#
# The question these answer is whether the bundling is stable. The pooled
# numbers say anti-pluralism is near-orthogonal to populism but strongly
# bundled with the economic right; if that is a recent, OECD-specific
# phenomenon (which is what the rising populism/anti-pluralism correlation in
# adhoc/vparty_corr_by_decade.R would suggest) then the pooled figure is an
# average over regimes that differ, and the confound is worse in exactly the
# cells the headline result leans on.
#
# FULL SAMPLE ONLY -- no restriction sweep is crossed with the cells. Ten cells
# times five restriction levels times two axes would be 100 subsamples of a
# 1,347-election build, most of them too thin to fit.
#
# Every cell gets its OWN bandwidth and binning, which is why these are ten
# separate panels assembled with patchwork rather than a facet_grid: a shared
# binning across cells with very different sample sizes would misrepresent all
# of them. The y axis IS shared within a figure, so the panels stay comparable.
# ==============================================================================

d_cells <- d_full |>
  filter(
    election_year >= CELL_YEAR_MIN,
    election_year <= CELL_YEAR_MAX
  ) |>
  mutate(
    group = oecd_group(country_text_id),
    decade = decade_label(election_year)
  )

n_dropped <- nrow(d_full) - nrow(d_cells)
cat(sprintf(
  "\nDecade x OECD cells: %d of %d elections fall in %d-%d (%d outside, pooled only)\n",
  nrow(d_cells), nrow(d_full), CELL_YEAR_MIN, CELL_YEAR_MAX, n_dropped
))

cell_keys <- expand_grid(
  group = levels(d_cells$group),
  decade = levels(d_cells$decade)
) |>
  mutate(
    group = factor(group, levels = levels(d_cells$group)),
    decade = factor(decade, levels = levels(d_cells$decade)),
    cell = paste(group, decade),
    n_cell = map2_int(group, decade, \(g, dd) {
      sum(d_cells$group == g & d_cells$decade == dd)
    })
  )

cell_df <- function(g, dd) {
  d_cells[d_cells$group == g & d_cells$decade == dd, , drop = FALSE]
}

# Usable observations for ONE outcome in ONE cell -- the number the fit
# actually rests on, and what CELL_MIN_N is compared against.
n_usable <- function(g, dd, var) sum(!is.na(cell_df(g, dd)[[var]]))

cat("\nCell sizes (elections, and usable observations per outcome):\n")
print(
  cell_keys |>
    select(group, decade, n_cell) |>
    mutate(!!!setNames(
      lapply(COMPARISON_SCORES, function(sc) {
        map2_int(cell_keys$group, cell_keys$decade,
                 \(g, dd) n_usable(g, dd, paste0("Zbin_", sc)))
      }),
      unname(PARTY_SCORE_LABELS[COMPARISON_SCORES])
    )) |>
    as.data.frame(),
  row.names = FALSE
)
cat(sprintf("(a cell is fitted when its usable count reaches %d)\n", CELL_MIN_N))

# ---- estimation --------------------------------------------------------------

cell_results <- pmap_dfr(
  cell_keys,
  function(group, decade, cell, n_cell) {
    dd <- cell_df(group, decade)
    expand_grid(form = names(FORMS), score = COMPARISON_SCORES) |>
      pmap_dfr(function(form, score) {
        var <- paste0(FORMS[[form]]$prefix, score)
        nu <- sum(!is.na(dd[[var]]))
        base <- tibble(
          group = group, decade = decade, cell = cell,
          n_cell = n_cell, n_usable = nu, fitted = nu >= CELL_MIN_N,
          form = form, score = score,
          score_label = unname(PARTY_SCORE_DISPLAY[score])
        )
        # A too-thin cell is still emitted, so it is visibly blank in the CSV
        # rather than silently absent from it.
        if (nu < CELL_MIN_N) {
          return(bind_cols(base, tibble(
            n = NA_integer_, est = NA_real_, se = NA_real_, pval = NA_real_,
            ci_lo = NA_real_, ci_hi = NA_real_, bw = NA_real_
          )))
        }
        bind_cols(base, estimate_cell(dd, var))
      })
  }
)

write_csv(cell_results, file.path(out_dir, "party_outcome_cells.csv"))

# ---- grids -------------------------------------------------------------------

# One file per form; one block per OECD group, so the two blocks get distinct
# gt row-group names (see make_block).
for (fm in names(FORMS)) {
  blocks <- list()
  for (g in levels(d_cells$group)) {
    sub <- cell_results |> filter(form == fm, group == g)
    if (nrow(sub) == 0) next
    decs <- levels(d_cells$decade)
    row_labels <- setNames(
      vapply(decs, function(dd) {
        k <- cell_keys |> filter(group == g, decade == dd)
        nf <- sum(sub$fitted[sub$decade == dd])
        sprintf(
          "%s (N = %d)%s", dd, k$n_cell,
          if (nf == 0) " -- all too thin" else ""
        )
      }, character(1)),
      decs
    )
    blocks[[length(blocks) + 1]] <- make_block(
      sub |> mutate(decade = as.character(decade)),
      "decade", row_labels, g, fm
    )
  }
  if (length(blocks) == 0) next

  color_grid_blocks(
    blocks,
    path = file.path(out_dir, sprintf("grid_cells_%s.html", fm)),
    title = sprintf(
      "Decade x OECD: RD on the winner's other party scores -- %s",
      FORMS[[fm]]$title
    ),
    subtitle = sprintf(
      "%s | Full sample within each cell, no restriction applied | %d-%d",
      subtitle_base, CELL_YEAR_MIN, CELL_YEAR_MAX
    ),
    note = paste(
      SIG_FOOTNOTE, unit_note_for(fm), REDUNDANCY_NOTE,
      sprintf(
        "Cells with fewer than %d elections are left blank rather than fitted.",
        CELL_MIN_N
      )
    )
  )
}

# ---- figures -----------------------------------------------------------------

# A stand-in so a thin cell keeps its slot and the 2 x 5 grid stays aligned,
# rather than the remaining panels reflowing and the decades falling out of
# line between the two rows.
blank_panel <- function(title, msg) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = msg, size = 2.6, colour = "grey45") +
    labs(title = title) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0))
}

for (fm in names(FORMS)) {
  spec <- FORMS[[fm]]
  for (sc in COMPARISON_SCORES) {
    var <- paste0(spec$prefix, sc)

    # One y range for all ten panels, so a tag moving between decades is
    # genuinely moving. Binary is already a probability on [0, 1]; for the
    # continuous form the 2nd-98th percentile of the score keeps the axis off
    # the tails without clipping the binned means.
    y_lim <- if (!is.null(spec$ylim)) {
      spec$ylim
    } else {
      q <- quantile(d_cells[[var]], c(0.02, 0.98), na.rm = TRUE)
      pad <- max(diff(q) * 0.12, 1e-6)
      unname(c(q[1] - pad, q[2] + pad))
    }

    n_dec <- nlevels(d_cells$decade)
    panels <- pmap(cell_keys, function(group, decade, cell, n_cell) {
      i <- which(levels(d_cells$decade) == decade)
      is_left <- i == 1
      is_bottom <- group == levels(d_cells$group)[nlevels(d_cells$group)]
      nu <- n_usable(group, decade, var)
      est <- cell_results |>
        filter(group == !!group, decade == !!decade, form == fm, score == sc)
      # The estimate and its interval go in the TITLE rather than inside the
      # panel. These cells are thin enough that their CI ribbons fill the
      # plotting area, so a reader left to eyeball the jump reads noise; and a
      # title is the one place guaranteed not to collide with the data.
      title <- sprintf(
        "%s \u00b7 %s (n = %d)%s",
        group, decade, nu,
        if (nrow(est) == 1 && !is.na(est$est)) {
          sprintf(
            "\nRD %s",
            fmt_est_ci(est$est, est$ci_lo, est$ci_hi, est$pval)
          )
        } else {
          ""
        }
      )

      # Axis titles only on the outer edge. Repeating the same two strings ten
      # times crowds the panels and, at this width, truncates them.
      trim_axes <- function(p) {
        p + labs(
          x = if (is_bottom) "Running variable (illiberal - other share, pp)" else NULL,
          y = if (is_left) spec$ylab else NULL
        )
      }

      if (nu < CELL_MIN_N) {
        return(blank_panel(
          title, sprintf("%d usable, need %d", nu, CELL_MIN_N)
        ))
      }
      p <- build_panel_plot(
        cell_df(group, decade),
        setNames(unname(PARTY_SCORE_LABELS[sc]), var),
        title, spec$ylab, y_lim = y_lim, show_legend = FALSE
      )
      if (is.null(p)) {
        return(blank_panel(title, "too few near the cutoff"))
      }
      trim_axes(p)
    })

    combined <- wrap_plots(panels, ncol = nlevels(d_cells$decade)) +
      plot_annotation(
        title = sprintf(
          "%s by decade and OECD -- %s",
          unname(PARTY_SCORE_DISPLAY[sc]),
          if (identical(fm, "binary")) {
            "is the winner the higher-scoring of the top 2?"
          } else {
            "the winner's own score"
          }
        ),
        subtitle = paste(
          strwrap(sprintf(paste(
            "Instrument: %s. Each panel is its own RD with its own bandwidth",
            "and binning; the y axis is shared across panels so they are",
            "comparable. Panel titles give the RD estimate and its robust 95%%",
            "CI. Shaded band = 95%% CI of the LOCAL-LINEAR FIT (conventional,",
            "not bias-corrected), which is a different thing and is why the",
            "two can disagree. Full sample within each cell.",
            "%d-%d; %d elections outside that window are in the pooled figures",
            "only."
          ), INSTRUMENT_DISPLAY[[ILLIBERALISM_VAR]],
          CELL_YEAR_MIN, CELL_YEAR_MAX, n_dropped), 130),
          collapse = "\n"
        )
      )
    ggsave(
      file.path(plots_dir, sprintf("cells_%s_%s.png", fm, PARTY_SCORE_LABELS[sc])),
      combined,
      width = 3.1 * nlevels(d_cells$decade) + 0.6,
      height = 2.9 * nlevels(d_cells$group) + 1.5,
      dpi = 150, limitsize = FALSE
    )
    cat(sprintf("Saved plots/cells_%s_%s.png\n", fm, PARTY_SCORE_LABELS[sc]))
  }
}

cat("\n=== Unrestricted sample ===\n")
print(as.data.frame(flat), row.names = FALSE)
message("\nParty-outcome RDs written to ", out_dir)
