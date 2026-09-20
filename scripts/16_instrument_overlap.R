# ==============================================================================
# How much do the candidate instrument definitions actually differ?
#
# Every alternative instrument (populism, economic left, anti-elitism, GAL-TAN)
# is supposed to pick out a different "illiberal" side of the top 2 than
# anti-pluralism does. If they mostly agree, then to-do 6's alternative
# instruments aren't really alternatives and their first stages aren't
# independent evidence. This script measures the agreement two ways.
#
#   1. UpSet plots over the sets {elections where the WINNER scored higher than
#      the loser on instrument s}, one set per instrument. Run for all
#      elections and for narrow elections only.
#
#      Winner-based is deliberate here and differs from the RD. The RD never
#      conditions on who won (that would select on treatment); but the question
#      "do these measures classify the same elections the same way" is a
#      descriptive measurement question, and the winner/loser split is the
#      natural, instrument-independent way to pose it.
#
#   2. Jaccard heatmaps at the top-2-finisher level, comparing anti-pluralism
#      deciles against populism deciles. Cell (i, j) is
#      |D_illib = i AND D_popul = j| / |D_illib = i OR D_popul = j|, so a
#      perfectly redundant pair of measures would light up the diagonal and
#      nothing else.
#
# "Narrow" is reported at two thresholds side by side, since the answer
# shouldn't hinge on where the line is drawn:
#   - a fixed +/-5 percentage-point margin
#   - the MSE-optimal bandwidth rdrobust picks for the main first stage, i.e.
#     exactly the elections the headline estimate is computed from
#
# Output: output/runs/_sweeps/instrument_overlap/
# ==============================================================================

library(tidyverse)
library(here)
library(gt)
library(patchwork)

source(here::here("scripts", "rdd_helpers.R"))
source(here::here("scripts", "vparty_helpers.R"))

data_dir <- here::here("data")

if (!exists("OVERLAP_INSTRUMENT_BUILD")) {
  OVERLAP_INSTRUMENT_BUILD <- "v2xpa_antiplural"
}
if (!exists("OVERLAP_WINDOW")) {
  OVERLAP_WINDOW <- 5
}
if (!exists("NARROW_MARGIN_PP")) {
  NARROW_MARGIN_PP <- 5
}
# Selects which build to read and is echoed into the sweep folder name. The
# instrument comparison itself doesn't depend on the treatment window, but the
# NARROW-election bandwidth does (it comes from the first stage), so the two
# conventions get separate folders rather than overwriting each other.
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- TRUE
}
build_suffix <- if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "" else "_exclyr"

out_dir <- sweep_dir(paste0("instrument_overlap", build_suffix))

build_path <- file.path(
  data_dir,
  "rdd_build",
  sprintf(
    "rdd_%s_w%d%s.rds",
    OVERLAP_INSTRUMENT_BUILD,
    OVERLAP_WINDOW,
    build_suffix
  )
)
parties_path <- sub("\\.rds$", "_parties.rds", build_path)
stopifnot(file.exists(build_path), file.exists(parties_path))

d <- readRDS(build_path)
parties <- readRDS(parties_path)
build_incl <- attr(d, "includes_election_year") %||% FALSE
if (!identical(build_incl, TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR)) {
  stop(
    "Build at ",
    build_path,
    " was made with ",
    "TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR = ",
    build_incl,
    " but this sweep asked for ",
    TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR,
    "."
  )
}

# ------------------------------------------------------------------------------
# Which instruments to include
#
# ep_galtan comes from the CHES bridge, which covers only European parties in
# recent decades: it has both top-2 scores for 53 of 1,347 elections here. Put
# on an UpSet alongside the others it would force the whole plot down to those
# 53 complete cases and make every intersection count meaningless. It is
# reported separately in the coverage table instead.
# ------------------------------------------------------------------------------

coverage <- map_dfr(ILLIBERALISM_VARS, function(v) {
  ok <- !is.na(d[[paste0(v, "__winner")]]) & !is.na(d[[paste0(v, "__loser")]])
  tibble(
    instrument = unname(INSTRUMENT_DISPLAY[v]),
    variable = v,
    n_elections_scored = sum(ok),
    pct_of_spine = round(100 * mean(ok), 1)
  )
})
write_csv(coverage, file.path(out_dir, "instrument_coverage.csv"))
save_table_html(
  coverage |>
    rename(
      Instrument = instrument,
      Variable = variable,
      `N elections with both top-2 scored` = n_elections_scored,
      `% of spine` = pct_of_spine
    ),
  file.path(out_dir, "instrument_coverage.html"),
  "Instrument coverage on the election spine",
  sprintf(
    "Spine = the %d elections scored under %s",
    nrow(d),
    INSTRUMENT_DISPLAY[[OVERLAP_INSTRUMENT_BUILD]]
  ),
  note = paste(
    "ep_galtan is the CHES bridge (European parties, recent decades only) and is",
    "excluded from the UpSet plots: including it would restrict every intersection",
    "to its own tiny complete-case sample."
  )
)

# ------------------------------------------------------------------------------
# Set membership + narrow-election thresholds
# ------------------------------------------------------------------------------

set_labels <- c(
  v2xpa_antiplural = "Anti-pluralism",
  v2xpa_popul = "Populism",
  v2pariglef_neg = "Economic left",
  v2paanteli = "Anti-elitism",
  ep_galtan = "GAL-TAN"
)

# Which instrument sets to draw an UpSet for. Each entry produces its own plot
# (one per sample definition), so the comparison set is an input rather than a
# hard-coded choice inside the script.
#
# Note the trade-off between entries: an UpSet needs COMPLETE CASES across every
# instrument in the set, since an election missing one score is neither in nor
# out of that set. Adding a sparsely-covered instrument therefore shrinks the
# whole plot's sample, not just its own bar -- ep_galtan (the CHES bridge,
# European parties in recent decades) has both top-2 scores for only 53 of 1,347
# elections, so the "all" set is drawn on those 53. See instrument_coverage.html.
if (!exists("UPSET_INSTRUMENT_SETS")) {
  UPSET_INSTRUMENT_SETS <- list(
    antiplural_vs_popul = c("v2xpa_antiplural", "v2xpa_popul"),
    core3 = c("v2xpa_antiplural", "v2xpa_popul", "v2pariglef_neg"),
    all = ILLIBERALISM_VARS
  )
}
stopifnot(
  all(unlist(UPSET_INSTRUMENT_SETS) %in% names(set_labels)),
  all(lengths(UPSET_INSTRUMENT_SETS) >= 2)
)

# One row per election; one logical column per instrument, TRUE when the actual
# election WINNER scored higher than the loser on that instrument.
#
# Built with mutate(), not select() + bind_cols(): bind_cols aligns by row
# position, so it was only correct while nothing filtered or reordered between
# the two sides -- a filter added to that pipeline would have misaligned the
# memberships silently. mutate() is immune.
build_memberships <- function(instruments) {
  cols <- unname(set_labels[instruments])
  out <- d |>
    select(
      election_id,
      election_type,
      country_text_id,
      election_year,
      running_var
    )
  for (v in instruments) {
    out[[unname(set_labels[v])]] <- d[[paste0(v, "__winner")]] >
      d[[paste0(v, "__loser")]]
  }
  n_before <- nrow(out)
  out <- out |> drop_na(all_of(cols))
  cat(sprintf(
    "  %-22s %d instruments, complete cases %d / %d elections\n",
    paste0("[", paste(instruments, collapse = ","), "]"),
    length(instruments),
    nrow(out),
    n_before
  ))
  out
}

# The bandwidth the headline first stage is actually computed over.
fs_fit <- safe_rdrobust(d$backsliding_Nyr, d$running_var)
mse_bw <- if (!is.null(fs_fit)) unname(fs_fit$bws["h", 1]) else NA_real_
cat(sprintf("MSE-optimal first-stage bandwidth: %.2f pp\n", mse_bw))

SAMPLES <- list(
  all = list(
    label = "All elections",
    fn = function(df) df
  ),
  narrow5 = list(
    label = sprintf("Narrow: |margin| <= %g pp", NARROW_MARGIN_PP),
    fn = function(df) df |> filter(abs(running_var) <= NARROW_MARGIN_PP)
  ),
  narrow_bw = list(
    label = sprintf("Narrow: |margin| <= MSE bandwidth (%.1f pp)", mse_bw),
    fn = function(df) df |> filter(abs(running_var) <= mse_bw)
  )
)

# ------------------------------------------------------------------------------
# UpSet plots
# ------------------------------------------------------------------------------

# An UpSet plot, built directly rather than via ComplexUpset: that package
# requires ggplot2 >= 3.5 (it registers the `axis.text.theta` theme element),
# while this project pins 3.4.4 for every other script. Three aligned panels --
# intersection sizes on top, the membership dot matrix below, set sizes on the
# left -- assembled with patchwork, which aligns the shared axes.
build_upset <- function(df, cols, title, subtitle) {
  combos <- df |>
    count(across(all_of(cols)), name = "n_elections") |>
    arrange(desc(n_elections)) |>
    mutate(combo = row_number())

  # Membership dot matrix: one row per set, one column per intersection.
  matrix_df <- combos |>
    select(combo, all_of(cols)) |>
    pivot_longer(all_of(cols), names_to = "set", values_to = "member") |>
    mutate(set = factor(set, levels = rev(cols)))

  # The connecting line runs between the first and last member set of each
  # intersection -- the visual cue that makes an UpSet readable at a glance.
  spans <- matrix_df |>
    filter(member) |>
    group_by(combo) |>
    summarise(
      lo = min(as.integer(set)),
      hi = max(as.integer(set)),
      .groups = "drop"
    ) |>
    filter(hi > lo)

  set_totals <- df |>
    summarise(across(all_of(cols), sum)) |>
    pivot_longer(everything(), names_to = "set", values_to = "n") |>
    mutate(set = factor(set, levels = rev(cols)))

  x_scale <- scale_x_continuous(
    breaks = combos$combo,
    limits = c(0.4, nrow(combos) + 0.6),
    expand = c(0, 0)
  )

  p_top <- ggplot(combos, aes(x = combo, y = n_elections)) +
    geom_col(fill = SERIES_COLORS[1], width = 0.68) +
    geom_text(
      aes(label = n_elections),
      vjust = -0.4,
      size = 2.5,
      colour = "grey25"
    ) +
    x_scale +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(y = "Elections", x = NULL) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "grey60"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  p_matrix <- ggplot(matrix_df, aes(x = combo, y = set)) +
    geom_point(colour = "grey88", size = 2.8) +
    geom_segment(
      data = spans,
      aes(x = combo, xend = combo, y = lo, yend = hi),
      inherit.aes = FALSE,
      colour = "grey30",
      linewidth = 0.5
    ) +
    geom_point(
      data = filter(matrix_df, member),
      colour = "grey20",
      size = 2.8
    ) +
    x_scale +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank()
    )

  p_sets <- ggplot(set_totals, aes(x = n, y = set)) +
    geom_col(fill = "grey60", width = 0.55) +
    geom_text(aes(label = n), hjust = -0.15, size = 2.4, colour = "grey25") +
    scale_x_reverse(expand = expansion(mult = c(0.28, 0))) +
    labs(x = "Set size", y = NULL) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = "grey60"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  (plot_spacer() + p_top + p_sets + p_matrix) +
    plot_layout(ncol = 2, widths = c(0.3, 1), heights = c(1, 0.75)) +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 7.5)
      )
    )
}

cat("\nBuilding memberships per instrument set:\n")
for (set_nm in names(UPSET_INSTRUMENT_SETS)) {
  instruments <- UPSET_INSTRUMENT_SETS[[set_nm]]
  set_cols <- unname(set_labels[instruments])
  memberships <- build_memberships(instruments)

  for (nm in names(SAMPLES)) {
    samp <- SAMPLES[[nm]]
    df <- samp$fn(memberships)
    if (nrow(df) < 10) {
      message(sprintf(
        "Skipping UpSet (only %d elections): %s / %s",
        nrow(df), set_nm, samp$label
      ))
      next
    }
    # Naming every instrument in the title overflows once there are more than
    # about three, and ggplot truncates rather than wraps. Past that, the title
    # states the count and the subtitle carries the names.
    plot_title <- if (length(set_cols) <= 3) {
      sprintf(
        "Instrument agreement (%s): %s (N = %d)",
        paste(set_cols, collapse = " vs "), samp$label, nrow(df)
      )
    } else {
      sprintf(
        "Instrument agreement across %d instruments: %s (N = %d)",
        length(set_cols), samp$label, nrow(df)
      )
    }

    p <- build_upset(
      df,
      set_cols,
      plot_title,
      paste(
        c(
          if (length(set_cols) > 3) {
            paste("Instruments:", paste(set_cols, collapse = ", "))
          },
          "Each column is one combination of instruments that agree the winner was the more extreme of the top 2.",
          # Only mention an instrument that is actually in this set.
          if ("Economic left" %in% set_cols) {
            "Economic left is the NEGATED v2pariglef, so its high end is economic LEFT."
          },
          if ("GAL-TAN" %in% set_cols) {
            sprintf(
              "GAL-TAN's sparse CHES coverage restricts this set to %d complete cases.",
              nrow(memberships)
            )
          }
        ),
        collapse = "\n"
      )
    )
    # Width scales with the number of intersection columns, which grows with the
    # set size (a 2-instrument set has at most 4, a 5-instrument set up to 32) --
    # but never below what the title needs, or ggplot truncates it mid-word.
    n_combos <- nrow(distinct(df[set_cols]))
    ggsave(
      file.path(out_dir, sprintf("upset_%s_%s.png", set_nm, nm)),
      p,
      width = max(8.5, 1.9 + 0.45 * n_combos),
      height = 2.6 + 0.55 * length(set_cols),
      dpi = 150
    )
    cat(sprintf("Saved upset_%s_%s.png\n", set_nm, nm))

    # The same information as a table, so the exact counts are readable and the
    # figure isn't the only record.
    counts <- df |>
      count(across(all_of(set_cols)), name = "n_elections") |>
      arrange(desc(n_elections)) |>
      mutate(across(all_of(set_cols), \(x) if_else(x, "yes", "no")))
    write_csv(
      counts,
      file.path(out_dir, sprintf("upset_%s_%s_counts.csv", set_nm, nm))
    )
  }
}

# ---- pairwise agreement ------------------------------------------------------
# The single number the UpSet is a decomposition of: for what share of
# elections do two instruments make the same winner-vs-loser call?
#
# Computed over EVERY instrument, and pairwise-complete rather than
# complete-across-all: each cell uses the elections where its own two
# instruments are both scored. Restricting to rows complete across all five --
# which is what reusing an UpSet set's memberships would do -- would compute the
# whole table on ep_galtan's 53 elections, including the cells galtan isn't in.
# The N behind each cell is reported alongside, since it now varies.
all_memberships <- d |> select(election_id)
for (v in ILLIBERALISM_VARS) {
  all_memberships[[unname(set_labels[v])]] <- d[[paste0(v, "__winner")]] >
    d[[paste0(v, "__loser")]]
}
all_cols <- unname(set_labels[ILLIBERALISM_VARS])

pairwise <- expand_grid(a = all_cols, b = all_cols) |>
  rowwise() |>
  mutate(
    n_pair = sum(!is.na(all_memberships[[a]]) & !is.na(all_memberships[[b]])),
    agree_pct = round(
      100 * mean(all_memberships[[a]] == all_memberships[[b]], na.rm = TRUE),
      1
    )
  ) |>
  ungroup()

to_mat <- function(col) {
  matrix(
    pairwise[[col]],
    nrow = length(all_cols),
    dimnames = list(all_cols, all_cols),
    byrow = TRUE
  )
}
agree_mat <- to_mat("agree_pct")
n_mat <- to_mat("n_pair")

write_csv(
  bind_rows(
    as_tibble(agree_mat, rownames = "instrument") |> mutate(stat = "agree_pct"),
    as_tibble(n_mat, rownames = "instrument") |> mutate(stat = "n_elections")
  ),
  file.path(out_dir, "pairwise_agreement.csv")
)

agree_tbl <- as_tibble(
  matrix(
    sprintf("%.1f%%  (N=%s)", agree_mat, format(n_mat, big.mark = ",", trim = TRUE)),
    nrow = length(all_cols),
    dimnames = list(all_cols, all_cols)
  ),
  rownames = "Instrument"
)
save_table_html(
  agree_tbl,
  file.path(out_dir, "pairwise_agreement.html"),
  "Pairwise agreement between instruments (% of elections)",
  "Each cell uses the elections where both of its instruments are scored (pairwise complete), so N varies across cells.",
  note = paste(
    "Share of elections on which two instruments make the same call about whether",
    "the winner scored higher than the loser. 50% is what two unrelated binary",
    "measures would give; 100% means the two instruments are interchangeable here."
  )
)

# ---- pairwise correlation ----------------------------------------------------
# The agreement table above is a binarised comparison: it only asks whether two
# instruments rank the winner above the loser, throwing away how far apart the
# scores are. This is the same question on the raw scores, at the PARTY level
# (one row per top-2 finisher) -- the grain at which the instruments are
# actually measured, and the same grain as the Jaccard heatmaps.
#
# Both Pearson and Spearman are reported. Spearman is the one to compare ACROSS
# pairs: the instruments are on different scales (v2xpa_* are [0,1] indices,
# v2pariglef_neg / v2paanteli / ep_galtan are expert scales), and a rank
# correlation is invariant to that, whereas Pearson is not invariant to the
# non-linear rescalings that differ between them.
#
# Pairwise-complete, like the agreement table and for the same reason: requiring
# rows complete across all five would compute every cell on ep_galtan's handful
# of parties. N is reported per cell since it varies a lot.
corr_vars <- ILLIBERALISM_VARS[ILLIBERALISM_VARS %in% unique(unlist(UPSET_INSTRUMENT_SETS))]
corr_labels <- unname(set_labels[corr_vars])

corr_cells <- expand_grid(a = corr_vars, b = corr_vars) |>
  rowwise() |>
  mutate(
    n_pair = sum(!is.na(parties[[a]]) & !is.na(parties[[b]])),
    pearson = if (n_pair >= 3) {
      cor(parties[[a]], parties[[b]], use = "pairwise.complete.obs")
    } else {
      NA_real_
    },
    spearman = if (n_pair >= 3) {
      cor(
        parties[[a]], parties[[b]],
        use = "pairwise.complete.obs", method = "spearman"
      )
    } else {
      NA_real_
    }
  ) |>
  ungroup()

corr_mat <- function(col) {
  matrix(
    corr_cells[[col]],
    nrow = length(corr_vars),
    dimnames = list(corr_labels, corr_labels),
    byrow = TRUE
  )
}
pearson_mat <- corr_mat("pearson")
spearman_mat <- corr_mat("spearman")
corr_n_mat <- corr_mat("n_pair")

write_csv(
  bind_rows(
    as_tibble(pearson_mat, rownames = "instrument") |> mutate(stat = "pearson"),
    as_tibble(spearman_mat, rownames = "instrument") |> mutate(stat = "spearman"),
    as_tibble(corr_n_mat, rownames = "instrument") |> mutate(stat = "n_parties")
  ),
  file.path(out_dir, "pairwise_correlation.csv")
)

corr_text <- matrix(
  ifelse(
    is.na(as.vector(pearson_mat)),
    "--",
    sprintf(
      "%.3f  (rho %.3f)  N=%s",
      as.vector(pearson_mat),
      as.vector(spearman_mat),
      format(as.vector(corr_n_mat), big.mark = ",", trim = TRUE)
    )
  ),
  nrow = length(corr_vars),
  dimnames = list(corr_labels, corr_labels)
)

corr_tbl <- as_tibble(corr_text, rownames = "Instrument")
gt_corr <- corr_tbl |>
  gt() |>
  tab_header(
    title = "Pairwise correlation between instrument scores, party level",
    subtitle = sprintf(
      "One observation per top-2 finisher (%s parties). Pearson, with Spearman in parentheses. Pairwise complete, so N varies across cells.",
      format(nrow(parties), big.mark = ",")
    )
  ) |>
  cols_align(align = "center", columns = -1) |>
  apply_table_style() |>
  tab_source_note(source_note = paste(
    "Compare SPEARMAN across pairs: the instruments sit on different scales, and",
    "a rank correlation is invariant to that while Pearson is not.",
    "Economic left is the NEGATED v2pariglef, so a positive correlation means",
    "more anti-pluralist / more populist parties are further LEFT economically.",
    "Cell colour = correlation, on a fixed -1 to +1 scale (green positive, red negative)."
  ))

# Fixed +/-1 domain, not data-driven: a correlation matrix has a natural scale,
# and letting the ramp stretch to the largest off-diagonal value would make weak
# correlations look strong.
for (j in seq_along(corr_labels)) {
  fills <- diverging_fill(pearson_mat[, j], max_abs = 1)
  for (i in seq_along(corr_labels)) {
    if (is.na(pearson_mat[i, j])) {
      next
    }
    gt_corr <- gt_corr |>
      tab_style(
        style = cell_fill(color = fills[i]),
        locations = cells_body(columns = corr_labels[j], rows = i)
      )
  }
}
gtsave(gt_corr, file.path(out_dir, "pairwise_correlation.html"))
cat(sprintf("Saved %s\n", file.path(out_dir, "pairwise_correlation.html")))

# ------------------------------------------------------------------------------
# Jaccard heatmaps, one pair of instruments at a time
#
# Unit is a top-2 finisher (a party in a specific election), the grain at which
# the instruments are actually measured.
#
# Bins are EQUAL WIDTH on [0, 1] -- ten fixed bins, [0, 0.1], (0.1, 0.2], ...
# The bin edges therefore mean the same thing in both dimensions and across
# every sample, so cells are comparable between heatmaps and the diagonal is a
# genuine "same score" diagonal. Quantile (decile) bins are available via
# JACCARD_BINS but are a different question: they equalise cell counts, which
# makes the diagonal "same RANK" rather than "same score", and they move with
# whichever sample is being cut.
# ------------------------------------------------------------------------------

if (!exists("JACCARD_PAIRS")) {
  JACCARD_PAIRS <- list(c("v2xpa_antiplural", "v2xpa_popul"))
}
if (!exists("JACCARD_BINS")) JACCARD_BINS <- "equal01" # "equal01" | "deciles"
if (!exists("JACCARD_N_BINS")) JACCARD_N_BINS <- 10
stopifnot(JACCARD_BINS %in% c("equal01", "deciles"))

parties_j_all <- parties |>
  left_join(
    d |> select(election_id, running_var),
    by = "election_id"
  )

# bin_breaks() and jaccard_matrix() are in vparty_helpers.R, shared with
# 18_vparty_jaccard_panels.R, which runs the same computation over raw V-Party
# rather than over the top-2 spine. They take the bin mode and count as
# arguments; this script's JACCARD_BINS / JACCARD_N_BINS are passed in at the
# call sites below.

JACCARD_SAMPLES <- list(
  all = list(label = "All top-2 finishers", fn = function(df) df),
  narrow5 = list(
    label = sprintf(
      "Top-2 finishers in elections with |margin| <= %g pp",
      NARROW_MARGIN_PP
    ),
    fn = function(df) df |> filter(abs(running_var) <= NARROW_MARGIN_PP)
  ),
  narrow_bw = list(
    label = sprintf(
      "Top-2 finishers in elections with |margin| <= %.1f pp",
      mse_bw
    ),
    fn = function(df) df |> filter(abs(running_var) <= mse_bw)
  )
)

for (pair in JACCARD_PAIRS) {
  var_a <- pair[1]
  var_b <- pair[2]
  # File names carry the instruments being compared, matching the UpSet naming
  # convention, so a directory of heatmaps is self-describing.
  pair_slug <- sprintf(
    "%s_vs_%s",
    INSTRUMENT_LABELS[[var_a]], INSTRUMENT_LABELS[[var_b]]
  )

  parties_pair <- parties_j_all |>
    filter(if_all(all_of(pair), \(x) !is.na(x)))

  brk_a <- bin_breaks(parties_pair[[var_a]], var_a, JACCARD_BINS, JACCARD_N_BINS)
  brk_b <- bin_breaks(parties_pair[[var_b]], var_b, JACCARD_BINS, JACCARD_N_BINS)
  n_a <- length(brk_a) - 1L
  n_b <- length(brk_b) - 1L
  # Label each bin by its interval, not by an opaque "D3" -- the whole point of
  # equal-width bins is that the edges are interpretable.
  lab_a <- sprintf("%.1f-%.1f", brk_a[-length(brk_a)], brk_a[-1])
  lab_b <- sprintf("%.1f-%.1f", brk_b[-length(brk_b)], brk_b[-1])

  parties_pair <- parties_pair |>
    mutate(
      bin_a = cut(.data[[var_a]], brk_a, include.lowest = TRUE, labels = FALSE),
      bin_b = cut(.data[[var_b]], brk_b, include.lowest = TRUE, labels = FALSE)
    )

for (nm in names(JACCARD_SAMPLES)) {
  js <- JACCARD_SAMPLES[[nm]]
  df <- js$fn(parties_pair)
  if (nrow(df) < 20) {
    message("Skipping Jaccard (too few observations): ", js$label)
    next
  }
  m <- jaccard_matrix(df, n_a, n_b, lab_a, lab_b)
  write_csv(
    as_tibble(m, rownames = sprintf("%s_bin", var_a)),
    file.path(out_dir, sprintf("jaccard_%s_%s.csv", pair_slug, nm))
  )

  long <- as_tibble(m, rownames = "row") |>
    pivot_longer(-row, names_to = "col", values_to = "jaccard") |>
    mutate(
      row_i = match(row, lab_a),
      col_i = match(col, lab_b)
    )

  # Sequential = ONE hue, light to dark. The quantity here is a magnitude
  # (similarity from 0 to 1) with no meaningful midpoint, so a diverging or
  # multi-hue ramp would invent a polarity the data doesn't have.
  p <- ggplot(long, aes(x = col_i, y = row_i, fill = jaccard)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(
      aes(label = if_else(is.na(jaccard), "", sprintf("%.2f", jaccard))),
      size = 2.5,
      colour = if_else(long$jaccard > 0.28, "white", "grey20")
    ) +
    scale_fill_gradient(
      low = "#eaf2f8",
      high = "#0072B2",
      na.value = "grey92",
      limits = c(0, NA),
      name = "Jaccard"
    ) +
    scale_x_continuous(
      breaks = seq_len(n_b), labels = lab_b, expand = c(0, 0)
    ) +
    scale_y_reverse(
      breaks = seq_len(n_a), labels = lab_a, expand = c(0, 0)
    ) +
    coord_fixed() +
    labs(
      title = sprintf(
        "%s vs %s: %s",
        INSTRUMENT_DISPLAY[[var_a]], INSTRUMENT_DISPLAY[[var_b]], js$label
      ),
      subtitle = paste(
        strwrap(
          sprintf(
            "Jaccard = |both| / |either|, over %d top-2 finishers. Bins are %s. A perfectly redundant pair would light only the diagonal.",
            nrow(df),
            if (JACCARD_BINS == "equal01") {
              sprintf("%d equal-width intervals on [0, 1]", JACCARD_N_BINS)
            } else {
              sprintf("%d quantiles of the pooled top-2 distribution", JACCARD_N_BINS)
            }
          ),
          width = 95
        ),
        collapse = "\n"
      ),
      x = sprintf("%s bin", INSTRUMENT_DISPLAY[[var_b]]),
      y = sprintf("%s bin", INSTRUMENT_DISPLAY[[var_a]])
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 7.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7)
    )

  ggsave(
    file.path(out_dir, sprintf("jaccard_%s_%s.png", pair_slug, nm)),
    p,
    width = 7.5,
    height = 6.5,
    dpi = 150
  )
  cat(sprintf("Saved jaccard_%s_%s.png\n", pair_slug, nm))
}
}

message("Instrument overlap written to ", out_dir)
