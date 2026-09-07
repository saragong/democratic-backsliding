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

data_dir <- here::here("data")

if (!exists("OVERLAP_INSTRUMENT_BUILD")) OVERLAP_INSTRUMENT_BUILD <- "v2xpa_antiplural"
if (!exists("OVERLAP_WINDOW")) OVERLAP_WINDOW <- 5
if (!exists("NARROW_MARGIN_PP")) NARROW_MARGIN_PP <- 5

out_dir <- sweep_dir("instrument_overlap")

build_path <- file.path(
  data_dir, "rdd_build",
  sprintf("rdd_%s_w%d.rds", OVERLAP_INSTRUMENT_BUILD, OVERLAP_WINDOW)
)
parties_path <- sub("\\.rds$", "_parties.rds", build_path)
stopifnot(file.exists(build_path), file.exists(parties_path))

d <- readRDS(build_path)
parties <- readRDS(parties_path)

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
      Instrument = instrument, Variable = variable,
      `N elections with both top-2 scored` = n_elections_scored,
      `% of spine` = pct_of_spine
    ),
  file.path(out_dir, "instrument_coverage.html"),
  "Instrument coverage on the election spine",
  sprintf(
    "Spine = the %d elections scored under %s",
    nrow(d), INSTRUMENT_DISPLAY[[OVERLAP_INSTRUMENT_BUILD]]
  ),
  note = paste(
    "ep_galtan is the CHES bridge (European parties, recent decades only) and is",
    "excluded from the UpSet plots: including it would restrict every intersection",
    "to its own tiny complete-case sample."
  )
)

UPSET_INSTRUMENTS <- coverage |>
  filter(n_elections_scored >= 0.5 * nrow(d)) |>
  pull(variable)
cat(sprintf(
  "UpSet instruments (>=50%% coverage): %s\n",
  paste(UPSET_INSTRUMENTS, collapse = ", ")
))

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

memberships <- d |>
  select(election_id, election_type, country_text_id, election_year, running_var) |>
  bind_cols(
    map_dfc(UPSET_INSTRUMENTS, function(v) {
      tibble(
        !!unname(set_labels[v]) := d[[paste0(v, "__winner")]] >
          d[[paste0(v, "__loser")]]
      )
    })
  )

set_cols <- unname(set_labels[UPSET_INSTRUMENTS])

# Complete cases only. An election missing one instrument's score can't be
# placed in or out of that set, so it has no well-defined intersection
# membership and would silently distort every bar.
n_before <- nrow(memberships)
memberships <- memberships |> drop_na(all_of(set_cols))
cat(sprintf(
  "Complete cases across %d instruments: %d / %d elections\n",
  length(set_cols), nrow(memberships), n_before
))

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
      lo = min(as.integer(set)), hi = max(as.integer(set)),
      .groups = "drop"
    ) |>
    filter(hi > lo)

  set_totals <- df |>
    summarise(across(all_of(cols), sum)) |>
    pivot_longer(everything(), names_to = "set", values_to = "n") |>
    mutate(set = factor(set, levels = rev(cols)))

  x_scale <- scale_x_continuous(
    breaks = combos$combo, limits = c(0.4, nrow(combos) + 0.6), expand = c(0, 0)
  )

  p_top <- ggplot(combos, aes(x = combo, y = n_elections)) +
    geom_col(fill = SERIES_COLORS[1], width = 0.68) +
    geom_text(aes(label = n_elections), vjust = -0.4, size = 2.5, colour = "grey25") +
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
      inherit.aes = FALSE, colour = "grey30", linewidth = 0.5
    ) +
    geom_point(
      data = filter(matrix_df, member),
      colour = "grey20", size = 2.8
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

for (nm in names(SAMPLES)) {
  samp <- SAMPLES[[nm]]
  df <- samp$fn(memberships)
  if (nrow(df) < 10) {
    message("Skipping UpSet (too few elections): ", samp$label)
    next
  }
  p <- build_upset(
    df, set_cols,
    sprintf("Instrument agreement: %s (N = %d)", samp$label, nrow(df)),
    paste0(
      "Each column is one combination of instruments that agree the winner was the more extreme of the top 2.\n",
      "Economic left is the NEGATED v2pariglef, so its high end is economic LEFT."
    )
  )
  ggsave(
    file.path(out_dir, sprintf("upset_%s.png", nm)),
    p, width = 9, height = 5.2, dpi = 150
  )
  cat(sprintf("Saved upset_%s.png\n", nm))

  # The same information as a table, so the exact counts are readable and the
  # figure isn't the only record.
  counts <- df |>
    count(across(all_of(set_cols)), name = "n_elections") |>
    arrange(desc(n_elections)) |>
    mutate(across(all_of(set_cols), \(x) if_else(x, "yes", "no")))
  write_csv(counts, file.path(out_dir, sprintf("upset_%s_counts.csv", nm)))
}

# ---- pairwise agreement ------------------------------------------------------
# The single number the UpSet is a decomposition of: for what share of
# elections do two instruments make the same winner-vs-loser call?
pairwise <- expand_grid(a = set_cols, b = set_cols) |>
  rowwise() |>
  mutate(
    agree_pct = round(100 * mean(memberships[[a]] == memberships[[b]]), 1)
  ) |>
  ungroup()

agree_mat <- matrix(
  pairwise$agree_pct,
  nrow = length(set_cols),
  dimnames = list(set_cols, set_cols),
  byrow = TRUE
)
write_csv(
  as_tibble(agree_mat, rownames = "instrument"),
  file.path(out_dir, "pairwise_agreement.csv")
)
save_table_html(
  as_tibble(agree_mat, rownames = "Instrument"),
  file.path(out_dir, "pairwise_agreement.html"),
  "Pairwise agreement between instruments (% of elections)",
  sprintf("N = %d elections, complete cases across all instruments", nrow(memberships)),
  note = paste(
    "Share of elections on which two instruments make the same call about whether",
    "the winner scored higher than the loser. 50% is what two unrelated binary",
    "measures would give; 100% means the two instruments are interchangeable here."
  )
)

# ------------------------------------------------------------------------------
# Jaccard heatmaps: anti-pluralism deciles x populism deciles
#
# Unit is a top-2 finisher (a party in a specific election), the grain at which
# the instruments are actually measured. Deciles are cut on the pooled
# distribution of ALL top-2 finishers, so the same bin boundaries apply to both
# samples and the two heatmaps are directly comparable.
# ------------------------------------------------------------------------------

JACCARD_PAIR <- c("v2xpa_antiplural", "v2xpa_popul")

parties_j <- parties |>
  left_join(
    d |> select(election_id, running_var),
    by = "election_id"
  ) |>
  filter(if_all(all_of(JACCARD_PAIR), \(x) !is.na(x)))

decile_breaks <- function(x) {
  unique(quantile(x, probs = seq(0, 1, 0.1), na.rm = TRUE))
}
brk_a <- decile_breaks(parties_j[[JACCARD_PAIR[1]]])
brk_b <- decile_breaks(parties_j[[JACCARD_PAIR[2]]])

parties_j <- parties_j |>
  mutate(
    dec_a = cut(.data[[JACCARD_PAIR[1]]], brk_a, include.lowest = TRUE, labels = FALSE),
    dec_b = cut(.data[[JACCARD_PAIR[2]]], brk_b, include.lowest = TRUE, labels = FALSE)
  )

jaccard_matrix <- function(df) {
  n_a <- max(df$dec_a, na.rm = TRUE)
  n_b <- max(df$dec_b, na.rm = TRUE)
  m <- matrix(NA_real_, nrow = n_a, ncol = n_b)
  for (i in seq_len(n_a)) {
    for (j in seq_len(n_b)) {
      in_a <- df$dec_a == i
      in_b <- df$dec_b == j
      union_n <- sum(in_a | in_b, na.rm = TRUE)
      m[i, j] <- if (union_n == 0) NA_real_ else sum(in_a & in_b, na.rm = TRUE) / union_n
    }
  }
  dimnames(m) <- list(paste0("D", seq_len(n_a)), paste0("D", seq_len(n_b)))
  m
}

JACCARD_SAMPLES <- list(
  all = list(label = "All top-2 finishers", fn = function(df) df),
  narrow5 = list(
    label = sprintf("Top-2 finishers in elections with |margin| <= %g pp", NARROW_MARGIN_PP),
    fn = function(df) df |> filter(abs(running_var) <= NARROW_MARGIN_PP)
  ),
  narrow_bw = list(
    label = sprintf("Top-2 finishers in elections with |margin| <= %.1f pp", mse_bw),
    fn = function(df) df |> filter(abs(running_var) <= mse_bw)
  )
)

for (nm in names(JACCARD_SAMPLES)) {
  js <- JACCARD_SAMPLES[[nm]]
  df <- js$fn(parties_j)
  if (nrow(df) < 20) {
    message("Skipping Jaccard (too few observations): ", js$label)
    next
  }
  m <- jaccard_matrix(df)
  write_csv(
    as_tibble(m, rownames = "antipluralism_decile"),
    file.path(out_dir, sprintf("jaccard_%s.csv", nm))
  )

  long <- as_tibble(m, rownames = "row") |>
    pivot_longer(-row, names_to = "col", values_to = "jaccard") |>
    mutate(
      row_i = as.integer(sub("D", "", row)),
      col_i = as.integer(sub("D", "", col))
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
      low = "#eaf2f8", high = "#0072B2", na.value = "grey92",
      limits = c(0, NA), name = "Jaccard"
    ) +
    scale_x_continuous(breaks = seq_len(ncol(m)), expand = c(0, 0)) +
    scale_y_reverse(breaks = seq_len(nrow(m)), expand = c(0, 0)) +
    coord_fixed() +
    labs(
      title = sprintf("Anti-pluralism vs populism deciles: %s", js$label),
      subtitle = sprintf(
        "Jaccard = |both| / |either|, over %d top-2 finishers. A perfectly redundant pair would light only the diagonal.",
        nrow(df)
      ),
      x = "Populism decile (v2xpa_popul)",
      y = "Anti-pluralism decile (v2xpa_antiplural)"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 7.5)
    )

  ggsave(
    file.path(out_dir, sprintf("jaccard_%s.png", nm)),
    p, width = 7, height = 6, dpi = 150
  )
  cat(sprintf("Saved jaccard_%s.png\n", nm))
}

message("Instrument overlap written to ", out_dir)
