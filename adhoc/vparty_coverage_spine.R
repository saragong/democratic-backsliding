# ==============================================================================
# V-Party score coverage along the election spine.
#
# Every RDD in this repo runs on the elections where BOTH of the top 2 parties
# carry a V-Party score -- 1,347 of the 3,096 elections that have a valid
# distinct-party top 2. The other 1,749 are dropped at Step 3 of
# 11_build_rdd_data.R and are invisible everywhere downstream. This asks what
# is being dropped and when, by plotting, per 5-year bin:
#
#   both     both top-2 parties have a v2xpa_antiplural score  -> in the RDD
#   one      exactly one does                                   -> dropped
#   (the headroom to 100% is "neither", also dropped)
#
# "one" is the interesting category: those elections have a real top 2 that
# matched parties_database, and the pair is broken by V-Party coverage alone,
# not by the election being unusable.
#
# OECD on top, non-OECD below, because coverage differs sharply between them
# and the pooled line is a mix of two different regimes.
#
#   Rscript --no-init-file adhoc/vparty_coverage_spine.R
#
# Output: output/runs/_sweeps/vparty_coverage_spine/
#           coverage_spine_<score>.png
#           coverage_spine_<score>.csv   one row per bin x group
# ==============================================================================

library(tidyverse)
library(here)

source(here::here("scripts", "rdd_helpers.R"))    # SERIES_COLORS, sweep_dir()
source(here::here("scripts", "vparty_helpers.R")) # oecd_group()

# ---- toggles -----------------------------------------------------------------

# The score whose coverage is being measured. Any of PARTY_SCORE_VARS: the
# spine carries all of them, so this needs no rebuild to change.
if (!exists("COVERAGE_SCORE")) COVERAGE_SCORE <- "v2xpa_antiplural"

# Which spine file, i.e. which (ELECTION_TYPE, SAMPLE_YEARS) the build ran under.
if (!exists("COVERAGE_SPINE")) COVERAGE_SPINE <- "both"

if (!exists("COVERAGE_BIN_WIDTH")) COVERAGE_BIN_WIDTH <- 5

# 1945 rather than the spine's true start of 1805. Before 1945 the bins hold 1
# to 22 elections each and a share computed on them is noise; from 1945 every
# bin holds at least 105. The five bins from 1945 to 1969 are kept even though
# they are structurally empty of scores, because seeing coverage switch on when
# V-Party's own coverage begins is the point of showing them.
if (!exists("COVERAGE_YEAR_MIN")) COVERAGE_YEAR_MIN <- 1945

stopifnot(COVERAGE_SCORE %in% PARTY_SCORE_VARS)

spine_path <- file.path(
  here::here("data", "rdd_build"),
  sprintf("top2_spine_%s.rds", COVERAGE_SPINE)
)
if (!file.exists(spine_path)) {
  stop(
    "No spine at ", spine_path, ". It is written by Step 2b of ",
    "11_build_rdd_data.R -- run that once (any window) to produce it.",
    call. = FALSE
  )
}

out_dir <- sweep_dir("vparty_coverage_spine")

# ---- one row per election ----------------------------------------------------

spine <- readRDS(spine_path)

# Two rows per election by construction (rank 1 and 2). Assert it rather than
# assume: a spine with a stray third row would silently turn "both scored" into
# a count of 2-out-of-3 and every share below would be wrong.
stopifnot(all(table(spine$election_id) == 2))

elections <- spine |>
  group_by(election_id) |>
  summarise(
    election_year = first(election_year),
    country_text_id = first(country_text_id),
    n_scored = sum(!is.na(.data[[COVERAGE_SCORE]])),
    .groups = "drop"
  ) |>
  filter(election_year >= COVERAGE_YEAR_MIN) |>
  mutate(
    group = factor(oecd_group(country_text_id), levels = c("OECD", "Non-OECD")),
    bin_start = (election_year %/% COVERAGE_BIN_WIDTH) * COVERAGE_BIN_WIDTH,
    coverage = factor(
      c("neither", "one", "both")[n_scored + 1],
      levels = c("both", "one", "neither")
    )
  )

year_max <- max(elections$election_year)

# Only claim a partial bin when there actually is one. Selecting with a logical
# that matches nothing yields NA, not NULL, so `%||%` would not catch it and
# the subtitle would read "The NA bin is partial".
partial_note <- local({
  starts <- sort(unique((elections$election_year %/% COVERAGE_BIN_WIDTH) * COVERAGE_BIN_WIDTH))
  cut <- starts[starts + COVERAGE_BIN_WIDTH - 1L > year_max]
  if (length(cut) == 0) {
    ""
  } else {
    sprintf(
      " The %d-%02d bin is partial (data end %d).",
      cut[1], (cut[1] + COVERAGE_BIN_WIDTH - 1L) %% 100, year_max
    )
  }
})

# ---- bin x group -------------------------------------------------------------

# complete() so a bin with no elections of one group is an explicit zero-N row
# rather than a silently missing bar, which would read as 0% coverage.
cells <- elections |>
  count(group, bin_start, coverage, name = "n") |>
  complete(
    group, bin_start = seq(
      min(elections$bin_start), max(elections$bin_start), by = COVERAGE_BIN_WIDTH
    ),
    coverage, fill = list(n = 0L)
  ) |>
  group_by(group, bin_start) |>
  mutate(n_elections = sum(n), share = if_else(n_elections > 0, n / n_elections, NA_real_)) |>
  ungroup() |>
  mutate(
    bin_end = bin_start + COVERAGE_BIN_WIDTH - 1L,
    # bin_label is a STRING ("1970-74"); bin_start next to it is the integer to
    # sort or filter on. read_csv() parses a label like "1970s" as the number
    # 1970 without complaint, which has bitten this repo twice, so the numeric
    # column always travels with the label.
    bin_label = sprintf("%d-%02d", bin_start, bin_end %% 100),
    partial = bin_end > year_max
  )

wide <- cells |>
  select(group, bin_start, bin_label, bin_end, partial, n_elections, coverage, n) |>
  pivot_wider(names_from = coverage, values_from = n, names_prefix = "n_") |>
  mutate(
    share_both = if_else(n_elections > 0, n_both / n_elections, NA_real_),
    share_one = if_else(n_elections > 0, n_one / n_elections, NA_real_),
    share_neither = if_else(n_elections > 0, n_neither / n_elections, NA_real_)
  ) |>
  arrange(group, bin_start)

write_csv(wide, file.path(out_dir, sprintf("coverage_spine_%s.csv", COVERAGE_SCORE)))

cat("\nCoverage by group (all bins pooled):\n")
print(
  elections |>
    count(group, coverage) |>
    group_by(group) |>
    mutate(share = sprintf("%.1f%%", 100 * n / sum(n))) |>
    ungroup() |>
    pivot_wider(names_from = coverage, values_from = c(n, share)) |>
    as.data.frame(),
  row.names = FALSE
)

# ---- figure ------------------------------------------------------------------

# Stacked bars, not lines: these are shares of a whole over discrete bins, so
# the bar top reads directly as "at least one party scored" and the white
# headroom to 100% is the elections V-Party cannot see at all. Only "both" and
# "one" are drawn -- "neither" is the headroom, named in the subtitle rather
# than given a third grey band that would hide the 100% reference.
plot_dat <- cells |>
  filter(coverage != "neither", n_elections > 0) |>
  mutate(coverage = factor(
    coverage,
    levels = c("one", "both"), # "both" drawn at the bottom of the stack
    labels = c(
      "Only one of the top 2 scored (dropped)",
      "Both top-2 parties scored (in the RDD)"
    )
  ))

# The count row sits below the axis line, where it could be mistaken for a
# second set of tick labels, so the leftmost bin of each panel carries an
# "n =" to say what the row is.
n_lab <- cells |>
  distinct(group, bin_start, bin_label, n_elections, partial) |>
  filter(n_elections > 0) |>
  group_by(group) |>
  mutate(label = if_else(
    bin_start == min(bin_start),
    paste0("n = ", n_elections), as.character(n_elections)
  )) |>
  ungroup()

p <- ggplot(plot_dat, aes(x = bin_start, y = share, fill = coverage)) +
  geom_col(width = COVERAGE_BIN_WIDTH * 0.85) +
  geom_text(
    data = n_lab, inherit.aes = FALSE,
    aes(x = bin_start, y = -0.055, label = label),
    size = 2.1, colour = "grey35"
  ) +
  facet_wrap(~group, ncol = 1) +
  scale_fill_manual(values = setNames(
    SERIES_COLORS[c(3, 1)],
    levels(plot_dat$coverage)
  ), name = NULL, guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, 0.25), limits = c(-0.09, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_continuous(
    breaks = sort(unique(cells$bin_start)),
    labels = function(b) sub("^\\d{2}", "", as.character(b)) |>
      (\(x) paste0("'", x))()
  ) +
  labs(
    title = sprintf(
      "V-Party coverage of the top 2: %s, by %d-year bin",
      PARTY_SCORE_DISPLAY[[COVERAGE_SCORE]], COVERAGE_BIN_WIDTH
    ),
    subtitle = paste(strwrap(sprintf(paste(
      "Denominator is every election in the spine with a valid distinct-party",
      "top 2 (%d elections from %d, both chambers). Bar height is the share",
      "with at least one top-2 party scored; the headroom to 100%% is the",
      "elections where NEITHER is scored. Only the dark band enters the RDD --",
      "the light band is elections broken by V-Party coverage alone, with a",
      "real top 2 that matched parties_database. V-Party's own coverage starts",
      "in 1970, which is why the earlier bins are empty; scores attach past",
      "2019 because a party's most recent score at or before the election year",
      "is carried forward. Small numbers under each bar are the elections in",
      "that bin. OECD is current membership applied to all years.%s"
    ), nrow(elections), COVERAGE_YEAR_MIN, partial_note), 128),
    collapse = "\n"),
    x = sprintf("Election year (%d-year bins, labelled by bin start)", COVERAGE_BIN_WIDTH),
    y = "Share of elections in the bin"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "bottom",
    legend.margin = margin(t = -2),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 7, colour = "grey25")
  )

ggsave(
  file.path(out_dir, sprintf("coverage_spine_%s.png", COVERAGE_SCORE)), p,
  width = 9.5, height = 6.4, dpi = 150
)

message("\nCoverage figure written to ", out_dir)
