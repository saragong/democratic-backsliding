# ==============================================================================
# Who are the illiberal parties? Named examples, on both sides of the cutoff.
#
# The RDD compares elections a highly anti-pluralist party narrowly WON against
# elections it narrowly LOST. That comparison is only persuasive if the score
# picks out parties a reader recognizes as illiberal, so this lists them by
# name.
#
# Three tables, because the top of the score distribution and the parties that
# actually identify the effect are not the same set:
#
#   examples_narrow        the identifying variation -- highest-scoring parties
#                          within +/- 5 pp of a tie, split into narrowly won
#                          and narrowly lost. These are the cases the estimate
#                          is built from.
#   examples_top_scorers   the top of the distribution regardless of margin.
#                          Almost all are hegemonic-party regimes winning by
#                          50-98 pp, which is exactly why they contribute
#                          nothing at the cutoff -- worth showing so the
#                          measure's top end is not mistaken for the sample.
#   examples_repeat        parties appearing in more than one narrow election,
#                          including ones that appear on BOTH sides.
#
#   Rscript --no-init-file adhoc/illiberal_party_examples.R
#
# Output: output/runs/_sweeps/illiberal_party_examples/
# ==============================================================================

library(tidyverse)
library(here)
library(gt)

source(here::here("scripts", "rdd_helpers.R"))
source(here::here("scripts", "vparty_helpers.R"))

# ---- toggles -----------------------------------------------------------------

if (!exists("EXAMPLES_INSTRUMENT")) EXAMPLES_INSTRUMENT <- "v2xpa_antiplural"
if (!exists("EXAMPLES_WINDOW")) EXAMPLES_WINDOW <- 5
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- FALSE
}

# "Narrow" for illustration only. The RDD's own bandwidth is about 18 pp, so
# this is far tighter than the estimation window -- it is picking the most
# vivid cases, not reproducing the sample. Stated on the table so the two are
# not confused.
if (!exists("EXAMPLES_NARROW_MARGIN")) EXAMPLES_NARROW_MARGIN <- 5

if (!exists("EXAMPLES_N")) EXAMPLES_N <- 15

build_suffix <- if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "" else "_exclyr"
stem <- sprintf(
  "rdd_%s_w%d%s", EXAMPLES_INSTRUMENT, EXAMPLES_WINDOW, build_suffix
)
build_path <- here::here("data", "rdd_build", paste0(stem, ".rds"))
parties_path <- here::here("data", "rdd_build", paste0(stem, "_parties.rds"))
if (!file.exists(build_path) || !file.exists(parties_path)) {
  stop("No build at ", build_path, ". Run 11_build_rdd_data.R first.", call. = FALSE)
}

out_dir <- sweep_dir("illiberal_party_examples")

# ---- party names -------------------------------------------------------------

d <- readRDS(build_path)
parties <- readRDS(parties_path)

# Name the party AS OF THE ELECTION, by the same rule 11_build_rdd_data.R's
# match_vparty() uses to pick the score: the most recent V-Party row at or
# before the election year. Taking the party's latest name instead would label
# a 1994 election with a name the party only adopted in 2015, which for
# renamed and merged parties is exactly the kind of quiet mismatch that makes
# a reader distrust the whole table.
vparty_names <- load_vparty_raw(
  c("v2paid", "year", "v2paenname", "v2pashname")
)

name_at <- parties |>
  distinct(vdem_id_1, election_year) |>
  filter(!is.na(vdem_id_1)) |>
  left_join(vparty_names, by = c("vdem_id_1" = "v2paid"),
            relationship = "many-to-many") |>
  filter(year <= election_year) |>
  group_by(vdem_id_1, election_year) |>
  slice_max(year, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(vdem_id_1, election_year, party_name = v2paenname, party_abbr = v2pashname)

# One row per election: the MORE anti-pluralist of the top 2 and its opponent.
# Ordered by the score, not by vote share, so "ap_" is always the illiberal
# side regardless of who won -- the same ordering the running variable uses.
pairs <- parties |>
  left_join(name_at, by = c("vdem_id_1", "election_year")) |>
  group_by(election_id) |>
  arrange(desc(.data[[EXAMPLES_INSTRUMENT]]), .by_group = TRUE) |>
  summarise(
    ap_name = first(party_name), ap_abbr = first(party_abbr),
    ap_score = first(.data[[EXAMPLES_INSTRUMENT]]),
    ap_share = first(final_share), ap_cand = first(candidate),
    op_name = last(party_name), op_abbr = last(party_abbr),
    op_score = last(.data[[EXAMPLES_INSTRUMENT]]),
    op_share = last(final_share), op_cand = last(candidate),
    .groups = "drop"
  )

popucut <- resolve_threshold_abs("popucut", "popucut", NULL)

ex <- d |>
  select(
    election_id, country_name, country_text_id, election_year, election_type,
    running_var, illiberal_score, other_score,
    backsliding_Nyr, Y_gdp_growth, Y_polyarchy
  ) |>
  left_join(pairs, by = "election_id") |>
  mutate(
    group = oecd_group(country_text_id),
    side = if_else(running_var > 0, "won", "lost"),
    # The headline restricted sample, flagged rather than filtered so the
    # examples stay drawn from the full build.
    in_popucut = illiberal_score > popucut & other_score <= popucut,
    label = if_else(
      is.na(ap_name), coalesce(ap_cand, "(unnamed)"),
      if_else(is.na(ap_abbr) | ap_abbr == "", ap_name,
              sprintf("%s (%s)", ap_name, ap_abbr))
    ),
    op_label = coalesce(op_name, op_cand, "(unnamed)")
  )

stopifnot(nrow(ex) == nrow(d))
write_csv(
  ex |> arrange(desc(ap_score)),
  file.path(out_dir, "illiberal_party_examples.csv")
)

# ---- shared formatting -------------------------------------------------------

TYPE_SHORT <- c(presidential = "Pres.", parliamentary = "Parl.")

fmt_rows <- function(x) {
  x |>
    transmute(
      Country = country_name,
      Year = election_year,
      Type = unname(TYPE_SHORT[election_type]),
      `More anti-pluralist party` = label,
      Score = round(ap_score, 3),
      Opponent = op_label,
      `Opp. score` = round(op_score, 3),
      `Margin (pp)` = round(running_var, 1),
      `ERT episode` = if_else(backsliding_Nyr == 1, "yes", ""),
      `Log GDP pc chg` = if_else(is.na(Y_gdp_growth), NA_character_,
                                 sprintf("%+.2f", Y_gdp_growth)),
      `Log polyarchy chg` = if_else(is.na(Y_polyarchy), NA_character_,
                                    sprintf("%+.3f", Y_polyarchy)),
      `PopuList 411` = if_else(in_popucut, "yes", "")
    )
}

OUTCOME_NOTE <- paste(
  "'ERT episode' is an autocratization episode starting within", EXAMPLES_WINDOW,
  "years of the election. Both outcome columns are LOG changes over the same",
  sprintf(
    "window, ln(value at election year + %d) - ln(value at election year - 1);",
    EXAMPLES_WINDOW
  ),
  "positive polyarchy means the index ROSE, i.e. no backsliding on that",
  "measure. These are single elections, shown to make the cases concrete --",
  "they are illustrations, not evidence, and any one of them can run the",
  "other way. 'PopuList 411' marks the elections in the restricted sample",
  sprintf(
    "where one top-2 party scores above %.4f and the other does not.", popucut
  ),
  "Where the opponent's score is also near the ceiling (Morocco 1993:",
  "0.985 vs 0.853) the sign of the running variable rests on a small",
  "difference between two parties both coded as anti-pluralist; SCORE_GAP_MIN",
  "in 12_rdd_analysis.R is the lever for excluding those."
)

save_grouped <- function(tbl, path, title, subtitle, note) {
  gt_tbl <- tbl |>
    gt(groupname_col = "panel") |>
    tab_header(title = title, subtitle = subtitle) |>
    sub_missing(missing_text = "--") |>
    # any_of(), not c(): the three tables below share this function but not
    # their column sets, and gt errors on a name that is not there.
    cols_align("left", columns = any_of(c(
      "Country", "More anti-pluralist party", "Opponent", "Party",
      "Elections (margin, pp)"
    ))) |>
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_row_groups()
    ) |>
    opt_row_striping() |>
    apply_table_style(font_size = 11) |>
    tab_source_note(note)
  gtsave(gt_tbl, path)
  cat(sprintf("Saved %s\n", path))
}

# ---- 1. the identifying variation -------------------------------------------

narrow <- ex |> filter(abs(running_var) <= EXAMPLES_NARROW_MARGIN)

narrow_tbl <- bind_rows(
  narrow |> filter(side == "won") |> arrange(desc(ap_score)) |> head(EXAMPLES_N) |>
    fmt_rows() |> mutate(panel = sprintf(
      "Anti-pluralist party narrowly WON (0 < margin <= %g pp)", EXAMPLES_NARROW_MARGIN
    )),
  narrow |> filter(side == "lost") |> arrange(desc(ap_score)) |> head(EXAMPLES_N) |>
    fmt_rows() |> mutate(panel = sprintf(
      "Anti-pluralist party narrowly LOST (-%g <= margin < 0 pp)", EXAMPLES_NARROW_MARGIN
    ))
)

save_grouped(
  narrow_tbl,
  file.path(out_dir, "examples_narrow.html"),
  sprintf(
    "Highest-scoring anti-pluralist parties in near-tied elections (top %d each side)",
    EXAMPLES_N
  ),
  sprintf(paste(
    "%s. Parties ordered by the score, so the left-hand party is the more",
    "anti-pluralist of the top 2 whether or not it won. Margin is its vote or",
    "seat share minus its opponent's. Note the +/- %g pp window is tighter",
    "than the RDD's own bandwidth of roughly 18 pp: these are the most vivid",
    "cases, not the estimation sample. Party names are as of the election year."
  ), INSTRUMENT_DISPLAY[[EXAMPLES_INSTRUMENT]], EXAMPLES_NARROW_MARGIN),
  OUTCOME_NOTE
)

# ---- 2. the top of the distribution -----------------------------------------

top_tbl <- ex |>
  arrange(desc(ap_score)) |>
  head(EXAMPLES_N * 2) |>
  fmt_rows() |>
  mutate(panel = "Highest scores in the sample, any margin")

save_grouped(
  top_tbl,
  file.path(out_dir, "examples_top_scorers.html"),
  sprintf("Highest %s scores overall, any margin", INSTRUMENT_DISPLAY[[EXAMPLES_INSTRUMENT]]),
  paste(
    "The ceiling of the score is occupied by hegemonic-party regimes winning",
    "by enormous margins, which is why they contribute nothing to an estimate",
    "taken at a near-tie. Shown so the top of the measure is not mistaken for",
    "the variation the design uses -- for that, see examples_narrow.html."
  ),
  OUTCOME_NOTE
)

# ---- 3. parties seen more than once in a near-tie ---------------------------

# A party that narrowly won once and narrowly lost another time is the design
# in miniature, so those are listed first.
repeats <- narrow |>
  filter(!is.na(ap_name)) |>
  group_by(country_name, ap_name, ap_abbr) |>
  filter(n() > 1) |>
  summarise(
    n_narrow = n(),
    n_won = sum(side == "won"), n_lost = sum(side == "lost"),
    both_sides = n_won > 0 & n_lost > 0,
    mean_score = mean(ap_score),
    years = paste(sprintf(
      "%d (%s%.1f)", election_year, if_else(running_var > 0, "+", ""), running_var
    ), collapse = ", "),
    .groups = "drop"
  ) |>
  arrange(desc(both_sides), desc(mean_score))

save_grouped(
  repeats |>
    transmute(
      Country = country_name,
      Party = if_else(is.na(ap_abbr) | ap_abbr == "", ap_name,
                      sprintf("%s (%s)", ap_name, ap_abbr)),
      `Mean score` = round(mean_score, 3),
      `Near-ties` = n_narrow, Won = n_won, Lost = n_lost,
      `Elections (margin, pp)` = years,
      panel = if_else(both_sides,
                      "Appears on BOTH sides of the cutoff",
                      "Appears more than once on one side")
    ),
  file.path(out_dir, "examples_repeat.html"),
  sprintf("Parties in more than one near-tied election (within %g pp)", EXAMPLES_NARROW_MARGIN),
  paste(
    "A party that narrowly won one election and narrowly lost another is the",
    "research design in miniature: the same party, the same country, a",
    "coin-flip apart. Those are listed first."
  ),
  sprintf(
    "Ordered by mean %s score within each block.",
    INSTRUMENT_DISPLAY[[EXAMPLES_INSTRUMENT]]
  )
)

cat(sprintf(
  "\n%d elections within %g pp: %d narrow wins, %d narrow losses; %d parties appear on both sides.\n",
  nrow(narrow), EXAMPLES_NARROW_MARGIN,
  sum(narrow$side == "won"), sum(narrow$side == "lost"), sum(repeats$both_sides)
))
message("\nExamples written to ", out_dir)
