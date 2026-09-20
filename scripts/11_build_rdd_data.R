# ==============================================================================
# Fuzzy RDD analysis dataset: illiberal-party elections -> democratic
# backsliding -> economic outcomes
#
# Running variable: vote/seat-share margin of the more-illiberal of the
# top-2 parties/candidates in an election (positive = illiberal side leads).
# Fuzzy treatment: whether a democratic-backsliding episode (ERT aut_ep)
# begins within BACKSLIDING_WINDOW_YEARS after the election. Outcomes are
# six combined_panel.rds variables, each measured over the identical
# (election_year - 1) -> (election_year + N) window as the treatment, so
# re-running with a different N moves every outcome's window consistently.
#
# See notes/rdd_plan.md for the full design and open scope notes (esp. the
# Bermeo classification's ~50% coverage limitation).
#
# Data:
#   data/elections_database/presidential_elections.dta
#   data/elections_database/parliamentary_elections.dta
#   data/elections_database/parties_database.dta
#   data/elections_database/CPD_V-Party_CSV_v2.zip
#   data/combined_panel.rds
#   data/episodes_bermeo.rds
#   output/ddcg_ert_miss_table.csv
#
# Output:
#   data/rdd_analysis_data.rds -- one row per election
# ==============================================================================

library(tidyverse)
library(haven)
library(here)
library(countrycode)
library(gt)

# Shared instrument/treatment metadata (ILLIBERALISM_VARS, INSTRUMENT_DISPLAY,
# ...) and table helpers, also used by scripts 12-16.
source(here::here("scripts", "rdd_helpers.R"))

data_dir <- here::here("data")
elections_dir <- file.path(data_dir, "elections_database")
out_dir <- here::here("output")

# ------------------------------------------------------------------------------
# Toggles
# ------------------------------------------------------------------------------

# Every toggle is set with `if (!exists(...))` so a driver script (14, 15) can
# source() this file into an environment that already defines some of them and
# have those overrides respected. Running the script standalone is unchanged --
# the defaults below apply.

# "presidential" | "parliamentary" | "both" -- "both" unions the two spines (see Step 1)
if (!exists("ELECTION_TYPE")) ELECTION_TYPE <- "both"
# years post-election to look for a backsliding start
if (!exists("BACKSLIDING_WINDOW_YEARS")) BACKSLIDING_WINDOW_YEARS <- 5
# Which party-level score defines "the illiberal one" of the top 2, and so the
# sign of the running variable. All five are V-Party party-year variables:
#   v2xpa_antiplural  anti-pluralism index          [0,1]
#   v2xpa_popul       populism index                [0,1]
#   ep_galtan         CHES-merged GAL-TAN score
#   v2pariglef_neg    economic left-right, NEGATED so economic LEFT scores high
#   v2paanteli        anti-elitism
# Note the differing scales: v2pariglef/v2paanteli are expert-scale variables,
# NOT [0,1] indices like the v2xpa_* ones, so any absolute cutoff on the score
# (12_rdd_analysis.R's ILLIBERAL_CUTOFF) has to be set per instrument -- see
# the quantile-based option there.
if (!exists("ILLIBERALISM_VAR")) ILLIBERALISM_VAR <- "v2xpa_antiplural"
# "ert" | "ddcg" | "llm" -- see Step 4b below
if (!exists("START_YEAR_SOURCE")) START_YEAR_SOURCE <- "ert"
# "all" | "ddcg_comparable" -- restrict elections to DDCG's coverage window (see
# DDCG_START/DDCG_END) for a robustness sample comparable to DDCG-based prior work
if (!exists("SAMPLE_YEARS")) SAMPLE_YEARS <- "all"

# Does the post-election treatment window open IN the election year, or the
# year after?
#
#   TRUE  (default) -- treatment window is [election_year, election_year + N]
#   FALSE            -- treatment window is (election_year, election_year + N]
#
# TRUE is the default because it is what actually aligns treatment with the
# outcomes. window_change() measures every outcome from (election_year - 1) to
# (election_year + N), so the outcome window has always spanned the election
# year; under FALSE the treatment window did not, and the two were a year out of
# step with each other.
#
# The argument for FALSE, which this script used to hardcode: ERT and DDCG both
# date events to a calendar year only, and elections are spread through the year
# (48% Jan-Jun, 52% Jul-Dec in this spine), so an episode dated to the election
# year may well have begun BEFORE the vote -- reverse causation rather than
# treatment. That ambiguity is real and unresolvable from annual data, but it
# costs 31 treated elections on the ERT arm (142 -> 173, +22%) and 38 on the
# union arm, which is a lot to give up on a first stage this underpowered.
# FALSE is kept so the earlier results stay reproducible.
#
# Whichever way this is set, the pre- and post-election windows stay a clean
# partition with no overlap and no gap: prior_backsliding covers the N years
# immediately before the treatment window opens.
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- TRUE
}
stopifnot(is.logical(TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR))

stopifnot(ELECTION_TYPE %in% c("presidential", "parliamentary", "both"))
stopifnot(ILLIBERALISM_VAR %in% ILLIBERALISM_VARS)
stopifnot(START_YEAR_SOURCE %in% c("ert", "ddcg", "llm"))
stopifnot(SAMPLE_YEARS %in% c("all", "ddcg_comparable"))

# DDCG (Acemoglu et al. 2019) coverage window -- same bounds used in
# 09_ddcg_overlap.R / 10_bermeo_classification.R for cross-source overlap
# comparisons, reused here both for SAMPLE_YEARS == "ddcg_comparable" and
# for Step 6's match-accounting diagnostic below.
DDCG_START <- 1960
DDCG_END <- 2010

# ------------------------------------------------------------------------------
# Step 1 -- load the election spine(s), reshape wide -> long, take top 2
# Presidential and parliamentary are loaded into the SAME schema (election_id,
# election_type, country_text_id, ..., party_id, final_share) so ELECTION_TYPE
# == "both" can just bind_rows() them and run every later step -- scoring,
# running variable, treatment/outcome merge, heterogeneity -- completely
# unmodified. Only the reshape (wide -> long) and any type-specific validity
# filtering happen inside the two loaders below.
# ------------------------------------------------------------------------------

# country_abb is NOT ISO3 for most countries (53% of codes in this data
# mismatch combined_panel's country_text_id -- e.g. "ROK" for South Korea,
# "GMY" for Germany, "FRN" for France -- it looks like a Correlates-of-War-
# style abbreviation scheme, not ISO3; only a coincidental few, like "AFG",
# happen to match both). Derive the real ISO3 code from country_cow instead,
# which resolves 69 of 73 originally-mismatched codes correctly. Of the
# remaining 4: Kiribati/Palau are likely just outside combined_panel's
# coverage, not a code problem. COW 340 (Serbia and Montenegro, country_abb
# "SCG") is genuinely ambiguous with no 1:1 modern ISO3 mapping -- left
# unmatched, same convention as this repo's DDCG wbcode crosswalk. COW 345
# and 678 are NOT actually ambiguous in this data specifically, even though
# countrycode() can't resolve them in general (COW 345 also covers
# historical Yugoslavia; COW 678 covers pre/post-unification Yemen) -- every
# row here already carries an unambiguous country_abb (spot-checked: 345 is
# "SRB"/Serbia for every row but one 1945 parliamentary election coded
# "YUG"/Yugoslavia; 678 is "YAR" continuously across both pre- and
# post-1990-unification years, all labeled Country == "Yemen", matching
# ISO3 YEM's own continuous coverage of that same lineage). Override with
# country_abb for these two specific cases rather than leaving real
# elections stranded on a package limitation.
derive_country_text_id <- function(country_cow, country_abb) {
  case_when(
    country_cow == 678 ~ "YEM",
    country_cow == 345 & country_abb == "SRB" ~ "SRB",
    TRUE ~ countrycode(country_cow, origin = "cown", destination = "iso3c")
  )
}

load_presidential_candidates <- function() {
  elections_raw <- read_dta(file.path(
    elections_dir,
    "presidential_elections.dta"
  ))

  elections_raw |>
    mutate(election_id = paste0("pres_", row_number())) |>
    pivot_longer(
      cols = matches("^(Candidate|Party_ID|Vote_Share1|Vote_Share2)_[0-9]+$"),
      names_to = c(".value", "slot"),
      names_pattern = "^(Candidate|Party_ID|Vote_Share1|Vote_Share2)_([0-9]+)$"
    ) |>
    filter(!is.na(Candidate), nzchar(trimws(Candidate))) |>
    transmute(
      election_id,
      election_type = "presidential",
      country_abb,
      country_text_id = derive_country_text_id(country_cow, country_abb),
      country_name = Country,
      election_year = Year,
      election_date = Date,
      flag_indirect,
      candidate = Candidate,
      # Empty string (not NA) means "ran without a coded party" (e.g.
      # independents like Karzai in Afghanistan 2004) -- normalize to NA now
      # so it fails the V-Dem match cleanly downstream instead of silently
      # not matching an empty-string key. Must NOT be filtered out here,
      # before ranking: dropping no-party candidates before determining the
      # top 2 would let 3rd/4th place candidates masquerade as the top 2
      # whenever the actual top-2 includes an independent.
      party_id = na_if(trimws(Party_ID), ""),
      # Final-round share: round 2 if the election went to a runoff, else round 1.
      final_share = coalesce(Vote_Share2, Vote_Share1)
    ) |>
    filter(!is.na(final_share))
}

# Unlike presidential (where the top-2 are always two distinct named
# candidates, independents included), a party-level "winner" can be
# genuinely ambiguous here: 76 elections have an exact tie for the top seat
# share across 2+ parties, and hundreds more report a top-2 seat share as a
# single bracket-notation coalition entry (e.g.
# "[PF4393/PF4386/PF4857]") -- an "Electoral Turnovers"-paper convention
# for jointly-contesting alliances, with no resolution appendix in this
# repo (see notes/rdd_plan.md's original scope note). For this first pass,
# rather than guess which single party should be credited, drop any
# election where the top seat share is tied, or where EITHER of the top-2
# (by share) is a bracket-notation coalition entry rather than one
# distinct, unambiguous party_id.
load_parliamentary_candidates <- function() {
  elections_raw <- read_dta(file.path(
    elections_dir,
    "parliamentary_elections.dta"
  ))

  long <- elections_raw |>
    mutate(election_id = paste0("parl_", row_number())) |>
    pivot_longer(
      cols = matches("^(Party_ID|Seat_Share)_[0-9]+$"),
      names_to = c(".value", "slot"),
      names_pattern = "^(Party_ID|Seat_Share)_([0-9]+)$"
    ) |>
    filter(!is.na(Seat_Share)) |>
    transmute(
      election_id,
      election_type = "parliamentary",
      country_abb,
      country_text_id = derive_country_text_id(country_cow, country_abb),
      country_name = Country,
      election_year = Year,
      election_date = Date,
      flag_indirect = NA_real_, # not applicable to parliamentary; kept for schema parity
      candidate = NA_character_,
      party_id = na_if(trimws(Party_ID), ""),
      final_share = Seat_Share
    ) |>
    filter(!is.na(final_share))

  # Check both the top-2 (by share, matching the generic top-2 selection
  # later on), not just the winner: a bracket-notation coalition entry can
  # land in either slot, and a coalition runner-up is just as unresolvable
  # as a coalition winner for identifying "which single party" occupies
  # that spot.
  top2_check <- long |>
    group_by(election_id) |>
    mutate(
      top_share = max(final_share),
      n_tied_for_first = sum(final_share == top_share)
    ) |>
    slice_max(final_share, n = 2, with_ties = FALSE) |>
    summarise(
      n_tied_for_first = first(n_tied_for_first),
      winner_party_id = party_id[which.max(final_share)],
      any_bracket = any(str_detect(party_id, "\\["), na.rm = TRUE),
      .groups = "drop"
    )

  distinct_winner_elections <- top2_check |>
    filter(n_tied_for_first == 1, !is.na(winner_party_id), !any_bracket) |>
    pull(election_id)

  n_before <- n_distinct(long$election_id)
  long <- long |> filter(election_id %in% distinct_winner_elections)
  cat(sprintf(
    "Step 1 (parliamentary): dropped %d / %d elections without a single distinct-party winner (tied top seat share, or a bracket-notation coalition entry in either top-2 slot)\n",
    n_before - n_distinct(long$election_id),
    n_before
  ))
  long
}

candidates_long <- bind_rows(
  if (ELECTION_TYPE %in% c("presidential", "both")) {
    load_presidential_candidates()
  },
  if (ELECTION_TYPE %in% c("parliamentary", "both")) {
    load_parliamentary_candidates()
  }
)

# Drop candidate-rows whose country code countrycode() couldn't resolve to a
# single ISO3 -- either genuinely unmatched (Kiribati/Palau, likely just
# outside combined_panel's coverage) or ambiguous (COW 340/345/678 =
# Serbia and Montenegro/Serbia/Yemen, whose COW-to-ISO3 mapping depends on
# a date countrycode() isn't given here, so it correctly refuses to guess).
# These can never merge onto combined_panel/episode_years, so drop them up
# front rather than letting them silently fail every downstream join.
n_rows_before_cc <- nrow(candidates_long)
candidates_long <- candidates_long |> filter(!is.na(country_text_id))
n_rows_dropped_cc <- n_rows_before_cc - nrow(candidates_long)
cat(sprintf(
  "Step 1: dropped %d / %d candidate-rows with an unresolved country code (country_cow -> ISO3 failed or ambiguous)\n",
  n_rows_dropped_cc,
  n_rows_before_cc
))

# SAMPLE_YEARS == "ddcg_comparable": restrict to elections falling inside
# DDCG's own coverage window, for a robustness sample directly comparable
# to prior DDCG-based work (rather than every election combined_panel
# happens to cover, which runs well past 2010).
if (SAMPLE_YEARS == "ddcg_comparable") {
  n_rows_before_sample <- nrow(candidates_long)
  candidates_long <- candidates_long |>
    filter(election_year >= DDCG_START, election_year <= DDCG_END)
  cat(sprintf(
    "Step 1: SAMPLE_YEARS = 'ddcg_comparable' -- dropped %d / %d candidate-rows outside %d-%d\n",
    n_rows_before_sample - nrow(candidates_long),
    n_rows_before_sample,
    DDCG_START,
    DDCG_END
  ))
}

n_elections_raw <- n_distinct(candidates_long$election_id)

# Rank by share over ALL candidates first -- filtering out no-party
# candidates before ranking would let a 3rd/4th place candidate masquerade
# as "top 2" whenever the true top-2 includes an independent (e.g. Karzai,
# Afghanistan 2004, ran with no coded party and would otherwise have been
# dropped before the ranking step, promoting Qanuni/Mohaqiq into his place).
top2 <- candidates_long |>
  group_by(election_id) |>
  slice_max(final_share, n = 2, with_ties = FALSE) |>
  mutate(rank = row_number()) |>
  ungroup() |>
  filter(rank <= 2)

n_elections_top2 <- top2 |>
  group_by(election_id) |>
  filter(n() == 2) |>
  ungroup() |>
  pull(election_id) |>
  n_distinct()

n_elections_top2_scoreable <- top2 |>
  group_by(election_id) |>
  filter(n() == 2, all(!is.na(party_id))) |>
  ungroup() |>
  pull(election_id) |>
  n_distinct()

n_by_type <- candidates_long |>
  group_by(election_type) |>
  summarise(n = n_distinct(election_id), .groups = "drop") |>
  summarise(
    label = paste(sprintf("%d %s", n, election_type), collapse = " + ")
  ) |>
  pull(label)

cat(sprintf(
  "Step 1: %d elections loaded (%s); %d have a valid top-2 by share; of those, %d have a party ID for both (elections where the top-2 includes an independent/no-party candidate are excluded here, not silently backfilled with the 3rd-place finisher)\n",
  n_elections_raw,
  n_by_type,
  n_elections_top2,
  n_elections_top2_scoreable
))

top2 <- top2 |>
  semi_join(
    top2 |> count(election_id) |> filter(n == 2),
    by = "election_id"
  )

# ------------------------------------------------------------------------------
# Step 2 -- identify the more-illiberal of the top 2
# Party_ID (Party Facts format, e.g. "PF7049") -> parties_database.party_id
# (same string format) -> vdem_id_1 -> V-Party's v2paid, matched to the
# closest available V-Party year for that party.
# ------------------------------------------------------------------------------

# inception/dissolution are "YYYY-MM-DD" strings with "00" placeholders for
# an unknown month/day (e.g. "1990-00-00"). Parse to real Dates -- comparing
# full dates (not just years) against the election's actual date matters:
# see PF7655 (South Korea) below, where two genealogy stages both "cover"
# 1997 at year granularity, but the true election date (18 Dec 1997) falls
# only within the later stage (inception 21 Nov 1997), not the earlier one
# (dissolution 21 Nov 1997) -- a same-day transition invisible at the
# year level. Unknown month/day defaults to the start (Jan 1) of a range
# for inception and the end (Dec 28, safe across all months) for
# dissolution, so a party isn't spuriously excluded from years it likely
# still covers. Genuinely empty strings get a far sentinel date instead of
# NA so date comparisons behave like -Inf/Inf.
parse_party_dates <- function(x, end_of_range = FALSE) {
  empty <- is.na(x) | !nzchar(x)
  yr <- as.integer(substr(x, 1, 4))
  mo <- as.integer(substr(x, 6, 7))
  dy <- as.integer(substr(x, 9, 10))
  mo[is.na(mo) | mo == 0] <- if (end_of_range) 12L else 1L
  dy[is.na(dy) | dy == 0] <- if (end_of_range) 28L else 1L
  out <- as.Date(ifelse(empty, NA, sprintf("%04d-%02d-%02d", yr, mo, dy)))
  out[empty] <- if (end_of_range) {
    as.Date("9999-12-31")
  } else {
    as.Date("1000-01-01")
  }
  out
}

parties_database <- read_dta(file.path(
  elections_dir,
  "parties_database.dta"
)) |>
  select(party_id, name, vdem_id_1, vdem_id_2, inception, dissolution) |>
  filter(!is.na(vdem_id_1)) |>
  # 20 rows carry a second, non-null vdem_id_2 -- an earlier predecessor/
  # alliance-partner V-Party coding for the same lineage stage (e.g.
  # PF36/N-VA: vdem_id_1=36 for 2010+, vdem_id_2=756 for a 2007 electoral
  # alliance; PF1508/SDSM-Macedonia: vdem_id_1=1508 covers 1998-2016,
  # vdem_id_2=2448 covers only 1994). Which id applies to a given election
  # isn't resolvable from parties_database alone, but vdem_id_1 alone is
  # still a perfectly good, unambiguous single-id match for elections
  # falling in ITS coverage window (as most do -- vdem_id_2 is typically a
  # short, chronologically distant earlier stage). Dropping these rows
  # entirely (rather than just discarding vdem_id_2) was overly costly: it
  # silently lost otherwise-cleanly-matchable elections like North
  # Macedonia's 1999/2004 (both resolve correctly via vdem_id_1 alone,
  # since vdem_id_2's one 1994 data point would never be the closer match
  # anyway). Discard vdem_id_2 and keep the row as an ordinary vdem_id_1
  # match instead of dropping it.
  select(-vdem_id_2) |>
  mutate(
    inception_date = parse_party_dates(inception, end_of_range = FALSE),
    dissolution_date = parse_party_dates(dissolution, end_of_range = TRUE)
  )

# 40 party_id values have more than one genealogy_rank (the party renamed/
# re-registered over time, each stage with its own vdem_id_1 and
# inception/dissolution date range) -- a plain join on party_id alone
# produces a many-to-many fan-out, so pick whichever stage's date range
# covers the election DATE (closest range if none covers it exactly).
# Matching must use the full date, not just the year: PF7655 (South Korea)
# has two stages whose ranges both "cover" 1997 at year granularity (DLP
# dissolves 1997-11-21, GNP's successor stage starts the same day), so a
# year-only comparison ties and arbitrarily picks the first row. Without
# an actual election date to break a multi-row tie like this, don't guess
# -- drop the match rather than falling back to a mid-year default.
match_party_lineage <- function(pid, target_date) {
  cand <- parties_database |> filter(party_id == pid)
  if (nrow(cand) == 0) {
    return(tibble(vdem_id_1 = NA_real_))
  }
  if (nrow(cand) == 1) {
    return(cand |> select(vdem_id_1))
  }
  if (is.na(target_date)) {
    return(tibble(vdem_id_1 = NA_real_))
  }
  covering <- cand |>
    filter(inception_date <= target_date, target_date <= dissolution_date)
  if (nrow(covering) >= 1) {
    return(covering |> slice(1) |> select(vdem_id_1))
  }
  cand |>
    mutate(
      dist = pmin(
        abs(as.numeric(inception_date - target_date)),
        abs(as.numeric(dissolution_date - target_date))
      )
    ) |>
    slice_min(dist, n = 1, with_ties = FALSE) |>
    select(vdem_id_1)
}

vparty_con <- unz(
  file.path(elections_dir, "CPD_V-Party_CSV_v2.zip"),
  "CPD_V-Party_CSV_v2/V-Dem-CPD-Party-V2.csv"
)
v_party <- read_csv(vparty_con, show_col_types = FALSE) |>
  select(
    v2paid,
    year,
    v2xpa_antiplural,
    v2xpa_popul,
    ep_galtan,
    v2pariglef,
    v2paanteli
  ) |>
  filter(!is.na(v2paid)) |>
  # v2pariglef runs right-positive in V-Party ("economic left-right position",
  # higher = further right). Every other candidate instrument here is oriented
  # so that HIGHER = the side hypothesized to erode democracy, and the ask is
  # to treat economic LEFT as the high end, so negate it. Carried as its own
  # column (rather than flipping v2pariglef in place) so the raw variable stays
  # available and the orientation is visible at every use site.
  mutate(v2pariglef_neg = -v2pariglef)

# Matching each top2 row's party_id to a V-Party score is a two-stage
# lookup: party_id -> parties_database row -> vdem_id_1 -> V-Party. The
# only ambiguity is in the first stage -- 40 party_ids have more than one
# genealogy_rank row (the party renamed/re-registered over time) --
# resolved above by match_party_lineage() using the election DATE against
# each row's inception/dissolution range, or dropped if that date is
# missing. (parties_database rows with a second vdem_id_2 are dropped
# entirely earlier, so there's no second-stage ambiguity to resolve here.)
lineage_matches <- map2_dfr(
  top2$party_id,
  top2$election_date,
  match_party_lineage
)
top2_with_vdem <- bind_cols(top2, lineage_matches)

n_matched_vdem <- sum(!is.na(top2_with_vdem$vdem_id_1))
cat(sprintf(
  "Step 2a: %d / %d top-2 candidate-rows matched to a V-Dem party ID via parties_database\n",
  n_matched_vdem,
  nrow(top2_with_vdem)
))

# V-Party is party-year grain, not necessarily covering every year -- pick
# the most recent row at or before the election year (never a later one:
# the illiberalism score is meant to characterize the party AS OF the
# election, not with the benefit of its post-election trajectory).
match_vparty <- function(vdem_id_1, target_year) {
  empty_result <- tibble(
    v2xpa_antiplural = NA_real_,
    v2xpa_popul = NA_real_,
    ep_galtan = NA_real_,
    v2pariglef_neg = NA_real_,
    v2paanteli = NA_real_
  )
  if (is.na(vdem_id_1)) {
    return(empty_result)
  }
  cand <- v_party |> filter(v2paid == vdem_id_1, year <= target_year)
  if (nrow(cand) == 0) {
    return(empty_result)
  }
  cand |>
    slice_max(year, n = 1, with_ties = FALSE) |>
    select(all_of(ILLIBERALISM_VARS))
}

vparty_scores <- map2_dfr(
  top2_with_vdem$vdem_id_1,
  top2_with_vdem$election_year,
  match_vparty
)

top2_scored <- bind_cols(top2_with_vdem, vparty_scores) |>
  mutate(illiberalism_score = .data[[ILLIBERALISM_VAR]])

n_matched_score <- sum(!is.na(top2_scored$illiberalism_score))
cat(sprintf(
  "Step 2b: %d / %d top-2 candidate-rows matched to a V-Party score (var = %s)\n",
  n_matched_score,
  nrow(top2_scored),
  ILLIBERALISM_VAR
))

# ------------------------------------------------------------------------------
# Step 3 -- running variable
# ------------------------------------------------------------------------------

top2_for_scoring <- top2_scored |>
  group_by(election_id) |>
  filter(sum(!is.na(illiberalism_score)) == 2) |> # both candidates must have a score
  ungroup()

# Carry EVERY candidate instrument's score for both top-2 members through to the
# election level, not just the active ILLIBERALISM_VAR. Script 16 needs all of
# them to compare how much the instrument definitions actually differ, and
# recomputing the party match there would duplicate ~200 lines of join logic.
# Ordered by SHARE (winner first), not by any score, so "__winner"/"__loser"
# mean the same thing for every instrument. Ties are impossible for
# parliamentary (dropped in Step 1) and vanishingly rare for presidential, but
# arrange() breaking a tie arbitrarily is harmless here since both members then
# carry the same share.
all_scores_wide <- top2_for_scoring |>
  group_by(election_id) |>
  arrange(desc(final_share), .by_group = TRUE) |>
  summarise(
    across(
      all_of(ILLIBERALISM_VARS),
      list(winner = ~ .x[1], loser = ~ .x[2]),
      .names = "{.col}__{.fn}"
    ),
    .groups = "drop"
  )

elections_scored <- top2_for_scoring |>
  group_by(election_id) |>
  arrange(desc(illiberalism_score), .by_group = TRUE) |>
  summarise(
    election_type = first(election_type),
    country_abb = first(country_abb), # raw source code, kept for display only
    country_text_id = first(country_text_id), # ISO3, used for every downstream match
    country_name = first(country_name),
    election_year = first(election_year),
    election_date = first(election_date),
    flag_indirect = first(flag_indirect),
    illiberal_share = first(final_share),
    other_share = last(final_share),
    illiberal_score = first(illiberalism_score),
    other_score = last(illiberalism_score),
    .groups = "drop"
  ) |>
  mutate(
    running_var = illiberal_share - other_share,
    score_gap = illiberal_score - other_score
  ) |>
  left_join(all_scores_wide, by = "election_id")

# score_gap_z: score_gap standardized against the variability of
# illiberalism_score itself, pooling BOTH top-2 members across EVERY
# scored election in that country (party-years), not the distribution of
# election-level gaps. This is equivalent to z-scoring each party's raw
# score against its country's own score distribution and taking the
# difference of those two z-scores for the top-2 (the country mean
# cancels out of that difference, leaving score_gap / sd_country_score).
# Standardizing against the GAP's own sample SD (the prior approach) has a
# mechanical ceiling: with a sample z-score, an n-election country can
# never exceed |z| ~ sqrt(n-1)-ish regardless of how large its actual
# gaps are (e.g. a 2-election country can never exceed |z| = 0.71), so
# raising a z-threshold silently excluded sparse-coverage countries no
# matter how genuinely high-contrast their races were. Pooling party-years
# instead uses ~2x the observations per country and isn't self-referential
# to the gap variable, so it doesn't have that hard ceiling.
# One row per COUNTRY, which is what the comment above describes and what the
# join below requires. Grouping by (country_text_id, party_id) -- as this did
# previously -- computed each PARTY's own within-party SD and then left_join()ed
# it on country_text_id alone, fanning every election out against every party in
# its country (1,347 elections became 6,602 rows) and attaching an arbitrary
# party's SD to each copy. Both the stated intent ("pooling BOTH top-2 members
# across every scored election in that country") and the ceiling argument that
# motivated it require the country-level pool.
country_score_sd <- top2_for_scoring |>
  group_by(country_text_id) |>
  summarise(
    n_party_years = sum(!is.na(illiberalism_score)),
    sd_country_score = sd(illiberalism_score, na.rm = TRUE),
    .groups = "drop"
  )

elections_scored <- elections_scored |>
  left_join(country_score_sd, by = "country_text_id", relationship = "many-to-one") |>
  mutate(
    score_gap_z = if_else(
      n_party_years > 1 & sd_country_score > 0,
      score_gap / sd_country_score,
      NA_real_
    )
  ) |>
  select(-n_party_years, -sd_country_score)

n_elections_scored <- nrow(elections_scored)
cat(sprintf(
  "Step 3: %d elections have a valid running variable (both top-2 candidates matched to a V-Party score)\n",
  n_elections_scored
))

# ------------------------------------------------------------------------------
# Step 4b -- backsliding start-year robustness toggle
# ERT's aut_ep_start_year is a Polyarchy-threshold-crossing artifact, not a
# real calendar date. START_YEAR_SOURCE swaps in a more precise alternative
# where available, always falling back to aut_ep_start_year otherwise.
# ------------------------------------------------------------------------------

episodes_bermeo <- readRDS(file.path(data_dir, "episodes_bermeo.rds"))

ddcg_years <- read_csv(
  file.path(out_dir, "ddcg_ert_miss_table.csv"),
  col_types = cols(ddcg_event_years = col_character(), .default = col_guess())
) |>
  select(country_name, ep_start, ep_end, ddcg_event_years) |>
  mutate(
    ddcg_start_year = map_dbl(ddcg_event_years, function(x) {
      if (is.na(x)) {
        return(NA_real_)
      }
      years <- as.numeric(str_extract_all(x, "\\d{4}")[[1]])
      if (length(years) == 0) NA_real_ else min(years)
    })
  ) |>
  select(country_name, ep_start, ep_end, ddcg_start_year)

# ERT backsliding episodes, use START_YEAR_SOURCE to override the ert episode start
# year with the DDCG or LLM coded start year for the episode
episode_years <- episodes_bermeo |>
  left_join(
    ddcg_years,
    by = c("country_name", "ep_start" = "ep_start", "ep_end" = "ep_end")
  ) |>
  mutate(
    llm_start_year = map_dbl(precise_start_date, function(x) {
      if (is.na(x)) {
        return(NA_real_)
      }
      years <- as.numeric(str_extract_all(x, "\\d{4}")[[1]])
      if (length(years) == 0) NA_real_ else min(years)
    }),
    start_year_used = case_when(
      START_YEAR_SOURCE == "ert" ~ ep_start,
      START_YEAR_SOURCE == "ddcg" ~ coalesce(ddcg_start_year, ep_start),
      START_YEAR_SOURCE == "llm" ~ coalesce(llm_start_year, ep_start)
    ),
    used_override = start_year_used != ep_start
  )

n_override <- sum(episode_years$used_override, na.rm = TRUE)
cat(sprintf(
  "Step 4b: START_YEAR_SOURCE = '%s' -- %d / %d episodes used an overridden start year (rest fell back to ERT's aut_ep_start_year)\n",
  START_YEAR_SOURCE,
  n_override,
  nrow(episode_years)
))

# ------------------------------------------------------------------------------
# Step 4 -- merge to backsliding treatment + all six outcomes, windows
# aligned to the identical (election_year - 1) -> (election_year + N) span.
# ------------------------------------------------------------------------------

combined_panel <- readRDS(file.path(data_dir, "combined_panel.rds"))

# Value of `var` for `country` in `yr`, or NA if that country-year isn't in
# the panel.
#
# This used to be a dplyr::filter() over all ~20k panel rows on every call.
# With 1,347 elections x ~19 (variable, year) lookups that was ~25k full table
# scans per build, and this script is now run 13+ times (once per instrument x
# window combination) -- so index the panel once into a hash keyed on
# "ISO3 year" and make each lookup O(1). Columns are pulled out of the data
# frame into a plain list of vectors up front too, since `[[` on a tibble is
# itself not free at this call count.
panel_row_index <- setNames(
  seq_len(nrow(combined_panel)),
  paste(combined_panel$country_text_id, combined_panel$year)
)
panel_cols <- as.list(combined_panel)

panel_value <- function(country, yr, var) {
  # Single-bracket, not [[ ]]: `[[` on a named vector ERRORS for a name that
  # isn't present, while `[` returns NA -- and a missing country-year is the
  # normal case here (elections run past the panel's coverage on both ends).
  i <- panel_row_index[paste(country, yr)]
  if (is.na(i)) NA_real_ else panel_cols[[var]][[i]]
}

# Change in `var` from (election_year - 1) to (election_year + N); if
# `compounded = TRUE`, instead compound the annual rates over that same span
# (used for inflation, which has no level column to difference).
window_change <- function(
  country,
  election_year,
  var,
  compounded = FALSE,
  log_transform = FALSE
) {
  pre_year <- election_year - 1
  post_year <- election_year + BACKSLIDING_WINDOW_YEARS
  if (!compounded) {
    pre_val <- panel_value(country, pre_year, var)
    post_val <- panel_value(country, post_year, var)
    return(post_val - pre_val) # standard log-difference approximation to percentage growth
  }
  yrs <- (pre_year + 1):post_year
  rates <- map_dbl(yrs, ~ panel_value(country, .x, var))
  if (any(is.na(rates))) {
    return(NA_real_)
  }
  growth_factor <- prod(1 + rates / 100)
  # log_transform: return cumulative log growth (sum of each year's
  # log(1 + rate/100)) instead of the raw cumulative % change. Hyperinflation
  # episodes (e.g. Angola 1992: cumulative inflation of 2,324,357%, Peru
  # 1985: 177,179%) are 4-5 orders of magnitude larger than a typical
  # observation (~40-80%) on the raw scale -- a handful of these points
  # dominate any local regression's slope/variance near the RD cutoff,
  # producing wildly unstable estimates. Logging compresses the right tail
  # the same way Y_gdp_growth's log-difference already does, at the cost of
  # changing units from "% change" to "log points" (still directly
  # interpretable: 0.05 ~ 5% for small values, same log-difference
  # convention used throughout this script).
  if (log_transform) log(growth_factor) else growth_factor - 1
}

# Every outcome is a change over the identical (election_year - 1) ->
# (election_year + N) window as the treatment, so re-running with a different N
# moves treatment and every outcome consistently.
#
# The three growth measures are all log-level differences of a log GDP per
# capita series, so they are directly comparable and can share one plot panel:
#   ln_gdp_pc      Penn World Table 11.0  (the original)
#   ln_gdp_pc_wb   World Bank             (outcomes.dta logGDPc_wb)
#   ln_gdp_pc_imf  IMF WEO                (log of outcomes.dta GDPc_imfweo)
build_outcomes <- function(country, election_year) {
  tibble(
    Y_gdp_growth = window_change(country, election_year, "ln_gdp_pc"),
    Y_gdp_growth_wb = window_change(country, election_year, "ln_gdp_pc_wb"),
    Y_gdp_growth_imf = window_change(country, election_year, "ln_gdp_pc_imf"),
    Y_inflation = window_change(
      country,
      election_year,
      "cpi_inflation",
      compounded = TRUE,
      log_transform = TRUE
    ),
    Y_unemployment = window_change(country, election_year, "unemployment_rate"),
    Y_trade_pct_gdp = window_change(country, election_year, "trade_pct_gdp"),
    Y_top10_share = window_change(country, election_year, "top10_share"),
    Y_gini_disp = window_change(country, election_year, "gini_disp"),
    Y_gini_mkt = window_change(country, election_year, "gini_mkt"),
    # Fiscal. deficit_pct_gdp is oriented in 01f so that positive = LARGER
    # deficit, matching the "higher = worse" direction of the other outcomes.
    Y_debt = window_change(country, election_year, "debt_pct_gdp"),
    Y_deficit = window_change(country, election_year, "deficit_pct_gdp"),
    # V-Dem institutional subcomponents. Note these run the OTHER way from the
    # economic outcomes: higher = more constrained executive = healthier
    # democracy, so a negative RD estimate is the "backsliding" sign here.
    Y_checks_balances = window_change(country, election_year, "checks_balances"),
    Y_jucon = window_change(country, election_year, "v2x_jucon"),
    Y_legcon = window_change(country, election_year, "v2xlg_legcon"),
    Y_hos_power = window_change(country, election_year, "hos_power_linear"),
    Y_hog_power = window_change(country, election_year, "hog_power_linear"),
    Y_hos_power_vdem = window_change(country, election_year, "hos_power_vdem"),
    Y_hog_power_vdem = window_change(country, election_year, "hog_power_vdem")
  )
}

# Backsliding treatment: does this country have an episode whose start_year_used
# falls in (election_year, election_year + N]? Also record prior_backsliding
# and the earliest triggering episode ID.
backsliding_for_election <- function(country, election_year) {
  # The single place the window convention is applied. Both the ERT arm and the
  # DDCG arm below read window_start, so the two sources can never drift onto
  # different conventions -- which would matter, since backsliding_union_Nyr is
  # their pmax() and any apparent power gain would then partly be the wider
  # window rather than the extra events.
  window_start <- if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) {
    election_year
  } else {
    election_year + 1
  }
  window_end <- election_year + BACKSLIDING_WINDOW_YEARS

  # ERT episodes that started in the post-election window
  post_matches <- episode_years |>
    filter(
      country_text_id == country,
      start_year_used >= window_start,
      start_year_used <= window_end
    ) |>
    arrange(start_year_used)
  # ERT episodes that started in the BACKSLIDING_WINDOW_YEARS immediately
  # before the treatment window opens -- a backsliding episode from decades
  # earlier isn't a relevant confound for this specific election. Anchored to
  # window_start rather than election_year so the pre- and post-windows abut
  # exactly, with no year counted twice or skipped, under either convention.
  prior_matches <- episode_years |>
    filter(
      country_text_id == country,
      start_year_used >= window_start - BACKSLIDING_WINDOW_YEARS,
      start_year_used < window_start
    )

  # DDCG (Acemoglu et al. 2019) reversal events in the same post-election
  # window. DDCG's panel only runs 1960-2010, so outside that span this is
  # always 0 -- ddcg_covered below records whether the window was even in
  # scope, so the union treatment's extra power isn't misread as coming from
  # elections DDCG never could have contributed to.
  ddcg_years_in_window <- window_start:window_end
  ddcg_hits <- vapply(
    ddcg_years_in_window,
    function(y) {
      v <- panel_value(country, y, "revevent")
      isTRUE(!is.na(v) && v == 1)
    },
    logical(1)
  )

  # Continuous treatment: how far V-Dem polyarchy FELL over the same
  # (election_year - 1) -> (election_year + N) window, negated so that, like
  # the binary treatments, higher = more backsliding.
  polyarchy_decline <- -window_change(country, election_year, "v2x_polyarchy")

  tibble(
    backsliding_Nyr = as.integer(nrow(post_matches) > 0),
    matched_aut_ep_id = if (nrow(post_matches) > 0) {
      post_matches$aut_ep_id[1]
    } else {
      NA_character_
    },
    prior_backsliding = as.integer(nrow(prior_matches) > 0),
    backsliding_ddcg_Nyr = as.integer(any(ddcg_hits)),
    ddcg_covered = as.integer(
      max(ddcg_years_in_window) >= DDCG_START &&
        min(ddcg_years_in_window) <= DDCG_END
    ),
    polyarchy_decline = polyarchy_decline
  )
}

cat(
  "Step 4: merging backsliding treatment + outcomes onto each election (this loops per election)...\n"
)

treatment_outcomes <- pmap_dfr(
  list(elections_scored$country_text_id, elections_scored$election_year),
  function(country, yr) {
    bind_cols(
      backsliding_for_election(country, yr),
      build_outcomes(country, yr)
    )
  }
)

elections_full <- bind_cols(elections_scored, treatment_outcomes) |>
  # The extended treatment (to-do 4): an ERT autocratization episode OR a DDCG
  # democratic reversal in the post-election window. Computed here rather than
  # made a build-time toggle so a single build serves all three treatment
  # definitions and 12_rdd_analysis.R can switch between them for free.
  mutate(
    backsliding_union_Nyr = pmax(backsliding_Nyr, backsliding_ddcg_Nyr),
    # The continuous polyarchy treatment, binarised: did polyarchy decline over
    # the window AT ALL, regardless of how far. Same quantity as
    # polyarchy_decline at a coarser level of measurement, so the pair isolates
    # what the magnitude information is worth -- the binary version discards it.
    # Strictly positive: an exactly unchanged index is not a decline.
    polyarchy_declined = as.integer(polyarchy_decline > 0)
  )

n_matched_panel <- sum(
  !is.na(elections_full$Y_gdp_growth) | elections_full$backsliding_Nyr == 1
)
cat(sprintf(
  "Step 4: %d treated elections (backsliding_Nyr==1); %d with prior_backsliding==1\n",
  sum(elections_full$backsliding_Nyr),
  sum(elections_full$prior_backsliding)
))
n_added <- sum(
  elections_full$backsliding_union_Nyr == 1 & elections_full$backsliding_Nyr == 0
)
cat(sprintf(
  "Step 4: DDCG reversals add %d treated elections on top of ERT (%d -> %d); %d elections have a window overlapping DDCG's %d-%d coverage at all\n",
  n_added,
  sum(elections_full$backsliding_Nyr),
  sum(elections_full$backsliding_union_Nyr),
  sum(elections_full$ddcg_covered),
  DDCG_START,
  DDCG_END
))
cat(sprintf(
  "Step 4: polyarchy_decline non-missing for %d elections (mean %.4f, sd %.4f); polyarchy_declined == 1 for %d\n",
  sum(!is.na(elections_full$polyarchy_decline)),
  mean(elections_full$polyarchy_decline, na.rm = TRUE),
  sd(elections_full$polyarchy_decline, na.rm = TRUE),
  sum(elections_full$polyarchy_declined, na.rm = TRUE)
))

# ------------------------------------------------------------------------------
# Step 5 -- heterogeneity covariates
# ------------------------------------------------------------------------------

elections_final <- elections_full |>
  left_join(
    episode_years |>
      select(matched_aut_ep_id = aut_ep_id, bermeo_category_primary),
    by = "matched_aut_ep_id"
  )

n_treated <- sum(elections_final$backsliding_Nyr)
n_bermeo_matched <- sum(
  elections_final$backsliding_Nyr == 1 &
    !is.na(elections_final$bermeo_category_primary)
)
cat(sprintf(
  "Step 5: of %d treated elections, %d (%.1f%%) have a Bermeo category attached (coverage is concentrated pre-2010 -- see notes/rdd_plan.md)\n",
  n_treated,
  n_bermeo_matched,
  if (n_treated > 0) 100 * n_bermeo_matched / n_treated else NA_real_
))

# ------------------------------------------------------------------------------
# Save + summary
# ------------------------------------------------------------------------------

# Cache one build per (instrument, window) pair. The treatment definitions and
# every outcome are all computed in a single pass, so this is the only axis the
# build actually varies along -- 12_rdd_analysis.R and the driver scripts read
# these files back rather than rebuilding.
build_dir <- file.path(data_dir, "rdd_build")
dir.create(build_dir, showWarnings = FALSE, recursive = TRUE)
# The window convention changes the treatment columns, so it has to be part of
# the build's identity or a TRUE build and a FALSE build would silently
# overwrite each other. Only the non-default (FALSE) convention gets a suffix,
# keeping the default filenames clean.
build_suffix <- if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "" else "_exclyr"
build_path <- file.path(
  build_dir,
  sprintf(
    "rdd_%s_w%d%s.rds",
    ILLIBERALISM_VAR, BACKSLIDING_WINDOW_YEARS, build_suffix
  )
)
# Record the convention on the object itself, so 12_rdd_analysis.R can read it
# back and label its run folder accordingly rather than having to be told.
attr(elections_final, "includes_election_year") <-
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR
saveRDS(elections_final, build_path)
message("Saved ", build_path)

# Party-level companion: the two top-2 finishers of every scored election, with
# every candidate instrument's score. 16_instrument_overlap.R needs this grain
# (a party-year observation, not an election) for its Jaccard heatmaps, and
# re-deriving it there would mean duplicating the whole Party Facts -> V-Dem ->
# V-Party matching chain above.
parties_path <- sub("\\.rds$", "_parties.rds", build_path)
saveRDS(
  top2_for_scoring |>
    select(
      election_id, election_type, country_text_id, country_name, election_year,
      party_id, candidate, final_share, all_of(ILLIBERALISM_VARS)
    ),
  parties_path
)
message("Saved ", parties_path)

# Also keep writing the historical default location, so anything still pointing
# at it (and a standalone run of 12_rdd_analysis.R) keeps working.
if (
  ILLIBERALISM_VAR == "v2xpa_antiplural" && BACKSLIDING_WINDOW_YEARS == 5 &&
    ELECTION_TYPE == "both" && START_YEAR_SOURCE == "ert" &&
    SAMPLE_YEARS == "all" && TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR
) {
  saveRDS(elections_final, file.path(data_dir, "rdd_analysis_data.rds"))
  message("Saved data/rdd_analysis_data.rds (default configuration)")
}

cat("\n=== Summary ===\n")
cat(sprintf("Election type: %s\n", ELECTION_TYPE))
cat(sprintf("Backsliding window: %d years\n", BACKSLIDING_WINDOW_YEARS))
cat(sprintf("Illiberalism variable: %s\n", ILLIBERALISM_VAR))
cat(sprintf("Start-year source: %s\n", START_YEAR_SOURCE))
cat(sprintf(
  "Treatment window: [%s, election_year + %d]\n",
  if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) {
    "election_year"
  } else {
    "election_year + 1"
  },
  BACKSLIDING_WINDOW_YEARS
))
cat(sprintf("N elections (final): %d\n", nrow(elections_final)))
cat(sprintf("N treated (backsliding_Nyr==1): %d\n", n_treated))
cat("\nNon-missing counts per outcome:\n")
outcome_cols <- grep("^Y_", names(elections_final), value = TRUE)
print(colSums(!is.na(elections_final[outcome_cols])))
cat(sprintf(
  "\nBermeo-match coverage (of treated): %.1f%%\n",
  100 * n_bermeo_matched / n_treated
))

# ------------------------------------------------------------------------------
# Step 6 -- diagnostic: ERT episode match accounting
# Of the 131 ERT autocratization episodes overlapping DDCG's 1960-2010
# window (same overlap definition as 09_ddcg_overlap.R / 10_bermeo_
# classification.R: ep_start <= 2010, ep_end >= 1960), only a fraction
# surface as a matched_aut_ep_id in elections_final. This walks every such
# episode through the pipeline (election exists? scored? backsliding_Nyr
# flagged? for THIS episode specifically?) to find exactly where it drops
# out, mirroring the miss-table pattern from scripts 07-10.
# ------------------------------------------------------------------------------

episodes_in_window <- episode_years #|>
#filter(ep_start <= DDCG_END, ep_end >= DDCG_START)

matched_ids <- unique(
  elections_full$matched_aut_ep_id[!is.na(elections_full$matched_aut_ep_id)]
)

# Mirrors backsliding_for_election()'s post_matches condition exactly, from the
# episode's point of view: for this episode to be a candidate treatment for an
# election, election_year must be in [start_year_used - N, start_year_used) --
# or [start_year_used - N, start_year_used] when the treatment window includes
# the election year, since an episode starting in the election year itself then
# counts.
classify_episode <- function(ep_id, country, start_year_used) {
  if (ep_id %in% matched_ids) {
    return("Matched")
  }

  win_lo <- start_year_used - BACKSLIDING_WINDOW_YEARS
  win_hi <- start_year_used

  candidate_elections <- candidates_long |>
    filter(
      country_text_id == country,
      election_year >= win_lo,
      if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) {
        election_year <= win_hi
      } else {
        election_year < win_hi
      }
    ) |>
    distinct(election_year)

  if (nrow(candidate_elections) == 0) {
    if (!(country %in% candidates_long$country_text_id)) {
      return("Country has no elections in the spine at all")
    }
    return("No election in the pre-episode window")
  }

  scored_in_window <- elections_scored |>
    filter(
      country_text_id == country,
      election_year %in% candidate_elections$election_year
    )

  if (nrow(scored_in_window) == 0) {
    return(
      "Election(s) in window never scored (top-2/party/V-Party match failed)"
    )
  }

  full_in_window <- elections_full |>
    filter(
      country_text_id == country,
      election_year %in% scored_in_window$election_year
    )

  if (any(full_in_window$backsliding_Nyr == 1)) {
    return("Shadowed by a closer episode in the same election's window")
  }

  "Unexplained (should not occur)"
}

ert_miss_table <- episodes_in_window |>
  rowwise() |>
  mutate(
    miss_reason = classify_episode(aut_ep_id, country_text_id, start_year_used)
  ) |>
  ungroup() |>
  select(
    country_name,
    aut_ep_id,
    ep_start,
    ep_end,
    start_year_used,
    outcome_label,
    prch_label,
    miss_reason
  ) |>
  arrange(ep_start)

n_episodes_window <- nrow(ert_miss_table)
n_matched_window <- sum(ert_miss_table$miss_reason == "Matched")

cat(sprintf(
  "\n=== Step 6: ERT episode match accounting (%d-%d overlap window) ===\n",
  DDCG_START,
  DDCG_END
))
cat(sprintf(
  "%d / %d episodes matched to an election (%.1f%%)\n",
  n_matched_window,
  n_episodes_window,
  100 * n_matched_window / n_episodes_window
))
cat("\nBreakdown of unmatched episodes by reason:\n")
print(
  ert_miss_table |>
    filter(miss_reason != "Matched") |>
    count(miss_reason, sort = TRUE)
)

# ---- gt miss table -----------------------------------------------------------

apply_table_style <- function(gt_tbl) {
  gt_tbl |>
    tab_options(
      table.font.size = px(11),
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

election_type_label <- if (ELECTION_TYPE == "both") {
  "presidential + parliamentary"
} else {
  ELECTION_TYPE
}

gt_ert_miss_table <- ert_miss_table |>
  gt() |>
  tab_header(
    title = paste0(
      "ERT autocratization episodes vs. the ",
      election_type_label,
      "-election RD spine"
    ),
    subtitle = paste0(
      "Episodes overlapping ",
      DDCG_START,
      "-",
      DDCG_END,
      "; window = ",
      BACKSLIDING_WINDOW_YEARS,
      " yr post-election"
    )
  ) |>
  cols_label(
    country_name = "Country",
    aut_ep_id = "Episode ID",
    ep_start = "Start",
    ep_end = "End",
    start_year_used = "Start yr used",
    outcome_label = "ERT outcome",
    prch_label = "From democracy",
    miss_reason = "Match status / miss reason"
  ) |>
  cols_align(
    align = "center",
    columns = c(ep_start, ep_end, start_year_used, prch_label)
  ) |>
  tab_style(
    style = list(cell_fill(color = "#c7e9c0"), cell_text(weight = "bold")),
    locations = cells_body(
      columns = miss_reason,
      rows = miss_reason == "Matched"
    )
  ) |>
  tab_style(
    style = cell_fill(color = "#eeeeee"),
    locations = cells_body(
      columns = miss_reason,
      rows = miss_reason != "Matched"
    )
  ) |>
  opt_row_striping() |>
  apply_table_style()

# Written into the per-build folder rather than a flat output/ file, so
# different instrument/window builds don't overwrite each other's accounting.
miss_dir <- file.path(out_dir, "runs", "_builds", sprintf(
  "%s_w%d%s",
  ILLIBERALISM_VAR,
  BACKSLIDING_WINDOW_YEARS,
  build_suffix
))
dir.create(miss_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(ert_miss_table, file.path(miss_dir, "ert_miss_table.csv"))
gtsave(gt_ert_miss_table, file.path(miss_dir, "ert_miss_table.html"))
message("Saved ", file.path(miss_dir, "ert_miss_table.csv/.html"))
