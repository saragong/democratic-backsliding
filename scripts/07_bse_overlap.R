# ==============================================================================
# BSE overlap analysis
# How do our ERT autocratization episodes relate to the treatment definitions
# in Boese-Schlosser & Eberhardt (CEPR 2025)?
#
# BSE define democratic breakdown using five families of indicator:
#   (A) ROW     — row_regch_event == -1 (ERT panel variable)
#   (B) ERT     — reg_trans == -1 (ERT panel variable)
#   (C) Polyarchy thresholds  — 5 variants of v2x_polyarchy crossing
#   (D) LibDem thresholds     — 5 variants of v2x_libdem crossing
#   (E) LibComp thresholds    — 5 variants of v2x_liberal crossing
#
# For (C)-(E), BSE construct thresholds as mean ± {1/4, 1/8} sd of the
# respective index over 1961-2023 (all countries). A country enters treatment
# in the year its index crosses from above to below the threshold.
#
# BSE sample: ~169 countries, 1999-2023.
#
# Outputs
#   output/bse_overlap_summary.csv  — one row per BSE definition
#   output/bse_novel_episodes.csv   — ERT outcome==2 episodes BSE cannot capture
#   output/bse_timing_gap.csv       — years from ERT ep start to BSE breakdown year
# ==============================================================================

library(tidyverse)
library(vdemdata)
library(gt)
library(here)

data_dir <- here::here("data")
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE)

BSE_START <- 1999
BSE_END <- 2023
MOMENT_WINDOW <- 1961:2023 # BSE: compute mean/SD over this window
EP_MATCH_TOL <- 1 # tolerance (years) when matching breakdown year to episode window

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

ert <- readRDS(file.path(data_dir, "ert_episodes.rds"))

vdem_panel <- vdemdata::vdem |>
  select(
    country_text_id,
    country_name,
    year,
    v2x_polyarchy,
    v2x_libdem,
    v2x_liberal,
    v2x_regime
  ) |>
  filter(!is.na(country_text_id), !is.na(year))

# ERT autocratization episodes originating from a democracy (aut_ep_prch == 1),
# comparable to BSE which only studies countries that were democracies pre-collapse.
aut_eps <- ert |>
  filter(aut_ep == 1, aut_ep_prch == 1) |>
  group_by(aut_ep_id) |>
  summarise(
    country_text_id = first(country_text_id),
    country_name = first(country_name),
    ep_start = first(aut_ep_start_year),
    ep_end = first(aut_ep_end_year),
    outcome_agg = first(aut_ep_outcome_agg),
    poly_start = v2x_polyarchy[year == first(aut_ep_start_year)][1],
    poly_end = v2x_polyarchy[year == first(aut_ep_end_year)][1],
    .groups = "drop"
  ) |>
  mutate(
    outcome_label = case_when(
      outcome_agg == 1 ~ "Breakdown (BSE territory)",
      outcome_agg == 2 ~ "Erosion only, no breakdown (novel territory)",
      outcome_agg == 3 ~ "Autocracy deepening",
      outcome_agg == 4 ~ "Uncertain / ongoing"
    )
  )

# ------------------------------------------------------------------------------
# Helper: threshold-crossing treatment events for a continuous V-Dem index
# Returns one row per (country, year, threshold_label) crossing event.
# ------------------------------------------------------------------------------

threshold_events <- function(
  panel,
  index_col,
  moment_window,
  bse_start,
  bse_end
) {
  stats <- panel |>
    filter(year %in% moment_window, !is.na(.data[[index_col]])) |>
    summarise(mu = mean(.data[[index_col]]), sigma = sd(.data[[index_col]]))

  thresholds <- tibble(
    label = c("M-1/4sd", "M-1/8sd", "Mean", "M+1/8sd", "M+1/4sd"),
    theta = stats$mu + stats$sigma * c(-1 / 4, -1 / 8, 0, 1 / 8, 1 / 4)
  )

  map_dfr(seq_len(nrow(thresholds)), function(i) {
    theta <- thresholds$theta[[i]]
    label <- thresholds$label[[i]]

    panel |>
      filter(!is.na(.data[[index_col]])) |>
      arrange(country_text_id, year) |>
      group_by(country_text_id) |>
      mutate(
        above = .data[[index_col]] >= theta,
        above_lag = lag(above),
        crossing = !above & above_lag # democracy -> autocracy
      ) |>
      ungroup() |>
      filter(crossing, year >= bse_start, year <= bse_end) |>
      transmute(
        country_text_id,
        country_name,
        year,
        index_name = index_col,
        threshold_label = label,
        theta,
        index_val = .data[[index_col]]
      )
  })
}

# ------------------------------------------------------------------------------
# Build all BSE treatment event sets (one row per breakdown event)
# ------------------------------------------------------------------------------

row_events <- ert |>
  filter(row_regch_event == -1, year >= BSE_START, year <= BSE_END) |>
  transmute(country_text_id, country_name, year, bse_definition = "ROW")

ert_events <- ert |>
  filter(reg_trans == -1, year >= BSE_START, year <= BSE_END) |>
  transmute(country_text_id, country_name, year, bse_definition = "ERT regime")

threshold_all <- bind_rows(
  threshold_events(
    vdem_panel,
    "v2x_polyarchy",
    MOMENT_WINDOW,
    BSE_START,
    BSE_END
  ),
  threshold_events(vdem_panel, "v2x_libdem", MOMENT_WINDOW, BSE_START, BSE_END),
  threshold_events(vdem_panel, "v2x_liberal", MOMENT_WINDOW, BSE_START, BSE_END)
) |>
  mutate(
    bse_definition = paste0(
      case_when(
        index_name == "v2x_polyarchy" ~ "Polyarchy",
        index_name == "v2x_libdem" ~ "LibDem",
        index_name == "v2x_liberal" ~ "LibComp"
      ),
      " (",
      threshold_label,
      ")"
    )
  ) |>
  select(country_text_id, country_name, year, bse_definition)

all_bse_events <- bind_rows(row_events, ert_events, threshold_all)

# Metadata lookup: family grouping, within-family strictness rank, and display label.
# strictness is NA for binary definitions; ranges from -2 (strictest threshold,
# lowest theta, fewest crossings) to +2 (loosest threshold, highest theta, most crossings).
# threshold_display gives a publication-friendly label used in the gt table.
FAMILY_ORDER <- c("ROW", "ERT regime", "Polyarchy", "LibDem", "LibComp")
THRESHOLD_CODES <- c("M-1/4sd", "M-1/8sd", "Mean", "M+1/8sd", "M+1/4sd")
THRESHOLD_DISPLAY <- c(
  "mean − σ/4",
  "mean − σ/8",
  "mean",
  "mean + σ/8",
  "mean + σ/4"
)

definition_metadata <- bind_rows(
  tibble(
    bse_definition = c("ROW", "ERT regime"),
    definition_family = c("ROW", "ERT regime"),
    threshold_display = c("ROW", "ERT regime"),
    strictness = NA_integer_
  ),
  tibble(
    bse_definition = paste0(
      rep(c("Polyarchy", "LibDem", "LibComp"), each = 5),
      " (",
      rep(THRESHOLD_CODES, 3),
      ")"
    ),
    definition_family = rep(c("Polyarchy", "LibDem", "LibComp"), each = 5),
    threshold_display = rep(THRESHOLD_DISPLAY, 3),
    strictness = rep(c(-2L, -1L, 0L, 1L, 2L), 3)
  )
)

# ------------------------------------------------------------------------------
# Match each BSE breakdown event to our ERT autocratization episodes.
#
# Within a single BSE definition, events are unique on (country_text_id, year)
# — a country can cross a threshold at most once per year. We exploit this by
# matching one definition at a time via group_modify, so the join left-side is
# never duplicated. This avoids a cross-definition fan-out and makes the
# many-to-many suppression unnecessary.
#
# For each definition's events we use join_by(between(...)) to match a
# breakdown year directly to ERT episode windows — same pattern as
# 02b_build_event_dataset.R. The right side (aut_eps) can still have multiple
# episodes per country, so slice_min keeps the closest one when that happens.
# ------------------------------------------------------------------------------

aut_eps_bounds <- aut_eps |>
  select(
    country_text_id,
    ep_start,
    ep_end,
    outcome_agg,
    outcome_label,
    poly_start,
    poly_end
  ) |>
  mutate(
    ep_start_match = ep_start - EP_MATCH_TOL,
    ep_end_match = ep_end + EP_MATCH_TOL
  )

match_one_definition <- function(events_df, ...) {
  events_df |>
    left_join(
      aut_eps_bounds,
      by = join_by(country_text_id, between(year, ep_start_match, ep_end_match))
    ) |>
    # between() join can still match one event to multiple overlapping episodes;
    # keep the one whose ep_end is closest to the breakdown year.
    group_by(country_text_id, year) |>
    slice_min(abs(ep_end - year), n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(-ep_start_match, -ep_end_match)
}

bse_matched <- all_bse_events |>
  group_by(bse_definition) |>
  group_modify(match_one_definition) |>
  ungroup()

# ------------------------------------------------------------------------------
# Summary table: one row per BSE definition
# ------------------------------------------------------------------------------

summary_table <- bse_matched |>
  group_by(bse_definition) |>
  summarise(
    n_breakdown_events = n(),
    n_countries = n_distinct(country_text_id),
    n_matched_ert_ep = sum(!is.na(outcome_agg)),
    # The next three columns use ERT's aut_ep_outcome_agg, which classifies
    # episodes against the fixed RoW binary threshold — not against BSE's
    # current definition. A "RoW breakdown" (outcome_agg==1) means the episode
    # crossed from democracy to autocracy on the RoW measure; "erosion only"
    # (outcome_agg==2) means it did not; "ongoing" (outcome_agg==4) means the
    # episode had not resolved when ERT was coded. The three sum to n_matched_ert_ep.
    n_ep_row_breakdown = sum(outcome_agg == 1, na.rm = TRUE),
    n_ep_erosion_only = sum(outcome_agg == 2, na.rm = TRUE),
    n_ep_ongoing = sum(outcome_agg == 4, na.rm = TRUE),
    n_no_ert_match = sum(is.na(outcome_agg)),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Reverse direction: ERT episodes → which BSE definitions capture them?
#
# For each BSE definition, we ask: of all ERT autocratization episodes
# (aut_ep_prch == 1) that overlap the BSE period, how many have at least one
# BSE breakdown event falling inside the episode window (± EP_MATCH_TOL)?
# Episodes with no such event are "invisible" to that definition.
#
# We use the same inequality join_by pattern as novel_episodes below, but
# with the roles swapped: ERT episode bounds are on the LEFT, BSE event year
# is on the RIGHT. group_modify processes one definition at a time so there
# is no cross-definition fan-out.
# ------------------------------------------------------------------------------

ert_eps_in_period <- aut_eps |>
  filter(ep_start <= BSE_END, ep_end >= BSE_START) |>
  mutate(
    ep_start_match = ep_start - EP_MATCH_TOL,
    ep_end_match = ep_end + EP_MATCH_TOL
  )

match_ert_to_definition <- function(events_df, ...) {
  ert_eps_in_period |>
    left_join(
      events_df |> select(country_text_id, bse_year = year),
      by = join_by(
        country_text_id,
        ep_start_match <= bse_year,
        ep_end_match >= bse_year
      )
    ) |>
    group_by(aut_ep_id, country_text_id, ep_start, ep_end, outcome_agg) |>
    summarise(has_bse_match = any(!is.na(bse_year)), .groups = "drop")
}

ert_coverage <- all_bse_events |>
  group_by(bse_definition) |>
  group_modify(match_ert_to_definition) |>
  ungroup()

ert_coverage_summary <- ert_coverage |>
  group_by(bse_definition) |>
  summarise(
    n_ert_ep_in_period = n(),
    n_ert_ep_captured = sum(has_bse_match),
    n_ert_ep_missed = sum(!has_bse_match),
    n_ert_missed_breakdown = sum(
      !has_bse_match & outcome_agg == 1,
      na.rm = TRUE
    ),
    n_ert_missed_erosion = sum(!has_bse_match & outcome_agg == 2, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Episode × BSE family miss matrix
# One row per ERT episode overlapping the BSE period; one column per BSE family.
# Cells contain the BSE year(s) within the episode window (comma-separated if
# multiple variants in the family fire in different years), or NA if not captured.
# "ERT regime" is renamed ERT_regime (no space) for safe use in gt row conditions.
# ------------------------------------------------------------------------------

FAMILY_COLS <- c("ROW", "ERT_regime", "Polyarchy", "LibDem", "LibComp")

# Attach family metadata to every BSE event row.
bse_with_family <- all_bse_events |>
  left_join(
    definition_metadata |> select(bse_definition, definition_family),
    by = "bse_definition"
  ) |>
  select(country_text_id, bse_year = year, definition_family)

# For each (episode, family): collect unique BSE years that fall within the window.
# Episodes with no BSE match in any family produce a definition_family=NA row after
# the left_join. Filtering those out before pivot_wider avoids a spurious NA column,
# but it means those episodes disappear from the pivot result entirely — pivot_wider
# cannot implicitly add rows, only fill missing *columns*. We fix this by left-joining
# the pivot result back onto the full ert_eps_in_period baseline, which restores all
# unmatched episodes with NA across all five family columns.
bse_years_by_family <- ert_eps_in_period |>
  left_join(
    bse_with_family,
    by = join_by(
      country_text_id,
      ep_start_match <= bse_year,
      ep_end_match >= bse_year
    ),
    relationship = "many-to-many"
  ) |>
  group_by(
    aut_ep_id,
    country_text_id,
    ep_start,
    ep_end,
    outcome_agg,
    definition_family
  ) |>
  summarise(
    bse_years = {
      yrs <- sort(unique(bse_year[!is.na(bse_year)]))
      if (length(yrs) == 0) NA_character_ else paste(yrs, collapse = ",")
    },
    .groups = "drop"
  ) |>
  filter(!is.na(definition_family)) |>
  pivot_wider(names_from = definition_family, values_from = bse_years) |>
  rename(ERT_regime = `ERT regime`)

ert_miss_matrix <- ert_eps_in_period |>
  select(aut_ep_id, country_text_id, ep_start, ep_end, outcome_agg) |>
  left_join(
    bse_years_by_family,
    by = c("aut_ep_id", "country_text_id", "ep_start", "ep_end", "outcome_agg")
  ) |>
  left_join(aut_eps |> select(aut_ep_id, country_name), by = "aut_ep_id") |>
  mutate(
    outcome_short = case_when(
      outcome_agg == 1 ~ "Breakdown",
      outcome_agg == 2 ~ "Erosion only",
      outcome_agg == 4 ~ "Ongoing"
    )
  ) |>
  arrange(ep_start) |>
  select(country_name, ep_start, ep_end, outcome_short, all_of(FAMILY_COLS))

summary_table <- summary_table |>
  left_join(ert_coverage_summary, by = "bse_definition") |>
  left_join(definition_metadata, by = "bse_definition") |>
  mutate(
    n_ert_episodes = n_ert_ep_in_period,
    n_bse_episodes = n_breakdown_events,
    n_overlap = n_ert_ep_captured, # ERT episodes with >= 1 BSE event in window
    pct_bse_in_ert = round(100 * n_matched_ert_ep / n_breakdown_events, 1),
    pct_bse_not_in_ert = round(100 * n_no_ert_match / n_breakdown_events, 1),
    pct_ert_in_bse = round(100 * n_ert_ep_captured / n_ert_ep_in_period, 1),
    pct_ert_not_in_bse = round(100 * n_ert_ep_missed / n_ert_ep_in_period, 1),
    definition_family = factor(definition_family, levels = FAMILY_ORDER)
  ) |>
  arrange(definition_family, strictness) |>
  select(
    definition_family,
    strictness,
    bse_definition,
    threshold_display,
    n_ert_episodes,
    n_bse_episodes,
    n_overlap,
    pct_bse_in_ert,
    pct_bse_not_in_ert,
    pct_ert_in_bse,
    pct_ert_not_in_bse,
    # detail columns kept for CSV
    n_countries,
    n_matched_ert_ep,
    n_no_ert_match,
    n_ert_ep_missed,
    n_ert_missed_breakdown,
    n_ert_missed_erosion,
    n_ep_row_breakdown,
    n_ep_erosion_only,
    n_ep_ongoing
  )

# ------------------------------------------------------------------------------
# Outcome==2 episodes: which BSE definitions (if any) capture them?
#
# outcome_agg==2 means the episode did not cross the RoW threshold, so the
# ROW and ERT regime definitions cannot capture it by construction. But the
# threshold-based definitions (Polyarchy, LibDem, LibComp) use different
# cutoffs and CAN match outcome_agg==2 episodes if the index dips below their
# threshold without crossing RoW.
#
# We check each outcome_agg==2 episode for temporal overlap with BSE treatment
# events (any definition) using join_by(between(...)), the same approach as
# match_one_definition above. Truly novel episodes are those with no BSE match
# under any definition.
# ------------------------------------------------------------------------------

ert_erosion_eps <- aut_eps |>
  filter(outcome_agg == 2, ep_start <= BSE_END, ep_end >= BSE_START) |>
  mutate(
    ep_start_match = ep_start - EP_MATCH_TOL,
    ep_end_match = ep_end + EP_MATCH_TOL
  )

novel_episodes <- ert_erosion_eps |>
  left_join(
    all_bse_events |> select(country_text_id, bse_year = year, bse_definition),
    # ep_start_match/ep_end_match are on the left (ert_erosion_eps);
    # bse_year is on the right (all_bse_events).
    # join_by(x >= y) means LEFT.x >= RIGHT.y, so the two conditions together
    # select BSE events whose year falls within the episode window.
    by = join_by(
      country_text_id,
      ep_start_match <= bse_year,
      ep_end_match >= bse_year
    )
  ) |>
  group_by(
    aut_ep_id,
    country_text_id,
    country_name,
    ep_start,
    ep_end,
    poly_start,
    poly_end
  ) |>
  summarise(
    n_bse_definitions_capturing = n_distinct(bse_definition[
      !is.na(bse_definition)
    ]),
    bse_definitions_capturing = paste(
      sort(unique(bse_definition[!is.na(bse_definition)])),
      collapse = "; "
    ),
    .groups = "drop"
  ) |>
  mutate(
    captured_by_any_bse = n_bse_definitions_capturing > 0,
    bse_definitions_capturing = if_else(
      captured_by_any_bse,
      bse_definitions_capturing,
      "None"
    )
  )

# ------------------------------------------------------------------------------
# Timing gap: years from ERT episode start to BSE breakdown year
# For episodes matched to outcome==1 (breakdown) only.
# ------------------------------------------------------------------------------

timing_gap <- bse_matched |>
  filter(
    bse_definition %in% c("ROW", "ERT regime"),
    outcome_agg == 1,
    !is.na(ep_start)
  ) |>
  mutate(gap_years = year - ep_start) |>
  select(
    bse_definition,
    country_name,
    ep_start,
    breakdown_year = year,
    ep_end,
    gap_years
  )

# ------------------------------------------------------------------------------
# Print results
# ------------------------------------------------------------------------------

cat(
  "\n=== Our ERT autocratization episodes (aut_ep_prch == 1, all years) ===\n"
)
aut_eps |>
  count(outcome_agg, outcome_label) |>
  arrange(outcome_agg) |>
  print()

cat("\n\n=== BSE <-> ERT overlap summary (1999-2023) ===\n")
cat(sprintf(
  "Baseline: %d ERT autocratization episodes (aut_ep_prch==1) overlapping the BSE period.\n",
  nrow(ert_eps_in_period)
))
cat(
  "n_ert_episodes     — ERT episodes overlapping 1999-2023 (same baseline for all rows)\n"
)
cat("n_bse_episodes     — BSE breakdown events identified by this definition\n")
cat(
  "n_overlap          — ERT episodes with >= 1 BSE event inside their window (+/- EP_MATCH_TOL)\n"
)
cat(
  "pct_bse_in_ert     — share of BSE events that land inside an ERT episode (%)\n"
)
cat("pct_bse_not_in_ert — share of BSE events with no ERT episode match (%)\n")
cat(
  "pct_ert_in_bse     — share of ERT episodes captured by this BSE definition (%)\n"
)
cat(
  "pct_ert_not_in_bse — share of ERT episodes invisible to this definition (%)\n\n"
)

summary_table |>
  select(
    definition_family,
    bse_definition,
    n_ert_episodes,
    n_bse_episodes,
    n_overlap,
    pct_bse_in_ert,
    pct_bse_not_in_ert,
    pct_ert_in_bse,
    pct_ert_not_in_bse
  ) |>
  print(n = Inf)

cat("\n\n=== ERT outcome==2 episodes overlapping 1999-2023: BSE coverage ===\n")
cat("(outcome_agg==2: erosion that never crossed the RoW threshold)\n")
cat("ROW and ERT regime definitions cannot capture these by construction.\n")
cat(
  "Threshold-based definitions (Polyarchy, LibDem, LibComp) may capture them\n"
)
cat("if the index dips below their cutoff during the episode window.\n\n")

cat(sprintf(
  "Total outcome==2 episodes overlapping 1999-2023: %d (%d countries)\n",
  nrow(novel_episodes),
  n_distinct(novel_episodes$country_text_id)
))
cat(sprintf(
  "Captured by at least one BSE threshold definition: %d\n",
  sum(novel_episodes$captured_by_any_bse)
))
cat(sprintf(
  "NOT captured by any BSE definition (truly novel): %d\n\n",
  sum(!novel_episodes$captured_by_any_bse)
))

novel_episodes |>
  arrange(captured_by_any_bse, ep_start) |>
  select(
    country_name,
    ep_start,
    ep_end,
    poly_start,
    poly_end,
    captured_by_any_bse,
    bse_definitions_capturing
  ) |>
  print(n = Inf)

cat(
  "\n\n=== Timing gap: years from ERT episode start to BSE breakdown year ===\n"
)
cat("(ROW and ERT regime binary treatments, outcome==1 episodes only)\n\n")

timing_gap |>
  group_by(bse_definition) |>
  summarise(
    n = n(),
    mean_gap = round(mean(gap_years), 1),
    median_gap = median(gap_years),
    min_gap = min(gap_years),
    max_gap = max(gap_years),
    pct_gap_gt5 = round(100 * mean(gap_years > 5), 1),
    .groups = "drop"
  ) |>
  print()

cat("\nCountry-level detail (ROW):\n")
timing_gap |>
  filter(bse_definition == "ROW") |>
  arrange(desc(gap_years)) |>
  print(n = Inf)

# ------------------------------------------------------------------------------
# Publication-ready gt table
# Rows grouped by definition_family; threshold_display as the row identifier.
# Two column spanners: precision (BSE->ERT) and recall (ERT->BSE).
# ------------------------------------------------------------------------------

# Shared table styling: smaller font and clean borders for browser copy-paste.
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
      table_body.hlines.color = "#cccccc",
      row_group.border.top.style = "solid",
      row_group.border.top.width = px(1),
      row_group.border.top.color = "#666666",
      row_group.border.bottom.style = "solid",
      row_group.border.bottom.width = px(0.5),
      row_group.border.bottom.color = "#666666"
    )
}

gt_overlap <- summary_table |>
  select(
    definition_family,
    threshold_display,
    n_ert_episodes,
    n_bse_episodes,
    n_overlap,
    pct_bse_in_ert,
    pct_bse_not_in_ert,
    pct_ert_in_bse,
    pct_ert_not_in_bse
  ) |>
  gt(groupname_col = "definition_family") |>
  tab_header(
    title = "Overlap between ERT autocratization episodes and BSE treatment definitions",
    subtitle = "1999–2023 | ERT: aut_ep_prch == 1 | matching window ± 1 year"
  ) |>
  tab_spanner(
    label = "BSE → ERT",
    columns = c(pct_bse_in_ert, pct_bse_not_in_ert)
  ) |>
  tab_spanner(
    label = "ERT → BSE",
    columns = c(pct_ert_in_bse, pct_ert_not_in_bse)
  ) |>
  cols_label(
    threshold_display = "Definition",
    n_ert_episodes = "ERT episodes",
    n_bse_episodes = "BSE events",
    n_overlap = "Overlap",
    pct_bse_in_ert = "In ERT (%)",
    pct_bse_not_in_ert = "Not in ERT (%)",
    pct_ert_in_bse = "In BSE (%)",
    pct_ert_not_in_bse = "Not in BSE (%)"
  ) |>
  fmt_number(
    columns = c(n_ert_episodes, n_bse_episodes, n_overlap),
    decimals = 0
  ) |>
  fmt_percent(
    columns = c(
      pct_bse_in_ert,
      pct_bse_not_in_ert,
      pct_ert_in_bse,
      pct_ert_not_in_bse
    ),
    scale_values = FALSE,
    decimals = 1
  ) |>
  cols_align(
    align = "center",
    columns = c(
      n_ert_episodes,
      n_bse_episodes,
      n_overlap,
      pct_bse_in_ert,
      pct_bse_not_in_ert,
      pct_ert_in_bse,
      pct_ert_not_in_bse
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) |>
  tab_footnote(
    footnote = "Overlap = number of ERT episodes with at least one BSE breakdown event within their window.",
    locations = cells_column_labels(columns = n_overlap)
  ) |>
  opt_row_striping() |>
  apply_table_style()

# ------------------------------------------------------------------------------
# Publication-ready gt table: episode × BSE family miss matrix
# ------------------------------------------------------------------------------

gt_miss_matrix <- ert_miss_matrix |>
  gt() |>
  tab_header(
    title = "ERT autocratization episodes overlapping 1999–2023: coverage by BSE definition family",
    subtitle = "Cells show the year(s) BSE identifies a breakdown within the episode window; — = not captured"
  ) |>
  tab_spanner(
    label = "BSE breakdown year(s) within episode window",
    columns = all_of(FAMILY_COLS)
  ) |>
  cols_label(
    country_name = "Country",
    ep_start = "Start",
    ep_end = "End",
    outcome_short = "ERT outcome",
    ROW = "ROW",
    ERT_regime = "ERT regime",
    Polyarchy = "Polyarchy",
    LibDem = "LibDem",
    LibComp = "LibComp"
  ) |>
  sub_missing(columns = all_of(FAMILY_COLS), missing_text = "—") |>
  cols_align(align = "center", columns = all_of(FAMILY_COLS)) |>
  cols_width(
    country_name ~ px(110),
    ep_start ~ px(42),
    ep_end ~ px(42),
    outcome_short ~ px(80),
    ROW ~ px(52),
    ERT_regime ~ px(62),
    Polyarchy ~ px(62),
    LibDem ~ px(62),
    LibComp ~ px(62)
  ) |>
  tab_footnote(
    footnote = "A family captures an episode if any of its threshold variants has a breakdown event within the episode window (± 1 year). Multiple years shown when variants within a family fire at different times.",
    locations = cells_column_spanners()
  ) |>
  opt_row_striping() |>
  apply_table_style()

# Color cells: green for captured (non-NA), light gray for missed (NA).
gt_miss_matrix <- reduce(
  FAMILY_COLS,
  function(tbl, col) {
    sym_col <- sym(col)
    tbl |>
      tab_style(
        style = list(cell_fill(color = "#c7e9c0"), cell_text(weight = "bold")),
        locations = cells_body(columns = all_of(col), rows = !is.na(!!sym_col))
      ) |>
      tab_style(
        style = cell_fill(color = "#eeeeee"),
        locations = cells_body(columns = all_of(col), rows = is.na(!!sym_col))
      )
  },
  .init = gt_miss_matrix
)

# Vertical divider between episode info columns and BSE family columns.
# Smaller font for the year cells so multi-year entries don't blow out column width.
gt_miss_matrix <- gt_miss_matrix |>
  tab_style(
    style = cell_borders(sides = "left", color = "#666666", weight = px(1.5)),
    locations = cells_body(columns = ROW)
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "#666666", weight = px(1.5)),
    locations = cells_column_labels(columns = ROW)
  ) |>
  tab_style(
    style = cell_text(size = px(9)),
    locations = cells_body(columns = all_of(FAMILY_COLS))
  )

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

write_csv(summary_table, file.path(out_dir, "bse_overlap_summary.csv"))
write_csv(ert_miss_matrix, file.path(out_dir, "bse_ert_miss_matrix.csv"))
write_csv(novel_episodes, file.path(out_dir, "bse_novel_episodes.csv"))
write_csv(timing_gap, file.path(out_dir, "bse_timing_gap.csv"))
gtsave(gt_overlap, file.path(out_dir, "bse_overlap_table.html"))
gtsave(gt_miss_matrix, file.path(out_dir, "bse_ert_miss_matrix.html"))

message("\nSaved outputs to ", out_dir)
