# ==============================================================================
# DDCG (Acemoglu et al. 2019) overlap analysis
# How well does DDCG's binary revevent flag capture our ERT autocratization
# episodes? This is the central motivating gap for the project: DDCG uses a
# single-year, dichotomous democracy indicator, so it should miss gradual
# erosion episodes that ERT records as autocratization but that never cross
# DDCG's binary threshold.
#
# DDCG sample:  184 countries, 1960-2010 (wbcode identifier).
# ERT sample:   183 countries, 1900-2025 (country_text_id identifier).
# Overlap analysis restricted to ERT episodes overlapping DDCG's 1960-2010
# coverage window. Unlike the Funke comparison, DDCG's country coverage is
# not a meaningful subset of ERT's, so no country-level restriction is applied.
#
# ERT episode base set: all aut_ep==1 episodes (not just aut_ep_prch==1),
# same scope as 08_funke_overlap.R. The miss table additionally reports
# whether each episode originated in a democracy (aut_ep_prch), since DDCG's
# dichotomized measure can only ever record a reversal for episodes that
# started in a democracy in the first place.
#
# Outputs:
#   output/ddcg_overlap_summary.csv  — single-row overlap-metrics table
#   output/ddcg_ert_miss_table.csv   — episode-level accounting table
#   output/ddcg_unmatched_ert.csv    — ERT episodes with no DDCG revevent
#   output/ddcg_unmatched_events.csv — DDCG revevent years with no ERT episode
#   output/ddcg_timing.csv           — timing gap for matched pairs
#   output/ddcg_overlap_table.html, ddcg_ert_miss_table.html,
#   ddcg_unmatched_events.html — gt versions of the above
# ==============================================================================

library(tidyverse)
library(gt)
library(here)

data_dir <- here::here("data")
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE)

DDCG_START <- 1960
DDCG_END <- 2010
EP_MATCH_TOL <- 1 # years of slack matching a DDCG revevent year to an ERT episode window

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

ert <- readRDS(file.path(data_dir, "ert_episodes.rds"))
ddcg <- readRDS(file.path(data_dir, "ddcg_panel.rds"))

# Crosswalk: DDCG wbcode -> ERT country_text_id (inverse of the standard
# ert_to_ddcg crosswalk used in 01a/02a/02b). Serbia/Montenegro excluded:
# DDCG's SER = pre-2006 Serbia & Montenegro union, not 1-to-1 mappable to
# ERT's separate SRB/MNE.
ert_to_ddcg <- c(
  COD = "ZAR",
  ROU = "ROM",
  SGP = "SIN",
  TWN = "TAW",
  ARE = "UAE"
)
ddcg_to_ert <- setNames(names(ert_to_ddcg), unname(ert_to_ddcg))

# ERT autocratization episodes, all aut_ep==1 (no aut_ep_prch filter — same
# scope as 08_funke_overlap.R). prch is retained so the miss table can show
# whether each episode originated in a democracy.
aut_eps <- ert |>
  filter(aut_ep == 1) |>
  group_by(aut_ep_id) |>
  summarise(
    country_text_id = first(country_text_id),
    country_name = first(country_name),
    ep_start = first(aut_ep_start_year),
    ep_end = first(aut_ep_end_year),
    outcome_agg = first(aut_ep_outcome_agg),
    prch = first(aut_ep_prch),
    poly_start = v2x_polyarchy[year == first(aut_ep_start_year)][1],
    poly_end = v2x_polyarchy[year == first(aut_ep_end_year)][1],
    .groups = "drop"
  ) |>
  mutate(
    outcome_label = case_when(
      outcome_agg == 1 ~ "Breakdown",
      outcome_agg == 2 ~ "Erosion only",
      outcome_agg == 3 ~ "Autocracy deepening",
      outcome_agg == 4 ~ "Uncertain / ongoing"
    ),
    prch_label = if_else(prch == 1, "Yes", "No")
  )

# DDCG reversal (autocratization) events, recoded onto ERT's country codes.
# A couple of DDCG country_name values (Côte d'Ivoire, São Tomé & Príncipe)
# carry mis-encoded bytes tagged as UTF-8, which crashes gtsave()'s HTML
# rendering downstream; re-encode from latin1 to fix.
fix_encoding <- function(x) {
  if_else(validUTF8(x), x, iconv(x, from = "latin1", to = "UTF-8"))
}

ddcg_events <- ddcg |>
  filter(revevent == 1) |>
  transmute(
    country_text_id = recode(wbcode, !!!ddcg_to_ert, .default = wbcode),
    country_name = fix_encoding(country_name),
    year
  )

# ------------------------------------------------------------------------------
# Restrict to ERT episodes overlapping DDCG's coverage window
# ------------------------------------------------------------------------------

aut_eps_in_period <- aut_eps |>
  filter(ep_start <= DDCG_END, ep_end >= DDCG_START) |>
  mutate(
    ep_start_match = ep_start - EP_MATCH_TOL,
    ep_end_match = ep_end + EP_MATCH_TOL
  )

ddcg_events_in_period <- ddcg_events |>
  filter(year >= DDCG_START, year <= DDCG_END)

# ------------------------------------------------------------------------------
# Match each DDCG revevent year to an ERT episode window (point event, same
# pattern as 07_bse_overlap.R's match_one_definition — a single "definition"
# here, so no group_modify by definition needed).
# ------------------------------------------------------------------------------

ddcg_matched <- ddcg_events_in_period |>
  left_join(
    aut_eps_in_period |>
      select(
        country_text_id,
        aut_ep_id,
        ep_start,
        ep_end,
        ep_start_match,
        ep_end_match,
        outcome_agg,
        outcome_label,
        prch,
        prch_label,
        poly_start,
        poly_end
      ),
    by = join_by(country_text_id, between(year, ep_start_match, ep_end_match))
  ) |>
  # between() join can still match one event to multiple overlapping episodes;
  # keep the one whose ep_end is closest to the revevent year.
  group_by(country_text_id, year) |>
  slice_min(abs(ep_end - year), n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(-ep_start_match, -ep_end_match)

# ------------------------------------------------------------------------------
# Episode-level miss table (ERT -> DDCG)
# left_join with join_by inequality silently drops left rows for episodes
# with zero DDCG matches — same fan-out issue fixed in 07/08. Fix: compute
# matches first, then left-join the result back to the full baseline.
# ------------------------------------------------------------------------------

ddcg_event_matches <- aut_eps_in_period |>
  left_join(
    ddcg_events_in_period |> select(country_text_id, ddcg_year = year),
    by = join_by(
      country_text_id,
      ep_start_match <= ddcg_year,
      ep_end_match >= ddcg_year
    )
  ) |>
  group_by(aut_ep_id) |>
  summarise(
    ddcg_event_years = {
      yrs <- sort(unique(ddcg_year[!is.na(ddcg_year)]))
      if (length(yrs) == 0) NA_character_ else paste(yrs, collapse = ", ")
    },
    .groups = "drop"
  )

ert_miss_table <- aut_eps_in_period |>
  select(aut_ep_id, country_name, ep_start, ep_end, outcome_agg, prch_label) |>
  left_join(ddcg_event_matches, by = "aut_ep_id") |>
  mutate(
    outcome_label = case_when(
      outcome_agg == 1 ~ "Democratic breakdown",
      outcome_agg == 2 ~ "Democratic erosion",
      outcome_agg == 3 ~ "Autocracy deepening",
      outcome_agg == 4 ~ "Ongoing"
    )
  ) |>
  arrange(ep_start) |>
  select(
    country_name,
    ep_start,
    ep_end,
    outcome_label,
    prch_label,
    ddcg_event_years
  )

# ------------------------------------------------------------------------------
# Summary table (one row, both directions)
# ------------------------------------------------------------------------------

n_ert_episodes <- nrow(aut_eps_in_period)
n_ddcg_events <- nrow(ddcg_events_in_period)
n_ddcg_matched <- sum(!is.na(ddcg_matched$aut_ep_id))
n_ert_matched <- sum(!is.na(ert_miss_table$ddcg_event_years))

pct_ddcg_in_ert <- round(100 * n_ddcg_matched / n_ddcg_events, 1)
pct_ert_in_ddcg <- round(100 * n_ert_matched / n_ert_episodes, 1)

summary_table <- tibble(
  n_ert_episodes = n_ert_episodes,
  n_ddcg_events = n_ddcg_events,
  n_overlap = n_ert_matched,
  pct_ddcg_in_ert = pct_ddcg_in_ert,
  pct_ddcg_not_in_ert = round(100 - pct_ddcg_in_ert, 1),
  pct_ert_in_ddcg = pct_ert_in_ddcg,
  pct_ert_not_in_ddcg = round(100 - pct_ert_in_ddcg, 1)
)

# ------------------------------------------------------------------------------
# Unmatched sets and timing gap
# ------------------------------------------------------------------------------

unmatched_ert <- ert_miss_table |> filter(is.na(ddcg_event_years))

unmatched_events <- ddcg_matched |>
  filter(is.na(aut_ep_id)) |>
  select(country_text_id, country_name, year)

timing <- ddcg_matched |>
  filter(!is.na(aut_ep_id)) |>
  mutate(gap_years = year - ep_start) |>
  select(
    country_name,
    ep_start,
    revevent_year = year,
    ep_end,
    outcome_label,
    prch_label,
    gap_years
  ) |>
  arrange(desc(gap_years))

# ------------------------------------------------------------------------------
# Print results
# ------------------------------------------------------------------------------

cat("=== Overview ===\n")
cat(sprintf(
  "DDCG revevent events (%d-%d):        %d\n",
  DDCG_START,
  DDCG_END,
  n_ddcg_events
))
cat(sprintf(
  "ERT aut_ep episodes overlapping %d-%d: %d\n",
  DDCG_START,
  DDCG_END,
  n_ert_episodes
))
cat(sprintf(
  "  of which originated in a democracy (prch==1): %d\n",
  sum(aut_eps_in_period$prch == 1, na.rm = TRUE)
))
cat(sprintf(
  "  of which already autocratic (prch==0):         %d\n",
  sum(aut_eps_in_period$prch == 0, na.rm = TRUE)
))

cat("\n\n=== DDCG <-> ERT overlap summary ===\n")
cat(sprintf(
  "n_overlap          — ERT episodes with >= 1 DDCG revevent inside their window (+/- %d yr)\n",
  EP_MATCH_TOL
))
cat(
  "pct_ddcg_in_ert    — share of DDCG revevent events that land inside an ERT episode (%)\n"
)
cat("pct_ert_in_ddcg    — share of ERT episodes captured by DDCG (%)\n\n")
print(summary_table)

cat("\n\n=== ERT episodes with NO matching DDCG revevent ===\n")
cat("(the core motivating gap: DDCG's binary flag misses these entirely)\n\n")
unmatched_ert |>
  count(outcome_label, prch_label) |>
  print()
cat("\n")
print(unmatched_ert, n = Inf)

cat("\n\n=== DDCG revevent events with NO matching ERT episode ===\n\n")
print(unmatched_events, n = Inf)

cat(
  "\n\n=== Timing gap: years from ERT episode start to DDCG revevent year ===\n"
)
timing |>
  summarise(
    n = n(),
    mean_gap = round(mean(gap_years), 1),
    median_gap = median(gap_years),
    min_gap = min(gap_years),
    max_gap = max(gap_years),
    pct_gap_gt5 = round(100 * mean(gap_years > 5), 1)
  ) |>
  print()
cat("\nCountry-level detail:\n")
print(timing, n = Inf)

# ------------------------------------------------------------------------------
# Publication-ready gt tables
# ------------------------------------------------------------------------------

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

# ---- gt summary table -------------------------------------------------------

gt_summary <- summary_table |>
  gt() |>
  tab_header(
    title = "ERT autocratization episodes vs. DDCG (Acemoglu et al. 2019) revevent",
    subtitle = paste0(
      "ERT episodes overlapping ",
      DDCG_START,
      "–",
      DDCG_END,
      "; matching window ± ",
      EP_MATCH_TOL,
      " year"
    )
  ) |>
  tab_spanner(
    label = "DDCG → ERT",
    columns = c(pct_ddcg_in_ert, pct_ddcg_not_in_ert)
  ) |>
  tab_spanner(
    label = "ERT → DDCG",
    columns = c(pct_ert_in_ddcg, pct_ert_not_in_ddcg)
  ) |>
  cols_label(
    n_ert_episodes = "ERT episodes",
    n_ddcg_events = "DDCG events",
    n_overlap = "Overlap",
    pct_ddcg_in_ert = "% in ERT",
    pct_ddcg_not_in_ert = "% not in ERT",
    pct_ert_in_ddcg = "% in DDCG",
    pct_ert_not_in_ddcg = "% not in DDCG"
  ) |>
  cols_align(align = "center", columns = everything()) |>
  cols_width(
    n_ert_episodes ~ px(85),
    n_ddcg_events ~ px(85),
    n_overlap ~ px(70),
    pct_ddcg_in_ert ~ px(75),
    pct_ddcg_not_in_ert ~ px(90),
    pct_ert_in_ddcg ~ px(75),
    pct_ert_not_in_ddcg ~ px(90)
  ) |>
  apply_table_style()

# ---- gt miss table -----------------------------------------------------------

gt_miss_table <- ert_miss_table |>
  gt() |>
  tab_header(
    title = paste0(
      "ERT autocratization episodes overlapping ",
      DDCG_START,
      "–",
      DDCG_END
    ),
    subtitle = "Cell shows DDCG revevent year(s) within the episode window; — = not captured"
  ) |>
  cols_label(
    country_name = "Country",
    ep_start = "Start",
    ep_end = "End",
    outcome_label = "ERT outcome",
    prch_label = "Originated in democracy",
    ddcg_event_years = "DDCG revevent year(s)"
  ) |>
  sub_missing(columns = ddcg_event_years, missing_text = "—") |>
  cols_align(align = "center", columns = c(prch_label, ddcg_event_years)) |>
  cols_width(
    country_name ~ px(120),
    ep_start ~ px(45),
    ep_end ~ px(45),
    outcome_label ~ px(90),
    prch_label ~ px(80),
    ddcg_event_years ~ px(95)
  ) |>
  tab_style(
    style = list(cell_fill(color = "#c7e9c0"), cell_text(weight = "bold")),
    locations = cells_body(
      columns = ddcg_event_years,
      rows = !is.na(ddcg_event_years)
    )
  ) |>
  tab_style(
    style = cell_fill(color = "#eeeeee"),
    locations = cells_body(
      columns = ddcg_event_years,
      rows = is.na(ddcg_event_years)
    )
  ) |>
  tab_style(
    style = cell_text(size = px(9)),
    locations = cells_body(columns = ddcg_event_years)
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "#666666", weight = px(1.5)),
    locations = cells_body(columns = ddcg_event_years)
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "#666666", weight = px(1.5)),
    locations = cells_column_labels(columns = ddcg_event_years)
  ) |>
  opt_row_striping() |>
  apply_table_style()

# ---- gt unmatched-events table ------------------------------------------------
# DDCG revevent years with no ERT autocratization episode within ± EP_MATCH_TOL.

gt_unmatched_events <- unmatched_events |>
  arrange(year) |>
  gt() |>
  tab_header(
    title = "DDCG revevent years with no matching ERT autocratization episode",
    subtitle = paste0(
      "DDCG events (",
      DDCG_START,
      "–",
      DDCG_END,
      ") with no ERT aut_ep episode within ± ",
      EP_MATCH_TOL,
      " year"
    )
  ) |>
  cols_label(
    country_text_id = "Code",
    country_name = "Country",
    year = "Revevent year"
  ) |>
  cols_align(align = "center", columns = c(country_text_id, year)) |>
  cols_width(
    country_text_id ~ px(60),
    country_name ~ px(140),
    year ~ px(90)
  ) |>
  opt_row_striping() |>
  apply_table_style()

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------

write_csv(summary_table, file.path(out_dir, "ddcg_overlap_summary.csv"))
write_csv(ert_miss_table, file.path(out_dir, "ddcg_ert_miss_table.csv"))
write_csv(unmatched_ert, file.path(out_dir, "ddcg_unmatched_ert.csv"))
write_csv(unmatched_events, file.path(out_dir, "ddcg_unmatched_events.csv"))
write_csv(timing, file.path(out_dir, "ddcg_timing.csv"))
gtsave(gt_summary, file.path(out_dir, "ddcg_overlap_table.html"))
gtsave(gt_miss_table, file.path(out_dir, "ddcg_ert_miss_table.html"))
gtsave(gt_unmatched_events, file.path(out_dir, "ddcg_unmatched_events.html"))

message("Saved outputs to ", out_dir)
