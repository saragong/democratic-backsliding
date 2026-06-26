# ==============================================================================
# Build event study panel
# One row per (episode, event_time), t = -15 to +15 relative to episode start.
# Outcomes are kept in levels and as changes relative to the t = -1 baseline.
# ==============================================================================

library(tidyverse)
library(here)

data_dir <- here::here("data")

# --- Toggle -------------------------------------------------------------------
# FALSE (a): keep all event-window observations even when windows overlap
#            across episodes in the same country
# TRUE  (b): truncate each episode's window symmetrically:
#              forward window clipped at the start of the next episode
#              backward window clipped at the end of the previous episode
#            "next/previous episode" = any episode (dem or aut) in the same country
TRIM_AT_NEXT_EPISODE <- TRUE

# --- Config -------------------------------------------------------------------
WINDOW <- -15L:15L
OUTCOMES <- c(
  "gdp_pc_growth",
  "cpi_inflation",
  "unemployment_rate",
  "trade_pct_gdp",
  "top10_share",
  "gini_disp"
)

# --- Load panel ---------------------------------------------------------------
panel <- readRDS(file.path(data_dir, "combined_panel.rds"))

outcomes_panel <- panel |>
  select(country_text_id, year, all_of(OUTCOMES))

# --- Build episode list -------------------------------------------------------
dem_eps <- panel |>
  filter(!is.na(dem_ep_id)) |>
  distinct(
    country_text_id,
    country_name,
    ep_id         = dem_ep_id,
    ep_start_year = dem_ep_start_year,
    ep_end_year   = dem_ep_end_year
  ) |>
  mutate(ep_type = "Democratization")

aut_eps <- panel |>
  filter(!is.na(aut_ep_id)) |>
  distinct(
    country_text_id,
    country_name,
    ep_id         = aut_ep_id,
    ep_start_year = aut_ep_start_year,
    ep_end_year   = aut_ep_end_year
  ) |>
  mutate(ep_type = "Autocratization")

episodes <- bind_rows(dem_eps, aut_eps) |>
  arrange(country_text_id, ep_start_year)

message(sprintf(
  "Episodes: %d democratization, %d autocratization",
  nrow(dem_eps),
  nrow(aut_eps)
))

# --- Event grid ---------------------------------------------------------------
event_grid <- episodes |>
  mutate(event_time = list(WINDOW)) |>
  unnest(event_time) |>
  mutate(year = ep_start_year + event_time)

# --- Option (b): truncate window at episode boundaries -----------------------
# Forward: clip at the start year of the next episode in the same country.
# Backward: clip at the end year of the previous episode in the same country.
# Both use all episodes (dem + aut) sorted by start year within country.
if (TRIM_AT_NEXT_EPISODE) {
  boundaries <- episodes |>
    group_by(country_text_id) |>
    arrange(ep_start_year) |>
    mutate(
      next_ep_start = lead(ep_start_year),
      prev_ep_end   = lag(ep_end_year)
    ) |>
    ungroup() |>
    select(country_text_id, ep_id, ep_type, next_ep_start, prev_ep_end)

  event_grid <- event_grid |>
    left_join(boundaries, by = c("country_text_id", "ep_id", "ep_type")) |>
    filter(is.na(next_ep_start) | year < next_ep_start) |>
    filter(is.na(prev_ep_end)   | year > prev_ep_end) |>
    select(-next_ep_start, -prev_ep_end)
}

# --- Join outcomes ------------------------------------------------------------
event_study <- event_grid |>
  left_join(outcomes_panel, by = c("country_text_id", "year"))

# --- Baseline demean (t = -1) -------------------------------------------------
baseline <- event_study |>
  filter(event_time == -1L) |>
  select(ep_id, ep_type, all_of(OUTCOMES)) |>
  rename_with(~ paste0(.x, "_base"), all_of(OUTCOMES))

event_study <- event_study |>
  left_join(baseline, by = c("ep_id", "ep_type"))

for (v in OUTCOMES) {
  event_study[[paste0(v, "_chg")]] <- event_study[[v]] -
    event_study[[paste0(v, "_base")]]
}

event_study <- event_study |>
  select(-ends_with("_base")) |>
  relocate(
    ep_type,
    ep_id,
    country_text_id,
    country_name,
    ep_start_year,
    event_time,
    year
  )

# --- Coverage summary ---------------------------------------------------------
message(sprintf(
  "Event study panel: %d rows, %d episodes, window [%d, +%d]",
  nrow(event_study),
  n_distinct(paste(event_study$ep_type, event_study$ep_id)),
  min(WINDOW),
  max(WINDOW)
))

event_study |>
  group_by(ep_type, event_time) |>
  summarise(
    n_episodes = n(),
    across(
      paste0(OUTCOMES, "_chg"),
      ~ mean(!is.na(.x)),
      .names = "pct_nonmissing_{.col}"
    ),
    .groups = "drop"
  ) |>
  filter(event_time %in% c(-15L, -5L, 0L, 5L, 10L, 15L)) |>
  print(width = 120)

# --- Save ---------------------------------------------------------------------
saveRDS(event_study, file.path(data_dir, "event_study_panel.rds"))
message("Saved event_study_panel.rds to ", data_dir)
