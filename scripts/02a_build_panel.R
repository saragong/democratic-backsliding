# ==============================================================================
# Build combined country-year panel
# ==============================================================================
# Joins all data sources into a single country-year panel keyed on ERT's
# country_text_id (which is standard ISO-3) and year.
#
# Country code notes:
#   ERT  country_text_id — standard ISO-3 (COD, ROU, SGP, TWN, ARE, ...)
#   DDCG wbcode          — non-standard for 5 countries (ZAR, ROM, SIN, TAW, UAE)
#   Outcomes iso3c       — standard ISO-3  → join directly on country_text_id
#   NED  country_abb     — COW codes for ~105 countries; ned_panel.rds carries
#                          country_iso3 (ERT-compatible) derived in 01c
#
# Key derived variable:
#   executive_entry_year — the year the current ruling party first took power
#     (most recent NED turnover of any election type at or before that year).
#     For each ERT episode, the value at the episode's start year identifies
#     when the leader who presided over the episode came to power.
#
# Output: combined_panel.rds
# ==============================================================================

library(tidyverse)
library(here)

source(here::here("scripts", "vdem_indices.R"))

data_dir <- here::here("data")

# --- Load all sources ---------------------------------------------------------

ert <- readRDS(file.path(data_dir, "ert_episodes.rds"))
ddcg <- readRDS(file.path(data_dir, "ddcg_panel.rds"))
pwt_gdp <- readRDS(file.path(data_dir, "pwt_gdp.rds"))
imf_cpi <- readRDS(file.path(data_dir, "imf_cpi.rds"))
ilo_unemp <- readRDS(file.path(data_dir, "ilo_unemployment.rds"))
wb_trade <- readRDS(file.path(data_dir, "wb_trade.rds"))
wid_ineq <- readRDS(file.path(data_dir, "wid_inequality.rds"))
swiid <- readRDS(file.path(data_dir, "swiid_gini.rds"))
ned_panel <- readRDS(file.path(data_dir, "ned_panel.rds"))
archigos_panel <- readRDS(file.path(data_dir, "archigos_panel.rds"))
et_outcomes <- readRDS(file.path(data_dir, "et_outcomes.rds"))
vdem_sub <- readRDS(file.path(data_dir, "vdem_subcomponents.rds"))

# --- Crosswalk (ERT → DDCG only; not needed for outcomes or NED) -------------

ert_to_ddcg <- c(
  COD = "ZAR",
  ROU = "ROM",
  SGP = "SIN",
  TWN = "TAW",
  ARE = "UAE"
)

# --- NED: derive last_turnover_year per country × election_type × year -------
# For each row, carry forward the year of the most recent turnover so far.
# This lets us look up "who came to power most recently" at any given year.

ned_enriched <- ned_panel |>
  group_by(country_abb, election_type) |>
  mutate(
    last_turnover_year = if_else(turnover == 1L, year, NA_integer_),
    last_turnover_party = if_else(turnover == 1L, winner_party, NA_character_)
  ) |>
  fill(last_turnover_year, last_turnover_party, .direction = "down") |>
  ungroup()

ned_pres <- ned_enriched |>
  filter(election_type == "presidential") |>
  slice_tail(n = 1, by = c(country_abb, year)) |> # keep last if two elections in same year
  select(
    country_iso3,
    year,
    ruling_party_pres = winner_party,
    ruling_candidate_pres = winner_candidate,
    election_year_pres = election_year,
    turnover_pres = turnover,
    last_turnover_year_pres = last_turnover_year
  )

ned_parl <- ned_enriched |>
  filter(election_type == "parliamentary") |>
  slice_tail(n = 1, by = c(country_abb, year)) |> # keep last if two elections in same year
  select(
    country_iso3,
    year,
    ruling_party_parl = winner_party,
    election_year_parl = election_year,
    turnover_parl = turnover,
    last_turnover_year_parl = last_turnover_year
  )

# --- ERT backbone -------------------------------------------------------------

ert_core <- ert |>
  select(
    country_text_id,
    country_name,
    year,
    v2x_polyarchy,
    reg_type,
    dem_ep,
    dem_ep_id,
    dem_ep_start_year,
    dem_ep_end_year,
    dem_ep_outcome_agg,
    dem_ep_prch,
    aut_ep,
    aut_ep_id,
    aut_ep_start_year,
    aut_ep_end_year,
    aut_ep_outcome_agg,
    aut_ep_prch
  )

# --- Join DDCG ----------------------------------------------------------------
# Recode ERT codes to DDCG wbcodes for the 5 non-standard countries, join,
# then drop the temporary wbcode column.

ddcg_core <- ddcg |>
  select(wbcode, year, dem, demevent, revevent) |>
  mutate(dem = as.numeric(dem))

panel <- ert_core |>
  mutate(
    wbcode = recode(country_text_id, !!!ert_to_ddcg, .default = country_text_id)
  ) |>
  left_join(ddcg_core, by = c("wbcode", "year")) |>
  select(-wbcode)

# --- Join economic outcomes (iso3c == country_text_id for standard ISO-3) ----

panel <- panel |>
  left_join(
    pwt_gdp |> select(iso3c, year, gdp_pc_growth, ln_gdp_pc),
    by = c("country_text_id" = "iso3c", "year")
  ) |>
  left_join(
    imf_cpi |> select(iso3c, year, cpi_inflation),
    by = c("country_text_id" = "iso3c", "year")
  ) |>
  left_join(
    ilo_unemp |> select(iso3c, year, unemployment_rate),
    by = c("country_text_id" = "iso3c", "year")
  ) |>
  left_join(
    wb_trade |> select(iso3c, year, trade_pct_gdp),
    by = c("country_text_id" = "iso3c", "year")
  ) |>
  left_join(
    wid_ineq |> select(iso3c, year, top10_share),
    by = c("country_text_id" = "iso3c", "year")
  ) |>
  left_join(
    swiid |>
      select(iso3c, year, gini_disp, gini_mkt) |>
      summarise(across(c(gini_disp, gini_mkt), mean, na.rm = TRUE),
                .by = c(iso3c, year)),
    by = c("country_text_id" = "iso3c", "year")
  )

# --- Join the additional outcomes built in 01f --------------------------------
# Both are already keyed on country_text_id/year (et_outcomes resolves
# outcomes.dta's labelled ISO3 `id`; vdem_subcomponents comes straight from
# V-Dem, which shares this key natively), so neither needs a crosswalk.
#
# ln_gdp_pc_wb / ln_gdp_pc_imf are alternative log GDP-per-capita series that
# difference exactly like PWT's ln_gdp_pc above, so the RDD growth outcome can
# be shown from three sources on one panel.

panel <- panel |>
  left_join(
    et_outcomes |>
      select(
        country_text_id,
        year,
        ln_gdp_pc_wb,
        ln_gdp_pc_imf,
        debt_pct_gdp,
        deficit_pct_gdp
      ),
    by = c("country_text_id", "year")
  ) |>
  left_join(
    vdem_sub |>
      select(
        country_text_id,
        year,
        v2x_jucon,
        v2xlg_legcon,
        checks_balances,
        hos_power_linear,
        hog_power_linear,
        hos_power_vdem,
        hog_power_vdem,
        # The high- and mid-level democracy indices, listed in
        # scripts/vdem_indices.R. all_of() so a name that 01f stopped
        # producing is an error here rather than a silently absent outcome.
        all_of(VDEM_INDEX_NEW_VARS)
      ),
    by = c("country_text_id", "year")
  )

# --- Join NED (on country_iso3, which maps COW codes to ERT country_text_id) --

panel <- panel |>
  left_join(ned_pres, by = c("country_text_id" = "country_iso3", "year")) |>
  left_join(ned_parl, by = c("country_text_id" = "country_iso3", "year"))

# --- Join Archigos (leader name + entry year, primary source 1875-2015) -------
# Archigos has explicit end-dates for every tenure so it avoids the NED
# forward-fill problem (e.g. Bonaparte 1848 filling all the way to 1965 in FR).
# For years after 2015 (Archigos cutoff) we fall back to NED.

panel <- panel |>
  left_join(archigos_panel, by = c("country_text_id" = "iso3c", "year"))

# --- Executive entry year -----------------------------------------------------
# Archigos start year for the leader currently in office.

panel <- panel |>
  mutate(
    executive_entry_year = archigos_start_year
  )

# --- Summary ------------------------------------------------------------------

n_ep_with_entry <- panel |>
  filter(dem_ep == 1 | aut_ep == 1) |>
  filter(!is.na(executive_entry_year)) |>
  distinct(
    country_text_id,
    ep_id = coalesce(as.character(dem_ep_id), as.character(aut_ep_id))
  ) |>
  nrow()

message(sprintf(
  "Combined panel: %d country-years, %d countries, %d-%d",
  nrow(panel),
  n_distinct(panel$country_text_id),
  min(panel$year),
  max(panel$year)
))
message(sprintf(
  "Episode-years with executive_entry_year: %d of %d (%.0f%%)",
  sum(
    !is.na(panel$executive_entry_year) &
      (panel$dem_ep == 1 | panel$aut_ep == 1),
    na.rm = TRUE
  ),
  sum(panel$dem_ep == 1 | panel$aut_ep == 1, na.rm = TRUE),
  100 *
    mean(
      !is.na(panel$executive_entry_year[panel$dem_ep == 1 | panel$aut_ep == 1]),
      na.rm = TRUE
    )
))

# --- Save ---------------------------------------------------------------------

saveRDS(panel, file.path(data_dir, "combined_panel.rds"))
message("Saved combined_panel.rds to ", data_dir)
