# ==============================================================================
# Load National Elections Database v2.0 (Marx, Pons & Rollet 2025)
# ==============================================================================
# 6,309 national elections, 1789–2023, 212 countries.
# Source: https://www.nationalelectionsdatabase.com/
# Codebook: nationalelectionsdatabase.com/codebook/codebook_elections_database_v2.pdf
#
# Presidential elections (1,409):
#   Candidates ranked 1–39 by electoral-college score, then second-round,
#   then first-round vote share. winner = candidate_1 / party_1.
#   winner_share = vote_share1_1 (first- or only-round share of winner).
#
# Parliamentary elections (4,900):
#   Parties ranked 1–74 by seat share. winner = party_1.
#   winner_share = seat_share_1 (note: codebook says seats_share_i but actual column is seat_share_i).
#
# country_abb: three-letter code, a mix of ISO-3 and COW alpha codes.
# country_cow: numeric COW code (the authoritative identifier in NED).
# country_iso3 is derived from country_cow via countrycode("cown","iso3c"),
#   with manual overrides for entities whose COW code is absent from the
#   countrycode package or that ERT encodes non-standardly (DDR, XKX, YMD, SML).
#
# Outputs:
#   ned_elections.rds  — one row per election, winner identified
#   ned_panel.rds      — country × election_type × year panel, winner
#                        forward-filled between elections
# ==============================================================================

library(tidyverse)
library(haven)
library(countrycode)
library(here)

data_dir <- here::here("data")

# --- Download -----------------------------------------------------------------

base_url <- "https://www.nationalelectionsdatabase.com/datasets"

for (lst in list(
  list(
    dest = file.path(data_dir, "ned_presidential_v2.dta"),
    file = "presidential_elections_v2.dta"
  ),
  list(
    dest = file.path(data_dir, "ned_parliamentary_v2.dta"),
    file = "parliamentary_elections_v2.dta"
  )
)) {
  if (!file.exists(lst$dest)) {
    message("Downloading ", lst$file, " ...")
    download.file(paste0(base_url, "/", lst$file), lst$dest, mode = "wb")
    message("Done.")
  } else {
    message(lst$file, " already present, skipping.")
  }
}

# --- Presidential elections ---------------------------------------------------

pres_raw <- read_dta(file.path(data_dir, "ned_presidential_v2.dta")) |>
  as_tibble() |>
  mutate(across(where(is.character), ~ na_if(as.character(.x), "")))

pres <- pres_raw |>
  filter(coalesce(flag_inconsequential, 0) == 0) |>
  transmute(
    country_abb = country_abb,
    country_cow = as.numeric(country_cow),
    year = as.integer(year),
    date = date,
    election_type = "presidential",
    winner_candidate = as.character(candidate_1),
    winner_party = as.character(party_1),
    winner_share = as.numeric(vote_share1_1),
    flag_coup = as.integer(coalesce(flag_coup, 0)),
    flag_indirect = as.integer(coalesce(flag_indirect, 0))
  ) |>
  arrange(country_abb, year)

message(sprintf(
  "Presidential: %d elections, %d countries",
  nrow(pres),
  n_distinct(pres$country_abb)
))

# --- Parliamentary elections --------------------------------------------------

parl_raw <- read_dta(file.path(data_dir, "ned_parliamentary_v2.dta")) |>
  as_tibble() |>
  mutate(across(where(is.character), ~ na_if(as.character(.x), "")))

parl <- parl_raw |>
  filter(
    coalesce(flag_inconsequential, 0) == 0,
    coalesce(flag_constituent, 0) == 0
  ) |>
  transmute(
    country_abb = country_abb,
    country_cow = as.numeric(country_cow),
    year = as.integer(year),
    date = date,
    election_type = "parliamentary",
    winner_candidate = NA_character_,
    winner_party = as.character(party_1),
    winner_share = as.numeric(seat_share_1),
    flag_coup = as.integer(coalesce(flag_coup, 0)),
    flag_indirect = 0L
  ) |>
  arrange(country_abb, year)

message(sprintf(
  "Parliamentary: %d elections, %d countries",
  nrow(parl),
  n_distinct(parl$country_abb)
))

# --- Combine and derive turnover ----------------------------------------------
# turnover = 1 when the winning party changes relative to the previous
# election of the same type in the same country.
# Note: party-name matching is imperfect (renames, coalitions); treat as
# approximate and supplement with manual review for key cases.

ned_elections <- bind_rows(pres, parl) |>
  arrange(country_abb, election_type, year) |>
  group_by(country_abb, election_type) |>
  mutate(
    turnover = as.integer(
      !is.na(winner_party) &
        !is.na(lag(winner_party)) &
        winner_party != lag(winner_party)
    ),
    election_year = 1L
  ) |>
  ungroup()

# --- Harmonise country codes to match ERT (country_text_id) -------------------
# Two-step: if country_abb is already a genuine ISO-3 code, use it directly.
# Otherwise convert via country_cow (numeric COW code) using countrycode.
# The exception: 5 country_abb values that ARE valid ISO-3 codes but NED uses
# them as COW alpha codes for different countries — force COW conversion for
# these (AUS=Austria, MAC=Macedonia, MNG=Montenegro, SLV=Slovenia, SWZ=Switzerland).
# Manual overrides handle 5 additional cases:
#   265 → DDR  (East Germany; ERT uses DDR)
#   347 → XKX  (Kosovo; ERT uses XKX)
#   680 → YMD  (South Yemen; ERT uses historical YMD)
#   SRB → SRB  (Serbia; COW 345 is shared with Yugoslavia)
#   SML → SML  (Somaliland; no COW code, ERT uses non-standard SML)
# Remaining NAs (Netherlands Antilles, Serbia & Montenegro, North Yemen,
# Yugoslavia) do not appear in ERT and are intentionally left unmatched.

NED_COWC_CONFLICTS <- c("AUS", "MAC", "MNG", "SLV", "SWZ")

ned_elections <- ned_elections |>
  mutate(
    country_iso3 = case_when(
      # country_abb is already the correct ISO-3 (and not a COW alpha collision)
      !is.na(countrycode(country_abb, "iso3c", "iso3c", warn = FALSE)) &
        !country_abb %in% NED_COWC_CONFLICTS ~
        countrycode(country_abb, "iso3c", "iso3c", warn = FALSE),
      # Use numeric COW code for all others
      TRUE ~ coalesce(
        countrycode(country_cow, "cown", "iso3c", warn = FALSE),
        case_when(
          country_cow == 265   ~ "DDR",
          country_cow == 347   ~ "XKX",
          country_cow == 680   ~ "YMD",
          country_abb == "SRB" ~ "SRB",
          country_abb == "SML" ~ "SML"
        )
      )
    )
  )

# --- Build country-year panel by forward-filling winner ----------------------
# Within each (country, election_type) group, expand to every year from first
# observed election through 2023, then carry winner forward.

ned_panel <- ned_elections |>
  group_by(country_abb, election_type) |>
  complete(year = min(year):2023L) |>
  fill(
    winner_party,
    winner_candidate,
    winner_share,
    flag_coup,
    flag_indirect,
    country_cow,
    country_iso3,
    .direction = "down"
  ) |>
  mutate(
    election_year = replace_na(election_year, 0L),
    turnover = replace_na(turnover, 0L)
  ) |>
  ungroup()

message(sprintf(
  "NED panel: %d country-election_type-years, %d countries",
  nrow(ned_panel),
  n_distinct(ned_panel$country_abb)
))

# --- Save ---------------------------------------------------------------------

saveRDS(ned_elections, file.path(data_dir, "ned_elections.rds"))
saveRDS(ned_panel, file.path(data_dir, "ned_panel.rds"))
message("Saved ned_elections.rds and ned_panel.rds to ", data_dir)
