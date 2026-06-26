# ==============================================================================
# Load Archigos v4.1 — Database of Political Leaders
# Goemans, Gleditsch & Chiozza (2009); data through December 31, 2015.
# Source: https://www.rochester.edu/college/faculty/hgoemans/data.htm
# ==============================================================================
# 3,409 leader-spells, ~188 countries, 1875-2015.
# Unit of analysis: leader-spell (one row per continuous tenure).
#
# Key variables:
#   ccode:     COW numeric country code
#   leader:    leader surname (or known name)
#   startdate: date entered office (character, "YYYY-MM-DD")
#   enddate:   date left office  (character, "YYYY-MM-DD")
#   entry:     "Regular" | "Irregular" | "Foreign Intervention"
#   exit:      "Regular" | "Irregular" | "Natural Death" | "Ill Health" |
#              "Still in Office" | "Other"
#
# Output: archigos_panel.rds
#   One row per (iso3c, year): the leader who held power at year-end
#   (latest entrant when two leaders shared a calendar year).
#   Also carries archigos_start_year for episode-shading in plot_helpers.R.
# ==============================================================================

library(tidyverse)
library(haven)
library(countrycode)
library(here)

data_dir <- here::here("data")

# --- Download -----------------------------------------------------------------

arch_path <- file.path(data_dir, "Archigos_4.1_stata14.dta")

if (!file.exists(arch_path)) {
  message("Downloading Archigos v4.1 ...")
  download.file(
    "https://www.rochester.edu/college/faculty/hgoemans/Archigos_4.1_stata14.dta",
    arch_path, mode = "wb"
  )
  message("Download complete.")
} else {
  message("Archigos already present, skipping download.")
}

# --- Read and clean -----------------------------------------------------------

archigos <- read_dta(arch_path) |>
  as_tibble() |>
  mutate(across(where(is.character), ~ na_if(as.character(.x), ""))) |>
  # Stata .dta files often embed Latin-1/Windows-1252 bytes in string fields.
  # Convert to valid UTF-8 so str_trunc / stringi never see invalid sequences.
  mutate(across(where(is.character), ~ iconv(.x, from = "latin1", to = "UTF-8", sub = ""))) |>
  transmute(
    ccode      = as.integer(ccode),
    leader     = as.character(leader),
    startdate  = as.Date(startdate),    # stored as "YYYY-MM-DD" character
    enddate    = as.Date(enddate),
    entry      = as.character(entry),
    exit       = as.character(exit)
  ) |>
  filter(!is.na(ccode), !is.na(startdate), !is.na(enddate)) |>
  mutate(
    start_year = as.integer(format(startdate, "%Y")),
    end_year   = as.integer(format(enddate,   "%Y"))
  ) |>
  filter(start_year <= end_year)        # drop any malformed rows

message(sprintf(
  "Archigos: %d leader-spells, %d countries, %d-%d",
  nrow(archigos),
  n_distinct(archigos$ccode),
  min(archigos$start_year),
  max(archigos$end_year)
))

# --- Expand to country-year ---------------------------------------------------
# Each spell covers [start_year, end_year] inclusive.
# When two leaders are present in the same year (mid-year handover), keep the
# one who entered latest — they held power at year-end.

archigos_panel <- archigos |>
  mutate(year = map2(start_year, end_year, \(s, e) seq.int(s, e))) |>
  unnest(year) |>
  mutate(year = as.integer(year)) |>
  arrange(ccode, year, startdate) |>
  slice_tail(n = 1, by = c(ccode, year)) |>
  mutate(iso3c = countrycode(ccode, "cown", "iso3c", warn = FALSE)) |>
  filter(!is.na(iso3c)) |>
  select(
    iso3c,
    year,
    archigos_leader     = leader,
    archigos_start_year = start_year,
    archigos_entry      = entry,
    archigos_exit       = exit
  )

message(sprintf(
  "Archigos panel: %d country-years, %d countries, %d-%d",
  nrow(archigos_panel),
  n_distinct(archigos_panel$iso3c),
  min(archigos_panel$year),
  max(archigos_panel$year)
))

# --- Save ---------------------------------------------------------------------

saveRDS(archigos_panel, file.path(data_dir, "archigos_panel.rds"))
message("Saved archigos_panel.rds to ", data_dir)
