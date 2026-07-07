# ==============================================================================
# Load and clean Funke, Schularick & Trebesch (2023) populist leaders data
#
# Source: updated panel (1900-2020) from Manuel Funke's website
#   https://sites.google.com/view/manuel-funke/data
# Cite: Funke, Schularick & Trebesch, "Populist Leaders and the Economy",
#   American Economic Review 113(12): 3249-3288, 2023.
#
# Key variables in raw panel:
#   pop   = 1 if a populist leader is in power (any variant)
#   lpop  = 1 if left-wing populist
#   rpop  = 1 if right-wing populist
#   independent = 1 if country is in the analytic sample for that year
#
# Outputs:
#   data/funke_panel.rds    — country-year panel, 60 countries, 1900-2020
#   data/funke_spells.rds   — one row per populist spell (entry/exit years)
# ==============================================================================

library(haven)
library(tidyverse)
library(here)

data_dir <- here::here("data")

# ------------------------------------------------------------------------------
# 1. Load raw panel
# ------------------------------------------------------------------------------

raw <- read_dta(file.path(data_dir, "funke_populists_2020.dta")) |>
  rename(pop_any = pop) |>      # rename to avoid shadowing dplyr::pop
  mutate(
    country_text_id = iso,      # Funke ISO3 codes match ERT country_text_id
    across(c(pop_any, lpop, rpop, independent), as.integer)
  )

# ------------------------------------------------------------------------------
# 2. Build clean country-year panel
# ------------------------------------------------------------------------------

funke_panel <- raw |>
  select(country, country_text_id, iso, year, independent,
         pop_any, lpop, rpop) |>
  arrange(country_text_id, year)

# ------------------------------------------------------------------------------
# 3. Extract populist spells (entry and exit years)
# A new spell begins whenever pop_any goes from 0 -> 1 (or at the start of
# the series if pop_any == 1 in the first observed year for that country).
# ------------------------------------------------------------------------------

funke_spells <- funke_panel |>
  arrange(country_text_id, year) |>
  group_by(country_text_id) |>
  mutate(
    prev_pop  = lag(pop_any, default = 0L),
    spell_start = pop_any == 1L & prev_pop == 0L,
    spell_end   = pop_any == 0L & prev_pop == 1L,
    spell_id    = cumsum(spell_start)
  ) |>
  ungroup() |>
  filter(pop_any == 1L) |>
  group_by(country_text_id, country, iso, spell_id) |>
  summarise(
    entry_year = min(year),
    exit_year  = max(year),       # last year in power (inclusive)
    lpop       = max(lpop),       # 1 if any year in spell was left-wing
    rpop       = max(rpop),
    variant    = case_when(
      max(lpop) == 1 & max(rpop) == 1 ~ "Mixed",
      max(lpop) == 1                  ~ "Left",
      max(rpop) == 1                  ~ "Right",
      TRUE                            ~ NA_character_
    ),
    .groups = "drop"
  ) |>
  select(-spell_id) |>
  arrange(country_text_id, entry_year)

# ------------------------------------------------------------------------------
# 4. Verify against AER Table 1 counts
# ------------------------------------------------------------------------------

cat("=== Funke populist spells (2020 panel) ===\n")
cat(sprintf("Countries in panel: %d\n", n_distinct(funke_panel$country_text_id)))
cat(sprintf("Total spells: %d  |  left: %d  |  right: %d\n",
            nrow(funke_spells),
            sum(funke_spells$variant == "Left",  na.rm = TRUE),
            sum(funke_spells$variant == "Right", na.rm = TRUE)))
cat(sprintf("Countries with at least one spell: %d\n",
            n_distinct(funke_spells$country_text_id)))

cat("\nPost-1990 spells:\n")
funke_spells |>
  filter(entry_year >= 1990) |>
  select(country, entry_year, exit_year, variant) |>
  arrange(entry_year) |>
  print(n = Inf)

# ------------------------------------------------------------------------------
# 5. Check country code alignment with ERT
# ------------------------------------------------------------------------------

ert_ids <- readRDS(file.path(data_dir, "ert_episodes.rds")) |>
  distinct(country_text_id)

unmatched <- funke_panel |>
  distinct(country_text_id) |>
  anti_join(ert_ids, by = "country_text_id")

if (nrow(unmatched) > 0) {
  cat("\nWARNING: Funke countries not found in ERT:\n")
  print(unmatched)
} else {
  cat("\nAll Funke country codes match ERT country_text_id.\n")
}

# ------------------------------------------------------------------------------
# 6. Save
# ------------------------------------------------------------------------------

saveRDS(funke_panel,  file.path(data_dir, "funke_panel.rds"))
saveRDS(funke_spells, file.path(data_dir, "funke_spells.rds"))
message("Saved funke_panel.rds and funke_spells.rds to ", data_dir)
