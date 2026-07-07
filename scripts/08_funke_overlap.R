# ==============================================================================
# Funke et al. (2023) overlap analysis
# How do populist leader spells relate to our ERT autocratization episodes?
#
# Key questions from the Jun 29 meeting:
#   (1) Are our ERT episodes the same events as when populist leaders take power?
#   (2) For episodes that DO overlap: which comes first — the populist entry
#       or the ERT episode start? (diagnostic for causal direction)
#   (3) Which ERT episodes have no Funke populist — i.e., are driven by
#       something other than a populist leader taking power?
#
# Funke sample: 60 specific countries, 1900-2020.
# ERT sample:   183 countries, 1900-2025.
# Overlap analysis restricted to (a) Funke's 60 countries and (b) 1900-2020
# (the period most relevant for our research design).
#
# Outputs:
#   output/funke_overlap_summary.csv    — episode-level match table
#   output/funke_timing.csv            — timing gap for matched pairs
#   output/funke_unmatched_ert.csv     — ERT episodes with no Funke populist
#   output/funke_unmatched_spells.csv  — Funke spells with no ERT episode
# ==============================================================================

library(tidyverse)
library(gt)
library(here)

data_dir <- here::here("data")
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE)

EP_MATCH_TOL <- 2 # years of slack in period-overlap matching

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

ert <- readRDS(file.path(data_dir, "ert_episodes.rds"))
spells <- readRDS(file.path(data_dir, "funke_spells.rds"))
panel <- readRDS(file.path(data_dir, "funke_panel.rds"))

funke_countries <- distinct(panel, country_text_id, country)

# ERT autocratization episodes originating from a democracy (aut_ep_prch == 1)
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
      outcome_agg == 1 ~ "Breakdown",
      outcome_agg == 2 ~ "Erosion only",
      outcome_agg == 4 ~ "Uncertain / ongoing"
    ),
    in_funke_sample = country_text_id %in% funke_countries$country_text_id
  )

# ------------------------------------------------------------------------------
# Match populist spells to ERT autocratization episodes
# Two periods overlap if: spell_start <= ep_end + TOL & spell_end >= ep_start - TOL
# ------------------------------------------------------------------------------

matched <- spells |>
  left_join(
    aut_eps |>
      select(
        country_text_id,
        aut_ep_id,
        ep_start,
        ep_end,
        outcome_agg,
        outcome_label,
        poly_start,
        poly_end
      ),
    by = "country_text_id",
    relationship = "many-to-many"
  ) |>
  filter(
    is.na(ep_start) |
      (entry_year <= ep_end + EP_MATCH_TOL &
        exit_year >= ep_start - EP_MATCH_TOL)
  ) |>
  mutate(
    has_ert_match = !is.na(ep_start),
    entry_before_ep = entry_year <= ep_start, # populist entered before episode
    years_populist_precedes = ep_start - entry_year # positive = populist came first
  )

# ------------------------------------------------------------------------------
# Summary 1: What fraction of Funke populist spells overlap with an ERT episode?
# ------------------------------------------------------------------------------

spell_summary <- matched |>
  group_by(country_text_id, country, entry_year, exit_year, variant) |>
  summarise(
    n_ert_matches = sum(has_ert_match),
    any_ert = any(has_ert_match),
    any_outcome_1 = any(outcome_agg == 1, na.rm = TRUE),
    any_outcome_2 = any(outcome_agg == 2, na.rm = TRUE),
    any_outcome_4 = any(outcome_agg == 4, na.rm = TRUE),
    matched_outcomes = paste(
      sort(unique(outcome_label[has_ert_match])),
      collapse = "; "
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Summary 2: What fraction of ERT episodes (in Funke countries) have a populist?
# ------------------------------------------------------------------------------

# For each ERT episode in Funke's 60-country sample, check for a matching spell
ep_summary <- aut_eps |>
  filter(in_funke_sample) |>
  left_join(
    spells |> select(country_text_id, entry_year, exit_year, variant),
    by = "country_text_id",
    relationship = "many-to-many"
  ) |>
  filter(
    is.na(entry_year) |
      (entry_year <= ep_end + EP_MATCH_TOL &
        exit_year >= ep_start - EP_MATCH_TOL)
  ) |>
  group_by(
    aut_ep_id,
    country_name,
    ep_start,
    ep_end,
    outcome_agg,
    outcome_label,
    poly_start,
    poly_end
  ) |>
  summarise(
    n_populist_spells = sum(!is.na(entry_year)),
    any_populist = any(!is.na(entry_year)),
    populist_variants = paste(
      sort(unique(variant[!is.na(variant)])),
      collapse = "; "
    ),
    earliest_entry_year = min(entry_year, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    earliest_entry_year = if_else(
      is.infinite(earliest_entry_year),
      NA_real_,
      earliest_entry_year
    ),
    years_populist_leads = ep_start - earliest_entry_year # positive = populist first
  )

# ------------------------------------------------------------------------------
# Period constants, gt helper, miss table, and summary table
# ------------------------------------------------------------------------------

FUNKE_START <- 1900
FUNKE_END <- 2020

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

# ---- Episode-level miss table (ERT → Funke) ---------------------------------
# left_join with join_by inequality silently drops left rows for countries that
# have spells but none overlapping the episode window — same fan-out issue fixed
# in 07 via the bse_years_by_family split. Fix: compute matches first, then
# left-join the result back to the full 101-row baseline to restore all episodes.

# Compute matches against all 101 ERT episodes (not just Funke's 60 countries).
# Episodes outside Funke's sample will get NA for funke_spell_years naturally
# since no spell row exists for those country_text_ids.
funke_spell_matches <- aut_eps |>
  filter(ep_start <= FUNKE_END) |>
  mutate(ep_end_cap = pmin(ep_end, FUNKE_END, na.rm = TRUE)) |>
  left_join(
    spells |>
      filter(entry_year >= FUNKE_START) |>
      select(country_text_id, entry_year, exit_year),
    by = join_by(
      country_text_id,
      ep_start <= exit_year,
      ep_end_cap >= entry_year
    )
  ) |>
  group_by(aut_ep_id) |>
  summarise(
    funke_spell_years = {
      e <- entry_year[!is.na(entry_year)]
      x <- exit_year[!is.na(entry_year)]
      if (length(e) == 0) {
        NA_character_
      } else {
        paste(sort(unique(paste0(e, "–", x))), collapse = ", ")
      }
    },
    .groups = "drop"
  )

ert_miss_table <- aut_eps |>
  filter(ep_start <= FUNKE_END) |>
  select(aut_ep_id, country_name, ep_start, ep_end, outcome_agg) |>
  left_join(funke_spell_matches, by = "aut_ep_id") |>
  mutate(
    outcome_short = case_when(
      outcome_agg == 1 ~ "Breakdown",
      outcome_agg == 2 ~ "Erosion only",
      outcome_agg == 4 ~ "Ongoing"
    )
  ) |>
  arrange(ep_start) |>
  select(country_name, ep_start, ep_end, outcome_short, funke_spell_years)

# ---- Funke → ERT direction (for summary counts) -----------------------------
# Same two-step fix as funke_spell_matches: compute matches first, then
# left-join back to the full spell list to restore non-matching spells.

spells_base <- spells |> filter(entry_year >= FUNKE_START)

ert_matches_for_spells <- spells_base |>
  left_join(
    aut_eps |>
      filter(ep_start <= FUNKE_END) |>
      mutate(ep_end_cap = pmin(ep_end, FUNKE_END, na.rm = TRUE)) |>
      select(country_text_id, aut_ep_id, ep_start, ep_end_cap),
    by = join_by(
      country_text_id,
      entry_year <= ep_end_cap,
      exit_year >= ep_start
    )
  ) |>
  group_by(country_text_id, country, entry_year, exit_year, variant) |>
  summarise(any_ert = any(!is.na(aut_ep_id)), .groups = "drop")

funke_ert_match <- spells_base |>
  left_join(
    ert_matches_for_spells,
    by = c("country_text_id", "country", "entry_year", "exit_year", "variant")
  ) |>
  mutate(any_ert = coalesce(any_ert, FALSE))

# ---- Summary table (one row, both directions) --------------------------------

n_ert_in_funke_count <- sum(!is.na(ert_miss_table$funke_spell_years))
n_funke_period <- nrow(funke_ert_match)
n_funke_in_ert <- sum(funke_ert_match$any_ert)
pct_funke_in_ert <- round(100 * mean(funke_ert_match$any_ert), 1)

summary_table <- tibble(
  n_ert_episodes = nrow(ert_miss_table),
  n_funke_spells = n_funke_period,
  n_overlap = n_ert_in_funke_count,
  pct_funke_in_ert = pct_funke_in_ert,
  pct_funke_not_in_ert = round(100 - pct_funke_in_ert, 1),
  pct_ert_in_funke = round(
    100 * n_ert_in_funke_count / nrow(ert_miss_table),
    1
  ),
  pct_ert_not_in_funke = round(
    100 * (nrow(ert_miss_table) - n_ert_in_funke_count) / nrow(ert_miss_table),
    1
  )
)

# ---- gt summary table -------------------------------------------------------

gt_summary <- summary_table |>
  gt() |>
  tab_header(
    title = "ERT autocratization episodes vs. Funke et al. (2023) populist spells",
    subtitle = paste0(
      "ERT episodes in Funke’s 60 countries with ep_start ≥ ",
      FUNKE_START,
      "; ",
      "Funke spells with entry year ≥ ",
      FUNKE_START,
      "."
    )
  ) |>
  tab_spanner(
    label = "Funke → ERT",
    columns = c(pct_funke_in_ert, pct_funke_not_in_ert)
  ) |>
  tab_spanner(
    label = "ERT → Funke",
    columns = c(pct_ert_in_funke, pct_ert_not_in_funke)
  ) |>
  cols_label(
    n_ert_episodes = "ERT episodes",
    n_funke_spells = "Funke spells",
    n_overlap = "Overlap",
    pct_funke_in_ert = "% in ERT",
    pct_funke_not_in_ert = "% not in ERT",
    pct_ert_in_funke = "% in Funke",
    pct_ert_not_in_funke = "% not in Funke"
  ) |>
  cols_align(align = "center", columns = everything()) |>
  cols_width(
    n_ert_episodes ~ px(85),
    n_funke_spells ~ px(85),
    n_overlap ~ px(70),
    pct_funke_in_ert ~ px(75),
    pct_funke_not_in_ert ~ px(90),
    pct_ert_in_funke ~ px(75),
    pct_ert_not_in_funke ~ px(90)
  ) |>
  apply_table_style()

# ---- gt miss table ----------------------------------------------------------

gt_miss_table <- ert_miss_table |>
  gt() |>
  tab_header(
    title = paste0(
      "ERT autocratization episodes in Funke’s 60 countries (post-",
      FUNKE_START,
      ")"
    ),
    subtitle = "Cell shows Funke populist spell year(s) overlapping the episode window (entry–exit); — = not captured"
  ) |>
  cols_label(
    country_name = "Country",
    ep_start = "Start",
    ep_end = "End",
    outcome_short = "ERT outcome",
    funke_spell_years = "Funke spell year(s)"
  ) |>
  sub_missing(columns = funke_spell_years, missing_text = "—") |>
  cols_align(align = "center", columns = funke_spell_years) |>
  cols_width(
    country_name ~ px(120),
    ep_start ~ px(45),
    ep_end ~ px(45),
    outcome_short ~ px(90),
    funke_spell_years ~ px(95)
  ) |>
  tab_style(
    style = list(cell_fill(color = "#c7e9c0"), cell_text(weight = "bold")),
    locations = cells_body(
      columns = funke_spell_years,
      rows = !is.na(funke_spell_years)
    )
  ) |>
  tab_style(
    style = cell_fill(color = "#eeeeee"),
    locations = cells_body(
      columns = funke_spell_years,
      rows = is.na(funke_spell_years)
    )
  ) |>
  tab_style(
    style = cell_text(size = px(9)),
    locations = cells_body(columns = funke_spell_years)
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "#666666", weight = px(1.5)),
    locations = cells_body(columns = funke_spell_years)
  ) |>
  tab_style(
    style = cell_borders(sides = "left", color = "#666666", weight = px(1.5)),
    locations = cells_column_labels(columns = funke_spell_years)
  ) |>
  opt_row_striping() |>
  apply_table_style()

# ------------------------------------------------------------------------------
# Summary 3: Timing — for matched pairs, who came first?
# (restricted to post-1900 for relevance)
# ------------------------------------------------------------------------------

timing <- matched |>
  filter(has_ert_match, entry_year >= 1900) |>
  select(
    country,
    entry_year,
    exit_year,
    variant,
    ep_start,
    ep_end,
    outcome_agg,
    outcome_label,
    entry_before_ep,
    years_populist_precedes
  ) |>
  arrange(years_populist_precedes)

# ------------------------------------------------------------------------------
# Print results
# ------------------------------------------------------------------------------

cat("=== Overview ===\n")
cat(sprintf(
  "Funke spells (all years):     %d across %d countries\n",
  nrow(spells),
  n_distinct(spells$country_text_id)
))
cat(sprintf(
  "Funke spells (post-1900):     %d\n",
  sum(spells$entry_year >= 1900)
))
cat(sprintf(
  "ERT aut episodes (aut_ep_prch==1, all years): %d\n",
  nrow(aut_eps)
))
cat(sprintf(
  "  of which in Funke's 60 countries:           %d\n",
  sum(aut_eps$in_funke_sample)
))
cat(sprintf(
  "  of which NOT in Funke's 60 countries:       %d\n",
  sum(!aut_eps$in_funke_sample)
))

cat(
  "\n\n=== (1) Funke populist spells → do they overlap with an ERT episode? ===\n"
)
cat("(All years)\n")
spell_summary |>
  count(any_ert) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  print()

cat("\nBreakdown by ERT outcome (among matched spells):\n")
spell_summary |>
  filter(any_ert) |>
  count(any_outcome_1, any_outcome_2) |>
  print()

cat("\nPost-1900 spells only:\n")
spell_summary |>
  filter(entry_year >= 1900) |>
  count(any_ert) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  print()

cat("\nFull match table (post-1900):\n")
spell_summary |>
  filter(entry_year >= 1900) |>
  select(country, entry_year, exit_year, variant, any_ert, matched_outcomes) |>
  arrange(entry_year) |>
  print(n = Inf)

cat(
  "\n\n=== (2) ERT episodes (in Funke countries) → do they have a populist spell? ===\n"
)
ep_summary |>
  count(outcome_label, any_populist) |>
  pivot_wider(
    names_from = any_populist,
    values_from = n,
    names_prefix = "populist_",
    values_fill = 0
  ) |>
  rename(no_populist = populist_FALSE, has_populist = populist_TRUE) |>
  mutate(
    pct_with_populist = round(
      100 * has_populist / (has_populist + no_populist),
      1
    )
  ) |>
  print()

cat("\nERT episodes with no Funke populist spell (in Funke's 60 countries):\n")
ep_summary |>
  filter(!any_populist) |>
  select(country_name, ep_start, ep_end, outcome_label, poly_start, poly_end) |>
  arrange(ep_start) |>
  print(n = Inf)

cat(
  "\n\n=== (3) Timing: who comes first — populist entry or ERT episode start? ===\n"
)
cat(
  "(Post-1900, matched pairs only; positive = populist entered BEFORE episode)\n\n"
)

timing |>
  group_by(outcome_label) |>
  summarise(
    n = n(),
    mean_lead = round(mean(years_populist_precedes, na.rm = TRUE), 1),
    median_lead = median(years_populist_precedes, na.rm = TRUE),
    pct_populist_first = round(100 * mean(entry_before_ep, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()

cat("\nCountry-level detail:\n")
timing |>
  select(
    country,
    variant,
    entry_year,
    ep_start,
    ep_end,
    outcome_label,
    years_populist_precedes
  ) |>
  arrange(desc(years_populist_precedes)) |>
  print(n = Inf)

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------

write_csv(spell_summary, file.path(out_dir, "funke_overlap_spell_summary.csv"))
write_csv(ep_summary, file.path(out_dir, "funke_overlap_ep_summary.csv"))
write_csv(timing, file.path(out_dir, "funke_timing.csv"))
write_csv(
  ep_summary |> filter(!any_populist),
  file.path(out_dir, "funke_unmatched_ert.csv")
)
write_csv(
  spell_summary |> filter(!any_ert),
  file.path(out_dir, "funke_unmatched_spells.csv")
)
write_csv(summary_table, file.path(out_dir, "funke_summary_table.csv"))
write_csv(ert_miss_table, file.path(out_dir, "funke_ert_miss_table.csv"))
gtsave(gt_summary, file.path(out_dir, "funke_summary_table.html"))
gtsave(gt_miss_table, file.path(out_dir, "funke_ert_miss_table.html"))

message("Saved outputs to ", out_dir)
