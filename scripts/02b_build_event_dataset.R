# ==============================================================================
# Build comparison event dataset: V-Dem ERT × Acemoglu et al. (2019)
# ==============================================================================

library(haven)
library(tidyverse)

# --- Setup --------------------------------------------------------------------

data_dir <- here::here("data")
TOLERANCE <- 1 # years of slack when matching a DDCG event year to an ERT window

# Country code crosswalk (applied in both directions)
ert_to_ddcg <- c(
  COD = "ZAR",
  ROU = "ROM",
  SGP = "SIN",
  TWN = "TAW",
  ARE = "UAE"
)
ddcg_to_ert <- setNames(names(ert_to_ddcg), unname(ert_to_ddcg))

# --- Load data ----------------------------------------------------------------

ddcg <- readRDS(file.path(data_dir, "ddcg_panel.rds"))
ert <- readRDS(file.path(data_dir, "ert_episodes.rds"))

ert_xw <- ert |>
  mutate(country_text_id = recode(country_text_id, !!!ert_to_ddcg))

# --- 1. DDCG events -----------------------------------------------------------

ddcg_events <- ddcg |>
  filter(demevent == 1 | revevent == 1) |>
  transmute(
    country_text_id = recode(wbcode, !!!ddcg_to_ert),
    country_name,
    year,
    direction = case_when(
      demevent == 1 ~ "democratization",
      revevent == 1 ~ "autocratization"
    )
  )

# --- 2. ERT episodes (one row per episode) ------------------------------------

extract_eps <- function(
  ert_df,
  ep_col,
  ep_id_col,
  start_col,
  end_col,
  outcome_col,
  direction_label
) {
  ert_df |>
    filter(.data[[ep_col]] == 1) |>
    group_by(.data[[ep_id_col]]) |>
    summarise(
      country_text_id = first(country_text_id),
      country_name = first(country_name),
      ep_id = first(.data[[ep_id_col]]),
      start_year = first(.data[[start_col]]),
      end_year = first(.data[[end_col]]),
      outcome_agg = first(.data[[outcome_col]]),
      reg_type = reg_type[year == first(.data[[start_col]])][1],
      polyarchy_start = v2x_polyarchy[year == first(.data[[start_col]])][1],
      polyarchy_end = v2x_polyarchy[year == first(.data[[end_col]])][1],
      .groups = "drop"
    ) |>
    mutate(direction = direction_label)
}

ert_eps <- bind_rows(
  extract_eps(
    ert_xw,
    "dem_ep",
    "dem_ep_id",
    "dem_ep_start_year",
    "dem_ep_end_year",
    "dem_ep_outcome_agg",
    "democratization"
  ),
  extract_eps(
    ert_xw,
    "aut_ep",
    "aut_ep_id",
    "aut_ep_start_year",
    "aut_ep_end_year",
    "aut_ep_outcome_agg",
    "autocratization"
  )
)

# --- 3. Match and build union via full outer join -----------------------------

ert_eps_w <- ert_eps |>
  mutate(match_start = start_year - TOLERANCE, match_end = end_year + TOLERANCE)

nrow(ert_eps_w)

# Full outer join: keeps matched pairs, DDCG-only rows, and ERT-only rows.
# country_name appears in both tables; coalesce after joining.
events <- full_join(
  ddcg_events,
  ert_eps_w |>
    select(
      country_text_id,
      direction,
      country_name,
      ep_id,
      start_year,
      end_year,
      outcome_agg,
      reg_type,
      polyarchy_start,
      polyarchy_end,
      match_start,
      match_end
    ),
  by = join_by(
    country_text_id,
    direction,
    between(year, match_start, match_end)
  )
) |>
  select(-match_start, -match_end) |>
  mutate(
    country_name = coalesce(country_name.x, country_name.y),
    in_ddcg = !is.na(year),
    in_ert = !is.na(ep_id),
    year = coalesce(year, start_year)
  ) |>
  select(-country_name.x, -country_name.y) |>
  arrange(country_text_id, year)

nrow(events)

# --- 5. Summary ---------------------------------------------------------------

print_summary <- function(df, label) {
  cat(sprintf("\n%s\n", label))
  df |>
    count(direction, in_ddcg, in_ert) |>
    mutate(
      source = case_when(
        in_ddcg & in_ert ~ "Both",
        in_ddcg & !in_ert ~ "DDCG only",
        !in_ddcg & in_ert ~ "ERT only"
      )
    ) |>
    select(direction, source, n) |>
    arrange(direction, source) |>
    pivot_wider(names_from = source, values_from = n, values_fill = 0) |>
    mutate(
      Total = rowSums(across(where(is.integer))),
      across(
        where(is.integer),
        \(x) sprintf("%d (%.0f%%)", x, 100 * x / Total),
        .names = "{.col}"
      )
    ) |>
    print()
}

cat("\n=== Union event dataset ===\n")
cat(sprintf("Total rows: %d\n", nrow(events)))
cat("(>1 row per DDCG event if it matched multiple ERT episodes)\n")

print_summary(events, "All years:")

ddcg_years <- unique(ddcg$year)
print_summary(
  events |> filter(year %in% ddcg_years),
  "DDCG window (1960-2010):"
)

saveRDS(events, file.path(data_dir, "events.rds"))
message("\nSaved events.rds to ", data_dir)
