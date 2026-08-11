# ==============================================================================
# Bermeo taxonomy classification
# Classifies each ERT autocratization episode by cause, using Nancy Bermeo's
# six-category taxonomy of how democracies backslide (open-ended coup,
# executive coup, election-day fraud, promissory coup, executive
# aggrandizement, strategic harassment/manipulation), plus "Other/unclear"
# for episodes the taxonomy doesn't cleanly cover (e.g. civil-war state
# collapse, decolonization-driven erosion).
#
# Classifications were produced via LLM web research (Wikipedia-first, with
# source URL + quote per episode) rather than DEED event-level data, which
# this project does not yet have access to (see project memory). Currently
# covers the 131 episodes overlapping DDCG's 1960-2010 window; the raw CSV
# is designed to be appended to as coverage expands to the remaining ~131
# episodes outside that window.
#
# Data:
#   data/bermeo_classifications_raw.csv - one row per classified episode
#
# Outputs:
#   data/episodes_bermeo.rds              - full aut_eps baseline + Bermeo cols
#   output/bermeo_episode_table.csv/.html - episode-level review table (gt)
#   output/bermeo_summary_table.csv/.html - category-frequency summary (gt)
#   output/episodes_full_accounting.csv   - Bermeo + BSE/Funke/DDCG miss tables
#   output/<dataset>_miss_table_bermeo.csv/.html
#     - gt_miss_table from 07/08/09 with Bermeo category + description
#       appended; toggle COMPARISON_DATASET below to "bse", "funke", or "ddcg"
# ==============================================================================

library(tidyverse)
library(gt)
library(here)

data_dir <- here::here("data")
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE)

# Toggle: which comparison dataset's miss table gets the Bermeo treatment.
COMPARISON_DATASET <- "ddcg" # one of "bse", "funke", "ddcg"

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

ert <- readRDS(file.path(data_dir, "ert_episodes.rds"))
bermeo_raw <- read_csv(
  file.path(data_dir, "bermeo_classifications_raw.csv"),
  show_col_types = FALSE
)

# ERT autocratization episodes, all aut_ep==1 (same scope as 08/09).
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
      outcome_agg == 1 ~ "Democratic breakdown",
      outcome_agg == 2 ~ "Democratic erosion",
      outcome_agg == 3 ~ "Autocracy deepening",
      outcome_agg == 4 ~ "Ongoing"
    ),
    prch_label = if_else(prch == 1, "Yes", "No")
  )

n_total_episodes <- nrow(aut_eps)
n_classified <- nrow(bermeo_raw)

# ------------------------------------------------------------------------------
# Join Bermeo classifications onto the full episode baseline
# ------------------------------------------------------------------------------

episodes_bermeo <- aut_eps |>
  left_join(bermeo_raw, by = "aut_ep_id") |>
  mutate(is_classified = !is.na(bermeo_category_primary))

saveRDS(episodes_bermeo, file.path(data_dir, "episodes_bermeo.rds"))

cat(sprintf(
  "Bermeo classification coverage: %d / %d episodes (%.1f%%)\n",
  n_classified,
  n_total_episodes,
  100 * n_classified / n_total_episodes
))

# ------------------------------------------------------------------------------
# Load the comparison dataset selected by COMPARISON_DATASET. Determines
# which miss-table columns are shown in gt_miss_bermeo below, and which
# comparison's match rate appears as a column in gt_bermeo_summary (i.e.
# "of episodes in category X, what share also appear in <dataset>?").
# ------------------------------------------------------------------------------

read_miss_table <- function(filename, rename_cols = character(0)) {
  path <- file.path(out_dir, filename)
  if (!file.exists(path)) {
    return(NULL)
  }
  df <- read_csv(path, show_col_types = FALSE)
  if (length(rename_cols) > 0) {
    df <- rename(df, !!!rename_cols)
  }
  df
}

if (COMPARISON_DATASET == "bse") {
  base_miss <- read_miss_table("bse_ert_miss_matrix.csv")
  stopifnot(!is.null(base_miss))
  match_cols <- c("ROW", "ERT_regime", "Polyarchy", "LibDem", "LibComp")
  match_labels <- c(
    ROW = "ROW",
    ERT_regime = "ERT regime",
    Polyarchy = "Polyarchy",
    LibDem = "LibDem",
    LibComp = "LibComp"
  )
  base_miss <- base_miss |>
    mutate(matched = if_any(all_of(match_cols), ~ !is.na(.)))
  title_txt <- "ERT autocratization episodes (BSE window) with Bermeo classification"
  subtitle_txt <- "Cells show BSE breakdown year(s) within the episode window; Bermeo category/description from LLM research"
  match_label_short <- "BSE"
  comparison_window <- "1999-2023"
} else if (COMPARISON_DATASET == "funke") {
  base_miss <- read_miss_table("funke_ert_miss_table.csv")
  stopifnot(!is.null(base_miss))
  match_cols <- c("funke_spell_years")
  match_labels <- c(funke_spell_years = "Funke spell year(s)")
  base_miss <- base_miss |> mutate(matched = !is.na(funke_spell_years))
  title_txt <- "ERT autocratization episodes (Funke window) with Bermeo classification"
  subtitle_txt <- "Cell shows Funke populist spell year(s) overlapping the episode window; Bermeo category/description from LLM research"
  match_label_short <- "Funke"
  comparison_window <- "1900-2020"
} else if (COMPARISON_DATASET == "ddcg") {
  base_miss <- read_miss_table("ddcg_ert_miss_table.csv")
  stopifnot(!is.null(base_miss))
  match_cols <- c("prch_label", "ddcg_event_years")
  match_labels <- c(
    prch_label = "Originated in democracy",
    ddcg_event_years = "DDCG revevent year(s)"
  )
  base_miss <- base_miss |> mutate(matched = !is.na(ddcg_event_years))
  title_txt <- "ERT autocratization episodes (DDCG window) with Bermeo classification"
  subtitle_txt <- "Cell shows DDCG revevent year(s) within the episode window; Bermeo category/description from LLM research"
  match_label_short <- "DDCG"
  comparison_window <- "1960-2010"
} else {
  stop("COMPARISON_DATASET must be one of 'bse', 'funke', 'ddcg'")
}

# ------------------------------------------------------------------------------
# Shared table styling (same helper as 07/08/09)
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

# ------------------------------------------------------------------------------
# Episode-level review table (classified episodes only)
# ------------------------------------------------------------------------------

bermeo_review <- episodes_bermeo |>
  filter(is_classified) |>
  arrange(ep_start) |>
  select(
    country_name,
    ep_start,
    ep_end,
    outcome_label,
    prch_label,
    bermeo_category_primary,
    bermeo_category_secondary,
    precise_start_date,
    precise_end_date,
    precipitating_event,
    consequential_events,
    description,
    confidence,
    source_url
  )

gt_bermeo_review <- bermeo_review |>
  gt() |>
  tab_header(
    title = "ERT autocratization episodes classified by Bermeo taxonomy",
    subtitle = sprintf(
      "%d of %d episodes classified (1960-2010 DDCG overlap window)",
      n_classified,
      n_total_episodes
    )
  ) |>
  cols_label(
    country_name = "Country",
    ep_start = "ERT start",
    ep_end = "ERT end",
    outcome_label = "ERT outcome",
    prch_label = "Originated in democracy",
    bermeo_category_primary = "Category",
    bermeo_category_secondary = "Secondary",
    precise_start_date = "Precise start",
    precise_end_date = "Precise end",
    precipitating_event = "Precipitating event (LLM)",
    consequential_events = "Consequential events (LLM)",
    description = "Description",
    confidence = "Confidence",
    source_url = "Source(s)"
  ) |>
  sub_missing(columns = bermeo_category_secondary, missing_text = "—") |>
  cols_width(
    country_name ~ px(110),
    ep_start ~ px(60),
    ep_end ~ px(60),
    outcome_label ~ px(90),
    prch_label ~ px(70),
    bermeo_category_primary ~ px(130),
    bermeo_category_secondary ~ px(120),
    precise_start_date ~ px(90),
    precise_end_date ~ px(90),
    precipitating_event ~ px(300),
    consequential_events ~ px(300),
    description ~ px(380),
    confidence ~ px(70),
    source_url ~ px(160)
  ) |>
  cols_align(
    align = "center",
    columns = c(ep_start, ep_end, prch_label, confidence)
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = confidence, rows = confidence == "High")
  ) |>
  tab_style(
    style = cell_text(color = "#999999"),
    locations = cells_body(columns = confidence, rows = confidence == "Low")
  ) |>
  opt_row_striping() |>
  apply_table_style()

# ------------------------------------------------------------------------------
# Category-frequency summary table
# ------------------------------------------------------------------------------

bermeo_review_matched <- bermeo_review |>
  left_join(
    base_miss |> select(country_name, ep_start, ep_end, matched),
    by = c("country_name", "ep_start", "ep_end")
  ) |>
  # Episodes outside the comparison dataset's scope (e.g. BSE's aut_ep_prch==1
  # + 1999-2023 restriction) have no row in base_miss at all, producing NA
  # here -- that correctly means "no match", not "unknown", so coalesce to FALSE.
  mutate(matched = coalesce(matched, FALSE))

bermeo_summary <- bermeo_review |>
  count(bermeo_category_primary, name = "n") |>
  mutate(pct_of_classified = round(100 * n / n_classified, 1)) |>
  left_join(
    bermeo_review |>
      count(bermeo_category_primary, prch_label) |>
      pivot_wider(
        names_from = prch_label,
        values_from = n,
        names_prefix = "prch_",
        values_fill = 0
      ),
    by = "bermeo_category_primary"
  ) |>
  left_join(
    bermeo_review |>
      count(bermeo_category_primary, outcome_label) |>
      pivot_wider(names_from = outcome_label, values_from = n, values_fill = 0),
    by = "bermeo_category_primary"
  ) |>
  left_join(
    bermeo_review_matched |>
      group_by(bermeo_category_primary) |>
      summarise(
        pct_matched_comparison = round(100 * mean(matched), 1),
        .groups = "drop"
      ),
    by = "bermeo_category_primary"
  ) |>
  arrange(desc(n))

outcome_cols <- intersect(
  c(
    "Democratic breakdown",
    "Democratic erosion",
    "Autocracy deepening",
    "Ongoing"
  ),
  names(bermeo_summary)
)

gt_bermeo_summary <- bermeo_summary |>
  gt() |>
  tab_header(
    title = "Bermeo category frequency across classified ERT episodes",
    subtitle = sprintf(
      "n = %d classified episodes; comparison dataset = %s (%s)",
      n_classified,
      match_label_short,
      comparison_window
    )
  ) |>
  tab_spanner(
    label = "Originated in democracy",
    columns = c(prch_Yes, prch_No)
  ) |>
  tab_spanner(label = "ERT outcome", columns = all_of(outcome_cols)) |>
  cols_label(
    bermeo_category_primary = "Bermeo category",
    n = "N",
    pct_of_classified = "% of classified",
    prch_Yes = "Yes",
    prch_No = "No",
    pct_matched_comparison = paste0(
      "% matched to ",
      match_label_short,
      " (",
      comparison_window,
      ")"
    )
  ) |>
  cols_align(align = "center", columns = -bermeo_category_primary) |>
  tab_style(
    style = cell_fill(color = "#dbeafe"),
    locations = cells_body(columns = pct_matched_comparison)
  ) |>
  tab_footnote(
    footnote = paste0(
      "Share of episodes in this category with a match in the ",
      match_label_short,
      " comparison dataset (",
      comparison_window,
      "; see script's COMPARISON_DATASET toggle)."
    ),
    locations = cells_column_labels(columns = pct_matched_comparison)
  ) |>
  opt_row_striping() |>
  apply_table_style()

# ------------------------------------------------------------------------------
# Combined reference table: Bermeo + BSE/Funke/DDCG miss tables
# Joined by (country_name, ep_start, ep_end), which is unique per episode.
# ------------------------------------------------------------------------------

bse_miss <- read_miss_table(
  "bse_ert_miss_matrix.csv",
  c(
    bse_outcome_short = "outcome_short",
    bse_ROW = "ROW",
    bse_ERT_regime = "ERT_regime",
    bse_Polyarchy = "Polyarchy",
    bse_LibDem = "LibDem",
    bse_LibComp = "LibComp"
  )
)
funke_miss <- read_miss_table(
  "funke_ert_miss_table.csv",
  c(funke_outcome_short = "outcome_short")
)
ddcg_miss <- read_miss_table(
  "ddcg_ert_miss_table.csv",
  c(ddcg_outcome_label = "outcome_label", ddcg_prch_label = "prch_label")
)

episodes_full_accounting <- episodes_bermeo

if (!is.null(bse_miss)) {
  episodes_full_accounting <- episodes_full_accounting |>
    left_join(bse_miss, by = c("country_name", "ep_start", "ep_end"))
}
if (!is.null(funke_miss)) {
  episodes_full_accounting <- episodes_full_accounting |>
    left_join(funke_miss, by = c("country_name", "ep_start", "ep_end"))
}
if (!is.null(ddcg_miss)) {
  episodes_full_accounting <- episodes_full_accounting |>
    left_join(ddcg_miss, by = c("country_name", "ep_start", "ep_end"))
}

write_csv(
  episodes_full_accounting,
  file.path(out_dir, "episodes_full_accounting.csv")
)

# ------------------------------------------------------------------------------
# Toggleable gt table: gt_miss_table (07/08/09 style) + Bermeo category/
# description appended. COMPARISON_DATASET picks which comparison's match
# columns to show; Bermeo columns are always appended.
# ------------------------------------------------------------------------------

bermeo_cols_for_join <- episodes_bermeo |>
  select(
    country_name,
    ep_start,
    ep_end,
    bermeo_category_primary,
    bermeo_category_secondary,
    precise_start_date,
    precipitating_event,
    precise_end_date,
    consequential_events,
    source_url,
    confidence
  )

# source_url is semicolon-separated when an episode cites multiple Wikipedia
# pages; render each as a short markdown link (full URLs are too long for a
# table cell) rather than printing the raw URL.
build_source_link <- function(urls) {
  if (is.na(urls)) {
    return(NA_character_)
  }
  parts <- trimws(strsplit(urls, ";")[[1]])
  if (length(parts) == 1) {
    sprintf("[Source](%s)", parts)
  } else {
    paste(sprintf("[Source %d](%s)", seq_along(parts), parts), collapse = ", ")
  }
}

miss_with_bermeo <- base_miss |>
  left_join(
    bermeo_cols_for_join,
    by = c("country_name", "ep_start", "ep_end")
  ) |>
  arrange(ep_start)

n_matched_bermeo <- sum(!is.na(miss_with_bermeo$bermeo_category_primary))
cat(sprintf(
  "\n%s miss table: %d episodes, %d with Bermeo classification\n",
  toupper(COMPARISON_DATASET),
  nrow(miss_with_bermeo),
  n_matched_bermeo
))

# ERT outcome column name differs by source table (outcome_short vs outcome_label).
outcome_col <- intersect(
  c("outcome_short", "outcome_label"),
  names(miss_with_bermeo)
)[1]

gt_miss_bermeo <- miss_with_bermeo |>
  mutate(
    precise_dates = if_else(
      is.na(precise_start_date) | is.na(precise_end_date),
      NA_character_,
      paste0(precise_start_date, " – ", precise_end_date)
    ),
    source_link = map_chr(source_url, build_source_link)
  ) |>
  select(
    country_name,
    ep_start,
    ep_end,
    all_of(outcome_col),
    all_of(match_cols),
    bermeo_category_primary,
    bermeo_category_secondary,
    precise_dates,
    precipitating_event,
    consequential_events,
    source_link,
    confidence
  ) |>
  gt() |>
  tab_header(title = title_txt, subtitle = subtitle_txt) |>
  cols_label(
    country_name = "Country",
    ep_start = "Start",
    ep_end = "End",
    bermeo_category_primary = "Bermeo category",
    bermeo_category_secondary = "Secondary",
    precise_dates = "Precise dates (LLM)",
    precipitating_event = "Precipitating event (LLM)",
    consequential_events = "Consequential events",
    source_link = "Source",
    confidence = "Confidence"
  ) |>
  cols_label(!!outcome_col := "ERT outcome") |>
  cols_label(!!!match_labels) |>
  fmt_markdown(columns = source_link) |>
  sub_missing(columns = all_of(match_cols), missing_text = "—") |>
  sub_missing(
    columns = c(
      bermeo_category_primary,
      bermeo_category_secondary,
      precise_dates,
      precipitating_event,
      consequential_events,
      source_link,
      confidence
    ),
    missing_text = "— not yet classified"
  ) |>
  cols_align(
    align = "center",
    columns = c(
      ep_start,
      ep_end,
      all_of(match_cols),
      precise_dates,
      source_link,
      confidence
    )
  ) |>
  cols_width(
    country_name ~ px(120),
    ep_start ~ px(45),
    ep_end ~ px(45),
    bermeo_category_primary ~ px(130),
    bermeo_category_secondary ~ px(120),
    precise_dates ~ px(170),
    precipitating_event ~ px(300),
    consequential_events ~ px(300),
    source_link ~ px(110),
    confidence ~ px(70)
  ) |>
  tab_style(
    style = list(cell_fill(color = "#dbeafe")),
    locations = cells_body(
      columns = bermeo_category_primary,
      rows = !is.na(bermeo_category_primary)
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = confidence, rows = confidence == "High")
  ) |>
  tab_style(
    style = cell_text(color = "#999999"),
    locations = cells_body(columns = confidence, rows = confidence == "Low")
  ) |>
  opt_row_striping() |>
  apply_table_style()

# Green/gray fill on the comparison-specific match column(s), same convention
# as 07/08/09's gt_miss_table.
for (col in match_cols) {
  sym_col <- sym(col)
  gt_miss_bermeo <- gt_miss_bermeo |>
    tab_style(
      style = list(cell_fill(color = "#c7e9c0")),
      locations = cells_body(
        columns = all_of(col),
        rows = !is.na(!!sym_col) & !(!!sym_col %in% c("Yes", "No"))
      )
    )
}

# ------------------------------------------------------------------------------
# Print console summary
# ------------------------------------------------------------------------------

cat("\n=== Bermeo category counts (classified episodes) ===\n")
print(bermeo_summary |> select(bermeo_category_primary, n, pct_of_classified))

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------

write_csv(bermeo_review, file.path(out_dir, "bermeo_episode_table.csv"))
write_csv(bermeo_summary, file.path(out_dir, "bermeo_summary_table.csv"))
write_csv(
  miss_with_bermeo,
  file.path(out_dir, paste0(COMPARISON_DATASET, "_miss_table_bermeo.csv"))
)
gtsave(gt_bermeo_review, file.path(out_dir, "bermeo_episode_table.html"))
gtsave(gt_bermeo_summary, file.path(out_dir, "bermeo_summary_table.html"))
gtsave(
  gt_miss_bermeo,
  file.path(out_dir, paste0(COMPARISON_DATASET, "_miss_table_bermeo.html"))
)

message("Saved outputs to ", out_dir)

# ------------------------------------------------------------------------------
# Optional: publish bermeo_summary and bermeo_review to a two-tab Google
# Sheet, formatted like the gt tables above -- pretty (Title Case) headers,
# Roboto font, centered text (left-aligned for long free-text columns),
# frozen bold/shaded header row, default gridlines hidden in favor of a
# deliberate header divider + column dividers + outer table border, and a
# blue highlight on the summary's "% matched" column. Off by default. To
# use: set PUBLISH_TO_SHEETS <- TRUE and run googlesheets4::gs4_auth() once
# interactively to authorize (the token caches after that). Re-running
# updates the same spreadsheet rather than creating a new one each time.
# ------------------------------------------------------------------------------

PUBLISH_TO_SHEETS <- TRUE

# Renames columns present in `labels` (a named character vector, old name ->
# pretty label) and leaves any other columns untouched.
apply_pretty_labels <- function(df, labels) {
  nm <- names(df)
  matched <- nm %in% names(labels)
  nm[matched] <- labels[nm[matched]]
  names(df) <- nm
  df
}

# bandedRanges (the mechanism behind zebra striping) reject a second
# addBanding call over a range that already has one -- since re-running
# this script re-applies formatting to the same sheet each time, look up
# and delete any bandedRangeId already on this sheet first.
get_existing_band_ids <- function(spreadsheet_id, sheet_id) {
  req <- gargle::request_build(
    method = "GET",
    path = "v4/spreadsheets/{spreadsheetId}",
    params = list(
      spreadsheetId = spreadsheet_id,
      fields = "sheets(properties.sheetId,bandedRanges.bandedRangeId)"
    ),
    token = googlesheets4::gs4_token(),
    base_url = "https://sheets.googleapis.com/"
  )
  resp <- gargle::request_make(req)
  meta <- gargle::response_process(resp)
  target <- Filter(function(s) s$properties$sheetId == sheet_id, meta$sheets)
  if (length(target) == 0 || is.null(target[[1]]$bandedRanges)) {
    return(integer(0))
  }
  vapply(target[[1]]$bandedRanges, function(b) b$bandedRangeId, integer(1))
}

# Sends one batchUpdate request applying gt-table-style formatting to a
# single sheet tab. Field masks are kept granular (dotted paths like
# "userEnteredFormat.textFormat.bold" rather than the bare "textFormat")
# so later requests only touch the specific property they mean to set,
# instead of clobbering sibling properties (font, alignment, ...) set by
# an earlier request over an overlapping range.
apply_sheet_formatting <- function(
  spreadsheet_id,
  sheet_id,
  df,
  col_widths = c(),
  left_align_columns = character(),
  highlight_columns = character()
) {
  col_pos <- setNames(seq_along(names(df)) - 1, names(df))
  n_rows <- nrow(df)
  n_cols <- ncol(df)
  gray <- list(red = 0.8, green = 0.8, blue = 0.8)
  black <- list(red = 0, green = 0, blue = 0)

  existing_band_ids <- get_existing_band_ids(spreadsheet_id, sheet_id)
  delete_band_requests <- lapply(existing_band_ids, function(id) {
    list(deleteBanding = list(bandedRangeId = id))
  })

  requests <- c(
    list(
      # Freeze header row; hide the sheet's default gridlines -- dividers
      # below are drawn deliberately instead.
      list(
        updateSheetProperties = list(
          properties = list(
            sheetId = sheet_id,
            gridProperties = list(frozenRowCount = 1, hideGridlines = TRUE)
          ),
          fields = "gridProperties.frozenRowCount,gridProperties.hideGridlines"
        )
      )
    ),
    delete_band_requests,
    list(
      # Zebra striping ("slightly shade every other row"), matching gt's
      # opt_row_striping(). Comes before the header/highlight requests
      # below so those more specific background colors win where they
      # overlap.
      list(
        addBanding = list(
          bandedRange = list(
            range = list(
              sheetId = sheet_id,
              startRowIndex = 0,
              endRowIndex = n_rows + 1,
              startColumnIndex = 0,
              endColumnIndex = n_cols
            ),
            rowProperties = list(
              headerColor = list(red = 0.90, green = 0.90, blue = 0.90),
              firstBandColor = list(red = 1, green = 1, blue = 1),
              secondBandColor = list(red = 0.96, green = 0.96, blue = 0.96)
            )
          )
        )
      ),
      # Base look for the whole table: a cleaner font, centered text, and
      # wrapping (harmless on short cells, needed for long-text columns).
      list(
        repeatCell = list(
          range = list(
            sheetId = sheet_id,
            startRowIndex = 0,
            endRowIndex = n_rows + 1,
            startColumnIndex = 0,
            endColumnIndex = n_cols
          ),
          cell = list(
            userEnteredFormat = list(
              textFormat = list(fontFamily = "Roboto", fontSize = 10),
              horizontalAlignment = "CENTER",
              verticalAlignment = "MIDDLE",
              wrapStrategy = "WRAP"
            )
          ),
          fields = paste(
            "userEnteredFormat.textFormat.fontFamily",
            "userEnteredFormat.textFormat.fontSize",
            "userEnteredFormat.horizontalAlignment",
            "userEnteredFormat.verticalAlignment",
            "userEnteredFormat.wrapStrategy",
            sep = ","
          )
        )
      ),
      # Bold, shaded header row with a divider under it (only touches
      # bold/background/border -- font and alignment set above are
      # untouched).
      list(
        repeatCell = list(
          range = list(
            sheetId = sheet_id,
            startRowIndex = 0,
            endRowIndex = 1,
            startColumnIndex = 0,
            endColumnIndex = n_cols
          ),
          cell = list(
            userEnteredFormat = list(
              backgroundColor = list(red = 0.90, green = 0.90, blue = 0.90),
              textFormat = list(bold = TRUE),
              borders = list(
                bottom = list(style = "SOLID_MEDIUM", color = black)
              )
            )
          ),
          fields = paste(
            "userEnteredFormat.textFormat.bold",
            "userEnteredFormat.backgroundColor",
            "userEnteredFormat.borders.bottom",
            sep = ","
          )
        )
      ),
      # Outer table border plus light vertical column dividers -- no
      # horizontal line between every row (that's what "hideGridlines"
      # above is meant to get rid of).
      list(
        updateBorders = list(
          range = list(
            sheetId = sheet_id,
            startRowIndex = 0,
            endRowIndex = n_rows + 1,
            startColumnIndex = 0,
            endColumnIndex = n_cols
          ),
          top = list(style = "SOLID_MEDIUM", color = black),
          bottom = list(style = "SOLID_MEDIUM", color = black),
          left = list(style = "SOLID", color = gray),
          right = list(style = "SOLID", color = gray),
          innerVertical = list(style = "SOLID", color = gray)
        )
      )
    )
  )

  # Left-align long free-text columns -- centered paragraphs read poorly.
  # Must come after the global center-align request above so it wins.
  for (nm in left_align_columns) {
    if (!nm %in% names(col_pos)) {
      next
    }
    idx <- col_pos[[nm]]
    requests <- c(
      requests,
      list(list(
        repeatCell = list(
          range = list(
            sheetId = sheet_id,
            startRowIndex = 0,
            endRowIndex = n_rows + 1,
            startColumnIndex = idx,
            endColumnIndex = idx + 1
          ),
          cell = list(userEnteredFormat = list(horizontalAlignment = "LEFT")),
          fields = "userEnteredFormat.horizontalAlignment"
        )
      ))
    )
  }

  # Explicit column widths -- gt-style fixed sizing sized to the pretty
  # (short) headers, rather than autofit.
  for (nm in names(col_widths)) {
    if (!nm %in% names(col_pos)) {
      next
    }
    idx <- col_pos[[nm]]
    requests <- c(
      requests,
      list(list(
        updateDimensionProperties = list(
          range = list(
            sheetId = sheet_id,
            dimension = "COLUMNS",
            startIndex = idx,
            endIndex = idx + 1
          ),
          properties = list(pixelSize = col_widths[[nm]]),
          fields = "pixelSize"
        )
      ))
    )
  }

  for (nm in highlight_columns) {
    if (!nm %in% names(col_pos)) {
      next
    }
    idx <- col_pos[[nm]]
    requests <- c(
      requests,
      list(list(
        repeatCell = list(
          range = list(
            sheetId = sheet_id,
            startRowIndex = 1,
            endRowIndex = n_rows + 1,
            startColumnIndex = idx,
            endColumnIndex = idx + 1
          ),
          cell = list(
            userEnteredFormat = list(
              backgroundColor = list(red = 0.86, green = 0.91, blue = 0.98)
            )
          ),
          fields = "userEnteredFormat.backgroundColor"
        )
      ))
    )
  }

  # Bypass request_generate()'s named-endpoint lookup (fragile across
  # googlesheets4 versions -- "sheets.spreadsheets.batchUpdate" wasn't
  # recognized in testing) and build the raw batchUpdate request via
  # gargle directly instead.
  req <- gargle::request_build(
    method = "POST",
    path = "v4/spreadsheets/{spreadsheetId}:batchUpdate",
    params = list(spreadsheetId = spreadsheet_id),
    body = list(requests = requests),
    token = googlesheets4::gs4_token(),
    base_url = "https://sheets.googleapis.com/"
  )
  resp <- gargle::request_make(req)
  gargle::response_process(resp)
}

if (PUBLISH_TO_SHEETS) {
  tryCatch(
    {
      library(googlesheets4)
      library(googledrive)

      sheet_title <- paste0("Bermeo summary (", COMPARISON_DATASET, ")")
      existing <- drive_find(sheet_title, type = "spreadsheet", n_max = 1)
      ss <- if (nrow(existing) > 0) existing$id[1] else gs4_create(sheet_title)
      spreadsheet_id <- as.character(as_sheets_id(ss))

      # gt's spanner labels ("Originated in democracy" / "Yes" / "No", "ERT
      # outcome" / "Breakdown" / ...) don't have a flat-header equivalent,
      # so these are kept short and self-explanatory rather than verbose.
      summary_labels <- c(
        bermeo_category_primary = "Bermeo Category",
        n = "N",
        pct_of_classified = "% Classified",
        prch_Yes = "Dem. Origin: Yes",
        prch_No = "Dem. Origin: No",
        "Democratic breakdown" = "Breakdown",
        "Democratic erosion" = "Erosion",
        "Autocracy deepening" = "Deepening",
        Ongoing = "Ongoing",
        pct_matched_comparison = paste0("% Matched (", match_label_short, ")")
      )
      bermeo_summary_pretty <- apply_pretty_labels(
        bermeo_summary,
        summary_labels
      )

      # Built from miss_with_bermeo (the same data behind gt_miss_bermeo)
      # rather than bermeo_review, so the Episodes tab has the exact same
      # columns as that gt table -- including the comparison-dataset match
      # column(s), which bermeo_review never had (it's the dataset-agnostic
      # classification review, not tied to COMPARISON_DATASET).
      build_source_formula <- function(urls) {
        if (is.na(urls)) {
          return(NA_character_)
        }
        parts <- trimws(strsplit(urls, ";")[[1]])
        if (length(parts) == 1) {
          sprintf('=HYPERLINK("%s","Source")', parts)
        } else {
          links <- sprintf(
            'HYPERLINK("%s","Source %d")',
            parts,
            seq_along(parts)
          )
          paste0("=", paste(links, collapse = '&", "&'))
        }
      }

      episodes_for_sheet <- miss_with_bermeo |>
        mutate(
          precise_dates = if_else(
            is.na(precise_start_date) | is.na(precise_end_date),
            NA_character_,
            paste0(precise_start_date, " - ", precise_end_date)
          ),
          source_formula = gs4_formula(map_chr(
            source_url,
            build_source_formula
          ))
        ) |>
        select(
          country_name,
          ep_start,
          ep_end,
          all_of(outcome_col),
          all_of(match_cols),
          bermeo_category_primary,
          bermeo_category_secondary,
          precise_dates,
          precipitating_event,
          consequential_events,
          source_formula,
          confidence
        )

      # Short Title Case labels for the comparison-specific match column(s)
      # -- separate from match_labels (used by gt_miss_bermeo's HTML) so
      # the Sheets header stays as concise as the rest of the pretty labels.
      match_pretty_labels <- switch(
        COMPARISON_DATASET,
        bse = c(
          ROW = "ROW",
          ERT_regime = "ERT Regime",
          Polyarchy = "Polyarchy",
          LibDem = "LibDem",
          LibComp = "LibComp"
        ),
        funke = c(funke_spell_years = "Funke Spell Year(s)"),
        ddcg = c(
          prch_label = "In Democracy",
          ddcg_event_years = "DDCG Match Year(s)"
        )
      )

      episode_labels <- c(
        country_name = "Country",
        ep_start = "ERT Start",
        ep_end = "ERT End",
        setNames("ERT Outcome", outcome_col),
        match_pretty_labels,
        bermeo_category_primary = "Bermeo Category",
        bermeo_category_secondary = "Secondary Category",
        precise_dates = "Precise Dates (LLM)",
        precipitating_event = "Precipitating Event (LLM)",
        consequential_events = "Consequential Events (LLM)",
        confidence = "Confidence",
        source_formula = "Source"
      )
      episodes_for_sheet_pretty <- apply_pretty_labels(
        episodes_for_sheet,
        episode_labels
      )

      sheet_write(bermeo_summary_pretty, ss = ss, sheet = "Summary")
      sheet_write(episodes_for_sheet_pretty, ss = ss, sheet = "Episodes")

      # A freshly created spreadsheet starts with an unused default tab;
      # drop it now that Summary/Episodes exist.
      props <- sheet_properties(ss)
      if ("Sheet1" %in% props$name) {
        try(sheet_delete(ss, sheet = "Sheet1"), silent = TRUE)
        props <- sheet_properties(ss)
      }
      summary_sheet_id <- props |> filter(name == "Summary") |> pull(id)
      episodes_sheet_id <- props |> filter(name == "Episodes") |> pull(id)

      summary_col_widths <- c(
        "Bermeo Category" = 150,
        N = 45,
        "% Classified" = 100,
        "Dem. Origin: Yes" = 110,
        "Dem. Origin: No" = 110,
        Breakdown = 95,
        Erosion = 85,
        Deepening = 95,
        Ongoing = 85,
        setNames(140, paste0("% Matched (", match_label_short, ")"))
      )
      apply_sheet_formatting(
        spreadsheet_id,
        summary_sheet_id,
        bermeo_summary_pretty,
        col_widths = summary_col_widths,
        highlight_columns = paste0("% Matched (", match_label_short, ")")
      )

      episodes_col_widths <- c(
        Country = 130,
        "ERT Start" = 85,
        "ERT End" = 85,
        "ERT Outcome" = 120,
        "In Democracy" = 100,
        "DDCG Match Year(s)" = 130,
        "Funke Spell Year(s)" = 140,
        ROW = 70,
        "ERT Regime" = 90,
        Polyarchy = 90,
        LibDem = 90,
        LibComp = 90,
        "Bermeo Category" = 150,
        "Secondary Category" = 150,
        "Precise Dates (LLM)" = 170,
        "Precipitating Event (LLM)" = 320,
        "Consequential Events (LLM)" = 320,
        Confidence = 90,
        Source = 90
      )
      apply_sheet_formatting(
        spreadsheet_id,
        episodes_sheet_id,
        episodes_for_sheet_pretty,
        col_widths = episodes_col_widths,
        left_align_columns = c(
          "Precipitating Event (LLM)",
          "Consequential Events (LLM)"
        )
      )

      message(
        "Published Bermeo summary + episodes to Google Sheet: ",
        gs4_get(ss)$spreadsheet_url
      )
    },
    error = function(e) {
      message(
        "Could not publish to Google Sheets (set PUBLISH_TO_SHEETS <- TRUE and run ",
        "googlesheets4::gs4_auth() once to authorize): ",
        conditionMessage(e)
      )
    }
  )
}


library(dplyr)
library(stringr)

d <- readr::read_csv("output/episodes_full_accounting.csv") |>
  filter(!is.na(ddcg_event_years), !is.na(precise_start_date)) |>
  mutate(
    ddcg_start = as.numeric(str_extract(
      as.character(ddcg_event_years),
      "\\d{4}"
    )),
    llm_start = as.numeric(str_extract(precise_start_date, "\\d{4}")),
    ddcg_diff = abs(ddcg_start - llm_start),
    ert_diff = abs(
      ep_start -
        llm_start
    ),
  )

d |>
  filter(!is.na(ddcg_diff), !is.na(ert_diff)) |>
  summarise(
    ddcg_closer = mean(ddcg_diff < ert_diff),
    ert_closer = mean(ddcg_diff > ert_diff),
    same = mean(ddcg_diff == ert_diff)
  )
