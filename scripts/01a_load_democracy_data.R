# ==============================================================================
# Load Acemoglu et al. (2019) and V-Dem ERT datasets
# ==============================================================================
# Acemoglu, Naidu, Restrepo, Robinson (2019) "Democracy Does Cause Growth"
# Journal of Political Economy 127(1): 47-100.
#
# V-Dem ERT: Episodes of Regime Transformation, built with default parameters.
# Package: https://github.com/vdeminstitute/ERT
# ==============================================================================

# --- Packages -----------------------------------------------------------------

pkgs <- c("haven", "tidyverse", "devtools")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(haven)
library(tidyverse)

if (!requireNamespace("ERT", quietly = TRUE)) {
  devtools::install_github("vdeminstitute/ERT")
}
library(ERT)

# --- Directories --------------------------------------------------------------

data_dir <- here::here("data")
acemoglu_dir <- file.path(data_dir, "acemoglu")
dir.create(acemoglu_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Acemoglu et al. (2019) replication files -----------------------------

rar_url <- "https://economics.mit.edu/sites/default/files/inline-files/replication_files_ddcg%20%281%29.rar"
rar_path <- file.path(acemoglu_dir, "replication_files_ddcg.rar")

# Download
if (!file.exists(rar_path)) {
  message("Downloading Acemoglu et al. replication files (~30 MB) ...")
  download.file(rar_url, rar_path, mode = "wb")
  message("Download complete.")
} else {
  message("RAR archive already present, skipping download.")
}

# Extract (requires unar: `brew install unar`)
dta_files <- list.files(
  acemoglu_dir,
  pattern = "\\.dta$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(dta_files) == 0) {
  unar_path <- Sys.which("unar")
  unrar_path <- Sys.which("unrar")

  if (nchar(unar_path) == 0 && nchar(unrar_path) == 0) {
    stop(
      "No RAR extraction tool found.\n",
      "Install one with:  brew install unar\n",
      "Then re-run this script."
    )
  }

  message("Extracting replication archive ...")
  if (nchar(unar_path) > 0) {
    system2("unar", args = c("-o", acemoglu_dir, "-f", rar_path))
  } else {
    system2("unrar", args = c("x", rar_path, acemoglu_dir))
  }

  dta_files <- list.files(
    acemoglu_dir,
    pattern = "\\.dta$",
    recursive = TRUE,
    full.names = TRUE
  )
}

message("Stata files found in replication archive:")
print(basename(dta_files))

# Load the main panel dataset.
# The DDCG replication archive contains multiple .dta files; the country-year
# panel used in the paper is DDCGdata.dta.
ddcg_path <- dta_files[grepl(
  "DDCGdata",
  basename(dta_files),
  ignore.case = TRUE
)]
if (length(ddcg_path) == 0) {
  # Fallback: list all and let the user pick
  stop(
    "Could not identify DDCGdata.dta automatically.\n",
    "Files found:\n",
    paste(basename(dta_files), collapse = "\n")
  )
}
ddcg_path <- ddcg_path[1]

ddcg <- read_dta(ddcg_path)
message(sprintf(
  "Loaded Acemoglu et al. panel: %d country-years, %d variables  [%s]",
  nrow(ddcg),
  ncol(ddcg),
  basename(ddcg_path)
))

# Quick look
cat("\n--- Acemoglu et al. (2019): column names ---\n")
print(names(ddcg))

cat("\n--- Coverage ---\n")
cat(sprintf("Years:     %d to %d\n", min(ddcg$year), max(ddcg$year)))
cat(sprintf("Countries: %d\n", n_distinct(ddcg$wbcode)))
cat(sprintf("Country-years: %d\n", nrow(ddcg)))
cat(sprintf(
  "Avg democracy across country-years: %f\n",
  mean(ddcg$dem, na.rm = T)
))

# Democracy indicator: the paper uses `dem` (Acemoglu-Naidu-Restrepo-Robinson
# dichotomous democracy measure, also called the ANRR indicator).
if ("dem" %in% names(ddcg)) {
  transitions <- ddcg %>%
    arrange(wbcode, year) %>%
    group_by(wbcode) %>%
    mutate(
      prev_dem = lag(dem),
      dem_change = dem - prev_dem
    ) %>%
    filter(!is.na(dem_change), dem_change != 0) %>%
    ungroup() %>%
    mutate(
      direction = if_else(dem_change > 0, "democratization", "autocratization")
    )

  cat("\n--- Acemoglu et al. transitions (dem indicator changes) ---\n")
  cat("Total transition events:", nrow(transitions), "\n")
  transitions |>
    count(direction) |>
    mutate(prop = proportions(n)) |>
    print()
} else {
  warning(
    "Expected column `dem` not found. Column names printed above; inspect and adjust."
  )
}

# sapply(ddcg, function(x) attr(x, 'label') %||% '')

# --- 2. V-Dem ERT (default parameters) ---------------------------------------
# Codebook: https://www.v-dem.net/documents/9/ert_codebook.pdf

message("\nLoading V-Dem ERT with default parameters ...")
message("  start_incl = 0.01  (min annual EDI change to trigger onset)")
message(
  "  cum_incl   = 0.10  (min cumulative EDI change for a manifest episode)"
)
message(
  "  year_turn  = 0.03  (annual reversal magnitude that terminates an episode)"
)
message("  cum_turn   = 0.10  (cumulative reversal over tolerance window)")
message(
  "  tolerance  = 5     (years of stasis/opposite movement before termination)"
)

ert <- get_eps() # uses ERT::vdem internally

message(sprintf(
  "Loaded V-Dem ERT: %d episode-rows, %d variables",
  nrow(ert),
  ncol(ert)
))

cat("\n--- V-Dem ERT: column names ---\n")
print(names(ert))

n_dem <- n_distinct(ert$dem_ep_id, na.rm = TRUE)
n_aut <- n_distinct(ert$aut_ep_id, na.rm = TRUE)
n_total <- n_dem + n_aut

cat("\n--- Distinct episodes ---\n")
cat(sprintf(
  "Democratization episodes: %d (%.1f%%)\n",
  n_dem,
  100 * n_dem / n_total
))
cat(sprintf(
  "Autocratization episodes: %d (%.1f%%)\n",
  n_aut,
  100 * n_aut / n_total
))

dem_outcome_labels <- c(
  "0" = "No democratization episode",
  "1" = "Democratic transition",
  "2" = "No democratic transition",
  "3" = "Deepened democracy",
  "4" = "Uncertain"
)

aut_outcome_labels <- c(
  "0" = "No autocratization episode",
  "1" = "Democratic breakdown",
  "2" = "No democratic breakdown",
  "3" = "Regressed autocracy",
  "4" = "Uncertain"
)

cat("\n--- Democratization episodes by outcome ---\n")
ert |>
  filter(dem_ep == 1) |>
  distinct(dem_ep_id, .keep_all = TRUE) |>
  count(dem_ep_outcome_agg) |>
  mutate(
    label = dem_outcome_labels[as.character(dem_ep_outcome_agg)],
    prop = proportions(n)
  ) |>
  print()

cat("\n--- Autocratization episodes by outcome ---\n")
ert |>
  filter(aut_ep == 1) |>
  distinct(aut_ep_id, .keep_all = TRUE) |>
  count(aut_ep_outcome_agg) |>
  mutate(
    label = aut_outcome_labels[as.character(aut_ep_outcome_agg)],
    prop = proportions(n)
  ) |>
  print()

cat("\n--- Coverage ---\n")
cat("Years:    ", range(ert$year, na.rm = TRUE), "\n")
cat("Countries:", n_distinct(ert$country_name), "\n")
cat(sprintf("Country-years: %d\n", nrow(ert)))
cat(sprintf(
  "Avg democracy across country-years: %f\n",
  mean(ert$reg_type, na.rm = T)
))

# --- 3. V-Dem ERT subset: DDCG countries and years --------------------------

# Crosswalk: ERT country_text_id -> DDCG wbcode, for countries that are the
# same entity but use different codes in each dataset.
# Serbia/Montenegro excluded: DDCG's SER = Serbia & Montenegro (pre-2006
# union); ERT's SRB and MNE are the successor states — not a 1-to-1 mapping.
ert_to_ddcg <- c(
  COD = "ZAR", # Democratic Republic of the Congo
  ROU = "ROM", # Romania
  SGP = "SIN", # Singapore
  TWN = "TAW", # Taiwan
  ARE = "UAE" # United Arab Emirates
)

country_coverage_report <- function(ert_df, ddcg_df) {
  ddcg_wbcodes <- unique(ddcg_df$wbcode)
  ert_text_ids <- unique(ert_df$country_text_id)
  ddcg_years <- unique(ddcg_df$year)

  unmatched_ert <- ert_df |>
    filter(year %in% ddcg_years) |>
    filter(!country_text_id %in% ddcg_wbcodes) |>
    distinct(country_text_id, country_name)

  unmatched_ddcg <- ddcg_df |>
    filter(!wbcode %in% ert_text_ids) |>
    distinct(wbcode, country_name)

  cat("ERT countries not matched to DDCG:\n")
  print(unmatched_ert)
  cat("\nDDCG countries not matched to ERT:\n")
  print(unmatched_ddcg)
}

# Apply crosswalk to ERT before subsetting
ert_xw <- ert |>
  mutate(country_text_id = recode(country_text_id, !!!ert_to_ddcg))

cat("\n--- Country coverage before crosswalk ---\n")
country_coverage_report(ert, ddcg)

cat("\n--- Country coverage after crosswalk ---\n")
country_coverage_report(ert_xw, ddcg)

ddcg_years <- unique(ddcg$year)
ddcg_wbcodes <- unique(ddcg$wbcode)

ert_sub <- ert_xw |>
  filter(year %in% ddcg_years)

cat("\n--- ERT subset (DDCG years only, after crosswalk) ---\n")
cat(sprintf("Years:     %d to %d\n", min(ert_sub$year), max(ert_sub$year)))
cat(sprintf("Countries: %d\n", n_distinct(ert_sub$country_text_id)))
cat(sprintf("Country-years: %d\n", nrow(ert_sub)))
cat(sprintf(
  "Avg democracy across country-years: %f\n",
  mean(ert_sub$reg_type, na.rm = T)
))
n_dem_sub <- n_distinct(ert_sub$dem_ep_id, na.rm = TRUE)
n_aut_sub <- n_distinct(ert_sub$aut_ep_id, na.rm = TRUE)
n_total_sub <- n_dem_sub + n_aut_sub

cat("\n--- Distinct episodes (ERT subset, DDCG years) ---\n")
cat(sprintf(
  "Democratization episodes: %d (%.1f%%)\n",
  n_dem_sub,
  100 * n_dem_sub / n_total_sub
))
cat(sprintf(
  "Autocratization episodes: %d (%.1f%%)\n",
  n_aut_sub,
  100 * n_aut_sub / n_total_sub
))

cat("\n--- Democratization episodes by outcome (ERT subset, DDCG years) ---\n")
ert_sub |>
  filter(dem_ep == 1) |>
  distinct(dem_ep_id, .keep_all = TRUE) |>
  count(dem_ep_outcome_agg) |>
  mutate(
    label = dem_outcome_labels[as.character(dem_ep_outcome_agg)],
    prop = proportions(n)
  ) |>
  print()

cat("\n--- Autocratization episodes by outcome (ERT subset, DDCG years) ---\n")
ert_sub |>
  filter(aut_ep == 1) |>
  distinct(aut_ep_id, .keep_all = TRUE) |>
  count(aut_ep_outcome_agg) |>
  mutate(
    label = aut_outcome_labels[as.character(aut_ep_outcome_agg)],
    prop = proportions(n)
  ) |>
  print()


# --- 4. Save cleaned objects for downstream use ------------------------------

saveRDS(ddcg, file.path(data_dir, "ddcg_panel.rds"))
saveRDS(ert, file.path(data_dir, "ert_episodes.rds"))
saveRDS(ert_sub, file.path(data_dir, "ert_sub.rds"))
message(
  "\nSaved ddcg_panel.rds, ert_episodes.rds, and ert_sub.rds to ",
  data_dir
)
