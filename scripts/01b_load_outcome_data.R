# ==============================================================================
# Load economic outcome data
# ==============================================================================
# Outcomes from Pon, Marx & Rollet (Electoral Turnovers):
#
# 1. GDP per capita growth — Penn World Tables 11.0
#    Variable: rgdpna — "Real GDP at constant 2017 national prices (in million
#    2017 USD)." *na variables are "primarily designed to capture growth over
#    time" (constant national prices isolate real volume changes from inflation).
#    Coverage: 185 countries, 1950-2023.
#    Codebook: https://www.rug.nl/ggdc/productivity/pwt/?lang=en
#    Download: https://dataverse.nl/api/access/datafile/554030  (Stata)
#
# 2. CPI inflation — IMF International Financial Statistics (via World Bank WDI)
#    Variable: FP.CPI.TOTL.ZG — "Inflation, consumer prices (annual %)."
#    Sourced from IMF IFS and data files.
#    Codebook: https://api.worldbank.org/v2/indicator/FP.CPI.TOTL.ZG?format=json
#
# 3. Unemployment rate — ILO modelled estimates (via World Bank WDI)
#    Variable: SL.UEM.TOTL.ZS — "Unemployment, total (% of total labor force)
#    (modelled ILO estimate)." Sourced from ILO ILOSTAT.
#    Codebook: https://api.worldbank.org/v2/indicator/SL.UEM.TOTL.ZS?format=json
#
# 4. Trade intensity — World Bank WDI
#    Variable: NE.TRD.GNFS.ZS — "Trade (% of GDP)" = (exports + imports of
#    goods and services) / GDP * 100.
#    Codebook: https://api.worldbank.org/v2/indicator/NE.TRD.GNFS.ZS?format=json
#
# 5. Top income shares — World Inequality Database (WID)
#    Variable: sptinc992j — share of pre-tax national income, equal-split adults.
#    Percentiles: p90p100 (top 10%), p99p100 (top 1%), p0p50 (bottom 50%).
#    Built from tax records + national accounts (distributional national accounts
#    methodology, Piketty/Saez/Zucman). Coverage: rich countries from ~1900,
#    most others from ~1980. Values stored as percent (0–100).
#    Codebook: https://wid.world/codes-dictionary/
#    Package: GitHub only — remotes::install_github("world-inequality-database/wid-r-tool")
#    Note: raw download is cached to data/wid_inequality_raw.rds.
#
# 6. Gini coefficients — Standardized World Income Inequality Database (SWIID)
#    Variables: gini_disp (disposable income, post-tax/transfer) and gini_mkt
#    (market income, pre-tax/transfer), both on 0–100 scale.
#    Harmonized across survey sources with multiple imputation; swiid_summary.csv
#    gives point estimates and is kept current in the repo (no versioned filename).
#    Reference: Solt (2020), Social Science Quarterly.
#    Coverage: ~180 countries, mostly 1960–present.
#    Download: https://github.com/fsolt/swiid (data/swiid_summary.csv)
#
# 1–4 replicate Pon, Marx & Rollet (Electoral Turnovers). 5–6 are added to
# support the inequality → backsliding research angle (reverse causality).
#
# Each source is saved separately as an RDS in data/. Country identifier is
# the ISO 3-letter code (iso3c) throughout, to facilitate later merging with
# the event dataset (which uses V-Dem / World Bank country codes).
# ==============================================================================

# --- Packages -----------------------------------------------------------------

pkgs <- c("tidyverse", "here", "haven", "WDI", "countrycode", "remotes")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
if (!requireNamespace("wid", quietly = TRUE)) {
  remotes::install_github("world-inequality-database/wid-r-tool")
}
library(tidyverse)
library(haven)
library(WDI)
library(wid)
library(countrycode)

# --- Directories --------------------------------------------------------------

data_dir <- here::here("data")

# --- 1. GDP per capita growth (Penn World Tables 11.0) ------------------------
# Released October 2025. No R package yet; downloaded directly from Dataverse.
#
# rgdpna: real GDP at constant 2017 national prices (millions 2017 USD)
# pop:    population (millions)
# -> gdp_pc in 2017 USD; gdp_pc_growth as annual percentage change.

pwt_path <- file.path(data_dir, "pwt110.dta")

if (!file.exists(pwt_path)) {
  message("Downloading Penn World Tables 11.0 (~15 MB) ...")
  download.file(
    "https://dataverse.nl/api/access/datafile/554030",
    pwt_path,
    mode = "wb"
  )
  message("Download complete.")
} else {
  message("PWT 11.0 already present, skipping download.")
}

pwt <- read_dta(pwt_path) |>
  as_tibble() |>
  select(country, iso3c = countrycode, year, rgdpna, pop) |>
  arrange(iso3c, year) |>
  group_by(iso3c) |>
  mutate(
    gdp_pc        = rgdpna / pop,
    gdp_pc_growth = (gdp_pc / lag(gdp_pc) - 1) * 100,
    ln_gdp_pc     = log(gdp_pc)
  ) |>
  ungroup() |>
  select(-rgdpna, -pop)

message(sprintf(
  "PWT: %d country-years, %d countries, %d-%d",
  nrow(pwt),
  n_distinct(pwt$iso3c),
  min(pwt$year),
  max(pwt$year)
))
message(sprintf(
  "  gdp_pc_growth: %d non-missing values",
  sum(!is.na(pwt$gdp_pc_growth))
))

saveRDS(pwt, file.path(data_dir, "pwt_gdp.rds"))
message("Saved pwt_gdp.rds\n")

# --- 2. CPI inflation + Unemployment + Trade intensity (World Bank WDI) -------
# FP.CPI.TOTL.ZG: Inflation, consumer prices (annual %)     — sourced from IMF IFS
# SL.UEM.TOTL.ZS: Unemployment, total (% of labour force)   — ILO modelled estimates
# NE.TRD.GNFS.ZS: Trade (% of GDP) = (exports + imports) / GDP * 100

message("Downloading World Bank WDI (CPI + Unemployment + Trade) ...")
wdi_raw <- WDI(
  indicator = c(
    cpi_inflation     = "FP.CPI.TOTL.ZG",
    unemployment_rate = "SL.UEM.TOTL.ZS",
    trade_pct_gdp     = "NE.TRD.GNFS.ZS"
  ),
  country = "all",
  start   = 1960,
  end     = 2025,
  extra   = TRUE
)

wdi <- wdi_raw |>
  as_tibble() |>
  filter(region != "Aggregates") |>
  select(country, iso3c, year, cpi_inflation, unemployment_rate, trade_pct_gdp) |>
  arrange(iso3c, year)

cpi   <- wdi |> select(country, iso3c, year, cpi_inflation)
unemp <- wdi |> select(country, iso3c, year, unemployment_rate)
trade <- wdi |> select(country, iso3c, year, trade_pct_gdp)

message(sprintf(
  "IMF CPI (via WDI): %d non-missing country-years, %d countries",
  sum(!is.na(cpi$cpi_inflation)),
  n_distinct(cpi$iso3c[!is.na(cpi$cpi_inflation)])
))
message(sprintf(
  "ILO unemployment (via WDI): %d non-missing country-years, %d countries",
  sum(!is.na(unemp$unemployment_rate)),
  n_distinct(unemp$iso3c[!is.na(unemp$unemployment_rate)])
))
message(sprintf(
  "WB Trade (via WDI): %d non-missing country-years, %d countries",
  sum(!is.na(trade$trade_pct_gdp)),
  n_distinct(trade$iso3c[!is.na(trade$trade_pct_gdp)])
))

saveRDS(cpi,   file.path(data_dir, "imf_cpi.rds"))
saveRDS(unemp, file.path(data_dir, "ilo_unemployment.rds"))
saveRDS(trade, file.path(data_dir, "wb_trade.rds"))
message("Saved imf_cpi.rds, ilo_unemployment.rds, and wb_trade.rds\n")

# --- 5. Top income shares (World Inequality Database) -------------------------
# sptinc992j: share of pre-tax national income, equal-split adults (ages=992, pop="j")
# p90p100 = top 10%; p99p100 = top 1%; p0p50 = bottom 50%.
# Values returned as fractions (0–1); converted to % below.
# Coverage: rich countries from ~1900, most others from ~1980. Sub-national
# entries (e.g. "US-NY") are dropped via countrycode conversion.
# Codebook: https://wid.world/codes-dictionary/

wid_cache <- file.path(data_dir, "wid_inequality_raw.rds")
if (!file.exists(wid_cache)) {
  message("Downloading WID income share data (all countries) ...")
  wid_raw <- download_wid(
    indicators = "sptinc",
    areas      = "all",
    perc       = c("p90p100", "p99p100", "p0p50"),
    ages       = 992,
    pop        = "j"
  )
  saveRDS(wid_raw, wid_cache)
  message("WID download complete.")
} else {
  message("WID raw cache found, skipping download.")
  wid_raw <- readRDS(wid_cache)
}

wid_ineq <- wid_raw |>
  as_tibble() |>
  mutate(
    iso3c = countrycode(country, "iso2c", "iso3c", warn = FALSE),
    var   = case_when(
      percentile == "p90p100" ~ "top10_share",
      percentile == "p99p100" ~ "top1_share",
      percentile == "p0p50"   ~ "bot50_share"
    ),
    value = value * 100   # fraction → percent
  ) |>
  filter(!is.na(iso3c), !is.na(var)) |>
  select(iso3c, year, var, value) |>
  pivot_wider(names_from = var, values_from = value) |>
  arrange(iso3c, year)

message(sprintf(
  "WID: %d country-years, %d countries, %d-%d",
  nrow(wid_ineq), n_distinct(wid_ineq$iso3c),
  min(wid_ineq$year), max(wid_ineq$year)
))
message(sprintf(
  "  top10_share: %d non-missing values",
  sum(!is.na(wid_ineq$top10_share))
))

saveRDS(wid_ineq, file.path(data_dir, "wid_inequality.rds"))
message("Saved wid_inequality.rds\n")

# --- 6. Gini coefficients (SWIID, current via swiid_summary.csv) -------------
# gini_disp: disposable income Gini (post-tax, post-transfer), 0–100 scale
# gini_mkt:  market income Gini (pre-tax, pre-transfer), 0–100 scale
# Reference: Solt (2020), SPI. Coverage: ~180 countries, mostly 1960–present.
# swiid_summary.csv is the point-estimates file kept current in the repo.
# Download: https://github.com/fsolt/swiid

swiid_path <- file.path(data_dir, "swiid_summary.csv")
if (!file.exists(swiid_path)) {
  message("Downloading SWIID summary ...")
  download.file(
    "https://raw.githubusercontent.com/fsolt/swiid/master/data/swiid_summary.csv",
    swiid_path, mode = "wb"
  )
  message("SWIID download complete.")
} else {
  message("SWIID already present, skipping download.")
}

swiid_gini <- read_csv(swiid_path, show_col_types = FALSE) |>
  mutate(iso3c = countrycode(country, "country.name", "iso3c", warn = FALSE)) |>
  filter(!is.na(iso3c)) |>
  select(iso3c, year, gini_disp, gini_mkt) |>
  arrange(iso3c, year)

message(sprintf(
  "SWIID: %d country-years, %d countries, %d-%d",
  nrow(swiid_gini), n_distinct(swiid_gini$iso3c),
  min(swiid_gini$year), max(swiid_gini$year)
))

saveRDS(swiid_gini, file.path(data_dir, "swiid_gini.rds"))
message("Saved swiid_gini.rds\n")

message("All outcome datasets saved to ", data_dir)
