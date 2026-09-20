# ==============================================================================
# Additional outcome variables for the RDD pipeline.
#
# Two sources, two output files:
#
#   1. data/elections_database/outcomes.dta (the Electoral Turnovers /
#      National Elections Database country-year outcomes panel, 254 MB,
#      859 columns) -> data/et_outcomes.rds
#        - alternative GDP-per-capita series (World Bank, IMF WEO) so growth
#          can be shown from three sources alongside the PWT series already in
#          combined_panel.rds
#        - government debt and deficit, both as a share of GDP
#      Read with col_select so only the handful of columns we need is parsed.
#
#   2. ERT::vdem (V-Dem v16, bundled with the ERT package -- the same V-Dem
#      release get_eps() builds the episodes from, so no version skew against
#      combined_panel.rds) -> data/vdem_subcomponents.rds
#        - judicial and legislative constraints on the executive, and their
#          mean (the "checks and balances" index used in the ET paper)
#        - HOS/HOG power indices, two constructions (see below)
#
# Output:
#   data/et_outcomes.rds          country_text_id x year
#   data/vdem_subcomponents.rds   country_text_id x year
# ==============================================================================

library(tidyverse)
library(haven)
library(here)
library(ERT)

source(here::here("scripts", "vdem_indices.R"))

data_dir <- here::here("data")
elections_dir <- file.path(data_dir, "elections_database")

report_coverage <- function(df, vars) {
  for (v in vars) {
    ok <- !is.na(df[[v]])
    if (!any(ok)) {
      cat(sprintf("  %-22s  EMPTY\n", v))
      next
    }
    cat(sprintf(
      "  %-22s n=%6d  %d-%d  %d countries\n",
      v,
      sum(ok),
      min(df$year[ok]),
      max(df$year[ok]),
      n_distinct(df$country_text_id[ok])
    ))
  }
}

# ------------------------------------------------------------------------------
# Part 1 -- Electoral Turnovers outcomes panel
# ------------------------------------------------------------------------------

# `id` is a haven_labelled integer whose VALUE LABELS are ISO3 codes
# ("ABW", "AFG", ...), not the integers themselves -- as_factor() recovers the
# code. It is non-missing for ~85% of rows; the rest are historical entities
# with no modern ISO3 and are dropped, since they could never join to
# combined_panel anyway.
#
# Variable choices, from the 859 available:
#   logGDPc_wb   "Log GDP per capita" (World Bank)  -- already logged, so it
#                differences exactly like combined_panel's ln_gdp_pc (PWT).
#   GDPc_imfweo  "GDP, constant prices, per capita, domestic currency" (IMF
#                WEO) -- a LEVEL in local currency, so log it here. Levels
#                aren't comparable across countries, but a within-country log
#                difference is, which is all the outcome uses.
#   bgen_gmd     "General government Debt as % of GDP" (Global Macro Database).
#                Chosen over b_imfweo (1980+ only) and b_wb (108 countries) for
#                coverage; over b_gmd because that one is a level, not a ratio.
#   dgen_gmd     "General government Deficit as % of GDP" (same source).
cat("Reading outcomes.dta (col_select -- only the columns we need)...\n")
et_raw <- read_dta(
  file.path(elections_dir, "outcomes.dta"),
  col_select = c(Country, Year, id, logGDPc_wb, GDPc_imfweo, bgen_gmd, dgen_gmd)
)

et_outcomes <- et_raw |>
  mutate(country_text_id = as.character(as_factor(id))) |>
  filter(!is.na(country_text_id), country_text_id != "") |>
  transmute(
    country_text_id,
    year = as.numeric(Year),
    ln_gdp_pc_wb = logGDPc_wb,
    # log of a local-currency constant-price level: only ever differenced
    # within country, so the units cancel.
    ln_gdp_pc_imf = if_else(GDPc_imfweo > 0, log(GDPc_imfweo), NA_real_),
    debt_pct_gdp = bgen_gmd,
    deficit_pct_gdp = dgen_gmd
  ) |>
  # A handful of (country, year) pairs appear twice where the source panel
  # carries two country_name spellings mapping to one ISO3. Collapse by mean,
  # the same convention 02a_build_panel.R already uses for SWIID.
  summarise(
    across(everything(), \(x) mean(x, na.rm = TRUE)),
    .by = c(country_text_id, year)
  ) |>
  mutate(across(everything(), \(x) if (is.numeric(x)) ifelse(is.nan(x), NA_real_, x) else x))

cat(sprintf(
  "\nET outcomes: %d country-years, %d countries\n",
  nrow(et_outcomes),
  n_distinct(et_outcomes$country_text_id)
))
report_coverage(
  et_outcomes,
  c("ln_gdp_pc_wb", "ln_gdp_pc_imf", "debt_pct_gdp", "deficit_pct_gdp")
)

# ---- deficit sign diagnostic -------------------------------------------------
# The Global Macro Database labels dgen_gmd "Deficit", while the sibling _get
# and _imfpfmh sources label the analogous column "primary balance" -- opposite
# sign conventions for the same quantity. Rather than trust the label, check
# against country-years that unambiguously ran large deficits (US 2009-2010,
# Greece 2009, Ireland 2010) and orient so that POSITIVE = LARGER DEFICIT,
# which is the direction that makes "backsliding worsens public finances" a
# positive coefficient, matching every other outcome's sign convention.
deficit_probe <- et_outcomes |>
  filter(
    (country_text_id == "USA" & year %in% 2009:2010) |
      (country_text_id == "GRC" & year == 2009) |
      (country_text_id == "IRL" & year == 2010)
  ) |>
  select(country_text_id, year, deficit_pct_gdp)

cat("\nDeficit sign diagnostic (known large-deficit country-years):\n")
print(as.data.frame(deficit_probe))

probe_mean <- mean(deficit_probe$deficit_pct_gdp, na.rm = TRUE)
if (is.finite(probe_mean) && probe_mean < 0) {
  cat(sprintf(
    "  -> mean %.2f is negative, so dgen_gmd is signed as a BALANCE; flipping so positive = larger deficit\n",
    probe_mean
  ))
  et_outcomes <- et_outcomes |> mutate(deficit_pct_gdp = -deficit_pct_gdp)
} else {
  cat(sprintf(
    "  -> mean %.2f is positive, so dgen_gmd is already signed as a DEFICIT; left as-is\n",
    probe_mean
  ))
}

saveRDS(et_outcomes, file.path(data_dir, "et_outcomes.rds"))
message("Saved data/et_outcomes.rds")

# ------------------------------------------------------------------------------
# Part 2 -- V-Dem subcomponents and executive-power indices
# ------------------------------------------------------------------------------

# --- checks and balances ---
# The Electoral Turnovers outcomes panel ships a pre-built checks_balances_vdem;
# it is exactly the simple mean of jucon_vdem and lg_legcon_vdem (verified
# against all 20,631 non-missing rows, max absolute deviation 3e-8). Rebuild it
# from V-Dem directly rather than carrying it over, so it comes from the same
# V-Dem release as everything else here.

# --- HOS/HOG power, construction 1 (primary): the Electoral Turnovers recipe ---
# "We considered that leaders hold power in different forms, following the
# variables from V-Dem: power to dissolve the legislature, appoint ministers,
# dismiss ministers, veto legislation, propose legislation. For each form of
# power, we normalized the V-Dem variable on the [0,1] segment, with 0 meaning
# least power and 1 most power. The means of these normalized variables give us
# indices reflecting the level of power of the HOS and HOG."
#
# Two things the paper's variable list gets wrong for the HOG side, both of
# which a literal transcription would hit silently:
#
#  (a) It names v2exdfdshg for "dissolve the legislature", v2exdfcbhg for
#      "appoint ministers" and v2exdfdmhg for "dismiss ministers". In V-Dem,
#      v2exdfcbhg and v2exdfdmhg DO NOT EXIST, and v2exdfdshg is "HOG dismisses
#      ministers in practice" -- NOT dissolution. The real tags are v2exdjdshg
#      (HOG dissolution) and v2exdjcbhg (HOG appoints cabinet). The mapping
#      below matches on CONCEPT, which is what the paper's prose describes.
#
#  (b) The propose-legislation items run BACKWARDS relative to the other eight.
#      For everything else 0 = "No" = no power; for v2exdfpphs/v2exdfpphg,
#      0 = "Yes, in all policy areas, including some exclusive domains" (most
#      power) and the top category is "cannot propose legislation" (none). They
#      are reverse-coded below; normalizing them like the rest would drag the
#      index the wrong way.
#
# `ord_max` is the top category from the V-Dem codebook, used as the
# normalization denominator (not the observed max, which could drift).
POWER_COMPONENTS <- tribble(
  ~concept,               ~hos_tag,        ~hos_max, ~hog_tag,        ~hog_max, ~reverse,
  "dissolve_legislature", "v2exdfdshs_ord", 3,       "v2exdjdshg_ord", 3,       FALSE,
  "appoint_ministers",    "v2exdfcbhs_ord", 4,       "v2exdjcbhg_ord", 2,       FALSE,
  "dismiss_ministers",    "v2exdfdmhs_ord", 3,       "v2exdfdshg_ord", 3,       FALSE,
  "veto_legislation",     "v2exdfvths_ord", 4,       "v2exdfvthg_ord", 4,       FALSE,
  "propose_legislation",  "v2exdfpphs_ord", 2,       "v2exdfpphg_ord", 2,       TRUE
)

normalize_component <- function(x, ord_max, reverse) {
  z <- x / ord_max
  if (reverse) 1 - z else z
}

vdem_raw <- ERT::vdem

needed <- c(
  "country_text_id", "year",
  "v2x_jucon", "v2xlg_legcon", "v2ex_hosw", "v2ex_hogw",
  POWER_COMPONENTS$hos_tag, POWER_COMPONENTS$hog_tag,
  # The high- and mid-level democracy indices used as RDD outcomes. The list
  # is in scripts/vdem_indices.R because 02a and rdd_helpers.R need the same
  # one. NEW_VARS excludes the three already reaching combined_panel.rds by
  # another route, which would otherwise arrive twice on the join.
  VDEM_INDEX_NEW_VARS
)
missing_vars <- setdiff(needed, names(vdem_raw))
if (length(missing_vars) > 0) {
  stop("Missing from ERT::vdem: ", paste(missing_vars, collapse = ", "))
}

vd <- vdem_raw[, needed]

hos_mat <- mapply(
  normalize_component,
  vd[POWER_COMPONENTS$hos_tag],
  POWER_COMPONENTS$hos_max,
  POWER_COMPONENTS$reverse
)
hog_mat <- mapply(
  normalize_component,
  vd[POWER_COMPONENTS$hog_tag],
  POWER_COMPONENTS$hog_max,
  POWER_COMPONENTS$reverse
)

# Mean over the components that are coded for that country-year. Requiring all
# five would throw away a lot of otherwise-usable rows, but a mean over one or
# two components is not the same index -- so require at least three and record
# the component count for auditing.
mean_if_enough <- function(m, min_components = 3) {
  n_ok <- rowSums(!is.na(m))
  out <- rowMeans(m, na.rm = TRUE)
  out[n_ok < min_components] <- NA_real_
  list(value = out, n_components = n_ok)
}

hos_idx <- mean_if_enough(hos_mat)
hog_idx <- mean_if_enough(hog_mat)

# --- HOS/HOG power, construction 2 (alternative) ---
# V-Dem's own v2ex_hosw / v2ex_hogw ("relative power of the HOS/HOG"), which
# take the four values {0, 0.25, 0.5, 1}. That spacing is V-Dem's, not linear;
# recode it to the evenly-spaced {0, 1/3, 2/3, 1} so the alternative is a
# genuinely LINEAR index, matching the spirit of ET's "_linear" suffix.
linear_recode_hosw <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x == 0 ~ 0,
    x == 0.25 ~ 1 / 3,
    x == 0.5 ~ 2 / 3,
    x == 1 ~ 1,
    # V-Dem occasionally carries an intermediate value from its own
    # aggregation; fall back to the identity rather than silently NA-ing it.
    TRUE ~ x
  )
}

vdem_sub <- tibble(
  country_text_id = vd$country_text_id,
  year = vd$year,
  v2x_jucon = vd$v2x_jucon,
  v2xlg_legcon = vd$v2xlg_legcon,
  checks_balances = rowMeans(
    cbind(vd$v2x_jucon, vd$v2xlg_legcon),
    na.rm = FALSE
  ),
  hos_power_linear = hos_idx$value,
  hog_power_linear = hog_idx$value,
  hos_power_n_components = hos_idx$n_components,
  hog_power_n_components = hog_idx$n_components,
  hos_power_vdem = linear_recode_hosw(vd$v2ex_hosw),
  hog_power_vdem = linear_recode_hosw(vd$v2ex_hogw)
) |>
  # The democracy indices pass through untouched -- unlike the power indices
  # above they need no reconstruction, they ARE V-Dem's published aggregates.
  bind_cols(vd[, VDEM_INDEX_NEW_VARS, drop = FALSE]) |>
  filter(!is.na(country_text_id), !is.na(year))

cat(sprintf(
  "\nV-Dem subcomponents: %d country-years, %d countries\n",
  nrow(vdem_sub),
  n_distinct(vdem_sub$country_text_id)
))
report_coverage(
  vdem_sub,
  c(
    "v2x_jucon", "v2xlg_legcon", "checks_balances",
    "hos_power_linear", "hog_power_linear",
    "hos_power_vdem", "hog_power_vdem"
  )
)

cat("\nV-Dem democracy indices (RDD outcomes):\n")
report_coverage(vdem_sub, VDEM_INDEX_NEW_VARS)

# ---- executive-power sanity check -------------------------------------------
# If the propose-legislation reverse-coding above were missed, these cases would
# come out muddled instead of sharp: a strong presidential system should sit
# high on hos_power_linear with hog_power_linear absent or low; a parliamentary
# system should invert that.
cat("\nExecutive-power sanity check (2015):\n")
print(
  vdem_sub |>
    filter(
      year == 2015,
      country_text_id %in% c("USA", "FRA", "DEU", "ITA", "GBR", "RUS", "TUR")
    ) |>
    select(
      country_text_id,
      hos_power_linear,
      hog_power_linear,
      hos_power_vdem,
      hog_power_vdem
    ) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    as.data.frame()
)

saveRDS(vdem_sub, file.path(data_dir, "vdem_subcomponents.rds"))
message("Saved data/vdem_subcomponents.rds")
