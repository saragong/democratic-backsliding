# ==============================================================================
# Cutpoint summary: every candidate "illiberal" threshold side by side.
#
# For each rule -- the PopuList cut carried over as a raw number, the same cut
# transferred at its global percentile, the within-country percentile cut
# (accuracy and Youden criteria), and a round 0.5 for reference -- report:
#
#   how well it separates PopuList's populists (from 01g_populist_threshold.R),
#   what share of all V-Party party-years it calls illiberal,
#   how many RDD elections have ONE top-2 party above it and the other not,
#     and how many of those are treated,
#   what raw anti-pluralism score it amounts to in a few reference countries,
#   and whether Front National clears it.
#
# Runs after 01g (the cuts) and 11 (the build, which carries the
# within-country percentile columns). Reads, never estimates.
#
# Output: output/sweeps/populist_threshold/
#           cutpoint_summary.csv / .html       one row per rule
#           implied_cut_by_country.csv / .html the raw anti-pluralism score each
#                                              within-country cut implies, by country
# ==============================================================================

library(tidyverse)
library(here)
library(gt)

source(here::here("scripts", "rdd_helpers.R"))
source(here::here("scripts", "vparty_helpers.R"))

if (!exists("SUMMARY_WINDOW")) SUMMARY_WINDOW <- DEFAULT_WINDOW
if (!exists("SUMMARY_TREATMENT")) SUMMARY_TREATMENT <- DEFAULT_TREATMENT
REFERENCE_COUNTRIES <- c("FRA", "HUN", "USA", "TUR", "IND", "BRA")
# Front National / Rassemblement National, the case the Sep 20 note raises: at
# 0.49 it sits below the 0.65 PopuList cut.
FN_VPARTY_ID <- 433

out_dir <- sweep_dir("populist_threshold")
thr <- readRDS(here::here("data", "populist_threshold.rds"))

vp <- load_vparty_raw(c("v2paid", "country_text_id", "year", "v2xpa_antiplural")) |>
  filter(!is.na(v2xpa_antiplural))
vp$pct_ctry <- country_pct(
  vp$v2xpa_antiplural, vp$country_text_id, vp$v2xpa_antiplural, vp$country_text_id
)

d <- load_build(DEFAULT_INSTRUMENT, SUMMARY_WINDOW)
stopifnot(all(c("illiberal_pct_ctry", "other_pct_ctry") %in% names(d)))

fn <- vp |> filter(v2paid == FN_VPARTY_ID) |> slice_max(year, n = 1)

acc_row <- function(method) {
  thr$cutpoints |> filter(method == !!method) |> slice(1)
}
acc_row_ctry <- function(method) {
  thr$cutpoints_ctry |> filter(method == !!method) |> slice(1)
}

RULES <- tribble(
  ~rule, ~spec, ~scale, ~cut, ~sens, ~spec_rate, ~acc,
  "PopuList cut, raw score (accuracy)", "popucut", "raw", thr$threshold_abs,
  acc_row("accuracy")$sensitivity, acc_row("accuracy")$specificity, acc_row("accuracy")$accuracy,
  "PopuList cut, global-percentile transfer", "popucut_pct", "raw", thr$threshold_pct,
  NA, NA, NA,
  "PopuList cut, within-country percentile (accuracy)", "popucut_ctry", "ctry_pct", thr$threshold_pct_ctry,
  acc_row_ctry("accuracy")$sensitivity, acc_row_ctry("accuracy")$specificity, acc_row_ctry("accuracy")$accuracy,
  "PopuList cut, within-country percentile (Youden)", "popucut_ctry_youden", "ctry_pct", thr$threshold_pct_ctry_youden,
  acc_row_ctry("youden")$sensitivity, acc_row_ctry("youden")$specificity, acc_row_ctry("youden")$accuracy,
  "Round reference cut", "0.5", "raw", 0.5,
  NA, NA, NA
)

# The raw anti-pluralism score a rule amounts to in one country: the cut itself
# on the raw scale, or the country's own quantile at the percentile cut.
implied_raw <- function(scale, cut, country) {
  if (scale == "raw") {
    return(cut)
  }
  x <- vp$v2xpa_antiplural[vp$country_text_id == country]
  if (length(x) < COUNTRY_PCT_MIN_N) {
    return(NA_real_)
  }
  unname(quantile(x, probs = cut, type = 1))
}

summary_tbl <- RULES |>
  rowwise() |>
  mutate(
    party_years_above = if (scale == "raw") {
      mean(vp$v2xpa_antiplural > cut)
    } else {
      mean(vp$pct_ctry > cut, na.rm = TRUE)
    },
    elections_pair = {
      iv <- d[[threshold_var("illiberal", scale)]]
      ov <- d[[threshold_var("other", scale)]]
      sum(iv > cut & ov <= cut, na.rm = TRUE)
    },
    treated_pair = {
      iv <- d[[threshold_var("illiberal", scale)]]
      ov <- d[[threshold_var("other", scale)]]
      sum((iv > cut & ov <= cut) & d[[SUMMARY_TREATMENT]] == 1, na.rm = TRUE)
    },
    fn_clears = if (scale == "raw") {
      fn$v2xpa_antiplural > cut
    } else {
      fn$pct_ctry > cut
    },
    !!!setNames(
      lapply(REFERENCE_COUNTRIES, function(ct) rlang::expr(implied_raw(scale, cut, !!ct))),
      paste0("cut_", REFERENCE_COUNTRIES)
    )
  ) |>
  ungroup()

write_csv(summary_tbl, file.path(out_dir, "cutpoint_summary.csv"))

fmt_cut <- function(scale, cut) {
  ifelse(scale == "raw", sprintf("%.3f (raw)", cut), sprintf("%.1fth pct (within country)", 100 * cut))
}

gt_summary <- summary_tbl |>
  transmute(
    Rule = rule,
    Spec = paste0("`", spec, "`"),
    Cut = fmt_cut(scale, cut),
    Sensitivity = sens, Specificity = spec_rate, Accuracy = acc,
    `Party-years above` = party_years_above,
    `Elections, one above one below` = elections_pair,
    `... treated` = treated_pair,
    `Front National clears it` = if_else(fn_clears, "yes", "no"),
    across(starts_with("cut_"))
  ) |>
  gt() |>
  fmt_markdown(columns = Spec) |>
  fmt_number(columns = c(Sensitivity, Specificity, Accuracy), decimals = 2) |>
  fmt_percent(columns = `Party-years above`, decimals = 1) |>
  fmt_number(columns = starts_with("cut_"), decimals = 3) |>
  sub_missing(missing_text = "--") |>
  tab_spanner(label = "PopuList fit", columns = c(Sensitivity, Specificity, Accuracy)) |>
  tab_spanner(
    label = sprintf("RDD sample, w%d", SUMMARY_WINDOW),
    columns = c(`Elections, one above one below`, `... treated`)
  ) |>
  tab_spanner(
    label = "Implied raw anti-pluralism cut",
    columns = starts_with("cut_")
  ) |>
  cols_label(.list = setNames(as.list(REFERENCE_COUNTRIES), paste0("cut_", REFERENCE_COUNTRIES))) |>
  tab_header(
    title = "Candidate illiberality thresholds",
    subtitle = sprintf(
      "Calibrated on populism against PopuList (%s), applied to anti-pluralism. N = %d scored RDD elections at w%d.",
      thr$coverage, nrow(d), SUMMARY_WINDOW
    )
  ) |>
  tab_source_note(paste(
    "Within-country percentiles rank a party-year among every V-Party party-year in its own country (all years);",
    "the implied raw cut is that country's anti-pluralism quantile at the percentile.",
    sprintf(
      "Front National: v2xpa_antiplural = %.3f in %d, the %.1fth percentile of French party-years.",
      fn$v2xpa_antiplural, fn$year, 100 * fn$pct_ctry
    ),
    "Party-years above: share of all scored V-Party party-years called illiberal.",
    "Treated: an ERT autocratization episode begins in the window."
  )) |>
  apply_table_style()
gtsave(gt_summary, file.path(out_dir, "cutpoint_summary.html"))

# ---- implied cut by country -------------------------------------------------

rdd_n <- d |> count(country_text_id, name = "rdd_elections")
by_country <- vp |>
  group_by(country_text_id) |>
  summarise(
    party_years = n(),
    median_score = median(v2xpa_antiplural),
    cut_ctry_accuracy = if (n() >= COUNTRY_PCT_MIN_N) {
      unname(quantile(v2xpa_antiplural, thr$threshold_pct_ctry, type = 1))
    } else {
      NA_real_
    },
    cut_ctry_youden = if (n() >= COUNTRY_PCT_MIN_N) {
      unname(quantile(v2xpa_antiplural, thr$threshold_pct_ctry_youden, type = 1))
    } else {
      NA_real_
    },
    .groups = "drop"
  ) |>
  left_join(rdd_n, by = "country_text_id") |>
  mutate(rdd_elections = coalesce(rdd_elections, 0L)) |>
  arrange(desc(rdd_elections), country_text_id)
write_csv(by_country, file.path(out_dir, "implied_cut_by_country.csv"))

gt_country <- by_country |>
  filter(rdd_elections > 0) |>
  gt() |>
  cols_label(
    country_text_id = "Country", party_years = "V-Party party-years",
    median_score = "Median anti-pluralism",
    cut_ctry_accuracy = sprintf("Cut, accuracy (%.1fth pct)", 100 * thr$threshold_pct_ctry),
    cut_ctry_youden = sprintf("Cut, Youden (%.1fth pct)", 100 * thr$threshold_pct_ctry_youden),
    rdd_elections = sprintf("RDD elections (w%d)", SUMMARY_WINDOW)
  ) |>
  fmt_number(columns = c(median_score, starts_with("cut_")), decimals = 3) |>
  sub_missing(missing_text = "--") |>
  tab_header(
    title = "The raw anti-pluralism score each within-country cut implies",
    subtitle = sprintf(
      "Countries in the RDD sample. The raw PopuList cut is %.3f everywhere, for comparison.",
      thr$threshold_abs
    )
  ) |>
  apply_table_style()
gtsave(gt_country, file.path(out_dir, "implied_cut_by_country.html"))

print(summary_tbl |> select(rule, cut, party_years_above, elections_pair, treated_pair, fn_clears, starts_with("cut_")), width = Inf)
message("Cutpoint summary written to ", out_dir)
