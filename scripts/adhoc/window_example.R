# ==============================================================================
# Ad hoc: inspect one election's outcome window, end to end.
#
# Pulls a single election out of an RDD build, re-derives its Y_* outcomes from
# the raw country-year panel by hand, checks them against what the build
# actually stored, and plots the underlying series with the window marked.
#
# Point of it: every Y_* outcome is produced by one function
# (window_change() in 11_build_rdd_data.R) with two different branches, and the
# units differ by branch. This makes the arithmetic visible for a case you know.
#
#   Rscript --no-init-file scripts/adhoc/window_example.R
#   # or, interactively:
#   COUNTRY <- "HUN"; YEAR <- 2010; source("scripts/adhoc/window_example.R")
# ==============================================================================

library(tidyverse)
library(patchwork)
library(here)
library(haven)

source(here::here("scripts", "rdd_helpers.R"))

# ---- what to look at ---------------------------------------------------------
if (!exists("COUNTRY")) {
  COUNTRY <- "PER"
}
if (!exists("YEAR")) {
  YEAR <- 2016
}
if (!exists("N")) {
  N <- 5
} # window length; must match the build below
if (!exists("BUILD")) {
  BUILD <- "rdd_v2xpa_antiplural_w5_exclyr.rds"
}

# Which instrument this build was scored on -- it decides which of the top 2
# counts as "the illiberal side", and therefore the sign of running_var.
ACTIVE_INSTRUMENT <- sub("^rdd_(.*)_w[0-9]+(_exclyr)?\\.rds$", "\\1", BUILD)

panel <- readRDS(here::here("data", "combined_panel.rds"))
build <- readRDS(here::here("data", "rdd_build", BUILD))
# The party-level companion 11_build_rdd_data.R writes alongside every build:
# two rows per election, one per top-2 finisher, carrying every instrument's
# score. Party NAMES aren't in it (only Party Facts ids), so they come from
# parties_database.dta below.
parties <- readRDS(here::here(
  "data",
  "rdd_build",
  sub("\\.rds$", "_parties.rds", BUILD)
))

pre_year <- YEAR - 1 # baseline: the last year before the window
post_year <- YEAR + N # end of the window

# A country-year can hold BOTH a presidential and a parliamentary election
# (Peru 2016, for one), and they are separate rows with different running
# variables. Set TYPE to disambiguate; without it, listing them and stopping is
# better than silently taking whichever comes first.
if (!exists("TYPE")) {
  TYPE <- NA_character_
}

row <- build |> filter(country_text_id == COUNTRY, election_year == YEAR)
if (!is.na(TYPE)) {
  row <- row |> filter(election_type == TYPE)
}
if (nrow(row) == 0) {
  stop(
    "No election for ",
    COUNTRY,
    " ",
    YEAR,
    if (!is.na(TYPE)) paste0(" (", TYPE, ")") else "",
    " in ",
    BUILD,
    ". Available: ",
    paste(
      build |>
        filter(country_text_id == COUNTRY) |>
        arrange(election_year) |>
        mutate(
          lab = sprintf("%d %s", election_year, substr(election_type, 1, 4))
        ) |>
        pull(lab),
      collapse = ", "
    )
  )
}
if (nrow(row) > 1) {
  stop(
    COUNTRY,
    " ",
    YEAR,
    " has ",
    nrow(row),
    " elections in this build (",
    paste(row$election_type, collapse = ", "),
    "). Set TYPE to pick one, e.g. TYPE <- \"",
    row$election_type[1],
    "\"."
  )
}

# A couple of years of padding either side, so you can see what the window cuts off.
win <- panel |>
  filter(
    country_text_id == COUNTRY,
    year >= pre_year - 3,
    year <= post_year + 3
  ) |>
  arrange(year)

cat(sprintf(
  "\n%s (%s) %s election %d | window: %d (baseline) -> %d\n",
  row$country_name,
  COUNTRY,
  row$election_type,
  YEAR,
  pre_year,
  post_year
))
cat(sprintf(
  "  running_var = %+.2f pp  (%s side led)   illiberal_score = %.3f vs %.3f\n",
  row$running_var,
  if (row$running_var > 0) "more-illiberal" else "less-illiberal",
  row$illiberal_score,
  row$other_score
))
cat(sprintf(
  "  backsliding_Nyr = %d | union = %d | polyarchy_decline = %+.4f | declined = %d\n",
  row$backsliding_Nyr,
  row$backsliding_union_Nyr,
  row$polyarchy_decline,
  row$polyarchy_declined
))

# ---- who the top 2 actually were --------------------------------------------
# party_id is a Party Facts code (e.g. "PF1204"); the readable name lives in
# parties_database.dta. A party_id can have several genealogy_rank rows (the
# party renamed or re-registered), so collapse the distinct names rather than
# silently picking one -- this is a display table, not the matching logic.
party_names <- read_dta(here::here(
  "data",
  "elections_database",
  "parties_database.dta"
)) |>
  select(party_id, name, short_name) |>
  summarise(
    party_name = paste(unique(na.omit(name)), collapse = " / "),
    short_name = paste(unique(na.omit(short_name)), collapse = " / "),
    .by = party_id
  )

top2 <- parties |>
  filter(election_id == row$election_id) |>
  left_join(party_names, by = "party_id") |>
  arrange(desc(final_share)) |>
  mutate(
    place = if_else(row_number() == 1, "WON", "lost"),
    # "the illiberal side" is the higher scorer on the ACTIVE instrument, which
    # is what defines the sign of running_var -- and it is NOT always the winner.
    side = if_else(
      .data[[ACTIVE_INSTRUMENT]] == max(.data[[ACTIVE_INSTRUMENT]]),
      "more illiberal",
      "less illiberal"
    )
  )

cat(sprintf(
  "\n  Top 2 by %s share, scored on %s:\n",
  if (row$election_type == "parliamentary") "seat" else "final-round vote",
  ACTIVE_INSTRUMENT
))
print(
  as.data.frame(
    top2 |>
      transmute(
        place,
        side,
        share = round(final_share, 2),
        party = if_else(
          is.na(party_name) | party_name == "",
          party_id,
          sprintf("%s (%s)", party_name, party_id)
        ),
        candidate = coalesce(candidate, "--"),
        across(all_of(ILLIBERALISM_VARS), \(x) round(x, 3))
      )
  ),
  row.names = FALSE
)
cat(sprintf(
  "  -> the more-illiberal side %s, so running_var = %+.2f\n\n",
  if (row$running_var > 0) "WON" else "LOST",
  row$running_var
))

# ---- re-derive the outcomes by hand ------------------------------------------
at <- function(v, y) {
  panel[[v]][match(paste(COUNTRY, y), paste(panel$country_text_id, panel$year))]
}

# Branch 1: endpoint difference. Used for every LEVEL variable -- so the units
# are whatever the level was in (percentage points of GDP, index points, ...).
endpoint_diff <- function(v) at(v, post_year) - at(v, pre_year)

# Branch 2: compound the annual RATES over the window, then log. Used only for
# inflation, which has no price-level column to difference. Units become log
# points of the price LEVEL -- NOT a change in the inflation rate.
compounded_log <- function(v) {
  rates <- vapply(
    seq(pre_year + 1, post_year),
    function(y) at(v, y),
    numeric(1)
  )
  if (any(is.na(rates))) {
    return(NA_real_)
  }
  log(prod(1 + rates / 100))
}

checks <- tribble(
  ~outcome            , ~units                      , ~hand                              ,
  "Y_gdp_growth"      , "log points (cumulative)"   , endpoint_diff("ln_gdp_pc")         ,
  "Y_inflation"       , "log points of price LEVEL" , compounded_log("cpi_inflation")    ,
  "Y_unemployment"    , "pp of labour force"        , endpoint_diff("unemployment_rate") ,
  "Y_trade_pct_gdp"   , "pp of GDP"                 , endpoint_diff("trade_pct_gdp")     ,
  "Y_debt"            , "pp of GDP (stock)"         , endpoint_diff("debt_pct_gdp")      ,
  "Y_deficit"         , "pp of GDP (annual flow)"   , endpoint_diff("deficit_pct_gdp")   ,
  "Y_gini_disp"       , "Gini points"               , endpoint_diff("gini_disp")         ,
  "Y_checks_balances" , "index points (0-1)"        , endpoint_diff("checks_balances")   ,
) |>
  mutate(
    stored = vapply(outcome, function(o) as.numeric(row[[o]]), numeric(1)),
    match = if_else(
      is.na(hand) & is.na(stored),
      "both NA",
      if_else(abs(hand - stored) < 1e-9, "ok", "MISMATCH")
    )
  )

cat("Re-derived by hand vs what the build stored:\n")
print(
  as.data.frame(checks |> mutate(across(c(hand, stored), \(x) round(x, 4)))),
  row.names = FALSE
)

# ---- plot the underlying series ----------------------------------------------
# cpi_inflation is a rate, so it gets an extra panel showing the cumulative
# price index the outcome is actually built from.
win <- win |>
  mutate(
    price_index = cumprod(if_else(
      year > pre_year & !is.na(cpi_inflation),
      1 + cpi_inflation / 100,
      1
    ))
  )

SERIES <- tribble(
  ~var                , ~label                              ,
  "v2x_polyarchy"     , "V-Dem polyarchy (index)"           ,
  "ln_gdp_pc"         , "log GDP per capita (PWT)"          ,
  "cpi_inflation"     , "CPI inflation (% per year)"        ,
  "price_index"       , "Price level (index, baseline = 1)" ,
  "debt_pct_gdp"      , "Debt (% of GDP)"                   ,
  "deficit_pct_gdp"   , "Deficit (% of GDP, annual)"        ,
  "unemployment_rate" , "Unemployment (%)"                  ,
  "checks_balances"   , "Checks & balances (0-1)"           ,
)

one_panel <- function(v, label) {
  dat <- win |> select(year, y = all_of(v)) |> filter(!is.na(y))
  if (nrow(dat) == 0) {
    return(NULL)
  }
  endpoints <- dat |> filter(year %in% c(pre_year, post_year))
  ggplot(dat, aes(year, y)) +
    # the window the outcome is measured over
    annotate(
      "rect",
      xmin = pre_year,
      xmax = post_year,
      ymin = -Inf,
      ymax = Inf,
      fill = "#0072B2",
      alpha = 0.07
    ) +
    geom_vline(
      xintercept = YEAR,
      linetype = "dashed",
      colour = "#D55E00",
      linewidth = 0.4
    ) +
    geom_line(colour = "grey35", linewidth = 0.5) +
    geom_point(size = 1, colour = "grey35") +
    # the two values the endpoint difference actually uses
    geom_point(data = endpoints, size = 2.4, colour = "#0072B2") +
    scale_x_continuous(breaks = seq(pre_year - 2, post_year + 2, 2)) +
    labs(title = label, x = NULL, y = NULL) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 8, face = "bold")
    )
}

# ---- who-fought-whom line for the subtitle ----------------------------------
# Plain sprintf in the default font. Column alignment isn't worth chasing here:
# the separators carry the structure, and forcing alignment needs either a
# monospace family (renders as Courier, looks wrong beside the panels) or a
# table grob (more machinery than one line of text deserves).
share_word <- if (row$election_type == "parliamentary") "seats" else "votes"

who <- top2 |>
  mutate(
    party_disp = if_else(
      is.na(party_name) | party_name == "",
      party_id,
      party_name
    ),
    label = if_else(
      is.na(candidate),
      party_disp,
      sprintf("%s, %s", candidate, party_disp)
    ),
    outcome_word = if_else(row_number() == 1, "WON", "LOST")
  )

subtitle_lines <- c(
  sprintf(
    "%s: %s \u00b7 %.1f%% of %s \u00b7 %s %.3f \u00b7 %s",
    who$outcome_word,
    who$label,
    who$final_share,
    share_word,
    ACTIVE_INSTRUMENT,
    who[[ACTIVE_INSTRUMENT]],
    who$side
  ),
  sprintf(
    "Margin %+.2f pp, so the more-illiberal side %s.",
    row$running_var,
    if (row$running_var > 0) "won" else "lost",
    ACTIVE_INSTRUMENT
  )
)

panels <- compact(map2(SERIES$var, SERIES$label, one_panel))
fig <- wrap_plots(panels, ncol = 2) +
  plot_annotation(
    title = sprintf(
      "%s %d %s election: the %d-year outcome window",
      row$country_name,
      YEAR,
      row$election_type,
      N
    ),
    subtitle = paste(subtitle_lines, collapse = "\n"),
    caption = sprintf(
      "Shaded = %d to %d, the span every outcome is measured over. Dashed = election year. Blue dots = the two endpoints an endpoint-difference outcome uses.\nInflation is the exception: it compounds the annual rates inside the window (see the price-level panel), it does not difference the rate.",
      pre_year,
      post_year
    ),
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(
        size = 8.5,
        colour = "grey20",
        lineheight = 1.4,
        margin = margin(t = 3, b = 8)
      ),
      plot.caption = element_text(size = 6.5, colour = "grey35", hjust = 0)
    )
  )

out <- here::here("output", "adhoc")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
path <- file.path(
  out,
  sprintf(
    "window_example_%s_%d_%s_w%d.png",
    COUNTRY,
    YEAR,
    substr(row$election_type, 1, 4),
    N
  )
)
ggsave(path, fig, width = 9, height = 9, dpi = 150)
cat(sprintf("\nSaved %s\n", path))
