# ==============================================================================
# Ad hoc: where does the variation in V-Party's illiberalism indices live?
#
# Companion to vparty_corr_by_decade.R, which plots the CORRELATION between
# anti-pluralism and populism under the same OECD x decade split and the same
# two fixed effects. This one drops the second index and asks a prior question
# about each index on its own: how much does it actually move, and along which
# margin?
#
# Two variance series per panel, both residual variances from a within
# transformation, so each is the part of the index's variation that survives
# one demeaning:
#   Election FE (across parties)  demeaned within election (country x year).
#                        The spread of the index ACROSS THE PARTIES CONTESTING
#                        ONE ELECTION -- how ideologically far apart a typical
#                        field is.
#   Country x party FE (across years)  demeaned within party. The spread of the
#                        index ACROSS YEARS FOR A GIVEN PARTY -- how much a
#                        party's own coding moves over time.
#
# Writes one figure per index (anti-pluralism, populism) plus a combined CSV.
#
# Reads the raw V-Party CSV straight out of the zip -- no election spine, no
# top-2 filter, no crosswalk. Self-contained on purpose (as the sibling adhoc
# scripts are), so the OECD list below is duplicated from vparty_corr_by_decade.R.
#
#   Rscript --no-init-file scripts/adhoc/vparty_variance_by_decade.R
# ==============================================================================

library(tidyverse)
library(here)

# ---- toggles -----------------------------------------------------------------
INDICES <- c(
  v2xpa_antiplural = "anti-pluralism",
  v2xpa_popul = "populism"
)

YEAR_MIN <- 1970 # V-Party reaches back to 1900, but coverage before 1970 is
YEAR_MAX <- 2019 # thin and lopsided; 2019 is the last year in v2.
MIN_CELL_N <- 30 # drop a decade x group cell thinner than this

# Current (2026) OECD membership, applied to the whole period -- so the split is
# "countries that ended up rich democracies", not membership as of each
# election. Kept identical to vparty_corr_by_decade.R.
OECD <- c(
  "AUS", "AUT", "BEL", "CAN", "CHL", "COL", "CRI", "CZE", "DNK", "EST",
  "FIN", "FRA", "DEU", "GRC", "HUN", "ISL", "IRL", "ISR", "ITA", "JPN",
  "KOR", "LVA", "LTU", "LUX", "MEX", "NLD", "NZL", "NOR", "POL", "PRT",
  "SVK", "SVN", "ESP", "SWE", "CHE", "TUR", "GBR", "USA"
)

# Which years a party is demeaned against, for the party-FE series only.
#   TRUE  demean within party x decade, so the cell's variance is built only
#         from movement INSIDE that decade. Honest to the x axis, but thin:
#         a party averages ~2 elections per decade here, so each group is
#         close to a single first difference.
#   FALSE demean within party across all of 1970-2019, i.e. spread around the
#         party's long-run mean, then attribute each party-year to its decade.
#         Roughly 5-11 observations per party and an order of magnitude larger
#         variances -- but a cell then reflects movement from outside it.
PARTY_FE_WITHIN_DECADE <- TRUE

# Which parties enter the party-FE series (this restriction does NOT touch the
# election-FE series or total_var, which stay on the full sample).
#   "balanced" keep only parties observed in EVERY election their country held
#             in that decade, so the within-party variance reflects genuine
#             movement rather than which parties happened to enter or exit
#             mid-decade. Costs 20-43% of party-years per cell.
#   "spanning" looser: present at the decade's first AND last election in that
#             country, but allowed to miss one in between.
#   "all"     no restriction (what the first version of this script did).
# A country-decade holding only ONE election qualifies trivially under either
# restriction, but such parties are singleton groups that contribute nothing to
# the variance either way -- see the n - G note below.
PARTY_PANEL <- "balanced"
stopifnot(PARTY_PANEL %in% c("balanced", "spanning", "all"))

# The two within variances differ by 10-100x under the default toggle, so a
# linear axis flattens the party-FE series onto zero and hides its shape.
LOG_Y <- FALSE

SERIES_COLORS <- c(
  "Election FE (spread across parties)" = "#c0562f",
  "Country x party FE (spread across years)" = "#3f8f6f"
)
SERIES_SHAPES <- c(15, 17)

# ---- data --------------------------------------------------------------------
vparty <- read_csv(
  unz(
    here::here("data", "elections_database", "CPD_V-Party_CSV_v2.zip"),
    "CPD_V-Party_CSV_v2/V-Dem-CPD-Party-V2.csv"
  ),
  show_col_types = FALSE
)

stopifnot(all(OECD %in% unique(vparty$country_text_id)))

# One historical_date per country-year in V-Party, so country x year IS the
# election. v2paid is already unique per country, so country x party is just
# party; the country term is kept for readability and is harmless.
d <- vparty |>
  select(
    country = country_text_id,
    party = v2paid,
    year,
    all_of(names(INDICES))
  ) |>
  # Both indices are required, not just the one being plotted, so the two
  # figures and the correlation script all rest on the same 6,304 party-years.
  filter(if_all(all_of(names(INDICES)), \(x) !is.na(x))) |>
  filter(year >= YEAR_MIN, year <= YEAR_MAX) |>
  mutate(
    decade = paste0(year %/% 10 * 10, "s"),
    group = if_else(country %in% OECD, "OECD", "Non-OECD")
  )

# ---- within variance ---------------------------------------------------------
# The residual variance from a within transformation, SSR / (n - G), NOT var()
# of the residuals. Demeaning by group costs one degree of freedom per group,
# and group sizes differ a lot across these cells (~4-6 parties per election vs
# ~2 elections per party-decade), so dividing by n-1 would understate the
# party-FE series far more than the election-FE one and make the two series
# non-comparable. n - G also disposes of singleton groups correctly: a party
# seen once contributes 0 to SSR and 1 to both n and G, so it neither inflates
# nor deflates the estimate.
within_var <- function(x, g) {
  n <- length(x)
  n_groups <- n_distinct(g)
  if (n <= n_groups) {
    return(NA_real_)
  }
  sum((x - ave(x, g))^2) / (n - n_groups)
}

# Balance flag: compare the election years a party is observed in, inside a
# decade, against the election years its country held in that same decade.
d <- d |>
  mutate(country_years = list(sort(unique(year))), .by = c(country, decade)) |>
  mutate(
    party_years = list(sort(unique(year))),
    .by = c(country, party, decade)
  ) |>
  mutate(
    party_keep = switch(
      PARTY_PANEL,
      balanced = map2_lgl(party_years, country_years, identical),
      spanning = map2_lgl(party_years, country_years, \(p, c) {
        min(p) == min(c) && max(p) == max(c)
      }),
      all = TRUE
    )
  ) |>
  select(-country_years, -party_years)

cat(sprintf(
  "Party-FE panel restriction %s: keeps %s of %s party-years (%.0f%%)\n",
  sQuote(PARTY_PANEL),
  format(sum(d$party_keep), big.mark = ","),
  format(nrow(d), big.mark = ","),
  100 * mean(d$party_keep)
))

party_group <- if (PARTY_FE_WITHIN_DECADE) {
  with(d, paste(country, party, decade))
} else {
  with(d, paste(country, party))
}

tab <- d |>
  mutate(election_group = paste(country, year), party_group = party_group) |>
  pivot_longer(
    all_of(names(INDICES)),
    names_to = "index_var",
    values_to = "x"
  ) |>
  summarise(
    n_party_years = n(),
    n_elections = n_distinct(election_group),
    n_parties = n_distinct(party_group),
    parties_per_election = n() / n_distinct(election_group),
    # the party-FE series runs on the restricted subset, so its own n and mean
    # group size are reported separately rather than shared with the row above
    n_party_years_kept = sum(party_keep),
    years_per_party = sum(party_keep) /
      n_distinct(party_group[party_keep]),
    total_var = var(x),
    `Election FE (spread across parties)` = within_var(x, election_group),
    `Country x party FE (spread across years)` = within_var(
      x[party_keep],
      party_group[party_keep]
    ),
    .by = c(group, decade, index_var)
  ) |>
  filter(n_party_years >= MIN_CELL_N) |>
  mutate(index = unname(INDICES[index_var]), .after = index_var) |>
  arrange(index_var, desc(group), decade)

cat(sprintf(
  "%s party-years, %s countries, %d-%d | party FE demeaned within %s\n\n",
  format(n_distinct(paste(d$country, d$party, d$year)), big.mark = ","),
  format(n_distinct(d$country), big.mark = ","),
  min(d$year), max(d$year),
  if (PARTY_FE_WITHIN_DECADE) "party x decade" else "party (full period)"
))
print(
  as.data.frame(tab |> mutate(across(where(is.double), \(x) round(x, 4)))),
  row.names = FALSE
)

# ---- plot --------------------------------------------------------------------
# Label declutter, as in vparty_corr_by_decade.R: stack a cell's labels in value
# order under a minimum gap, and keep them clear of the point markers too.
# Distances are in whatever units the y axis is drawn in, so on a log axis the
# stacking happens in log space and is mapped back at the end.
LABEL_PAD <- 0.035 # label-to-point clearance, as a share of the axis range
LABEL_SEP <- 0.050 # minimum gap between two labels in the same cell

stack_labels <- function(v, span, pad = LABEL_PAD, sep = LABEL_SEP,
                         floor = -Inf) {
  pad <- pad * span
  sep <- sep * span
  pos <- numeric(length(v))
  prev <- Inf
  for (i in order(v, decreasing = TRUE)) {
    cand <- if (is.infinite(prev)) v[i] + pad else min(v[i] - pad, prev - sep)
    if (cand < floor) {
      # On a linear axis the smallest variances sit a hair above 0, so the
      # usual "put it below the point" would push the label off the bottom of
      # a panel that starts at zero. Flip it above the point instead.
      cand <- v[i] + pad
      repeat {
        hit <- which(abs(cand - v) < pad - 1e-9)
        if (!length(hit)) break
        cand <- max(max(v[hit]) + pad, cand + pad / 2)
      }
    } else {
      repeat {
        hit <- which(abs(cand - v) < pad - 1e-9)
        if (!length(hit)) break
        cand <- min(min(v[hit]) - pad, cand - pad / 2)
      }
    }
    pos[i] <- cand
    prev <- cand
  }
  pos
}

make_fig <- function(var_name) {
  label <- unname(INDICES[var_name])

  long <- tab |>
    filter(index_var == var_name) |>
    pivot_longer(
      all_of(names(SERIES_COLORS)),
      names_to = "series",
      values_to = "v"
    ) |>
    mutate(
      series = factor(series, levels = names(SERIES_COLORS)),
      group = factor(group, levels = c("OECD", "Non-OECD"))
    )

  scale_v <- if (LOG_Y) log10(long$v) else long$v
  span <- diff(range(scale_v))

  long <- long |>
    mutate(
      scale_v = if (LOG_Y) log10(v) else v,
      label_pos = stack_labels(
        scale_v,
        span,
        floor = if (LOG_Y) -Inf else 0
      ),
      label_y = if (LOG_Y) 10^label_pos else label_pos,
      .by = c(group, decade)
    )

  ggplot(long, aes(decade, v, colour = series, group = series)) +
    geom_line(linewidth = 0.7) +
    geom_point(aes(shape = series), size = 2.2) +
    geom_text(
      aes(y = label_y, label = sprintf("%.4f", v)),
      size = 2.7,
      show.legend = FALSE
    ) +
    facet_wrap(~group) +
    scale_colour_manual(values = SERIES_COLORS, name = NULL) +
    scale_shape_manual(values = SERIES_SHAPES, name = NULL) +
    (if (LOG_Y) {
      scale_y_log10(
        expand = expansion(mult = 0.12),
        labels = \(x) sprintf("%g", x)
      )
    } else {
      # A variance has a meaningful zero, so the linear axis starts there
      # rather than at the smallest observed value.
      scale_y_continuous(
        limits = c(0, NA),
        expand = expansion(mult = c(0.02, 0.14))
      )
    }) +
    labs(
      title = sprintf("Spread of %s by decade", label),
      subtitle = sprintf(
        "Within variance of V-Party %s (residual variance, SSR/(n-G))%s\n%s",
        var_name,
        if (LOG_Y) " · log scale" else "",
        switch(
          PARTY_PANEL,
          balanced = "Party series restricted to parties observed at every election their country held in the decade",
          spanning = "Party series restricted to parties present at the decade's first and last election in their country",
          all = "Party series uses all parties, including those observed only once in the decade"
        )
      ),
      x = "Decade",
      y = "Within variance"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 8, colour = "grey30")
    )
}

out <- here::here("output", "adhoc")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

for (v in names(INDICES)) {
  path <- file.path(out, sprintf("vparty_variance_%s_by_decade.png", v))
  ggsave(path, make_fig(v), width = 10, height = 4.8, dpi = 150)
  cat(sprintf("\nSaved %s", path))
}

csv_path <- file.path(out, "vparty_variance_by_decade.csv")
write_csv(tab, csv_path)
cat(sprintf("\nSaved %s\n", csv_path))
