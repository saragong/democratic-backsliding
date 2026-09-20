# ==============================================================================
# Ad hoc: does anti-pluralism track populism the same way everywhere, always?
#
# The pooled correlation between V-Party's two illiberalism indices (see
# vparty_raw_correlation.R) hides two splits that move it a lot: OECD vs the
# rest, and decade. Rich democracies show a positive and rising association --
# by the 2010s a populist party is also, reliably, an anti-pluralist one --
# while outside the OECD the association is negative early on and roughly zero
# after 1990.
#
# Three series per panel, each a Pearson r on the same party-years, differing
# only in what gets demeaned out first:
#   Raw correlation      no demeaning. Across all party-years in the cell.
#   Election fixed effects  demeaned within election (country x year), i.e. the
#                        WITHIN-election correlation. Asks whether the party
#                        that is more populist than its rivals in a given race
#                        is also the more anti-pluralist one, netting out any
#                        common shift in how a whole election's field is coded.
#   Country x party fixed effects  demeaned within party. Asks whether a party
#                        that drifts more populist than its OWN average also
#                        drifts more anti-pluralist -- a within-party,
#                        over-time correlation, with every cross-party
#                        difference netted out. (v2paid is already unique per
#                        country, so country x party and party are the same
#                        cell here; the country term is harmless.)
#
# Note the FE ladder is NOT nested: election FE holds the race fixed and
# compares parties, party FE holds the party fixed and compares races. They
# answer different questions and can easily point opposite ways.
#
# Reads the raw V-Party CSV straight out of the zip -- no election spine, no
# top-2 filter, no crosswalk -- so this is a property of the indices, not of
# this project's sample construction. The read, the OECD list and the coverage
# bounds come from scripts/vparty_helpers.R.
#
#   Rscript --no-init-file scripts/adhoc/vparty_corr_by_decade.R
# ==============================================================================

library(tidyverse)
library(here)

source(here::here("scripts", "vparty_helpers.R"))

# ---- toggles -----------------------------------------------------------------
A <- "v2xpa_antiplural" # "illiberalism" axis
B <- "v2xpa_popul" # "populism" axis

# Coverage bounds live in vparty_helpers.R with the reason (pre-1970 V-Party
# is thin and lopsided; 2019 is the last year in v2).
YEAR_MIN <- VPARTY_YEAR_MIN
YEAR_MAX <- VPARTY_YEAR_MAX
MIN_CELL_N <- 30 # drop a decade x group cell thinner than this

# OECD membership (and the reason it is applied to the whole period rather
# than by accession date) is in vparty_helpers.R, shared with the other
# raw-V-Party scripts.
OECD <- OECD_ISO3

# Which parties enter the party-FE series (this restriction does NOT touch the
# raw or election-FE series, which stay on the full sample). Mirrors the toggle
# of the same name in vparty_variance_by_decade.R.
#   "balanced" keep only parties observed in EVERY election their country held
#             in that decade, so the within-party correlation reflects genuine
#             movement rather than which parties happened to enter or exit
#             mid-decade.
#   "spanning" looser: present at the decade's first AND last election in that
#             country, but allowed to miss one in between.
#   "all"     no restriction.
PARTY_PANEL <- "balanced"
stopifnot(PARTY_PANEL %in% c("balanced", "spanning", "all"))

SERIES_COLORS <- c(
  "Raw correlation" = "#2c6e91",
  "Election fixed effects" = "#c0562f",
  "Country x party fixed effects" = "#3f8f6f"
)
SERIES_SHAPES <- c(16, 15, 17)

# Party FE only has bite where a party appears more than once in the cell.
#   TRUE  demean within party x decade, so the fixed effect is nested inside
#         the cell being reported (as election FE already is) and the series is
#         a clean within-cell statistic -- but a party seen once in a decade
#         contributes two zeros and drops out.
#   FALSE demean within party over the whole 1970-2019 sample, which keeps
#         those parties by scoring them against their long-run mean -- at the
#         cost of the demeaning using data from outside the decade plotted.
PARTY_FE_WITHIN_DECADE <- TRUE

# ---- data --------------------------------------------------------------------
vparty <- load_vparty_raw()

stopifnot(all(OECD %in% unique(vparty$country_text_id)))

# One historical_date per country-year in V-Party, so country x year IS the
# election -- no need to carry the date around as a separate key.
d <- vparty |>
  transmute(
    country = country_text_id,
    party = v2paid,
    year,
    a = .data[[A]],
    b = .data[[B]]
  ) |>
  filter(!is.na(a), !is.na(b), year >= YEAR_MIN, year <= YEAR_MAX) |>
  mutate(
    decade = paste0(year %/% 10 * 10, "s"),
    group = if_else(country %in% OECD, "OECD", "Non-OECD")
  )

cat(sprintf(
  "%s party-years, %s countries, %d-%d (%s OECD / %s non-OECD)\n\n",
  format(nrow(d), big.mark = ","),
  format(n_distinct(d$country), big.mark = ","),
  min(d$year),
  max(d$year),
  format(sum(d$group == "OECD"), big.mark = ","),
  format(sum(d$group == "Non-OECD"), big.mark = ",")
))

# Balance flag: compare the election years a party is observed in, inside a
# decade, against the election years its country held in that same decade. The
# flag is constant within a party-decade, so under PARTY_FE_WITHIN_DECADE it
# makes no difference whether the subset is taken before or after demeaning.
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
  "Party-FE panel restriction %s: keeps %s of %s party-years (%.0f%%)\n\n",
  sQuote(PARTY_PANEL),
  format(sum(d$party_keep), big.mark = ","),
  format(nrow(d), big.mark = ","),
  100 * mean(d$party_keep)
))

# Demeaning. A group seen once contributes two zeros and so drops out of the
# correlation on its own -- no filtering needed for that -- but it does mean
# the three series rest on different effective samples, which n_eff_* records.
# Demeaning happens on whatever subset is passed in, so the party series is
# centred on the means of the RETAINED party-years only.
within_cor <- function(x, y, g, keep = TRUE) {
  x <- x[keep]
  y <- y[keep]
  g <- g[keep]
  if (length(x) < 3) {
    return(NA_real_)
  }
  rx <- x - ave(x, g)
  ry <- y - ave(y, g)
  if (sd(rx) == 0 || sd(ry) == 0) {
    return(NA_real_)
  }
  cor(rx, ry)
}

party_group <- if (PARTY_FE_WITHIN_DECADE) {
  with(d, paste(country, party, decade))
} else {
  with(d, paste(country, party))
}

tab <- d |>
  # assigned ungrouped: both are full-length vectors keyed to the rows of d
  mutate(election_group = paste(country, year), party_group = party_group) |>
  mutate(n_election = n(), .by = election_group) |>
  # a party-year only counts toward its party group if it is RETAINED
  mutate(n_party = sum(party_keep), .by = party_group) |>
  summarise(
    n_party_years = n(),
    n_elections = n_distinct(election_group),
    n_countries = n_distinct(country),
    n_parties = n_distinct(party),
    # party-years that actually carry within-group variation
    n_eff_election = sum(n_election > 1),
    n_eff_party = sum(party_keep & n_party > 1),
    n_party_years_kept = sum(party_keep),
    `Raw correlation` = cor(a, b),
    `Election fixed effects` = within_cor(a, b, election_group),
    `Country x party fixed effects` = within_cor(a, b, party_group, party_keep),
    .by = c(group, decade)
  ) |>
  filter(n_party_years >= MIN_CELL_N) |>
  arrange(desc(group), decade)

print(
  as.data.frame(tab |> mutate(across(where(is.numeric), \(x) round(x, 4)))),
  row.names = FALSE
)

long <- tab |>
  pivot_longer(
    all_of(names(SERIES_COLORS)),
    names_to = "series",
    values_to = "r"
  ) |>
  mutate(
    series = factor(series, levels = names(SERIES_COLORS)),
    # OECD panel on the left, as in the original.
    group = factor(group, levels = c("OECD", "Non-OECD"))
  )

# ---- plot --------------------------------------------------------------------
# Labels are placed by hand (no ggrepel here). With three series several cells
# hold near-ties -- non-OECD 2010s spans 0.015 across all three -- so anchoring
# each label to its own point overlaps. Instead, stack a cell's labels in value
# order under a minimum gap: the top value sits just above its point, each next
# one below the previous, sliding away from its own point only as far as the
# tie forces. Reading order still matches line order, and colour carries the
# mapping.
LABEL_PAD <- 0.045 # label-to-point clearance, in correlation units
LABEL_SEP <- 0.062 # minimum gap between two labels in the same cell

stack_labels <- function(r, pad = LABEL_PAD, sep = LABEL_SEP) {
  pos <- numeric(length(r))
  prev <- Inf
  for (i in order(r, decreasing = TRUE)) {
    cand <- if (is.infinite(prev)) r[i] + pad else min(r[i] - pad, prev - sep)
    # A label also has to clear every POINT MARKER in the cell, not just the
    # labels above it: in OECD 1980s the election-FE label lands exactly on the
    # raw point otherwise. Drop below whichever marker it collides with and
    # re-check, since that drop can land it on another.
    # min(r[hit]) - pad sits exactly `pad` from that marker, which floating
    # point can round back into the collision set -- hence the tolerance and
    # the forced strict decrease, without which this repeat never exits.
    repeat {
      hit <- which(abs(cand - r) < pad - 1e-9)
      if (!length(hit)) {
        break
      }
      cand <- min(min(r[hit]) - pad, cand - pad / 2)
    }
    pos[i] <- cand
    prev <- cand
  }
  pos
}

long <- long |>
  mutate(label_y = stack_labels(r), .by = c(group, decade))

fig <- ggplot(long, aes(decade, r, colour = series, group = series)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(shape = series), size = 2.2) +
  geom_text(
    aes(y = label_y, label = sprintf("%.3f", r)),
    size = 2.9,
    show.legend = FALSE
  ) +
  facet_wrap(~group) +
  scale_colour_manual(values = SERIES_COLORS, name = NULL) +
  scale_shape_manual(values = SERIES_SHAPES, name = NULL) +
  scale_y_continuous(
    limits = c(-0.65, 0.65),
    breaks = round(seq(-0.6, 0.6, 0.2), 1),
    labels = \(x) sprintf("%.1f", x)
  ) +
  labs(
    title = "Populism and illiberalism by decade",
    x = "Decade",
    y = "Pearson correlation",
    caption = switch(
      PARTY_PANEL,
      balanced = "Country x party series restricted to parties observed at every election their country held in the decade; the other two series use all party-years.",
      spanning = "Country x party series restricted to parties present at the decade's first and last election in their country; the other two series use all party-years.",
      all = NULL
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 12),
    plot.caption = element_text(size = 6.5, colour = "grey35", hjust = 0)
  )

out <- here::here("output", "adhoc")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

png_path <- file.path(out, "vparty_corr_by_decade.png")
csv_path <- file.path(out, "vparty_corr_by_decade.csv")
ggsave(png_path, fig, width = 10, height = 4.6, dpi = 150)
write_csv(tab, csv_path)
cat(sprintf("\nSaved %s\nSaved %s\n", png_path, csv_path))
