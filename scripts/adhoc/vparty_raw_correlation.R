# ==============================================================================
# Ad hoc: anti-pluralism vs populism in the RAW V-Party dataset.
#
# Deliberately reads the V-Party CSV straight out of the zip -- no election
# spine, no top-2 filter, no Party Facts crosswalk, no nearest-year matching.
# Every party-year V-Dem coded, so the correlation is a property of the two
# indices themselves rather than of this project's sample construction.
#
# Compare against the party-level correlation on the analysis sample
# (16_instrument_overlap.R's pairwise_correlation.html), which is computed on
# top-2 finishers only.
#
#   Rscript --no-init-file scripts/adhoc/vparty_raw_correlation.R
# ==============================================================================

library(tidyverse)
library(patchwork)
library(here)

A <- "v2xpa_antiplural"
B <- "v2xpa_popul"

vparty <- read_csv(
  unz(
    here::here("data", "elections_database", "CPD_V-Party_CSV_v2.zip"),
    "CPD_V-Party_CSV_v2/V-Dem-CPD-Party-V2.csv"
  ),
  show_col_types = FALSE
)

cat(sprintf(
  "Raw V-Party: %s rows (party-years), %s parties, %s countries, %d-%d\n",
  format(nrow(vparty), big.mark = ","),
  format(n_distinct(vparty$v2paid), big.mark = ","),
  format(n_distinct(vparty$country_text_id), big.mark = ","),
  min(vparty$year), max(vparty$year)
))
cat(sprintf(
  "  %s non-missing: %s | %s non-missing: %s | both: %s\n\n",
  A, format(sum(!is.na(vparty[[A]])), big.mark = ","),
  B, format(sum(!is.na(vparty[[B]])), big.mark = ","),
  format(sum(!is.na(vparty[[A]]) & !is.na(vparty[[B]])), big.mark = ",")
))

d <- vparty |>
  select(
    v2paid, party = v2paenname, country = country_text_id, year,
    a = all_of(A), b = all_of(B)
  ) |>
  filter(!is.na(a), !is.na(b))

pear <- cor(d$a, d$b)
spear <- cor(d$a, d$b, method = "spearman")
cat(sprintf(
  "ALL %s party-years:  Pearson r = %.4f   Spearman rho = %.4f   R^2 = %.4f\n",
  format(nrow(d), big.mark = ","), pear, spear, pear^2
))

# The unit of observation matters: V-Party repeats a party across election
# years, so party-years are not independent. Collapsing to one row per party
# (its mean score) shows whether the correlation is driven by within-party
# repetition.
by_party <- d |>
  summarise(a = mean(a), b = mean(b), .by = c(v2paid, party))
cat(sprintf(
  "Collapsed to %s distinct parties (mean score): Pearson r = %.4f   Spearman rho = %.4f\n",
  format(nrow(by_party), big.mark = ","),
  cor(by_party$a, by_party$b),
  cor(by_party$a, by_party$b, method = "spearman")
))

# For contrast: the same pair on this project's analysis sample.
parties_path <- here::here(
  "data", "rdd_build", "rdd_v2xpa_antiplural_w5_exclyr_parties.rds"
)
if (file.exists(parties_path)) {
  p <- readRDS(parties_path) |> filter(!is.na(.data[[A]]), !is.na(.data[[B]]))
  cat(sprintf(
    "Analysis sample (top-2 finishers only, %s rows): Pearson r = %.4f   Spearman rho = %.4f\n",
    format(nrow(p), big.mark = ","),
    cor(p[[A]], p[[B]]),
    cor(p[[A]], p[[B]], method = "spearman")
  ))
}

# ---- scatter -----------------------------------------------------------------
# Two views of the same data because 11k points on a [0,1] x [0,1] square
# overplot badly: raw points with heavy transparency show the individual
# observations and the empty regions, the binned version shows where the mass
# actually is. Blue line = OLS; dashed grey = the 45-degree line the two indices
# would follow if they measured the same thing on the same scale.
base <- list(
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey55", linewidth = 0.4),
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "#0072B2", fill = "#0072B2", alpha = 0.2, linewidth = 0.7),
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)),
  labs(x = sprintf("Populism (%s)", B), y = sprintf("Anti-pluralism (%s)", A)),
  theme_bw(base_size = 9),
  theme(panel.grid.minor = element_blank())
)

p_pts <- ggplot(d, aes(b, a)) +
  geom_point(alpha = 0.10, size = 0.7, colour = "grey25") +
  base +
  labs(title = sprintf("Every party-year (n = %s)", format(nrow(d), big.mark = ",")))

p_bin <- ggplot(d, aes(b, a)) +
  geom_bin2d(bins = 40) +
  scale_fill_gradient(low = "#eaf2f8", high = "#08519c", name = "party-years") +
  base +
  labs(title = "Same data, binned by density")

fig <- (p_pts | p_bin) +
  plot_annotation(
    title = "V-Party anti-pluralism vs populism: raw, unfiltered",
    subtitle = sprintf(
      "All %s party-years with both indices coded · Pearson r = %.3f · Spearman rho = %.3f · R² = %.3f\nDashed line = 45° (what perfect agreement would look like); blue = OLS fit.",
      format(nrow(d), big.mark = ","), pear, spear, pear^2
    ),
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 8, colour = "grey25")
    )
  )

out <- here::here("output", "adhoc")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
path <- file.path(out, "vparty_raw_antiplural_vs_popul.png")
ggsave(path, fig, width = 11, height = 5.6, dpi = 150)
cat(sprintf("\nSaved %s\n", path))

write_csv(d, file.path(out, "vparty_raw_antiplural_vs_popul.csv"))
cat(sprintf("Saved %s\n", file.path(out, "vparty_raw_antiplural_vs_popul.csv")))
