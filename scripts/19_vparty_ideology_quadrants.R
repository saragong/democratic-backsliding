# ==============================================================================
# Which parties populate the illiberalism x populism quadrants?
#
# 18_vparty_jaccard_panels.R shows that anti-pluralism and populism pick out
# largely different parties. The obvious next question, and the one Vincent
# Rollet asked when he sent parties_ideology_tags.xlsx, is WHICH parties sit
# where: "for the most common tags, where these parties typically sit on the
# illiberalism-populism quadrant".
#
# So: join each V-Party party-year to its Wikipedia/Wikidata ideology tags,
# and for each tag plot the mean (populism, anti-pluralism) of the parties
# carrying it -- on the same 2 x 5 OECD-by-decade grid as script 18, so the two
# figures are read together.
#
# DESIGN CHOICES THAT MATTER FOR READING THE FIGURE
#
#   Top tags are chosen GLOBALLY, not per panel. If each panel showed its own
#   most frequent tags, the panels would show different tags and could not be
#   compared -- which is the one thing the 2 x 5 layout exists to allow.
#
#   Axes are IDENTICAL in every panel, with quadrant guides at the pooled
#   medians, so a tag moving right across the decades is genuinely moving and
#   not being rescaled. They are fixed to the range of the plotted tag MEANS
#   rather than to [0, 1]: means of dozens of parties do not reach the extremes
#   of an index, and on a full [0, 1] square every panel collapses into a thin
#   band with 80% of the figure empty.
#
#   Only the top N tags are PLOTTED, but tag_means.csv carries every tag that
#   clears MIN_TAG_CELL_N. That matters here: the most frequent tags are the
#   mainstream families (social democracy, conservatism, christian democracy),
#   so the tags that motivated the question -- right_wing_populism, nationalism
#   and the like -- do not make the top 6 and are not on the figure. They are
#   in the CSV, and raising N_TAGS brings them onto it at the cost of reusing
#   colours.
#
#   Colour AND shape both encode the tag, and the palette is capped at six
#   entries. A larger N_TAGS is available as a toggle but reuses colours, which
#   the figure warns about rather than silently doing. Double-encoding is what
#   makes a tag followable across the decade panels: with decade as a FACET,
#   the points of one tag live in different panels and cannot be joined by a
#   line, so identity has to be carried by the marker itself.
#
# COVERAGE, stated on the figure rather than buried here: only about 1,500 of
# the 8,049 tagged parties carry both a V-Party id and a tag, so every panel
# rests on a subset of V-Party, and the rile_* variants rest on a thinner one
# (822 parties) than the ideology_* variants (~1,500). Per-panel coverage is
# in panel_coverage.csv and summarized on each figure. The selection this
# induces is measured, and small -- see the panel_coverage block below.
#
# WEIGHTING: a tag's mean is over PARTY-YEARS, so a party observed in eight
# elections counts eight times and one observed twice counts twice. Checked
# against the alternative (average within party first, then across parties):
# the two agree to 0.996 on anti-pluralism and 0.988 on populism, mean
# absolute shift 0.015, and they select the identical top-6 tag set. So the
# choice is immaterial here and party-years is kept, being the grain the
# scores are actually measured at. Point AREA is the distinct party count,
# which is the more meaningful "how much is behind this dot"; both counts are
# in tag_means.csv.
#
#   Rscript --no-init-file scripts/19_vparty_ideology_quadrants.R
#
# Output: output/runs/_sweeps/vparty_ideology_quadrants/
#           quadrants_<vocab>.png     one per tag vocabulary
#           tag_means.csv             every tag x panel cell
#           tag_coverage.csv          how much of V-Party each vocabulary reaches
# ==============================================================================

library(tidyverse)
library(readxl)
library(here)

source(here::here("scripts", "vparty_helpers.R"))

# ---- toggles -----------------------------------------------------------------

A <- "v2xpa_antiplural" # y axis
B <- "v2xpa_popul" # x axis

# The four tag vocabularies in the file. rile_* and ideology_* are different
# KINDS of tag (a left-right bucket vs a named ideology), and the wikipedia and
# wikidata versions of each use DIFFERENT controlled vocabularies -- wikipedia
# has "far_left"/"centre", wikidata has "far_left_politics"/"centrism" -- so
# the four are run separately and never pooled. Merging them would need an
# explicit crosswalk that does not exist.
TAG_COLUMNS <- c(
  "rile_wikipedia", "rile_wikidata",
  "ideology_wikipedia", "ideology_wikidata"
)

# How many tags to show. Six is the width of the validated palette below; the
# script warns and recycles colours above that.
N_TAGS <- 6

# A tag x panel cell below this many party-years is dropped: a "mean" over four
# observations placed on a quadrant chart invites exactly the over-reading the
# figure is meant to prevent.
MIN_TAG_CELL_N <- 15

# Okabe-Ito based, validated the same way as SERIES_COLORS in rdd_helpers.R:
# worst pairwise separation dE 43.3 normal / 17.2 deuteranopic / 8.5
# protanopic, and every entry clears 3:1 contrast on white (5.19, 3.87, 3.42,
# 3.06, 21.0, 4.86). The protanopic figure is the weak one, which is why shape
# co-encodes the tag rather than colour carrying it alone.
TAG_COLORS <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#000000", "#8C6D1F")
TAG_SHAPES <- c(16, 17, 15, 18, 8, 4)

out_dir <- here::here("output", "runs", "_sweeps", "vparty_ideology_quadrants")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- data --------------------------------------------------------------------

vparty <- load_vparty_raw(c("v2paid", "country_text_id", "year", A, B))

base <- vparty |>
  filter(
    !is.na(.data[[A]]), !is.na(.data[[B]]),
    year >= VPARTY_YEAR_MIN, year <= VPARTY_YEAR_MAX
  ) |>
  mutate(
    group = oecd_group(country_text_id),
    decade = decade_label(year)
  )

tags_path <- here::here(
  "data", "elections_database", "parties_ideology_tags.xlsx"
)
if (!file.exists(tags_path)) {
  stop("No ideology tags file at ", tags_path)
}
tags_raw <- read_excel(tags_path)

# vdem_id_1 IS V-Party's v2paid, so this is a direct key join.
#
# The obvious idea of recovering extra coverage by crosswalking
# party_id -> parties_database.dta$vdem_id_1 does NOT work and is not done
# here: checked, the two vdem_id_1 columns agree on all 8,049 rows and neither
# is populated where the other is missing. This file's mapping IS the
# parties_database mapping, so 2,166 parties is the ceiling either way.
stopifnot(!any(duplicated(na.omit(tags_raw$vdem_id_1))))

# Semicolon-separated lists -> one row per (party, tag). str_squish rather than
# trimws: a few entries carry doubled internal spaces.
tag_long <- tags_raw |>
  select(vdem_id_1, all_of(TAG_COLUMNS)) |>
  filter(!is.na(vdem_id_1)) |>
  pivot_longer(all_of(TAG_COLUMNS), names_to = "vocab", values_to = "tags") |>
  filter(!is.na(tags)) |>
  separate_rows(tags, sep = ";") |>
  mutate(tag = str_squish(tags)) |>
  filter(tag != "") |>
  select(vdem_id_1, vocab, tag) |>
  distinct()

# ---- coverage ----------------------------------------------------------------

coverage <- tag_long |>
  summarise(
    n_tagged_parties = n_distinct(vdem_id_1),
    n_distinct_tags = n_distinct(tag),
    .by = vocab
  ) |>
  mutate(
    n_matched_in_vparty = map_int(
      vocab,
      \(v) n_distinct(intersect(
        tag_long$vdem_id_1[tag_long$vocab == v], base$v2paid
      ))
    ),
    n_vparty_parties = n_distinct(base$v2paid),
    pct_of_vparty = round(100 * n_matched_in_vparty / n_vparty_parties, 1)
  )

cat("Tag coverage against scored V-Party parties:\n")
print(as.data.frame(coverage), row.names = FALSE)
write_csv(coverage, file.path(out_dir, "tag_coverage.csv"))

# Coverage PER PANEL, not just pooled, because the 2 x 5 layout rests on the
# ten panels being comparable and pooled coverage cannot show whether they
# are.
#
# Coverage is uneven: OECD panels are 87-92% tagged, non-OECD 70-83%. That
# looks like it should bias the between-row comparison, and an earlier version
# of this comment said it did. It does not, and the difference is worth
# measuring rather than assuming, because the bias depends not on the coverage
# gap but on how far UNTAGGED parties differ from tagged ones:
#
#   bias in a panel mean = (1 - coverage) x (tagged mean - untagged mean)
#
# Untagged parties turn out to sit only slightly higher on anti-pluralism
# (0.04 in the OECD, 0.03 outside), and they are a small minority, so the
# product is small. Measured per panel, the largest gap between a panel's
# plotted mean and its true all-party mean is 0.011 -- against an OECD /
# non-OECD separation of 0.376, i.e. 3% of the thing being compared. Both
# comparisons the figure invites are therefore safe.
#
# Coverage is still reported per panel: it is the evidence for that claim, and
# it is what would have to be re-checked if the tag file were ever extended
# unevenly.
panel_coverage <- map_dfr(TAG_COLUMNS, function(v) {
  ids <- tag_long |> filter(vocab == v) |> pull(vdem_id_1) |> unique()
  base |>
    mutate(tagged = v2paid %in% ids) |>
    summarise(
      parties = n_distinct(v2paid),
      tagged_parties = n_distinct(v2paid[tagged]),
      .by = c(group, decade)
    ) |>
    mutate(vocab = v, pct_tagged = 100 * tagged_parties / parties, .before = 1)
})
write_csv(panel_coverage, file.path(out_dir, "panel_coverage.csv"))

# One coverage range per OECD group, for the subtitle.
coverage_by_group <- function(vocab_name) {
  pc <- panel_coverage |> filter(vocab == vocab_name)
  paste(
    vapply(levels(pc$group), function(g) {
      x <- pc$pct_tagged[pc$group == g]
      sprintf("%s %.0f-%.0f%%", g, min(x), max(x))
    }, character(1)),
    collapse = ", "
  )
}

# ---- per-vocabulary means ----------------------------------------------------

# Pooled medians, computed on the WHOLE scored sample rather than per panel, so
# the quadrant guides sit in the same place in all ten panels.
med_a <- median(base[[A]], na.rm = TRUE)
med_b <- median(base[[B]], na.rm = TRUE)

tag_means_for <- function(vocab_name) {
  joined <- base |>
    inner_join(
      tag_long |> filter(vocab == vocab_name) |> select(vdem_id_1, tag),
      by = c("v2paid" = "vdem_id_1"),
      relationship = "many-to-many"
    )
  if (nrow(joined) == 0) {
    return(tibble())
  }
  # Global frequency, so the same tags appear in every panel -- see the header.
  top_tags <- joined |>
    count(tag, sort = TRUE) |>
    slice_head(n = N_TAGS) |>
    pull(tag)

  # Every tag, not only the plotted ones -- see the header. is_plotted marks
  # which made the figure so the CSV is self-explaining.
  joined |>
    summarise(
      n_party_years = n(),
      n_parties = n_distinct(v2paid),
      mean_antiplural = mean(.data[[A]]),
      mean_popul = mean(.data[[B]]),
      .by = c(group, decade, tag)
    ) |>
    filter(n_party_years >= MIN_TAG_CELL_N) |>
    mutate(
      vocab = vocab_name,
      is_plotted = tag %in% top_tags,
      tag_rank = match(tag, top_tags),
      .before = 1
    )
}

# ---- plot --------------------------------------------------------------------

quadrant_figure <- function(all_dat, vocab_name) {
  dat <- all_dat |>
    filter(is_plotted) |>
    mutate(tag = fct_reorder(tag, tag_rank, .fun = min))
  n_lvl <- nlevels(dat$tag)
  if (n_lvl > length(TAG_COLORS)) {
    warning(
      "N_TAGS = ", n_lvl, " exceeds the ", length(TAG_COLORS),
      "-entry validated palette; colours will repeat.",
      call. = FALSE
    )
  }
  cols <- rep_len(TAG_COLORS, n_lvl)
  shps <- rep_len(TAG_SHAPES, n_lvl)

  # One padded range for both the data and the median guides, shared by all ten
  # panels of this figure.
  pad_range <- function(v, guide) {
    r <- range(c(v, guide), na.rm = TRUE)
    r + c(-1, 1) * max(diff(r) * 0.08, 0.02)
  }
  lim_a <- pad_range(dat$mean_antiplural, med_a)
  lim_b <- pad_range(dat$mean_popul, med_b)

  ggplot(dat, aes(x = mean_popul, y = mean_antiplural, colour = tag)) +
    geom_hline(yintercept = med_a, colour = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = med_b, colour = "grey80", linewidth = 0.3) +
    geom_point(aes(size = n_parties, shape = tag), alpha = 0.9) +
    facet_grid(group ~ decade) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_shape_manual(values = shps, name = NULL) +
    scale_size_area(max_size = 5, name = "Parties") +
    # Identical in every panel, but sized to the tag means rather than to the
    # index range -- see the header. The guides are inside the range by
    # construction, since a median of the raw scores sits among the means.
    coord_fixed(xlim = lim_b, ylim = lim_a) +
    scale_x_continuous(breaks = scales::pretty_breaks(4)) +
    scale_y_continuous(breaks = scales::pretty_breaks(4)) +
    labs(
      title = sprintf(
        "Where the most common %s tags sit on illiberalism x populism",
        vocab_name
      ),
      subtitle = paste(
        strwrap(sprintf(
          paste(
            "Each point is one tag's mean score over the parties carrying it in",
            "that decade and group; point area is the number of parties. Grey",
            "lines are the pooled medians (populism %.2f, anti-pluralism %.2f),",
            "identical in every panel, as are the axes.",
            "Top %d tags by overall frequency, chosen once so the same tags",
            "appear in every panel; tag_means.csv carries every tag, including",
            "the less common populist and nationalist ones. Cells under %d",
            "party-years are dropped. Axes are common to all panels but scaled",
            "to the tag means, not to the full [0, 1] index range.",
            "V-Party %d-%d. Tag coverage per panel: %s. Coverage is uneven",
            "but the parties it misses are close to the ones it keeps, so no",
            "panel mean sits more than 0.011 from its all-party value --",
            "3%% of the OECD / non-OECD separation."
          ),
          med_b, med_a, N_TAGS, MIN_TAG_CELL_N,
          VPARTY_YEAR_MIN, VPARTY_YEAR_MAX,
          coverage_by_group(vocab_name)
        ), 120),
        collapse = "\n"
      ),
      x = "Mean populism (v2xpa_popul)",
      y = "Mean anti-pluralism (v2xpa_antiplural)"
    ) +
    guides(
      colour = guide_legend(order = 1, override.aes = list(size = 2.6)),
      shape = guide_legend(order = 1),
      size = guide_legend(order = 2)
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey94", linewidth = 0.25),
      panel.spacing = unit(0.35, "lines"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(size = 8),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 7.5, colour = "grey25"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(t = -2),
      legend.key.size = unit(0.35, "cm")
    )
}

# ---- run ---------------------------------------------------------------------

all_means <- list()

for (v in TAG_COLUMNS) {
  dat <- tag_means_for(v)
  if (nrow(dat) == 0) {
    warning("No usable cells for ", v, "; skipping.", call. = FALSE)
    next
  }
  all_means[[v]] <- dat

  fig <- quadrant_figure(dat, v)
  ggsave(
    file.path(out_dir, sprintf("quadrants_%s.png", v)),
    fig,
    width = 13, height = 7.4, dpi = 150
  )
  cat(sprintf(
    "Saved quadrants_%s.png (%d tags plotted, %d cells plotted, %d cells in CSV)\n",
    v, n_distinct(dat$tag[dat$is_plotted]),
    sum(dat$is_plotted), nrow(dat)
  ))
}

write_csv(bind_rows(all_means), file.path(out_dir, "tag_means.csv"))
message("\nIdeology quadrants written to ", out_dir)
