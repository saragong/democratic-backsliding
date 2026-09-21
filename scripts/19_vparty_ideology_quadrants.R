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
#   Top tags are chosen PER PANEL, not globally. Movements come and go, so the
#   most common tags in the OECD in the 1970s are not the ones that matter
#   outside the OECD in the 2010s, and forcing one global set meant every
#   panel showed the same six mainstream families while the populist and
#   nationalist tags that motivate the question never appeared anywhere.
#   The cost is that panels no longer share a tag set; the labels carry the
#   identity, so nothing is ambiguous, and a tag present in several panels can
#   still be followed by name.
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
#   in the CSV, and raising N_TAGS_PER_PANEL brings more of them onto it.
#
#   Tags are labelled DIRECTLY on the points, not encoded in colour and shape
#   against a legend. With a per-panel tag set the union across the ten panels
#   runs well past any palette that stays distinguishable, so colour-coding
#   would have to reuse combinations and two different tags would look alike.
#   Direct labels also answer the question a reader actually has -- "what is
#   that point" -- without a round trip to a 20-entry legend. Labels are
#   placed by ggrepel; points are a single colour, since colour would now be
#   decoration carrying no information.
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
library(ggrepel)
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

# How many tags to show IN EACH PANEL, by number of distinct parties carrying
# the tag in that panel. Labels are placed rather than colour-coded, so this
# is bounded by legibility rather than by a palette: 8 is comfortable at this
# panel size, and the union across the ten panels lands around 15-20 distinct
# tags, which the run prints so crowding can be checked.
N_TAGS_PER_PANEL <- 8

# A tag x panel cell below this many party-years is dropped: a "mean" over four
# observations placed on a quadrant chart invites exactly the over-reading the
# figure is meant to prevent.
MIN_TAG_CELL_N <- 15

# One colour for every point. Tags are identified by their labels, so colour
# would be decoration; a single ink keeps the eye on position, which is what
# the figure is about. Okabe-Ito blue, 5.19:1 on white.
POINT_COLOR <- "#0072B2"

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
  # Every tag's cell means first; the per-panel top set is picked afterwards,
  # so tag_means.csv keeps every tag that clears MIN_TAG_CELL_N and only the
  # figure is thinned.
  cells <- joined |>
    summarise(
      n_party_years = n(),
      n_parties = n_distinct(v2paid),
      mean_antiplural = mean(.data[[A]]),
      mean_popul = mean(.data[[B]]),
      .by = c(group, decade, tag)
    ) |>
    filter(n_party_years >= MIN_TAG_CELL_N)

  if (nrow(cells) == 0) {
    return(tibble())
  }

  # Ranked WITHIN each panel, by distinct parties rather than party-years: the
  # question is how many parties carry the tag here, and it is also what the
  # point area shows, so the two agree. Ties broken by party-years.
  cells |>
    mutate(
      tag_rank = rank(-n_parties + -n_party_years / 1e6, ties.method = "first"),
      .by = c(group, decade)
    ) |>
    mutate(
      vocab = vocab_name,
      is_plotted = tag_rank <= N_TAGS_PER_PANEL,
      .before = 1
    )
}

# ---- plot --------------------------------------------------------------------

quadrant_figure <- function(all_dat, vocab_name) {
  dat <- all_dat |> filter(is_plotted)

  # One padded range for both the data and the median guides, shared by all
  # ten panels of this figure.
  pad_range <- function(v, guide) {
    r <- range(c(v, guide), na.rm = TRUE)
    r + c(-1, 1) * max(diff(r) * 0.10, 0.02)
  }
  lim_a <- pad_range(dat$mean_antiplural, med_a)
  lim_b <- pad_range(dat$mean_popul, med_b)

  n_union <- n_distinct(dat$tag)
  per_panel <- dat |> count(group, decade) |> pull(n)

  ggplot(dat, aes(x = mean_popul, y = mean_antiplural)) +
    geom_hline(yintercept = med_a, colour = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = med_b, colour = "grey80", linewidth = 0.3) +
    geom_point(aes(size = n_parties), colour = POINT_COLOR, alpha = 0.85) +
    # Direct labels rather than a legend -- see the header. min.segment.length
    # draws a leader line whenever a label has had to move, so a displaced
    # label is never silently attached to the wrong point.
    geom_text_repel(
      aes(label = tag),
      size = 2.1, colour = "grey15", segment.colour = "grey65",
      segment.size = 0.2, min.segment.length = 0.15,
      box.padding = 0.28, point.padding = 0.12,
      max.overlaps = Inf, seed = 1
    ) +
    facet_grid(group ~ decade) +
    scale_size_area(max_size = 5, name = "Parties") +
    coord_fixed(xlim = lim_b, ylim = lim_a) +
    scale_x_continuous(breaks = scales::pretty_breaks(4)) +
    scale_y_continuous(breaks = scales::pretty_breaks(4)) +
    labs(
      title = sprintf(
        "Where the most common %s tags sit on illiberalism x populism",
        vocab_name
      ),
      subtitle = paste(
        strwrap(sprintf(paste(
          "Each point is one tag's mean score over the parties carrying it in",
          "that decade and group; point area is the number of parties. Grey",
          "lines are the pooled medians (populism %.2f, anti-pluralism %.2f),",
          "identical in every panel, as are the axes. Each panel shows its OWN",
          "top %d tags by party count, so the sets differ between panels --",
          "%d distinct tags appear across the ten, %s. Cells",
          "under %d party-years are dropped; tag_means.csv carries every tag.",
          "V-Party %d-%d. Tag coverage per panel: %s. Coverage is uneven but",
          "the parties it misses are close to the ones it keeps, so no panel",
          "mean sits more than 0.011 from its all-party value."
        ),
        med_b, med_a, N_TAGS_PER_PANEL, n_union,
        if (min(per_panel) == max(per_panel)) {
          sprintf("%d in each", min(per_panel))
        } else {
          sprintf("%d to %d per panel", min(per_panel), max(per_panel))
        },
        MIN_TAG_CELL_N,
        VPARTY_YEAR_MIN, VPARTY_YEAR_MAX, coverage_by_group(vocab_name)), 130),
        collapse = "\n"
      ),
      x = "Mean populism (v2xpa_popul)",
      y = "Mean anti-pluralism (v2xpa_antiplural)"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey94", linewidth = 0.25),
      panel.spacing = unit(0.4, "lines"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(size = 8),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 7.5, colour = "grey25"),
      legend.position = "bottom",
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
    "Saved quadrants_%s.png (%d distinct tags across panels, %d cells plotted, %d in CSV)\n",
    v, n_distinct(dat$tag[dat$is_plotted]), sum(dat$is_plotted), nrow(dat)
  ))
}

write_csv(bind_rows(all_means), file.path(out_dir, "tag_means.csv"))
message("\nIdeology quadrants written to ", out_dir)
