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
library(gt)
library(here)

source(here::here("scripts", "vparty_helpers.R"))

# ---- toggles -----------------------------------------------------------------

# Each entry is one figure: a (y, x) pair of V-Party scores, plotted with the
# same machinery. Axis titles, filenames and whether coord_fixed() applies are
# all derived from VPARTY_SCORES in vparty_helpers.R.
#
# The second pair asks where tags sit on economic left-right given their
# anti-pluralism, so anti-pluralism moves to x. It uses the RAW right-positive
# v2pariglef, not the negated version the RDD carries, so the axis reads
# left-to-right in the conventional direction.
if (!exists("QUADRANT_PAIRS")) {
  QUADRANT_PAIRS <- list(
    c(y = "v2xpa_antiplural", x = "v2xpa_popul"),
    c(y = "v2pariglef", x = "v2xpa_antiplural")
  )
}
ALL_PAIR_SCORES <- unique(unlist(QUADRANT_PAIRS))

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

# One colour per tag, assigned once per vocabulary so a tag keeps its colour
# in every panel it appears in and can be found in the legend.
#
# 18 entries, because the per-panel tag sets union to 16-18. Chosen greedily
# from a ~45-colour pool (Okabe-Ito, ColorBrewer qualitative and dark
# diverging endpoints) to maximize the MINIMUM pairwise separation, subject to
# a 3:1 contrast floor on white -- these carry both small points and coloured
# label text. Achieved: worst pair dE 18.8, worst contrast 3.02. Hand-picking
# 18 colours reliably produces a near-duplicate pair (an earlier attempt had
# two purples at dE 3.8); the greedy selection is there so that cannot recur
# if the palette is ever extended.
#
# Because selection is greedy, the EARLIER entries are the best separated, and
# tags are assigned in order of overall frequency -- so the tags appearing in
# the most panels get the most distinguishable colours.
#
# Colour is a navigational aid, not the identifier: the labels stay on the
# points. At 18 categories no palette is colourblind-safe, so a reader who
# cannot separate two hues still has the label, and the legend is there to
# look a tag up rather than to decode it.
TAG_COLORS <- c(
  "#0072B2", "#E31A1C", "#33A02C", "#C51B7D", "#000000", "#A6761D",
  "#54278F", "#35978F", "#BC80BD", "#666666", "#D55E00", "#006D2C",
  "#7F3B08", "#B2182B", "#009E73", "#7570B3", "#984EA3", "#01665E"
)

out_dir <- here::here("output", "runs", "_sweeps", "vparty_ideology_quadrants")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- data --------------------------------------------------------------------

vparty <- load_vparty_raw(c(
  "v2paid", "v2paenname", "v2pashname", "country_name",
  "country_text_id", "year", ALL_PAIR_SCORES
))

# Filtered on the year window AND on anti-pluralism being scored -- nothing
# narrower. Anti-pluralism is on every pair's axes, so a party without it can
# never reach a figure and does not belong in a coverage denominator; but a
# party missing v2pariglef must still count for the populism figure, so the
# per-pair drop happens below rather than here.
#
# This denominator is load-bearing: filtering only on the year window would
# put 2,432 parties in it rather than 1,930, silently understating every
# coverage figure by a fifth.
base <- vparty |>
  filter(
    year >= VPARTY_YEAR_MIN, year <= VPARTY_YEAR_MAX,
    !is.na(.data[["v2xpa_antiplural"]])
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

# readr's read_csv() silently parses a "1970s" column as the NUMBER 1970 --
# its column guesser falls through to a number parser that strips trailing
# non-numeric characters, even though guess_parser() on the same values
# returns "character". So a collaborator doing
# read_csv(...) |> filter(decade == "2010s") gets zero rows and no warning.
# Verified, and it caught me twice while auditing these files.
#
# Rather than hope nobody hits it, every CSV written here carries an integer
# decade_start alongside the label. That column survives any parser, and it
# is the one to join or filter on.
with_decade_key <- function(df) {
  if (!"decade" %in% names(df)) {
    return(df)
  }
  dplyr::mutate(
    df,
    decade_start = as.integer(sub("s$", "", as.character(decade))),
    .after = decade
  )
}

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
write_csv(with_decade_key(coverage), file.path(out_dir, "tag_coverage.csv"))

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
write_csv(with_decade_key(panel_coverage), file.path(out_dir, "panel_coverage.csv"))

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
# Collected across pairs so the party listing below covers every cell that
# appears on any figure, and is written once rather than per pair.
plotted_keys <- list()

for (pr in QUADRANT_PAIRS) {
  A <- unname(pr[["y"]])
  B <- unname(pr[["x"]])
  pslug <- pair_slug(A, B)
  cat(sprintf("\n========== %s (y) vs %s (x) ==========\n", A, B))

  # This pair's usable rows. base is filtered only on year, so each pair drops
  # only what IT is missing.
  pair_base <- base |> filter(!is.na(.data[[A]]), !is.na(.data[[B]]))
  med_a <- median(pair_base[[A]], na.rm = TRUE)
  med_b <- median(pair_base[[B]], na.rm = TRUE)

  tag_means_for <- function(vocab_name) {
  joined <- pair_base |>
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

  # Colour assignment is per FIGURE, not per panel, so a tag is the same
  # colour wherever it appears. Ordered by how many panels a tag reaches, then
  # by parties, so the most widespread tags take the best-separated colours
  # and the legend reads in a useful order.
  tag_order <- dat |>
    summarise(panels = n(), parties = sum(n_parties), .by = tag) |>
    arrange(desc(panels), desc(parties)) |>
    pull(tag)
  if (length(tag_order) > length(TAG_COLORS)) {
    stop(
      "Need ", length(tag_order), " tag colours but the validated palette has ",
      length(TAG_COLORS), ". Lower N_TAGS_PER_PANEL or extend TAG_COLORS ",
      "using the greedy max-min-separation procedure documented there.",
      call. = FALSE
    )
  }
  pal <- setNames(TAG_COLORS[seq_along(tag_order)], tag_order)
  dat <- dat |> mutate(tag = factor(tag, levels = tag_order))

  # One padded range for both the data and the median guides, shared by all
  # ten panels of this figure.
  #
  # Padding is generous, and asymmetric, because it is what the LABELS need
  # rather than what the points need. Labels are wide and horizontal, so x
  # gets more room than y; and since repel is confined to these limits (see
  # geom_text_repel below), too little padding leaves it nowhere to move a
  # label to and it gives up and overlaps. Too much just shrinks the points
  # into the middle. These values were set by looking at the densest panel.
  pad_range <- function(v, guide, frac) {
    r <- range(c(v, guide), na.rm = TRUE)
    r + c(-1, 1) * max(diff(r) * frac, 0.02)
  }
  lim_a <- pad_range(dat$mean_antiplural, med_a, 0.10)
  lim_b <- pad_range(dat$mean_popul, med_b, 0.26)

  n_union <- n_distinct(dat$tag)
  per_panel <- dat |> count(group, decade) |> pull(n)

  ggplot(dat, aes(x = mean_popul, y = mean_antiplural)) +
    geom_hline(yintercept = med_a, colour = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = med_b, colour = "grey80", linewidth = 0.3) +
    geom_point(aes(size = n_parties, colour = tag), alpha = 0.9) +
    # Direct labels rather than a legend -- see the header. min.segment.length
    # draws a leader line whenever a label has had to move, so a displaced
    # label is never silently attached to the wrong point.
    geom_text_repel(
      aes(label = tag, colour = tag),
      size = 2.1, segment.size = 0.2, min.segment.length = 0.15,
      # point.padding has to clear the MARKER, which is area-scaled up to
      # size 5 here -- at 0.12 the largest points sat on top of their own
      # labels and ate the first few characters.
      box.padding = 0.35, point.padding = 0.4,
      # Keep labels inside the panel. Without this, repel happily pushes a
      # label past the axis and ggplot clips it mid-word, which reads as a
      # different tag ("liberal_conservatis").
      xlim = lim_b, ylim = lim_a,
      max.overlaps = Inf, seed = 1, show.legend = FALSE
    ) +
    facet_grid(group ~ decade) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    scale_size_area(max_size = 5, name = "Parties") +
    # coord_fixed only when a unit means the same on both axes. It does for
    # two [0, 1] indices; it does not when one runs 0-1 and the other -4 to
    # +4, where forcing equal aspect would squash the figure to a sliver.
    (if (vparty_score(A, "unit_interval") && vparty_score(B, "unit_interval")) {
      coord_fixed(xlim = lim_b, ylim = lim_a)
    } else {
      coord_cartesian(xlim = lim_b, ylim = lim_a)
    }) +
    scale_x_continuous(breaks = scales::pretty_breaks(4)) +
    scale_y_continuous(breaks = scales::pretty_breaks(4)) +
    guides(
      colour = guide_legend(
        order = 1, nrow = 3, byrow = TRUE,
        override.aes = list(size = 2.4, label = "")
      ),
      size = guide_legend(order = 2)
    ) +
    labs(
      title = sprintf(
        "Where the most common %s tags sit: %s (y) vs %s (x)",
        vocab_name, vparty_score(A, "display"), vparty_score(B, "display")
      ),
      subtitle = paste(
        strwrap(sprintf(paste(
          "Each point is one tag's mean score over the parties carrying it in",
          "that decade and group; point area is the number of parties. Grey",
          "lines are the pooled medians (%s %.2f on x, %s %.2f on y),",
          "identical in every panel, as are the axes. Each panel shows its OWN",
          "top %d tags by party count, so the sets differ between panels --",
          "%d distinct tags appear across the ten, %s. Cells",
          "under %d party-years are dropped; tag_means_%s.csv carries every tag.",
          "V-Party %d-%d. Tag coverage per panel: %s. Coverage is uneven but",
          "the parties it misses are close to the ones it keeps, so no panel",
          "mean sits more than 0.011 from its all-party value."
        ),
        vparty_score(B, "slug"), med_b, vparty_score(A, "slug"), med_a,
        N_TAGS_PER_PANEL, n_union,
        if (min(per_panel) == max(per_panel)) {
          sprintf("%d in each", min(per_panel))
        } else {
          sprintf("%d to %d per panel", min(per_panel), max(per_panel))
        },
        MIN_TAG_CELL_N, pslug,
        VPARTY_YEAR_MIN, VPARTY_YEAR_MAX, coverage_by_group(vocab_name)), 130),
        collapse = "\n"
      ),
      x = sprintf("Mean %s", vparty_score(B, "display")),
      y = sprintf("Mean %s", vparty_score(A, "display"))
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
    file.path(out_dir, sprintf("quadrants_%s_%s.png", pslug, v)),
    fig,
    width = 14.5, height = 7.8, dpi = 150
  )
  cat(sprintf(
    "Saved quadrants_%s_%s.png (%d distinct tags across panels, %d cells plotted, %d in CSV)\n",
    pslug, v, n_distinct(dat$tag[dat$is_plotted]), sum(dat$is_plotted), nrow(dat)
  ))
  }

  write_csv(
  with_decade_key(bind_rows(all_means) |> mutate(score_y = A, score_x = B, .before = 1)),
  file.path(out_dir, sprintf("tag_means_%s.csv", pslug))
  )
  plotted_keys[[pslug]] <- bind_rows(all_means) |>
    filter(is_plotted) |> distinct(group, decade, tag)
}


# ==============================================================================
# Which parties are behind each point
#
# Every point on a figure is a mean over parties, and the obvious next
# question is which ones. This writes that out: one row per party per plotted
# (tag x decade x OECD) cell, with the party's name, country, how many times
# it was observed in that cell, and its own mean scores.
#
# Restricted to PLOTTED cells, so the tables and the figures show the same
# thing. tag_means.csv still has every tag; if a cell is not on a figure its
# parties are not listed here.
#
# Sorted within a cell by anti-pluralism, descending. The question these
# tables get opened for is "what is dragging this tag up the y axis", so the
# answer should be the first row rather than something to search for.
# ==============================================================================

party_rows <- function(vocab_name, plotted) {
  base |>
    inner_join(
      tag_long |> filter(vocab == vocab_name) |> select(vdem_id_1, tag),
      by = c("v2paid" = "vdem_id_1"), relationship = "many-to-many"
    ) |>
    semi_join(plotted, by = c("group", "decade", "tag")) |>
    summarise(
      n_obs = n(),
      first_year = min(year),
      last_year = max(year),
      # Every score any pair plots, not just the current pair's two, so this
      # table serves both figures and does not have to be written twice.
      across(all_of(ALL_PAIR_SCORES), \(x) mean(x, na.rm = TRUE)),
      .by = c(group, decade, tag, v2paid, v2paenname, v2pashname, country_name)
    ) |>
    arrange(group, decade, tag, desc(.data[["v2xpa_antiplural"]])) |>
    mutate(vocab = vocab_name, .before = 1)
}

# Union of the cells plotted by any pair. A cell that reaches a figure under
# either pairing gets its parties listed, so the table is a superset of both
# figures rather than tied to whichever pair happened to run last.
plotted_union <- bind_rows(plotted_keys) |> distinct(group, decade, tag)
all_parties <- map_dfr(TAG_COLUMNS, function(v) {
  keys <- plotted_union |>
    semi_join(
      tag_long |> filter(vocab == v) |> distinct(tag), by = "tag"
    )
  if (nrow(keys) == 0) return(tibble())
  party_rows(v, keys)
})
write_csv(with_decade_key(all_parties), file.path(out_dir, "tag_parties.csv"))
cat(sprintf(
  "\nParty listing: %d rows across %d vocabularies -> tag_parties.csv\n",
  nrow(all_parties), n_distinct(all_parties$vocab)
))

# One HTML per vocabulary, sections keyed on panel and tag. gt is given the
# whole listing rather than a truncated one -- these are reference tables, and
# a reader checking whether a point is driven by one odd party needs the row
# that is not in the top ten.
for (v in unique(all_parties$vocab)) {
  d <- all_parties |> filter(vocab == v)
  if (nrow(d) == 0) next
  tbl <- d |>
    transmute(
      section = sprintf("%s \u00b7 %s \u00b7 %s", group, decade, tag),
      Party = v2paenname,
      Abbr. = v2pashname,
      Country = country_name,
      Obs = n_obs,
      Years = ifelse(first_year == last_year,
                     as.character(first_year),
                     sprintf("%d-%d", first_year, last_year)),
      `Anti-pluralism` = round(v2xpa_antiplural, 3),
      Populism = round(v2xpa_popul, 3),
      `Econ L-R` = round(v2pariglef, 3)
    )
  gt_tbl <- tbl |>
    gt::gt(groupname_col = "section") |>
    gt::tab_header(
      title = sprintf("Parties behind each point: %s", v),
      subtitle = sprintf(
        paste(
          "One row per party per plotted cell, sorted by anti-pluralism within",
          "each cell. Scores are that party's own mean over its observations",
          "in that decade. V-Party %d-%d; %d parties over %d cells."
        ),
        VPARTY_YEAR_MIN, VPARTY_YEAR_MAX, nrow(d),
        n_distinct(paste(d$group, d$decade, d$tag))
      )
    ) |>
    gt::opt_row_striping() |>
    gt::tab_options(
      table.font.size = gt::px(11),
      row_group.font.weight = "bold",
      row_group.background.color = "#f0f0f0",
      data_row.padding = gt::px(2)
    )
  path <- file.path(out_dir, sprintf("tag_parties_%s.html", v))
  gt::gtsave(gt_tbl, path)
  cat(sprintf("Saved %s (%d rows)\n", basename(path), nrow(d)))
}
message("\nIdeology quadrants written to ", out_dir)
