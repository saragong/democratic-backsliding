# ==============================================================================
# Anti-pluralism vs populism over the FULL V-Party dataset, by decade and OECD.
#
# 16_instrument_overlap.R already draws this heatmap, but only for the 2,694
# top-2 finishers in this project's election spine. That sample is small, it is
# selected on having reached a top 2 AND on having a Party Facts -> V-Party
# match, and it pools 50 years. So when it shows the two indices disagreeing
# more than they agree (Pearson 0.115), it is fair to ask how much of that is
# the sample.
#
# This script asks the same question of the 11,898-party-year source dataset --
# no election spine, no top-2 filter, no crosswalk -- and splits it the two
# ways the correlation is known to move (Vincent Rollet's email, and
# adhoc/vparty_corr_by_decade.R): OECD vs the rest, and decade.
#
# The output is deliberately one 2 x 5 figure rather than ten: rows are
# OECD / non-OECD, columns are the 1970s to the 2010s, and every panel shares
# one fill scale, so the eye can read the change over time straight across.
# Ten separately-scaled heatmaps would make any comparison between them an
# artefact of their own colour ranges.
#
# Two universes, because the sample question cuts both ways:
#   all_parties  every scored V-Party party-year. What the indices look like.
#   ever_top2    only parties that reached a top 2 somewhere in the election
#                spine. Isolates the SELECTION from the rest of the spine's
#                construction: if the top-2 heatmap differs from all_parties
#                but matches this, the difference is who reaches a top 2, not
#                the crosswalk or the pooling.
#
#   Rscript --no-init-file scripts/18_vparty_jaccard_panels.R
#
# Output: output/runs/_sweeps/vparty_jaccard_panels/
#           jaccard_panels_<universe>.png
#           jaccard_cells_<universe>.csv     every cell of every panel
#           panel_summary.csv                N and diagonal mass per panel
# ==============================================================================

library(tidyverse)
library(here)

source(here::here("scripts", "vparty_helpers.R"))

# ---- toggles -----------------------------------------------------------------

# Each entry is one figure: a (y, x) pair of V-Party scores. Both are plotted
# with the SAME machinery; only the axes, the bin mode and the filenames
# differ, all derived from VPARTY_SCORES in vparty_helpers.R.
#
# The second pair puts anti-pluralism on x against economic left-right on y.
# Note the orientation differs from the first pair deliberately: in
# antiplural_vs_popul the question is "do these two indices agree", and
# anti-pluralism is the subject on y; in econlr_vs_antiplural the question is
# "where do parties sit on left-right given their anti-pluralism", so
# anti-pluralism becomes the x conditioner.
if (!exists("JACCARD_PANEL_PAIRS")) {
  JACCARD_PANEL_PAIRS <- list(
    c(y = "v2xpa_antiplural", x = "v2xpa_popul"),
    c(y = "v2pariglef", x = "v2xpa_antiplural")
  )
}

N_BINS <- 10

# A panel below this many party-years is reported but not drawn: a 10x10 grid
# over 30 observations is mostly empty cells and a few 1.0s, which reads as a
# strong signal and is noise.
MIN_CELL_N <- 100

# Where the build's top-2 parties come from, for the ever_top2 universe. Any
# window works -- the top-2 membership does not depend on the outcome window --
# so this just picks the one that exists.
TOP2_BUILD <- here::here(
  "data", "rdd_build", "rdd_v2xpa_antiplural_w5_exclyr_parties.rds"
)

out_dir <- here::here("output", "runs", "_sweeps", "vparty_jaccard_panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- data --------------------------------------------------------------------

for (pr in JACCARD_PANEL_PAIRS) {
  A <- unname(pr[["y"]])
  B <- unname(pr[["x"]])
  bin_mode <- pair_bin_mode(A, B)
  pslug <- pair_slug(A, B)
  cat(sprintf(
    "\n========== %s (y) vs %s (x)  [%s bins] ==========\n",
    A, B, bin_mode
  ))

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

  cat(sprintf(
    "V-Party %d-%d, both indices scored: %d party-years, %d parties, %d countries\n",
    VPARTY_YEAR_MIN, VPARTY_YEAR_MAX,
    nrow(base), n_distinct(base$v2paid), n_distinct(base$country_text_id)
  ))

  # The top-2 universe. _parties.rds carries vdem_id_1, which IS V-Party's
  # v2paid, so this is a direct key join -- no name matching.
  if (!file.exists(TOP2_BUILD)) {
    stop(
      "No party-level build at ", TOP2_BUILD,
      ". Run 11_build_rdd_data.R first (the ever_top2 universe needs it)."
    )
  }
  top2_ids <- readRDS(TOP2_BUILD)$vdem_id_1 |> unique()
  top2_ids <- top2_ids[!is.na(top2_ids)]
  cat(sprintf(
    "Parties reaching a top 2 at least once: %d (%d of them scored here)\n",
    length(top2_ids), n_distinct(base$v2paid[base$v2paid %in% top2_ids])
  ))

  UNIVERSES <- list(
    all_parties = list(
      label = "All V-Party parties",
      fn = function(df) df
    ),
    ever_top2 = list(
      label = "Parties that reached a top 2 in the election spine",
      fn = function(df) df |> filter(v2paid %in% top2_ids)
    )
  )

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

  # ---- cells -------------------------------------------------------------------

  # Breaks are computed ONCE over the whole pair, not per panel, so a cell
# means the same thing in every panel. Mode comes from pair_bin_mode().
  brk_a <- bin_breaks(base[[A]], A, bin_mode, N_BINS)
  brk_b <- bin_breaks(base[[B]], B, bin_mode, N_BINS)

  # Every (group, decade) combination, including any that turn out empty, so a
  # missing cell shows as an empty panel rather than silently collapsing the grid
  # and shifting the decades out of line between the two rows.
  cell_keys <- expand_grid(
    group = levels(base$group),
    decade = levels(base$decade)
  ) |>
    mutate(group = factor(group, levels = levels(base$group)),
           decade = factor(decade, levels = levels(base$decade)))

  build_cells <- function(df) {
    pmap_dfr(cell_keys, function(group, decade) {
      cell <- df |> filter(.data$group == !!group, .data$decade == !!decade)
      if (nrow(cell) == 0) {
        return(tibble())
      }
      jaccard_long(cell, A, B, brk_a, brk_b) |>
        mutate(group = group, decade = decade, .before = 1)
    })
  }

  # Diagonal mass: the mean Jaccard on the leading diagonal, over the mean off
  # it. One number per panel for "how concentrated is this heatmap", so the
  # visual reading can be checked against something.
  #
  # IT IS NOT A MEASURE OF ASSOCIATION, and the difference matters for reading
  # the figure. Jaccard is |both| / |either| on RAW counts, so a cell is bright
  # wherever the two marginals both pile up, whether or not the indices are
  # related. Across these ten panels diag_ratio correlates 0.78 with the share
  # of parties in the lowest anti-pluralism bin -- almost as strongly as it
  # correlates with the actual association (0.71 with Pearson). Two thirds of
  # OECD parties sit in that lowest bin, which is most of why the OECD (1,1)
  # corner glows in every decade.
  #
  # So the correlation is carried alongside, per panel, and shown on the panel
  # strip. It is the marginal-free statement, and the two genuinely disagree:
  # across OECD decades Pearson rises monotonically 0.22 -> 0.53 while
  # diag_ratio wanders 1.62, 1.77, 1.78, 1.14, 1.55. Read the correlation for
  # "are these indices related here"; read the heatmap for "where do the parties
  # actually sit".
  #
  # PEARSON is what goes on the strip, not Spearman, so these numbers are the
  # same quantity already in circulation on this project: Vincent Rollet's email
  # figures (0.43 recent/large, 0.21 all-periods 16-country, 0.14 full sample),
  # Sara's "the Pearson correlation is 0.14" on full V-Party, and the "Raw
  # correlation" column of adhoc/vparty_corr_by_decade.R. Verified: this
  # script's per-panel Pearson and N reproduce that script's ten rows exactly
  # (0.2154, 0.1113, 0.4020, 0.4458, 0.5262 / -0.3810, -0.2686, -0.0397,
  # 0.0031, -0.0151), which is also a check that the two build the same sample.
  # Spearman is kept in the CSV, since the Jaccard bins are themselves ranks and
  # it is the natural companion there.
  summarise_cells <- function(cells, df) {
    assoc <- df |>
      summarise(
        pearson = cor(.data[[A]], .data[[B]]),
        spearman = cor(.data[[A]], .data[[B]], method = "spearman"),
        .by = c(group, decade)
      )
    cells |>
      summarise(
        n_party_years = first(n_pairs),
        diag_mean = mean(jaccard[row_i == col_i], na.rm = TRUE),
        offdiag_mean = mean(jaccard[row_i != col_i], na.rm = TRUE),
        .by = c(group, decade)
      ) |>
      mutate(diag_ratio = diag_mean / offdiag_mean) |>
      left_join(assoc, by = c("group", "decade"))
  }

  # ---- plot --------------------------------------------------------------------

  # Single-hue sequential, matching 16_instrument_overlap.R so the two figures
  # are read the same way. ONE scale across all ten panels: the whole point is
  # the left-to-right comparison, and a per-panel scale would make every decade
  # look equally concentrated.
  panel_figure <- function(cells, summary, universe_label, fill_max) {
    thin <- summary |> filter(n_party_years < MIN_CELL_N)
    drawn <- cells |>
      semi_join(
        summary |> filter(n_party_years >= MIN_CELL_N),
        by = c("group", "decade")
      )

    # Each panel states the N it rests on, so one built from 200 party-years is
    # not read as the equal of one built from 1,000.
    #
    # facet_wrap on a combined label rather than facet_grid(group ~ decade),
    # which would give tidier shared strips but has nowhere to put the N: the two
    # rows have DIFFERENT Ns for the same decade, so a shared column strip could
    # only ever show one of them. Annotating inside the panel is not an option
    # either -- the tiles cover the full [0.5, 10.5] range on both axes, so any
    # in-range label sits on top of data. Levels are ordered OECD-then-Non-OECD
    # by decade, which with ncol = 5 reproduces exactly the 2 x 5 layout.
    strip_levels <- expand_grid(
      group = levels(summary$group),
      decade = levels(summary$decade)
    ) |>
      left_join(summary, by = c("group", "decade")) |>
      mutate(strip = sprintf(
        "%s \u00b7 %s \u00b7 n = %s \u00b7 r = %s",
        group, decade,
        ifelse(is.na(n_party_years), "0",
               format(n_party_years, big.mark = ",")),
        ifelse(is.na(pearson), "--", sprintf("%.2f", pearson))
      )) |>
      pull(strip)

    drawn <- drawn |>
      left_join(
        summary |> mutate(strip = sprintf(
          "%s \u00b7 %s \u00b7 n = %s \u00b7 r = %s",
          group, decade, format(n_party_years, big.mark = ","),
          sprintf("%.2f", pearson)
        )) |> select(group, decade, strip),
        by = c("group", "decade")
      ) |>
      mutate(strip = factor(strip, levels = strip_levels))

    ggplot(drawn, aes(x = col_i, y = row_i, fill = jaccard)) +
      geom_tile(colour = "white", linewidth = 0.25) +
      facet_wrap(~strip, ncol = nlevels(summary$decade), drop = FALSE) +
      scale_fill_gradient(
        low = "#eaf2f8", high = "#0072B2", na.value = "grey93",
        limits = c(0, fill_max), name = "Jaccard"
      ) +
      # Bin INDICES on both axes rather than the "0.0-0.1" labels: ten text
      # labels per axis times ten panels is unreadable, and the indices carry the
      # same information once the axis title says what they are.
      #
      # y runs UPWARD (bin 1 at the bottom, bin 10 at the top), not in matrix
      # order. 16_instrument_overlap.R reverses it because it is drawing a
      # matrix; here the figure sits next to 19_vparty_ideology_quadrants.R,
      # which plots anti-pluralism on a conventional axis, and flipping between
      # the two invites reading one of them upside down. The leading diagonal
      # therefore runs bottom-left to top-right, as on a scatter.
      scale_x_continuous(breaks = c(1, 5, 10), expand = c(0, 0)) +
      scale_y_continuous(breaks = c(1, 5, 10), expand = c(0, 0)) +
      coord_fixed() +
      labs(
        title = sprintf(
        "%s (y) vs %s (x): %s", vparty_score(A, "display"),
        vparty_score(B, "display"), universe_label),
        subtitle = paste(
          strwrap(sprintf(
            paste(
              "Jaccard = the count in both bins over the count in either bin,",
              "across %s %s.",
              "A perfectly redundant pair would light only the diagonal.",
              "One shared colour scale across all panels. r is the Pearson",
              "correlation in that panel -- the same quantity as the raw",
              "correlation series in adhoc/vparty_corr_by_decade.R. Read IT for",
              "whether the two indices are related, because Jaccard runs on raw",
              "counts and is bright wherever both marginals pile up, related or",
              "not. V-Party %d-%d.%s"
            ),
            N_BINS,
            if (bin_mode == "equal01") {
              "equal-width bins of each index on [0, 1], so the diagonal is a same-SCORE diagonal"
            } else {
              paste(
                "quantile bins of each index, used because at least one axis",
                "is an expert scale rather than a [0, 1] index -- so the",
                "diagonal here is a same-RANK diagonal, not a same-score one"
              )
            },
            VPARTY_YEAR_MIN, VPARTY_YEAR_MAX,
            if (nrow(thin) > 0) {
              sprintf(
                " %d panel(s) with fewer than %d party-years are left blank.",
                nrow(thin), MIN_CELL_N
              )
            } else {
              ""
            }
          ), 115),
          collapse = "\n"
        ),
        x = sprintf("%s bin (1 = lowest, %d = highest)", vparty_score(B, "display"), N_BINS),
        y = sprintf("%s bin (1 = lowest, %d = highest)", vparty_score(A, "display"), N_BINS)
      ) +
      theme_bw(base_size = 9) +
      theme(
        panel.grid = element_blank(),
        panel.spacing = unit(0.35, "lines"),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(size = 8),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "grey25"),
        legend.key.height = unit(1.1, "cm"),
        legend.key.width = unit(0.35, "cm")
      )
  }

  # ---- run ---------------------------------------------------------------------

  all_summaries <- list()

  # Computed for every universe BEFORE anything is drawn, so the fill scale can
  # be shared across figures as well as across panels.
  all_cells <- list()
  for (u in names(UNIVERSES)) {
    cells <- build_cells(UNIVERSES[[u]]$fn(base))
    if (nrow(cells) == 0) {
      warning("No cells for universe ", u, "; skipping.", call. = FALSE)
      next
    }
    all_cells[[u]] <- cells
  }
  stopifnot(length(all_cells) > 0)
  fill_max <- max(bind_rows(all_cells)$jaccard, na.rm = TRUE)
  cat(sprintf("\nShared fill scale: 0 to %.3f\n", fill_max))

  for (u in names(all_cells)) {
    spec <- UNIVERSES[[u]]
    cells <- all_cells[[u]]
    summary <- summarise_cells(cells, spec$fn(base))
    all_summaries[[u]] <- summary |>
      mutate(universe = u, score_y = A, score_x = B, .before = 1)

    write_csv(with_decade_key(cells), file.path(out_dir, sprintf("jaccard_cells_%s_%s.csv", pslug, u)))

    fig <- panel_figure(cells, summary, spec$label, fill_max)
    ggsave(
      file.path(out_dir, sprintf("jaccard_panels_%s_%s.png", pslug, u)),
      fig,
      width = 13, height = 6.2, dpi = 150
    )
    cat(sprintf("Saved jaccard_panels_%s_%s.png\n", pslug, u))

    cat(sprintf("\n=== %s ===\n", spec$label))
    print(
      summary |>
        mutate(across(where(is.numeric), \(x) round(x, 4))) |>
        as.data.frame(),
      row.names = FALSE
    )
  }

  write_csv(with_decade_key(bind_rows(all_summaries)), file.path(out_dir, sprintf("panel_summary_%s.csv", pslug)))
}

message("\nV-Party Jaccard panels written to ", out_dir)
