# ==============================================================================
# Shared helpers for analyses over the RAW V-Party dataset.
#
# Distinct from rdd_helpers.R on purpose. rdd_helpers.R is about the RDD's
# election spine -- runs, slugs, thresholds, rdrobust. Everything here is about
# V-Party as a party-year dataset in its own right: the file, the OECD split,
# the decade split, and the Jaccard binning. Scripts 18 and 19 answer
# measurement questions about the indices and never touch a build, so they
# source this and not rdd_helpers.R.
#
# 16_instrument_overlap.R sources BOTH: it does the same Jaccard computation on
# the top-2 spine rather than on raw V-Party.
#
# Sourced, not run. Assumes tidyverse + here are attached by the caller.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. The dataset
# ------------------------------------------------------------------------------

# V-Party v2, read straight out of the zip -- no election spine, no top-2
# filter, no crosswalk. Anything computed on this is a property of the indices,
# not of this project's sample construction, which is the whole point of the
# raw-V-Party analyses.
#
# 11_build_rdd_data.R reads the same file but immediately select()s the five
# instrument columns and joins to the spine, so it is not a substitute.
load_vparty_raw <- function(cols = NULL) {
  path <- here::here("data", "elections_database", "CPD_V-Party_CSV_v2.zip")
  if (!file.exists(path)) {
    stop("No V-Party zip at ", path, call. = FALSE)
  }
  d <- readr::read_csv(
    unz(path, "CPD_V-Party_CSV_v2/V-Dem-CPD-Party-V2.csv"),
    show_col_types = FALSE
  )
  if (!is.null(cols)) {
    missing <- setdiff(cols, names(d))
    if (length(missing) > 0) {
      stop(
        "Not in V-Party: ", paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    d <- d[, cols, drop = FALSE]
  }
  d
}

# ------------------------------------------------------------------------------
# 2. Party scores: slug, display name, and scale type
#
# Scripts 18 and 19 plot one PAIR of scores against each other and need three
# things per variable: a filesystem-safe slug for the output name, a display
# name for the axis, and whether it is a [0, 1] index.
#
# The scale type is not decoration. It decides two things that silently break
# otherwise: Jaccard binning (equal-width bins are only meaningful on [0, 1];
# bin_breaks() errors on an expert scale) and coord_fixed() (equal aspect is
# only meaningful when a unit means the same on both axes, which it does not
# when one runs 0-1 and the other -4 to +4).
#
# rdd_helpers.R has a similar registry, but 18 and 19 deliberately do not
# source it -- they never touch a build -- so this is the V-Party-side copy.
# It carries the RAW v2pariglef, right-positive, not the negated version the
# RDD uses: on a descriptive scatter "left-right" should read left-to-right.
VPARTY_SCORES <- list(
  v2xpa_antiplural = list(
    slug = "antiplural",
    display = "Anti-pluralism (v2xpa_antiplural)",
    unit_interval = TRUE
  ),
  v2xpa_popul = list(
    slug = "popul",
    display = "Populism (v2xpa_popul)",
    unit_interval = TRUE
  ),
  v2pariglef = list(
    slug = "econlr",
    display = "Economic left-right (v2pariglef; higher = RIGHT)",
    unit_interval = FALSE
  ),
  v2paanteli = list(
    slug = "anteli",
    display = "Anti-elitism (v2paanteli)",
    unit_interval = FALSE
  ),
  v2paminor = list(
    slug = "minor",
    display = "Minority rights (v2paminor; higher = MORE supportive)",
    unit_interval = FALSE
  )
)

vparty_score <- function(var, field) {
  if (is.null(VPARTY_SCORES[[var]])) {
    stop(
      "No entry for '", var, "' in VPARTY_SCORES (scripts/vparty_helpers.R). ",
      "Add one giving its slug, display name and whether it is a [0, 1] index.",
      call. = FALSE
    )
  }
  VPARTY_SCORES[[var]][[field]]
}

# Equal-width bins mean "same score" and are only defined on [0, 1]; quantile
# bins mean "same rank" and work on any scale. A pair is binned by score only
# when BOTH axes allow it, so the diagonal means one thing rather than two.
pair_bin_mode <- function(var_a, var_b) {
  if (vparty_score(var_a, "unit_interval") && vparty_score(var_b, "unit_interval")) {
    "equal01"
  } else {
    "deciles"
  }
}

pair_slug <- function(var_a, var_b) {
  sprintf("%s_vs_%s", vparty_score(var_a, "slug"), vparty_score(var_b, "slug"))
}

# ------------------------------------------------------------------------------
# 2. The two splits
# ------------------------------------------------------------------------------

# Current (2026) OECD membership, applied to the WHOLE period -- so the split is
# "countries that ended up rich democracies", not membership as of each
# election. Switching to date-of-accession membership would empty the 1970s
# non-OECD-then cells of exactly the countries that make the OECD line.
#
# Previously duplicated verbatim in adhoc/vparty_corr_by_decade.R and
# adhoc/vparty_variance_by_decade.R, whose own header flagged the duplication.
OECD_ISO3 <- c(
  "AUS", "AUT", "BEL", "CAN", "CHL", "COL", "CRI", "CZE", "DNK", "EST",
  "FIN", "FRA", "DEU", "GRC", "HUN", "ISL", "IRL", "ISR", "ITA", "JPN",
  "KOR", "LVA", "LTU", "LUX", "MEX", "NLD", "NZL", "NOR", "POL", "PRT",
  "SVK", "SVN", "ESP", "SWE", "CHE", "TUR", "GBR", "USA"
)

# OECD first, so it is the TOP row of a facet_grid and the left panel of a
# facet_wrap without every caller re-stating the level order.
oecd_group <- function(country_text_id) {
  factor(
    ifelse(country_text_id %in% OECD_ISO3, "OECD", "Non-OECD"),
    levels = c("OECD", "Non-OECD")
  )
}

# V-Party reaches back to 1900, but coverage before 1970 is thin and lopsided,
# and 2019 is the last year in v2. That window is also exactly five decades,
# which is what the 2 x 5 decade panels need.
VPARTY_YEAR_MIN <- 1970
VPARTY_YEAR_MAX <- 2019

# "1970s", "1980s", ... as an ordered factor, so decades read left to right in
# a facet_grid without relying on alphabetical luck.
decade_label <- function(year) {
  dec <- (year %/% 10) * 10
  lvls <- paste0(seq(
    (VPARTY_YEAR_MIN %/% 10) * 10,
    (VPARTY_YEAR_MAX %/% 10) * 10,
    by = 10
  ), "s")
  factor(paste0(dec, "s"), levels = lvls)
}

# ------------------------------------------------------------------------------
# 3. Jaccard binning
#
# Moved here from 16_instrument_overlap.R so the top-2 heatmaps and the raw
# V-Party decade panels compute the same statistic on the same bins. The two
# globals 16 used (JACCARD_BINS, JACCARD_N_BINS) are arguments now -- a shared
# function that silently reads its caller's globals is worse than a duplicate.
# ------------------------------------------------------------------------------

# Bin edges for one variable.
#
# mode = "equal01": ten fixed bins on [0, 1], [0, 0.1], (0.1, 0.2], ... The
#   edges mean the same thing in both dimensions and across every sample, so
#   cells are comparable between heatmaps and the diagonal is a genuine
#   "same score" diagonal.
# mode = "deciles": quantile bins. A different question -- they equalise cell
#   counts, which makes the diagonal "same RANK" rather than "same score", and
#   they move with whichever sample is being cut.
#
# Equal-width bins are only meaningful for a variable actually scaled to [0, 1].
# v2xpa_antiplural and v2xpa_popul are; v2pariglef_neg (about -1.9 to 3.8),
# v2paanteli (-2.4 to 4.4) and ep_galtan (4.5 to 9.4) are not, and silently
# binning them onto [0, 1] would drop nearly every observation out of range.
bin_breaks <- function(x, var, mode = c("equal01", "deciles"), n_bins = 10) {
  mode <- match.arg(mode)
  if (mode == "deciles") {
    return(unique(quantile(
      x,
      probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE
    )))
  }
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng))) {
    stop("bin_breaks: ", var, " has no non-missing values.", call. = FALSE)
  }
  if (rng[1] < 0 || rng[2] > 1) {
    stop(
      "bin mode 'equal01' needs a variable on [0, 1], but ", var,
      " ranges ", sprintf("%.2f to %.2f", rng[1], rng[2]),
      ". Use mode = 'deciles' for expert-scale variables.",
      call. = FALSE
    )
  }
  seq(0, 1, length.out = n_bins + 1)
}

# Cell (i, j) = |bin_a == i AND bin_b == j| / |bin_a == i OR bin_b == j|.
# A perfectly redundant pair of measures would light the diagonal and nothing
# else. Bins are fixed up front, so an EMPTY bin still gets its row/column and
# renders as NA rather than silently shifting the grid.
jaccard_matrix <- function(df, n_a, n_b, lab_a, lab_b) {
  m <- matrix(NA_real_, nrow = n_a, ncol = n_b)
  for (i in seq_len(n_a)) {
    for (j in seq_len(n_b)) {
      in_a <- df$bin_a == i
      in_b <- df$bin_b == j
      union_n <- sum(in_a | in_b, na.rm = TRUE)
      m[i, j] <- if (union_n == 0) {
        NA_real_
      } else {
        sum(in_a & in_b, na.rm = TRUE) / union_n
      }
    }
  }
  dimnames(m) <- list(lab_a, lab_b)
  m
}

# Long-format Jaccard cells for one (var_a, var_b) pair on one data frame.
# Returns row/col bin indices, their labels and the Jaccard value -- the shape
# both ggplot and the CSVs want. Callers that need the matrix itself still call
# jaccard_matrix() directly.
jaccard_long <- function(df, var_a, var_b, breaks_a, breaks_b) {
  lab_a <- sprintf("%.1f-%.1f", head(breaks_a, -1), tail(breaks_a, -1))
  lab_b <- sprintf("%.1f-%.1f", head(breaks_b, -1), tail(breaks_b, -1))
  binned <- data.frame(
    bin_a = cut(df[[var_a]], breaks_a, include.lowest = TRUE, labels = FALSE),
    bin_b = cut(df[[var_b]], breaks_b, include.lowest = TRUE, labels = FALSE)
  )
  binned <- binned[!is.na(binned$bin_a) & !is.na(binned$bin_b), , drop = FALSE]
  m <- jaccard_matrix(binned, length(lab_a), length(lab_b), lab_a, lab_b)
  tibble::tibble(
    row_i = rep(seq_along(lab_a), times = length(lab_b)),
    col_i = rep(seq_along(lab_b), each = length(lab_a)),
    row_lab = factor(rep(lab_a, times = length(lab_b)), levels = lab_a),
    col_lab = factor(rep(lab_b, each = length(lab_a)), levels = lab_b),
    jaccard = as.vector(m),
    n_pairs = nrow(binned)
  )
}
