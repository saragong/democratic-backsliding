# ==============================================================================
# The PopuList 411, one row per election.
#
# The headline restricted sample: elections where one of the top 2 parties
# scores above the PopuList-calibrated illiberality cut and the other does not.
# This lists every one of them by name, with the pre-election level and the
# 1-, 5- and 10-year change of both outcomes, so the sample can be read
# directly instead of through a coefficient.
#
# Winner and loser are by VOTE SHARE, not by score -- so the illiberal party is
# sometimes the winner and sometimes the loser, which is the whole point of the
# design. The gap column is signed winner minus loser:
#
#   POSITIVE  the WINNER is the more illiberal of the two  -> shaded RED
#   NEGATIVE  the winner is the LESS illiberal of the two  -> shaded GREEN
#
# which needs diverging_fill()'s argument negated, since that ramp runs red for
# negative and green for positive.
#
# The script ends by re-running the headline RDD on nothing but the columns
# this table displays, and checking it against the canonical run folder. See
# "sanity check" at the bottom.
#
#   Rscript --no-init-file adhoc/popucut_411_table.R
#
# Output: output/runs/_sweeps/popucut_411_table/
#           popucut_411.html   shaded, grouped by country
#           popucut_411.csv    same numbers, unformatted
#           sanity_check_rdd.csv
#         adhoc/popucut_411.csv  the same frame, next to this script
# ==============================================================================

library(tidyverse)
library(here)
library(gt)
library(rdrobust)

source(here::here("scripts", "rdd_helpers.R"))
source(here::here("scripts", "vparty_helpers.R"))

# ---- toggles -----------------------------------------------------------------

if (!exists("T411_INSTRUMENT")) T411_INSTRUMENT <- "v2xpa_antiplural"
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- FALSE
}
# The horizons shown as change columns. The sample itself does not depend on
# them -- the same 411 elections appear at every window.
if (!exists("T411_HORIZONS")) T411_HORIZONS <- c(1, 5, 10)
# The window whose build defines the sample. Any of them would give the same
# 411; 5 is the headline.
if (!exists("T411_SAMPLE_WINDOW")) T411_SAMPLE_WINDOW <- 5

build_suffix <- if (TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR) "" else "_exclyr"
build_path <- function(n) {
  here::here("data", "rdd_build", sprintf(
    "rdd_%s_w%d%s.rds", T411_INSTRUMENT, n, build_suffix
  ))
}
parties_path <- here::here("data", "rdd_build", sprintf(
  "rdd_%s_w%d%s_parties.rds", T411_INSTRUMENT, T411_SAMPLE_WINDOW, build_suffix
))

needed <- unique(c(T411_SAMPLE_WINDOW, T411_HORIZONS))
missing <- needed[!file.exists(vapply(needed, build_path, character(1)))]
if (length(missing) > 0 || !file.exists(parties_path)) {
  stop(
    "Missing build(s) for window(s) ", paste(missing, collapse = ", "),
    ". Run 11_build_rdd_data.R at those windows first.", call. = FALSE
  )
}

out_dir <- sweep_dir("popucut_411_table")

# "Anti-pluralism (v2xpa_antiplural)" is right for a subtitle and far too long
# for a column header repeated twice. Strip the parenthetical -- this yields a
# usable short name for every instrument in the registry ("Populism",
# "Anti-elitism", "Minority rights"), so it does not need a second lookup
# table that could drift out of step with the first.
SCORE_SHORT <- sub(" \\(.*$", "", PARTY_SCORE_DISPLAY[[T411_INSTRUMENT]])

# ---- the sample --------------------------------------------------------------

cut_abs <- resolve_threshold_abs("popucut", "popucut", NULL)

d <- readRDS(build_path(T411_SAMPLE_WINDOW))
sample_ids <- d |>
  filter(illiberal_score > cut_abs, other_score <= cut_abs) |>
  pull(election_id)

cat(sprintf(
  "PopuList cut = %.4f; %d of %d elections have one top-2 party above it and one at or below.\n",
  cut_abs, length(sample_ids), nrow(d)
))

# ---- winner and loser, by vote share ----------------------------------------

parties <- readRDS(parties_path) |>
  filter(election_id %in% sample_ids) |>
  left_join(vparty_names_at(readRDS(parties_path)),
            by = c("vdem_id_1", "election_year"))

wl <- parties |>
  group_by(election_id) |>
  arrange(desc(final_share), .by_group = TRUE) |>
  summarise(
    win_party = party_label(first(party_name), first(party_abbr), first(candidate)),
    win_share = first(final_share),
    win_score = first(.data[[T411_INSTRUMENT]]),
    lose_party = party_label(last(party_name), last(party_abbr), last(candidate)),
    lose_share = last(final_share),
    lose_score = last(.data[[T411_INSTRUMENT]]),
    .groups = "drop"
  ) |>
  mutate(
    score_diff = win_score - lose_score,
    vote_margin = win_share - lose_share,
    # The RDD's running variable is the margin of the ILLIBERAL party, which
    # is not the margin of the winner: vote_margin is always positive, and the
    # sign comes from which side the illiberal party was on, which is exactly
    # what score_diff records. Carried in the exported data so the CSV alone
    # is enough to re-run the design; the sanity check at the bottom asserts
    # it against the build's own running_var.
    running_var = vote_margin * sign(score_diff)
  )

# By construction of the sample the two scores straddle the cut, so they can
# never tie -- which also means none of the 22 tied-score elections, or the 4
# where the same party fills both top-2 slots, can reach this table.
stopifnot(all(wl$score_diff != 0), nrow(wl) == length(sample_ids))

# ---- levels and changes ------------------------------------------------------

# Read the outcome LEVELS straight from the country-year panel rather than
# inverting the build's log changes: the table shows a level in the
# pre-election year, and a level has to come from somewhere anyway. The
# reconciliation check below is what ties this back to the analysis.
panel <- readRDS(here::here("data", "combined_panel.rds")) |>
  select(country_text_id, year, ln_gdp_pc, v2x_polyarchy)

base <- d |>
  filter(election_id %in% sample_ids) |>
  select(election_id, country_name, country_text_id, election_year, election_type)

level_at <- function(offset) {
  base |>
    mutate(year = election_year + offset) |>
    left_join(panel, by = c("country_text_id", "year")) |>
    select(election_id, ln_gdp_pc, v2x_polyarchy)
}

pre <- level_at(-1L) |>
  rename(ln_gdp_pre = ln_gdp_pc, polyarchy_pre = v2x_polyarchy)

changes <- reduce(
  map(T411_HORIZONS, function(h) {
    level_at(h) |>
      left_join(pre, by = "election_id") |>
      transmute(
        election_id,
        # exp(delta ln) - 1: a percentage change, which reads far more
        # naturally next to a dollar level than log points do.
        !!sprintf("gdp_chg_%dy", h) := exp(ln_gdp_pc - ln_gdp_pre) - 1,
        # Polyarchy in INDEX POINTS, not log points. The build's Y_polyarchy is
        # a log change (see the reconciliation below), but the index is bounded
        # [0,1] and a reader comparing it to the level in the adjacent column
        # wants the same units.
        !!sprintf("poly_chg_%dy", h) := v2x_polyarchy - polyarchy_pre
      )
  }),
  left_join, by = "election_id"
)

tbl_dat <- base |>
  left_join(wl, by = "election_id") |>
  left_join(pre, by = "election_id") |>
  left_join(changes, by = "election_id") |>
  mutate(gdp_pc_pre = exp(ln_gdp_pre)) |>
  arrange(country_name, election_year)

# RECONCILIATION. The changes above are computed here, from the panel; the
# build computes its own over the same span inside window_change(). They must
# agree, or this table is describing a different window from the one every
# estimate in the repo uses. Checked at each horizon on the rows where both
# are observed, rather than asserted in a comment.
for (h in T411_HORIZONS) {
  b <- readRDS(build_path(h)) |>
    filter(election_id %in% sample_ids) |>
    select(election_id, Y_gdp_growth, Y_polyarchy)
  chk <- tbl_dat |>
    select(election_id, g = !!sprintf("gdp_chg_%dy", h)) |>
    left_join(b, by = "election_id") |>
    filter(!is.na(g), !is.na(Y_gdp_growth))
  stopifnot(nrow(chk) > 0)
  worst <- max(abs(log1p(chk$g) - chk$Y_gdp_growth))
  cat(sprintf(
    "  w%-2d reconciled against the build on %d elections (max |diff| in log points: %.2e)\n",
    h, nrow(chk), worst
  ))
  stopifnot(worst < 1e-8)
}

# Two destinations on purpose. The copy in the run folder is the companion to
# the HTML and travels with it; the copy in adhoc/ is the dataset itself, next
# to the script that builds it, where it can be picked up without knowing the
# sweep-folder convention. Same frame, written once from one object, so they
# cannot disagree.
csv_paths <- c(
  file.path(out_dir, "popucut_411.csv"),
  here::here("adhoc", "popucut_411.csv")
)
walk(csv_paths, function(f) {
  write_csv(tbl_dat, f)
  cat(sprintf("Saved %s\n", f))
})

# ---- the table ---------------------------------------------------------------

gdp_cols <- sprintf("gdp_chg_%dy", T411_HORIZONS)
poly_cols <- sprintf("poly_chg_%dy", T411_HORIZONS)

# Shading domains. One symmetric domain PER BLOCK, shared across that block's
# horizons, so a colour means the same magnitude at 1, 5 and 10 years and the
# eye can follow a row across. Capped at the 95th percentile of |value|:
# hyperinflation-era collapses would otherwise take the whole scale and leave
# every ordinary election white. diverging_fill() clamps, so the capped cases
# render at full intensity rather than dropping out.
cap <- function(cols) {
  v <- abs(unlist(tbl_dat[cols], use.names = FALSE))
  v <- v[is.finite(v)]
  if (length(v) == 0) 1 else unname(quantile(v, 0.95))
}
gdp_cap <- cap(gdp_cols)
poly_cap <- cap(poly_cols)
diff_cap <- max(abs(tbl_dat$score_diff), na.rm = TRUE)

# One diverging ramp everywhere, so red always means the direction a reader of
# this paper cares about: worse growth, falling polyarchy, and an illiberal
# party winning. `flip` negates the input for the gap column, whose raw sign
# runs the other way (positive = illiberal winner, which must be RED).
shade <- function(gt_tbl, cols, max_abs, flip = FALSE) {
  data_color(
    gt_tbl, columns = all_of(cols),
    fn = function(x) diverging_fill(if (flip) -x else x, max_abs = max_abs)
  )
}

disp <- tbl_dat |>
  transmute(
    Country = country_name,
    Year = election_year,
    Type = recode(election_type, presidential = "Pres.", parliamentary = "Parl."),
    # Who ran and how the vote fell first, then how illiberal each side was.
    # The three score columns sit together in their own spanner rather than
    # one inside each of Winner and Loser: the reader wants to compare the two
    # scores and their gap side by side, and it is the only arrangement that
    # puts the margin ahead of them without splitting a spanner.
    `Winning party` = win_party, `Vote %` = win_share,
    `Losing party` = lose_party, `Vote %.` = lose_share,
    `Margin (pp)` = vote_margin,
    `Score` = win_score, `Score.` = lose_score, `Score diff` = score_diff,
    `GDP pc (pre)` = gdp_pc_pre,
    !!!setNames(map(gdp_cols, ~ tbl_dat[[.x]]), sprintf("GDP %dy", T411_HORIZONS)),
    `Polyarchy (pre)` = polyarchy_pre,
    !!!setNames(map(poly_cols, ~ tbl_dat[[.x]]), sprintf("Poly %dy", T411_HORIZONS))
  )

gdp_disp <- sprintf("GDP %dy", T411_HORIZONS)
poly_disp <- sprintf("Poly %dy", T411_HORIZONS)

gt_tbl <- disp |>
  gt(groupname_col = "Country") |>
  tab_header(
    title = sprintf(
      "The PopuList %d: one top-2 party illiberal, the other not",
      length(sample_ids)
    ),
    subtitle = html(sprintf(paste(
      "Both party scores are <b>%s</b>, V-Party's illiberalism measure, on",
      "[0,&thinsp;1] with higher = more illiberal; the sample is elections",
      "where one top-2 party scores above %.4f and the other does not (the",
      "accuracy-maximizing PopuList cut, calibrated in",
      "01g_populist_threshold.R). Winner and loser are by vote or seat share,",
      "so the illiberal party is sometimes which ever one won:",
      "<b>a positive gap means the WINNER is the more illiberal of the two",
      "(shaded <span style=\"color:#b2182b\">red</span>); a negative gap means",
      "the winner is the LESS illiberal",
      "(shaded <span style=\"color:#1a9850\">green</span>)</b>. Changes are",
      "measured from the year before the election, the same span every",
      "estimate in the repo uses: GDP per capita as a percentage change,",
      "polyarchy in index points."
    ), INSTRUMENT_DISPLAY[[T411_INSTRUMENT]], cut_abs))
  ) |>
  tab_spanner("Winner", columns = c("Winning party", "Vote %")) |>
  tab_spanner("Loser", columns = c("Losing party", "Vote %.")) |>
  tab_spanner(
    html(sprintf("%s score", SCORE_SHORT)),
    columns = c("Score", "Score.", "Score diff")
  ) |>
  tab_spanner("GDP pc change", columns = all_of(gdp_disp)) |>
  tab_spanner("Polyarchy change", columns = all_of(poly_disp)) |>
  cols_label(
    `Winning party` = "Party", `Losing party` = "Party", `Vote %.` = "Vote %",
    `Score` = "Winner", `Score.` = "Loser",
    `Score diff` = html(
      "Gap<br><span style=\"font-weight:normal\">(winner &minus; loser)</span>"
    ),
    !!!setNames(as.list(paste0(T411_HORIZONS, "y")), gdp_disp),
    !!!setNames(as.list(paste0(T411_HORIZONS, "y")), poly_disp)
  ) |>
  fmt_number(columns = c("Vote %", "Vote %.", "Margin (pp)"), decimals = 1) |>
  fmt_number(columns = c("Score", "Score.", "Score diff", "Polyarchy (pre)"), decimals = 3) |>
  fmt_number(columns = "GDP pc (pre)", decimals = 0, use_seps = TRUE) |>
  fmt_percent(columns = all_of(gdp_disp), decimals = 0) |>
  fmt_number(columns = all_of(poly_disp), decimals = 3, force_sign = TRUE) |>
  sub_missing(missing_text = "--") |>
  shade(gdp_disp, gdp_cap) |>
  shade(poly_disp, poly_cap) |>
  shade("Score diff", diff_cap, flip = TRUE) |>
  cols_align("left", columns = c("Winning party", "Losing party")) |>
  tab_style(cell_text(weight = "bold"), locations = cells_row_groups()) |>
  apply_table_style(font_size = 10) |>
  tab_source_note(sprintf(paste(
    "Shading is a diverging scale centred at zero -- red below, green above --",
    "with one symmetric domain per block so a colour means the same magnitude",
    "at 1, 5 and 10 years. The domain is capped at the 95th percentile of",
    "|value| (GDP %.0f%%, polyarchy %.3f) so that a handful of collapses do",
    "not wash out the rest; capped cells render at full intensity. The %s gap",
    "is shaded on the same red-green scale over its own full range, with the",
    "sign reversed so that RED is an election the more illiberal party WON and",
    "GREEN one it lost -- red therefore means the same direction everywhere in",
    "the table: worse growth, falling polyarchy, an illiberal party winning.",
    "The deepest cells either way are the elections where the two top-2",
    "parties are furthest apart. %d elections, %d countries, %d-%d; the",
    "illiberal party won %d and lost %d. Blank cells are years the outcome",
    "panel does not cover."
  ), 100 * gdp_cap, poly_cap, SCORE_SHORT,
  nrow(disp), n_distinct(disp$Country),
  min(disp$Year), max(disp$Year),
  sum(tbl_dat$score_diff > 0), sum(tbl_dat$score_diff < 0)))

gtsave(gt_tbl, file.path(out_dir, "popucut_411.html"))
cat(sprintf("Saved %s\n", file.path(out_dir, "popucut_411.html")))

cat(sprintf(
  "\n%d elections, %d countries, %d-%d. Illiberal party WON %d, LOST %d.\n",
  nrow(tbl_dat), n_distinct(tbl_dat$country_name),
  min(tbl_dat$election_year), max(tbl_dat$election_year),
  sum(tbl_dat$score_diff > 0), sum(tbl_dat$score_diff < 0)
))

# ==============================================================================
# Sanity check: does this table reproduce the headline RDD?
#
# The table is assembled by a different route from the analysis -- its sample
# is filtered here, its winner/loser ordering is by vote share, and its
# outcome columns are recomputed from combined_panel.rds. If all of that is
# right, the headline estimate has to fall out of the displayed columns alone.
#
# The running variable is tbl_dat$running_var, built up where wl is assembled
# as margin x sign(gap) -- from the two columns this table displays, never
# read off the build. That is the point: a sign error in the gap column would
# surface here as an estimate with the wrong sign instead of passing silently.
# Step 1 below is what closes that loop.
# ==============================================================================

cat("\n", strrep("-", 78), "\nSANITY CHECK: the headline RDD, rebuilt from this table's own columns\n",
    strrep("-", 78), "\n", sep = "")

rv <- tbl_dat$running_var

# Step 1: the reconstruction must reproduce the build's running variable.
# Vote shares are carried at full precision in the parties file, so this is an
# equality check, not an approximate one.
rv_build <- d |>
  filter(election_id %in% sample_ids) |>
  arrange(match(election_id, tbl_dat$election_id)) |>
  pull(running_var)
stopifnot(identical(tbl_dat$election_id, d$election_id[match(tbl_dat$election_id, d$election_id)]))
max_rv_gap <- max(abs(rv - rv_build))
cat(sprintf(
  "  running variable rebuilt from Margin x sign(gap): max |diff| = %.2e over %d elections\n",
  max_rv_gap, length(rv)
))
stopifnot(max_rv_gap < 1e-9)
cat(sprintf(
  "  treated side (illiberal party won) = %d, control side = %d\n",
  sum(rv > 0), sum(rv < 0)
))

# Step 2: the RDD itself, on log1p of the displayed percentage change, which
# is the log change the analysis uses.
canon_path <- function(h) {
  file.path(RUNS_ROOT, run_slug(list(
    instrument = T411_INSTRUMENT, window = h, treatment = "backsliding_Nyr",
    score_gap_min = -Inf, illiberal_cutoff = cut_abs, other_cutoff_max = cut_abs,
    incl_election_year = TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR
  )), "rdd_results.csv")
}

# The treatment is the one input the table does not display; it is pulled from
# the build purely so the fuzzy arm can be reproduced too. It must come from
# the build AT HORIZON h, not from the sample build: backsliding_Nyr asks
# whether an episode starts within h years, so it is a different variable at
# every horizon. Reusing the w5 treatment against the w1 and w10 outcomes made
# the fuzzy bandwidths disagree while the reduced form matched exactly --
# which is what a window-mismatched treatment looks like.
trt_at <- function(h) {
  readRDS(build_path(h)) |>
    filter(election_id %in% sample_ids) |>
    arrange(match(election_id, tbl_dat$election_id)) |>
    pull(backsliding_Nyr)
}

check <- map_dfr(T411_HORIZONS, function(h) {
  y <- log1p(tbl_dat[[sprintf("gdp_chg_%dy", h)]])
  trt <- trt_at(h)
  rf <- extract_rd(safe_rdrobust(y, rv))
  late <- extract_rd(safe_rdrobust(y, rv, fuzzy = trt))
  row <- tibble(
    window = h, n_table = rf$N,
    est_table = rf$coef, se_table = rf$se, pval_table = rf$pval,
    # BOTH bandwidths, because rdd_results.csv's "bandwidth" column is the
    # LATE's, not the reduced form's (12_rdd_analysis.R line 382). Comparing
    # the reduced form's against it would fail on a naming difference and look
    # like a substantive mismatch -- which is exactly what happened when this
    # check was first written.
    bw_rf_table = rf$bw, bw_late_table = late$bw,
    late_table = late$coef
  )
  cp <- canon_path(h)
  if (!file.exists(cp)) {
    return(mutate(row, n_run = NA_integer_, est_run = NA_real_,
                  se_run = NA_real_, bw_run = NA_real_, late_run = NA_real_,
                  source = "run folder missing"))
  }
  canon <- read_csv(cp, show_col_types = FALSE) |> filter(outcome == "Y_gdp_growth")
  mutate(
    row, n_run = canon$n_outcome, est_run = canon$rd_estimate,
    se_run = canon$rd_se, bw_run = canon$bandwidth,
    late_run = canon$late_estimate, source = basename(dirname(cp))
  )
})

write_csv(check, file.path(out_dir, "sanity_check_rdd.csv"))

cat("\n  Reduced-form RD on GDP per capita, outcome measured from the year before the election:\n\n")
print(
  check |>
    transmute(
      window,
      `N (table)` = n_table, `N (run)` = n_run,
      `est (table)` = round(est_table, 5), `est (run)` = round(est_run, 5),
      `se (table)` = round(se_table, 5), `se (run)` = round(se_run, 5),
      `LATE bw (table)` = round(bw_late_table, 3), `bw col (run)` = round(bw_run, 3),
      `RF bw (table)` = round(bw_rf_table, 3)
    ) |>
    as.data.frame(),
  row.names = FALSE
)
cat(paste(
  "\n  Note: rdd_results.csv's 'bandwidth' column holds the fuzzy LATE's",
  "bandwidth, not the\n  reduced form's, so it is checked against LATE bw.",
  "The reduced form's own bandwidth\n  is shown last and is not recorded in",
  "the run folder at all.\n"
))

matched <- check |> filter(!is.na(est_run))
if (nrow(matched) == 0) {
  warning(
    "No canonical run folder to check against. Run 14_window_sweep.R with ",
    "SWEEP_ILLIBERAL_CUTOFF = SWEEP_OTHER_CUTOFF_MAX = \"popucut\" first.",
    call. = FALSE
  )
} else {
  worst <- with(matched, max(c(
    abs(est_table - est_run), abs(se_table - se_run),
    abs(late_table - late_run), abs(bw_late_table - bw_run)
  )))
  stopifnot(all(matched$n_table == matched$n_run))
  stopifnot(worst < 1e-8)
  cat(sprintf(
    "\n  PASS -- matches %s at %d horizon(s); largest discrepancy %.2e.\n",
    paste(unique(matched$source), collapse = ", "), nrow(matched), worst
  ))
}

message("\nTable written to ", out_dir)
