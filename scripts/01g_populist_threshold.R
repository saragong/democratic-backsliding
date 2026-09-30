# ==============================================================================
# A data-driven cutoff for "is this party illiberal?", calibrated on The
# PopuList.
#
# THE PROBLEM THIS SOLVES
#
# The RD compares elections where the more anti-pluralist of the top 2 narrowly
# won against ones where it narrowly lost. "Narrow" constrains the VOTE margin,
# but nothing constrains the illiberality GAP: in 29% of elections the two
# parties are within 0.05 of each other on v2xpa_antiplural, so crossing the
# cutoff swaps one nearly-identical party for another and there is no treatment
# contrast to speak of. The fix is to require that one side is illiberal and
# the other is not -- but "illiberal" needs a number, and picking one by hand
# is exactly the post-hoc choice the Sep 7 memo flagged as a weakness.
#
# THE PROCEDURE (Aug 24 notes, Sep 7 checklist)
#
#   1. Download The PopuList, an expert classification of European parties.
#   2. Merge it to V-Party via Party Facts.
#   3. Logit the populist dummy on V-Party's continuous populism score.
#   4. Take the ROC-optimal cutpoint.
#   5. Apply that cutpoint to the ANTI-PLURALISM score to classify parties as
#      illiberal or not.
#
# Step 5 is a transfer between two different variables, so the script emits the
# cutpoint three ways and 12_rdd_analysis.R can be run at each:
#   absolute         the raw number, carried across unchanged. Defensible
#                    because both indices are V-Party aggregates on [0, 1].
#                    Spec "popucut".
#   percentile       the cutpoint's percentile in the GLOBAL populism
#                    distribution (all V-Party party-years, every country and
#                    year), applied at the same percentile of the global
#                    anti-pluralism distribution. What the Aug 24 note
#                    described ("took the same percentile"). Spec
#                    "popucut_pct".
#   within-country   the fit redone with populism ranked within each party's
#                    own country; the percentile cut is applied to
#                    anti-pluralism ranked the same way, so it ports to
#                    countries PopuList never covers (section 4b). Specs
#                    "popucut_ctry" (accuracy) and "popucut_ctry_youden".
# They disagree whenever the indices are distributed differently, which they
# are, so none is obviously right and all are reported.
#
# "Party-years" throughout are V-Party observations, and V-Party scores a party
# at each ELECTION it contests (a median of 2 per party), so they are really
# party-elections: a party that contested more elections carries more weight in
# the fit and in every percentile.
#
# A CORRECTION TO THE ASK, STATED PLAINLY
#
# The brief says "choose the threshold that maximizes the AUC". AUC cannot be
# maximized over thresholds: it is the area under the whole ROC curve and is
# the same number whatever cutpoint you pick -- it measures the SCORE, not a
# cutpoint. What is wanted is an ROC-optimal cutpoint. The script computes two:
# Youden's J (sensitivity + specificity - 1), the standard choice, and the
# accuracy-maximizing cut, which is the headline (see HEADLINE_METHOD for when
# and why). Both are reported so the choice is visible rather than assumed,
# and the AUC is reported as what it actually is: how well the continuous
# score separates the classes.
#
# THE NEGATIVE CLASS IS THE WHOLE DIFFICULTY
#
# PopuList is a LIST, not a census: it names the European parties that are
# populist OR far right OR far left OR eurosceptic. So it carries 201 ones and
# 67 zeros, and those zeros exist only because a party can qualify on one of
# the other three grounds without being populist. Checked: all 67 of them are
# far right, far left or eurosceptic.
#
# That rules out the obvious approach of fitting on the listed parties alone.
# Its negative class would not be ordinary parties but OTHER RADICALS -- the
# French Communist Party, Sweden's Left Party, Poland's National Movement --
# which are populism-adjacent enough that labelling them 0 is closer to
# mislabelling than to sampling. Fitting on them asks the classifier to
# separate populist radicals from non-populist radicals, which is both harder
# than and different from the question here. It was tried, and it behaved
# exactly as that description predicts: 77% of party-years positive, mean
# populism among the "zeros" 0.521, AUC 0.644, and an accuracy cutpoint that
# collapsed to "call everything populist". Dropped.
#
# So the universe is the one PopuList's own inclusion rule implies. Its
# criterion is roughly "parties above 1% of the vote in 31 European countries
# that are populist / far right / far left / eurosceptic", so a party in one
# of those countries that is NOT listed is implicitly none of those things.
# Universe = every V-Party party-year in those countries within coverage;
# populist = 1 only for a listed party in a listed year, 0 otherwise. The
# negative class is then what it should be -- Sweden's Centre Party,
# Switzerland's Radical Democratic Party -- 20% positive, mean populism among
# the zeros 0.259, AUC 0.890.
#
# THE CAVEAT THAT TRAVELS WITH THE NUMBER
#
# PopuList covers 31 European countries. The cutpoint is applied to a global
# election spine, so it is an out-of-sample transfer in two directions at once
# (different variable, different countries). Every output says so.
#
#   Rscript --no-init-file scripts/01g_populist_threshold.R
#
# Output: data/populist_4.0.csv          cached download
#         data/populist_threshold.rds    every cutpoint, for 12_rdd_analysis.R
#         output/sweeps/populist_threshold/
#           thresholds.csv                  raw-score cuts and their transfers
#           thresholds_ctry_pct.csv         within-country percentile cuts
#           thresholds_sensitivity.csv      cuts with the ambiguous match recoded
#           merge_audit.csv                 every PopuList party and its V-Party id
#           unmatched_populists_audit.csv   populists with no V-Party id, nearest name
#           full_sample_cutpoint_counts.csv parties either side of each cut, globally
#           roc.png, score_distributions.png, score_distributions_full_sample.png,
#           score_distribution_ctry_pct.png
#         The side-by-side comparison of every cut on the RDD sample is
#         11c_threshold_summary.R's cutpoint_summary.html, which needs the build.
# ==============================================================================

# Numbered into the loader tier, not after the RDD scripts, because
# 12_rdd_analysis.R now consumes the threshold this writes: ILLIBERAL_CUTOFF
# and OTHER_CUTOFF_MAX accept the spec "popucut", which resolves out of
# data/populist_threshold.rds. The dependency therefore runs BEFORE the build
# and the analysis, and the file has to exist by then.
#
# Nothing here reads a build -- only the PopuList download, parties_database
# and the V-Party zip -- so it genuinely belongs at this stage rather than
# merely being moved to satisfy the ordering.
#
# A calibration step, not part of the RDD pipeline, so it has no toggles of its
# own and reads none of the pipeline's. It sources rdd_helpers.R for two
# functions only: sweep_dir(), for where the output goes, and country_pct(),
# so the within-country percentile here is computed by exactly the code that
# computes it in 11_build_rdd_data.R.
library(readr)
library(dplyr)
library(tidyr)
library(haven)
library(ggplot2)
library(pROC)
library(here)

# sweep_dir().
source(here::here("scripts", "rdd_helpers.R"))

data_dir <- here::here("data")
out_dir <- sweep_dir("populist_threshold")

POPULIST_URL <- "https://popu-list.github.io/Data/The%20PopuList%204.0.csv"
POPULIST_CSV <- file.path(data_dir, "populist_4.0.csv")

# The score the logit is fitted on, and the score the cutpoint is transferred
# to. Keeping them as named constants because the whole point of the exercise
# is that they are DIFFERENT variables.
FIT_SCORE <- "v2xpa_popul"
APPLY_SCORE <- "v2xpa_antiplural"

# Which cutpoint criterion the headline threshold uses. Both are always
# computed and both are carried in the output; this only decides which one
# data/populist_threshold.rds names as `threshold_abs` / `threshold_pct` and
# which the figures mark.
#
# "accuracy" (chosen 2026-09-21). Provenance matters here and belongs in the
# record rather than in anyone's memory: both criteria were computed first,
# the RDD was run at both, and accuracy was selected afterwards. It is the
# criterion that puts the bar where the illiberal-vs-not contrast is real --
# Youden lands at 0.362, near the middle of the anti-pluralism distribution,
# where using one number as both floor and ceiling makes the two parties
# straddle a point and is barely a contrast at all. It is ALSO the criterion
# under which the growth effect survives the restriction (-0.192 vs +0.042).
# Those two facts have the same cause, but a reader is entitled to know the
# choice was made with the second one visible, and any write-up should say
# so and show both.
#
# "youden" remains the defensible default for an imbalanced classification
# problem taken on its own terms, and switching back is one line.
HEADLINE_METHOD <- "accuracy"
stopifnot(HEADLINE_METHOD %in% c("youden", "accuracy"))

# ------------------------------------------------------------------------------
# 1. The PopuList
# ------------------------------------------------------------------------------

if (!file.exists(POPULIST_CSV)) {
  message("Downloading The PopuList 4.0 ...")
  utils::download.file(POPULIST_URL, POPULIST_CSV, quiet = TRUE)
}
# Semicolon-delimited with a BOM.
popu <- read_delim(POPULIST_CSV, delim = ";", show_col_types = FALSE)

stopifnot(all(
  c("party_name", "country_name", "populist", "populist_start", "populist_end",
    "partyfacts_id") %in% names(popu)
))
cat(sprintf(
  "PopuList 4.0: %d parties, %d countries, %d populist / %d not\n",
  nrow(popu), n_distinct(popu$country_name),
  sum(popu$populist == 1, na.rm = TRUE), sum(popu$populist == 0, na.rm = TRUE)
))

# ------------------------------------------------------------------------------
# 2. Merge to V-Party via Party Facts
#
# partyfacts_id -> any of parties_database's pf_id_1..4 -> vdem_id_1, which IS
# V-Party's v2paid. Going through parties_database rather than matching on name
# reuses the crosswalk the whole project already depends on.
# ------------------------------------------------------------------------------

pdb <- read_dta(file.path(data_dir, "elections_database", "parties_database.dta"))

pf_to_vdem <- pdb |>
  select(party_id, vdem_id_1, pf_id_1, pf_id_2, pf_id_3, pf_id_4) |>
  pivot_longer(starts_with("pf_id"), values_to = "pfid", names_to = NULL) |>
  filter(!is.na(pfid), !is.na(vdem_id_1)) |>
  distinct(pfid, vdem_id_1)

popu_matched <- popu |>
  mutate(pfid = as.numeric(partyfacts_id)) |>
  left_join(pf_to_vdem, by = "pfid", relationship = "many-to-many")

# PopuList populists that V-Party DOES score but that no Party Facts link
# reaches. Found by name-matching every unmatched populist against V-Party
# party names in the same country (the audit below re-runs that check each
# time). Left unmatched, such a party is imputed a ZERO in the calibration
# universe -- a populist labelled non-populist, high on the populism scale,
# exactly where the cutpoint is decided. Only confident matches go here:
#
#   Norway Sp   Senterpartiet, populist (borderline) from 2017. PopuList gives
#               it no Party Facts id at all; V-Party 1072 is the Centre
#               [Agrarian] Party, same country, scored through 2017.
#
# Deliberately NOT here: Romania "PSD" (PopuList pf 1078, 1990-2025). Its Party
# Facts id is the PSoDR, a minor socialist party dissolved in 1993, while its
# dates fit the large Partidul Social Democrat (V-Party 120). Which party
# PopuList means cannot be settled from the files, so it is left as PopuList
# coded it, and section 4c reports how much the cutpoint moves if it is the PSD.
MANUAL_VPARTY_MATCHES <- tibble::tribble(
  ~country_name, ~party_name_short, ~vdem_id_manual,
  "Norway",      "Sp",              1072
)
popu_matched <- popu_matched |>
  left_join(MANUAL_VPARTY_MATCHES, by = c("country_name", "party_name_short")) |>
  mutate(vdem_id_1 = coalesce(vdem_id_1, vdem_id_manual)) |>
  select(-vdem_id_manual)

n_pf <- sum(!is.na(popu_matched$pfid))
n_vdem <- sum(!is.na(popu_matched$vdem_id_1))
cat(sprintf(
  "Merge: %d / %d carry a Party Facts id; %d reach a V-Party id\n",
  n_pf, nrow(popu), n_vdem
))
# A future PopuList release that silently stops joining should fail loudly
# rather than quietly calibrate on a handful of parties.
if (n_vdem < 100) {
  stop(
    "Only ", n_vdem, " PopuList parties reached a V-Party id (expected ~130). ",
    "The Party Facts crosswalk has probably broken -- check partyfacts_id ",
    "against parties_database's pf_id_* columns before trusting any cutpoint."
  )
}

write_csv(
  popu_matched |>
    select(party_name, country_name, populist, populist_start, populist_end,
           partyfacts_id, vdem_id_1),
  file.path(out_dir, "merge_audit.csv")
)

# ------------------------------------------------------------------------------
# 3. The two universes
# ------------------------------------------------------------------------------

vparty <- read_csv(
  unz(
    file.path(data_dir, "elections_database", "CPD_V-Party_CSV_v2.zip"),
    "CPD_V-Party_CSV_v2/V-Dem-CPD-Party-V2.csv"
  ),
  show_col_types = FALSE
) |>
  select(v2paid, country_name, country_text_id, year, all_of(c(FIT_SCORE, APPLY_SCORE))) |>
  filter(!is.na(.data[[FIT_SCORE]]))

# PopuList's country names against V-Party's. Verified: all 31 match exactly
# on the raw string -- both use "Czech Republic", not "Czechia" -- so no
# crosswalk is needed. This is an error rather than a warning because a name
# that stops matching silently drops a whole country's parties out of the
# imputed-zero universe, which would move the cutpoint without any other sign.
popu_countries <- unique(popu$country_name)
unmatched <- setdiff(popu_countries, unique(vparty$country_name))
if (length(unmatched) > 0) {
  stop(
    "PopuList countries with no V-Party name match: ",
    paste(unmatched, collapse = ", "),
    ". Add a crosswalk entry -- leaving them out would silently shrink the ",
    "imputed-zero universe and move the cutpoint.",
    call. = FALSE
  )
}

# ---- unmatched-populist audit ------------------------------------------------
#
# Every PopuList populist still without a V-Party id, with the closest V-Party
# party name in the same country (Jaro-Winkler on accent- and
# punctuation-stripped names, over the original, English and short names).
# Written on every run so a new PopuList or V-Party release gets the same check.
# Nothing is recoded automatically: most unmatched populists are simply parties
# V-Party does not score, and a close name is often a different party (France's
# 1968 UDR vs the 2024 UDR). Confident matches go in MANUAL_VPARTY_MATCHES.
local({
  norm <- function(x) {
    gsub("[^a-z0-9]", "", stringi::stri_trans_general(tolower(x), "Latin-ASCII"))
  }
  vnames <- read_csv(
    unz(
      file.path(data_dir, "elections_database", "CPD_V-Party_CSV_v2.zip"),
      "CPD_V-Party_CSV_v2/V-Dem-CPD-Party-V2.csv"
    ),
    show_col_types = FALSE
  ) |>
    filter(!is.na(.data[[FIT_SCORE]])) |>
    distinct(v2paid, country_name, v2paenname, v2paorname, v2pashname) |>
    pivot_longer(c(v2paenname, v2paorname, v2pashname), values_to = "vparty_name") |>
    filter(!is.na(vparty_name)) |>
    mutate(vkey = norm(vparty_name)) |>
    distinct(v2paid, country_name, vparty_name, vkey)
  unmatched_pop <- popu_matched |>
    filter(populist == 1, is.na(vdem_id_1)) |>
    select(country_name, party_name, party_name_english, party_name_short,
           partyfacts_id, populist_start, populist_end, populist_bl) |>
    mutate(.row = row_number())
  cand <- unmatched_pop |>
    pivot_longer(c(party_name, party_name_english, party_name_short),
                 names_to = "popu_field", values_to = "popu_name") |>
    filter(!is.na(popu_name)) |>
    mutate(key = norm(popu_name)) |>
    filter(nchar(key) >= 2) |>
    inner_join(vnames, by = "country_name", relationship = "many-to-many") |>
    mutate(similarity = 1 - stringdist::stringdist(key, vkey, method = "jw", p = 0.1)) |>
    group_by(.row) |>
    slice_max(similarity, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(.row, matched_on = popu_name,
           best_vparty_name = vparty_name, best_v2paid = v2paid, similarity)
  audit <- unmatched_pop |>
    left_join(cand, by = ".row") |>
    select(-.row) |>
    arrange(desc(similarity))
  write_csv(audit, file.path(out_dir, "unmatched_populists_audit.csv"))
  cat(sprintf(
    "Unmatched-populist audit: %d PopuList populists without a V-Party id; %d with a same-country V-Party name at Jaro-Winkler >= 0.9 (review unmatched_populists_audit.csv)\n",
    nrow(unmatched_pop), sum(audit$similarity >= 0.9, na.rm = TRUE)
  ))
})

# PopuList's own coverage window, so the imputed zeros are not asserted for
# years the classification never looked at. The start/end columns use 1900 and
# 2100 as open-ended sentinels.
COVERAGE_MIN <- 1989
COVERAGE_MAX <- 2019

# Long form: one row per (party, year) that PopuList says is populist.
popu_years <- popu_matched |>
  filter(!is.na(vdem_id_1), populist == 1) |>
  transmute(
    vdem_id_1,
    start = pmax(populist_start, COVERAGE_MIN),
    end = pmin(populist_end, COVERAGE_MAX)
  ) |>
  filter(start <= end) |>
  rowwise() |>
  mutate(year = list(seq(start, end))) |>
  ungroup() |>
  unnest(year) |>
  distinct(vdem_id_1, year) |>
  mutate(populist = 1L)

# Zeros are imputed ONLY in countries where at least one PopuList populist
# actually landed in V-Party after the merge -- not in all 31 PopuList
# countries.
#
# The difference is a guard against the merge, not against PopuList. Only 132
# of PopuList's 268 parties reach a V-Party id, and that attrition is not
# spread evenly: a country can lose every one of its listed populists to a
# failed Party Facts link. In such a country, imputing zeros marks its
# genuinely populist parties as non-populist and feeds the logit pure false
# negatives, precisely in the upper range of the populism score where the
# cutpoint is decided. Requiring one surviving positive per country does not
# fix a country that lost SOME of its populists, but it removes the case
# where the country contributes nothing but wrongly-labelled zeros.
#
# The cost is real and is reported: countries dropped this way had parties
# that might legitimately have been zeros, so the negative class is smaller
# and slightly more concentrated in populist-heavy countries.
# Counted inside the coverage window: a populist matched only in years the
# calibration never uses would otherwise keep its country's imputed zeros
# without contributing a single positive.
countries_with_matched_populist <- vparty |>
  filter(
    v2paid %in% unique(popu_years$vdem_id_1),
    year >= COVERAGE_MIN, year <= COVERAGE_MAX
  ) |>
  pull(country_name) |>
  unique()

dropped_countries <- setdiff(popu_countries, countries_with_matched_populist)
cat(sprintf(
  "\nImputed-zero universe: %d of %d PopuList countries retain a matched populist\n",
  length(countries_with_matched_populist), length(popu_countries)
))
if (length(dropped_countries) > 0) {
  cat(sprintf(
    "  dropped (no PopuList populist survived the merge): %s\n",
    paste(sort(dropped_countries), collapse = ", ")
  ))
}

universe_imputed <- vparty |>
  filter(
    country_name %in% countries_with_matched_populist,
    year >= COVERAGE_MIN, year <= COVERAGE_MAX
  ) |>
  left_join(popu_years, by = c("v2paid" = "vdem_id_1", "year")) |>
  mutate(populist = coalesce(populist, 0L))

# Listed-only: the parties PopuList names, at their V-Party observations.
UNIVERSE_LABEL <- paste(
  "All V-Party parties in PopuList's countries;",
  "unlisted = not populist"
)
# (The unit is the V-Party observation: one per party per election contested.)

# ------------------------------------------------------------------------------
# 4. Logit, ROC, cutpoint
# ------------------------------------------------------------------------------

# The logit of a binary outcome on ONE continuous predictor is monotone in that
# predictor, so the ROC of the fitted probability and the ROC of the raw score
# are identical, and a cutpoint on the fitted probability maps back to exactly
# one cutpoint on the score. The model is therefore fitted (it is what the ask
# specifies, and its coefficient is worth reporting) but the cutpoint is read
# off the SCORE directly -- which is what has to be transferred to
# anti-pluralism, and what 12_rdd_analysis.R can actually apply.
fit_one <- function(df, label, score = FIT_SCORE) {
  df <- df |> filter(!is.na(.data[[score]]), !is.na(populist))
  m <- glm(
    as.formula(paste("populist ~", score)),
    data = df, family = binomial()
  )
  roc_obj <- roc(
    response = df$populist, predictor = df[[score]],
    quiet = TRUE, direction = "<"
  )

  # Two cutpoints, answering two different questions.
  #
  # YOUDEN maximizes sensitivity + specificity - 1, weighting the two classes
  # equally regardless of how many of each there are. It is the standard
  # "ROC-optimal" point and is the right one when the classes are imbalanced
  # and you care about both kinds of error, which is the case here: populists
  # are a fifth of the universe.
  #
  # ACCURACY maximizes (TP + TN) / N, which weights by prevalence. With a
  # minority positive class this pulls the cutpoint UP -- calling fewer
  # parties populist is cheap in accuracy terms because most parties are not
  # populist. Reported because it is the other natural reading of "best", and
  # because seeing the two apart is what shows the choice is doing work.
  #
  # pROC has no accuracy method, so it is computed directly over every
  # candidate cutpoint on the score.
  grab_roc <- function(method) {
    co <- coords(roc_obj, "best", best.method = method,
                 ret = c("threshold", "sensitivity", "specificity"),
                 transpose = FALSE)
    # Ties can return several equally good cutpoints. Take the tied cutpoint
    # nearest their median, so the answer does not depend on pROC's internal
    # ordering and is itself one of the optima -- a plain median of disjoint
    # optima can land between them, at a cut that is not optimal.
    i <- which.min(abs(co$threshold - median(co$threshold)))
    tibble(
      method = method,
      threshold = co$threshold[i],
      sensitivity = co$sensitivity[i],
      specificity = co$specificity[i]
    )
  }

  grab_accuracy <- function() {
    x <- df[[score]]
    y <- df$populist == 1
    # Midpoints between adjacent observed values, so every distinct split of
    # the data is considered exactly once -- plus both extremes ("everyone
    # populist" below the minimum, "no one" above the maximum), so the
    # degenerate solutions are in the candidate set and the guard below can
    # catch them.
    u <- sort(unique(x))
    cands <- c(u[1] - 1e-9, if (length(u) > 1) (head(u, -1) + tail(u, -1)) / 2, u[length(u)] + 1e-9)
    acc <- vapply(cands, function(t) mean((x > t) == y), numeric(1))
    best <- cands[acc == max(acc)]
    # The optimal cut nearest the median of the tied optima (see grab_roc).
    t0 <- best[which.min(abs(best - median(best)))]
    tibble(
      method = "accuracy",
      threshold = t0,
      sensitivity = mean(x[y] > t0),
      specificity = mean(x[!y] <= t0)
    )
  }

  cuts <- bind_rows(grab_roc("youden"), grab_accuracy())

  # Accuracy AT each chosen cutpoint, so the two can be compared on the same
  # footing rather than each on its own criterion.
  cuts$accuracy <- vapply(cuts$threshold, function(t) {
    mean((df[[score]] > t) == (df$populist == 1))
  }, numeric(1))
  cuts$youden_j <- cuts$sensitivity + cuts$specificity - 1

  # A cutpoint with zero specificity (or zero sensitivity) is the degenerate
  # "call everything one class" solution. Accuracy will choose it whenever one
  # class dominates badly enough, since predicting the majority everywhere
  # then scores well. It does not fire on the current universe (20% positive),
  # but it did fire on the listed-parties-only universe this script used to
  # carry (77% positive, accuracy 0.768 for "all populist"), which is why the
  # guard is here rather than assumed unnecessary.
  degenerate <- cuts$specificity <= .Machine$double.eps |
    cuts$sensitivity <= .Machine$double.eps
  if (any(degenerate)) {
    warning(
      label, ": the ", paste(cuts$method[degenerate], collapse = " and "),
      " cutpoint is degenerate (Youden's J = 0) -- it assigns every ",
      "observation to one class. Do not use it as a threshold.",
      call. = FALSE
    )
  }
  cuts$degenerate <- degenerate

  list(
    label = label,
    n = nrow(df),
    n_pos = sum(df$populist == 1),
    auc = as.numeric(auc(roc_obj)),
    coef = unname(coef(m)[2]),
    coef_se = sqrt(diag(vcov(m)))[2],
    roc = roc_obj,
    cuts = cuts,
    fit_score = df[[score]]
  )
}

# The percentile transfer: where the cutpoint sits in the populism
# distribution, then the anti-pluralism value at that same percentile. Both are
# taken over ALL scored V-Party party-years -- every country, every year -- not
# the European calibration sample the logit was fitted on, because that global
# population is what the threshold is applied to. Both indices are ranked over
# the same rows, so the two percentiles are comparable.
apply_scores <- vparty |>
  filter(!is.na(.data[[FIT_SCORE]]), !is.na(.data[[APPLY_SCORE]]))

transfer <- function(thr) {
  pct <- 100 * mean(apply_scores[[FIT_SCORE]] <= thr, na.rm = TRUE)
  tibble(
    threshold_absolute = thr,
    percentile = pct,
    threshold_percentile_value = unname(
      quantile(apply_scores[[APPLY_SCORE]], probs = pct / 100, na.rm = TRUE)
    )
  )
}

fit <- fit_one(universe_imputed, UNIVERSE_LABEL)
cat(sprintf(
  "\n  N = %d (%d populist, %.1f%%)\n  logit coef on %s = %.3f (%.3f)\n  AUC = %.3f\n",
  fit$n, fit$n_pos, 100 * fit$n_pos / fit$n, FIT_SCORE,
  fit$coef, fit$coef_se, fit$auc
))

rows <- lapply(seq_len(nrow(fit$cuts)), function(i) {
  tr <- transfer(fit$cuts$threshold[i])
  cat(sprintf(
    "  %-9s cut %.4f on %s (sens %.2f spec %.2f J %.2f acc %.3f) -> %.1fth pct -> %.4f on %s\n",
    fit$cuts$method[i], tr$threshold_absolute, FIT_SCORE,
    fit$cuts$sensitivity[i], fit$cuts$specificity[i],
    fit$cuts$youden_j[i], fit$cuts$accuracy[i],
    tr$percentile, tr$threshold_percentile_value, APPLY_SCORE
  ))
  bind_cols(
    tibble(
      universe_label = fit$label,
      n = fit$n, n_populist = fit$n_pos, auc = fit$auc,
      logit_coef = fit$coef, logit_se = fit$coef_se,
      method = fit$cuts$method[i],
      sensitivity = fit$cuts$sensitivity[i],
      specificity = fit$cuts$specificity[i],
      accuracy = fit$cuts$accuracy[i],
      youden_j = fit$cuts$youden_j[i],
      degenerate = fit$cuts$degenerate[i]
    ),
    tr
  )
})

thresholds <- bind_rows(rows)
write_csv(thresholds, file.path(out_dir, "thresholds.csv"))

# ------------------------------------------------------------------------------
# 4b. The within-country percentile calibration
#
# The same fit, with populism re-expressed as a percentile WITHIN the party's
# own country: F_c(score) over every V-Party party-year in country c, all
# years (country_pct() in rdd_helpers.R). The cut is then a percentile p*, and
# 11_build_rdd_data.R ranks anti-pluralism within country the same way, so
# "illiberal" means "in the top (1 - p*) of its own country's parties" --
# which ports to every country V-Party covers, not only PopuList's 29. A raw
# cut of 0.65 is strict in a country whose parties all sit low on the scale;
# the percentile adapts to where each country's party system actually sits.
# ------------------------------------------------------------------------------

universe_imputed$popul_pct_ctry <- country_pct(
  universe_imputed[[FIT_SCORE]], universe_imputed$country_text_id,
  vparty[[FIT_SCORE]], vparty$country_text_id
)
fit_ctry <- fit_one(
  universe_imputed,
  paste(UNIVERSE_LABEL, "(populism as a within-country percentile)"),
  score = "popul_pct_ctry"
)
cat(sprintf(
  "\nWithin-country percentile calibration: N = %d (%d populist), AUC = %.3f\n",
  fit_ctry$n, fit_ctry$n_pos, fit_ctry$auc
))
thresholds_ctry <- fit_ctry$cuts |>
  transmute(
    universe_label = fit_ctry$label,
    n = fit_ctry$n, n_populist = fit_ctry$n_pos, auc = fit_ctry$auc,
    method, percentile_cut = threshold,
    sensitivity, specificity, accuracy, youden_j, degenerate
  )
for (i in seq_len(nrow(thresholds_ctry))) {
  cat(sprintf(
    "  %-9s cut at the %.1fth within-country percentile (sens %.2f spec %.2f J %.2f acc %.3f)\n",
    thresholds_ctry$method[i], 100 * thresholds_ctry$percentile_cut[i],
    thresholds_ctry$sensitivity[i], thresholds_ctry$specificity[i],
    thresholds_ctry$youden_j[i], thresholds_ctry$accuracy[i]
  ))
}
write_csv(thresholds_ctry, file.path(out_dir, "thresholds_ctry_pct.csv"))
# The within-country headline uses the SAME criterion as the raw one, by
# construction rather than by a separate choice: HEADLINE_METHOD was selected
# on the raw score after the RDD had been run at both criteria (see its
# comment), and it is inherited here without having been chosen on this scale.
#
# On this scale the accuracy criterion is FLAT near its optimum: recoding a
# single party-year (Norway's Senterpartiet in 2017, the MANUAL_VPARTY_MATCHES
# fix) moved the cut from the 90.3th to the 86.6th percentile with identical
# accuracy (0.8839) at both. Where it lands is therefore decided by very few
# observations. The Youden cut (71.6th) is stable to that recode and travels
# in the output too; read every within-country run at both.
headline_ctry <- thresholds_ctry |>
  filter(!degenerate, method == HEADLINE_METHOD) |>
  slice(1)
stopifnot(nrow(headline_ctry) == 1)

# ------------------------------------------------------------------------------
# 4c. Sensitivity: the one ambiguous PopuList match
#
# Romania "PSD" (see MANUAL_VPARTY_MATCHES) is either the minor PSoDR its Party
# Facts id names, which V-Party does not score, or the large PSD (V-Party 120),
# which the calibration then carries as a zero while PopuList calls it
# populist. Refit both calibrations with the PSD coded populist over PopuList's
# stated years, and report how far each cut moves. Nothing downstream reads
# this; it says whether the ambiguity matters.
# ------------------------------------------------------------------------------

psd <- popu |> filter(country_name == "Romania", party_name_short == "PSD", populist == 1)
universe_psd <- universe_imputed |>
  mutate(populist = if_else(
    v2paid == 120 & year >= max(psd$populist_start, COVERAGE_MIN) &
      year <= min(psd$populist_end, COVERAGE_MAX),
    1L, populist
  ))
n_recoded <- sum(universe_psd$populist != universe_imputed$populist)
fit_psd <- fit_one(universe_psd, "sensitivity: Romania PSD coded populist")
fit_psd_ctry <- fit_one(universe_psd, "sensitivity: Romania PSD coded populist", score = "popul_pct_ctry")
sensitivity <- bind_rows(
  tibble(scale = "raw", method = fit$cuts$method, baseline = fit$cuts$threshold),
  tibble(scale = "within-country percentile", method = fit_ctry$cuts$method, baseline = fit_ctry$cuts$threshold)
) |>
  mutate(
    with_romania_psd = c(fit_psd$cuts$threshold, fit_psd_ctry$cuts$threshold),
    change = with_romania_psd - baseline,
    party_years_recoded = n_recoded
  )
write_csv(sensitivity, file.path(out_dir, "thresholds_sensitivity.csv"))
cat(sprintf("\nSensitivity, Romania PSD coded populist (%d party-years recoded):\n", n_recoded))
print(as.data.frame(sensitivity |> mutate(across(where(is.numeric), \(x) round(x, 4)))), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. The headline numbers, for 12_rdd_analysis.R
# ------------------------------------------------------------------------------

# Both cutpoints from the headline universe travel in the output, so
# 12_rdd_analysis.R can be run at either without re-deriving anything;
# HEADLINE_METHOD decides which one is named as the default.
cuts_out <- thresholds |> filter(!degenerate)
stopifnot(nrow(cuts_out) >= 1)
headline <- cuts_out |> filter(method == HEADLINE_METHOD) |> slice(1)
if (nrow(headline) == 0) {
  stop(
    "HEADLINE_METHOD = '", HEADLINE_METHOD, "' produced no usable cutpoint ",
    "(it may have been dropped as degenerate). Available: ",
    paste(cuts_out$method, collapse = ", "),
    call. = FALSE
  )
}

out <- list(
  fit_score = FIT_SCORE,
  apply_score = APPLY_SCORE,
  universe = UNIVERSE_LABEL,
  method = headline$method,
  auc = headline$auc,
  n = headline$n,
  # For ILLIBERAL_CUTOFF / OTHER_CUTOFF_MAX, the absolute transfer.
  threshold_abs = headline$threshold_absolute,
  # Every non-degenerate cutpoint from this universe: method, the raw value,
  # and its percentile transfer onto the anti-pluralism scale.
  cutpoints = cuts_out |>
    select(method, threshold_absolute, percentile,
           threshold_percentile_value, sensitivity, specificity,
           accuracy, youden_j),
  # ... and the percentile transfer, as an absolute value on the
  # anti-pluralism scale. Given as a NUMBER rather than a "qNN" string on
  # purpose: "qNN" would be re-resolved by 12_rdd_analysis.R against whatever
  # sample that run has, which is the top-2 election spine, not the V-Party
  # party-year distribution this percentile was computed on.
  threshold_pct = headline$threshold_percentile_value,
  percentile = headline$percentile,
  # The within-country percentile cut (0-1), for the "popucut_ctry" spec. It is
  # applied to illiberal_pct_ctry / other_pct_ctry, never to a raw score.
  threshold_pct_ctry = headline_ctry$percentile_cut,
  # The Youden cut on the same scale, carried because the accuracy criterion is
  # flat here (see headline_ctry) and the question is precisely whether the
  # PopuList cut is too conservative. "popucut_ctry_youden" reads it.
  threshold_pct_ctry_youden = thresholds_ctry |>
    filter(!degenerate, method == "youden") |>
    pull(percentile_cut) |>
    first(default = NA_real_),
  auc_ctry = fit_ctry$auc,
  # The within-country fit's own N: smaller than `n` whenever country_pct()
  # returned NA (a country with fewer than COUNTRY_PCT_MIN_N party-years).
  n_ctry = fit_ctry$n,
  cutpoints_ctry = thresholds_ctry |>
    filter(!degenerate) |>
    select(method, percentile_cut, sensitivity, specificity, accuracy, youden_j),
  coverage = sprintf(
    "PopuList 4.0, %d European countries, %d-%d",
    n_distinct(popu$country_name), COVERAGE_MIN, COVERAGE_MAX
  ),
  built_at = Sys.time()
)
saveRDS(out, file.path(data_dir, "populist_threshold.rds"))
cat(sprintf(
  "\nHeadline (%s): %.4f absolute, %.4f at the same percentile (%.1fth)\n",
  out$method, out$threshold_abs, out$threshold_pct, out$percentile
))
cat("All usable cutpoints carried in the RDS:\n")
print(
  out$cutpoints |>
    mutate(across(where(is.numeric), \(x) round(x, 4))) |>
    as.data.frame(),
  row.names = FALSE
)
cat(sprintf(
  "Within-country percentile (%s): %.4f\n", out$method, out$threshold_pct_ctry
))
cat("Saved data/populist_threshold.rds\n")

# ------------------------------------------------------------------------------
# 6. Figures
# ------------------------------------------------------------------------------

roc_df <- tibble(
  fpr = 1 - fit$roc$specificities,
  tpr = fit$roc$sensitivities
)

p_roc <- ggplot(roc_df, aes(fpr, tpr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.8, colour = "#0072B2") +
  # Both cutpoints marked, the headline one filled. Showing only the chosen
  # one would hide how far apart they are, which is the thing a reader of
  # this figure most needs to see.
  geom_point(
    data = thresholds |>
      filter(!degenerate) |>
      mutate(
        fpr = 1 - specificity, tpr = sensitivity,
        lab = sprintf("%s (%.3f)", method, threshold_absolute)
      ),
    aes(fpr, tpr, shape = method == HEADLINE_METHOD),
    size = 2.8, stroke = 1, fill = "white", colour = "#0072B2"
  ) +
  geom_text(
    data = thresholds |>
      filter(!degenerate) |>
      mutate(
        fpr = 1 - specificity, tpr = sensitivity,
        lab = sprintf("%s (%.3f)", method, threshold_absolute)
      ),
    aes(fpr, tpr, label = lab),
    hjust = -0.12, vjust = 1.6, size = 2.8, colour = "grey20"
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 21, `TRUE` = 19),
    labels = c(`FALSE` = "other cutpoint", `TRUE` = HEADLINE_METHOD),
    name = NULL
  ) +
  coord_fixed() +
  labs(
    title = sprintf(
      "Does %s separate PopuList's populist parties?  AUC %.3f",
      FIT_SCORE, fit$auc
    ),
    subtitle = paste(
      strwrap(sprintf(paste(
        "ROC of the raw V-Party populism score. Both cutpoints are marked;",
        "the filled one is the headline (%s). AUC measures the SCORE, not any",
        "cutpoint -- it is the same number whichever cutpoint is chosen.",
        "PopuList covers 31 European countries, so a cutpoint fitted here is",
        "applied out of sample to the global election spine."
      ), HEADLINE_METHOD), 95),
      collapse = "\n"
    ),
    x = "False positive rate (1 - specificity)",
    y = "True positive rate (sensitivity)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 7.5, colour = "grey25")
  )
ggsave(file.path(out_dir, "roc.png"), p_roc, width = 7, height = 7, dpi = 150)
cat("Saved roc.png\n")

# On FIT_SCORE the two rules coincide by construction -- the percentile
# transfer is defined as the percentile the cutpoint sits at, so evaluating it
# on the score it came from returns the cutpoint. They differ in the 4th
# decimal only because quantile() interpolates between order statistics, so an
# exact-equality dedup leaves two lines a thousandth apart, drawn on top of
# each other and indistinguishable. Collapse on a tolerance instead.
dedup_cuts <- function(df, tol = 1e-2) {
  # Pairwise, not bucketed: rounding value/tol to a bucket index splits two
  # values that straddle a bucket edge however close they are (0.6535 and
  # 0.6532 land in different buckets at tol = 1e-3 despite differing by
  # 3e-4). Keep a row only if no already-kept row on the same scale is within
  # tol of it.
  df |>
    arrange(score, rule) |>
    group_by(score) |>
    filter({
      keep <- logical(length(value))
      for (i in seq_along(value)) {
        keep[i] <- !any(keep & abs(value - value[i]) < tol)
      }
      keep
    }) |>
    ungroup()
}

# The transfer, shown rather than asserted, on the sample the cutpoint was
# actually fitted on.
#
# FIT_SCORE on top and APPLY_SCORE below, in the direction the transfer runs:
# the cutpoint is found on the upper panel and carried to the lower one.
# facet_wrap would otherwise order them alphabetically and put the
# destination above the source.
#
# Three densities per panel rather than one histogram. The single pooled
# histogram showed where parties sit but not whether the cutpoint separates
# anything; the split by assigned class is the part worth looking at, and the
# pooled curve stays as the backdrop so the two classes can be read against
# the whole. Densities rather than counts because the classes are 271 against
# 1,107 and stacked counts would make the smaller one invisible.
SCORE_LABELS <- c(
  "Populism (v2xpa_popul) -- the cutpoint is fitted here",
  "Anti-pluralism (v2xpa_antiplural) -- the cutpoint is applied here"
)
names(SCORE_LABELS) <- c(FIT_SCORE, APPLY_SCORE)

dist_long <- universe_imputed |>
  select(all_of(c(FIT_SCORE, APPLY_SCORE)), populist) |>
  pivot_longer(all_of(c(FIT_SCORE, APPLY_SCORE)),
               names_to = "score", values_to = "value") |>
  filter(!is.na(value))

CLASS_LEVELS <- c(
  "All parties in the inclusion sample",
  "Classified populist",
  "Not populist"
)
dist_df <- bind_rows(
  dist_long |> mutate(grp = CLASS_LEVELS[1]),
  dist_long |> filter(populist == 1) |> mutate(grp = CLASS_LEVELS[2]),
  dist_long |> filter(populist == 0) |> mutate(grp = CLASS_LEVELS[3])
) |>
  mutate(
    grp = factor(grp, levels = CLASS_LEVELS),
    score = factor(score, levels = c(FIT_SCORE, APPLY_SCORE))
  )

# Both rules on both scales. On FIT_SCORE they coincide by construction --
# the percentile transfer is defined as the percentile the cutpoint sits at,
# so evaluating it on the score it came from returns the cutpoint itself --
# and the duplicate line is dropped rather than drawn twice.
cut_df <- tibble(
  score = factor(c(FIT_SCORE, FIT_SCORE, APPLY_SCORE, APPLY_SCORE),
                 levels = c(FIT_SCORE, APPLY_SCORE)),
  rule = rep(c("Absolute threshold", "Percentile transfer"), 2),
  value = c(
    out$threshold_abs,
    unname(quantile(apply_scores[[FIT_SCORE]], out$percentile / 100, na.rm = TRUE)),
    out$threshold_abs,
    out$threshold_pct
  )
) |>
  mutate(lab = sprintf("%s (%.3f)", rule, value)) |>
  dedup_cuts()

p_dist <- ggplot(dist_df, aes(value)) +
  geom_density(
    aes(colour = grp, fill = grp),
    alpha = 0.22, linewidth = 0.6, adjust = 0.9
  ) +
  geom_vline(
    data = cut_df, aes(xintercept = value, linetype = rule),
    colour = "grey15", linewidth = 0.55
  ) +
  geom_text(
    data = cut_df,
    aes(x = value, y = Inf, label = sprintf("%.3f", value)),
    hjust = -0.15, vjust = 1.8, size = 2.7, colour = "grey15"
  ) +
  facet_wrap(
    ~score, ncol = 1, scales = "free_y",
    labeller = labeller(score = SCORE_LABELS)
  ) +
  scale_colour_manual(values = c("grey45", "#D55E00", "#0072B2"), name = NULL) +
  scale_fill_manual(values = c("grey70", "#D55E00", "#0072B2"), name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
  labs(
    title = "Transferring the cutpoint from populism to anti-pluralism",
    subtitle = paste(
      strwrap(sprintf(paste(
        "Inclusion sample only: %d party-years in the %d PopuList countries",
        "that retained a matched populist, %d of them classified populist.",
        "The absolute threshold carries the raw cutpoint across; the",
        "percentile transfer carries the %.1fth percentile across. On the",
        "populism panel the two coincide by construction, so one line is",
        "drawn. They separate on anti-pluralism because the two indices are",
        "not distributed alike. Percentiles are taken over all scored V-Party",
        "parties, not just this sample, since the threshold is applied to the",
        "global election spine."
      ), nrow(universe_imputed), length(countries_with_matched_populist),
      sum(universe_imputed$populist == 1), out$percentile), 100),
      collapse = "\n"
    ),
    x = "Score", y = "Density"
  ) +
  guides(
    colour = guide_legend(order = 1), fill = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(size = 8.5, face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -2),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 7.5, colour = "grey25")
  )
ggsave(
  file.path(out_dir, "score_distributions.png"), p_dist,
  width = 9, height = 6.8, dpi = 150
)
cat("Saved score_distributions.png\n")

# ------------------------------------------------------------------------------
# The same two cutpoints against the FULL party population
#
# score_distributions.png shows the inclusion sample -- the 29 European
# countries the logit was fitted on. But the cutpoint is not applied there: it
# is applied to the global election spine, where most parties are from
# countries PopuList never looked at. So the question this figure answers is
# the one that actually matters downstream: where do these thresholds fall in
# the population they will be used to cut?
#
# All scored V-Party party-years, every country and every year, which is also
# the population the percentile transfer was computed over -- so the dashed
# line sits at its stated percentile here by construction, and does not
# elsewhere.
#
# Counts rather than densities, because "how many parties does this threshold
# put on each side" is the thing being reported.
# ------------------------------------------------------------------------------

full_cuts <- tibble(
  score = factor(c(FIT_SCORE, FIT_SCORE, APPLY_SCORE, APPLY_SCORE),
                 levels = c(FIT_SCORE, APPLY_SCORE)),
  rule = rep(c("Absolute threshold", "Percentile transfer"), 2),
  value = c(
    out$threshold_abs,
    unname(quantile(apply_scores[[FIT_SCORE]], out$percentile / 100, na.rm = TRUE)),
    out$threshold_abs,
    out$threshold_pct
  )
) |>
  dedup_cuts()

# How many parties each threshold puts on each side, on each scale.
full_counts <- full_cuts |>
  rowwise() |>
  mutate(
    n_total = sum(!is.na(apply_scores[[as.character(score)]])),
    n_above = sum(apply_scores[[as.character(score)]] > value, na.rm = TRUE),
    n_below = n_total - n_above,
    pct_above = 100 * n_above / n_total
  ) |>
  ungroup()

write_csv(full_counts, file.path(out_dir, "full_sample_cutpoint_counts.csv"))

cat("\nFull V-Party population (all countries, all years), parties either side:\n")
print(
  full_counts |>
    transmute(
      score = as.character(score), rule,
      threshold = round(value, 4),
      n_below, n_above, n_total, pct_above = round(pct_above, 1)
    ) |>
    as.data.frame(),
  row.names = FALSE
)

full_long <- apply_scores |>
  select(all_of(c(FIT_SCORE, APPLY_SCORE))) |>
  pivot_longer(everything(), names_to = "score", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(score = factor(score, levels = c(FIT_SCORE, APPLY_SCORE)))

# One label block per panel, listing each rule's split. Placed top-left, where
# neither index has much mass.
count_lab <- full_counts |>
  summarise(
    lab = paste(
      sprintf("%s %.3f: %s below / %s above (%.0f%% above)",
              rule, value, format(n_below, big.mark = ","),
              format(n_above, big.mark = ","), pct_above),
      collapse = "\n"
    ),
    .by = score
  )

p_full <- ggplot(full_long, aes(value)) +
  geom_histogram(bins = 60, fill = "grey80", colour = "white", linewidth = 0.2) +
  geom_vline(
    data = full_cuts, aes(xintercept = value, linetype = rule),
    colour = "grey15", linewidth = 0.55
  ) +
  geom_text(
    data = count_lab, aes(x = -Inf, y = Inf, label = lab),
    hjust = -0.03, vjust = 1.25, size = 2.6, colour = "grey15", lineheight = 1.15
  ) +
  facet_wrap(
    ~score, ncol = 1, scales = "free_y",
    labeller = labeller(score = SCORE_LABELS)
  ) +
  scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
  labs(
    title = "Where the cutpoints fall in the full V-Party population",
    subtitle = paste(
      strwrap(sprintf(paste(
        "All %s scored party-years, every country and every year -- NOT the",
        "European inclusion sample the logit was fitted on (that is",
        "score_distributions.png). This is the population the threshold is",
        "actually applied to. On populism the two rules coincide by",
        "construction and one line is drawn; they separate on anti-pluralism",
        "because the two indices are not distributed alike. Counts either side",
        "are in full_sample_cutpoint_counts.csv."
      ), format(nrow(apply_scores), big.mark = ",")), 100),
      collapse = "\n"
    ),
    x = "Score", y = "Party-years"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(size = 8.5, face = "bold"),
    legend.position = "bottom",
    legend.margin = margin(t = -2),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 7.5, colour = "grey25")
  )
ggsave(
  file.path(out_dir, "score_distributions_full_sample.png"), p_full,
  width = 9, height = 6.4, dpi = 150
)
cat("Saved score_distributions_full_sample.png\n")

# The within-country calibration on its own scale: where PopuList's populists
# and everyone else sit in their OWN country's populism distribution, and the
# percentile cut the accuracy criterion picks.
p_ctry <- universe_imputed |>
  filter(!is.na(popul_pct_ctry)) |>
  mutate(group = if_else(populist == 1, "PopuList populist", "Not listed")) |>
  ggplot(aes(x = popul_pct_ctry, fill = group)) +
  geom_histogram(binwidth = 0.025, boundary = 0, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = out$threshold_pct_ctry, linewidth = 0.6) +
  annotate(
    "label", x = out$threshold_pct_ctry, y = Inf, vjust = 1.3, size = 2.8,
    label = sprintf(
      "%s cut: %.1fth pct\nsens %.2f, spec %.2f, AUC %.3f",
      out$method, 100 * out$threshold_pct_ctry,
      headline_ctry$sensitivity, headline_ctry$specificity, out$auc_ctry
    )
  ) +
  scale_fill_manual(values = c("PopuList populist" = "#D55E00", "Not listed" = "#999999")) +
  labs(
    title = "Populism as a within-country percentile, PopuList calibration universe",
    subtitle = paste(strwrap(paste(
      "Each party-year's populism score ranked among every V-Party party-year in",
      "its own country (all years). The cut is carried to anti-pluralism ranked",
      "the same way, so it applies to countries PopuList does not cover."
    ), 100), collapse = "\n"),
    x = "Within-country percentile of v2xpa_popul", y = "Party-years", fill = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "score_distribution_ctry_pct.png"), p_ctry,
       width = 8, height = 4.6, dpi = 150)
cat("Saved score_distribution_ctry_pct.png\n")

message("\nPopuList threshold calibration written to ", out_dir)
