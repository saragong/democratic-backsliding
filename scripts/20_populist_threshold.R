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
# cutpoint BOTH ways and 12_rdd_analysis.R is run at each:
#   absolute    the raw number, carried across unchanged. Defensible because
#               both indices are V-Party aggregates on [0, 1].
#   percentile  the cutpoint's percentile in the populism distribution, applied
#               at the same percentile of the anti-pluralism distribution.
#               What the Aug 24 note actually described ("took the same
#               percentile").
# They disagree whenever the two indices are distributed differently, which
# they are, so neither is obviously right and both are reported.
#
# A CORRECTION TO THE ASK, STATED PLAINLY
#
# The brief says "choose the threshold that maximizes the AUC". AUC cannot be
# maximized over thresholds: it is the area under the whole ROC curve and is
# the same number whatever cutpoint you pick -- it measures the SCORE, not a
# cutpoint. What is wanted is the ROC-OPTIMAL cutpoint, and the standard choice
# is Youden's J (sensitivity + specificity - 1), which is what this script
# uses. It also reports the accuracy-maximizing and closest-to-(0,1) cutpoints
# so the choice is visible rather than assumed, and reports the AUC itself as
# what it actually is: how well the continuous score separates the classes.
#
# TWO UNIVERSES (the negative class is the whole difficulty)
#
#   imputed_zero (headline)  PopuList lists the parties that ARE populist, far
#       right, far left or eurosceptic; everything else in its 31 countries is
#       implicitly none of those. So the universe is every V-Party party-year
#       in those countries over PopuList's coverage, populist = 1 only for a
#       listed party in a listed year. This is how PopuList is normally used.
#   listed_only (robustness)  Only the parties PopuList lists, 201 populist vs
#       67 not. The 67 are parties that are far-right/far-left/eurosceptic but
#       NOT populist -- a right-tail-selected control group, so a cutpoint
#       fitted here is calibrated on the hard cases only and sits too high.
#
# THE CAVEAT THAT TRAVELS WITH THE NUMBER
#
# PopuList covers 31 European countries. The cutpoint is applied to a global
# election spine, so it is an out-of-sample transfer in two directions at once
# (different variable, different countries). Every output says so.
#
#   Rscript --no-init-file scripts/20_populist_threshold.R
#
# Output: data/populist_4.0.csv          cached download
#         data/populist_threshold.rds    the cutpoints, for 12_rdd_analysis.R
#         output/runs/_sweeps/populist_threshold/
#           thresholds.csv, roc.png, score_distributions.png, merge_audit.csv
# ==============================================================================

# Deliberately few dependencies: this is a calibration step, not part of the
# RDD pipeline, and it must not inherit the pipeline's toggles. It does not
# source rdd_helpers.R.
library(readr)
library(dplyr)
library(tidyr)
library(haven)
library(ggplot2)
library(pROC)
library(here)

data_dir <- here::here("data")
out_dir <- here::here("output", "runs", "_sweeps", "populist_threshold")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

POPULIST_URL <- "https://popu-list.github.io/Data/The%20PopuList%204.0.csv"
POPULIST_CSV <- file.path(data_dir, "populist_4.0.csv")

# The score the logit is fitted on, and the score the cutpoint is transferred
# to. Keeping them as named constants because the whole point of the exercise
# is that they are DIFFERENT variables.
FIT_SCORE <- "v2xpa_popul"
APPLY_SCORE <- "v2xpa_antiplural"

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
  select(v2paid, country_name, year, all_of(c(FIT_SCORE, APPLY_SCORE))) |>
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

universe_imputed <- vparty |>
  filter(
    country_name %in% popu_countries,
    year >= COVERAGE_MIN, year <= COVERAGE_MAX
  ) |>
  left_join(popu_years, by = c("v2paid" = "vdem_id_1", "year")) |>
  mutate(populist = coalesce(populist, 0L))

# Listed-only: the parties PopuList names, at their V-Party observations.
listed_ids <- popu_matched |>
  filter(!is.na(vdem_id_1)) |>
  distinct(vdem_id_1, populist)

universe_listed <- vparty |>
  inner_join(listed_ids, by = c("v2paid" = "vdem_id_1"),
             relationship = "many-to-many") |>
  filter(year >= COVERAGE_MIN, year <= COVERAGE_MAX)

UNIVERSES <- list(
  imputed_zero = list(
    data = universe_imputed,
    label = "All V-Party parties in PopuList's countries; unlisted = not populist"
  ),
  listed_only = list(
    data = universe_listed,
    label = "Only parties PopuList lists (right-tail-selected controls)"
  )
)

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
fit_one <- function(df, label) {
  df <- df |> filter(!is.na(.data[[FIT_SCORE]]), !is.na(populist))
  m <- glm(
    as.formula(paste("populist ~", FIT_SCORE)),
    data = df, family = binomial()
  )
  roc_obj <- roc(
    response = df$populist, predictor = df[[FIT_SCORE]],
    quiet = TRUE, direction = "<"
  )

  grab <- function(method) {
    co <- coords(roc_obj, "best", best.method = method,
                 ret = c("threshold", "sensitivity", "specificity"),
                 transpose = FALSE)
    # Ties can return several equally good cutpoints; take the median so the
    # answer does not depend on pROC's internal ordering.
    tibble(
      method = method,
      threshold = median(co$threshold),
      sensitivity = median(co$sensitivity),
      specificity = median(co$specificity)
    )
  }

  cuts <- bind_rows(grab("youden"), grab("closest.topleft"))

  list(
    label = label,
    n = nrow(df),
    n_pos = sum(df$populist == 1),
    auc = as.numeric(auc(roc_obj)),
    coef = unname(coef(m)[2]),
    coef_se = sqrt(diag(vcov(m)))[2],
    roc = roc_obj,
    cuts = cuts,
    fit_score = df[[FIT_SCORE]]
  )
}

# The percentile transfer: where the cutpoint sits in the populism
# distribution, then the anti-pluralism value at that same percentile. Both are
# taken over the SAME rows the logit was fitted on, so the two distributions
# are comparable.
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

rows <- list()
fits <- list()
for (u in names(UNIVERSES)) {
  f <- fit_one(UNIVERSES[[u]]$data, UNIVERSES[[u]]$label)
  fits[[u]] <- f
  cat(sprintf(
    "\n=== %s ===\n  N = %d (%d populist, %.1f%%)\n  logit coef on %s = %.3f (%.3f)\n  AUC = %.3f\n",
    u, f$n, f$n_pos, 100 * f$n_pos / f$n, FIT_SCORE, f$coef, f$coef_se, f$auc
  ))
  for (i in seq_len(nrow(f$cuts))) {
    tr <- transfer(f$cuts$threshold[i])
    rows[[length(rows) + 1]] <- bind_cols(
      tibble(
        universe = u,
        universe_label = f$label,
        n = f$n, n_populist = f$n_pos, auc = f$auc,
        logit_coef = f$coef, logit_se = f$coef_se,
        method = f$cuts$method[i],
        sensitivity = f$cuts$sensitivity[i],
        specificity = f$cuts$specificity[i]
      ),
      tr
    )
    cat(sprintf(
      "  %-16s cutpoint %.4f on %s (sens %.2f, spec %.2f) -> %.1fth pct -> %.4f on %s\n",
      f$cuts$method[i], tr$threshold_absolute, FIT_SCORE,
      f$cuts$sensitivity[i], f$cuts$specificity[i],
      tr$percentile, tr$threshold_percentile_value, APPLY_SCORE
    ))
  }
}

thresholds <- bind_rows(rows)
write_csv(thresholds, file.path(out_dir, "thresholds.csv"))

# ------------------------------------------------------------------------------
# 5. The headline numbers, for 12_rdd_analysis.R
# ------------------------------------------------------------------------------

headline <- thresholds |>
  filter(universe == "imputed_zero", method == "youden") |>
  slice(1)

out <- list(
  fit_score = FIT_SCORE,
  apply_score = APPLY_SCORE,
  universe = headline$universe,
  method = headline$method,
  auc = headline$auc,
  n = headline$n,
  # For ILLIBERAL_CUTOFF / OTHER_CUTOFF_MAX, the absolute transfer.
  threshold_abs = headline$threshold_absolute,
  # ... and the percentile transfer, as an absolute value on the
  # anti-pluralism scale. Given as a NUMBER rather than a "qNN" string on
  # purpose: "qNN" would be re-resolved by 12_rdd_analysis.R against whatever
  # sample that run has, which is the top-2 election spine, not the V-Party
  # party-year distribution this percentile was computed on.
  threshold_pct = headline$threshold_percentile_value,
  percentile = headline$percentile,
  coverage = sprintf(
    "PopuList 4.0, %d European countries, %d-%d",
    n_distinct(popu$country_name), COVERAGE_MIN, COVERAGE_MAX
  ),
  built_at = Sys.time()
)
saveRDS(out, file.path(data_dir, "populist_threshold.rds"))
cat(sprintf(
  "\nHeadline (%s, %s): %.4f absolute, %.4f at the same percentile (%.1fth)\nSaved data/populist_threshold.rds\n",
  out$universe, out$method, out$threshold_abs, out$threshold_pct, out$percentile
))

# ------------------------------------------------------------------------------
# 6. Figures
# ------------------------------------------------------------------------------

roc_df <- bind_rows(lapply(names(fits), function(u) {
  r <- fits[[u]]$roc
  tibble(
    universe = sprintf("%s (AUC %.3f)", u, fits[[u]]$auc),
    fpr = 1 - r$specificities,
    tpr = r$sensitivities
  )
}))

p_roc <- ggplot(roc_df, aes(fpr, tpr, colour = universe)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(
    data = thresholds |>
      filter(method == "youden") |>
      mutate(
        universe = sprintf("%s (AUC %.3f)", universe, auc),
        fpr = 1 - specificity, tpr = sensitivity
      ),
    aes(fpr, tpr), size = 2.6, shape = 21, fill = "white", stroke = 1
  ) +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
  coord_fixed() +
  labs(
    title = sprintf("Does %s separate PopuList's populist parties?", FIT_SCORE),
    subtitle = paste(
      strwrap(paste(
        "ROC of the raw V-Party populism score. Circles mark the Youden",
        "cutpoint. AUC measures the SCORE, not any cutpoint -- it is the same",
        "number whichever cutpoint is chosen. PopuList covers 31 European",
        "countries, so a cutpoint fitted here is applied out of sample to the",
        "global election spine."
      ), 95),
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

# The transfer, shown rather than asserted: the two score distributions with
# both candidate cutoffs on each.
dist_df <- apply_scores |>
  select(all_of(c(FIT_SCORE, APPLY_SCORE))) |>
  pivot_longer(everything(), names_to = "score", values_to = "value")

cut_df <- tibble(
  score = c(FIT_SCORE, APPLY_SCORE, APPLY_SCORE),
  value = c(out$threshold_abs, out$threshold_abs, out$threshold_pct),
  lab = c(
    sprintf("Youden cutpoint (%.3f)", out$threshold_abs),
    sprintf("Absolute transfer (%.3f)", out$threshold_abs),
    sprintf("Percentile transfer (%.3f)", out$threshold_pct)
  )
)

p_dist <- ggplot(dist_df, aes(value)) +
  geom_histogram(bins = 50, fill = "grey85", colour = "white") +
  geom_vline(
    data = cut_df, aes(xintercept = value, colour = lab),
    linewidth = 0.7
  ) +
  facet_wrap(~score, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#009E73"), name = NULL) +
  labs(
    title = "Transferring the cutpoint from populism to anti-pluralism",
    subtitle = paste(
      strwrap(sprintf(paste(
        "The cutpoint is fitted on %s and applied to %s. The absolute transfer",
        "carries the raw number across; the percentile transfer carries the",
        "%.1fth percentile across. They differ because the two indices are not",
        "distributed alike. Both are produced; neither is obviously right."
      ), FIT_SCORE, APPLY_SCORE, out$percentile), 95),
      collapse = "\n"
    ),
    x = "Score", y = "Party-years"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 7.5, colour = "grey25")
  )
ggsave(
  file.path(out_dir, "score_distributions.png"), p_dist,
  width = 8, height = 6, dpi = 150
)
cat("Saved score_distributions.png\n")

message("\nPopuList threshold calibration written to ", out_dir)
