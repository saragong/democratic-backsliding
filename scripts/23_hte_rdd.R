# ==============================================================================
# Heterogeneous effects of a narrow anti-pluralist victory
# (Calonico, Cattaneo, Farrell, Palomba & Titiunik 2025, "Treatment Effect
# Heterogeneity in Regression Discontinuity Designs"; the rdhte package).
#
# The estimator is their eq. 2.5: the local linear RD with a pre-election
# covariate W entering as a level, interacted with treatment T, with the
# running variable X, and with X*T:
#
#   Y = a + T*tau + X*b- + X*T*b+ + W*l + W*T*d + W*X*g- + W*X*T*g+
#
# so the effect at the cutoff is tau + d*W. d is the heterogeneity: how the
# effect of a narrow win moves with W.
#
# DESIGN. Sharp, reduced form: T = the more anti-pluralist party wins
# (running_var > 0), which is the "narrow illiberal victory" treatment, and the
# outcome is PWT GDP-per-capita growth since the pre-election year. rdhte has
# no fuzzy arm, and the first stage is weak enough that the reduced form is the
# object worth decomposing anyway.
#
# W is every pre-election covariate 11b_build_covariates.R writes, one at a
# time:
#   party   the more- and less-illiberal parties' V-Party scores (anti-
#           pluralism, populism, left-right, anti-elitism, GAL-TAN, minority
#           rights) and their gaps
#   country levels at t-1 and 5-year pre-trends of every outcome series
#   design  OECD, decade, presidential, prior backsliding
# Continuous W is centred at its mean among the CLOSE elections -- weighted by
# the triangular kernel at the pooled fit's bandwidth -- and scaled by its SD
# among those same elections (both written to the CSV as w_mean_close /
# w_sd_close). So tau is the effect at the W of the typical close election, which
# is what the pooled RD estimates, and lines up with it; d is per standard
# deviation. (Centred at the full-sample mean instead, tau would be the effect
# at a W the close elections do not have on average, and would not match the
# pooled effect whenever d is non-zero.) Binary W stays
# 0/1: d is the change in the effect from W = 0 to W = 1. rdhte treats any 0/1
# covariate as two groups, so it is fitted that way with bw.joint = TRUE -- one
# bandwidth for both groups, which makes the group difference exactly eq. 2.5's
# d for a dummy (with separate bandwidths it would not be). Decade is a
# factor too: one effect per decade (dated t-1, like every W), plus a Wald test
# that they are all equal. A decade with fewer than HTE_MIN_GROUP elections, or
# fewer than HTE_MIN_GROUP_H inside its bandwidth, is left out of that fit and
# named in the CSV's note. A small joint
# model puts the party-score gaps and the growth pre-trend in together.
#
# ESTIMATOR. rdhte's bias correction differs from rdrobust's, so its pooled
# effect is NOT the headline number from 12_rdd_analysis.R. rdhte fixes the
# bias-correction bandwidth b equal to the main bandwidth h and uses an HC3
# variance; rdrobust (the headline) estimates b separately -- about twice h here
# -- with a nearest-neighbour variance. The conventional estimates agree (full
# sample, w5: -0.047 vs -0.046); the bias-corrected ones do not (-0.078 vs
# -0.053), and rdhte's SE is larger (0.035 vs 0.027). Both are valid robust
# bias-corrected estimates. Every row here uses rdhte, so rows are comparable
# with rdhte's pooled effect, not with the headline; both pooled numbers are
# written to the CSV and the figure. (rdrobust at b = h with vce = "hc3"
# reproduces rdhte's pooled estimate exactly.)
#
# The sample is any spec 12_rdd_analysis.R can run -- prepare_rdd_sample() cuts
# it identically -- and the output lands in that spec's folder.
#
# Output: output/runs/<spec>/hte/
#           hte_w<NN>_<outcome>.csv        one row per W (and per group / slope)
#           hte_w<NN>_<outcome>_forest_<family>.png  one figure per family:
#             party_scores, country_levels, country_pretrends, design,
#             decade, joint
# ==============================================================================

library(tidyverse)
library(here)
library(rdhte)

source(here::here("scripts", "rdd_helpers.R"))

if (!exists("HTE_WINDOW")) HTE_WINDOW <- DEFAULT_WINDOW
if (!exists("HTE_OUTCOME")) HTE_OUTCOME <- "Y_gdp_growth"
# The samples to run, each a list of sample_opts() arguments.
if (!exists("HTE_SAMPLES")) {
  HTE_SAMPLES <- list(
    full = list(),
    popucut = list(illiberal_cutoff = "popucut", other_cutoff_max = "popucut")
  )
}
# Add the outcome's local projection (11b) as an efficiency covariate, per
# Calonico et al. (2019). Off by default so the heterogeneity estimates are
# the plain eq. 2.5 ones.
if (!exists("HTE_COVS_EFF")) HTE_COVS_EFF <- FALSE
# Cluster standard errors by country (rdhte's CR1). Off by default, matching
# the rdrobust runs; countries contribute up to 40 elections, so worth a check.
if (!exists("HTE_CLUSTER")) HTE_CLUSTER <- FALSE
# A W with fewer usable elections than this, or a group smaller than this
# within the bandwidth, is reported but not estimated.
if (!exists("HTE_MIN_N")) HTE_MIN_N <- 100
# A continuous W with more than this share of its elections at ONE value is
# reported but not estimated. It carries almost no variation to identify a
# slope from, and standardising it blows the few other values up into huge
# z-scores: suffrage (98% at one value) gave a CI of about -3 to +5 that
# flattened every other row of the forest plot, and the V-Dem executive-power
# pre-trends (93-94%) returned no standard errors at all. 0.9 sits in a clean
# gap: the next most concentrated W is at 0.85.
if (!exists("HTE_MAX_MODE_SHARE")) HTE_MAX_MODE_SHARE <- 0.9

# Categorical W estimated as one effect per level rather than as a slope.
FACTOR_W <- c("W_decade")
# A level of a categorical W with fewer elections than this is dropped from
# its fit: the 1960s decade, for instance, is only the 17 elections held in
# 1970 (t-1 = 1969), 7 of them within 20 points of the cutoff.
if (!exists("HTE_MIN_GROUP")) HTE_MIN_GROUP <- 30
# ... and a level with fewer elections than this INSIDE its own bandwidth is
# dropped and the fit redone once. Each level gets its own bandwidth, so a
# level can clear HTE_MIN_GROUP overall and still be nearly empty near the
# cutoff (the 1970s in the PopuList sample: 38 elections, 4 in-bandwidth), and
# one degenerate level leaves the whole fit without standard errors.
if (!exists("HTE_MIN_GROUP_H")) HTE_MIN_GROUP_H <- 20

W_FAMILY <- function(w) {
  case_when(
    str_detect(w, "^W_(ill|oth|gap)_") ~ "Party scores",
    str_detect(w, "^W_lev_") ~ "Country history: levels at t-1",
    str_detect(w, "^W_pre5_") ~ "Country history: 5-year pre-trends",
    TRUE ~ "Design"
  )
}

# The party scores in words, with their direction wherever it is not "higher =
# more of the thing named". "More-illiberal" is the top-2 party with the higher
# anti-pluralism score, whose narrow win is the treatment; the gap is always
# that party's score MINUS the less-illiberal party's, so it is signed.
SCORE_WORDS <- c(
  v2xpa_antiplural = "anti-pluralism",
  v2xpa_popul = "populism",
  v2pariglef_neg = "economic-left position (higher = further left)",
  v2paanteli = "anti-elitism",
  ep_galtan = "GAL-TAN (higher = more TAN)",
  v2paminor = "minority-rights support (higher = more supportive)"
)

w_label <- function(w) {
  out <- w
  m <- str_match(w, "^W_(ill|oth|gap)_(.+)$")
  i <- !is.na(m[, 1])
  words <- unname(SCORE_WORDS[m[i, 3]])
  first_word <- sub(" \\(.*$", "", words)
  out[i] <- case_when(
    m[i, 2] == "ill" ~ paste0("More-illiberal party's ", words),
    m[i, 2] == "oth" ~ paste0("Less-illiberal party's ", words),
    TRUE ~ paste0(
      str_to_sentence(first_word), " gap: more- minus less-illiberal party",
      if_else(grepl("(", words, fixed = TRUE), sub("^[^(]*", " ", words), "")
    )
  )
  m <- str_match(w, "^W_(lev|pre5)_(Y_.+)$")
  i <- !is.na(m[, 1])
  out[i] <- paste0(
    c(lev = "Level at t-1: ", pre5 = "5-yr pre-trend: ")[m[i, 2]],
    vapply(m[i, 3], outcome_full_label, character(1))
  )
  out <- recode(
    out,
    W_oecd = "OECD member", W_presidential = "Presidential election",
    W_prior_backsliding = "ERT episode starting in t-5..t-1", W_decade = "Decade"
  )
  out
}

is_binary <- function(x) all(x[!is.na(x)] %in% c(0, 1)) && length(unique(x[!is.na(x)])) == 2

# The mean and SD of W among close elections: weighted by the triangular
# kernel at the pooled bandwidth, i.e. over the elections the RD effect is
# about. W is centred and scaled by these, so a row reads "at the typical close
# election's W (mean), one SD more (sd) changes the effect by d".
close_stats <- function(v, x, h) {
  k <- pmax(0, 1 - abs(x) / h)
  m <- stats::weighted.mean(v, k)
  c(mean = m, sd = sqrt(sum(k * (v - m)^2) / sum(k)))
}
center_at_cutoff <- function(v, x, h) {
  st <- close_stats(v, x, h)
  (v - st[["mean"]]) / st[["sd"]]
}

# The same fit at HALF the data-driven bandwidth, as a standing check that a
# heterogeneity term comes from close elections: at h ~ 20 vote-margin points
# the fits still use races decided by 15-20 points, and a pattern that lives
# there (the polyarchy pre-trend: +0.10 at h = 20.7, ~0 at h = 10) is not
# heterogeneity in the close-election effect. Returned as *_half columns keyed
# by term.
#
# HC1 (CR1 with clustering) rather than rdhte's default HC3: rdhte 0.2.0's HC3
# path mis-sizes its leverage correction when h is supplied (it warns "longer
# object length is not a multiple of shorter object length"). HC1 gives the
# same estimates and SEs within a few percent where both run.
half_bw_rows <- function(y, x, covs, ce, cl, binary, h_main) {
  h_half <- h_main / 2
  f <- tryCatch(
    rdhte(y, x, covs.hte = covs, covs.eff = ce, cluster = cl, h = h_half,
          vce = if (is.null(cl)) "hc1" else "cr1", bw.joint = binary),
    error = function(e) NULL
  )
  if (is.null(f)) {
    return(NULL)
  }
  out <- tibble(
    term = if (binary) paste0("effect at W = ", f$W.lev) else c("effect at cutoff-mean W", "slope per SD of W"),
    estimate_bc_half = unname(f$Estimate.bc), se_rb_half = unname(f$se.rb),
    ci_lo_half = f$ci.rb[, 1], ci_hi_half = f$ci.rb[, 2], pval_half = unname(f$pv.rb),
    h_half = h_half, n_h_half = rowSums(f$Nh)
  )
  if (binary && nrow(out) == 2) {
    v <- f$vcov
    dlt <- out$estimate_bc_half[2] - out$estimate_bc_half[1]
    se <- sqrt(v[1, 1] + v[2, 2] - 2 * v[1, 2])
    z <- qnorm(0.975)
    out <- bind_rows(out, tibble(
      term = "change from W = 0 to 1", estimate_bc_half = dlt, se_rb_half = se,
      ci_lo_half = dlt - z * se, ci_hi_half = dlt + z * se,
      pval_half = 2 * pnorm(-abs(dlt / se)), h_half = h_half, n_h_half = sum(f$Nh)
    ))
  }
  out
}

fit_one_w <- function(dd, w, covs_eff, cluster, h_pooled) {
  wv <- dd[[w]]
  ok <- !is.na(wv) & !is.na(dd[[HTE_OUTCOME]]) &
    (if (is.null(covs_eff)) TRUE else !is.na(covs_eff))
  if (sum(ok) < HTE_MIN_N) {
    return(tibble(w = w, term = NA_character_, note = sprintf("only %d usable", sum(ok))))
  }
  y <- dd[[HTE_OUTCOME]][ok]
  x <- dd$running_var[ok]
  ce <- if (is.null(covs_eff)) NULL else covs_eff[ok]
  cl <- if (cluster) dd$country_text_id[ok] else NULL
  is_factor <- w %in% FACTOR_W
  dropped_note <- NA_character_
  if (is_factor) {
    counts <- table(wv[ok])
    small <- names(counts)[counts < HTE_MIN_GROUP]
    if (length(small) > 0) {
      dropped_note <- sprintf(
        "dropped levels with < %d elections: %s",
        HTE_MIN_GROUP, paste(sprintf("%s (%d)", small, counts[small]), collapse = ", ")
      )
      keep <- !(as.character(wv[ok]) %in% small)
      ok[which(ok)[!keep]] <- FALSE
      y <- dd[[HTE_OUTCOME]][ok]
      x <- dd$running_var[ok]
      ce <- if (is.null(covs_eff)) NULL else covs_eff[ok]
      cl <- if (cluster) dd$country_text_id[ok] else NULL
    }
  }
  binary <- is_binary(wv[ok])
  if (!binary && !is_factor) {
    mode_share <- max(table(wv[ok])) / sum(ok)
    if (mode_share > HTE_MAX_MODE_SHARE) {
      return(tibble(w = w, term = NA_character_, note = sprintf(
        "not estimated: %.0f%% of elections at one value (limit %.0f%%)",
        100 * mode_share, 100 * HTE_MAX_MODE_SHARE
      )))
    }
  }
  st <- if (binary || is_factor) c(mean = NA_real_, sd = NA_real_) else close_stats(wv[ok], x, h_pooled)
  covs <- if (binary || is_factor) factor(wv[ok]) else (wv[ok] - st[["mean"]]) / st[["sd"]]
  fit <- tryCatch(
    rdhte(y, x, covs.hte = covs, covs.eff = ce, cluster = cl, bw.joint = binary),
    error = function(e) e
  )
  if (is_factor && !inherits(fit, "error")) {
    thin <- fit$W.lev[rowSums(fit$Nh) < HTE_MIN_GROUP_H]
    if (length(thin) > 0) {
      dropped_note <- paste(na.omit(c(
        dropped_note,
        sprintf(
          "dropped levels with < %d elections in their bandwidth: %s",
          HTE_MIN_GROUP_H,
          paste(sprintf("%s (%d)", thin, rowSums(fit$Nh)[fit$W.lev %in% thin]), collapse = ", ")
        )
      )), collapse = "; ")
      keep <- !(as.character(covs) %in% thin)
      y <- y[keep]; x <- x[keep]
      ce <- if (is.null(ce)) NULL else ce[keep]
      cl <- if (is.null(cl)) NULL else cl[keep]
      covs <- droplevels(covs[keep])
      ok[which(ok)[!keep]] <- FALSE
      fit <- tryCatch(
        rdhte(y, x, covs.hte = covs, covs.eff = ce, cluster = cl),
        error = function(e) e
      )
    }
  }
  if (inherits(fit, "error")) {
    return(tibble(w = w, term = NA_character_, note = conditionMessage(fit)))
  }
  base <- tibble(
    w = w, binary = binary, n = sum(ok),
    term = if (is_factor) {
      paste0("group ", fit$W.lev)
    } else if (binary) {
      paste0("effect at W = ", fit$W.lev)
    } else {
      c("effect at cutoff-mean W", "slope per SD of W")
    },
    estimate = unname(fit$Estimate), estimate_bc = unname(fit$Estimate.bc),
    se_rb = unname(fit$se.rb), ci_lo = fit$ci.rb[, 1], ci_hi = fit$ci.rb[, 2],
    pval = unname(fit$pv.rb),
    h = fit$h[, 1], n_h = rowSums(fit$Nh),
    # For continuous W: its mean and SD among close elections, in W's own
    # units -- the point tau is evaluated at, and the size of "one SD" in d.
    w_mean_close = st[["mean"]], w_sd_close = st[["sd"]],
    note = dropped_note
  )
  if (binary && nrow(base) == 2) {
    # The heterogeneity for a binary W: the change in the effect from W = 0 to
    # W = 1, from rdhte's own group-level covariance matrix. With one joint
    # bandwidth this is eq. 2.5's d for a 0/1 W.
    v <- fit$vcov
    diff <- base$estimate_bc[2] - base$estimate_bc[1]
    se <- sqrt(v[1, 1] + v[2, 2] - 2 * v[1, 2])
    z <- qnorm(0.975)
    base <- bind_rows(base, tibble(
      w = w, binary = TRUE, n = sum(ok), term = "change from W = 0 to 1",
      estimate = base$estimate[2] - base$estimate[1], estimate_bc = diff,
      se_rb = se, ci_lo = diff - z * se, ci_hi = diff + z * se,
      pval = 2 * pnorm(-abs(diff / se)), note = NA_character_
    ))
  }
  # Half-bandwidth check. Not for decade: its levels each have their own
  # bandwidth, and at half width most fall below HTE_MIN_GROUP_H.
  if (!is_factor) {
    half <- half_bw_rows(y, x, covs, ce, cl, binary, fit$h[1, 1])
    if (!is.null(half)) base <- left_join(base, half, by = "term")
  }
  if (is_factor && nrow(base) > 2) {
    # Heterogeneity across more than two groups: a Wald test that every group
    # effect is equal, contrasting each with the first, on rdhte's own
    # group-level covariance matrix of the bias-corrected estimates.
    b <- base$estimate_bc
    k <- length(b)
    C <- cbind(-1, diag(k - 1))
    cb <- C %*% b
    stat <- drop(t(cb) %*% solve(C %*% fit$vcov %*% t(C)) %*% cb)
    base <- bind_rows(base, tibble(
      w = w, binary = FALSE, n = sum(ok),
      term = sprintf("equality across %d groups (chi2, %d df)", k, k - 1),
      estimate = stat, pval = pchisq(stat, df = k - 1, lower.tail = FALSE),
      note = dropped_note
    ))
  }
  base
}

for (samp_name in names(HTE_SAMPLES)) {
  opts <- do.call(sample_opts, c(list(window = HTE_WINDOW), HTE_SAMPLES[[samp_name]]))
  cat(sprintf("\n==== sample '%s' ====\n", samp_name))
  prep <- prepare_rdd_sample(opts)
  # Only what this script uses: the W covariates, and the outcome's Z when it
  # is the efficiency covariate. Joining all of cv would collide with the Z
  # columns prepare_rdd_sample() has already joined under rd_covariates = "lp".
  cv <- load_covars(opts$instrument, opts$window, opts$incl, opts$placebo)
  cv <- cv[, c("election_id", grep("^W_", names(cv), value = TRUE),
               if (HTE_COVS_EFF) paste0("Z_", HTE_OUTCOME))]
  missing_ids <- setdiff(prep$d$election_id, cv$election_id)
  if (length(missing_ids) > 0) {
    stop(length(missing_ids), " sample elections have no covariates; re-run 11b_build_covariates.R.",
         call. = FALSE)
  }
  dd <- left_join(prep$d, cv, by = "election_id", relationship = "one-to-one")
  covs_eff <- if (HTE_COVS_EFF) dd[[paste0("Z_", HTE_OUTCOME)]] else NULL
  if (HTE_COVS_EFF) stopifnot(!is.null(covs_eff))

  w_vars <- grep("^W_", names(cv), value = TRUE)
  # Covariates that are constant in this sample (e.g. a restriction fixes
  # them) carry no heterogeneity to estimate.
  w_vars <- w_vars[vapply(w_vars, function(w) length(unique(na.omit(dd[[w]]))) > 1, logical(1))]

  # Pooled RD first, as the reference every heterogeneity estimate is read
  # against.
  ok <- !is.na(dd[[HTE_OUTCOME]])
  pooled <- rdhte(
    dd[[HTE_OUTCOME]][ok], dd$running_var[ok],
    covs.eff = if (is.null(covs_eff)) NULL else covs_eff[ok],
    cluster = if (HTE_CLUSTER) dd$country_text_id[ok] else NULL
  )
  # The headline estimator on the same sample, for reference only: see
  # ESTIMATOR in the header for why the two pooled numbers differ.
  headline <- extract_rd(safe_rdrobust(
    dd[[HTE_OUTCOME]], dd$running_var,
    covs = if (is.null(covs_eff)) NULL else covs_eff
  ))
  cat(sprintf(
    "Pooled sharp RD on %s: %.4f (bc %.4f, robust SE %.4f), h = %.2f, N = %d\n  (rdrobust, the headline estimator, on the same sample: bc %.4f, robust SE %.4f -- rdhte's bias-correction bandwidth is h, rdrobust's is estimated)\n",
    HTE_OUTCOME, pooled$Estimate, pooled$Estimate.bc, pooled$se.rb, pooled$h[1, 1], sum(ok),
    headline$coef, headline$se
  ))

  res <- map_dfr(w_vars, function(w) fit_one_w(dd, w, covs_eff, HTE_CLUSTER, pooled$h[1, 1]))

  # A joint model: the three party-score gaps that the "bundled
  # characteristics" question is about, plus the growth pre-trend.
  joint_vars <- intersect(
    c("W_gap_v2xpa_antiplural", "W_gap_v2xpa_popul", "W_gap_v2pariglef_neg", "W_pre5_Y_gdp_growth"),
    w_vars
  )
  z_eff <- if (HTE_COVS_EFF) paste0("Z_", HTE_OUTCOME)
  jd <- as.data.frame(dd[complete.cases(dd[, c(HTE_OUTCOME, joint_vars, z_eff)]), ])
  joint_st <- lapply(jd[joint_vars], close_stats, x = jd$running_var, h = pooled$h[1, 1])
  jd[joint_vars] <- lapply(jd[joint_vars], center_at_cutoff, x = jd$running_var, h = pooled$h[1, 1])
  # rdhte takes several continuous W only as a one-sided formula STRING looked
  # up in `data`, with y, x (and cluster) then given as column names too; a
  # matrix or data frame of W errors inside the package.
  jfit <- do.call(rdhte, c(
    list(
      y = as.name(HTE_OUTCOME), x = as.name("running_var"),
      covs.hte = paste("~", paste(joint_vars, collapse = " + ")),
      data = jd
    ),
    if (HTE_CLUSTER) list(cluster = as.name("country_text_id")),
    # The same efficiency covariate as every other fit, so the covs_eff column
    # is true of the joint rows too.
    if (HTE_COVS_EFF) list(covs.eff = as.name(z_eff))
  ))
  jfit_half <- tryCatch(do.call(rdhte, c(
    list(
      y = as.name(HTE_OUTCOME), x = as.name("running_var"),
      covs.hte = paste("~", paste(joint_vars, collapse = " + ")),
      data = jd, h = jfit$h[1, 1] / 2,
      vce = if (HTE_CLUSTER) "cr1" else "hc1"
    ),
    if (HTE_CLUSTER) list(cluster = as.name("country_text_id")),
    if (HTE_COVS_EFF) list(covs.eff = as.name(z_eff))
  )), error = function(e) NULL)
  joint <- tibble(
    w = c("(joint) intercept", paste0("(joint) ", joint_vars)), binary = FALSE, n = nrow(jd),
    term = c("effect at cutoff-mean W", rep("slope per SD of W", length(joint_vars))),
    estimate = unname(jfit$Estimate), estimate_bc = unname(jfit$Estimate.bc),
    se_rb = unname(jfit$se.rb), ci_lo = jfit$ci.rb[, 1], ci_hi = jfit$ci.rb[, 2],
    pval = unname(jfit$pv.rb), h = jfit$h[1, 1], n_h = sum(jfit$Nh[1, ]),
    w_mean_close = c(NA, vapply(joint_st, `[[`, numeric(1), "mean")),
    w_sd_close = c(NA, vapply(joint_st, `[[`, numeric(1), "sd")),
    note = NA_character_
  )
  if (!is.null(jfit_half)) {
    joint <- joint |> mutate(
      estimate_bc_half = unname(jfit_half$Estimate.bc), se_rb_half = unname(jfit_half$se.rb),
      ci_lo_half = jfit_half$ci.rb[, 1], ci_hi_half = jfit_half$ci.rb[, 2],
      pval_half = unname(jfit_half$pv.rb), h_half = jfit$h[1, 1] / 2,
      n_h_half = sum(jfit_half$Nh[1, ])
    )
  }

  res <- bind_rows(res, joint) |>
    mutate(
      family = if_else(str_starts(w, "\\(joint\\)"), "Joint model", W_FAMILY(w)),
      label = if_else(
        str_starts(w, "\\(joint\\)"),
        paste0("Joint: ", w_label(str_remove(w, "^\\(joint\\) "))),
        w_label(w)
      ),
      sample = samp_name, spec = spec_slug(prep$cfg), window = HTE_WINDOW,
      outcome = HTE_OUTCOME,
      pooled_estimate = pooled$Estimate, pooled_estimate_bc = pooled$Estimate.bc,
      pooled_se_rb = pooled$se.rb,
      pooled_rdrobust_estimate_bc = headline$coef,
      pooled_rdrobust_se_rb = headline$se,
      estimator_note = paste(
        "All rows use rdhte, whose bias-correction bandwidth equals h (vce HC3);",
        "the headline rdrobust estimate (pooled_rdrobust_*) estimates it separately (vce NN),",
        "so its bias-corrected effect differs. Compare rows with pooled_estimate_bc."
      ),
      covs_eff = if (HTE_COVS_EFF) paste0("Z_", HTE_OUTCOME) else "none",
      cluster = if (HTE_CLUSTER) "country" else "none",
      .before = 1
    )

  out_dir <- spec_subdir(prep$cfg, "hte")
  stem <- sprintf(
    "hte_w%02d_%s%s%s", HTE_WINDOW, HTE_OUTCOME,
    if (HTE_COVS_EFF) "_covseff" else "", if (HTE_CLUSTER) "_cluster" else ""
  )
  write_csv(res, file.path(out_dir, paste0(stem, ".csv")))

  # ---- forest plot: the heterogeneity terms only ----------------------------
  plot_dat <- res |>
    filter(
      term %in% c("slope per SD of W", "change from W = 0 to 1") |
        (w %in% FACTOR_W & str_starts(term, "group ")),
      !is.na(estimate_bc)
    ) |>
    mutate(
      # A categorical W is shown as its per-level EFFECTS, in its own panel,
      # since there is no single heterogeneity term to plot.
      family = if_else(w %in% FACTOR_W, "Decade: effect in each", family),
      # The figure title says levels or pre-trends, so the row need not.
      label = if_else(
        str_detect(w, "^W_(lev|pre5)_"),
        str_remove(label, "^(Level at t-1|5-yr pre-trend): "),
        label
      ),
      label = case_when(
        w %in% FACTOR_W ~ paste0(str_remove(term, "^group "), "s"),
        # Name on the first line(s), close-election mean and SD beneath it, so
        # the labels do not squeeze the data panel.
        !is.na(w_sd_close) ~ sprintf(
          "%s\n[close elections: mean %s, SD %s]", str_wrap(label, 60),
          formatC(signif(w_mean_close, 2), format = "fg", digits = 2),
          formatC(signif(w_sd_close, 2), format = "fg", digits = 2)
        ),
        TRUE ~ str_wrap(label, 60)
      ),
      family = factor(family, levels = c(
        "Party scores", "Country history: levels at t-1",
        "Country history: 5-year pre-trends", "Design",
        "Decade: effect in each", "Joint model"
      )),
      sig = pval < 0.05,
      label = fct_reorder(label, estimate_bc)
    )
  # One figure per family: a single sheet of ~100 rows was unreadable. Each is
  # sized to its own rows and shares the subtitle, so any one can stand alone.
  subtitle <- paste(strwrap(sprintf(
    paste(
      "Sample: %s. w%d. Each row is its own rdhte fit (Calonico et al. 2025, eq. 2.5), with the effect at the cutoff linear in W.",
      "Continuous W is centred and scaled by its mean and SD among close elections (kernel-weighted at the pooled bandwidth; shown",
      "in brackets, in W's units), so a row reads: at the typical close election's W, one SD more changes the effect by the estimate;",
      "binary W as 0/1, the change in the effect from W = 0 to W = 1 (one bandwidth). Pooled effect %.3f (robust SE %.3f).",
      "That is rdhte's estimate, which sets the bias-correction bandwidth equal to h; the headline rdrobust estimate on this",
      "sample, which estimates that bandwidth separately, is %.3f (%.3f). Rows are comparable with the rdhte number.",
      "Bias-corrected estimates, robust 95%% CIs; orange = p < 0.05, not adjusted for testing many W. Hollow points: the same fit at half",
      "the data-driven bandwidth (HC1 variance), i.e. closer races only -- a pattern that vanishes there is not about close elections.%s"
    ),
    prep$label, HTE_WINDOW, pooled$Estimate.bc, pooled$se.rb,
    headline$coef, headline$se,
    if (HTE_CLUSTER) " SEs clustered by country." else ""
  ), 170), collapse = "\n")
  FAMILY_FILES <- c(
    "Party scores" = "party_scores",
    "Country history: levels at t-1" = "country_levels",
    "Country history: 5-year pre-trends" = "country_pretrends",
    "Design" = "design",
    "Decade: effect in each" = "decade",
    "Joint model" = "joint"
  )
  # The per-family figures replace the single combined one.
  unlink(file.path(out_dir, paste0(stem, "_forest.png")))
  for (fam in intersect(levels(plot_dat$family), unique(as.character(plot_dat$family)))) {
    fd <- plot_dat |> filter(family == fam) |> mutate(label = fct_drop(label))
    # Two estimates per row: the data-driven bandwidth (filled) and half of
    # it (hollow, just below), so a pattern that disappears in the closest
    # races is visible on the figure itself. Decade has no half-bandwidth fit.
    levs <- levels(fd$label)
    has_half <- "estimate_bc_half" %in% names(fd)
    both <- bind_rows(
      fd |> transmute(y = as.numeric(label), est = estimate_bc, lo = ci_lo, hi = ci_hi,
                      sig = sig, bw = "data-driven h"),
      if (has_half) {
        fd |> transmute(y = as.numeric(label) - 0.3, est = estimate_bc_half, lo = ci_lo_half,
                        hi = ci_hi_half, sig = pval_half < 0.05, bw = "h / 2")
      }
    ) |> filter(!is.na(est))
    p <- ggplot(both, aes(x = est, y = y, colour = sig)) +
      geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
      geom_errorbarh(aes(xmin = lo, xmax = hi, linetype = bw), height = 0, linewidth = 0.7) +
      geom_point(aes(shape = bw), size = 2.2, fill = "white") +
      scale_shape_manual(values = c("data-driven h" = 16, "h / 2" = 21), name = NULL) +
      scale_linetype_manual(values = c("data-driven h" = "solid", "h / 2" = "22"), name = NULL) +
      scale_y_continuous(breaks = seq_along(levs), labels = levs, expand = expansion(add = 0.6)) +
      scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey35"), guide = "none", na.value = "grey35") +
      labs(
        title = sprintf(
          "%s. Heterogeneity in the effect of a narrow anti-pluralist win on %s",
          fam, outcome_full_label(HTE_OUTCOME)
        ),
        subtitle = subtitle,
        # The decade panel shows each decade's own effect, not a change in it.
        x = if (fam == "Decade: effect in each") {
          "RD effect in each decade (log points)"
        } else {
          "Change in the RD effect (log points)"
        },
        y = NULL
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 13),
        plot.subtitle = element_text(size = 9),
        axis.text.y = element_text(size = 10.5),
        # Anchored to the whole figure, not the panel: the y labels are long,
        # and a panel-anchored title runs off the right edge.
        plot.title.position = "plot",
        legend.position = "top", legend.justification = "left"
      )
    ggsave(
      file.path(out_dir, sprintf("%s_forest_%s.png", stem, FAMILY_FILES[[fam]])), p,
      width = 12, height = 3.2 + 0.42 * nrow(fd), dpi = 150, limitsize = FALSE
    )
  }
  cat(sprintf("Saved %s (%d W)\n", file.path(out_dir, paste0(stem, ".csv")), length(w_vars)))
  print(
    plot_dat |>
      arrange(pval) |>
      head(8) |>
      transmute(label = substr(as.character(label), 1, 60), term, estimate_bc = round(estimate_bc, 4),
                se_rb = round(se_rb, 4), pval = round(pval, 3), n),
    width = Inf
  )
}
