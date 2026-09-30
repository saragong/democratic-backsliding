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
# Continuous W is standardised within the estimation sample, so tau is the
# effect at the sample mean and d is per standard deviation. Binary W is a
# factor, giving one effect per group and their difference. Decade is a
# factor too: one effect per decade (dated t-1, like every W), plus a Wald test
# that they are all equal. A decade with fewer than HTE_MIN_GROUP elections, or
# fewer than HTE_MIN_GROUP_H inside its bandwidth, is left out of that fit and
# named in the CSV's note. A small joint
# model puts the party-score gaps and the growth pre-trend in together.
#
# The sample is any spec 12_rdd_analysis.R can run -- prepare_rdd_sample() cuts
# it identically -- and the output lands in that spec's folder.
#
# Output: output/runs/<spec>/hte/
#           hte_w<NN>_<outcome>.csv        one row per W (and per group / slope)
#           hte_w<NN>_<outcome>_forest.png heterogeneity slopes and group gaps
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
    str_detect(w, "^W_(lev|pre5)_") ~ "Country history (pre-election)",
    TRUE ~ "Design"
  )
}

w_label <- function(w) {
  score_lab <- function(s) unname(PARTY_SCORE_LABELS[s] %||% s)
  out <- w
  m <- str_match(w, "^W_(ill|oth|gap)_(.+)$")
  i <- !is.na(m[, 1])
  out[i] <- paste0(
    c(ill = "More-illiberal party: ", oth = "Less-illiberal party: ", gap = "Gap (more - less): ")[m[i, 2]],
    vapply(m[i, 3], score_lab, character(1))
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
    W_prior_backsliding = "ERT episode in the prior window", W_decade = "Decade"
  )
  out
}

is_binary <- function(x) all(x[!is.na(x)] %in% c(0, 1)) && length(unique(x[!is.na(x)])) == 2

fit_one_w <- function(dd, w, covs_eff, cluster) {
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
  covs <- if (binary || is_factor) factor(wv[ok]) else as.numeric(scale(wv[ok]))
  fit <- tryCatch(
    rdhte(y, x, covs.hte = covs, covs.eff = ce, cluster = cl),
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
  grouped <- binary || is_factor
  base <- tibble(
    w = w, binary = binary, n = sum(ok),
    term = if (grouped) paste0("group ", fit$W.lev) else c("effect at mean W", "slope per SD of W"),
    estimate = unname(fit$Estimate), estimate_bc = unname(fit$Estimate.bc),
    se_rb = unname(fit$se.rb), ci_lo = fit$ci.rb[, 1], ci_hi = fit$ci.rb[, 2],
    pval = unname(fit$pv.rb),
    h = fit$h[, 1], n_h = rowSums(fit$Nh),
    note = dropped_note
  )
  if (binary && nrow(base) == 2) {
    # The heterogeneity for a binary W is the difference between the two
    # group effects, from rdhte's own group-level covariance matrix.
    v <- fit$vcov
    diff <- base$estimate_bc[2] - base$estimate_bc[1]
    se <- sqrt(v[1, 1] + v[2, 2] - 2 * v[1, 2])
    z <- qnorm(0.975)
    base <- bind_rows(base, tibble(
      w = w, binary = TRUE, n = sum(ok), term = "difference (1 - 0)",
      estimate = base$estimate[2] - base$estimate[1], estimate_bc = diff,
      se_rb = se, ci_lo = diff - z * se, ci_hi = diff + z * se,
      pval = 2 * pnorm(-abs(diff / se)), note = NA_character_
    ))
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
  cat(sprintf(
    "Pooled sharp RD on %s: %.4f (bc %.4f, robust SE %.4f), h = %.2f, N = %d\n",
    HTE_OUTCOME, pooled$Estimate, pooled$Estimate.bc, pooled$se.rb, pooled$h[1, 1], sum(ok)
  ))

  res <- map_dfr(w_vars, function(w) fit_one_w(dd, w, covs_eff, HTE_CLUSTER))

  # A joint model: the three party-score gaps that the "bundled
  # characteristics" question is about, plus the growth pre-trend.
  joint_vars <- intersect(
    c("W_gap_v2xpa_antiplural", "W_gap_v2xpa_popul", "W_gap_v2pariglef_neg", "W_pre5_Y_gdp_growth"),
    w_vars
  )
  z_eff <- if (HTE_COVS_EFF) paste0("Z_", HTE_OUTCOME)
  jd <- as.data.frame(dd[complete.cases(dd[, c(HTE_OUTCOME, joint_vars, z_eff)]), ])
  jd[joint_vars] <- lapply(jd[joint_vars], function(v) as.numeric(scale(v)))
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
  joint <- tibble(
    w = c("(joint) intercept", paste0("(joint) ", joint_vars)), binary = FALSE, n = nrow(jd),
    term = c("effect at mean W", rep("slope per SD of W", length(joint_vars))),
    estimate = unname(jfit$Estimate), estimate_bc = unname(jfit$Estimate.bc),
    se_rb = unname(jfit$se.rb), ci_lo = jfit$ci.rb[, 1], ci_hi = jfit$ci.rb[, 2],
    pval = unname(jfit$pv.rb), h = jfit$h[1, 1], n_h = sum(jfit$Nh[1, ]),
    note = NA_character_
  )

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
      term %in% c("slope per SD of W", "difference (1 - 0)") |
        (w %in% FACTOR_W & str_starts(term, "group ")),
      !is.na(estimate_bc)
    ) |>
    mutate(
      # A categorical W is shown as its per-level EFFECTS, in its own panel,
      # since there is no single heterogeneity term to plot.
      family = if_else(w %in% FACTOR_W, "Decade: effect in each", family),
      label = if_else(
        w %in% FACTOR_W,
        paste0(str_remove(term, "^group "), "s"),
        label
      ),
      family = factor(family, levels = c(
        "Party scores", "Country history (pre-election)", "Design",
        "Decade: effect in each", "Joint model"
      )),
      sig = pval < 0.05,
      label = fct_reorder(label, estimate_bc)
    )
  p <- ggplot(plot_dat, aes(x = estimate_bc, y = label, colour = sig)) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0, linewidth = 0.5) +
    geom_point(size = 1.4) +
    facet_grid(family ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey35"), guide = "none") +
    labs(
      title = sprintf("Heterogeneity in the effect of a narrow anti-pluralist win on %s", outcome_full_label(HTE_OUTCOME)),
      subtitle = paste(strwrap(sprintf(
        paste(
          "Sample: %s. w%d. Each row is its own rdhte fit (Calonico et al. 2025, eq. 2.5):",
          "continuous W standardised, so the estimate is the change in the RD effect per SD of W;",
          "binary W, the difference between the two group effects. Pooled effect %.3f (robust SE %.3f).",
          "Bias-corrected estimates, robust 95%% CIs; orange = p < 0.05.%s"
        ),
        prep$label, HTE_WINDOW, pooled$Estimate.bc, pooled$se.rb,
        if (HTE_CLUSTER) " SEs clustered by country." else ""
      ), 130), collapse = "\n"),
      x = "Heterogeneity in the RD effect (log points)", y = NULL
    ) +
    theme_bw(base_size = 8) +
    theme(
      strip.text.y = element_text(angle = 0, size = 7),
      panel.grid.minor = element_blank(),
      plot.subtitle = element_text(size = 7)
    )
  ggsave(
    file.path(out_dir, paste0(stem, "_forest.png")), p,
    width = 10, height = 2 + 0.13 * nrow(plot_dat), dpi = 150, limitsize = FALSE
  )
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
