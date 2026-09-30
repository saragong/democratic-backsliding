# ==============================================================================
# Pre-election covariates for the RD, one file per build.
#
# Two kinds, for two different uses:
#
#   Z_<outcome>  covariate ADJUSTMENT (Calonico, Cattaneo, Farrell & Titiunik
#                2019, "Regression Discontinuity Designs Using Covariates").
#                For each outcome and window N, a local projection of that
#                outcome's N-year change on pre-election information only,
#                fitted on the country-year panel. It enters the RD linearly
#                and is never interacted with treatment -- the point is
#                precision, not heterogeneity.
#
#   W_*          heterogeneity covariates (Calonico, Cattaneo, Farrell,
#                Palomba & Titiunik 2025). The more- and less-illiberal
#                parties' V-Party scores and their gaps, and the country's
#                pre-election levels and 5-year trends of every outcome
#                series, plus OECD, decade and election type.
#
# THE LOCAL PROJECTION. For outcome Y built from panel series v (see
# OUTCOME_SOURCES in rdd_helpers.R) and window N, over every country-year t:
#
#   LHS_t = v_{t+N} - v_{t-1}            the outcome's own definition,
#           (or sum_{k=0..N} log(1 + r_{t+k}/100) for a compounded rate)
#
#   on:     v_{t-1}; its annual changes at t-1, t-2, t-3; its N-year
#           pre-trend v_{t-1} - v_{t-1-N}; PWT log GDP per capita at t-1 and
#           its growth at t-1, t-2; polyarchy at t-1; decade dummies.
#
# Every regressor is dated t-1 or earlier, so for an election in year t the
# projection uses nothing the election could have caused. A missing regressor
# is set to 0 with its own missing-value dummy, so one gappy series does not
# make Z missing for the whole row.
#
# LEAVE-COUNTRY-OUT. Each election's Z comes from coefficients estimated on
# every OTHER country's country-years. Fitted on the full panel, the
# projection would have partly seen each election's own realised outcome (its
# country-year is in the fitting sample, LHS and all), and Z would carry some
# of Y's post-election variation into the RD as a "pre-election" covariate.
# Computed by downdating one X'X and X'y, not by refitting.
#
# Data:   data/combined_panel.rds, data/rdd_build/rdd_<instr>_w<N><sfx>.rds
# Output: data/rdd_build/covars_<instr>_w<N><sfx>.rds   (election_id, Z_*, W_*)
#         output/builds/<instr>_w<N><sfx>/lp_fit.csv    out-of-sample fit of each Z
# ==============================================================================

library(tidyverse)
library(here)

source(here::here("scripts", "rdd_helpers.R"))
source(here::here("scripts", "vparty_helpers.R"))

if (!exists("ILLIBERALISM_VAR")) ILLIBERALISM_VAR <- DEFAULT_INSTRUMENT
if (!exists("COVAR_WINDOWS")) COVAR_WINDOWS <- DEFAULT_WINDOWS
if (!exists("TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR")) {
  TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR <- DEFAULT_INCL_ELECTION_YEAR
}
# The projection targets the POST-election window. A placebo build's outcomes
# are pre-election changes, which the regressors here would overlap.
PLACEBO_PRE_WINDOW <- FALSE

# Years before this contribute nothing to the fit: the series used here are
# essentially absent, and the decade dummies would be fitted on a handful of
# rows.
LP_YEAR_MIN <- 1950

# ------------------------------------------------------------------------------
# The panel on a complete country x year grid, so a lag of k rows is a lag of
# k years.
# ------------------------------------------------------------------------------

source_vars <- unique(c(unname(OUTCOME_SOURCES), "ln_gdp_pc", "v2x_polyarchy"))
panel <- readRDS(here::here("data", "combined_panel.rds")) |>
  select(country_text_id, year, all_of(source_vars)) |>
  complete(country_text_id, year = full_seq(year, 1)) |>
  arrange(country_text_id, year)

by_country <- function(x, f) ave(x, panel$country_text_id, FUN = f)
lag_k <- function(x, k) by_country(x, function(z) dplyr::lag(z, k))
lead_k <- function(x, k) by_country(x, function(z) dplyr::lead(z, k))

# The series a change is taken in: the level for differenced outcomes, the
# per-year log growth factor for compounded rates.
series_for <- function(outcome) {
  v <- panel[[OUTCOME_SOURCES[[outcome]]]]
  if (outcome %in% OUTCOME_COMPOUNDED) log(1 + v / 100) else v
}

# LHS for window N, exactly as window_change() defines the outcome.
lp_lhs <- function(outcome, N) {
  s <- series_for(outcome)
  if (outcome %in% OUTCOME_COMPOUNDED) {
    Reduce(`+`, lapply(0:N, function(k) lead_k(s, k)))
  } else {
    lead_k(s, N) - lag_k(s, 1)
  }
}

# Own-history regressors for window N, all dated t-1 or earlier.
own_history <- function(outcome, N) {
  s <- series_for(outcome)
  if (outcome %in% OUTCOME_COMPOUNDED) {
    list(
      own_lev = lag_k(s, 1),
      own_d1 = lag_k(s, 2),
      own_d2 = lag_k(s, 3),
      own_pre = Reduce(`+`, lapply(1:N, function(k) lag_k(s, k)))
    )
  } else {
    list(
      own_lev = lag_k(s, 1),
      own_d1 = lag_k(s, 1) - lag_k(s, 2),
      own_d2 = lag_k(s, 2) - lag_k(s, 3),
      own_d3 = lag_k(s, 3) - lag_k(s, 4),
      own_pre = lag_k(s, 1) - lag_k(s, 1 + N)
    )
  }
}

common <- list(
  gdp_lev = lag_k(panel$ln_gdp_pc, 1),
  gdp_d1 = lag_k(panel$ln_gdp_pc, 1) - lag_k(panel$ln_gdp_pc, 2),
  gdp_d2 = lag_k(panel$ln_gdp_pc, 2) - lag_k(panel$ln_gdp_pc, 3),
  poly_lev = lag_k(panel$v2x_polyarchy, 1)
)
decade <- factor((panel$year %/% 10) * 10)

# Design matrix: regressors zero-filled with a missing dummy each, an
# intercept, decade dummies. The GDP and polyarchy controls are dropped where
# they would duplicate the outcome's own history.
design <- function(outcome, N) {
  regs <- own_history(outcome, N)
  src <- OUTCOME_SOURCES[[outcome]]
  keep_common <- common
  if (src == "ln_gdp_pc") keep_common <- keep_common[c("poly_lev")]
  if (src == "v2x_polyarchy") keep_common$poly_lev <- NULL
  regs <- c(regs, keep_common)
  cols <- list()
  for (nm in names(regs)) {
    x <- regs[[nm]]
    miss <- is.na(x)
    x[miss] <- 0
    cols[[nm]] <- x
    if (any(miss)) cols[[paste0(nm, "_miss")]] <- as.numeric(miss)
  }
  X <- cbind(1, do.call(cbind, cols), stats::model.matrix(~ decade)[, -1, drop = FALSE])
  colnames(X)[1] <- "(Intercept)"
  X
}

# ------------------------------------------------------------------------------
# Leave-country-out projection
# ------------------------------------------------------------------------------

# Coefficients for every country at once, by downdating the full-sample normal
# equations: beta_{-c} = (X'X - X_c'X_c)^{-1} (X'y - X_c'y_c). Columns that are
# aliased in the full fitting sample (a missing dummy that never varies there)
# are dropped first, and a country whose removal still leaves the system
# singular falls back to a pseudo-inverse.
lco_predict <- function(X, y, fit_rows, pred_rows, country) {
  fit <- fit_rows & !is.na(y)
  qx <- qr(X[fit, , drop = FALSE])
  keep_cols <- qx$pivot[seq_len(qx$rank)]
  X <- X[, keep_cols, drop = FALSE]
  Xf <- X[fit, , drop = FALSE]
  yf <- y[fit]
  cf <- country[fit]
  XtX <- crossprod(Xf)
  Xty <- crossprod(Xf, yf)
  out <- rep(NA_real_, nrow(X))
  for (ct in unique(country[pred_rows])) {
    ic <- cf == ct
    A <- XtX - crossprod(Xf[ic, , drop = FALSE])
    b <- Xty - crossprod(Xf[ic, , drop = FALSE], yf[ic])
    beta <- tryCatch(solve(A, b), error = function(e) MASS::ginv(A) %*% b)
    ip <- which(pred_rows & country == ct)
    out[ip] <- drop(X[ip, , drop = FALSE] %*% beta)
  }
  out
}

# ------------------------------------------------------------------------------
# Heterogeneity covariates (window-independent country history, plus the
# party scores from the build)
# ------------------------------------------------------------------------------

w_country <- list()
for (oc in names(OUTCOME_SOURCES)) {
  s <- series_for(oc)
  w_country[[paste0("W_lev_", oc)]] <- lag_k(s, 1)
  w_country[[paste0("W_pre5_", oc)]] <- if (oc %in% OUTCOME_COMPOUNDED) {
    Reduce(`+`, lapply(1:5, function(k) lag_k(s, k)))
  } else {
    lag_k(s, 1) - lag_k(s, 6)
  }
}
w_country <- bind_cols(
  panel |> select(country_text_id, year),
  as_tibble(w_country)
)
# Several outcomes share a source series (the V-Dem indices each have their
# own, but the three GDP series and the executive-power pairs do not), so drop
# exact duplicate columns.
w_country <- w_country[, !duplicated(as.list(w_country))]

# ------------------------------------------------------------------------------
# Per window
# ------------------------------------------------------------------------------

fit_rows <- panel$year >= LP_YEAR_MIN

for (N in COVAR_WINDOWS) {
  d <- load_build(ILLIBERALISM_VAR, N, TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR, FALSE)
  cat(sprintf("[w=%d] %d elections\n", N, nrow(d)))
  row_of <- match(
    paste(d$country_text_id, d$election_year),
    paste(panel$country_text_id, panel$year)
  )
  if (anyNA(row_of)) {
    stop(sum(is.na(row_of)), " elections have no country-year in the panel.", call. = FALSE)
  }
  pred_rows <- seq_len(nrow(panel)) %in% row_of

  Z <- list()
  fit_stats <- list()
  for (oc in names(OUTCOME_SOURCES)) {
    y <- lp_lhs(oc, N)
    # The registry has to reproduce the build's own outcome at every election,
    # or Z is projecting a different variable from the one the RD estimates.
    y_e <- y[row_of]
    ok <- isTRUE(all.equal(y_e, d[[oc]], tolerance = 1e-8, check.attributes = FALSE))
    if (!ok) {
      stop(
        "[w=", N, "] ", oc, " recomputed from OUTCOME_SOURCES does not match the build ",
        "(", sum(xor(is.na(y_e), is.na(d[[oc]]))), " NA mismatches). Fix the registry ",
        "in rdd_helpers.R to agree with build_outcomes() in 11_build_rdd_data.R.",
        call. = FALSE
      )
    }
    X <- design(oc, N)
    z <- lco_predict(X, y, fit_rows, pred_rows, panel$country_text_id)[row_of]
    Z[[paste0("Z_", oc)]] <- z
    yy <- d[[oc]]
    has <- !is.na(yy) & !is.na(z)
    fit_stats[[oc]] <- tibble(
      window = N, outcome = oc, n_fit = sum(fit_rows & !is.na(y)),
      n_elections = sum(has),
      oos_r2 = 1 - sum((yy[has] - z[has])^2) / sum((yy[has] - mean(yy[has]))^2),
      cor = cor(yy[has], z[has])
    )
  }
  fit_stats <- bind_rows(fit_stats)

  sides <- illiberal_side_scores(d)
  w_party <- list()
  for (sc in PARTY_SCORE_VARS) {
    if (!paste0("ill_", sc) %in% names(sides)) next
    w_party[[paste0("W_ill_", sc)]] <- sides[[paste0("ill_", sc)]]
    w_party[[paste0("W_oth_", sc)]] <- sides[[paste0("oth_", sc)]]
    w_party[[paste0("W_gap_", sc)]] <- sides[[paste0("ill_", sc)]] - sides[[paste0("oth_", sc)]]
  }

  covars <- bind_cols(
    tibble(election_id = d$election_id),
    as_tibble(Z),
    as_tibble(w_party),
    tibble(
      W_prior_backsliding = as.numeric(d$prior_backsliding),
      W_oecd = as.numeric(oecd_group(d$country_text_id) == "OECD"),
      W_decade = (d$election_year %/% 10) * 10,
      W_presidential = as.numeric(d$election_type == "presidential")
    ),
    w_country[row_of, setdiff(names(w_country), c("country_text_id", "year"))]
  )
  attr(covars, "lp_fit") <- fit_stats
  path <- covars_path(ILLIBERALISM_VAR, N, TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR, FALSE)
  saveRDS(covars, path)

  diag_dir <- file.path(
    BUILDS_OUT_ROOT,
    sprintf("%s_w%d%s", ILLIBERALISM_VAR, N, build_suffix(TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR, FALSE))
  )
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(fit_stats, file.path(diag_dir, "lp_fit.csv"))
  cat(sprintf(
    "[w=%d] saved %s (%d Z, %d W). Out-of-sample R^2 of Z: GDP (PWT) %.3f, inflation %.3f, polyarchy %.3f; median over outcomes %.3f\n",
    N, basename(path), length(Z), ncol(covars) - 1 - length(Z),
    fit_stats$oos_r2[fit_stats$outcome == "Y_gdp_growth"],
    fit_stats$oos_r2[fit_stats$outcome == "Y_inflation"],
    fit_stats$oos_r2[fit_stats$outcome == "Y_polyarchy"],
    median(fit_stats$oos_r2, na.rm = TRUE)
  ))
}
