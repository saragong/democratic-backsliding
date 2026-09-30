# ==============================================================================
# Project-wide defaults for the fuzzy-RDD pipeline (scripts 11-22, adhoc/).
#
# Each default is defined HERE, once. Scripts still guard their toggles with
# `if (!exists("X")) X <- DEFAULT_X`, so a driver that sys.source()s a child
# into an environment with the toggle already set still overrides it; what
# this file removes is a dozen scripts each hard-coding their own copy of the
# fallback and drifting apart. (They had: running 12_rdd_analysis.R alone used
# score_gap_z >= 0 and illiberal_score > 0.6, while every driver's baseline was
# the unrestricted sample.)
#
# Driver-specific names (SWEEP_INSTRUMENT, CELL_INSTRUMENT, SPLIT_INSTRUMENT,
# ...) are kept as they are; they default from these same constants.
#
# Sourced by rdd_helpers.R, not run.
# ==============================================================================

# Which V-Party score picks the "illiberal" member of the top 2.
DEFAULT_INSTRUMENT <- "v2xpa_antiplural"

# Post-election window N, in years, for scripts that run a single window.
DEFAULT_WINDOW <- 5

# Every window a sweep over N runs.
DEFAULT_WINDOWS <- 1:10

# The treatment the fuzzy RD instruments for, and the set a window sweep runs.
DEFAULT_TREATMENT <- "backsliding_Nyr"
DEFAULT_SWEEP_TREATMENTS <- c(
  "backsliding_Nyr",
  "backsliding_union_Nyr",
  "polyarchy_decline"
)

# Sample restrictions. -Inf is the no-op for the two floors and +Inf the no-op
# for the ceiling, so the default is the full scored sample.
DEFAULT_SCORE_GAP_MIN <- -Inf
DEFAULT_ILLIBERAL_CUTOFF <- -Inf
DEFAULT_OTHER_CUTOFF_MAX <- Inf

# The treatment window opens the year AFTER the election. Project convention;
# do not change. See 11_build_rdd_data.R for the argument either way.
DEFAULT_INCL_ELECTION_YEAR <- FALSE

# The pre-election placebo window. Off except when a placebo is asked for.
DEFAULT_PLACEBO <- FALSE

# Left-right gap between the two top-2 parties, |v2pariglef| difference on
# V-Party's expert scale. A floor (keep pairs at least this far apart, for the
# left/right split) and a ceiling (keep pairs at most this far apart, the
# "no meaningful left-right difference" placebo). The no-ops keep every pair.
DEFAULT_LR_GAP_MIN <- -Inf
DEFAULT_LR_GAP_MAX <- Inf

# Keep only pairs on opposite sides of the left-right centre: one top-2 party
# with v2pariglef < 0 and the other > 0.
DEFAULT_LR_STRADDLE <- FALSE

# Drop elections whose [election_year, election_year + N] window overlaps a
# Funke, Schularick & Trebesch populist-leader spell.
DEFAULT_EXCLUDE_FUNKE <- FALSE

# Covariate adjustment for the RD (Calonico, Cattaneo, Farrell & Titiunik
# 2019): "none", or "lp" for the leave-country-out local projection of each
# outcome built by 11b_build_covariates.R, entered linearly via covs=.
DEFAULT_RD_COVARIATES <- "none"
