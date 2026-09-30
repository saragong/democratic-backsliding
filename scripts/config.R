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
