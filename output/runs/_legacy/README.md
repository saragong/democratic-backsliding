# Legacy RDD output

Everything in this folder predates the run-folder layout introduced with
`scripts/rdd_helpers.R`. It was produced by hand-editing `SCORE_GAP_MIN` /
`ILLIBERAL_CUTOFF` at the top of `12_rdd_analysis.R` and re-running, with each
run distinguished only by a filename suffix.

It is kept for reference only. Current output lives in `output/runs/<slug>/`,
one folder per configuration, indexed by `output/runs/manifest.csv`.

Two caveats if you compare these numbers against a current run:

1. The committed `11_build_rdd_data.R` that produced some of these files had a
   join bug (`country_score_sd` grouped by `(country_text_id, party_id)` but
   joined on `country_text_id` alone), which fanned each election out against
   every party in its country. Files written while that bug was live are based
   on an inflated 6,602-row sample rather than the correct 1,347 elections.
2. The outcome set has since grown from 7 variables to 18.

3. The treatment window now opens in the election year by default
   (`TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR = TRUE`), where this output used
   the year after. That is worth 31 treated elections on the ERT arm
   (142 -> 173) and it moves the estimates.

All three are reproducible: set `TREATMENT_WINDOW_INCLUDES_ELECTION_YEAR` to
`FALSE` and the unrestricted baseline first stage comes back as
`-0.025 (0.058)`, bw 21.08, and the `illiberal_score > 0.6` run as
`0.181 (0.103)*`, bw 20.86 -- both exact. Those runs land in folders with an
`_exclyr` suffix.
