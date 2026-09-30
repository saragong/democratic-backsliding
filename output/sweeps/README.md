# Sweeps

Output that compares ACROSS specs, or is not an RD run at all. A sweep over a
single spec's windows and treatments lives in that spec's folder instead, as
`output/runs/<spec>/pooled/`; see `output/runs/README.md`.

Each folder is created by `sweep_dir()` in `scripts/rdd_helpers.R`, and most
record what they held fixed and what they varied in a `sweep_config.csv`.
`<suffix>` is the build suffix: `_exclyr` under the project's window
convention, plus `_pre` for the pre-election placebo.

| Folder | Script | What it holds |
|---|---|---|
| `restriction_grid_<instr>_w<N>_trt-<treatment><suffix>/` | 13 | to-do 1: marginal tables and colour-coded pair grids over the five sample-restriction axes |
| `sample_composition_<instr>_w<N><suffix>/` | 13 | who each restriction axis drops, independent of treatment and outcome |
| `alt_specs_<sample><suffix>/` | 15 | to-dos 4-6: instrument x treatment-definition grids, plus what the DDCG extension actually adds |
| `instrument_overlap<suffix>/` | 16 | to-do 7: UpSet plots and Jaccard heatmaps measuring how much the instruments really differ |
| `party_outcome_rdd_<instr>_w<N><suffix>/` | 17 | 1a: the same RD with the winner's OTHER party scores as the outcome — is a narrow anti-pluralist victory also a populist / left / minority-hostile one? |
| `vparty_jaccard_panels/` | 18 | 1c: the anti-pluralism x populism Jaccard heatmap over the full V-Party dataset, as a 2 x 5 OECD-by-decade grid |
| `vparty_ideology_quadrants/` | 19 | 1d: where the most common Wikipedia/Wikidata ideology tags sit on the illiberalism x populism plane, same 2 x 5 grid |
| `populist_threshold/` | 01g | 2b: the PopuList-calibrated cutoff for "illiberal", and the ROC it comes from |
| `cell_rdd_<instr><suffix>/` | 21 | reduced-form RD on the decade x OECD grid, every outcome, w1-10, full and PopuList-restricted samples |

`output/builds/<instrument>_w<N><suffix>/` (written by `11_build_rdd_data.R`)
holds each build's ERT-episode match accounting: which episodes the election
spine can and cannot reach, and why.
