# ==============================================================================
# Event study plots: mean change in economic outcomes around episode starts
# Estimator: OLS with event-time dummies, SEs clustered by country (fixest)
# Baseline: t = -1 (normalized to 0 by construction in event_study_panel.rds)
# Output: figures/event_study_democratization.pdf
#         figures/event_study_autocratization.pdf
# ==============================================================================

library(tidyverse)
library(fixest)
library(patchwork)
library(here)

data_dir <- here::here("data")
fig_dir  <- here::here("figures")
dir.create(fig_dir, showWarnings = FALSE)

event_study <- readRDS(file.path(data_dir, "event_study_panel.rds"))

# --- Config -------------------------------------------------------------------

OUTCOMES <- c(
  "gdp_pc_growth_chg",
  "cpi_inflation_chg",
  "unemployment_rate_chg",
  "trade_pct_gdp_chg",
  "top10_share_chg",
  "gini_disp_chg"
)

OUTCOME_LABELS <- c(
  gdp_pc_growth_chg     = "GDP growth (pp vs t=-1)",
  cpi_inflation_chg     = "Inflation rate (pp vs t=-1)",
  unemployment_rate_chg = "Unemployment rate (pp vs t=-1)",
  trade_pct_gdp_chg     = "Trade/GDP (pp vs t=-1)",
  top10_share_chg       = "Top 10% share (pp vs t=-1)",
  gini_disp_chg         = "Gini, disposable (pp vs t=-1)"
)

EP_COLORS <- c(
  "Democratization" = "#2ca25f",
  "Autocratization" = "#e34a33"
)

# --- Estimate -----------------------------------------------------------------
# OLS of outcome_chg on event-time dummies, reference = t=-1, clustered by country.
# Returns a tidy tibble with one row per event_time including the t=-1 zero row.

estimate_es <- function(outcome_var, ep_type_val) {
  d <- event_study |>
    filter(ep_type == ep_type_val, !is.na(.data[[outcome_var]]))

  if (n_distinct(d$country_text_id) < 5) return(NULL)

  fit <- feols(
    as.formula(paste0(outcome_var, " ~ i(event_time, ref = -1)")),
    cluster = ~country_text_id,
    data    = d,
    warn    = FALSE,
    notes   = FALSE
  )

  tidy_res <- broom::tidy(fit, conf.int = TRUE, conf.level = 0.95) |>
    mutate(
      event_time = as.integer(str_extract(term, "-?\\d+")),
      ep_type    = ep_type_val,
      outcome    = outcome_var
    ) |>
    select(ep_type, outcome, event_time, estimate, conf.low, conf.high)

  # Re-attach the t=-1 reference row (zero by construction)
  bind_rows(
    tibble(ep_type = ep_type_val, outcome = outcome_var,
           event_time = -1L, estimate = 0, conf.low = 0, conf.high = 0),
    tidy_res
  ) |>
    arrange(event_time)
}

results <- expand_grid(
  outcome = OUTCOMES,
  ep_type = c("Democratization", "Autocratization")
) |>
  pmap(\(outcome, ep_type) estimate_es(outcome, ep_type)) |>
  bind_rows()

# --- Single outcome panel -----------------------------------------------------

plot_es_panel <- function(outcome_var, ep_type_val, show_x = FALSE) {
  d   <- results |> filter(outcome == outcome_var, ep_type == ep_type_val)
  col <- EP_COLORS[[ep_type_val]]

  n_ep <- event_study |>
    filter(ep_type == ep_type_val, event_time == 0L,
           !is.na(.data[[outcome_var]])) |>
    nrow()

  ggplot(d, aes(x = event_time, y = estimate)) +
    # Shade the post-episode region
    annotate("rect", xmin = -0.5, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = col, alpha = 0.04) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                fill = col, alpha = 0.20) +
    geom_line(color = col, linewidth = 0.7, na.rm = TRUE) +
    geom_point(color = col, size = 1.0, na.rm = TRUE) +
    geom_vline(xintercept = -0.5, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    scale_x_continuous(breaks = seq(-15, 15, by = 5)) +
    labs(
      y       = OUTCOME_LABELS[[outcome_var]],
      x       = if (show_x) "Years relative to episode start" else NULL,
      caption = if (show_x)
        sprintf("n = %d episodes with outcome data at t = 0", n_ep)
      else NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.title.x     = if (show_x) element_text(size = 8) else element_blank(),
      axis.text.x      = if (show_x) element_text(size = 8) else element_blank(),
      axis.ticks.x     = if (show_x) element_line() else element_blank(),
      axis.title.y     = element_text(angle = 0, hjust = 1, vjust = 0.5,
                                      size = 7, color = "grey30"),
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(size = 6, color = "grey50", hjust = 0)
    )
}

# --- Full figure for one episode type -----------------------------------------

make_es_figure <- function(ep_type_val) {
  panels <- imap(OUTCOMES, \(v, i)
    plot_es_panel(v, ep_type_val, show_x = (i == length(OUTCOMES)))
  )

  wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title    = sprintf("Economic outcomes around %s episodes", tolower(ep_type_val)),
      subtitle = "OLS; mean change relative to t=-1; 95% CI clustered by country",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey40")
      )
    )
}

# --- Save ---------------------------------------------------------------------

for (ep in c("Democratization", "Autocratization")) {
  path <- file.path(fig_dir, sprintf("event_study_%s.pdf", tolower(ep)))
  ggsave(path, make_es_figure(ep), width = 7, height = 11)
  message("Saved ", path)
}
