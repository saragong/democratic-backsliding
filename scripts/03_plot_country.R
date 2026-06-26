# ==============================================================================
# Country profile plot: democracy series + economic outcomes, one country.
# Top panel: V-Dem EDI (grey), ERT episodes (colored), DDCG dem (dashed).
# Lower panels: GDP growth, inflation, unemployment, trade, top-10% share, Gini.
# Data source: combined_panel.rds (built by 02a_build_panel.R).
# ==============================================================================

library(tidyverse)
library(patchwork)

source(here::here("scripts", "plot_helpers.R"))

data_dir <- here::here("data")
fig_dir <- here::here("figures")
dir.create(fig_dir, showWarnings = FALSE)

COUNTRY <- "USA" # ERT country_text_id (standard ISO-3)

# --- Load and filter ----------------------------------------------------------

panel <- readRDS(file.path(data_dir, "combined_panel.rds"))
panel_c <- panel |> filter(country_text_id == COUNTRY) |> arrange(year)

stopifnot(nrow(panel_c) > 0)

# --- Plot and save ------------------------------------------------------------

out_path <- file.path(
  fig_dir,
  sprintf("country_profile_%s.pdf", tolower(COUNTRY))
)
ggsave(out_path, make_country_plot(panel_c), width = 11, height = 10)
message("Saved ", out_path)
