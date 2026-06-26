# ==============================================================================
# Country profile plots: democracy series + economic outcomes — all countries.
# One page per country, sorted alphabetically. See plot_helpers.R for layout.
# Data source: combined_panel.rds (built by 02a_build_panel.R).
# ==============================================================================

library(tidyverse)
library(patchwork)

source(here::here("scripts", "plot_helpers.R"))

data_dir <- here::here("data")
fig_dir <- here::here("figures")
dir.create(fig_dir, showWarnings = FALSE)

# --- Load data ----------------------------------------------------------------

panel <- readRDS(file.path(data_dir, "combined_panel.rds"))

# --- Country list: ERT countries with at least one episode, sorted by name ---

countries <- panel |>
  filter(dem_ep == 1 | aut_ep == 1) |>
  distinct(country_text_id, country_name) |>
  arrange(country_name)

message(sprintf("Generating plots for %d countries ...", nrow(countries)))

# --- Main loop ----------------------------------------------------------------

out_path <- file.path(fig_dir, "country_profiles.pdf")
pdf(out_path, width = 11, height = 10, onefile = TRUE)

for (i in seq_len(nrow(countries))) {
  code <- countries$country_text_id[i]
  panel_c <- panel |> filter(country_text_id == code) |> arrange(year)

  print(make_country_plot(panel_c))

  if (i %% 20 == 0 || i == nrow(countries)) {
    message(sprintf(
      "  %d / %d  (%s)",
      i,
      nrow(countries),
      countries$country_name[i]
    ))
  }
}

dev.off()
message("Saved ", out_path)
