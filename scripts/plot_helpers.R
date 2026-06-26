# ==============================================================================
# Shared plot helpers: sourced by 03_plot_country.R and 04_plot_all_countries.R
# ==============================================================================

# Country code crosswalks (ERT <-> DDCG wbcode <-> ISO 3-letter)
ert_to_ddcg <- c(
  COD = "ZAR",
  ROU = "ROM",
  SGP = "SIN",
  TWN = "TAW",
  ARE = "UAE"
)
ddcg_to_iso3 <- c(
  ZAR = "COD",
  ROM = "ROU",
  SIN = "SGP",
  TAW = "TWN",
  UAE = "ARE"
)

EP_COLORS <- c("Democratization" = "#2ca25f", "Autocratization" = "#e34a33")

make_ep_bands <- function(panel_c) {
  # xmin: when the leader in power at episode start took office (falls back to
  # the episode start year itself when NED executive_entry_year is unavailable).
  # xmax: last year of the episode (unchanged).
  bind_rows(
    panel_c |>
      filter(dem_ep == 1) |>
      arrange(year) |>
      group_by(ep_id = dem_ep_id) |>
      summarise(
        xmin = coalesce(first(executive_entry_year), first(year)),
        xmax = max(year),
        ep_type = "Democratization",
        .groups = "drop"
      ),
    panel_c |>
      filter(aut_ep == 1) |>
      arrange(year) |>
      group_by(ep_id = aut_ep_id) |>
      summarise(
        xmin = coalesce(first(executive_entry_year), first(year)),
        xmax = max(year),
        ep_type = "Autocratization",
        .groups = "drop"
      )
  )
}

# Assign one qualitative color per unique Archigos leader that actually appears in
# the rug.
# Order is chronological (first year of appearance) so the legend reads
# left-to-right in time order within each row.
assign_leader_colors <- function(panel_c) {
  parties <- panel_c |>
    arrange(year) |>
    mutate(
      leader = archigos_leader
    ) |>
    filter(!is.na(leader)) |>
    distinct(leader) |>
    pull(leader)
  if (length(parties) == 0) {
    return(setNames(character(0), character(0)))
  }
  setNames(scales::hue_pal()(length(parties)), parties)
}

# Thin horizontal legend showing one colored square + name per ruling leader.
# Uses the same x coordinate system and axis structure as all other panels
# so patchwork aligns it correctly with the data panels above and below.
make_ruler_legend <- function(leader_colors, x_lim, has_ddcg = FALSE) {
  sec <- if (has_ddcg) dup_axis(labels = NULL, name = NULL) else waiver()

  base <- ggplot() +
    scale_x_continuous(breaks = seq(1900, x_lim[2], by = 5)) +
    scale_y_continuous(
      name = NULL,
      limits = c(0, 1),
      breaks = NULL,
      sec.axis = sec
    ) +
    coord_cartesian(xlim = x_lim) +
    theme_bw(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      panel.grid = element_blank()
    )

  if (length(leader_colors) == 0) {
    return(base)
  }

  n <- length(leader_colors)
  n_cols <- min(n, 7L)
  n_rows <- ceiling(n / n_cols)

  year_span <- x_lim[2] - x_lim[1]
  step <- year_span / n_cols

  ld <- tibble(
    leader = str_trunc(names(leader_colors), 18, ellipsis = "..."),
    color = unname(leader_colors),
    col = (seq_len(n) - 1L) %% n_cols,
    row = (seq_len(n) - 1L) %/% n_cols
  ) |>
    mutate(
      x_sq = x_lim[1] + (col + 0.3) * step,
      x_lbl = x_lim[1] + (col + 0.5) * step,
      y = 1 - (row + 0.5) / n_rows
    )

  base +
    geom_point(
      data = ld,
      aes(x = x_sq, y = y, color = I(color)),
      shape = 15,
      size = 3.5,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = ld,
      aes(x = x_lbl, y = y, label = leader),
      hjust = 0,
      vjust = 0.5,
      size = 2.1,
      color = "grey20",
      inherit.aes = FALSE
    )
}

# has_ddcg: when TRUE, adds a blank right axis to match p_dem's right margin,
# keeping x-axes aligned across all panels in the patchwork stack.
make_outcome_panel <- function(
  panel_c,
  y_var,
  y_label,
  ep_bands,
  x_lim,
  has_ddcg = FALSE,
  show_x = FALSE
) {
  ggplot() +
    geom_rect(
      data = ep_bands,
      aes(xmin = xmin, xmax = xmax, fill = ep_type),
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.12,
      inherit.aes = FALSE
    ) +
    geom_line(
      data = panel_c |> filter(!is.na(.data[[y_var]])),
      mapping = aes(x = year, y = .data[[y_var]]),
      color = "grey30",
      linewidth = 0.5
    ) +
    scale_fill_manual(values = EP_COLORS, guide = "none") +
    scale_x_continuous(breaks = seq(1900, x_lim[2], by = 5)) +
    scale_y_continuous(
      name = y_label,
      n.breaks = 3,
      sec.axis = if (has_ddcg) {
        dup_axis(labels = NULL, name = NULL)
      } else {
        waiver()
      }
    ) +
    coord_cartesian(xlim = x_lim) +
    theme_bw(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = if (show_x) {
        element_text(angle = 45, hjust = 1, size = 8)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.title.y = element_text(
        angle = 0,
        hjust = 1,
        vjust = 0.5,
        size = 7,
        color = "grey40"
      ),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_blank(),
      panel.grid.minor = element_blank()
    )
}

make_country_plot <- function(panel_c) {
  has_ddcg <- any(!is.na(panel_c$dem))
  n_dem_eps <- n_distinct(panel_c$dem_ep_id[panel_c$dem_ep == 1], na.rm = TRUE)
  n_aut_eps <- n_distinct(panel_c$aut_ep_id[panel_c$aut_ep == 1], na.rm = TRUE)
  cname <- panel_c$country_name[1]
  code <- panel_c$country_text_id[1]
  ddcg_c <- panel_c |> filter(!is.na(dem)) |> select(year, dem)
  ep_bands <- make_ep_bands(panel_c)
  x_lim <- c(min(panel_c$year), as.integer(format(Sys.Date(), "%Y")))

  leader_colors <- assign_leader_colors(panel_c)

  # Contiguous leader/leader spells for the rug: one row per spell so that
  # geom_rect borders appear only at actual transitions, not between every year.
  # Use Archigos leader only
  rug_periods <- panel_c |>
    arrange(year) |>
    mutate(leader = archigos_leader) |>
    filter(!is.na(leader)) |>
    mutate(
      new_spell = is.na(lag(leader)) | leader != lag(leader),
      spell_id = cumsum(new_spell)
    ) |>
    group_by(spell_id, leader) |>
    summarise(xmin = min(year) - 0.5, xmax = max(year) + 0.5, .groups = "drop")

  # Episode shading and ruler rug both use fill, so they share one combined scale
  all_fill <- c(EP_COLORS, leader_colors)

  df_dem_ep <- panel_c |>
    filter(dem_ep == 1) |>
    mutate(ep_type = "Democratization", ep_id = dem_ep_id)
  df_aut_ep <- panel_c |>
    filter(aut_ep == 1) |>
    mutate(ep_type = "Autocratization", ep_id = aut_ep_id)

  p_dem <- ggplot() +
    # Episode shading — capped at y=[0,1] so the rug area below stays unshaded
    geom_rect(
      data = ep_bands,
      aes(xmin = xmin, xmax = xmax, fill = ep_type),
      ymin = 0,
      ymax = 1,
      alpha = 0.12,
      inherit.aes = FALSE
    ) +
    # Ruling-leader rug: one rect per contiguous spell, black border at transitions
    geom_rect(
      data = rug_periods,
      aes(xmin = xmin, xmax = xmax, fill = leader),
      ymin = -0.08,
      ymax = 0,
      color = "black",
      linewidth = 0.2,
      inherit.aes = FALSE
    ) +
    geom_line(
      data = panel_c,
      aes(x = year, y = v2x_polyarchy),
      color = "grey65",
      linewidth = 0.55
    ) +
    geom_line(
      data = df_dem_ep,
      aes(x = year, y = v2x_polyarchy, group = ep_id, color = ep_type),
      linewidth = 1.8
    ) +
    geom_line(
      data = df_aut_ep,
      aes(x = year, y = v2x_polyarchy, group = ep_id, color = ep_type),
      linewidth = 1.8
    ) +
    {
      if (has_ddcg) {
        geom_step(
          data = ddcg_c,
          aes(x = year, y = dem),
          color = "black",
          linetype = "dashed",
          linewidth = 0.7,
          na.rm = TRUE
        )
      }
    } +
    scale_fill_manual(values = all_fill, guide = "none") +
    scale_color_manual(name = "ERT episode", values = EP_COLORS) +
    # No explicit limits so the rug tiles below 0 are not clipped by the scale;
    # coord_cartesian controls the visible range. Breaks stay at 0–1.
    scale_y_continuous(
      name = "V-Dem EDI",
      breaks = seq(0, 1, by = 0.25),
      sec.axis = if (has_ddcg) {
        sec_axis(~., name = "DDCG dem", breaks = c(0, 1))
      } else {
        waiver()
      }
    ) +
    scale_x_continuous(breaks = seq(1900, x_lim[2], by = 5)) +
    coord_cartesian(xlim = x_lim, ylim = c(-0.09, 1)) +
    labs(
      title = sprintf(
        "%s (%s)  --  democratization: %d  |  autocratization: %d",
        cname,
        code,
        n_dem_eps,
        n_aut_eps
      ),
      caption = paste(
        "Solid grey = V-Dem EDI; colored segments = ERT episodes; rug = ruling leader (see legend).",
        if (has_ddcg) {
          sprintf(
            "Dashed = DDCG dem (right axis; %d-%d).",
            min(ddcg_c$year),
            max(ddcg_c$year)
          )
        } else {
          "Country not in DDCG."
        }
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.caption = element_text(size = 7, color = "grey50", hjust = 0),
      axis.title.y.left = element_text(
        angle = 0,
        hjust = 1,
        vjust = 0.5,
        color = "grey40",
        size = 7
      ),
      axis.title.y.right = element_text(
        angle = 0,
        hjust = 0,
        vjust = 0.5,
        color = "black",
        size = 7
      ),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  n_legend_rows <- max(1L, ceiling(length(leader_colors) / 7L))
  legend_height <- 0.35 * n_legend_rows

  p_legend <- make_ruler_legend(leader_colors, x_lim, has_ddcg)

  p_gdp <- make_outcome_panel(
    panel_c,
    "gdp_pc_growth",
    "GDP growth",
    ep_bands,
    x_lim,
    has_ddcg
  )
  p_cpi <- make_outcome_panel(
    panel_c,
    "cpi_inflation",
    "Inflation",
    ep_bands,
    x_lim,
    has_ddcg
  )
  p_unemp <- make_outcome_panel(
    panel_c,
    "unemployment_rate",
    "Unemp.",
    ep_bands,
    x_lim,
    has_ddcg
  )
  p_trade <- make_outcome_panel(
    panel_c,
    "trade_pct_gdp",
    "Trade/GDP",
    ep_bands,
    x_lim,
    has_ddcg
  )
  p_top10 <- make_outcome_panel(
    panel_c,
    "top10_share",
    "Top 10% Share",
    ep_bands,
    x_lim,
    has_ddcg
  )
  p_gini <- make_outcome_panel(
    panel_c,
    "gini_disp",
    "Gini",
    ep_bands,
    x_lim,
    has_ddcg,
    show_x = TRUE
  )

  (p_dem / p_legend / p_gdp / p_cpi / p_unemp / p_trade / p_top10 / p_gini) +
    plot_layout(heights = c(2.5, legend_height, 1, 1, 1, 1, 1, 1))
}
