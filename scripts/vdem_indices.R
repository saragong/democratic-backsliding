# ==============================================================================
# The V-Dem democracy indices used as RDD outcomes, defined once.
#
# Pure data: no libraries, no side effects, safe to source from anywhere. It
# exists because the same list is needed at three layers that must not depend
# on each other:
#
#   01f_load_extra_outcomes.R  pulls these columns out of ERT::vdem
#   02a_build_panel.R          joins them into combined_panel.rds
#   rdd_helpers.R              labels them and groups them into plot panels
#                              (11_build_rdd_data.R then differences each one
#                              over the window to make its Y_ column)
#
# Keeping the list in the loader would make the analysis layer guess at what is
# available; keeping it in rdd_helpers.R would make a data loader depend on the
# RDD's own helper. So it lives on its own, and adding an index is one line
# here plus one panel entry in rdd_helpers.R.
#
# GROUPING follows V-Dem's OWN taxonomy -- the five high-level indices, then
# the mid-level indices that aggregate into each of them. That is a principled
# grouping rather than a convenience one, and it happens to keep every panel at
# or under five series, which is the most build_panel_plot()'s palette carries.
#
# DIRECTION: every index here runs HIGHER = MORE DEMOCRATIC, the opposite of
# the economic outcomes. A NEGATIVE RD estimate is the "backsliding" sign.
# ==============================================================================

# Y_ column name -> the V-Dem variable it is the window-change of.
# v2x_polyarchy is deliberately absent from the loader's pull: it already
# reaches combined_panel.rds through the ERT episode data (01a), so pulling it
# again in 01f would create a duplicate column on the join.
VDEM_INDEX_VARS <- c(
  # --- high level: the five headline democracy indices ---
  Y_polyarchy  = "v2x_polyarchy",
  Y_libdem     = "v2x_libdem",
  Y_partipdem  = "v2x_partipdem",
  Y_delibdem   = "v2x_delibdem",
  Y_egaldem    = "v2x_egaldem",
  # --- mid level: electoral (aggregate into polyarchy) ---
  Y_elecoff    = "v2x_elecoff",
  Y_frefair    = "v2xel_frefair",
  Y_frassoc    = "v2x_frassoc_thick",
  Y_suffr      = "v2x_suffr",
  Y_freexp     = "v2x_freexp_altinf",
  # --- mid level: liberal ---
  Y_cl_rol     = "v2xcl_rol",
  Y_jucon      = "v2x_jucon",
  Y_legcon     = "v2xlg_legcon",
  # --- mid level: participatory ---
  Y_cspart     = "v2x_cspart",
  Y_dd         = "v2xdd_dd",
  Y_locelec    = "v2xel_locelec",
  Y_regelec    = "v2xel_regelec",
  # --- mid level: deliberative and egalitarian ---
  Y_delib      = "v2xdl_delib",
  Y_eqprotec   = "v2xeg_eqprotec",
  Y_eqaccess   = "v2xeg_eqaccess",
  Y_eqdr       = "v2xeg_eqdr"
)

# Already in combined_panel.rds before this registry existed, from other
# sources, so the loader must NOT pull them again:
#   v2x_polyarchy            via 01a (ERT episode data)
#   v2x_jucon, v2xlg_legcon  via 01f's original executive-constraints block
VDEM_INDEX_ALREADY_IN_PANEL <- c(
  "v2x_polyarchy", "v2x_jucon", "v2xlg_legcon"
)

# What 01f must newly pull from ERT::vdem and 02a must newly join.
VDEM_INDEX_NEW_VARS <- setdiff(
  unname(VDEM_INDEX_VARS), VDEM_INDEX_ALREADY_IN_PANEL
)
