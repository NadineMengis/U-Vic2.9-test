diago.f
size.h
param.h
pconst.h
stdunits.h
calendar.h
ctmb.h
coord.h
cprnts.h
cregin.h
diag.h
docnam.h
emode.h
grdvar.h
iounit.h
levind.h
mw.h
scalar.h
stab.h
state.h
switch.h
tmngr.h
vmixc.h
cembm.h
ctavg.h
csbc.h
diaga.h
npzd.h
dens.h
#if defined O_mom
# if defined O_embm
# endif
# if defined O_tracer_averages
# endif
# if defined O_save_carbon_carbonate_chem
# endif
# if defined O_save_carbon_totals
# endif
# if defined O_mom_tbt
# endif
# if defined O_diagnostic_surf_height && defined O_stream_function
#  if !defined O_rigid_lid_surface_pressure && !defined O_implicit_free_surface
#  endif
# endif
# if defined O_tracer_averages
# endif
# if defined O_time_step_monitor
#  if defined O_save_carbon_totals
#  else
#  endif
# endif
# if defined O_stability_tests
# endif
# if defined O_energy_analysis
# endif
# if defined O_term_balances
# endif
# if defined O_gyre_components
# endif
# if defined O_meridional_overturning
# endif
# if defined O_show_external_mode
#  if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
#   if defined O_rigid_lid_surface_pressure
#   endif
#   if defined O_implicit_free_surface
#   endif
#   if defined O_rigid_lid_surface_pressure
#   endif
#   if defined O_implicit_free_surface
#   endif
#  endif
#  if defined O_stream_function
#  endif
# endif
# if defined O_meridional_tracer_budget
# endif
# if defined O_show_zonal_mean_of_sbc
# endif
# if defined O_time_averages
# endif
# if defined O_xbts
# endif
# if defined O_mom_tbt
# endif
# if defined O_save_carbon_carbonate_chem
# endif
#endif
#if defined O_mom && defined O_time_step_monitor
# if defined O_npzd
# endif
# if defined O_tai_otsf
# endif
# if defined O_tai_otsf_from_averages
# endif
# if defined O_tai_slh_from_averages
# endif
# if defined O_tai_otsf && !defined O_tai_otsf_from_averages
# endif
# if defined O_tai_slh && !defined O_tai_slh_from_averages
# endif
# if defined O_carbon
#  if defined O_carbon_14
#  endif
# endif
# if defined O_npzd_alk
# endif
# if defined O_npzd_o2
# endif
# if defined O_npzd
#  if defined O_npzd_nitrogen
#  endif
# endif
# if defined O_cfcs_data_transient
# endif
# if defined O_save_carbon_carbonate_chem
# endif
# if defined O_tai_otsf_from_averages
# else
# endif
# if defined O_tai_slh_from_averages
# else
# endif
#endif
#if defined O_mom && defined O_tai_slh
#endif
