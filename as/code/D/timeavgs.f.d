timeavgs.f
size.h
stdunits.h
param.h
pconst.h
levind.h
cregin.h
timeavgs.h
npzd.h
coord.h
docnam.h
csbc.h
emode.h
grdvar.h
iounit.h
scalar.h
switch.h
tmngr.h
calendar.h
cembm.h
atm.h
#if defined O_mom
# if defined O_time_averages
#  if defined O_npzd
#  endif
#  if defined O_embm
#  endif
#  if defined O_embm && !defined O_save_fluxes_with_virtual
#  else
#  endif
#  if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
#  endif
#  if defined O_stream_function
#  endif
#  if defined O_save_convection
#  endif
#  if defined O_save_carbon_carbonate_chem
#  endif
#  if defined O_gent_mcwilliams
#  endif
#  if defined O_carbon && defined O_carbon_14
#  endif
#  if defined O_save_kv
#  endif
#  if defined O_save_npzd
#   if defined O_npzd_nitrogen
#   endif
#   if defined O_npzd_nitrogen
#   endif
#  endif
#  if defined O_save_convection
#  endif
#  if defined O_save_carbon_carbonate_chem
#  endif
#  if defined O_gent_mcwilliams
#  endif
#  if defined O_carbon && defined O_carbon_14
#  endif
#  if defined O_save_kv
#  endif
#  if defined O_save_npzd
#   if defined O_npzd_nitrogen
#   endif
#  endif
#  if defined O_save_convection
#  endif
#  if defined O_save_carbon_carbonate_chem
#  endif
#  if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
#  endif
#  if defined O_stream_function
#  endif
#  if defined O_gent_mcwilliams
#  endif
#  if defined O_carbon && defined O_carbon_14
#  endif
#  if defined O_save_kv
#  endif
#  if defined O_save_npzd
#   if defined O_npzd_nitrogen
#   endif
#  endif
#  if defined O_save_convection
#  endif
#  if defined O_save_carbon_carbonate_chem
#  endif
#  if defined O_save_convection
#  endif
#  if defined O_save_carbon_carbonate_chem
#  endif
#  if defined O_save_convection
#  endif
#  if defined O_save_carbon_carbonate_chem
#  endif
#  if defined O_gent_mcwilliams
#  endif
#  if defined O_save_npzd
#   if defined O_npzd_nitrogen
#   endif
#  endif
#  if defined O_save_convection
#  endif
#  if defined O_save_carbon_carbonate_chem
#  endif
# endif
#endif
