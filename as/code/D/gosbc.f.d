gosbc.f
size.h
param.h
pconst.h
stdunits.h
calendar.h
csbc.h
grdvar.h
tmngr.h
switch.h
cembm.h
atm.h
mw.h
ice.h
mtlm.h
levind.h
sed.h
#if defined O_mom || defined O_embm
# if defined O_embm
# endif
# if defined O_mom
# endif
# if defined O_ice
#  if defined O_ice_cpts
#  endif
# endif
# if defined O_mtlm
# endif
# if defined O_fwa 
# endif
# if defined O_sed && !defined O_sed_uncoupled
# endif
# if defined O_sulphate_data_transient
# endif
# if defined O_embm_read_sflx || defined O_embm_write_sflx
# endif
# if defined O_mom && defined O_embm
#  if defined O_plume
#   if !defined O_plume_brine
#   endif
#  endif
#  if defined O_convect_brine
#  endif
#  if defined O_plume
#  endif
#  if defined O_convect_brine
#  endif
#  if defined O_ice_evp || defined O_embm_awind
#  endif
#  if defined O_convect_brine
#  endif
# endif
# if defined O_embm
# endif
# if defined O_mtlm
#  if defined O_carbon
#  endif
#  if defined O_mtlm && defined O_carbon
#  endif
#  if defined O_carbon
#  endif
# else
# endif
# if defined O_embm
#  if defined O_sed && !defined O_sed_uncoupled
#   if defined O_global_sums
#   endif
#   if defined O_carbon
#    if defined O_global_sums
#    endif
#   endif
#   if defined O_npzd_alk
#   endif
#  endif
# endif
# if defined O_sealev_data_transient && defined O_sealev_salinity
# endif
# if defined O_mom && defined O_embm
#  if defined O_fwa
#   if defined O_fwa_precip
#   endif
#   if defined O_fwa_precip
#   endif
#   if defined O_fwa_compevap
#   endif
#   if defined O_fwa_compevap
#   endif
#  endif
#  if defined O_carbon || defined O_npzd_alk || defined O_npzd_o2 || defined O_npzd || defined O_cfcs_data_transient
#   if defined O_carbon
#    if defined O_carbon_14
#    endif
#   endif
#   if defined O_npzd_alk
#   endif
#   if defined O_npzd_o2
#   endif
#   if defined O_npzd
#    if !defined O_npzd_no_vflux
#    endif
#    if defined O_npzd_nitrogen
#     if !defined O_npzd_no_vflux
#     endif
#    endif
#   endif
#   if defined O_cfcs_data_transient
#   endif
#   if defined O_carbon
#    if defined O_carbon_14
#    endif
#   endif
#   if defined O_npzd_alk
#   endif
#   if defined O_npzd_o2
#   endif
#   if defined O_npzd
#    if !defined O_npzd_no_vflux
#    endif
#    if defined O_npzd_nitrogen
#     if !defined O_npzd_no_vflux
#     endif
#    endif
#   endif
#   if defined O_cfcs_data_transient
#   endif
#  endif
# endif
# if defined O_embm_read_sflx
# elif defined O_embm_write_sflx
# endif
#endif
# if defined O_sulphate_data_transient
# endif
#if defined O_embm_read_sflx
#endif
#if defined O_embm_write_sflx
#endif
