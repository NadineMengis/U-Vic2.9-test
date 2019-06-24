global_sums.f
#if defined O_global_sums || defined O_co2emit_diag
# if defined O_embm
# endif
# if defined O_ice && defined O_embm
#  if defined O_ice_cpts
#  endif
# endif
# if defined O_mom
#  if defined O_npzd
#  endif
# endif
# if defined O_mtlm
# endif
# if defined O_sed
# endif
# if defined O_ism
# endif
# if !defined O_embm
# endif
# if !defined O_sed
# endif
# if defined O_mtlm
#  if defined O_mtlm_segday
#  else
#  endif
# endif
# if defined O_mom
# endif
# if defined O_embm
#  if defined O_carbon_co2_2d
#  else
#  endif
# endif
# if defined O_ice && defined O_embm
#  if defined O_ism
#  endif          
#  if defined O_ice_cpts
#  endif
# endif
# if defined O_mtlm && defined O_embm
# endif
# if defined O_embm
#  if defined O_mtlm
#   if defined O_carbon
#   endif
#  endif
# endif
# if defined O_mom
#  if defined O_carbon
#   if defined O_npzd
#    if defined O_npzd_nitrogen
#    endif
#   endif
#  endif
#  if defined O_sealev_data && defined O_sealev_salinity
#  endif
# endif
#  if defined O_ism
#  endif
# if defined O_embm
# endif
# if defined O_ice && defined O_embm
#  if defined O_ice_cpts
#  else
#  endif
# endif
# if defined O_ism      
# endif
# if defined O_mom
# endif
# if defined O_sed
# endif
# if defined O_co2emit_diag
# endif
# if defined O_global_sums
# endif
# if defined O_co2emit_diag
# endif
# if defined O_global_sums
# endif
# if defined O_global_sums
#  if defined O_embm
#  else
#  endif
# endif
#endif
