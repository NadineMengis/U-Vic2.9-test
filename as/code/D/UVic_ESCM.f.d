UVic_ESCM.f
size.h
param.h
pconst.h
stdunits.h
coord.h
csbc.h
iounit.h
emode.h
levind.h
scalar.h
switch.h
tmngr.h
cembm.h
atm.h
mw.h
ism.h
calendar.h
accel.h
cnep.h
cprnts.h
diag.h
fwa.h
hmixc.h
insolation.h
isopyc.h
mtlm.h
npzd.h
sed.h
stab.h
veg.h
vmixc.h
#if !defined O_ism_standalone
#if defined O_embm
#endif
#if defined O_mom
#endif
#if defined O_ism
#endif      
#if defined O_mom
# if defined O_sed
# endif
#endif
#if defined O_embm
# if defined O_ism
# endif
# if defined O_mtlm
# endif
#endif
#if defined O_mtlm
# if defined O_mtlm_segday
# else
# endif
#endif
#if !defined O_embm && defined O_mom
# if !defined O_replacst
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
# if !defined O_replacst
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
#endif
#if defined O_global_sums || defined O_co2emit_diag
#endif
#if defined O_embm
# if !defined O_mom
#  if defined O_mtlm
#  endif
# endif
#else
#endif
#if defined O_mtlm
#endif
# if defined O_ism	
# endif
#if defined O_mom
# if defined O_ism && !defined O_constant_ice	  	  
# endif
# if defined O_embm
# endif
# if defined O_mtlm
# endif
# if defined O_sed
# endif
#endif
#if defined O_global_sums || defined O_co2emit_diag
#endif
#if defined O_global_sums || defined O_co2emit_diag
#endif
#if defined O_ism      
#endif
#if defined O_sealev_salinity && O_mom
#endif
# if defined O_even_fluxes
#endif
#if defined O_embm && !defined O_mom && !defined O_replacst
#endif
#if defined O_redi_diffusion || defined O_gent_mcwilliams
# if !defined O_isopycmix
# endif
#endif
#if defined O_mom
#endif
#if defined O_embm
#endif
#if defined O_ism
#endif
#if defined O_mtlm
#endif
#if defined O_sed
#endif
#if defined O_mom && defined O_embm && defined O_restorst
#endif
#if defined O_mom && defined O_embm && defined O_replacst
#endif
#if defined O_restorst && defined O_replacst
#endif
#if defined O_mom
#endif
#if defined O_embm_awind && defined O_embm
#endif
#if defined O_embm
# if defined O_carbon_co2_2d
# endif
#endif
#if defined O_shortwave
#endif
#if defined O_ice_evp
#endif
#if defined O_carbon
# if defined O_carbon_14
# endif
#endif
#if defined O_npzd_alk
#endif
#if defined O_npzd_o2
#endif
#if defined O_npzd
# if !defined O_npzd_no_vflux
# endif
# if defined O_npzd_nitrogen
#  if !defined O_npzd_no_vflux
#  endif
# endif
#endif
#if defined O_cfcs_data_transient
#endif
#if defined O_mtlm
#endif
#if defined O_mtlm && defined O_carbon
#endif
#if defined O_sed
# if defined O_carbon
# endif
# if defined O_npzd_alk
# endif
#endif
#if defined O_use_bias_SAT_correction
#endif
#if defined O_data_sat
#endif
#if defined O_mom
# if defined O_carbon
#  if defined O_carbon_14
#  endif
# endif
# if defined O_cfcs_data_transient
# endif
# if defined O_npzd_o2
# endif
# if defined O_npzd_alk
# endif
# if defined O_npzd
#  if defined O_npzd_nitrogen
#  endif
# endif
# if defined O_carbon && defined O_carbon_14
# endif
# if defined O_npzd
#  if defined O_carbon
#  endif
#  if defined O_npzd_o2
#  endif
#  if defined O_npzd_alk
#  endif
#  if defined O_npzd_nitrogen
#  endif
# endif
#endif
#if defined O_embm
# if defined O_carbon && defined O_carbon_co2_2d
# endif
#endif
#if !defined O_mom || !defined O_time_step_monitor
#endif             
#if !defined O_mom || !defined O_gent_mcwilliams
#endif       
#if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion
#else
#endif
# if defined O_bryan_lewis_vertical
# endif
# if defined O_bryan_lewis_vertical
# endif
# if defined O_implicitvmix
# else
# endif
# if defined O_embm_snow_transient
# endif
# if defined O_embm
# else
# endif
# if defined O_mom
#  if defined O_isopycmix
#  endif
# else
# endif
# if defined O_carbon_14_coupled
# endif
#endif
