csbc.h
#if defined O_embm_awind
#endif
#if defined O_carbon_co2_2d
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
#endif
#if defined O_npzd_nitrogen
# if !defined O_npzd_no_vflux
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
#if defined O_plume
#endif
#if defined O_mtlm
#endif
#if defined O_global_sums || defined O_co2emit_diag
#endif
