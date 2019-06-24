icectl.f
O/ice_modules.o
switch.h
#if defined O_ism
# if defined O_ism_uncoupled
# else
# endif
# if defined CO2INTER 
# endif
# if defined O_ism
#  if defined O_ism_pdd
#  endif
# endif
# if defined LINEE && defined RESTARTE
#  include <fort.26x>
# endif
# if defined SUN
# endif
# if !defined O_ism_uncoupled
# endif
# if defined SUN
# endif
# if defined NESTING
# else
#  if !defined O_ism_standalone
#  elseif defined O_ism_standalone
#  endif      
# endif
# if defined NHA || defined CARB
# else
# endif
# if defined CO2INTER
# endif
# if defined EISANTA
# else
# endif
# if defined O_ism_data_geoflux
# else
# endif
# if defined HEINO 
# endif
# if defined LINEE && defined RESTARTE
#  if defined WEDGE4KM
#  elif defined WEDGE8KM
#  elif defined WEDGE12KM
#  elif defined WEDGE16KM
#  endif
#  endif
# if !defined O_ism_uncoupled
# endif
# if defined O_ism_uncoupled
# endif
# if defined O_ism_glacial_transient
# endif      
# if !defined O_ism_uncoupled
# endif
# if !defined O_ism_uncoupled
# if !defined O_ism_standalone
# else 
# endif                  
# endif
# if defined EIS1 || defined EIS2
# elif defined EISA
# elif defined EISLINE && ! defined TRANSECTA
# elif defined EISLINE && defined TRANSECTA
# endif
# if defined O_ism_uncoupled
#  if defined CO2INTER
#  endif
#  if defined GCMMATRIX
#   if defined CO2INTER
#   endif
#  elif defined SUBICE
#  else
#   if defined CO2INTER
#   endif
#  endif
#  if defined CO2INTER
#  endif
# else
#  if defined O_ism_glacial_transient
#  endif
#  if defined fade2glac
# endif  
# endif
# if defined BEDROCK
# endif
# if defined NOMOVEW
# else
# endif
# if defined SEDIMENT
# endif
# if defined NESTING
# endif
# if defined EISLINE
# elif defined EISANTA || defined NHA || defined CARB || defined HEINO
# endif
# if defined HEINO 
# endif
# if defined HOM 
# endif
# if defined O_ism_uncoupled      
# else
# endif
# if !defined O_ism_standalone
# else
# endif
# if !defined O_ism_uncoupled
# endif 
# if defined O_ism_uncoupled
# else
# endif
# if defined SUN || defined CRAY || defined IGUANA
# elif defined LION
# elif defined NCOM
# endif
# if defined EISLINE
#  if defined LINEB
#   if defined SIDEDRAG
#   else
#   endif
#  elif defined LINEC
#  elif defined LINED
#  elif defined LINEE
#   if defined SEALEVCHANGE
#    if defined SEALEVFALL
#    elif defined SEALEVRISE
#    endif
#   elif defined RESTARTE
#   elif defined INITEE
#   else
#   endif
#  elif defined LINEF
#  elif defined LINEG
#  elif defined TRANSECTA
#   if defined SIDEDRAG
#   endif
#  else
#  endif
# elif defined EISANTA
#  if defined FORCEPLEIST
#  elif defined FORCELGM
#  elif defined FORCEMIS31
#  elif defined FORCEEO
#  else
#  endif
# elif defined NHA || defined CARB
# else
# endif
# if defined EISANTA
#  if defined SUN || defined CRAY || defined IGUANA
#  elif defined LION
#  endif
#  if defined SUN || defined CRAY || defined IGUANA
#  elif defined LION
#  endif
#  if defined SUN || defined CRAY || defined IGUANA
#  elif defined LION
#  endif
#  if defined SUN || defined CRAY || defined IGUANA
#  elif defined LION
#  endif
#  if defined FORCEEO
#  else
#  endif
#  if defined SUN || defined CRAY || defined IGUANA
#  elif defined LION
#  elif defined NCOM
#  endif
# endif
#else
#endif
