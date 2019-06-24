iceshow.f
O/ice_modules.o
netcdf.inc
#if defined EISA
#endif
#if defined EISVAR20 || defined EISVAR40
#endif
#if defined HOMA || defined HOMB || defined HOMC || defined HOMD || defined HOME
#else
#endif
#if defined EISANTA
#endif
#if defined HEINO
#endif
#if defined SCHOOFGL
#endif
#if defined HOMA || defined HOMB || defined HOMC || defined HOMD || defined HOMF
#if defined HOMF
#endif
#endif
#if defined HEINO
#endif
#if defined HOMA || defined HOMB || defined HOMC || defined HOMD
#endif
#if defined HOMA || defined HOMC
#endif
#if defined LONLAT
#else
#endif
#if defined LONLAT
#else
#endif
#if defined EISLINE && defined LINEB
#elif defined EISLINE && defined LINED
#elif defined EISLINE && defined LINEE
#elif defined EISLINE && defined LINEF
#elif defined EISLINE && defined LINEG
#elif defined EISLINE && defined TRANSECTA
#else
#endif
#if defined EISLINE && defined LINEE
#endif
#if defined EISLINE && defined LINEE
#endif
#if defined EISLINE && defined LINEE
#else
#endif
#if defined NETCDF
#include <netcdf.inc>
# if !defined O_ism
# endif     
#else
# if !defined O_ism 
# endif
#endif
#if defined O_ism
#endif
#if defined HEINO
#endif
#if defined NETCDF 
#else
#if defined SGI || (defined SUN && defined F77)
#else
#endif
#endif
#if defined NETCDF 
#if defined O_ism
#endif
#else
#endif
#if defined NETCDF 
#endif
#if defined NETCDF 
#else
#endif
#if defined HEINO
#endif
#if defined HOM
#if defined HOMC 
#endif
#if defined HOMD 
#endif
#endif
