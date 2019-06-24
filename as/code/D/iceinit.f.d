iceinit.f
O/ice_modules.o
O/sgmodule.o
pconst.h
switch.h
#if defined LONLAT
#elif defined STEREO
#if defined EISLINE && defined TRANSECTA
#else
#if defined NESTING
#else
#endif
#if defined HEINO || defined HOM || defined EISLINE
#endif
#if defined HEINO || defined HOM || defined EISLINE
#endif
#if defined HEINO || defined HOM || defined EISLINE
#endif
#if defined HEINO || defined HOM || defined EISLINE
#endif
#if defined EISLINE 
#else
#endif
#if defined NESTING
#endif
#endif
#else
#endif
#if defined HEINO || defined HOM
#else
#endif
#if defined EISLINE      
#if defined LINEB 
#elif defined LINEC
#elif defined LINED
#elif defined LINEE 
#if defined INITEE
#else
#endif
#elif defined LINEF 
#elif defined LINEG
#elif defined TRANSECTA && defined IDEAL
#endif
#endif
#if defined BEDMAP
#elif defined EISMINT2
#elif defined NHA || defined CARB
#endif
#if defined O_ism      
#endif      
#if defined HEINO 
#elif defined HOM
#elif (EISLINE) && ! ( defined LINEB || defined LINEF || defined LINEG || defined TRANSECTA )
#elif (EISLINE) &&   ( defined LINEB || defined LINEF || defined LINEG || defined TRANSECTA )
#else
#endif
#if defined EIS1 || defined EIS2 
#elif defined EIS3 || defined EIS4 
#elif defined EISA || defined EISB || defined EISC || defined EISD
#elif defined EISLINE &&  ( defined LINEB || defined LINEC || defined LINED || defined LINEE || defined LINEG || (defined TRANSECTA && defined IDEAL) )
#elif defined EISLINE &&  defined LINEF
#elif defined BEDMAP
#if (defined EISLINE && defined TRANSECTA) || (defined EISANTA && defined NESTING) 
#elif defined EISANTA && ! defined NESTING
#endif
#if defined icefree_init
#endif
#if defined icefree_init
#endif
#if defined icefree_init
#endif    
#elif defined EISMINT2
#elif defined NHA
#elif defined CARB
#elif defined HEINO
#elif defined HOMA || defined HOMB || defined HOMC||defined HOMD
#if defined HOMA || defined HOMB
#elif defined HOMC || defined HOMD
#endif
#if defined HOMA
#elif defined HOMB
#endif
#elif defined HOME
#elif defined HOMF
#else
#endif
#if defined EISANTA
#endif
#if defined NESTING
#endif
#if defined O_ism      
#endif      
#if ! defined NESTING
# if !defined O_ism
# elif defined O_ism	   
# endif	  
#endif
#if defined HEINO
#if defined HEINOSMALL
#else
#endif
#endif
#if defined BEDMAP &&  defined EISLINE && defined TRANSECTA
#endif
#if defined BEDMAP &&  ((defined EISLINE && defined TRANSECTA) || (defined EISANTA && defined NESTING))
#endif
#if defined BEDMAP &&  (defined EISANTA && ! defined NESTING)
#endif
#if defined BEDMAP
#if defined SUN || defined CRAY || defined IGUANA
#elif defined LION 
#elif defined CONDOR 
#elif defined SGI
#elif defined NCOM
#elif defined O_ism
#endif
#endif
#if defined EISMINT2
#if defined SUN || defined CRAY || defined IGUANA
#elif defined LION 
#elif defined CONDOR 
#elif defined SGI 
#endif
#endif
#if defined NHA
#if defined SUN || defined CRAY  || defined IGUANA
#elif defined LION 
#elif defined CONDOR 
#elif defined SGI 
#endif
#endif
#if defined CARB
#if defined SUN || defined CRAY  || defined IGUANA
#elif defined LION 
#elif defined CONDOR 
#elif defined SGI 
#endif
#endif
#if !defined O_ism_standalone
#endif
#if defined O_isnan
#endif      
#if defined O_isnan      
#endif      
#if defined O_isnan       
#endif      
#if defined O_isnan       
#endif      
#if defined O_isnan       
#endif      
#if defined O_isnan       
#endif      
#if ! defined NESTING
#endif
#if defined NESTING
#endif
