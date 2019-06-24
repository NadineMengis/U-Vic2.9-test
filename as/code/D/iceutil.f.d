iceutil.f
#if defined SGI || defined SUN
#if defined SUN && defined F90
#else
#endif
#elif defined AIX
#endif
#if ! defined CRAY
#if defined AIX
#elif defined SUN
#else
#endif
#endif
#if defined SUN
#endif
#if defined SUN || defined CRAY || (IGUANA)
#elif defined CONDOR
#elif defined SGI
#endif
