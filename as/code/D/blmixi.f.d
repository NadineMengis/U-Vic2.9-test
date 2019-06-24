blmixi.f
size.h
param.h
pconst.h
stdunits.h
coord.h
hmixc.h
vmixc.h
#if defined O_mom
# if defined O_bryan_lewis_vertical || defined O_bryan_lewis_horizontal
#  if defined O_bryan_lewis_vertical
#  endif
#  if defined O_bryan_lewis_horizontal
#  else
#  endif
# endif
#endif
