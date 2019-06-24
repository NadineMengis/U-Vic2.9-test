hmixc.f
size.h
param.h
pconst.h
stdunits.h
grdvar.h
hmixc.h
mw.h
scalar.h
switch.h
#if defined O_mom
# if defined O_full_tensor
# endif
# if defined O_consthmix
#  if defined O_biharmonic
#  else
#  endif
# endif
# if defined O_consthmix
#  if defined O_bryan_lewis_horizontal
#  else
#   if defined O_biharmonic
#   else
#   endif
#   if defined O_full_tensor
#   endif
#  endif
# else
#  if defined O_smagnlmix && !defined O_consthmix
#  endif
# endif
#endif
