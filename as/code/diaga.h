! source file: /net/mare/home1/eby/as/2.9.old3/source/mom/diaga.h
!====================== include file "diaga.h" =========================

!     variables used for computing diagnostics:

!     totalk   = total number of levels involved in convection
!     vdepth   = ventilation depth (cm)
!     pe       = potential energy lost due to explicit convection (g/s2)

      real totalk, vdepth, pe
      common /cdiaga_r/ totalk(imt,jsmw:jemw), vdepth(imt,jsmw:jemw)
      common /cdiaga_r/ pe(imt,jsmw:jemw)
