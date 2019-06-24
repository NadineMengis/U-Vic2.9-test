! source file: /net/mare/home1/eby/as/2.9.old3/source/common/levind.h
!====================== include file "levind.h" ========================

!     vertical level indicators which define model geometry & bottom
!     topography:

!     kmt = number of vertical boxes over "t" points
!     kmu = number of vertical boxes over "u,v" points

      integer kmt, kmu
      common /levind/ kmt(imt,jmt), kmu(imt,jmt)
