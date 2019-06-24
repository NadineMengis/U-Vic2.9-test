! source file: /net/mare/home1/eby/as/ism/gisbc.F
      subroutine gisbc

!=======================================================================
!     Accumulation boundary conditions for ice sheet model
!=======================================================================

      use subgrid

      implicit none

      include "size.h"
      include "ism.h"
      include "param.h"
      include "atm.h"
      include "ice.h"
      include "csbc.h"
      include "tmngr.h"
      include "cembm.h"
      integer i,j,n
      integer ntrec
      data ntrec /0/
      save ntrec
      character(120) g_st

      integer monprev,monnext
      real modtime,wtnext,wtprev

      do j=2,jmtm1
        do i=2,imtm1
          !Get subgrid ice temperature boundary condition
	  do n=strtiE(i,j),endiE(i,j)
            !Ice temperature
            ismticeiE(n) = ismticeiE(n) + ticeiE(n)
	  enddo
	  ismtice(i,j) = ismtice(i,j) + ticel(i,j)
        enddo
      enddo

      !counter
      iac=iac+1

      ntrec=ntrec+1

      return
      end
