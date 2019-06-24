! source file: /net/mare/home1/eby/as/updates/gosbc.F
      subroutine gosbc (is, ie, js, je)

!=======================================================================
!     calculate the average fluxes for next ocean time step
!=======================================================================

      implicit none

      integer ie, is, je, js, i, j, nc
      integer iem1, isp1, jem1, jsp1, k

      real f1, f1a, f1l, fh, fs, fwcflx, fwaflx, time
      real area, tarea, tsflx, rsocn, tmp, fgs

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "csbc.h"
      include "grdvar.h"
      include "tmngr.h"
      include "switch.h"
      include "cembm.h"
      include "atm.h"
      include "mw.h"
      include "ice.h"
      include "mtlm.h"
      include "levind.h"
      include "sed.h"

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

      f1 = 1./atatm
      fh = 2.389e-8/atatm
      fs = -socn/atatm
!-----------------------------------------------------------------------
!     calculate average net fluxes. convert heat flux to cal/cm**2/s
!     from kW/m**2 and fresh water flux (cm/s) to an apparent salt
!     flux (g/cm**2/s) using global ocean average salinity, socn
!-----------------------------------------------------------------------

      do j=2,jmtm1
        do i=2,imtm1
          if (tmsk(i,j) .gt. 0.0) then
            sbc(i,j,ihflx) = sbc(i,j,ihflx) + fh*fluxo(i,j,isat)
!           add virtual fluxes of salinity
            sbc(i,j,isflx) = sbc(i,j,isflx) + fs*fluxo(i,j,ishum)

          else

            sbc(i,j,ihflx) = 0.
            sbc(i,j,isflx) = 0.
          endif
          if (umsk(i,j) .ge. 0.5) then
            sbc(i,j,itaux) = f1*fluxo(i,j,nat+1)
            sbc(i,j,itauy) = f1*fluxo(i,j,nat+2)
          else
            sbc(i,j,itaux) = 0.
            sbc(i,j,itauy) = 0.
          endif
        enddo
      enddo

      call setbcx (sbc(1,1,ihflx), imt, jmt)
      call setbcx (sbc(1,1,isflx), imt, jmt)
      call setbcx (sbc(1,1,itaux), imt, jmt)
      call setbcx (sbc(1,1,itauy), imt, jmt)

!-----------------------------------------------------------------------
!     update boundary conditions from the land model
!     do this now instead of in gasbc so fields can be written out
!-----------------------------------------------------------------------

      f1l = 0.
      f1a = 0.
      if (atatm .ne. 0.) f1a = 1.0/atatm
      if (atlnd .ne. 0.) f1l = 1.0/atlnd
      do j=2,jmtm1
        do i=2,imtm1
          if (land_map(i,j) .ne. 0) then
            sbc(i,j,iro) = sbc(i,j,iro)*f1l
            sbc(i,j,iscal) = sbc(i,j,iscal)*f1l
            sbc(i,j,ievap) = sbc(i,j,ievap)*f1l
            sbc(i,j,ilwr) = sbc(i,j,ilwr)*f1l
            sbc(i,j,isens) = sbc(i,j,isens)*f1l
            sbc(i,j,inpp) =  sbc(i,j,inpp)*f1l
            sbc(i,j,isr) =  sbc(i,j,isr)*f1l
            sbc(i,j,iburn) =  sbc(i,j,iburn)*f1l
          else
            sbc(i,j,iro) = sbc(i,j,iro)*f1a
            sbc(i,j,ievap) = 0.
            sbc(i,j,ilwr) = 0.
            sbc(i,j,isens) = 0.
            sbc(i,j,inpp) = 0.
            sbc(i,j,isr) = 0.
            sbc(i,j,iburn) = 0.
          endif
        enddo
      enddo
      call setbcx (sbc(1,1,iscal), imt, jmt)
      call setbcx (sbc(1,1,iscao), imt, jmt)
      call setbcx (sbc(1,1,ievap), imt, jmt)
      call setbcx (sbc(1,1,ilwr), imt, jmt)
      call setbcx (sbc(1,1,isens), imt, jmt)
      call setbcx (sbc(1,1,inpp), imt, jmt)
      call setbcx (sbc(1,1,isr), imt, jmt)
      call setbcx (sbc(1,1,iro), imt, jmt)

!-----------------------------------------------------------------------
!     zero diagnostic for river discharge and call river model
!-----------------------------------------------------------------------
      disch(:,:) = 0.
      call rivmodel

!-----------------------------------------------------------------------
!     apply CaCO3 weathering flux proportional to discharge
!-----------------------------------------------------------------------
      if (weathflx .gt. 0) then
        tmp = 0.
        do j=2,jmtm1
          do i=2,imtm1
            if (tmsk(i,j) .gt. 0.00001) then
              tmp = tmp + disch(i,j)*dxt(i)*dyt(j)*cst(j)
            endif
          enddo
        enddo
        tmp = weathflx/tmp
        do j=2,jmtm1
          do i=2,imtm1
            if (tmsk(i,j) .gt. 0.0) then
              sbc(i,j,idicflx) = sbc(i,j,idicflx) + disch(i,j)*tmp
              sbc(i,j,ialkflx) = sbc(i,j,ialkflx) + disch(i,j)*2.*tmp
            endif
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     add normalized virtual fluxes to other tracers
!-----------------------------------------------------------------------
      tarea = 0.
      tsflx = 0.
      rsocn = 1./socn
      do j=2,jmtm1
        do i=2,imtm1
          if (tmsk(i,j) .gt. 0.00001) then
            area = dxt(i)*dyt(j)*cst(j)
            tarea = tarea + area
            tsflx = tsflx + sbc(i,j,isflx)*area
          endif
        enddo
      enddo
      tsflx = tsflx/tarea
      do j=2,jmtm1
        do i=2,imtm1
          if (tmsk(i,j) .gt. 0.00001) then
            tmp = (sbc(i,j,isflx) - tsflx)*rsocn
            vflux(i,j) = tmp
            sbc(i,j,idicflx) = sbc(i,j,idicflx) + gaost(idic)*tmp
            sbc(i,j,ialkflx) = sbc(i,j,ialkflx) + gaost(ialk)*tmp
            sbc(i,j,io2flx) = sbc(i,j,io2flx) + gaost(io2)*tmp
            sbc(i,j,ipo4flx) = sbc(i,j,ipo4flx) + gaost(ipo4)*tmp
            sbc(i,j,iphytflx) = sbc(i,j,iphytflx) + gaost(iphyt)*tmp
            sbc(i,j,izoopflx) = sbc(i,j,izoopflx) + gaost(izoop)*tmp
            sbc(i,j,idetrflx) = sbc(i,j,idetrflx) + gaost(idetr)*tmp
            sbc(i,j,ino3flx) = sbc(i,j,ino3flx) + gaost(ino3)*tmp
            sbc(i,j,idiazflx) = sbc(i,j,idiazflx) + gaost(idiaz)*tmp
          endif
        enddo
      enddo
      call setbcx (sbc(1,1,idicflx), imt, jmt)
      call setbcx (sbc(1,1,ialkflx), imt, jmt)
      call setbcx (sbc(1,1,io2flx), imt, jmt)
      call setbcx (sbc(1,1,ipo4flx), imt, jmt)
      call setbcx (sbc(1,1,iphytflx), imt, jmt)
      call setbcx (sbc(1,1,izoopflx), imt, jmt)
      call setbcx (sbc(1,1,idetrflx), imt, jmt)
      call setbcx (sbc(1,1,ino3flx), imt, jmt)
      call setbcx (sbc(1,1,idiazflx), imt, jmt)

      return
      end

