! source file: /net/mare/home1/eby/as/updates/rivmodel.F
      subroutine rivmodel

!=======================================================================
!     river model for energy-moisture balance model
!     calculates river runoff from precipitation over land
!=======================================================================

      implicit none

      integer i, j, n, m

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "riv.h"
      include "atm.h"
      include "cembm.h"
      include "grdvar.h"
      include "ism.h"
!-----------------------------------------------------------------------
!     sum (area integral) the precip over each basin
!-----------------------------------------------------------------------
      psum(:) = 0.0
      do j=2,jmtm1
        do i=2,imtm1
          n = nriv(i,j)
          if (n .ne. 0) then
            psum(n) = sbc(i,j,iro)*dxt(i)*cst(j)*dyt(j) + psum(n)
          endif
        enddo
      enddo
      hsum(:) = 0.0

      do n=1,nb
        if (isb(n) .gt. 0) then
          !add ice melt to basin runoff flux
          psum(n) = psum(n) + isflxm(n)
          !remove flux from ice sheet basin runoff 'bank'
          dicevol(n)=dicevol(n)-isflxm(n)/0.910!cm3 ice/s
     &                                   /1.e6 !m3 ice/s
     &                                   *atatm!m3 ice
          !accumulate total change in ice volume
          disvol = disvol - isflxm(n)/0.910!cm3 ice/s
     &                                   /1.e6 !m3 ice/s
     &                                   *atatm

          !add heat to basin heat runoff flux
          hsum(n) = hsum(n) + isflxh(n)
          !remove heat from ice sheet heat 'bank'
          diceheat(n)=diceheat(n)-isflxh(n)*atatm
        endif
      enddo
!-----------------------------------------------------------------------
!     calculate discharge for each discharge point in each basin
!     and add it to the fresh water flux
!-----------------------------------------------------------------------

      do j=2,jmtm1
        do i=2,imtm1
          !set n to discharge number of cell
          n = ndis(i,j)
          if (n .ne. 0) then
            !set discharge to fraction of the total value for that basin
            !discharge fraction based on wdar
            disch(i,j) = psum(n)*wdar(i,j)
            sbc(i,j,isflx) = sbc(i,j,isflx) - socn*psum(n)*wdar(i,j)
            sbc(i,j,ihflx) = sbc(i,j,ihflx) + 0.2389*hsum(n)*wdar(i,j)
          endif
        enddo
      enddo
      call setbcx (sbc(1,1,isflx), imt, jmt)

      return
      end

      subroutine rivinit

!=======================================================================
!     read river model basin and discharge points
!=======================================================================

      implicit none

      character(120) :: fname, new_file_name

      logical exists

      integer i, ios, iou, j, k, m, n, ib(10), ic(10)

      real wt

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "atm.h"
      include "riv.h"
      include "grdvar.h"

      character(1) :: line(imt)

      real tmpij(imtm2,jmtm2)

      wdar(:,:) = 0.
      nriv(:,:) = 0
      ndis(:,:) = 0
      where (tmsk(:,:) .lt. 1.) nriv(:,:) = 1
!     read the netcdf river file.
      fname = new_file_name ("rivers.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Warning => ", trim(fname), " does not exist."
      else
        ib(:) = 1
        ic(:) = imtm2
        ic(2) = jmtm2
        call openfile (fname, iou)
        call getvara ('rivers', iou, imtm2*jmtm2, ib, ic, tmpij
     &, c1, c0)
        wdar(2:imtm1,2:jmtm1) = tmpij(1:imtm2,1:jmtm2)
        call embmbc (wdar)
        nriv(:,:) = wdar(:,:)
        call getvara ('discharge', iou, imtm2*jmtm2, ib, ic, tmpij
     &, c1, c0)
        wdar(2:imtm1,2:jmtm1) = tmpij(1:imtm2,1:jmtm2)
        call embmbc (wdar)
        ndis(:,:) = wdar(:,:)
        call getvara ('weights', iou, imtm2*jmtm2, ib, ic, tmpij
     &, c1, c0)
        wdar(2:imtm1,2:jmtm1) = tmpij(1:imtm2,1:jmtm2)
        call embmbc (wdar)
      endif
      call expand_basins
      call init_discharge

      return
      end

      subroutine expand_basins

!=======================================================================
!     expand basins out over the ocean in case of land boundary charges
!     or areas without discharge points
!=======================================================================

      implicit none

      integer i, j, n, nmax

      logical done

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "riv.h"

      done = .false.
      n = 0
      nmax = max (imt,jmt) + 1
      do while (.not. done .and. n .le. nmax)
        nrfill(:,:) = nriv(:,:)
        done = .true.
        do j=2,jmtm1
          do i=2,imtm1
            if (nrfill(i,j) .eq. 0) then
              if (nriv(i-1,j+1) .ne. 0) then
                nrfill(i,j) = nriv(i-1,j+1)
              elseif (nriv(i,j+1) .ne. 0) then
                nrfill(i,j) = nriv(i,j+1)
              elseif (nriv(i+1,j+1) .ne. 0) then
                nrfill(i,j) = nriv(i+1,j+1)
              elseif (nriv(i-1,j) .ne. 0) then
                nrfill(i,j) = nriv(i-1,j)
              elseif (nriv(i+1,j) .ne. 0) then
                nrfill(i,j) = nriv(i+1,j)
              elseif (nriv(i-1,j-1) .ne. 0) then
                nrfill(i,j) = nriv(i-1,j-1)
              elseif (nriv(i,j-1) .ne. 0) then
                nrfill(i,j) = nriv(i,j-1)
              elseif (nriv(i+1,j-1) .ne. 0) then
                nrfill(i,j) = nriv(i+1,j-1)
              endif
              if (nrfill(i,j) .eq. 0) done = .false.
            endif
          enddo
          nrfill(1,j) = nrfill(imtm1,j)
          nrfill(imt,j) = nrfill(2,j)
        enddo
        nriv(:,:) = nrfill(:,:)
        n = n + 1
      enddo
      if (n .ge. nmax) stop 'rivinit nmax'

      return
      end

      subroutine init_discharge

!=======================================================================
!     sets river model basin and discharge points
!=======================================================================

      implicit none

      character(120) :: fname, new_file_name

      integer i, ios, iou, j, k, m, n

      logical done

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "atm.h"
      include "riv.h"
      include "grdvar.h"

      character(1) :: line(imt)
      real wp(maxnb), wt(maxnb)

      nriv(:,:) = nrfill(:,:)

      call basin_perimeter

!     check each basin has discharge
      done = .true.
      do n=1,nb
        if (nbp(n) .ne. 0 .and. ndp(n) .eq. 0) done = .false.
      enddo

!     if no discharge points expand other basins to fill this one
      if (.not. done) then
        do j=2,jmtm1
          do i=2,imtm1
            if (nriv(i,j) .ne. 0) then
              n = nriv(i,j)
              if (nbp(n) .ne. 0 .and. ndp(n) .eq. 0) nriv(i,j) = 0
            endif
          enddo
        enddo
        call expand_basins
        call basin_perimeter
      endif

!     find total specified weight for discharge points for each basin
      wt(:) = 0.
      do j=2,jmtm1
        do i=2,imtm1
          if (wdar(i,j) .ne. 0.0 .and. ndis(i,j) .ne. 0.) then
            wt(ndis(i,j)) = wt(ndis(i,j)) + wdar(i,j)
          endif
        enddo
      enddo

!     calculate discharge weight to perimeter points
      do n=1,nb
        wp(n) = 0.
        if (nbp(n) .ne. 0) then
          if (wt(n) .lt. 1. .and. ndp(n) .gt. 0) then
            wp(n) = (1. - wt(n))/ndp(n)
            wt(n) = 1.
          endif
        endif

      enddo

!     calculate normalized weight over area for all discharge points
      do j=2,jmtm1
        do i=2,imtm1
          n = ndis(i,j)
          if (n .ne. 0) then
            if (nbp(n) .ne. 0) then
              if (wt(n)*dxt(i)*cst(j)*dyt(j) .eq. 0.) stop 'rivmodel wt'
	      wdar(i,j) = (wdar(i,j)+wp(n))/(wt(n)*dxt(i)*cst(j)*dyt(j))
            endif
          endif
        enddo
      enddo

      return
      end

      subroutine basin_perimeter

!=======================================================================
!     sets basin perimeters
!=======================================================================

      implicit none

      integer i, j, n, nmax

      logical done

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "atm.h"
      include "riv.h"

!     find perimeter points
      nb = 0
      nbp(:) = 0
      ndp(:) = 0
      do j=2,jmtm1
        do i=2,imtm1
          if (nriv(i,j) .gt. nb) nb = nriv(i,j)
          if (nb .gt. maxnb) stop "basin_perimeter maxnb"
          if (tmsk(i,j) .eq. 1.0) then ! if full ocean
            if (ndis(i,j) .ne. 0) then ! if discharge point (set via netcdf)
              ndp(ndis(i,j)) = ndp(ndis(i,j)) + 1 !increase # of discharge points... what if this is in middle of ocean?
            elseif (tmsk(i-1,j+1) .lt. 1. .or. tmsk(i,j+1) .lt. 1.
     &        .or. tmsk(i+1,j+1) .lt. 1. .or. tmsk(i-1,j) .lt. 1.
     &        .or. tmsk(i+1,j) .lt. 1. .or. tmsk(i-1,j-1) .lt. 1.
     &        .or. tmsk(i,j-1) .lt. 1. .or. tmsk(i+1,j-1) .lt. 1.)
     &        then
              ndis(i,j) = nriv(i,j)
              ndp(ndis(i,j)) = ndp(ndis(i,j)) + 1
            endif
            nriv(i,j) = 0
          else
            if (nriv(i,j) .eq. 0) stop "basin_perimeter zero nriv"
            nbp(nriv(i,j)) = nbp(nriv(i,j)) + 1
            ndis(i,j) = 0
          endif
        enddo
      enddo
      nriv(1,:) = 0
      nriv(imt,:) = 0
      nriv(:,1) = 0
      nriv(:,jmt) = 0
      ndis(1,:) = 0
      ndis(imt,:) = 0
      ndis(:,1) = 0
      ndis(:,jmt) = 0

      return
      end
