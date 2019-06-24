! source file: /net/mare/home1/eby/as/2.9.old3/source/common/setcom.F
      subroutine setcom (is, ie, js, je)

!=======================================================================
!     set up everything which must be done only once per run
!=======================================================================

      implicit none

      character (120) :: fname, new_file_name, logfile

      integer ie, is, je, js, i, indp, iormsk, iou, jrow, k, ke, ks, n

      logical exists

      real fx, fxa, fxb, tarea, uarea, ocnp, ej, ek

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "scalar.h"
      include "stdunits.h"
      include "coord.h"
      include "cregin.h"
      include "grdvar.h"
      include "index.h"
      include "iounit.h"
      include "levind.h"
      include "csbc.h"
      include "cpolar.h"

!     set latitudes used in filtering of tracer and velocity fields

      rjfrst = -87.3
      rjft0 = -67.5
      rjft1 = -69.3
      rjft2 = 69.3
      rjfu0 = -68.4
      rjfu1 = -70.2
      rjfu2 = 70.2

!     default average surface salinity
      socn = 0.03475

!     compute sin and cos values for vector correction before filtering

      fx =  1.0e-10
      fxa = dxt(1)/radius

      do i=2,imtm1
        fxb = fxa*float(i-2)
        spsin(i) = sin(fxb)
        spcos(i) = cos(fxb)
        if (abs(spsin(i)) .lt. fx) spsin(i) = c0
        if (abs(spcos(i)) .lt. fx) spcos(i) = c0
      enddo

      spsin(1) = c0
      spcos(1) = c0
      spsin(imt) = c0
      spcos(imt) = c0

!     set up model indices for filtering high latitudes

      jfrst  = indp (rjfrst, yt, jmt)
      jft0   = indp (rjft0, yt, jmt)
      jft1   = indp (rjft1, yt, jmt)
      jft2   = indp (rjft2, yt, jmt)
      jfu0   = indp (rjfu0, yu, jmt)
      jfu1   = indp (rjfu1, yu, jmt)
      jfu2   = indp (rjfu2, yu, jmt)
      jskpt  = jft2-jft1
      jskpu  = jfu2-jfu1
      njtbft = (jft1-jfrst+1)+(jmtm1-jft2+1)
      njtbfu = (jfu1-jfrst+1)+(jmtm1-jfu2+1)
      if (njtbft .gt. jmtfil .or. njtbfu .gt. jmtfil) then
        write (stdout,9599) njtbft, njtbfu
        write (stderr,9599) njtbft, njtbfu
        stop '>ocean 1'
      endif
9551  format (/' ==== start and end indices for Fourier filtering ====')
9552  format (' only 11 sets of indices fit across the page.',
     &       '  others will not be printed.'/)
9553  format (/,' == filtering indices for t,s ==')
9554  format (/,' == filtering indices for u,v ==')
9555  format (/,' == filtering indices for stream function ==')
9599  format (/,' error => jmtfil must be >= max(njtbft,njtbfu)',
     &        /,'          njtbft=',i8,' njtbfu=',i8)

!-----------------------------------------------------------------------
!     compute surface area and volume of ocean ("t" cells and "u" cells)
!     (note that areas are defined at each level)
!-----------------------------------------------------------------------

      do k=1,km
        tcella(k) = c0
        ucella(k) = c0
      enddo
      ocnp   = c0
      tcellv = c0
      ucellv = c0

!     this comment directive turns off autotasking on the YMP
!     for the following loop

      tgarea(:,:) = 0.
      ugarea(:,:) = 0.
      do jrow=2,jmtm1
        do i=2,imtm1
          tarea = cst(jrow)*dxt(i)*dyt(jrow)
          uarea = csu(jrow)*dxu(i)*dyu(jrow)
          tgarea(i,jrow) = tarea
          ugarea(i,jrow) = uarea
          if (kmt(i,jrow) .gt. 0) then
            do k=1,kmt(i,jrow)
              tcella(k) = tcella(k) + tarea
            enddo
            tcellv = tcellv + tarea*zw(kmt(i,jrow))
            ocnp   = ocnp + float(kmt(i,jrow))
          endif
          if (kmu(i,jrow) .gt. 0) then
            do k=1,kmu(i,jrow)
              ucella(k) = ucella(k) + uarea
            enddo
            ucellv = ucellv + uarea*zw(kmu(i,jrow))
          endif
        enddo
      enddo

      write (stdout,9341) tcella(1), tcellv, ucella(1), ucellv
9341  format (//,'  Global ocean statistics:'
     &,/,'  the total ocean surface area (t cells) =',1pe15.8,'cm**2'
     &,/,'  the total ocean volume (t cells)       =',1pe15.8,'cm**3'
     &,/,'  the total ocean surface area (u cells) =',1pe15.8,'cm**2'
     &,/,'  the total ocean volume (u cells)       =',1pe15.8,'cm**3')

!---------------------------------------------------------------------
!     set the horizontal regional masks and names to be used when
!     computing averages on the "t" grid in subroutine "region.F".
!     also set the vertical regional masks and names for use (along
!     with the horizontal ones) in term balance calculations for
!     tracers & momentum. For term balance calculations the number
!     of masks is the product of horizontal & vertical regions.
!---------------------------------------------------------------------

      fname = new_file_name ("region_masks")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
!     read in horizontal & vertical regional masks ("mskhr" & "mskvr")
!     and names ("hregnm" & "vregnm") on unit iormsk

        call getunit (iormsk, trim(fname)
     &,               'formatted sequential rewind')
        call reg1st (iormsk, .true., .false., .false., .false., .true.)
        call relunit (iormsk)
      else

!     set up the horizontal regions:
!     specify "mskhr" and "hregnm" values rather than reading them in
!     from a file (arbitrarily defaulted below for 5 zonal bands).
!     note: there must be "nhreg" calls ... one for each horizontal
!     region. The form is:
!     call sethr (region number, starting lon, ending lon, starting lat,
!                 ending lat)
!     where the lon and lat limits are for the edges of the "t" cells
!     "sethr" will fit the region using the nearest model grid points

        do jrow=1,jmt
          do i=1,imt
            mskhr(i,jrow) = 0
          enddo
        enddo

        js = 1
        ej = 0.0
        do n=1, nhreg
          ej = ej + float(jmtm1)/float(nhreg)
          je = min(jmtm1, nint(ej))
          call sethr (n, xu(1), xu(imtm1), yu(js), yu(je))
          js = je
        enddo

        do jrow=1,jmt
          mskhr(1,jrow)   = 0
          mskhr(imt,jrow) = 0
        enddo

!     set up the vertical regions:
!     specify "mskvr" and "vregnm" values. also the starting & ending
!     levels for the vertical regions (arbitrarily defaulted for two
!     vertical regions). Note: there must be "nvreg" calls... one for
!     each vertical region.
!     The form is:
!     call setvr (region number, starting depth, ending depth)
!     where the depth limits are for the edges of the "t" cells
!     "setvr" will fit the region using the nearest model grid points

        ek = 0.0
        do n=1, nvreg
          ek = ek + float(km)/float(nvreg)
          ke = min(km, nint(ek))
          if (n .eq. 1) then
            call setvr (n, 0.0, zw(ke))
          else
            call setvr (n, zw(ks), zw(ke))
          endif
          ks = ke
        enddo

        do k=1,km
          mskvr(k) = 0
        enddo
        do n=1,nvreg
          ks = llvreg(n,1)
          ke = llvreg(n,2)
          do k=ks,ke
            mskvr(k) = n
          enddo
        enddo

!     print out the regional mask to a ascii file

        call getunit (iou, trim(fname), 'f s r')
        call reg1st (iou, .true., .false., .false., .true., .false.)
        call relunit (iou)
      endif

      return
      end

      subroutine sethr (nr, xstart, xend, ystart, yend)

!=======================================================================
!     discretizes the horizontal region to nearest model grid points

!     nr     = the horizontal region number
!     xstart = starting longitude at edge of "t" box region
!     xend   = ending longitude at edge of  "t" box region
!     ystart = starting latitude at edge of  "t" box region
!     yend   = ending latitude at edge of "t" box region
!=======================================================================

      implicit none

      integer jer, jsr, indp, ier, isr, nr, j, i, ker

      real xsrl, xerl, ysrl, yerl, xstart, xend, ystart, yend

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cregin.h"
      include "levind.h"

!     find the nearest "t" box indices within the region

      jsr = min(indp (ystart, yu, jmt)+1, jmt)
      jer = indp (yend, yu, jmt)

      isr = min(indp (xstart, xu, imt)+1, imt)
      ier = indp (xend, xu, imt)

!     define "edges" of "t" box region

      if (isr .eq. 1) then
        xsrl = xu(1) - dxudeg(1)
      else
        xsrl = xu(isr-1)
      endif
      xerl = xu(ier)
      if (jsr .eq. 1) then
        ysrl = yu(1) - dyudeg(1)
      else
        ysrl = yu(jsr-1)
      endif
      yerl = yu(jer)

      write (hregnm(nr),9000) xsrl, xerl, ysrl, yerl
      write (stdout,*) ' Defining horizontal region # ',nr
     &, ' as "t" cells within ', hregnm(nr)
      if (isr .gt. ier) then
        write (stdout,*) ' Error: isr=',isr,' >  ier=',ier,' in sethr'
        stop '=>sethr'
      endif
      if (jsr .gt. jer) then
        write (stdout,*) ' Error: jsr=',jsr,' >  jer=',jer, 'in sethr'
        stop '=>sethr'
      endif
      do j=jsr,jer
        do i=isr,ier
          if (kmt(i,j) .gt. 0)  then
            mskhr(i,j) = nr
          endif
        enddo
      enddo
9000  format ('lon: ',f5.1,' => ',f5.1,'  lat: ',f5.1,' => ',f5.1)

      return
      end

      subroutine setvr (nr, zstart, zend)

!=======================================================================
!     discretizes the vertical region to nearest model grid points

!     nr     = the vertical region number
!     zstart = starting depth at edge of "t" box region in cm.
!     zend   = ending depth at edge of "t" box region in cm.
!=======================================================================

      implicit none

      integer ker, ksr, nr, indp

      real tmin, tmax, reflat, zstart, ztopb, zend

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cregin.h"

!     find the nearest "t" box indices within the region

      if (zstart .lt. p5*zw(1)) then
        ksr = 1
        ztopb = 0.0
      else
        ksr = min(indp (zstart, zw, km)+1, km)
        ztopb = zw(ksr-1)*0.01
      endif
      ker = indp (zend, zw, km)
      llvreg(nr,1) = ksr
      llvreg(nr,2) = ker

      write (vregnm(nr),9000) ztopb, zw(ker)*0.01
      write (stdout,*) ' Defining vertical region # ',nr
     &,  ' as "t" cells within ',  vregnm(nr)
      if (ksr .gt. ker) then
        write (stdout,*) ' Error: ksr=',ksr,' >  ker=',ker, 'in setvr'
        stop '=>setvr'
      endif
9000  format (' dpt:',f6.1, '=>',f6.1, 'm')

      return
      end
