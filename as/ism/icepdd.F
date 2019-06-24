c***********************
#if defined O_ism && defined O_ism_pdd
c***********************

c-----------------------------------------------------------------------

      subroutine climate_uvic (budgsnow, budgrain, budgevap, budgmelt,
     *                         tsurf, tsurfi,
     *                         h, hs, maskwater,
     *                         tempmon, precmon, topogmon,
     *                         timeice)

c     Calls snow/ice "model" (surfice_pdd) years to get annual ice mass 
c     balance forcing, using monthly means T and P supplied from
c     climate model. All on ice model grid.
c     Returned:
c       budgsnow (snowfall,     annual mean, m/y)
c       budgrain (rainfall,     annual mean, m/y)
c       budgevap (evaporation,  annual mean, m/y)
c       budgmelt (surface melt, annual mean, m/y)
c       tsurf    (surface air temp, annual mean, deg K)
c       tsurfi   (surface ice temp, annual mean, deg K)
c     Supplied:
c       tempmon (monthly mean surface air temp, deg C)
c       precmon (monthly mean precip, mm/day)
c       topogmon (elevations above sea level of climate model, m)
c       h (ice thickness, m)
c       hs (surface elevation, m)
c       maskwater (0=land, 1=ocean, 2+=inland lake)
c       timeice (current time of run, yrs, -ve=BP)
c       ifrest (-1=cold start, 0=initial (re)start, 1=true restart)

      use comicephys
      use comicegrid

      dimension 
     *  budgsnow(nx,ny),     budgrain(nx,ny),
     *  budgevap(nx,ny),     budgmelt(nx,ny),
     *  tsurf(nx,ny),        tsurfi(nx,ny),
     *  h(nx,ny),            hs(nx,ny),           maskwater(nx,ny),
     *  tempmon(nx,ny,nmon), precmon(nx,ny,nmon), topogmon(nx,ny)

      
      logical first, firstprint 
      data first, firstprint /.true., .true./
      save first, firstprint

c     local 2-D ice-grid arrays:
      dimension
     *  tag(nx,ny),      prg(nx,ny), 
     *  pmemr(nx,ny),    arain(nx,ny),   asnow(nx,ny),
     *  amelt(nx,ny),    ablat(nx,ny),   apdd(nx,ny),
     *  hsnow(nx,ny),    iftodo(nx,ny),
     *  indl(nx,ny),     weil(nx,ny),
     *  indb(nx,ny),     weib(nx,ny),
     *  hsprevcalc(nx,ny)

c     save since may only be set for a subset of points 
c     (whose elev changed enough since last call)
      save pmemr, arain, asnow, amelt, ablat, apdd, hsnow

      save indl, weil, indb, weib, hsprevcalc, timeiceprevall
      data timeiceprevall /0./       ! just to avoid indef below

c     local 1-D ice-strip arrays:
      dimension 
     *  toa(nvec), 
     *  pra(nvec),    raa(nvec),    sna(nvec),
     *  ta(nvec),
     *  tice(nvec),   tsno(nvec),
     *  dmelta(nvec), pdd(nvec),
c       annual accumulators:
     *  araina(nvec), asnowa(nvec), amelta(nvec),
     *  taira(nvec),  ticea(nvec),  pdda(nvec),
c       annual quantities (and hsnowa multiyear):
     *  ablata(nvec), pmemra(nvec), hsnowa(nvec),
     *  iland(nvec,2)
      print*, 'Error: save statements in climate_uvic not'
      stop 'multi-ice-sheet compatible yet.'

      call zero (budgsnow, nx*ny)
      call zero (budgrain, nx*ny)
      call zero (budgevap, nx*ny)
      call zero (budgmelt, nx*ny)
      
      dtimed = 5.*86400.   ! for degree day method, several-day timestep
      dtforc = 20.*86400. 

      nyear  = 2     ! more than 1 year to allow for spinup due to hsnow
 
      nspy   = max (1, nint((86400.*365  )/dtimed))
      nspf   = max (1, nint((dtforc      )/dtimed))

c        Set points to be calculated this call (iftodo, npoi, iland)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (first .or. nint(abs(timeice-timeiceprevall)).ge.1000) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c      if first call this run or every 1000 years, do entire domain:
        npoi = 0
        do j=1,ny
          do i=1,nx
            iftodo(i,j) = 1
            npoi = npoi + 1
            iland(npoi,1) = i
            iland(npoi,2) = j
          enddo
        enddo
        call scopy (nx*ny, hs, 1, hsprevcalc, 1)
        timeiceprevall = timeice
        first = .false.

c~~~~~~~~~
      else
c~~~~~~~~~

c       do points whose sfc elev has changed > 50 m since last calc:
        npoi = 0
        do j=1,ny
          do i=1,nx
            if (abs(hs(i,j)-hsprevcalc(i,j)).gt.50.) then
              iftodo(i,j) = 1
              npoi = npoi + 1
              iland(npoi,1) = i
              iland(npoi,2) = j
              hsprevcalc(i,j) = hs(i,j)
            else
              iftodo(i,j) = 0
            endif
          enddo
        enddo
        if (npoi.eq.0) return

c~~~~~~~~~~
      endif
c~~~~~~~~~~

c        Zero initial snow thickness (carried over from year to year)

      do ip=1,npoi
        hsnowa(ip) = 0.
      enddo

c        Loop over years

c=====================
      do iyear=1,nyear
c=====================

c          Zero annual accumulators

        do ip=1,npoi
          araina(ip) = 0.
          asnowa(ip) = 0.
          amelta(ip) = 0.
          taira (ip) = 0.
          ticea (ip) = 0.
          pdda  (ip) = 0.
        enddo

c          Loop over timesteps through one year

c======================
        do nstep=1,nspy
c======================

          isecyr = nint((nstep-0.5)*dtimed)           ! for degree-days
          isecdy = mod(isecyr,86400)
          calday = 1. + isecyr/86400.

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          if (nstep.eq.1 .or. mod(nstep,nspf).eq.0) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c              Interpolate current surface forcing from monthly mean
c              to current time

            call intdrive (calday, tag, tempmon, nx, ny, nmon)
            call intdrive (calday, prg, precmon, nx, ny, nmon)

c              Change units to MKS

            do j=1,ny
              do i=1,nx 
                prg(i,j) = prg(i,j)/86400.      ! mm/day to kg/m2/s
                tag(i,j) = tag(i,j) + tmelt     ! deg C to K
              enddo
            enddo

c              Transfer meteorological forcing to current 1-D strip

c-----------------------
            do ip=1,npoi
c-----------------------

              i = iland(ip,1)
              j = iland(ip,2)

              toa(ip) = topogmon(i,j)
              pra(ip) = prg(i,j)
              ta(ip)  = tag(i,j)

c                Apply elevation corrections to sfc air temperature ta,
c                and set rainfall vs snowfall

              zdelev = hs(i,j) - toa(ip)

              zt1 = ta(ip)
              zt2 = ta(ip) - rlapse*zdelev
              ta(ip) = zt2

              zpr = pra(ip)
              if (ta(ip).gt.tmelt) then
                raa(ip) = zpr
                sna(ip) = 0.
              else
                raa(ip) = 0.
                sna(ip) = zpr 
              endif

c----------------
            enddo
c----------------

c>>>>>>>>>>>>>>
          endif
c>>>>>>>>>>>>>>

          call surfice_pdd
     *      ( dtimed,
     *        ta,      raa,    sna,
     *        dmelta,  hsnowa, 
     *        tice,    tsno,   pdd, 
     *        npoi)

c            Calculate diagnostic min/max surface temperatures
c            and test for blowup

          ztmin =  1.e20
          ztmax = -1.e20
          do ip=1,npoi
            if (tice(ip).lt.ztmin .or. tsno(ip).lt.ztmin) then
              ztmin = min (tice(ip), tsno(ip))
              ipmin = ip
            endif
            if (tice(ip).gt.ztmax .or. tsno(ip).gt.ztmax) then
              ztmax = max (tice(ip), tsno(ip))
              ipmax = ip
            endif
          enddo

          if (ztmin.lt.tmelt-120.) then
             imin = iland(ipmin,1)
             jmin = iland(ipmin,2)
             write (ioterm,273) iyear, nstep, isecyr/86400.,
     *         'tmin',ztmin-tmelt, 
     *         imin, jmin, alond(imin,jmin), alatd(imin,jmin)
  273        format('*** Warning: Extreme surface temperature(s)',
     *              ' at iyear=',i3,'   nstep=',i6,'   day=',f8.2
     *             /'  ',a,'=',f10.2,'  at i,j=',2i4,'  lon,lat=',2f8.2)
c            call flush(ioterm)    ! doesn't work on linux/iguana
          endif
          if (ztmax.gt.tmelt+90.) then
             imax = iland(ipmax,1)
             jmax = iland(ipmax,2)
             write (ioterm,273) iyear, nstep, isecyr/86400.,
     *         'max', ztmax-tmelt, 
     *         imax, jmax, alond(imax,jmax), alatd(imax,jmax)
c            call flush(ioterm)    ! doesn't work on linux/iguana
          endif

c            Increment annual snowfall, rainfall, melt (kg/m2/s)

          do ip=1,npoi
            araina(ip) = araina(ip) + dtimed*raa(ip)
            asnowa(ip) = asnowa(ip) + dtimed*sna(ip)
            amelta(ip) = amelta(ip) + dtimed*dmelta(ip)
            taira(ip)  = taira(ip)  + ta(ip)/nspy
            ticea(ip)  = ticea(ip)  + tice(ip)/nspy
            pdda(ip)   = pdda(ip)   + pdd(ip)
          enddo

c============
        enddo ! timestep loop
c============

c==========
      enddo ! year loop
c==========

c        Calculate annual budgets for last year, with correction for
c        refreezing. Assumes rain passes through snow with no
c        interaction.

      do ip=1,npoi

c          (1) No runoff if melt/snow < 0.95
c          (2) Fraction runs off if 0.95 < melt/snow < 1.05
c          (3) All runs off if 1.05 < melt/snow, or h=0 (no icesheet)

c            For degree-days model

c       ratsno = 2.7
c       ratice = 7.2
        ratsno = 3.5
        ratice = 7.0

        zpdd = pdda(ip)

        rsw = 0.25/1.0                                      ! rhos/rhow

c       melt 1.-rsw fraction of snow that will saturate snowpack:
        zsuper = min (asnowa(ip)*(1.-rsw), ratsno*zpdd)
        zpdd = zpdd - zsuper/ratsno
        zmelt = zsuper + ratice*zpdd

c       melt embedded snow in the saturated snowpack:
        zrunsno = min (asnowa(ip), (ratice/rsw)*zpdd)
        zpdd = zpdd - zrunsno/(ratice/rsw)

c       melt ice below snow:
        zrunice = ratice*zpdd

        ablata(ip) = zrunsno + zrunice                     ! passed back
        amelta(ip) = zmelt                                 ! diagnostic
        pmemra(ip) = asnowa(ip) - ablata(ip)               ! diagnostic

      enddo

c        Save budget quantities in full ice-grid arrays
c        (so can be displayed by printmap...otherwise could
c        do final calc of budg[snow,rain,ablat] below over 1:npoi,
c        and avoid 2-D arrays for arain, etc)

      do ip = 1,npoi
        i = iland(ip,1)
        j = iland(ip,2)
        arain(i,j) = araina(ip)
        asnow(i,j) = asnowa(ip)
        amelt(i,j) = amelta(ip)
        ablat(i,j) = ablata(ip)
        pmemr(i,j) = pmemra(ip)
        apdd(i,j)  = pdda(ip)
        hsnow(i,j) = hsnowa(ip)
c       returned to icectl:
        tsurf(i,j) = taira(ip)
        tsurfi(i,j)= ticea(ip)
      enddo

c        Transfer to budg* for ice model, 
c        converting flux units from kg/m2/yr to m/yr 
          
      do j=1,ny
        do i=1,nx
          budgsnow(i,j) = asnow(i,j)*.001             ! kg/m2/yr to m/yr
          budgrain(i,j) = arain(i,j)*.001             ! kg/m2/yr to m/yr
          budgmelt(i,j) = ablat(i,j)*.001             ! kg/m2/yr to m/yr
          budgevap(i,j) = 0.

c         ensure zero or slightly negative budget over open ocean:
          if (maskwater(i,j).eq.1) then
            zf = max(0., min (1.,  (h(i,j)-2.5)/2.5))              ! 777
            budgsnow(i,j) = zf * budgsnow(i,j)  + (1.-zf)*(0   )   ! 777
            budgrain(i,j) = zf * budgrain(i,j)  + (1.-zf)*(0.  )   ! 777
            budgevap(i,j) = zf * budgevap (i,j) + (1.-zf)*(0.  )   ! 777
            budgmelt(i,j) = zf * budgmelt(i,j)  + (1.-zf)*(0.02)   ! 777
          endif
        enddo
      enddo

c        Print out maps, if not a restart. 
c        Ok first call since whole domain has been calculated

#if ! defined O_ism
      if (firstprint .and. ifrest.ne.1) then
        write (iuout2d, '(/"*** METEO FIELDS (FIRST CALL) ***"/)') 
        call printmap (timeice, pmemr,'pmemr, mm/yr', 50., h, iuout2d,0)
        call printmap (timeice, arain,'arain, mm/yr', 50., h, iuout2d,0)
        call printmap (timeice, asnow,'asnow, mm/yr', 50., h, iuout2d,0)
        call printmap (timeice, amelt,'amelt, mm/yr', 50., h, iuout2d,0)
        call printmap (timeice, ablat,'ablat, mm/yr', 50., h, iuout2d,0)
        call printmap (timeice, apdd, 'apdd,  C*day', 50., h, iuout2d,0)
        call printmap (timeice, hsnow,'hsnow, m',      .1, h, iuout2d,0)
        write (iuout2d,'(/"*** END OF METEO FIELDS (FIRST CALL) ***"/)')
        call flush(iuout)
        firstprint = .false.
      endif
#endif

      return
      end

c----------------------------------------------------------------------

      subroutine intdrive (calday, field, fieldmon, nx, ny, nmon)

c    Lineary interpolates fieldmon(..,..,nmon) in time between the
c    midmonth points to current calday, and sets current array field.

c calday    = day number within current year (1. to 366.) (supplied)
c field     = 2D agcm current-time array (returned)
c fieldmon  = 2D agcm field for all mid months (supplied)

      dimension field(nx,ny), fieldmon(nx,ny,nmon)

      dimension iday(0:13)
      data iday / -15.,  16.,  45.,  75., 105., 136., 166.,
     *                  197., 228., 258., 289., 319., 350.,381. /

c        Find interpolation indices for iday  (nb: calday runs from
c        1. to 366. from 00:00 GMT Jan 1st to 24:00 GMT Dec 31st)

      do im=1,nmon+1
        if (calday.le.iday(im)) then
          im1 = im-1
          im2 = im
          wm1 = (iday(im2)-calday) / (iday(im2)-iday(im1))
          if (im1.eq.0) im1 = nmon
          if (im2.eq.nmon+1) im2 = 1
          goto 102
        endif
      enddo
  102 continue

c        Interpolate to current time

      do j=1,ny
        do i=1,nx
          field(i,j) = (      wm1 *fieldmon(i,j,im1)
     *                  + (1.-wm1)*fieldmon(i,j,im2) )
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine surfice_pdd
     *      ( dt,
     *        ta,      raa,     sna,
     *        dmelta,  hsnowa,  
     *        tice,    tsno,    pdd,
     *        npoi)
 
c        Simple one-level slab "snow/ice" model 
 
      use comicephys
      use comicegrid
 
      dimension
     *  ta(nvec),      raa(nvec),    sna(nvec),
     *  dmelta(nvec),  hsnowa(nvec), 
     *  tice(nvec),    tsno(nvec),   pdd(nvec)

      dimension qdd(nvec)

      parameter (rhosno = 250.)     ! snow density, kg/m3

      logical surf3, surf4
      parameter (surf3=.true.)      ! degree day 
      parameter (surf4=.false.)     ! degree day,  parameterized diurnal
c     parameter (surf3=.false.)     ! degree day 
c     parameter (surf4=.true.)      ! degree day,  parameterized diurnal
 
c=================
      do ip=1,npoi
c=================

        zs = max (0., min (1., hsnowa(ip)/.15))  !fractional snow area

c----------------------
        if (surf3) then
c----------------------
c         pdd(ip) = max (ta(ip)-(tmelt+5.), 0.) * (dt/86400.)
c         zfac    = 3.0*(zs*1.0  + (1.-zs)*1.5)

          pdd(ip) = max (ta(ip)-tmelt, 0.) * (dt/86400.)
          zfac    = zs*3.5  + (1.-zs)*8.16                 !Ritz,CD,1997

          dmelta(ip) = zfac*pdd(ip) / dt

c---------------------------
        else if (surf4) then
c---------------------------
c         use mean ta for degree-day:
          call degree_day (dt, ta(ip), pdd(ip), qdd(ip), tmelt, 1)

          zpr = raa(ip) + sna(ip)
          raa(ip) = (1.-qdd(ip)) * zpr
          sna(ip) =     qdd(ip)  * zpr

c         zfac = zs*2.0  + (1.-zs)*3.0
c         zfac = zs*3.5  + (1.-zs)*7.0
          zfac = zs*3.0  + (1.-zs)*8.0
c         zfac = zs*4.5  + (1.-zs)*9.0
c         zfac = zs*10.  + (1.-zs)*20.
          dmelta(ip) = zfac*pdd(ip) / dt

c------------
        endif
c------------

        tice(ip) = min (ta(ip), tmelt)    ! needed for therm

c          Increment snow thickness
 
        hsnowa(ip) = hsnowa(ip) + sna(ip)*dt/rhosno
        if (hsnowa(ip).gt.0.) hsnowa(ip) = hsnowa(ip) 
c    *                        + (raa(ip) - dmelta(ip))*dt/rhosno
     *                        + (         -dmelta(ip))*dt/rhosno
        hsnowa(ip) = max (hsnowa(ip), 0.)
        tsno(ip) = tice(ip)

c==========
      enddo 
c==========
 
      return
      end

c----------------------------------------------------------------------

      subroutine degree_day (dt, ta, pdd, qdd, t0, np)

c        Degree-day parameterization, for np points, involving  
c        error function for normally-distributed diurnal/synoptic 
c        variability. Degree-day increments returned in pdd,
c        and snowfall/precip ratio in qdd.
c        t0 is zero-pdd temperature (deg K) (supplied).

      use comicephys
      use comicegrid

      dimension ta(np), pdd(np), qdd(np)

      parameter (nerf=3000, xerf=3., derf = xerf/nerf)
      dimension erf(0:nerf)
      save erf, sdt2, sdp2, pisq
c     parameter (sdt=5.5, sdp=4.5)
      parameter (sdt=5.0, sdp=5.0)
c     parameter (sdt=4.0, sdp=4.0)

      logical first
      data first /.true./
      save first

      if (first) then
        erf (0) = 0.
        do i=1,nerf
          if (i.eq.1) then
            zx  = 0.25*derf
            zdx = 0.5 *derf
          else
            zx  = (i-1)*derf
            zdx = derf
          endif
          erf(i)  = erf(i-1)  +   exp(-zx*zx)*zdx
        enddo
        sdt2 = sdt*sqrt(2.)
        sdp2 = sdp*sqrt(2.)
        pisq = sqrt(pi)
        first = .false.
      endif

      do ip=1,np

c       positive degree days:
        ztm = ta(ip) - t0
        zx = abs(ztm/sdt2)
        ierf = max (1, min (nerf, nint(zx/derf + 0.5) ))
        zpdd = 0.5*(sdt2/pisq) * exp(-(ztm/sdt2)**2)
     *       + (ztm/pisq) * (0.5*pisq + sign(1.,ztm)*erf(ierf))
        if (zpdd.lt..001) zpdd = 0.      !protect against roundoff error
        pdd(ip) = zpdd*(dt/86400.)

c       snowfall:precip ratio (use tmelt for cutoff value):
        ztp = ta(ip) - tmelt
        zx = abs(ztp/sdp2)
        ierf = max (1, min (nerf, nint(zx/derf + 0.5) ))
        zaop = (1./pisq) * (0.5*pisq - sign(1.,ztp)*erf(ierf))
        if (zaop.lt..001) zaop = 0.      !protect against roundoff error
        if (zaop.gt..999) zaop = 1.      !protect against roundoff error
        qdd(ip) = zaop

      enddo

      return
      end

c-----------------------------------------------------------------------

c*****
#else
c*****
      subroutine climate_uvic
      return
      end
c*****
#endif
c*****
