c----------------------------------------------------------------------- 


c preprocessor defines:

c LONLAT      = lon,lat grid (must have one of LONLAT,STEREO)
c STEREO      = polar stereographic grid
c LINSHELF    = linear shelf rheology (still non-linear in sheet)
c NOARRHEN    = no dependence of ice rheol coeff on T (subr arrhenius)
c NOBASET     = no dependence of basal sliding coeff on basal temp
c NOWATER     = no non-ocean water anywhere (open lakes or sub-ice)
c NOMOVET     = no advection or diffusion of temps (icetherm,vdif,hdif)
c NOMOVEW     = no advection of liquid water (movewater)
c SEDIMENT    = deforming sediment model
c SEDTRACK    = seiment tracking (for strata) (if SEDIMENT)
c BEDROCK     = bedrock model (locally relaxing asthenosphere, and...
c IFLITH      = elastic lithospheric flexure (if BEDROCK)

c NOSHELF     = no floating ice/shelf dyn, just SIA (w/o -DSCHOOFGL)
c SCHOOFGL    = Schoof (2007,JGR) imposed velocities at grounding line 

c GCMMATRIX   = surface budget from GCM meteo files with matrix method
c EULER       = simple Euler (upstream) advection of tracers (ice temp)
c               in icetherm, not parabolic upstream.
c WRAPAROUND  = E-W (dateline) wraparound in iceolddyyn and icetherm  
c               (not yet in icedyn)
c NOLONG      = new dynamics, but no longitudinal stretching (SAI)
c NUMREC,SLAP,SUPERLU,WATSON,MKL,GAUSS = sparse matrix solvers
c SUN         = Sun FPE hanlding
c---

c EIS1-4      = Eismint I shelf tests 
c EISA-D      = Eismint I sheet tests (and EISVAR[20,40] grid sizes) 
c EISANTA     = 2-D Antarctica 
c NHA         = 2-D Northern Hemisphere (using LONLAT and GCMMATRIX)
c CARB        = 2-D global (Permo-Carb) (using LONLAT and GCMMATRIX)

c---

c BEDMAP      = Bedmap Antarc data (with EISANTA or EISLINE/TRANSECTA)  
c IDEAL       = piecewise profile, not BEDMAP (with EISLINE/TRANSECTA)
c EISMINT2    = Eismint II Antarc data (with EISANTA)

c FORCEPLEIST  = with EISANTA, Specmap (sealev) and Epica (dtantann). 
c FORCELGM     = with EISANTA, sealev=-125 and dtantann=-10.
c FORCEMIS31   = with EISANTA, MIS31 forcing from Lisiecki+Raymo 
c                (sealev) and Laskar orbit/insol (dtantann, dtantjan).
c FORCEEO      = with EISANTA, Eocene-Oligocene transition
c BASACCUM     = with EISANTA, Arthern modern obs Antarc accumulation
c                (budgbas).

c---

c HEINO       = ISMIP/HEINO 2-D ice-sheet + Heinrich Event tests 
c ..._T1,T2,B1,B2,S1,S2,S3 = perturbs to temps, budg, basal sliding 
c HEINOSMALL  = use subset of domain

c---

c HOM[A-F]    = ISMIP/HOM 2-D tests

c---

c EISLINE & ...  = 1-D flowline (mostly marine, WAIS-like):

c EISLINE, LINEB = Quaternary sawtooth cycles (weirun: 0=modern, 1=lgm)

c EISLINE, LINEC = Cenozoic "orbital" cycles (weirun: 0=warm, 1=cold)

c EISLINE, LINED = basic mountain/ocean profile, quarrying+transport

c EISLINE, LINEE = wedge simulations, time indep or dep
c                  Either INITEE, RESTARTE, SEALEVCHANGE (time indep),
c                  or none of those (time dependent).
c   INITEE         = no wedge, time indep
c   RESTARTE       = with wedge, time indep
c   WEDGEnnKM      = wedge extent, nn = 4,8,12,16 (km), time indep
c   SEALEVCHANGE,RISE,FALL = instantaneous s.l. change, time indep

c EISLINE, LINEF = grounding line behavior, divide at l.h.edge
c                    grounding line (dzoom, etc, in comicegrid.h)
c   lhs_uh (comments) = specified h and u at l.h.s. domain

c EISLINE, LINEG = sub-glacial lakes and hydrology  

c EISLINE, TRANSECTA = WAIS Weddell-Ross Sea transect, with BEDMAP 

c O_ism = sheetshelf/uvic model coupling

c---

c-----------------------------------------------------------------------

#if defined O_ism

# if defined O_ism_uncoupled
      program sheetshelf
      
      use comicephys
      use comicegrid
      use comicesparse
      use comicesheet
      use ice_sheet_storage, only: strtis,endis          
# else
      subroutine ism(isn)
      
      use comicephys
      use comicegrid
      use comicesparse
      use comicesheet
      use ice_sheet_storage, only: strtis,endis 
            
      include "switch.h"   !contains eorun, init
# endif

# if defined CO2INTER 
      dimension
     *  tempmon(nx,ny,nmon), precmon(nx,ny,nmon)
      logical firstco2
      save firstco2
      data firstco2 /.true./
# endif
# if defined O_ism
#  if defined O_ism_pdd
      !note: must make these allocatable if using them
      real tempmon(nx,ny,nmon),precmon(nx,ny,nmon),topogmon(nx,ny)
#  endif
      character (120) fname      
# endif
      integer isn,ntrec
      save ntrec

      logical ifexist

# if defined LINEE && defined RESTARTE
#  include <fort.26x>
# endif

# if defined SUN
      external hand1, hand2, hand3
# endif

c-------------------------
c Start of executable code
c-------------------------

c#########################
# if !defined O_ism_uncoupled
c     start initialization       
      if (firstcall) then
# endif
c#########################

c        Set floating-point exception handler for sun (solaris) 

# if defined SUN
      ier = ieee_handler ('set', 'division', hand1)
c     ier = ieee_handler ('set', 'overflow', hand2)
c     ier = ieee_handler ('set', 'invalid',  hand3)
      ier = ieee_handler ('set', 'overflow', hand1)
      ier = ieee_handler ('set', 'invalid',  hand1)
# endif

c        Delete stub file that will signify normal end (not abend)

      open (iuokend, file='okend', status='unknown') 
      close (iuokend, status='delete') 
      print*, '------------------------------------------------'
      print*, 'Initializing', isname
      print*, '------------------------------------------------'
      print*, 'nyearstart=', nyearstart      
      print*, 'nyearres=', nyearres        
      print*, 'ifrest=', ifrest      
      print*, 'dtimeice=',dtimeice         
      print*, 'ice_accelerator=',ice_accelerator  
      print*, 'nyearhis=', nyearhis        
      print*, 'nyearout1d=', nyearout1d      
      print*, 'nyearout2d=', nyearout2d      
      print*, 'nyeartab=', nyeartab        
      print*, 'nyearplot1d=', nyearplot1d     

      nyearend   = 10000000 !set to longer than plausible coupled run
      nloopendall= 10000000 !ditto... should be removed eventually, 
                              !but used during writing of diagnostics 
                              !within loop
      if (isn==endis) ntrec=1 !Jer diagnostic

      call izero (hislist, 1000)
      
      nhislist = 0
      do m=1000,1,-1
        if (hislist(m).ne.0) then
          nhislist = m
          goto 100  
        endif
      enddo
  100 continue

c        Set remaining run variables 

      dtimeclim_gcm = 200.    ! must divide into 1000 yr (icegcm)
      dtimeclim_ebm = 200.    ! to run climate_ebm (surface budgets)  
      dtimerun_ebm = 2000.    ! to run ebm model
# if defined NESTING
      dtimetherm   = dtimeice
# else
#  if !defined O_ism_standalone
      dtimetherm   = dtimeice !Jer
#  elseif defined O_ism_standalone
      dtimetherm   = 5. !Jer
#  endif      
# endif
      dtimehyd     = dtimeice
# if defined NHA || defined CARB
      dtimebed     = 50.
# else
      dtimebed     = 10.  !Jer
# endif
      dtimesed     = dtimeice
      dtimeocn     = dtimeice
      nyearsedbud  = nyeartab
# if defined CO2INTER
      dtimeco2 = 200.
# endif

c        Initialize horizontal and vertical grids, 
c        and physical constants in commons

      call initgrid (.true., xglu)

c        Set geoflux (used in initphys below, and in icetherm/vdif).  geoflux_unif read in from namelist
      if (isname == 'Antarctica') then
      
        do j=1,ny
          do i=1,nx
# if defined EISANTA
c         higher geothermal flux under all WAIS:
            if ( alatd(i,j).gt. -86. .and.
     *           (alond(i,j).gt.170. .or. alond(i,j).lt.-30.) ) then
              geoflux(i,j) = geoflux_wais
            else
              geoflux(i,j) = geoflux_unif
            endif
# else
            geoflux(i,j) = geoflux_unif
# endif
          enddo
        enddo
      elseif (isname == 'Greenland') then 
# if defined O_ism_data_geoflux
        !overwrite prescribed geoflux value with spatial data from file
        fname='bheatflx'
	call bamberaggreg (fname, geoflux)       
	geoflux = geoflux*31556926. !Jer: convert to W/year presumably	
# else
	geoflux(:,:)=geoflux_unif
# endif
      endif      

c        Initialize physical variables, or read restart file 
      !if first ice sheet through, set iceinit to false.  Will switch to true if either ice sheet is initialized from data.
      if (isn==strtis) iceinit=.false.
c~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (ifrest.eq.-1) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~
        iceinit=.true.
        nwrit = 0     ! current number of writes on history file

        nloopstart = 1
        timeice    = nyearstart

        call initphys (h, hb, t, 
     *                 sedim,  sedimeq, sedimun, tsed, wsed, 
     *                 topbed, topbedeq, tbed, equiload, geoflux,
     *                 hw, tw, maskpres, sedpres, 
     *                 ua, va, sealev_notused, timeice)

c       following are used in icetherm before first set in icedyn
        call zero (ub, (nxp+1)*(nyp+1))
        call zero (vb, (nxp+1)*(nyp+1))
        call zero (dfu,(nxp+1)*(nyp+1))
        call zero (dfv,(nxp+1)*(nyp+1))
        call zero (heatb,nx*ny)
        call zero (heats,nx*ny*nsed)       ! set in sedflow
        call zero (basefrml,nx*ny)         ! set in vdif
        call zero (basefrml,nx*ny)         ! set in vdif

        rco2inter = 4.
c~~~~~~~~~
      else
c~~~~~~~~~

c---
# if defined HEINO 
        call readheino
# endif
        call readres  (h, hs, hb, t, 
     *                 baseperc, basefrml, heatb,
     *                 sedim, sedimeq, sedimun, tsed, wsed, heats, 
     *                 topbed, topbedeq, tbed, equiload,
     *                 hw, tw, maskpres, sedpres,
     *                 ua, va, ub, vb, dfu, dfv,
     *                 sealev_notused, rco2inter, nwrit, timeicein)

        if (ifrest.eq.0) rco2inter = 4.

# if defined LINEE && defined RESTARTE
      xglwedge1 = 1.e20
      xglwedge2 = 1.e20
      do j=1,ny
        do i=1,nx
          if (hw_wedge(i).gt.0.) then
            if (xglwedge1.eq. 1.e20) xglwedge1 = xh(i,j)-0.5*dx(i,j)
            if (xh(i,j)-xglwedge1 .gt. 0. .and.
#  if defined WEDGE4KM
     *          xh(i,j)-xglwedge1 .lt. 4.e3) then     !  4 km
#  elif defined WEDGE8KM
     *          xh(i,j)-xglwedge1 .lt. 8.e3) then     !  8 km
#  elif defined WEDGE12KM
     *          xh(i,j)-xglwedge1 .lt. 12.e3) then    ! 12 km
#  elif defined WEDGE16KM
     *          xh(i,j)-xglwedge1 .lt. 16.e3) then    ! 16 km
#  endif
              topbedeq(i,j) = hswedge(i) - hwedge(i) + .01
              topbed(i,j) = topbedeq(i,j)
              hb(i,j) = topbedeq(i,j)
              hw(i,j) = 0.
            else
              if (xglwedge2.eq. 1.e20) then
                xglwedge2 = xh(i-1,j)
                zglwedge2 = topbedeq(i-1,j)
              endif
              topbedeq(i,j) = max (topbedeq(i,j),
     *                             zglwedge2 - .01*(xh(i,j)-xglwedge2))
              topbed(i,j) = topbedeq(i,j)
              hb(i,j) = topbedeq(i,j)
              hw(i,j) = hs(i,j) - h(i,j) - topbed(i,j)
            endif
          endif
        enddo
      enddo
#  endif

        if (ifrest.eq.1) then
          nyearendin = nint(timeicein)+isyear0
	  if (isyear0 .ne. 0.)
     &	    print*,'Init time changed from',timeicein,' to',nyearendin
          if (nyearend.le.nyearendin) then
            write (6,'(a/a,i8,a,i8)')
     *      '*** Error: end date is earlier than restart-file end date',
     *      '    nyearend=',nyearend, '  nyearendin=',nyearendin
            stop
          endif
          nloopstart = nint((nyearendin-nyearstart)/dtimeice) + 1
          timeice    = nyearendin
        else
          nwrit = 0
          nloopstart = 1
          timeice    = nyearstart
        endif
# if !defined O_ism_uncoupled
        nwrit = 0 !overwrite nwrit to force creation of a new history
                  !file during each coupled model initialization 
# endif

c~~~~~~~~~~
      endif  ! ifrest  
c~~~~~~~~~~

c        Set initial maskwater, indlake, etc (not on restart file),
c        and adjust hw and hs (adjustpres)

      call findwater (maskwater, indlake, npoilake, nlake,
     *                h, hb, hw, sealev, timeice)

      call adjustpres (maskwater, indlake, npoilake, nlake,
     *                 h, hb, hw, hs, sealev)

      call scopy_i (nx*ny, maskwater, 1, maskdisplay, 1)  ! for printmap

c^^^^^^^^^^^^^^^^^^^^^^^
# if defined O_ism_uncoupled
c^^^^^^^^^^^^^^^^^^^^^^^
      nloopend    = nint((nyearend-nyearstart)/dtimeice)
      nloopendall = nint((nyearendall-nyearstart)/dtimeice)
      
c^^^^^
# endif
c^^^^^            
# if defined O_ism_glacial_transient
      !Jer: load Lisiecki/Raymo (2005) weight curve
      call set_uvic_lr_clim 
      call uvic_lr_clim(zclim,sealev)
# endif      
          
c#######################

# if !defined O_ism_uncoupled
c     end initialization
      firstcall = .false.
! !       
      return
      endif
# endif
c#######################
      facice = 0.
      facorb = 0.
      facco2 = 0.
cvvvvvvvvvvvvvvvvvvvvvvv
# if !defined O_ism_uncoupled
cvvvvvvvvvvvvvvvvvvvvvvv
# if !defined O_ism_standalone
      if (mod(real(ice_accelerator),dtimeice).gt.1.e-15)
     &  stop '->icectl: non-integer ice sheet loop counter' 
      !set number of ice sheet timesteps to carry out
      nloopend = nloopstart - 1 
     &+ ice_accelerator/dtimeice  
# else 
      !override timeice from input restarts
      timeice=-100000.
      print*, 'Warning: u r on Lisiecki/Raymo time now;'
      print*, 'Timeice:', timeice      
      !set number of ice sheet timesteps to carry out
      nloopstart=1
      nloopend = nloopstart + abs(timeice)/dtimeice - 1    
# endif                  
cvvvvv
# endif
cvvvvv
c        Overall time loop
      print*, 'ISM running ',isname    
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do iloop=nloopstart,nloopend
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c          Impose constraints for particular experiments

# if defined EIS1 || defined EIS2
        do j=1,ny
          h(1,j) = 500.*(1. + 0.2*sin(2.*pi*timeice/200.))
        enddo
# elif defined EISA
        if (nx.gt.1) then
          do j=1,ny
            h(1, j) = 0.
            h(nx,j) = 0.
          enddo
        endif
        if (ny.gt.1) then
          do i=1,nx
            h(i, 1) = 0.
            h(i,ny) = 0.
          enddo
        endif
# elif defined EISLINE && ! defined TRANSECTA
        if (nx.gt.ny) then
          do j=1,ny
            h(nx,j) = 0.  ! comment out for 999_1pt
          enddo
        else
          do i=1,nx
            h(i,ny) = 0.  ! comment out for 999_1pt
          enddo
        endif
# elif defined EISLINE && defined TRANSECTA
        do j=1,ny
c         h(1,j)  = 0.     ! comment out for central divide
          h(nx,j) = 0.
        enddo
# endif

        call scopy (nx*ny, sedim, 1, sedimold, 1)


c++++++++++++++++++++++++
# if defined O_ism_uncoupled
c++++++++++++++++++++++++

c          Set long-term weighting (weirun, or sea level, dt's, co2,
c          and orbit for EISANTA)
        call longterm (weirun, sealev, dtantann, dtantjan, dtseas, 
     *                 rco2, ecc, obl, prec,
     *                 timeice, nyearstart, nyearend, iutsin)
#  if defined CO2INTER
        rco2 = rco2inter
#  endif

c          Set climate (budg[snow,rain,evap,melt], tsurf[,i,o])
#  if defined GCMMATRIX
        n = max (1, nint (dtimeclim_gcm/dtimeice))
        if (mod (iloop, n).eq.0 .or. firstclim) then
          call climate_gcm (budgsnow, budgrain, budgevap, budgmelt,
     *                      tsurf, tsurfi,
     *                      h, hs, maskwater,
     *                      sealev, dtantann, dtantjan, dtseas,
     *                      rco2, ecc, obl, prec,
     *                      facice, facorb, facco2,
#   if defined CO2INTER
     *                      tempmon, precmon,
#   endif
     *                      timeice, weirun, ifrest, firstclim)
          firstclim = .false.
        endif
#  elif defined SUBICE
        n = max (1, nint (dtimeclim_ebm/dtimeice))
        if (mod (iloop, n).eq.0 .or. firstclim) then
          call climate_ebm (budgsnow, budgrain, budgevap, budgmelt,
     *                      tsurf, tsurfi,
     *                      h, hs, topbedeq, maskwater,
     *                      sealev, dtantann, dtantjan, dtseas,
     *                      rco2, ecc, obl, prec, 
     *                      timeice, iloop, dtimeice, idtimerun_ebm)
          firstclim = .false.
        endif
#  else
        call climate_simp (budgsnow, budgrain, budgevap, budgmelt,
     *                     tsurf, tsurfi,
     *                     h, hs, maskwater,
     *                     sealev, dtantann, dtantjan, dtseas,
     *                     rco2, ecc, obl, prec,
#   if defined CO2INTER
     *                     tempmon, precmon,
#   endif
     *                     timeice, weirun)
#  endif

#  if defined CO2INTER
        n = max (1, nint (dtimeco2/dtimeice))
        if (mod (iloop, n).eq.0 .or. firstco2) then
          call doco2 (tempmon, precmon, h, maskwater,
     *                rco2inter, ecc, obl, prec, timeice, dtimeco2)
          firstco2 = .false.
        endif
#  endif

c+++++
# else
c+++++
c     for uvic, used in calculation of ocean melt in subr ocean:
      dtantann = 0.
      dtantjan = 0.
      dtseas   = 20. ! not used 
      rco2     = 1.
#  if defined O_ism_glacial_transient
        call uvic_lr_clim(zclim,sealev)
#  endif
#  if defined fade2glac
      fyri=0.
      fyre=100.
      zclim = max((fyre-timeice),0.)/(fyre-fyri)
      sealev = -125.* (1.-zclim)                 
# endif  
c+++++
# endif
c+++++
c          Step lithospheric and asthenosheric bedrock response  
        
# if defined BEDROCK
        n = max (1, nint (dtimebed/dtimeice))
        if (mod (iloop, n).eq.0) then
          call bedrock (h, hb, hw, topbed, topbedeq, sedim, sedimeq,
     *                  equiload, dtimebed)
c         adjust hw and hs for change in bedrock elev (and so hb):
          call adjustpres (maskwater, indlake, npoilake, nlake,
     *                      h, hb, hw, hs, sealev)
        endif
# endif

c          Step liquid water (ocean, open lakes, sub-ice hydrol/lakes)
        n = max (1, nint (dtimehyd/dtimeice))
        if (mod (iloop, n).eq.0) then

# if defined NOMOVEW
          call zero (heatw,nx*ny)
          call zero (uw, (nxp+1)*(nyp+1))
          call zero (vw, (nxp+1)*(nyp+1))
# else
          call movewater (maskwater, h, hb, hw, tw, 
     *                    heatw, baseperc, basefrml, 
     *                    uw, vw, dtimehyd)
# endif

c         reset ocean/lake locations and indices 
          call findwater (maskwater, indlake, npoilake, nlake,
     *                    h, hb, hw, sealev, timeice)

c         adjust hw and hs for change in hw:
          call adjustpres (maskwater, indlake, npoilake, nlake,
     *                     h, hb, hw, hs, sealev)
        endif
c          Step temperatures t,tw,tsed,tbed, set basefrml,wsed

        n = max (1, nint (dtimetherm/dtimeice))
	!don't run icetherm if not starting from IC: instead, using read-in values for t, baseperc, basefrml.
	if (ifrest .ne. -1) firsttherm=.false.
        if (mod (iloop, n).eq.0 .or. firsttherm) then
          call icetherm (h, hs, t,  hw, tw, 
     *                   budgsnow, budgrain, budgevap, budgmelt, 
     *                   baseperc, basefrml,
     *                   tsurf, tsurfi,
     *                   arhap,  s2a0,
     *                   heati, heatb, heatw, geoflux,
     *                   w, wa,
     *                   sedim, tsed, wsed, heats, tbed,
     *                   u, v, ub, vb, dfu, dfv,
     *                   timeice, dtimetherm)
          firsttherm = .false.
        endif
c          Set crh[u,v] terms for basal sliding, used in icedyn 
        call basecoef (h, hb, hw, t, tsurfi, sedpres, sealev, 
     *                 crhu, crhv, powbu, powbv, 
     *                 fracgu, fracgv, hgu, hgv)

c          Sub-glacial sediment deformation. Step sediment thickness, 
c          calculate frictional heating for therm/vdif

# if defined SEDIMENT
        nloopsed = nint (dtimeice/min(dtimesed,dtimeice))
        dt = dtimeice/nloopsed
        do iloopsed = 1,nloopsed
          call sedflow (h, hb, topbed, topbedeq, sedim, sedimun, 
     *                  disturb, t,
     *                  dfu, dfv, crhu, crhv, fsedu, fsedv,
     *                  maskwater, quarryrate,
     *                  heats, heatb, dt, 1)
        enddo
c          Move sediment within ocean and lakes

        call sedocean (sedim, sedimun, h, hb, topbed, topbedeq,
     *                 maskwater)

c          Adjust hw and hs (for change in sedim, hence hb):

        call adjustpres (maskwater, indlake, npoilake, nlake, 
     *                   h, hb, hw, hs, sealev)
# endif

c          Set oceanmelt and calvrate 

        n = max (1, nint (dtimeocn/dtimeice))
        if (mod (iloop, n).eq.0 .or. firstocn) then
          call ocean (zclim,h, hb, hw, maskwater, maskpres,
     *                oceanmelt, calvrate, sealev, 
     *                dtantann, dtantjan, rco2,
     *                timeice, weirun, dtimeocn, firstocn, iloop)
          firstocn = .false.
        endif

c          Set net ice budget for convenient passing to icedyn, output.
c          Nb: budgmelt and budgrain have contributed to baseperc
c          in icetherm/vdif (= budgmelt+budgrain if NOMOVET or if h=0).

        do j=1,ny
          do i=1,nx
            budgall(i,j) = budgsnow(i,j) - budgevap(i,j)
     *                   + budgrain(i,j) - baseperc(i,j) + basefrml(i,j)
     *                   - oceanmelt(i,j) -calvrate(i,j)	     
     
# if defined NESTING
c           so h around edges is unchanged in icedyn: 
            if (i.eq.1 .or. i.eq.nx .or. j.eq.1 .or. j.eq.ny)
     *        budgall(i,j) = 0.
# endif
          enddo
        enddo
	
!         if (iloop==nloopend) then
!           fname = 'budgsnow_'// trim(isname)
!           call jerncdf_ice_timedep(budgsnow,fname,fname,ntrec)  
!         endif
c          Step ice dynamics, update ice thickness h and sfc elev hs,
c          set heatb,ub,vb for icetherm) 

        call mascalc (h, toth0, tota0)
        call icedyn (h, hs, hb, hw, sedim, t,
     *               arhap, s2a0, heatb, budgall,
     *               maskh, dfu, dfv,
     *               ua, va, ui, vi, ub, vb, uadv, vadv, 
     *               hu, hv, uw, vw, masku, maskv, 
     *               crhu, crhv, fsedu, fsedv,
     *               powbu, powbv, fracgu, fracgv, 
     *               muind, mvind, thetau, thetav, hgu, hgv,
     *               fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *               numh, itera, iterc, totflow, totneg,
     *               maskwater, indlake, npoilake, nlake, sealev,
     *               timeice, dtimeice, ifrest)
        call mascalc (h, toth, tota)

        timeice = timeice + dtimeice

c          Write ascii output, 1-D profiles (1-D or 2-D domains)
! 	call iceshow1d (h, hs, hb, t,
!      *  		basefrml, oceanmelt, budgall,
!      *  		tsurf, tsurfi, heati, heatb,
!      *  		w, wa, maskh,
!      *  		sedim, sedimun, tsed, wsed, heats, 
!      *  		quarryrate,
!      *  		topbed, topbedeq, tbed,
!      *  		hw, tw, maskwater,
!      *  		u, v, ua, va, ui, vi, ub, vb,
!      *  		uadv, vadv, hu, hv, dfu, dfv,
!      *  		muind, mvind, thetau, thetav,
!      *  		crhu, crhv, fsedu, fsedv, 
!      *  		numh, itera, iterc,
!      *  		toth0, tota0, toth, tota, totflow, totneg,
!      *  		sealev, dtantann, dtantjan, dtseas, rco2,
!      *  		timeice, dtimeice, weirun,
!      *  		iloop, nloopendall, nyearout1d)

c          Write ascii output, 2-D maps  (2-D domains only - 
c          does nothing if nyearout2d = 0, as set for 1-D domains)

!        call iceshow2d (h, hs, hb, t, 
!     *                  budgall, budgsnow, budgrain, budgmelt,
!     *                  baseperc, basefrml, oceanmelt, calvrate, 
!     *                  tsurf, tsurfi, heati, heatb,
!     *                  w, wa, 
!     *                  sedim,
!     *                  quarryrate,
!     *                  topbed, topbedeq,
!     *                  hw, tw, maskwater,
!     *                  u, v, ua, va, ub, vb,
!     *                  uadv, vadv, hu, hv, dfu, dfv,
!     *                  thetau, thetav, hgu, hgv,
!     *                  fluxgrdu, fluxgrdv, fluxschu, fluxschv,
!     *                  crhu, crhv, 
!     *                  numh, itera, iterc,
!     *                  toth0, tota0, toth, tota, totflow, totneg,
!     *                  sealev, timeice, dtimeice, weirun,
!     *                  iloop, nloopendall, nyearout2d)


! # if defined SEDIMENT
! c       Accumulate and write (if time) sediment budget quantities:
!         call sedbudg (sedim, timeice, dtimeice, iloop, 
!      *                nloopstart, nloopendall, nyearsedbud)
! 
! c       Track dated strata in sediment, write to output file (if time):
! #  if defined SEDTRACK
!         call sedtracking_eul (h, hs, hb, 
!      *                        sedim, sedimold, 
!      *                        timeice, dtimeice, weirun, sealev, 
!      *                        nyearstart, iloop, nloopendall)
! #  endif
! # endif

c       1D-flowline plotting output (1-D or 2-D domains):

!        call iceplot1d (h, hs, hb, t,  
!     *                  budgall, basefrml, 
!     *                  oceanmelt, calvrate,
!     *                  tsurf, tsurfi, heati, heatb,
!     *                  w, wa, maskh,
!     *                  sedim, sedimun, tsed, wsed, heats, quarryrate,
!     *                  topbed, topbedeq, tbed,
!     *                  hw, tw,
!     *                  u, v, ua, va, ui, vi, ub, vb, 
!     *                  uadv, vadv, hu, hv, 
!     *                  masku, maskv, crhu, crhv, 
!     *                  fracgu, fracgv, thetau, thetav,
!     *                  sealev, dtantann, dtantjan, dtseas, rco2, 
!     *                  timeice, dtimeice, weirun, 
!     *                  iloop, nloopendall, nyearplot1d)

c       tabular output:

# if defined EISLINE

!        call icetab1d (h, hs, hb, t,  
!     *                 budgall, basefrml, 
!     *                 oceanmelt, calvrate,
!     *                 tsurf, tsurfi, heati, heatb,
!     *                 w, wa, maskh,
!     *                 sedim, tsed, wsed, heats, 
!     *                 topbed, topbedeq, tbed,
!     *                 hw, tw, maskwater,
!     *                 u, v, ua, va, ui, vi, ub, vb, 
!     *                 uadv, vadv, hu, hv, 
!     *                 masku, maskv, muind, mvind, crhu, crhv, 
!     *                 fracgu, fracgv, thetau, thetav, 
!     *                 fluxgrdu, fluxgrdv, fluxschu, fluxschv,
!     *                 sealev, dtantann, dtantjan,  dtseas, rco2,
!     *                 timeice, dtimeice, weirun,
!     *                 iloop, nloopendall, nyeartab)

# elif defined EISANTA || defined NHA || defined CARB || defined HEINO

!        call icetab2d (h, hw, maskwater,
!     *                 sealev, dtantann, dtantjan, dtseas, rco2,
!     *                 ecc, obl, prec, facice, facorb, facco2,
!     *                 timeice, dtimeice, weirun,
!     *                 iloop, nloopendall, nyeartab, ifrest)

# endif

# if defined HEINO 
!      call iceshowheino (h, hs, hb, t, heatb, u, v, ub, vb, 
!     *                   timeice, dtimeice, iloop, nyearendall) 
# endif

# if defined HOM 
!      call iceshowhom (h, hs, hb, u, v, ub, vb, wa, dfu, dfv,
!     *                 timeice, dtimeice, iloop, nloopendall)
# endif

c          Write to  history file

        call writehis (h, hs, hb, t,  
     *                 budgall, basefrml, 
     *                 oceanmelt, calvrate, 
     *                 tsurf, tsurfi, heati, heatb,
     *                 w, wa, sedim, tsed, wsed, heats, 
     *                 topbed, topbedeq, tbed,
     *                 hw, tw, maskwater, 
     *                 u, v,  ua, va, ub, vb, uadv, vadv, 
     *                 timeice, dtimeice, nwrit,
     *                 iloop, nloopendall, nyearhis, hislist, nhislist, 
     *                 ifrest)

c          Write restart file
        ismrflag=.false.
        if(
     *	  ((nyearres.ne.0) 
     *	     .and. 
     *    ((mod(abs(timeice)+0.5*dtimeice,max(float(nyearres),dtimeice))
     *        .lt.dtimeice)))
     *       .or. 
     *    ((eorun) .and. iloop.eq.nloopend) 
# if defined O_ism_uncoupled      
     *         .or. iloop.eq.nloopend) then
# else
     *	   ) then
# endif
# if !defined O_ism_standalone
          ismrflag=.true.
# else
          call writeres
# endif
        endif
c>>>>>>>>>>
      enddo
c>>>>>>>>>>

cvvvvvvvvvvvvvvvvvvvvvvv
# if !defined O_ism_uncoupled
cvvvvvvvvvvvvvvvvvvvvvvv      
      !ready loop counter for next ice sheet timestep
      nloopstart=nloopend+1   
cvvvvv
# endif 
cvvvvv
      if(isn==endis) ntrec=ntrec+1    
c        Write file for nested runs using this large-scale run

!# if ! defined NESTING
!      call writenest (h, hb, t, ua, va)
!# endif
# if defined O_ism_uncoupled
c     create stub file to signify normal end:
      open (iuokend, file='okend', status='new') 
      stop
      end
# else
      return
      end subroutine ism
# endif

c-----------------------------------------------------------------------

      subroutine readbas (alond, alatd, accum, nx, ny, iu)

c        Reads Arthern et al (BAS) Antarctic snow accumulation
c        data file. Transfers to model grid by nearest neighbor
c        using lon and lat.

      dimension alond(nx,ny), alatd(nx,ny), accum(nx,ny)

c     input data:
      parameter (ninmax = 150000)
      dimension alonin(ninmax), alatin(ninmax), accumin(ninmax),
     *          xin(ninmax), yin(ninmax), zin(ninmax)

      parameter (tmelt = 273.16)
      parameter (pi = 3.14159265358979)

      character*200 cfilin

# if defined SUN || defined CRAY || defined IGUANA
      cfilin = '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Data_bas/'
     *      // 'amsr_accumulation_map.txt'
# elif defined LION
       cfilin = '/home2/dxp21/Sheetshelf/Data_bas/'
     *       // 'amsr_accumulation_map.txt'
# elif defined NCOM
       cfilin = '/home/geosc/tbellato/Sheetshelf/Data_bas'
     *       // 'amsr_accumulation_map.txt'
# endif

c        Read BAS input file

      open (iu, file=cfilin, status='old')
      do  iskip=1,21
        read (iu,*)
      enddo
      nin = 0
  100 continue
      read (iu,*,end=200,err=200) zlat, zlon, izx, zy, zacc, zerr
      nin = nin + 1  
      if (nin.gt.ninmax) then
        write (6,'(a,2i6)') 'Error readbas, nin,ninmax=',nin,ninmax
        stop
      endif
      alatin(nin) = zlat
      alonin(nin) = zlon
      accumin(nin) = zacc
      xin(nin) = cos(alatin(nin)*pi/180.)*cos(alonin(nin)*pi/180.)
      yin(nin) = cos(alatin(nin)*pi/180.)*sin(alonin(nin)*pi/180.)
      zin(nin) = sin(alatin(nin)*pi/180.)
c     write (6,'(i6,3f10.2)') nin, zlat, zlon, zacc
      goto 100
  200 close (iu)
c     write (6,'(a,i6)') 'nin=',nin

c--------

c         Interpolate to model grid (nearest neighbor)

      do j=1,ny
        write (6,'(a,i3)') 'doing j=',j
        do i=1,nx
           x1 = cos(alatd(i,j)*pi/180.)*cos(alond(i,j)*pi/180.)
           y1 = cos(alatd(i,j)*pi/180.)*sin(alond(i,j)*pi/180.)
           z1 = sin(alatd(i,j)*pi/180.)
           zdistmin = 1.e20

           do k=1,nin
             zdist = 
c    *               sqrt ( 
     *                  (x1-xin(k))**2 + (y1-yin(k))**2 + (z1-zin(k))**2
c    *                    )
c            zdist = acos (x1*xin(k) + y1*yin(k) + z1*zin(k)) ! slower
             if (zdist.lt.zdistmin) then
               zdistmin = zdist
               kmin = k
             endif
           enddo

           if (nint(accumin(kmin)).eq.-999) then 
c            if null bas value, use crude parameterization vs lat 
             zhs = 0.  ! sea level
             zts = tmelt + 34.46 - .00914*zhs
     *                           - .68775*abs(alatd(i,j))
             accum(i,j) = 1.5 * (2.**((zts-tmelt)/10.))           ! m/yr
           else
             accum(i,j) = .001*accumin(kmin)          ! kg/m2/yr to m/yr
           endif
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine longterm (weirun, sealev, dtantann, dtantjan, dtseas, 
     *                     rco2, ecc, obl, prec,
     *                     timeice, nyearstart, nyearend, iutsin)

      parameter (pi     = 3.14159265358979)

c        Sets long-term variables (sealev, dtant*, dtseas, rco2, 
c        orbit, if EISANTA),
c        or long-term weighting (weirun, sealev, if not EISANTA)

      dtantann = 0.
      dtantjan = 0.
      dtseas   = 20.
      ecc  = 0.016706             ! modern
      obl  = 23.4377 * pi/180.    ! modern
      prec = 77.0613 * pi/180.    ! modern

c====================
# if defined EISLINE
c====================

c---
c     constant:
      weirun = 0.
c     weirun = 0.5
c     weirun = 1.
      sealev = 0. - 130.*weirun
c---
c     alternating steps:
c     if ( mod ( int((timeice-nyearstart)*1.e-3/20.), 2 ).eq.1 ) then
c       weirun = 0.
c     else
c       weirun = 1.
c     endif
c     sealev = 0. - 130.*weirun
c---
c     Eismint g.l. steps:
c     sealev0 =    0.
c     sealev0 = -100.
c     sealev0 =  100.
c     dsealev = -125.
c     dsealev =    0.
c     dsealev =  125.
c     kyear = nint(timeice-nyearstart)
c     if (kyear.le.50000) then
c       weirun = 0. 
c     else if (kyear.le.100000) then
c       weirun = min((kyear-50000)/1000., 1.)
c     else 
c       weirun = 1.- min((kyear-100000)/1000.,1.)
c     endif
c     sealev = sealev0 + dsealev *weirun
c---
#  if defined LINEB
c---
c     Quaternary sawtooth cycles (weirun: 0=modern, 1=lgm):
      zper = 120000.
      zsub = 105000.
      zt = mod (timeice-nyearstart, zper)
      if (zt.lt.zsub) then 
        weirun = zt/zsub
      else
        weirun = 1. - (zt-zsub)/(zper-zsub)
      endif

c     gradual s.l. lowering, then gradual rising:
c     ztm = timeice - nyearstart
c     if(ztm.le.125000.) then
c       weirun = max (0., min (1.,    (ztm         )/125000.))
c     else
c       weirun = max (0., min (1., 1.-(ztm -125000.)/125000.))
c     endif

c     sealev = 0. - 130.*weirun  ! old
c     sealev = 0.                ! old

#   if defined SIDEDRAG
c     sealev = 180. -   50.*weirun ! with sidedrag
      sealev = 250. -  125.*weirun ! with sidedrag
#   else
      sealev =  30. -   95.*weirun ! no   sidedrag
#   endif

c     reverse slope:
c     ztm = timeice - nyearstart
c     weirun = max (0., min (1., (ztm - 50000.)/50000.))
c# if defined SIDEDRAG
c     sealev = 180. - 50.*weirun  ! with sidedrag, sl falling
c     sealev = 130. + 50.*weirun  ! with sidedrag, sl rising
c     sealev = 185. - 60.*weirun  ! with sidedrag, sl sawtooth
c# else
c     sealev = -55. -  10.*weirun ! no sidedrag, sl falling
c     sealev =  20. +  10.*weirun ! no sidedrag, sl rising
c     sealev =  35. - 105.*weirun ! no sidedrag, sl sawtooth
c# endif

c---
#  elif defined LINEC
c---
c     Cenozoic "orbital" cycles (weirun: 0=warm, 1=cold):

c     weirun = 0.5 
c     sealev = 0. - 130.*weirun

c     zt = timeice - nyearstart 
c     weirun = 0.5 + 0.5*sin(2.*3.1415927*zt/100000.)
c     sealev = 100. - 200.*weirun

c     ztm = mod (timeice-nyearstart, 200000.)
c     weirun = min ( (ztm - 0.)/5000., (105000.-ztm)/5000. )
c     weirun = max (0., min (1., weirun))
c     sealev = 100. - 200.*weirun

      ztm = mod (timeice-nyearstart, 100000.)
      weirun = max ( (90000.-ztm)/5000., (ztm-95000.)/5000. )
      weirun = 0.5* max (0., min (1., weirun))
c     weirun = 0.5
      sealev = 300.*(0.5-weirun) +  (-130.)*weirun

c     if (timeice-nyearstart.lt.9900000.) then
c       sealev = -65.
c     else    
c       sealev = 150.
c     endif   
c     weirun = 0.5

c---
#  elif defined LINED
c---
c     basic mountain/ocean profile, just quarrying+transport

      weirun = 0.5
      ztm = mod (timeice-nyearstart, 1000000.)
      if (ztm.gt.900000.) weirun = 0.
      sealev = 0.

c     weirun = 0.5
c     if (timeice-nyearstart.gt.9900000.) weirun = 0.
c     sealev = 0.

c---
#  elif defined LINEE
c---
c     Richard Alley sediment wedge WAIS tests
      weirun = 0.5
#   if defined SEALEVCHANGE
#    if defined SEALEVFALL
      sealev = -40. - SEALEVFALL
#    elif defined SEALEVRISE
c     sealev = -40. + 30   ! + 0 5 10 20 30
      sealev = -40. + SEALEVRISE
#    endif
#   elif defined RESTARTE
c     sealev = 0.
      sealev = -40. ! -35
#   elif defined INITEE
      zwei =  max (0., min (1., (timeice-1000.)/1000.))
ca    sealev = -320.*(1.-zwei) - 0.*zwei   ! a
c     sealev = -500.*(1.-zwei) - 0.*zwei
      sealev = -500.*(1.-zwei) - 40.*zwei  ! -35
#   else
      zwei =  max (0., min (1., (timeice-1000.)/(nyearend-1000.)))
      sealev = -320.*(1.-zwei) + 100.*    zwei   ! b,c
c     sealev = -320.*(   zwei) + 100.*(1.-zwei)  ! d,e
c     sealev = min (sealev, -75.)                 ! 777   
c     sealev = -75. + 00.*min (1., timeice/2000.) ! 777     ! 49.9 km
c     sealev = -75. +  5.*min (1., timeice/2000.) ! 777     ! 
c     sealev = -75. + 10.*min (1., timeice/2000.) ! 777     ! 20.9 km
c     sealev = -75. + 25.*min (1., timeice/2000.) ! 777     ! 10.3 km
#   endif
c---
#  elif defined LINEF
c---
      weirun = 0.5

c     zwei =  max (0., min (1., (timeice-20000.)/100.))
c     sealev = -500.*(1.-zwei) - 0.*zwei
c     sealev = -501.*(1.-zwei) - 200.*zwei

c     sealev = -310.  ! on slope
c     sealev = -320.  ! on shelf
c     sealev = -330.  ! to rh edge

c     zwei =  max (0., min (1., (timeice - 0.)/200000.))
c     sealev = -300.*(1.-zwei) - 330.*zwei

c     if (timeice.lt.50000.) then
c       zwei =  max (0., min (1., (timeice - 0.)/50000.))
c       sealev = -220.
c     else if (timeice.lt.100000.) then
c       zwei =  max (0., min (1., (timeice - 50000.)/50000.))
c       sealev = -220.*(1.-zwei) - 233.*zwei ! 240
c     else if (timeice.lt.150000.) then
c       zwei =  max (0., min (1., (timeice - 100000.)/50000.))
c       sealev = -233.
c     else if (timeice.lt.250000.) then
c       zwei =  max (0., min (1., (timeice - 100000.)/50000.))
c       sealev = -233.
c     else
c       zwei =  max (0., min (1., (timeice - 150000.)/50000.))
c       sealev = -233.*(1.-zwei) - 130.*zwei
c     endif

c     sealev = -150.
c     sealev = -150. + 1.*(timeice-40000)/2000.

c     sealev = 0. 
ca    sealev = -100. 
cb    sealev = 0. 
cc    sealev = -200. 

c     sealev = -100.
c     if (timeice.gt.100000.) sealev = -150.

      sealev = 0. 
c     if (timeice.gt.100000. .and. timeice.le.200000.) sealev = -125.
c     if (timeice.gt.100000. .and. timeice.le.200000.) sealev =  125.

c     schoof a,b:
c     if (timeice.gt.100000.) sealev = -125.
c     if (timeice.gt.200000.) sealev =  100.

c     schoof c:
c     zk = timeice/1000.
c     if (zk.lt.200.) sealev = 150. - 350.*zk/200.
c     if (zk.gt.200. .and. zk.lt.300.) sealev = -200.
c     if (zk.gt.300.) sealev = -200. + 350.*(zk-300.)/200.

      if  (timeice.lt.50000.) then
        sealev = 100.   
      else if  (timeice.lt.100000.) then
        sealev = -125.
      else
        sealev = 100.
      endif

c---
#  elif defined LINEG
c---

      weirun = 0.5
      sealev = 0.

c---
#  elif defined TRANSECTA
c---

      weirun = 0.5
c     sealev =    0. ! a

#   if defined SIDEDRAG
      if (timeice.lt.50000.) then
        sealev = -125.  ! lgm
c       sealev = 0.     ! mis31
      else if (timeice.lt.100000.) then
        sealev = 0.     ! mod
      else
        sealev = 0.     ! mis31
c       sealev = -125.  ! lgm
      endif
#   endif

c---
#  else
c---

      weirun = 0.5
      sealev = 0.

c---
#  endif
c---

c======================
# elif defined EISANTA
c======================

c     weirun set in ocean (diagnostic) 
c     dtantann = annual mean insol at 80 S - modern * 0.1 deg C/W/m2  
c     dtantjan = January     insol at 80 S - modern * 0.1 deg C/W/m2  
c     rco2     = atmos CO2 level / preindustrial  

#  if defined FORCEPLEIST
c     call specmap  (sealev,   timeice, iutsin)  
cc    call bassinot (sealev,   timeice, iutsin)
c     call epica    (dtantann, timeice, iutsin)  
      call lisiecki (sealev, timeice, iutsin)  
      call laskar (dtantann, dtantjan, dtseas, ecc, obl, prec, 
     *             timeice, iutsin)
      rco2 = 1.
#  elif defined FORCELGM
      sealev = -125.
      dtantann = -10.
      dtantjan = 0.
      dtseas = 20.
      rco2 = 1.
#  elif defined FORCEMIS31
      call lisiecki (sealev, timeice, iutsin)  
      call laskar (dtantann, dtantjan, dtseas, ecc, obl, prec, 
     *             timeice, iutsin)
      rco2 = 1.
#  elif defined FORCEEO
      sealev = 0.
      call laskar (dtantann, dtantjan, dtseas, ecc, obl, prec,
     *             timeice, iutsin)
      rco2 = max (4. - 2.*(timeice-nyearstart)/1.e6, 1.)

c     synthetic cycles:
c     dtantann = 40. - 30.*timeice/1.e6
c     dtantann = 0.
c     dtantjan =       7.5*cos(2.*pi*timeice/20000.)
c     dtseas   = 20. + 15.*cos(2.*pi*timeice/20000.)
#  else
      sealev = 0.  
      dtantann = 0.
      dtantjan = 0.
      dtseas = 20.
      rco2 = 1.
#  endif

c====================================
# elif defined NHA || defined CARB
c====================================

      weirun = 0.5
      sealev = 0.

c====
# else
c====

      weirun = 0.5
c     sea level unchanged from initial value     

c=====
# endif
c=====

      return
      end

c----------------------------------------------------------------------

c!!!!!!!!!!!!!!!!!!!!
# if defined EISANTA
c!!!!!!!!!!!!!!!!!!!!

      subroutine specmap (sealev, timeice, iu)

c        If first call, read SPECMAP d18O record for last 782 ka, 
c        convert to sea level.
c        Interpolate to model time (timeice), return sealev.

      character*200 cfilin
      parameter (nspec=782)
      dimension timespec(0:nspec), o18spec(0:nspec), sealevspec(0:nspec)
      logical firstspec
      data firstspec /.true./
      save timespec, o18spec, sealevspec, firstspec

#  if defined SUN || defined CRAY || defined IGUANA
      cfilin = '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Data_specmap/'
     *      // 'specmap.017'
#  elif defined LION
      cfilin = '/home2/dxp21/Sheetshelf/Data_specmap/'
     *      // 'specmap.017'
#  endif

c        If first call, read data file
      
c------------------------
      if (firstspec) then
c------------------------
        write (6,'(a)') 'Reading Specmap file'
        open (iu, file=cfilin, status='old')
        read (iu,'(/)') 
        do i=0,nspec
          read (iu,*,err=200,end=200) timespec(i), o18spec(i)
        enddo
        close (iu)

c          Convert to sea level, and convert time to (-) years BP  

        o18mod = o18spec(0)
        o18lgm = o18spec(19)
        do i=0,nspec
          zweil = (o18spec(i) - o18mod) / (o18lgm-o18mod) 
          sealevspec(i) = 0.*(1.-zweil) - 125.*zweil
          timespec(i) = -1000.*timespec(i)
c         write (7,'(i6,2f20.2)') i, timespec(i), sealevspec(i)
        enddo

        firstspec = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.ge.timespec(0)) then
        sealev = sealevspec(0)
        return
      else if (timeice.le.timespec(nspec)) then
        sealev = sealevspec(nspec)
        return
      else
        do i=0,nspec-1
          if (timeice.le.timespec(i) .and. timeice.ge.timespec(i+1))then
            zwei = (timeice-timespec(i+1)) / (timespec(i)-timespec(i+1))
            sealev = zwei*sealevspec(i) + (1.-zwei)*sealevspec(i+1)
            return
          endif
        enddo
      endif
      write (6,'(a)') "Error in specmap - shouldn't get here"
      stop

  200 write (6,'(a,i6)') 'Error 200 in reading Specmap input file: i=',i
      stop

      end

c-----------------------------------------------------------------------

      subroutine bassinot (sealev, timeice, iu)

c        If first call, read Bassinot(1994) d18O record for last 890 ka,
c        convert to sea level.
c        Interpolate to model time (timeice), return sealev.

      character*200 cfilin
      parameter (nbas=443)
      dimension timebas(nbas), o18bas(nbas), sealevbas(nbas)
      logical firstbas
      data firstbas /.true./
      save timebas, o18bas, sealevbas, firstbas

#  if defined SUN || defined CRAY || defined IGUANA
      cfilin = '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Data_bassinot/'
     *      // 'bassinot.final.dat'
#  elif defined LION
      cfilin = '/home2/dxp21/Sheetshelf/Data_bassinot/'
     *      // 'bassinot.final.dat'
#  endif

c        If first call, read data file
      
c------------------------
      if (firstbas) then
c------------------------
        write (6,'(a)') 'Reading Bassinot file'
        open (iu, file=cfilin, status='old')
        do i=1,nbas
          read (iu,*,err=200,end=200) timebas(i), o18bas(i)
        enddo
        close (iu)

c          Convert to sea level, and convert time to (-) years BP  

        o18mod = o18bas(1)
        o18lgm = o18bas(6)
        do i=1,nbas
          zweil = (o18bas(i) - o18mod) / (o18lgm-o18mod) 
          sealevbas(i) = 0.*(1.-zweil) - 125.*zweil
          timebas(i) = -1000.*timebas(i)
c         write (7,'(i6,2f20.2)') i, timebas(i), sealevbas(i)
        enddo

        firstbas = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.ge.timebas(1)) then
        sealev = sealevbas(1)
        return
      else if (timeice.le.timebas(nbas)) then
        sealev = sealevbas(nbas)
        return
      else
        do i=1,nbas-1
          if (timeice.le.timebas(i) .and. timeice.ge.timebas(i+1))then
            zwei = (timeice-timebas(i+1)) / (timebas(i)-timebas(i+1))
            sealev = zwei*sealevbas(i) + (1.-zwei)*sealevbas(i+1)
            return
          endif
        enddo
      endif
      write (6,'(a)') "Error in basssinot - shouldn't get here"
      stop

  200 write(6,'(a,i6)') 'Error 200 in reading Bassinot input file: i=',i
      stop

      end

c-----------------------------------------------------------------------

      subroutine epica (dtantann, timeice, iu)

c        If first call, read EPICA (Dome C) dD record for last ~800 ka, 
c        bin into 1000 yr intervals. 
c        Interpolate to model time (timeice), return dtantann.

      character*200 cfilin
      parameter (nepin=5788)
      dimension timein(nepin), dtin(nepin)

      parameter (nepmax=1600)
      dimension timeepica(nepmax), dtepica(nepmax)

      logical firstepica
      data firstepica /.true./
      save timeepica, dtepica, nepica, firstepica

#  if defined SUN || defined CRAY || defined IGUANA
      cfilin = '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Data_epica/'
     *      // 'edc3deuttemp2007.txt'
#  elif defined LION
      cfilin = '/home2/dxp21/Sheetshelf/Data_epica/'
     *      // 'edc3deuttemp2007.txt'
#  endif

c        If first call, read data file, bin into uniform intervals
      
c-------------------------
      if (firstepica) then
c-------------------------
        write (6,'(a)') 'Reading Epica file'
        open (iu, file=cfilin, status='old')
        read (iu,'()') 
        do i=1,nepin
          read (iu,*,err=200,end=200) 
     *      ibag, skip1, timein(i), dum2, dtin(i)
          if (ibag.ne.i+12) then
            write (6,'(2a,2i6)') 'Error reading Epica input file.',
     *                          ' ibag, i =', ibag, i
            stop
          endif
        enddo
        close (iu)

c          Bin into 500 yr intervals, make timeepica -ve (yrs BP).
c          Also shift dt's so most recent (dtepica(1)) = 0. 

        ie = 0
        do m=1,1600
          timea = (m-1.)*500.
          timeb = m*500.

          dtbin = 0.
          nbin = 0  
          do i=1,nepin 
            if (timein(i).gt.timea .and. timein(i).le.timeb) then
              dtbin = dtbin + dtin(i)
              nbin = nbin + 1
            endif
          enddo
          if (nbin.gt.0) then
            ie = ie + 1
            timeepica(ie) = -0.5*(timea+timeb) 
            dtepica(ie) = dtbin/nbin
          endif
        enddo

        nepica = ie

        do i=1,nepica
          dtepica(i) = dtepica(i) - dtepica(1)
c         write (7,'(i6,2f20.2)') i, timeepica(i), dtepica(i)
        enddo

        firstepica = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.ge.timeepica(1)) then
        dtantann = dtepica(1)
        return
      else if (timeice.le.timeepica(nepica)) then
        dtantann = dtepica(nepica)
        return
      else
        do i=1,nepica
          if (timeice.le.timeepica(i) .and. timeice.ge.timeepica(i+1))
     *      then
            zwei = (timeice     -timeepica(i+1))
     *           / (timeepica(i)-timeepica(i+1))
            dtantann = zwei*dtepica(i) + (1.-zwei)*dtepica(i+1)
            return
          endif
        enddo
      endif
      write (6,'(a)') "Error in Epica - shouldn't get here"
      stop

  200 write (6,'(a,i6)') 'Error 200 in reading Epica input file: i=',i
      stop

      end

c-----------------------------------------------------------------------

      subroutine lisiecki (sealev, timeice, iu)

c        If first call, read Lisiecki+Raymo d18O record for last 5.32Ma,
c        convert to sea level.
c        Interpolate to model time (timeice), return sealev.

      character*200 cfilin
      parameter (nlis=2114)
      dimension timelis(0:nlis), o18lis(0:nlis), sealevlis(0:nlis)
      logical firstlis
      data firstlis /.true./
      save timelis, o18lis, sealevlis, firstlis

#  if defined SUN || defined CRAY || defined IGUANA
      cfilin = '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Data_lisiecki/'
     *      // 'LR04stack.txt'
#  elif defined LION
      cfilin = '/home2/dxp21/Sheetshelf/Data_lisiecki/'
     *      // 'LR04stack.txt'
#  endif

c        If first call, read data file
      
c------------------------
      if (firstlis) then
c------------------------
        write (6,'(a)') 'Reading Lisiecki+Raymo file'
        open (iu, file=cfilin, status='old')
        read (iu,'(////)') 
        do i=0,nlis
          read (iu,*,err=200,end=200) timelis(i), o18lis(i), skip
        enddo
        close (iu)

c          Convert to sea level, and convert time to (-) years BP  

        o18mod = o18lis(0)
        o18lgm = o18lis(18)
        do i=0,nlis
          zweil = (o18lis(i) - o18mod) / (o18lgm-o18mod) 
          sealevlis(i) = 0.*(1.-zweil) - 125.*zweil
          timelis(i) = -1000.*timelis(i)
c         write (7,'(i6,3f20.2)') i, timelis(i), o18lis(i), sealevlis(i)
        enddo

        firstlis = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.ge.timelis(0)) then
        sealev = sealevlis(0)
        return
      else if (timeice.le.timelis(nlis)) then
        sealev = sealevlis(nlis)
        return
      else
        do i=0,nlis-1
          if (timeice.le.timelis(i) .and. timeice.ge.timelis(i+1))then
            zwei = (timeice-timelis(i+1)) / (timelis(i)-timelis(i+1))
            sealev = zwei*sealevlis(i) + (1.-zwei)*sealevlis(i+1)
            return
          endif
        enddo
      endif
      write (6,'(a)') "Error in lisiecki - shouldn't get here"
      stop

  200 write (6,'(a,i6)') 
     *  'Error 200 in reading Lisiecki+Raymo input file: i=',i
      stop

      end

c-----------------------------------------------------------------------

      subroutine laskar (dtantann, dtantjan, dtseas, ecc, obl, prec, 
     *                   timeice, iu)

c        If first call, read Laskar (2004) file for orbital elements
c        from 51 Ma to present, and calculates annual and January
c        insolation for latitude alatlask.
c        Interpolate to model time (timeice), return dtant[ann,jan] 
c        and dtseas.

      character*200 cfilin
#  if defined FORCEEO
      parameter (nlask=51000) ! Laskar file has last 51 Myr
#  else
      parameter (nlask=2000)  ! just last 2 Myr
#  endif
      dimension timelask(0:nlask), 
     *           ecclask(0:nlask), obllask(0:nlask), preclask(0:nlask)
      logical firstlaskread, firstlask
      data firstlaskread, firstlask /.true., .true./
      save firstlaskread, firstlask, 
     *     timelask, ecclask,  obllask, preclask, 
     *     fluxann0, fluxjan0, fluxjul0
      parameter (pi = 3.14159265358979)
      dimension cosq(1), fraq(1), cosq24(1), fraq24(1), alatq(1)
      parameter (solcon = 1367.)
      parameter (alatlask = -80.)
      parameter (nflux=365)
      parameter (dtdf_ann  = 0.1)         ! d[~sfc temperature]/d[insol]
      parameter (dtdf_seas = 0.1)         ! d[~sfc temperature]/d[insol]

#  if defined SUN || defined CRAY || defined IGUANA
      cfilin = '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Data_laskar/'
     *      // 'INSOLN.LA2004.BTL.ASC'
#  elif defined LION
      cfilin = '/home2/dxp21/Sheetshelf/Data_laskar/'
     *      // 'INSOLN.LA2004.BTL.ASC'
#  elif defined NCOM
      cfilin = '/home/geosc/tbellato/Sheetshelf/Data_laskar/'
     *      // 'INSOLN.LA2004.BTL.ASC'
#  endif

c        If first call, read data file
      
c----------------------------
      if (firstlaskread) then
c----------------------------
        write (6,'(a)') 'Reading Laskar orbit file'
        open (iu, file=cfilin, status='old')

        do i=0,nlask
          read (iu,*,err=200,end=200) 
     *      timelask(i), ecclask(i), obllask(i), preclask(i)

c         1000's yrs to yrs (-ve past):
          timelask(i) = 1.e3*timelask(i)  

c         convert laskar prec (NHAE->peri) to Genesis (peri->NHVE):
          preclask(i) = mod (pi - preclask(i) + 2.*pi, 2.*pi)
        enddo

        close (iu)
        firstlaskread = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.ge.timelask(0)) then
        ecc_cur  = ecclask(0)
        obl_cur  = obllask(0)
        prec_cur = preclask(0)
        goto 100
      else if (timeice.le.timelask(nlask)) then
        ecc_cur  = ecclask(nlask)
        obl_cur  = obllask(nlask)
        prec_cur = preclask(nlask)
        goto 100
      else
        do i=0,nlask-1
          if (timeice.le.timelask(i) .and. timeice.ge.timelask(i+1))then
            zwei = (timeice-timelask(i+1)) / (timelask(i)-timelask(i+1))
            ecc_cur  = zwei*ecclask(i) + (1.-zwei)*ecclask(i+1)
            obl_cur  = zwei*obllask(i) + (1.-zwei)*obllask(i+1)

            if (preclask(i+1)-preclask(i).lt.-pi) then
              zpm = preclask(i)
              zpp = preclask(i+1) + 2.*pi
            else if (preclask(i+1)-preclask(i).gt.pi) then
              zpm = preclask(i) + 2.*pi
              zpp = preclask(i+1)
            else
              zpm = preclask(i)
              zpp = preclask(i+1)
            endif
            if (abs(zpp-zpm).gt.pi) then
              write (6,'(/2a,2f8.2, 3x, 2f8.2, 3x, 2i6)') 
     *           ' Error (laskar):',
     *           ' preclask(i:i+1), zp[m,p], i,i+1 =', 
     *           preclask(i), preclask(i+1), zpm,zpp, i,i+1
              stop
            endif
            prec_cur = zwei*zpm + (1.-zwei)*zpp
            prec_cur = mod (prec_cur + 2.*pi, 2.*pi)

            goto 100
          endif
        enddo
      endif
      write (6,'(a)') "Error in laskar - shouldn't get here"
      stop
  100 continue

c777  ecc_cur = .05 
c777  prec_cur = mod ( 0.5*pi - (2.*pi*timeice/40.e3) + 2.*pi, 2.*pi )
c777  obl_cur = (  0.5*(22.2+24.5) 
c777 *           + 0.5*(24.5-22.2)*cos(2.*pi*timeice/40.e3) ) * pi/180.

c        If first call, do flux calculations (zencal) twice,
c        first time for modern to store flux[ann,jan,jul]0.

      if (firstlask) then
        nloopf = 2
        firstlask = .false.
      else
        nloopf = 1
      endif

c~~~~~~~~~~~~~~~~~~~~~~~
      do iloopf=1,nloopf
c~~~~~~~~~~~~~~~~~~~~~~~

        if (nloopf.eq.2 .and. iloopf.eq.1) then
          ecc  = ecclask(0)
          obl  = obllask(0)
          prec = preclask(0)
          time = 0.
        else
          ecc  = ecc_cur
          obl  = obl_cur
          prec = prec_cur
          time = timeice
        endif

        zs = sin (0.5*prec) / sqrt ((1.+ecc)/(1.-ecc))
        zc = cos (0.5*prec)
        ze = 2. * atan2 (zs,zc)
        vern = ze - ecc * sin(ze)

c          Calculate daily mean insols at latitude alatlask,
c          and save averages for annual, January and July

        fluxann = 0.
        fluxjan = 0.
        fluxjul = 0.
        nann = 0
        njan = 0
        njul = 0

        isecdy = 86400/2
        dtq    = 86400.
        alatq(1) = alatlask*pi/180.
        do m=1,nflux
          isecyr = nint (86400.*365.*(m-.5)/nflux)
          call zencal (nint(time), isecyr, isecdy, dtq,
     *                 ecc, obl, prec, 
     *                 vern, dist, eccf,
     *                 cosq, fraq, cosq24, fraq24, alatq, 1)
          zflux = cosq(1)*fraq(1)*eccf*solcon
          fluxann = fluxann + zflux 
          nann = nann + 1

c         either calendar January:
c         if (isecyr.lt.86400*31) then  
c         or Dec 21 to Jan 20, as in Laskar program:
          if (isecyr.lt.86400*20 .or. isecyr.gt.86400*(365-11)) then
            fluxjan = fluxjan + zflux
            njan = njan + 1
          endif
c         either calendar July:
c         if (isecyr.gt.86400*181 .and. isecyr.lt.86400*212) then  
c         or Jun 21 to Jul 20:
          if (isecyr.gt.86400*171 .and. isecyr.gt.86400*202) then
            fluxjul = fluxjul + zflux
            njul = njul + 1
          endif
        enddo
        fluxann = fluxann / nann
        fluxjan = fluxjan / njan
        fluxjul = fluxjul / njul

        if (nloopf.eq.2 .and. iloopf.eq.1) then
          fluxann0 = fluxann
          fluxjan0 = fluxjan
          fluxjul0 = fluxjul
        endif

c~~~~~~~~~~
      enddo
c~~~~~~~~~~

c        Change to (differences from modern insolation)*(dT/dinsol)

      dtantann = ( fluxann - fluxann0 ) * dtdf_ann
      dtantjan = ( fluxjan - fluxjan0 ) * dtdf_seas
      dtseas   = ( fluxjan - fluxjul  ) * dtdf_seas 

c     write (7,'(f10.1, f10.6, 6f10.2)') 
c    *  time*1.e-3, 
c    *  ecc, obl*180./pi, prec*180./pi,
c    *  fluxann, fluxjan, fluxjul, 
c    *  (fluxann-fluxann0)*dtdf_ann, 
c    *  (fluxjan-fluxjan0)*dtdf_seas,
c    *  (fluxjul-fluxjul0)*dtdf_seas

      return

  200 write (6,'(a,i6)') 'Error 200 in reading Laskar input file: i=',i
      stop

      end
      
c!!!!!
# endif
c!!!!!
#else
      subroutine ism

      return
      end
#endif
