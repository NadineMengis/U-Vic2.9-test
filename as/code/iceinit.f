! source file: /net/mare/home1/eby/as/ism/iceinit.F
c-----------------------------------------------------------------------

      subroutine initgrid (firstgrid, xglu)

c     Initializes horizontal and vertical grids

      use comicephys
      use comicegrid

      logical firstgrid

      dimension dzeta10(10)
      data dzeta10
     *  /.05, .10, .15, .15, .20, .15, .10, .05, .03, .02/

c     dimension dzsed3(3)
c     data dzsed3
c    *  /.10, .45, .45/

      dimension dzbed6(6)
      data dzbed6
     *  /10., 30., 100., 300., 600., 1000./
      character(120) g_st
      dimension temp(nx,ny)
c        Set horizontal grids

c===================
c=====================

c*******************************************
c****
      if (isname=='Antarctica') then
        zlambda = 1. + cos((stdparallel+90.)*pi/180.)
        xoffa = 0.
        yoffa = 0.
      elseif (isname=='Greenland') then
        zlambda = 1. + cos((90.-stdparallel)*pi/180.)
        !input grid lower left corner x offset from NP (m):
	!see:http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
	lleftx = -800000.
        !input grid lower left corner y offset from NP (m)
	!see: http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland
	llefty = -3400000.
        xoffa = lleftx+dx0*nx/2
	yoffa = llefty+dy0*ny/2
      else
        stop 'no stereographic ice grid initialized: iceinit'
      endif

      totarea = 0.
      do j=1,ny
        do i=1,nx
          zx = dx0*(i-0.5-0.5*nx) + xoffa
          zy = dy0*(j-0.5-0.5*ny) + yoffa

          zr = sqrt (zx**2 + zy**2) / radius
          zscale = 2. / ( zlambda * (1. + (zr/zlambda)**2) )
          dx(i,j) = dx0 * zscale
          dy(i,j) = dy0 * zscale
          darea(i,j) = dx(i,j)*dy(i,j)
          totarea = totarea + darea(i,j)

          if (isname=='Antarctica') then

            if (zx.eq.0. and. zy.eq.0.) then
              alond(i,j) = 0.
            else
	      alond(i,j) = (0.5*pi - atan2(zy,zx)) * 180./pi           !SP
            endif
            alond(i,j) = mod (alond(i,j) + 720., 360.)
            if (alond(i,j).gt.180.) alond(i,j) = alond(i,j) - 360.

            alatd(i,j) = -acos (2.*(zr/zlambda) / (1.+(zr/zlambda)**2))
     *                   * 180./pi 	  	                       !SP
            temp(i,j) = zscale

          elseif (isname=='Greenland') then

	    alond(i,j) = -(atan(zx/zy))*180./pi -39.
            alatd(i,j) = acos (2.*(zr/zlambda) / (1.+(zr/zlambda)**2))
     *                   * 180./pi 	  	                       !SP
            temp(i,j) = zscale

          endif

        enddo
      enddo
      do j=0,nyp
        do i=0,nxp
          zr = sqrt (   (dx0*(i    -0.5*nx) + xoffa)**2
     *                + (dy0*(j-0.5-0.5*ny) + yoffa)**2 ) / radius
          zscale = 2. / ( zlambda * (1. + (zr/zlambda)**2) )
          dxu(i,j) = dx0 * zscale
          dyu(i,j) = dy0 * zscale

          zr = sqrt (   (dx0*(i-0.5-0.5*nx) + xoffa)**2
     *                + (dy0*(j    -0.5*ny) + yoffa)**2 ) / radius
          zscale = 2. / ( zlambda * (1. + (zr/zlambda)**2) )
          dxv(i,j) = dx0 * zscale
          dyv(i,j) = dy0 * zscale
        enddo
      enddo

      do j=0,ny
        do i=0,nx
          zr = sqrt (   (dx0*(i    -0.5*nx) + xoffa)**2
     *                + (dy0*(j    -0.5*ny) + yoffa)**2 ) / radius
          zscale = 2. / ( zlambda * (1. + (zr/zlambda)**2) )
          dxc(i,j) = dx0 * zscale
          dyc(i,j) = dy0 * zscale
        enddo
      enddo

      do j=1,ny
c       zx = 0.5*dxu((nx+1)/2,j) * (1-mod(nx,2))
        zx = 0.5*dx0             * (1-mod(nx,2))
        do i=nx/2+1,nx
          xh(i,j) = zx + xoffa
          xh(nx+1-i,j) = -zx + xoffa
c         zx = zx + dxu(i-1,j)
          zx = zx + dx0
        enddo
      enddo
      do i=1,nx
c       zy = 0.5*dyv(i,(ny+1)/2) * (1-mod(ny,2))
        zy = 0.5*dy0             * (1-mod(ny,2))
        do j=ny/2+1,ny
          yh(i,j) = zy + yoffa
          yh(i,ny+1-j) = -zy + yoffa
c         zy = zy + dyv(i,j-1)
          zy = zy + dy0
        enddo
      enddo

c*****
c*****

c====
c=====

c        Set vertical ice grid

      if (nlev.eq.10) then
        call scopy (nlev, dzeta10, 1, dzeta, 1)
      else if (nlev.eq.30) then
        do k=1,10
          dzeta(k) = 0.1/10.
        enddo
        do k=11,20
          dzeta(k) = 0.8/10.
        enddo
        do k=21,30
          dzeta(k) = 0.1/10.
        enddo
      else
        do k=1,nlev
          dzeta(k) = 1./nlev
        enddo
c       write (ioterm,*) 'Error: unknown nlev=',nlev
c       stop
      endif

      zeta(0) = 0.
      zeta(1) = 0.5*dzeta(1)
      do k=2,nlev
        zeta(k) = zeta(k-1) + 0.5*(dzeta(k-1)+dzeta(k))
      enddo
      zeta(nlevp) = 1.

      do k=0,nlev
        dzetah(k) = zeta(k+1) - zeta(k)
      enddo

      zetah(0) = 0.
      do k=1,nlev
        zetah(k) = zetah(k-1) + dzeta(k)
      enddo

c        Set vertical sediment grid

c     if (nsed.eq.3) then
c       call scopy (nsed, dzsed3, 1, dzsed, 1)
c     else
        do k=1,nsed
          dzsed(k) = 1./nsed
        enddo
c     endif

      zsedm(0) = 0.
      zsed(0) = 0.
      do k=1,nsed
        zsed(k)  = zsedm(k-1) + 0.5*dzsed(k)
        zsedm(k) = zsedm(k-1) +     dzsed(k)
      enddo
      if (abs(zsedm(nsed)-1.).gt.1.e-8) then
        write(ioterm,*) '*** Error: sum(dzsed) = ',zsedm(nsed)
        stop
      else
        zsedm(nsed) = 1.
        zsed(nsedp) = 1.
      endif

c        Set vertical bedrock grid

      if (nbed.eq.6) then
        call scopy (nbed, dzbed6, 1, dzbed, 1)
      else if (nbed.eq.1) then
        dzbed(1) = 30.  !   10.
      else
        do k=1,nbed
          dzbed(k) = 300.
        enddo
      endif

      zbedm(0) = 0.
      zbed(0) = 0.
      bedthick = 0.
      do k=1,nbed
        zbedm(k) = bedthick
        zbed(k)  = bedthick + dzbed(k)*0.5
        bedthick = bedthick + dzbed(k)
      enddo
      zbedm(nbed) = bedthick
      zbed(nbedp) = bedthick
      return
      end

c-----------------------------------------------------------------------

      subroutine initphys (h, hb, t,
     *                     sedim, sedimeq, sedimun, tsed, wsed,
     *                     topbed, topbedeq, tbed, equiload, geoflux,
     *                     hw, tw, maskpres, sedpres,
     *                     ua, va, sealev, timeice)

c     Initializes prognostic model variables

      use comicephys
      use comicegrid
      use comicebed
c     passed:
      dimension
     *  h(nx,ny),            hb(nx,ny),
     *  t(nx,ny,0:nlevp),
     *  sedim(nx,ny),        sedimeq(nx,ny),     sedimun(nx,ny),
     *  tsed(nx,ny,nsed),    wsed(nx,ny,nsed),
     *  topbed(nx,ny),       topbedeq(nx,ny),    tbed(nx,ny,nbed),
     *  equiload(nx,ny),     geoflux(nx,ny),
     *  hw(nx,ny),           tw(nx,ny),
     *  maskpres(nx,ny),     sedpres(nx,ny),
     *  ua(0:nxp,0:nyp),     va(0:nxp,0:nyp)

      dimension htmp(nx,ny), hbtmp(nx,ny), hwtmp(nx,ny)

      dimension bathy(nx,ny), bedelev(nx,ny), groundbed(nx,ny),
     *          thicice(nx,ny), surface(nx,ny), water(nx,ny)
      character(120) g_st,fname
      integer iou,ntrec,ib(10),ic(10),ln
      real tsur(nx,ny)
      logical exists
c set sealev just for initial settings below, may be changed in first
c timestep of run (icectl)
      sealev = 0.

      call izero (maskpres, nx*ny)
      call zero  (sedpres, nx*ny)

c        Ice thickness, sed thickness, elevations:

c===================================
c=====================

      if (isname == 'Antarctica') then

c        Read Bedmap files, aggregate to model grid
c        (assumes all have same nulval).

        call bedmapaggreg ('bathy',     bathy,     nulval)
        call bedmapaggreg ('bedelev',   bedelev,   nulval)
        call bedmapaggreg ('groundbed', groundbed, nulval)
        call bedmapaggreg ('icethic',   thicice,   nulval)
        call bedmapaggreg ('surface',   surface,   nulval)
        call bedmapaggreg ('water',     water,     nulval)

        do j=1,ny
          do i=1,nx

c         set bed elevations everywhere:
            if (nint(bedelev(i,j)).eq.nulval) then
              topbed(i,j) = -4000.                        ! lats > -60
            else
              topbed(i,j) = bedelev(i,j)
            endif
	    !Fill in deep hole around Berkner Island and offshore, E. Ant. Pen. and under Amery to avoid crazy melt there
	   if (i.ge.95.and.i.le.110.and.j.ge.158.and.j.le.174)
     & 	   topbed(i,j) = max(topbed(i,j),-1000.)
	   if (i.ge.67.and.i.le.84.and.j.ge.145.and.j.le.163)
     & 	   topbed(i,j) = max(topbed(i,j),-1000.)
	   if (i.ge.222.and.i.le.233.and.j.ge.172.and.j.le.183)
     & 	   topbed(i,j) = max(topbed(i,j),-1000.)

c           set ice thicknesses everywhere:
            if (nint(thicice(i,j)).eq.nulval) then
              h(i,j) = 0.
            else
              h(i,j) = thicice(i,j)
            endif

c         set topbedeq, allowing for current local equilibrated
c         isostatic depression by grounded ice:
            if ( alatd(i,j).gt.-60.5        .or.
     *           nint(bathy(i,j)).ne.nulval .or.
     *           nint(water(i,j)).ne.nulval ) then
c           maskpres is modern ice shelf distrib, used in ocean:
              if (h(i,j).gt.0.) then
                maskpres(i,j) = 2
              else
                maskpres(i,j) = 1
              endif
c             bathymetry or floating ice:

              topbed(i,j) = min (sealev-2., topbed(i,j))  ! ocn bed < sl
              topbedeq(i,j) = topbed(i,j)
            else
              maskpres(i,j) = 0
              topbedeq(i,j) = topbed(i,j) + (rhoice/rhobed)*h(i,j)
            endif

c         no ice, start from rebounded equil (comment out for modobs):
            sedimeq(i,j) = 0.

          enddo
        enddo

      elseif (isname == 'Greenland') then
         fname='topg'
         call bamberaggreg (fname, topbed)
         fname='thk'
         call bamberaggreg (fname, h)

	 !all ocean
	 maskpres(:,:) = 1
	 !mask grounded ice (note: 2= ice shelves, none exist)
	 where (h .gt. 0.) maskpres = 0
         topbedeq(:,:) = topbed(:,:) + (rhoice/rhobed)*h(:,:)
         sedimeq(:,:) = 0.

      endif
!       g_st = trim(isname)//'topbed'
!       call jerncdf_ice_snapshot(topbed,g_st)

c=======================
c=====

c        Set sedim,sedimun,hb

      do j=1,ny
        do i=1,nx
          sedim(i,j)   = sedimeq(i,j)
          sedimun(i,j) = sedim(i,j)
          hb(i,j) = topbed(i,j) + sedim(i,j)
        enddo
      enddo

c        Set equilibrium (no-ice) load for subr bedrock, using hydwater
c        to find location of ocean and lakes for that (temporary)
c        topography.
c        Also set sedpres is for sedim distrib, used in basewater.

      do j=1,ny
        do i=1,nx
          htmp(i,j) = 0.
          hbtmp(i,j) = topbedeq(i,j) + sedimeq(i,j)
        enddo
      enddo
      call initwater (htmp, hbtmp, hwtmp, sealev)
      do j=1,ny
        do i=1,nx
          equiload(i,j) = rholiq*hwtmp(i,j) + rhosed*sedimeq(i,j)
          if (hwtmp(i,j).gt.0.) then
            sedpres(i,j) = 1.
          else
            sedpres(i,j) = 0.
          endif
        enddo
      enddo

c        Use ocean/lake search algorithm initwater
c        to set initial water hw (no sub-ice lakes)

      call initwater (h, hb, hw, sealev)

c        Ice, sed, bed temperatures, sed liquid content, tracers:

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      fname = 'data/'//trim(isname)//'tsurf.nc'
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, fname,' doesnt exist, kid.'
	stop '-> iceinit.F'
      elseif (exists) then
        call openfile (trim(fname), iou)
        ntrec=1
        tsur(:,:)=0.
        ib(1) = 1
        ic(1) = nx
        ib(2) = 1
        ic(2) = ny
        ln = ic(1)*ic(2)
        call getvara ('output', iou, ln, ib, ic, tsur, 1., 0.)
        !call closefile(iou)
      endif
c----------------
      do j=1,ny
        do i=1,nx
c----------------

c            Initialize ice and sediment temperatures
c            where significant ice or sed, assuming vertical linear
c            diffusive profiles and basal ice at pressure melting point

c         leave ice temps unchanged (from readnest) if nested run:
          dtdz = geoflux(i,j) / (2.1*31556926)        !conduc ~ in vdif2
          do k=0,nlevp
	  !Jer: set internal ice temperatures to linear mix between surface temperature data at surface
	  !and 5C below pressure melting point at base.
	  ttop=tsur(i,j)
	  if (h(i,j) .gt. 0.) then
	    ttop=tsur(i,j)
	    tbot=tmelt-dtmdh*h(i,j)-5.
	    t(i,j,k)=ttop*(1.-zeta(k))+tbot*zeta(k)
	  else
	    t(i,j,k)=tmelt
	  endif

          enddo

c            Sediment temperatures and liquid content of pores (wsed,
c            0-1, rest is ice, assumed saturated, set in vdif):

          do k=1,nsed
            if (sedim(i,j).gt.0.01) then
              tsed(i,j,k) = tmelt-dtmdh*h(i,j)
              wsed(i,j,k) = 1.
            else
              tsed(i,j,k) = tmelt
              wsed(i,j,k) = 0.
            endif
          enddo

c            Bedrock temperatures:

          dtdz = geoflux(i,j) / (3.3*31556926)       !conduc as in vdif2
          do k=1,nbed
            tbed(i,j,k) = tmelt - dtmdh*h(i,j) + dtdz*zbed(k)
          enddo

c            Tracers (1st tracer = temperature)

!Jer: comment out addition tracer initializations
!           do n=ntraca+1,ntrace
!             do k=0,nlevp
!               t(i,j,k) = 0.
!             enddo
!           enddo

c           Water temperature (set to pressure melt point)

          tw(i,j) = tmelt - dtmdh*h(i,j)

c------------
        enddo
      enddo
c------------

      return
      end

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c)))))

c-----------------------------------------------------------------------

c))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c)))))

c-----------------------------------------------------------------------

c((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
c((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      subroutine bedmapaggreg (cname, arro, nulval)

c     Reads Bedmap data file cname.asc, aggregates to model resolution,
c     returns field in arro. Also returns "no-data" value nulval.

      use comicegrid

      character*(*) cname
      dimension arro(nx,ny)

      parameter (nximax=1500, nyimax=1500)
      dimension arri(nximax,nyimax)

      dimension narro(nx,ny), narroall(nx,ny)
      character*1 carro(nx,ny)

      character*80 cfilout
      nxo=nx
      nyo=ny
c        Read Bedmap input file into local array arri

      call bedmapread (cname, arri, nximax, nyimax,
     *                 nxi, nyi, ddi, nulval, iubed)

c        Aggregate to output resol. Loop over Bedmap points, calculate
c        the model grid box it's in (via zx,zy, distances from SP),
c        assign mean value of all those points if there's >= 1/2 good
c        data points out of all in the box.

      do jo=1,nyo
        do io=1,nxo
          arro(io,jo) = 0.
          narro(io,jo) = 0
          narroall(io,jo) = 0
        enddo
      enddo

c        zx[nr], zy[nr] are +/- distances from origin (SP):
c        rotate input 90 deg clockwise: so set zx to zynr, zy to -zxnr

      do j=1,nyi
        do i=1,nxi
          zxnr = ddi*(i-0.5-0.5*nxi)
          zynr = ddi*(j-0.5-0.5*nyi)

          zx =  zynr
          zy = -zxnr

          io = nint(zx/dd0 + 0.5*(1-mod(nxo,2))) + (nxo+1)/2
          jo = nint(zy/dd0 + 0.5*(1-mod(nyo,2))) + (nyo+1)/2

          if (io.ge.1.and.io.le.nxo .and. jo.ge.1.and.jo.le.nyo) then
            if (nint(arri(i,j)).ne.nulval) then
              arro(io,jo) = arro(io,jo) + arri(i,j)
              narro(io,jo) = narro(io,jo) + 1
            else
              narroall(io,jo) = narroall(io,jo) + 1
            endif
          endif
        enddo
      enddo

      do jo=1,nyo
        do io=1,nxo
          if (narro(io,jo).eq.0) then
            arro(io,jo) = nulval
          else if (narro(io,jo).ge.(narroall(io,jo)+1)/2) then
            arro(io,jo) = arro(io,jo)/narro(io,jo)
          else
            arro(io,jo) = nulval
          endif
        enddo
      enddo

c        Write scratch output file

c----------------------------
      if (iubedout.ne.0) then
c----------------------------

      zscale = 100.
      if (cname.eq.'water') zscale = 200.

      do jo=1,nyo
        do io=1,nxo
          if (nint(arro(io,jo)).eq.nulval) then
            carro(io,jo) = ' '
          else
            zv = arro(io,jo)/zscale
            izv = nint(zv)
            if (zv.eq.0.) then
              carro(io,jo) = '.'

            else if (izv.eq.0 .and. zv.gt.0.) then
              carro(io,jo) = '+'
            else if (izv.ge.1 .and. izv.le.26) then
              carro(io,jo) = char (ichar('A') + izv-1)
            else if (izv.gt.26) then
              carro(io,jo) = '*'

            else if (izv.eq.0 .and. zv.lt.0.) then
              carro(io,jo) = '-'
            else if (izv.le.-1 .and. izv.ge.-26) then
              carro(io,jo) = char (ichar('a') - izv-1)
            else if (izv.lt.-26) then
              carro(io,jo) = '#'
            endif
          endif
        enddo
      enddo

      cfilout = cname(1:lenchr(cname)) // '.out'
      open (iubedout, file=cfilout, status='unknown')
      close (iubedout, status='delete')
      open (iubedout, file=cfilout, status='new', form='formatted')

      write (iubedout,'(a)') cname
      write (iubedout,'(6x,200i5)') (io,io=5,nxo,5)
      do jo=nyo,1,-1
        write (iubedout,'(i4,2x,1000a1)') jo, (carro(io,jo),io=1,nxo)
      enddo

      close (iubedout)

c----------
      endif
c----------

      return
      end

      subroutine bamberaggreg (cname, arro)

c     Reads Bamber-sourced data, aggregates to model resolution,
c     returns field in arro.

      use comicegrid

      character(120) cname,fname
      dimension arro(nx,ny)
      parameter (nxi=300, nyi=560, ddi=5000.)
      dimension arri(nxi,nyi)
      dimension narro(nx,ny)
      integer ib(10), ic(10), iou
      logical inqvardef, exists

      nxo=nx
      nyo=ny
c        Read Bamber input file into local array arri
      fname = 'data/Greenland_5km_v6.nc'
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, trim(fname),' doesnt exist, kid.'
	stop '-> iceinit.F'
      elseif (exists) then
        call openfile (trim(fname), iou)
        exists = inqvardef(trim(cname), iou)
        ib(:) = 1
        ic(:) = 1
        ic(1) = nxi
        ic(2) = nyi
        if (.not. exists) then
          print*, 'Ice init field', trim(cname), ' doesnt exist.'
	  stop '-> iceinit.F'
        elseif (exists) then
	  call getvara (trim(cname), iou, nxi*nyi, ib, ic, arri, 1., 0.)
        endif
      endif

c        Aggregate to output resol. Loop over Bedmap points, calculate
c        the model grid box it's in (via zx,zy, distances from SP),
c        assign mean value of all those points if there's >= 1/2 good
c        data points out of all in the box.

      do jo=1,ny
        do io=1,nx
          arro(io,jo) = 0.
          narro(io,jo) = 0
        enddo
      enddo

      do j=1,nyi
        do i=1,nxi
          zx = ddi*(i-0.5)
          zy = ddi*(j-0.5)

          io = nint(zx/dd0+0.5)
          jo = nint(zy/dd0+0.5)

          if (io.ge.1.and.io.le.nx .and. jo.ge.1.and.jo.le.ny) then
              arro(io,jo) = arro(io,jo) + arri(i,j)
              narro(io,jo) = narro(io,jo) + 1
          endif
        enddo
      enddo
      do jo=1,ny
        do io=1,nx
          if (narro(io,jo).gt.0) then
            arro(io,jo) = arro(io,jo)/real(narro(io,jo))
	  endif
        enddo
      enddo

      return
      end
c(((((
c(((((

c-----------------------------------------------------------------------

      subroutine bedmapread (cnamin, arri, nximax, nyimax,
     *                       nxi, nyi, ddi, nulval, iu)

c     Basic read of one Bedmap data file (cnamin.asc) into temp array
c     arri, and also reads and returns nxi,nyi, ddi, nulval

      character*(*) cnamin
      dimension arri(nximax,nyimax)

      character*20 cname, cdum
      character*80 cfilin

      cname = cnamin

      cfilin = './data/' //
     *      cname(1:lenchr(cname)) // '.asc'

      open (iu, file=cfilin, status='old', form='formatted', err=500)

c        Read header lines, input resol, etc

      read (iu,*) cdum, nyi
      read (iu,*) cdum, nxi
      read (iu,*)
      read (iu,*)
      read (iu,*) cdum,ddi
      read (iu,*) cdum,nulval

      write (6,'(2a)')
     *  ' bedmapread: reading ', cfilin(1:lenchr(cfilin))
      write (6,'(a,2i6,f10.3,i6)')
     *  '   nxi, nyi, ddi, nulval =', nxi,nyi,ddi,nulval

c        Read Bedmap data

      do i=1,nxi
        read (iu,*) (arri(i,j),j=1,nyi)
      enddo

      close (iu)

      return

  500 write (6,'(/2a)') '*** Error (bedmapread) opening file ',
     *  cfilin(1:lenchr(cfilin))
      stop
      end

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

      subroutine readres (h, hs, hb, t,
     *                    baseperc, basefrml, heatb,
     *                    sedim, sedimeq, sedimun, tsed, wsed, heats,
     *                    topbed, topbedeq, tbed, equiload,
     *                    hw, tw, maskpres, sedpres,
     *                    ua, va, ub, vb, dfu, dfv,
     *                    sealev, rco2inter, nwrit, timeicein)

c     Read restart file

      use comicephys
      use comicegrid
      use comicesparse
      use subgrid, only: rtag

c     h grid:
      dimension
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny),
     *  t(nx,ny,0:nlevp),
     *  baseperc(nx,ny),     basefrml(nx,ny),    heatb(nx,ny),
     *  sedim(nx,ny),        sedimeq(nx,ny),     sedimun(nx,ny),
     *  tsed(nx,ny,nsed),    wsed(nx,ny,nsed),
     *  heats(nx,ny,nsed),
     *  topbed(nx,ny),       topbedeq(nx,ny),    tbed(nx,ny,nbed),
     *  equiload(nx,ny),
     *  hw(nx,ny),           tw(nx,ny),
     *  maskpres(nx,ny),     sedpres(nx,ny)

        dimension tracer(nx,ny,0:nlevp,1)

c     u grid:
      dimension
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),       dfv(0:nxp,0:nyp)

      include "pconst.h"
      integer iou,ntrec,ib(10),ic(10)
      real tmp
      real, allocatable :: tmpij(:,:),tmpijk(:,:,:)

      character(120) fname,g_st
      logical exists

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      fname = 'data/'//trim(isname)//'_restart.nc'
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, fname,' doesnt exist, kid.'
	stop '-> iceinit.F'
      endif
      write (ioterm,'(/a,a)')
     *  'Reading restart file ',fname(1:lenchr(fname))
      call openfile (fname, iou)
      ntrec = 1

!-----------------------------------------------------------------------
!     read 1d data
!-----------------------------------------------------------------------
      !tmp=versresin
      call getvars ('versres', iou, ntrec, tmp, c1, c0)
      versresin=tmp
      call getvars ('timeice', iou, ntrec, tmp, c1, c0)
      timeicein=tmp
      call getvars ('nx', iou, ntrec, tmp, c1, c0)
      nxin=tmp
      call getvars ('ny', iou, ntrec, tmp, c1, c0)
      nyin=tmp
      call getvars ('nlev', iou, ntrec, tmp, c1, c0)
      nlevin=tmp
      call getvars ('nsed', iou, ntrec, tmp, c1, c0)
      nsedin=tmp
      call getvars ('nbed', iou, ntrec, tmp, c1, c0)
      nbedin=tmp
      call getvars ('nwrit', iou, ntrec, tmp, c1, c0)
      nwrit=tmp
      call getvars ('sealev', iou, ntrec, tmp, c1, c0)
      sealev=tmp
      call getvars ('rtag', iou, ntrec, tmp, c1, c0)
      if (tmp .eq. rtag) then
        print*, 'Ice sheet and climate appear to be consistent.'
      else
        print*, 'Error: ice restart not consistent with climate.'
 	print*, 'ice tag=',tmp,'climate tag=',rtag
        stop '->subroutine readres (iceinit.F)'
      endif
      write (ioterm,'(a, i10)') '  timeicein  =', nint(timeicein)
      write (ioterm,'(a,2i8)') '  nxin, nyin =', nxin, nyin
      write (ioterm,'(a,2i8)') '  nx  , ny   =', nx,   ny
      write (ioterm,'(a,f8.2)') '   versres =', versresin
      if (       versresin.ne.1.
     *     .and. versresin.ne.2.
     *     .and. versresin.ne.3.
     *   ) then
        write(ioterm,*)
     *    '*** Error: restart file version # must be 1, 2 or 3'
        stop
      endif
      if (nxin.ne.nx .or. nyin.ne.ny) then
         write (ioterm,'(/a/a,2i6/a,2i6)')
     *      ' *** Error (readres): input ice grid not consistent.'
        stop
      endif

      if ( nlev.ne.nlevin .or. nsed.ne.nsedin .or. nbed.ne.nbedin ) then
        write (ioterm,'(a,4(/a,i4,a,i4))')
     *                              '*** Error: dimension mismatch:',
     *                              '    nlev=',nlev,'  nlevin=',nlevin,
     *                              '    nsed=',nsed,'  nsedin=',nsedin,
     *                              '    nbed=',nbed,'  nbedin=',nbedin
        stop
      endif

!-----------------------------------------------------------------------
!     read 2d data (x,y)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = nx
      ib(2) = 1
      ic(2) = ny
      ib(3) = 1
      ic(3) = ntrec
      ln = ic(1)*ic(2)
      allocate (tmpij(nx,ny))
      call getvara ('h', iou, ln, ib, ic, tmpij, c1, c0)
      h = tmpij
      call getvara ('hs', iou, ln, ib, ic, tmpij, c1, c0)
      hs = tmpij
      call getvara ('hb', iou, ln, ib, ic, tmpij, c1, c0)
      hb = tmpij
      call getvara ('baseperc', iou, ln, ib, ic, tmpij, c1, c0)
      baseperc = tmpij
      call getvara ('basefrml', iou, ln, ib, ic, tmpij, c1, c0)
      basefrml = tmpij
      call getvara ('heatb', iou, ln, ib, ic, tmpij, c1, c0)
      heatb = tmpij
      call getvara ('sedim', iou, ln, ib, ic, tmpij, c1, c0)
      sedim = tmpij
      call getvara ('sedimeq', iou, ln, ib, ic, tmpij, c1, c0)
      sedimeq = tmpij
      call getvara ('sedimun', iou, ln, ib, ic, tmpij, c1, c0)
      sedimun = tmpij
      call getvara ('topbed', iou, ln, ib, ic, tmpij, c1, c0)
      topbed = tmpij
      call getvara ('topbedeq', iou, ln, ib, ic, tmpij, c1, c0)
      topbedeq = tmpij
      call getvara ('equiload', iou, ln, ib, ic, tmpij, c1, c0)
      equiload = tmpij
      call getvara ('hw', iou, ln, ib, ic, tmpij, c1, c0)
      hw = tmpij
      call getvara ('tw', iou, ln, ib, ic, tmpij, c1, c0)
      tw = tmpij
      call getvara ('maskpres', iou, ln, ib, ic, tmpij, c1, c0)
      maskpres = tmpij
      call getvara ('sedpres', iou, ln, ib, ic, tmpij, c1, c0)
      sedpres = tmpij
       call getvara ('is2sg', iou, ln, ib, ic, tmpij, c1, c0)
       is2sg = nint(tmpij)
       call getvara ('iceiind', iou, ln, ib, ic, tmpij, c1, c0)
       iceiind = nint(tmpij)
       call getvara ('icejind', iou, ln, ib, ic, tmpij, c1, c0)
       icejind = nint(tmpij)

      deallocate(tmpij)

      ib(1) = 1
      ic(1) = nxp
      ib(2) = 1
      ic(2) = nyp
      ib(3) = 1
      ic(3) = ntrec
      ln = ic(1)*ic(2)
      allocate (tmpij(nxp+1,nyp+1))
      call getvara ('ua', iou, ln, ib, ic, tmpij, c1, c0)
      where (abs(tmpij) .gt. 1.e30) tmpij = 0.
      ua=tmpij
      call getvara ('va', iou, ln, ib, ic, tmpij, c1, c0)
      where (abs(tmpij) .gt. 1.e30) tmpij = 0.
      va=tmpij
      call getvara ('ub', iou, ln, ib, ic, tmpij, c1, c0)
      where (abs(tmpij) .gt. 1.e30) tmpij = 0.
      ub=tmpij
      call getvara ('vb', iou, ln, ib, ic, tmpij, c1, c0)
      where (abs(tmpij) .gt. 1.e30) tmpij = 0.
      vb=tmpij
      call getvara ('dfu', iou, ln, ib, ic, tmpij, c1, c0)
      where (abs(tmpij) .gt. 1.e30) tmpij = 0.
      dfu=tmpij
      call getvara ('dfv', iou, ln, ib, ic, tmpij, c1, c0)
      where (abs(tmpij) .gt. 1.e30) tmpij = 0.
      dfv=tmpij
      deallocate(tmpij)

!-----------------------------------------------------------------------
!     read 3d data (nx,ny,nlev)
!-----------------------------------------------------------------------

      ib(1) = 1
      ic(1) = nx
      ib(2) = 1
      ic(2) = ny
      ib(3) = 1
      ic(3) = nlevp+1
      ib(4) = 1
      ic(4) = ntrec
      ln = ic(1)*ic(2)*ic(3)
      allocate(tmpijk(nx,ny,nlevp+1))
      call getvara ('t', iou, ln, ib, ic, tmpijk, c1, c0)
      tracer(:,:,:,1) = tmpijk
      t(:,:,:) = tracer(:,:,:,1)
      deallocate(tmpijk)

      ib(1) = 1
      ic(1) = nx
      ib(2) = 1
      ic(2) = ny
      ib(3) = 1
      ic(3) = nbed
      ib(4) = 1
      ic(4) = ntrec

      ln = ic(1)*ic(2)*ic(3)
      allocate(tmpijk(nx,ny,nbed))
      tmpijk = tbed
      call getvara ('tbed', iou, ln, ib, ic, tmpijk, c1, c0)
      tbed = tmpijk
      deallocate(tmpijk)

      ib(1) = 1
      ic(1) = nx
      ib(2) = 1
      ic(2) = ny
      ib(3) = 1
      ic(3) = nsed
      ib(4) = 1
      ic(4) = ntrec

      ln = ic(1)*ic(2)*ic(3)
      allocate(tmpijk(nx,ny,nsed))
      call getvara ('tsed', iou, ln, ib, ic, tmpijk, c1, c0)
      tsed = tmpijk
      call getvara ('wsed', iou, ln, ib, ic, tmpijk, c1, c0)
      wsed = tmpijk
      call getvara ('heats', iou, ln, ib, ic, tmpijk, c1, c0)
      heats = tmpijk
      deallocate(tmpijk)

      !if (nint(versresin).ge.3) then
      !  read (iures,'(e25.15)') rco2inter
      !else
      !  rco2inter = 1.
      !endif

c     in common (used in sedtracking_eul):
      !call readin  (sedtrack, 1,nx,1,ny,   ntr+1,         1, 1, 0)
      !read (iures,'(e25.15)') timebot
      !call readini (itrtop,   1,nx,1,ny,   1,             1, 0, 0)
      !read (iures,'(e25.15)') versresend
      !close (iures)

      !if (abs(versresin-versresend).gt.1.e-6) then
      !  write (ioterm,'(a,f8.2)') '   versresend=', versresend
      !  write (ioterm,'(a)')'*** Error: versresin and versresend differ'
      !  stop
      !endif

      !call closefile (iou)

      return
      end

c-----------------------------------------------------------------------

      subroutine readin (arr, nx1, nx2, ny1, ny2, nvert,
     *                   ifhgrid, ifreal, iffirst)

       use comicegrid

      entry     readini (iarr,nx1, nx2, ny1, ny2, nvert,
     *                   ifhgrid, ifreal, iffirst)

c     Read in one field (all horiz slices) from restart file.
c     If different resolution, process one horiz slice at a time:
c       if ifhgrid=1, call interph(i) to transfer to current h-grid.
c          if (ifreal=1), bilinear interp for arr,
c          else nearest neighbor for iarr.
c       else set to zero (u,v grids).

           dimension  arr(nx1:nx2, ny1:ny2, nvert),
     *          iarr(nx1:nx2, ny1:ny2, nvert)

      dimension
     *  work(ninmx,ninmy), iwork(ninmx,ninmy)

c        If same resolution, just read into arr or iarr, return
c        (nxin,nyin, xhin,yhin are in common, have been read from
c        restart file in readsres)

      if (nxin.eq.nx .and. nyin.eq.ny) then
        do k=1,nvert
          if (ifreal.eq.1) then
            read (iures,'(e25.15)') (( arr(i,j,k),i=nx1,nx2), j=ny1,ny2)
          else
            read (iures,'(i10)') ((iarr(i,j,k),i=nx1,nx2), j=ny1,ny2)
          endif
        enddo
        return
      endif
      if (ifhgrid.eq.1) then
        nxin1 = 1
        nxin2 = nxin
        nyin1 = 1
        nyin2 = nyin
      else
        nxin1 = 0
        nxin2 = nxin+1
        nyin1 = 0
        nyin2 = nyin+1
      endif

c        Read in each horizontal slice for this field from
c        restart file into (i)work, transfer to (i)arr

c-----------------
      do k=1,nvert
c-----------------
        if (ifreal.eq.1) then
          read (iures) (( work(i,j),i=nxin1,nxin2), j=nyin1,nyin2)
        else
          read (iures) ((iwork(i,j),i=nxin1,nxin2), j=nyin1,nyin2)
        endif

c          Interpolate each slice to current grid
c          (xhin,yhin,nxin,nyin  are in common)

        if (ifhgrid.eq.1) then

          if (ifreal.eq.1) then
c           bilinear interp:
            call interph (arr(1,1,k), work, xhin, yhin, nxin, nyin,
     *                    ifreal, iffirst)
          else
c           nearest neighbor:
            call interphi (iarr(1,1,k), iwork, xhin, yhin, nxin, nyin,
     *                     ifreal, iffirst)
          endif

        else

c         for u,v grid variables, just set to zero:
          if (ifreal.eq.1) then
            call zero (arr(nx1,ny1,k), (nx2-nx1+1)*(ny2-ny1+1))
          else
            call izero (iarr(nx1,ny1,k), (nx2-nx1+1)*(ny2-ny1+1))
          endif

        endif
c----------
      enddo
c----------

      return
      end

c-----------------------------------------------------------------------

      subroutine interph  (arr, arro,  xho,yho, nxo,nyo, ifreal,iffirst)

      use comicegrid

      entry      interphi (iarr,iarro, xho,yho, nxo,nyo, ifreal,iffirst)

c        Horizontal interpolation or nearest neigbor for one horizontal
c        slice, on h-grid, from arro(nxo,nyo) with coords xho(nxo,nyo),
c        yho(nxo,nyo), to arr(nx,ny).
c        If iffirst=1, compute and save bilinedar interpolation
c        indices (indhx,indhy) and weights (weih).
c        If ifreal=1, arrays are real, do bilinear interp.
c        Else, arrays are integer, do nearest neighbor.

      dimension  arr(nx,ny),  arro(nxo,nyo),
     *           iarr(nx,ny), iarro(nxo,nyo)

c        If iffirst=1, find horiz-interp indices and weights (saved)

c---------------------------
      if (iffirst.eq.1) then
c---------------------------
        do j=1,ny
          do i=1,nx
            call findindex (xh(i,j), yh(i,j),
     *                      indhx(1,i,j), indhy(1,i,j), weih(1,i,j),
     *                      xho, yho, 1,nxo, 1,nyo)
          enddo
        enddo
c----------
      endif
c----------

c================
      do j=1,ny
        do i=1,nx
c================
          if (ifreal.eq.1) then
c           bilinear interp, real:
            arr(i,j) = weih(1,i,j)*arro(indhx(1,i,j),indhy(1,i,j))
     *               + weih(2,i,j)*arro(indhx(2,i,j),indhy(2,i,j))
     *               + weih(3,i,j)*arro(indhx(3,i,j),indhy(3,i,j))
     *               + weih(4,i,j)*arro(indhx(4,i,j),indhy(4,i,j))

          else
c           nearest neighbor, integer:
            weimax =
     *        max (weih(1,i,j),weih(2,i,j),weih(3,i,j),weih(4,i,j))
            if (weimax.eq.weih(1,i,j)) then
              iarr(i,j) = iarro(indhx(1,i,j),indhy(1,i,j))
            else if (weimax.eq.weih(2,i,j)) then
              iarr(i,j) = iarro(indhx(2,i,j),indhy(2,i,j))
            else if (weimax.eq.weih(3,i,j)) then
              iarr(i,j) = iarro(indhx(3,i,j),indhy(3,i,j))
            else if (weimax.eq.weih(4,i,j)) then
              iarr(i,j) = iarro(indhx(4,i,j),indhy(4,i,j))
            endif

          endif
c============
        enddo
      enddo
c============

      return
      end

c-----------------------------------------------------------------------

      subroutine findindex (x, y, indx, indy, wei,
     *                      xa, ya, nx1, nx2, ny1, ny2)

c     Finds bilinear interpolation indices (indx,indy) and weights (wei)
c     for location (x,y) in coordinate-arrays xa, ya

      dimension indx(4), indy(4), wei(4),
     *          xa(nx1:nx2,ny1:ny2), ya(nx1:nx2,ny1:ny2)

c       For now, just do for 1-D in x (ny=1)

      call zero (wei, 4)

      if (ny1.eq.ny2) then

        j = ny1
        indy(1) = j
        indy(2) = j
        indy(3) = j
        indy(4) = j

        if (x.le.xa(nx1,j)) then
          indx(1) = nx1
          indx(2) = nx1
          indx(3) = nx1
          indx(4) = nx1
          wei(1)  = 1.
        else if (x.ge.xa(nx2,j)) then
          indx(1) = nx2
          indx(2) = nx2
          indx(3) = nx2
          indx(4) = nx2
          wei(1)  = 1.
        else
          do i=nx1+1,nx2
            if (x.lt.xa(i,j)) then
              indx(1) = i-1
              wei(1)  =  (xa(i,j)-x) / (xa(i,j)-xa(i-1,j))
              indx(2) = i
              wei(2)  = 1. - wei(1)
              indx(3) = indx(1)
              indx(4) = indx(1)
              goto 100
            endif
          enddo
          write (ioterm,*) "findindex: shouldn't get here"
          stop
  100     continue
        endif

      else

        write (ioterm,*) 'findindex: not ready for ny1 < ny2'
        stop

      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine writeres

c     Write restart file
      use comicesheet, only:
     *h, hs, hb, t,
     *baseperc, basefrml, heatb,
     *sedim, sedimeq, sedimun, tsed, wsed, heats,
     *topbed, topbedeq, tbed, equiload,
     *hw, tw, maskpres, sedpres,
     *ua, va, ub, vb, dfu, dfv,
     *sealev, rco2inter, nwrit, timeice

      use comicephys
      use comicegrid
      use comicesparse
      use subgrid, only: rtag

      integer id_time,id_nx,id_ny,id_nxp,id_nyp,id_nxp2,id_nyp2
      integer id_nlevp,id_1d
      integer iou, i,j, n, it(10), iu(10), iv(10), ib(10), ic(10), ntrec
      real tmp
      include "pconst.h"
      include "switch.h"

      real, allocatable :: tmpij(:,:),tmpijk(:,:,:)
      integer ln
      character(120) :: fname,fname2
      character cyr*80
      character cmd*200

      integer numrests

      versres=3.
      ntrec = 1
      c1e20=1.e20

      numrests=1
      if (eorun) numrests=2

      do n=1,numrests
        if (n==1) then
          write (cyr, '(i10)') nint(abs(timeice))
          do i=1,lenchr(cyr)
            if (cyr(1:1).eq.' ') cyr = cyr(2:)
          enddo
          if (timeice.lt.0.) cyr = cyr(1:lenchr(cyr)) // 'm'
          fname=trim(isname)//'_restart'//trim(cyr)//'.nc'
	  print*, 'writing ice restart to', fname
        endif
        if (n==2) then
          fname='data/'//trim(isname)//'_restart.nc'
	  print*, 'writing ice restart to', fname
        endif
!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
        call openfile (fname, iou)
!-----------------------------------------------------------------------
!     start definitions
!-----------------------------------------------------------------------
        call redef (iou)
!-----------------------------------------------------------------------
!     define dimensions
!-----------------------------------------------------------------------
        call defdim ('time', iou, 0, id_time)
        call defdim ('nx',   iou, nx, id_nx)
        call defdim ('ny',   iou, ny, id_ny)
        call defdim ('nxp',  iou, nxp, id_nxp)
        call defdim ('nyp',  iou, nyp, id_nyp)
        call defdim ('nxp',  iou, nx+2, id_nxp2)
        call defdim ('nyp',  iou, ny+2, id_nyp2)
        call defdim ('nlevp', iou, nlevp+1, id_nlevp)
        call defdim ('nbed', iou, nbed, id_nbed)
        call defdim ('nsed', iou, nsed, id_nsed)
        call defdim ('oned', iou, 1, id_1d)
!-----------------------------------------------------------------------
!     define 1d data
!-----------------------------------------------------------------------
        it(1) = id_1d
        call defvar ('timeice', iou, 1, it, c0, c0, 'T', 'D'
     &  , 'timeice', '', '')
        call defvar ('versres', iou, 1, it, c0, c0, 'T', 'D'
     &  , 'versres', '', '')
        call defvar ('nx', iou, 1, it, c0, c0, 'T', 'D'
     &  , 'nx', '', '')
        call defvar ('ny', iou, 1, it, c0, c0, 'T', 'D'
     &  , 'ny', '', '')
        call defvar ('nlev', iou, 1, it, c0, c0, 'T', 'D'
     &  , 'nlev', '', '')
        call defvar ('nsed', iou, 1, it, c0, c0, 'T', 'D'
     &  , 'nsed', '', '')
        call defvar ('nbed', iou, 1, it, c0, c0, 'T', 'D'
     &  , 'nbed', '', '')
        call defvar ('nwrit', iou, 1, it, c0, c0, 'T', 'D'
     &  , 'nwrit', '', '')
        call defvar ('sealev', iou, 1, it, c0, c0, 'T', 'D'
     &  , 'sealev', '', '')
        call defvar ('rtag', iou, 1, it, c0, c0, 'T', 'D'
     &  , 'sealev', '', '')

!-----------------------------------------------------------------------
!     define 2d data (nx,ny)
!-----------------------------------------------------------------------
        it(1) = id_nx
        iu(1) = id_nxp
        iv(1) = id_nxp2
        it(2) = id_ny
        iu(2) = id_nyp
        iv(2) = id_nxp2
        it(3) = id_time
        iu(3) = id_time
        iv(3) = id_time

        call defvar ('xh', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'xh', ' ', ' ')
        call defvar ('yh', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'yh', ' ', ' ')
        call defvar ('h', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'h', ' ', ' ')
        call defvar ('hs', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'hs', ' ', ' ')
        call defvar ('hb', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'hb', ' ', ' ')
        call defvar ('baseperc', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'baseperc', ' ', ' ')
        call defvar ('basefrml', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'basefrml', ' ', ' ')
        call defvar ('heatb', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'heatb', ' ', ' ')
        call defvar ('sedim', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'sedim', ' ', ' ')
        call defvar ('sedimeq', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'sedimeq', ' ', ' ')
        call defvar ('sedimun', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'sedimun', ' ', ' ')
        call defvar ('topbed', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'topbed', ' ', ' ')
        call defvar ('topbedeq', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'topbedeq', ' ', ' ')
        call defvar ('equiload', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'equiload', ' ', ' ')
        call defvar ('hw', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'hw', ' ', ' ')
        call defvar ('tw', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'tw', ' ', ' ')
        call defvar ('maskpres', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'maskpres', ' ', ' ')
        call defvar ('sedpres', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'sedpres', ' ', ' ')
        call defvar ('is2sg', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'is2sg', ' ', ' ')
        call defvar ('iceiind', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'iceiind', ' ', ' ')
        call defvar ('icejind', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &  , 'icejind', ' ', ' ')
        call defvar ('ua', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &  , 'ua', ' ', ' ')
        call defvar ('va', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &  , 'va', ' ', ' ')
        call defvar ('ub', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &  , 'ub', ' ', ' ')
        call defvar ('vb', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &  , 'vb', ' ', ' ')
        call defvar ('dfu', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &  , 'dfu', ' ', ' ')
        call defvar ('dfv', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &  , 'dfv', ' ', ' ')

!-----------------------------------------------------------------------
!     define3d data (nx,ny,nlev)
!-----------------------------------------------------------------------
        it(1) = id_nx
        it(2) = id_ny
        it(3) = id_nlevp
        it(4) = id_time
        call defvar ('t', iou, 4, it,-c1e20, c1e20, ' ', 'D',
     &  't', ' ', ' ')

        it(1) = id_nx
        it(2) = id_ny
        it(3) = id_nbed
        it(4) = id_time
        call defvar ('tbed', iou, 4, it,-c1e20, c1e20, ' ', 'D',
     &  'tbed', ' ', ' ')

        it(1) = id_nx
        it(2) = id_ny
        it(3) = id_nsed
        it(4) = id_time
        call defvar ('tsed', iou, 4, it,-c1e20, c1e20, ' ', 'D',
     &  'tsed', ' ', ' ')
        call defvar ('wsed', iou, 4, it,-c1e20, c1e20, ' ', 'D',
     &  'wsed', ' ', ' ')
        call defvar ('heats', iou, 4, it,-c1e20, c1e20, ' ', 'D',
     &  'heats', ' ', ' ')

        call enddef (iou)

!-----------------------------------------------------------------------
!     write 1d data
!-----------------------------------------------------------------------
        ntrec=1

        tmp=timeice
        call putvars ('timeice', iou, ntrec, tmp, c1, c0)
        tmp=versres
        call putvars ('versres', iou, ntrec, tmp, c1, c0)
        tmp=nx
        call putvars ('nx', iou, ntrec, tmp, c1, c0)
        tmp=ny
        call putvars ('ny', iou, ntrec, tmp, c1, c0)
        tmp=nlev
        call putvars ('nlev', iou, ntrec, tmp, c1, c0)
        tmp=nsed
        call putvars ('nsed', iou, ntrec, tmp, c1, c0)
        tmp=nbed
        call putvars ('nbed', iou, ntrec, tmp, c1, c0)
        tmp=nwrit
        call putvars ('nwrit', iou, ntrec, tmp, c1, c0)
        tmp=sealev
        call putvars ('sealev', iou, ntrec, tmp, c1, c0)
        tmp=rtag
        call putvars ('rtag', iou, ntrec, tmp, c1, c0)

!-----------------------------------------------------------------------
!     write 2d data
!-----------------------------------------------------------------------

        ln=nx*ny
        ib(:) = 1
        ic(:) = 1
        ic(1) = nx
        ic(2) = ny
        ic(3) = ntrec
        allocate(tmpij(nx,ny))
        tmpij = h
        call putvara ('h', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = hs
        call putvara ('hs', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = hb
        call putvara ('hb', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = baseperc
        call putvara ('baseperc', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = basefrml
        call putvara ('basefrml', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = heatb
        call putvara ('heatb', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = sedim
        call putvara ('sedim', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = sedimeq
        call putvara ('sedimeq', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = sedimun
        call putvara ('sedimun', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = topbed
        call putvara ('topbed', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = topbedeq
        call putvara ('topbedeq', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = equiload
        call putvara ('equiload', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = hw
        call putvara ('hw', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = tw
        call putvara ('tw', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = maskpres
        call putvara ('maskpres', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = sedpres
        call putvara ('sedpres', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = real(is2sg)
        call putvara ('is2sg', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = real(iceiind)
        call putvara ('iceiind', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = real(icejind)
        call putvara ('icejind', iou, ln, ib, ic, tmpij, c1, c0)
        deallocate(tmpij)

        ln=nxp*nyp
        ib(:) = 1
        ic(:) = 1
        ic(1) = nxp
        ic(2) = nyp
        ic(3) = ntrec
        allocate(tmpij(nxp,nyp))
        tmpij = ua
      !set nans to 0 if isnan function available (not on kikou compiler)
        do i=1,nxp
          do j=1,nyp
            if (isnan(tmpij(i,j))) tmpij(i,j)= 0.
          enddo
        enddo
        call putvara ('ua', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = va
        do i=1,nxp
          do j=1,nyp
            if (isnan(tmpij(i,j))) tmpij(i,j)= 0.
          enddo
        enddo
        call putvara ('va', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = ub
        do i=1,nxp
          do j=1,nyp
            if (isnan(tmpij(i,j))) tmpij(i,j)= 0.
          enddo
        enddo
        call putvara ('ub', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = vb
        do i=1,nxp
          do j=1,nyp
            if (isnan(tmpij(i,j))) tmpij(i,j)= 0.
          enddo
        enddo
        call putvara ('vb', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = dfu
        do i=1,nxp
          do j=1,nyp
            if (isnan(tmpij(i,j))) tmpij(i,j)= 0.
          enddo
        enddo
        call putvara ('dfu', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij = dfv
        do i=1,nxp
          do j=1,nyp
            if (isnan(tmpij(i,j))) tmpij(i,j)= 0.
          enddo
        enddo
        call putvara ('dfv', iou, ln, ib, ic, tmpij, c1, c0)
        deallocate(tmpij)

!-----------------------------------------------------------------------
!     write 3d data
!-----------------------------------------------------------------------

        ln=nx*ny*(nlevp+1)
        ib(:) = 1
        ic(:) = 1
        ic(1) = nx
        ic(2) = ny
        ic(3) = nlevp+1
        ic(4) = ntrec

        allocate(tmpijk(nx,ny,nlevp+1))
      !t equiv to tracer(...,1)
        tmpijk = t
        call putvara ('t', iou, ln, ib, ic, tmpijk, c1, c0)
        deallocate(tmpijk)

        ln=nx*ny*nbed
        ib(:) = 1
        ic(:) = 1
        ic(1) = nx
        ic(2) = ny
        ic(3) = nbed
        ic(4) = ntrec

        allocate(tmpijk(nx,ny,nbed))
        tmpijk = tbed
        call putvara ('tbed', iou, ln, ib, ic, tmpijk, c1, c0)
        deallocate(tmpijk)

        ln=nx*ny*nsed
        ib(:) = 1
        ic(:) = 1
        ic(1) = nx
        ic(2) = ny
        ic(3) = nsed
        ic(4) = ntrec

        allocate(tmpijk(nx,ny,nsed))
        tmpijk = tsed
        call putvara ('tsed', iou, ln, ib, ic, tmpijk, c1, c0)
        tmpijk = wsed
        call putvara ('wsed', iou, ln, ib, ic, tmpijk, c1, c0)
        tmpijk = heats
        call putvara ('heats', iou, ln, ib, ic, tmpijk, c1, c0)
        deallocate(tmpijk)

        !call closefile (iou)

        !Not used in coupled climate runs
        !write (iures,'(e25.15)') rco2inter

c     in common (used in sedtracking_eul):
        !call writeout  (sedtrack, 1,nx,1,ny,   ntr+1,           1)
        !write (iures,'(e25.15)') timebot
        !call writeouti (itrtop,   1,nx,1,ny,   1,               0)
      enddo
! # if defined 1
!       if (eorun) then
! !-----------------------------------------------------------------------
! !          copy ice sheet restart to data directory
! !-----------------------------------------------------------------------
!       fname2= trim(isname)//'_restart.nc'
!       print*, 'Overwriting ice sheet restart file in /data'
!       cmd = 'cp -p ' // trim(fname) // ' data/' // trim(fname2)
!       print*, 'cmd=',cmd
!       ier = ishell (cmd)
!       endif
! # endif
      return
      end

c-----------------------------------------------------------------------

      subroutine writeout  (arr, nx1, nx2, ny1, ny2, nvert, ifreal)

      use comicegrid

      entry      writeouti (iarr,nx1, nx2, ny1, ny2, nvert, ifreal)

c     Writes one field (all horiz slices) to restart file

      dimension  arr(nx1:nx2, ny1:ny2, nvert),
     *          iarr(nx1:nx2, ny1:ny2, nvert)

      do k=1,nvert
        if (ifreal.eq.1) then
          write (iures,'(e25.15)') (( arr(i,j,k),i=nx1,nx2), j=ny1,ny2)
        else
          write (iures,'(i10)') ((iarr(i,j,k),i=nx1,nx2), j=ny1,ny2)
        endif
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine writenest (h, hb, t, ua, va)

c     Writes fields needed for subsequent nested run, for initializing
c     and specifying edge b.c's. Called at end of a continental-scale
c     run. File will be read by readnest at start of a later nested run.

      use comicephys
      use comicegrid

      dimension
     *  h(nx,ny),        hb(nx,ny),      t(nx,ny,0:nlevp),
     *  ua(0:nxp,0:nyp), va(0:nxp,0:nyp)

      write (iunest,'(2i6)') nx,ny,nlev
      write (iunest,'(e25.15)') h
      write (iunest,'(e25.15)') hb
      do k=0,nlevp
        write (iunest,'(e25.15)') ((t(i,j,k),i=1,nx),j=1,ny)
      enddo
      write (iunest,'(e25.15)') ua
      write (iunest,'(e25.15)') va

      return
      end

c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

