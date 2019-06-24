! source file: /net/mare/home1/eby/as/ism/iceocean.F
      subroutine ocean (zclim,h, hb, hw, maskwater, maskpres,
     *                  oceanmelt, calvrate, sealev,
     *                  dtantann, dtantjan, rco2,
     *                  timeice, weirun, dt, firstocn, iloop)

c        Set basal melting (oceanmelt) for floating ice, and calving
c        (calvrate) at floating ice edge.

      use comicephys
      use comicegrid

      dimension
     *  h(nx,ny),          hb(nx,ny),        hw(nx,ny),
     *  maskwater(nx,ny),  maskpres(nx,ny),
     *  oceanmelt(nx,ny),  calvrate(nx,ny)

      logical firstocn

c     for calving at ice shelf edge:
      dimension ish(4), jsh(4)
      save ish, jsh
c     E-W/N-S:
      data ish /-1, 1,  0, 0 /
      data jsh / 0, 0, -1, 1 /
      parameter (nlook = 4)
c     parameter (calvr =10.0)                   ! calving rate, m/yr
c     parameter (calvr = 2.0)                   ! calving rate, m/yr
      parameter (calvr = 0.)

c     dimension work(nx,ny)
      parameter (dtarc = 0.)   ! yr-interval to call oceanarc (0=always)
c     parameter (dtarc = 20.)

c     ocean melt for 3 zones (ocmelt[a-c]) and 3 climates (1:3):
c     3 zones: a = protected by arc, s = unprotected shelf, b = deep ocn
c     3 climates: 1 = lgm, 2 = mod, 3 = hot.
      dimension ocmelta(3), ocmelts(3), ocmeltb(3)
      save ocmelta, ocmelts, ocmeltb
c     data ocmelta / 0., .3, 1./     ! mis31 a,b
c     data ocmelts / 0., 2., 6./     ! mis31 a,b
c     data ocmeltb / 2., 5., 6./     ! mis31 a,b
ccur  data ocmelta / 0., .3, 2./     ! mis31 c
      data ocmelta / 0., .1, 2./     ! mis31 c

c777  data ocmelts / 0., 2.,10./     ! mis31 c
ccur  data ocmelts / 0., 3.,10./     ! mis31 c
      data ocmelts / 0., 5.,10./     ! mis31 c
c777  data ocmelts / 0., 4.,10./     ! mis31 c

      data ocmeltb / 2., 5.,10./     ! mis31 c

      if (isname == 'Greenland') then
      !Set oceanmelt very high in Greenland to avoid any ice shelf form'n
        do j=1,ny
          do i=1,nx
            if (maskwater(i,j).eq.1) then
              oceanmelt(i,j) = 100.
            else
              oceanmelt(i,j) = 0.
            endif
          enddo
        enddo
        return
      endif

c--------------------
c------------------------------------------

c        Calculate degrees of arc with straight-line access to deep ocn

      n = max (1, nint (dtarc/dt))
      if (mod (iloop, n).eq.0 .or. firstocn) then
        call oceanarc (h, maskwater)
c       call printmap (timeice,h,  'ice thickness (m)', 150., h, 19,0)
c       call printmap (timeice,h,  'ice thickness_10 (m)', 10., h, 19,0)
c       do j=1,ny
c         do i=1,nx
c           work(i,j) = 10.*maskwater(i,j)
c         enddo
c       enddo
c       call printmap (timeice,work,'maskwater*10',1.,h,19,0)
c       call printmap (timeice,arcocn,  'arcocn (deg)', 10., h, 19, 0)
c       call printmap (timeice,distocn, 'distocn (m)', dx0, h, 19, 0)
c       call printmap (timeice,distgl,  'distgl  (m)', dx0, h, 19, 0)
c       stop
      endif

c        Set current oceanmelt. First set weights zw[lgm,mod,hot]
c        depending on "climate" (zclim = f(dtant,...)) to set ocean
c        melts for 3 zones (zoca,s,b). Then for each ocean point,
c        calculate weights for each zone (zws,zwb) depending on
c        arc to open ocean (arcocn), bathymetry, distocn, distgl.

      if (zclim.lt.1.) then
c8888   zwmod = max (0., (zclim-.3)/.7)
        zwmod = max (0., (zclim-.4)/.6)  ! 8888
c       zwmod = max (0., min (1., (zclim-.4)/.5))  ! 8888a
        zwlgm = 1.-zwmod
        zwhot = 0.
      else
        zwhot = zclim - 1.
        zwmod = 1.-zwhot
        zwlgm = 0.
      endif
      weirun = zclim  ! diagnostic only

c     set ocean melt for 3 zones.
c     zoca = protected by arc, zocs = unprotected shelf, zocb =deep ocn.
c     zw[lgm,mod,hot] = weights for current climate (lgm, modern, hot)
      zoca = zwlgm*ocmelta(1) + zwmod*ocmelta(2) + zwhot*ocmelta(3)
      zocs = zwlgm*ocmelts(1) + zwmod*ocmelts(2) + zwhot*ocmelts(3)
      zocb = zwlgm*ocmeltb(1) + zwmod*ocmeltb(2) + zwhot*ocmeltb(3)

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then
c           set zone weights (zws for zones a vs s, zwb for zone b):
            zws = max (0., min (1., (arcocn(i,j)-80.)/30.)) ! 80, 90...
            zwb = max (0., min (1., (sealev-hb(i,j)-1400.)/200.))

c           more melting N. of Ross Island, W. Ross Sea:.
c           if (alond(i,j).gt.160. .and. alatd(i,j).gt.-78.) then
c             zf = min(1., ((180.-alond(i,j))/5.)*((alatd(i,j)+78.)/2.))
c             zws = (1.-zf)*zws + zf*1.
c           endif

            zws = zws * exp (-distocn(i,j)/100.e3)   ! cone

c           set ocean melt depending on zone:
            oceanmelt(i,j) =   (1.-zwb) * ((1.-zws)*zoca + zws*zocs)
     *                       +     zwb  * zocb

c           or: still protected if deep and protected (problem with
c           excessive ice extent at LGM!)
!             if (i==1.and.j==1)print*, 'NOTE: using Berkner Island patch'
!             oceanmelt(i,j) =       zws  * ((1.-zwb)*zocs + zwb*zocb)
!      *                       + (1.-zws) * zoca

            oceanmelt(i,j) = max (oceanmelt(i,j), 0.)

ctwo:
c           oceanmelt(i,j) = oceanmelt(i,j) *
c    *        max (exp(-distocn(i,j)/100.e3), exp(-distgl(i,j)/100.e3))

          else
            oceanmelt(i,j) = 0.
          endif
        enddo
      enddo

c        Set oceanmelt for grounded ice points adjacent to and
c        physically in contact with ocean

      do j=1,ny
        do i=1,nx
          if ( maskwater(i,j).ne.1 .and.
c    *         (hb(i,j).lt.sealev.and.h(i,j).gt.0.)
     *         (hb(i,j).lt.sealev)
     *       ) then
            nprox = 0
            oprox = 0.
            do look=1,nlook
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (maskwater(ii,jj).eq.1) then
                nprox = nprox + 1
                oprox = oprox + oceanmelt(ii,jj)
              endif
            enddo
            oceanmelt(i,j) = oprox / max (nprox,1)
          endif
        enddo
      enddo

c------------------------------------
c-----

c        Calving at ice edge. Only if adjacent (E-W,N-S) to open ocean,
c        and depending on thickness of thickest adjacent ice pt (not
c        local thickness, to allow advance of ice front if thin for
c        several grid pts near edge). Also compute for current
c        open-ocean pts adjacent to ice ("hypothetical")

      call zero (calvrate, nx*ny)

c     if (calvr.ne.0.) then
c       do j=1,ny
c         do i=1,nx
cc          if (maskwater(i,j).eq.1 .and. h(i,j).gt.0) then
c           if (maskwater(i,j).eq.1) then         ! include hypothetical
c             heff = h(i,j)
c             ncalv = 0
c             do look=1,nlook
c               ii = max (1, min (nx, i + ish(look)))
c               jj = max (1, min (ny, j + jsh(look)))
c               heff = max (heff, h(ii,jj))
c               if (maskwater(ii,jj).eq.1 .and. h(ii,jj).eq.0.)
c    *            ncalv = ncalv + 1
c             enddo
cc            if (ncalv.gt.0 .and. heff.gt.0.) !hypoth just 1pt adjacent
c             if (ncalv.gt.0)                  !hypoth for all ocean
c    *          calvrate(i,j) = ncalv * ((max(1.-(heff/50.),0.))**0.3)
c    *                          * calvr
c           endif
c         enddo
c       enddo
c     endif

c     if (calvr.ne.0.) then
c       do j=1,ny
c         do i=1,nx
c           if (maskwater(i,j).eq.1 .and. h(i,j).lt.100.) then
c             ncalv = 0
c             nland = 0
c             do look=1,nlook
c               ii = max (1, min (nx, i + ish(look)))
c               jj = max (1, min (ny, j + jsh(look)))
c               if (maskwater(ii,jj).eq.1 .and. h(ii,jj).eq.0.)
c    *            ncalv = ncalv + 1
c               if (maskwater(ii,jj).eq.0) nland = nland + 1
c             enddo
c             if (nland.eq.0 .and. ncalv.gt.0) calvrate(i,j) = calvr
c           endif
c         enddo
c       enddo
c     endif

      return
      end

c-----------------------------------------------------------------------

      subroutine oceanarc (h, maskwater)

c       For each ocean point, computes degrees of arc (0 to 360) with
c       straight-line all-ocean paths to open ocean (domain edge).

c       For inland lakes, set arcocn=1 (so full diffusion in oceantoc)

      use comicephys
      use comicegrid

      dimension h(nx,ny), maskwater(nx,ny)
      parameter (narc=72)   ! number of directions (every  5 deg)
c     parameter (narc=36)   ! number of directions (every 10 deg)

      dimension angarc(narc), tanarc(narc)
      logical firstarc
      data firstarc /.true./
      save angarc, tanarc, firstarc

      parameter (nlook = 9)
      dimension ish(nlook), jsh(nlook)
      save ish, jsh
      data ish /0, -1, -1, -1,  0,  1, 1, 1, 0/
      data jsh /0,  1,  0, -1, -1, -1, 0, 1, 1/

      if (firstarc) then
        do m=1,narc
          angarc(m) = -pi + (m-0.5)*(2.*pi)/narc
          tanarc(m) = tan(angarc(m))
        enddo
        firstarc = .false.
      endif

      do jo=1,ny
        do io=1,nx
          arcocn(io,jo) = 0.
          distocn(io,jo) = 10000.e3
          distgl(io,jo)  = 10000.e3

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          if (maskwater(io,jo).eq.1) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c           don't do for ocean pts with no adjacent h or land, for speed
            ifdo = 0
            do look=1,nlook
              ii = max (1, min (nx, io + ish(look)))
              jj = max (1, min (ny, jo + jsh(look)))
              if (h(ii,jj).gt.0. .or. maskwater(ii,jj).ne.1) ifdo = 1
c             if (h(ii,jj).gt.0.) ifdo = 1
            enddo
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            if (ifdo.eq.1) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c~~~~~~~~~~~~~~~~~~~~~~~
             do m=1,narc
c~~~~~~~~~~~~~~~~~~~~~~~
               ang = angarc(m)
               tang = tanarc(m)

c                 Either go along one cell in i direction, and use
c                 nearest j to current direction (ifi=1),or v.v.
c                 for j and i(ifi=0):

               if (abs(ang).ge.0.75*pi) then
                 ifi = 1
                 idir = -1
               else if (ang.ge.-0.75*pi .and. ang.le.-0.25*pi) then
                 ifi = 0
                 idir = -1
               else if (abs(ang).le.0.25*pi) then
                 ifi = 1
                 idir = 1
               else
                 ifi = 0
                 idir = 1
               endif

               if (ifi.eq.1) then
                 i = io
  100            i = i + idir
                 if (i.lt.1.or.i.gt.nx) goto 200
                 j = nint ((i-io)*tang) + jo
                 if (j.lt.1.or.j.gt.ny) goto 200
c                if (maskwater(i,j).ne.1) goto 300
                 if (maskwater(i,j).ne.1) then
                   zdist = dx0 * sqrt(float((i-io)**2+(j-jo)**2)) !polst
                   distgl(io,jo) = min (distgl(io,jo), zdist)
                   goto 300
                 endif
                 if (h(i,j).eq.0.) then
                   zdist = dx0 * sqrt(float((i-io)**2+(j-jo)**2)) !polst
                   distocn(io,jo) = min (distocn(io,jo), zdist)
                 endif
                 go to 100
               else
                 j = jo
  150            j = j + idir
                 if (j.lt.1.or.j.gt.ny) goto 200
                 i = nint ((j-jo)/tang) + io
                 if (i.lt.1.or.i.gt.nx) goto 200
c                if (maskwater(i,j).ne.1) goto 300
                 if (maskwater(i,j).ne.1) then
                   zdist = dx0 * sqrt(float((i-io)**2+(j-jo)**2)) !polst
                   distgl(io,jo) = min (distgl(io,jo), zdist)
                   goto 300
                 endif
                 if (h(i,j).eq.0.) then
                   zdist = dx0 * sqrt(float((i-io)**2+(j-jo)**2)) !polst
                   distocn(io,jo) = min (distocn(io,jo), zdist)
                 endif
                 go to 150
               endif
  200          arcocn(io,jo) = arcocn(io,jo) + 360./narc
  300          continue

c~~~~~~~~~~~~~~~~~
             enddo
c~~~~~~~~~~~~~~~~~
c>>>>>>>>>>>>>>>
            else
c>>>>>>>>>>>>>>>
              arcocn(io,jo) = 360.
              distocn(io,jo)= 0.
              distgl(io,jo) = 0.
c>>>>>>>>>>>>>>>>
            endif
c>>>>>>>>>>>>>>>>
c>>>>>>>>>>>>>>
          endif
c>>>>>>>>>>>>>>
        enddo
      enddo

      return
      end
