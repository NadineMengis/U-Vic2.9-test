! source file: /net/mare/home1/eby/as/ism/icetherm.F
c-----------------------------------------------------------------------

      subroutine icetherm (h, hs, t, hw, tw,
     *                     budgsnow, budgrain, budgevap, budgmelt,
     *                     baseperc, basefrml,
     *                     tsurf, tsurfi,
     *                     arhap, s2a0,
     *                     heati, heatb, heatw, geoflux,
     *                     w, wa,
     *                     sedim, tsed, wsed, heats, tbed,
     *                     u, v, ub, vb, dfu, dfv,
     *                     timeice, dt)

c        Does temperature (and tracer) calculations. Steps ice
c        temperatures (t), sediment temperatures (tsed) and
c        bed temperatures through dt. Also sets wsed, basefrml, w
c        and wa (all on h-grid).
c
c        t,tsed,wsed,tbed,tracer,w,basefrml,heati,heatb,heatw,heats
c        are on the main h-grid, as are tsurf[i],hs,h,budg*,base*,sedim.
c        ub,u,uh are on u-grid (staggered 1/2 box to right)
c        vb,v,vh are on v-grid (staggered 1/2 box to top)
c        (ie, an Arakawa "C" grid).
c
c        t,u,v,w(=dzeta/dt),wa(=dz/dt),uh,vh are on mid-levels
c        (0:nlevp). Extra levels k=0 (sfc) and k=nlevp (base) are
c        set for completeness, and are not prognostic.
c
c        Vertical integrals for u, uh, and heati (involving arh) are
c        done following Ritz et al.(1997), using:
c          s1a (at mid levels, for u,v)
c          s2a (at mid levels, for uh,vh and hence w)
c          s3a (integral over each layer, for heati).
c        Top-level s2a0 is set for icedyn to compute sheet-flow
c        dynamics
c
c        Calcs for u, uh use driving force dfu, dfv modified by
c        shelf eqns left-hand sides (=basal stress), passed from icedyn.
c
c        nb: ice temps t are equivalenced to tracer(*,*,*,1).

      use comicephys
      use comicegrid

c     passed:
      dimension
     *  h(nx,ny),            hs(nx,ny),
     *  t(nx,ny,0:nlevp),
     *  hw(nx,ny),           tw(nx,ny),
     *  budgsnow(nx,ny),     budgrain(nx,ny),
     *  budgevap(nx,ny),     budgmelt(nx,ny),
     *  baseperc(nx,ny),     basefrml(nx,ny),
     *  tsurf(nx,ny),        tsurfi(nx,ny),
     *  arhap(nx,ny),        s2a0(nx,ny),
     *  heati(nx,ny,nlev),   heatb(nx,ny),
     *  heatw(nx,ny),        geoflux(nx,ny),
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp),
     *  sedim(nx,ny),        tsed(nx,ny,nsed),    wsed(nx,ny,nsed),
     *  heats(nx,ny,nsed),
     *  tbed(nx,ny,nbed),
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),    dfv(0:nxp,0:nyp)

c     local:
      dimension
     *  arh(nx,ny,nlev),
     *  s1a(nx,ny,0:nlevp),      s2a(nx,ny,0:nlevp),    s3a(nx,ny,nlev),
     *  uh(0:nxp,0:nyp,0:nlevp), vh(0:nxp,0:nyp,0:nlevp),
     *  tro(nx,ny,0:nlevp),      dhdt(nx,ny)

      call zero (u,  (nxp+1)*(nyp+1)*(nlevp+1))
      call zero (v,  (nxp+1)*(nyp+1)*(nlevp+1))
      call zero (uh, (nxp+1)*(nyp+1)*(nlevp+1))
      call zero (vh, (nxp+1)*(nyp+1)*(nlevp+1))
      call zero (heati, nx*ny*nlev)
      call zero (dhdt, nx*ny)
      call zero (w, nx*ny*(nlevp+1))
      call zero (wa, nx*ny*(nlevp+1))

c       (Re)set true tracers to zero where no ice

!       do j=1,ny
!         do i=1,nx
!           if (h(i,j).eq.0.) then
!             do n=ntraca+1,ntrace
!               do k=0,nlevp
!                 t(i,j,k) = 0.
!               enddo
!             enddo
!           endif
!         enddo
!       enddo

c        Prescribe boundary tracer values

c     if (timeice.ge.-100200.) then
!       if (.true.) then
!
!         do n=ntraca+1,ntrace
!
!           if (n.eq.ntraca+1) then
!             do j=1,ny
!               do i=1,nx
!                 if (h(i,j).gt.0.) then
!                   t(i,j,0) = 10.
!                   t(i,j,1) = 10.
! c                 t(i,j,nlev, n) = 20.
! c                 t(i,j,nlevp,n) = 20.
!                 endif
!               enddo
!             enddo
!
!           else if (n.eq.ntraca+2) then
!             do j=1,ny
!               if (j.eq.(ny+1)/2) then
!                 do k=0,nlevp
!                   do i=1,nx
!                     if (h(i,j).gt.0.) t(i,j,k) = 10.
!                   enddo
!                 enddo
!               endif
!             enddo
!
!           endif
!
!         enddo
!
!       endif

c-----------------------------------------------------------------
c        Vertical integrals of ice rheology-temperature dependence
c        (on h-grid)
c-----------------------------------------------------------------

c        Calculate ice flow coefficients arh (for satable)
c        and arhap (for shelf flow in icedyn)

      call arrhenius (h, t, arh, arhap, tsurf)
c
c        Calculate s1a, s2a, s3a (vertical integral terms) for
c        computation of u,uh,v,vh,heati

      call satable (arh, s1a, s2a, s3a)

c        Copy top level s2a0 for icedyn

      call scopy (nx*ny, s2a(1,1,0), 1, s2a0, 1)

c---------------------------------------------
c        Calculate eastward u, uh (on u-grid).
c        u and uh at layer mid levels.
c---------------------------------------------

      do j = 1,ny
        do i = 1,nx-1
          zh  = 0.5*(h(i,j) + h(i+1,j))
          if (zh.gt.0.) then
            if (.not. (hw(i,j).gt.hwcut. and. hw(i+1,j).gt.hwcut) ) then
              zdfv2 = 0.25*(dfv(i,j)+dfv(i+1,j)+dfv(i,j-1)+dfv(i+1,j-1))
              zdf   = sqrt (dfu(i,j)**2 + zdfv2**2)
              ziu   = 2. * (zdf**powi) * zh * dfu(i,j)/max(zdf,1.e-10)
              ziu2  = ziu * zh
              do k=0,nlevp
                zs1a = 0.5*(s1a(i,j,k) + s1a(i+1,j,k))
                zs2a = 0.5*(s2a(i,j,k) + s2a(i+1,j,k))
                u (i,j,k) = ub(i,j) - ziu*zs1a
                uh(i,j,k) = ub(i,j)*zh*(1.-zeta(k)) - ziu2*zs2a
              enddo
            else
              do k=0,nlevp
                u (i,j,k) = ub(i,j)
                uh(i,j,k) = ub(i,j)*zh*(1.-zeta(k))
              enddo
            endif
          endif
        enddo
      enddo

c----------------------------------------------
c        Calculate northward v, vh (on v-grid).
c        v and vh at layer mid levels.
c----------------------------------------------

      do j = 1,ny-1
        do i = 1,nx
          zh  = 0.5*(h(i,j) + h(i,j+1))
          if (zh.gt.0.) then
            if (.not. (hw(i,j).gt.hwcut. and. hw(i,j+1).gt.hwcut) ) then
              zdfu2 = 0.25*(dfu(i-1,j)+dfu(i,j)+dfu(i-1,j+1)+dfu(i,j+1))
              zdf   = sqrt (zdfu2**2 + dfv(i,j)**2)
              ziv   = 2. * (zdf**powi) * zh * dfv(i,j)/max(zdf,1.e-10)
              ziv2  = ziv * zh
              do k=0,nlevp
                zs1a = 0.5*(s1a(i,j,k) + s1a(i,j+1,k))
                zs2a = 0.5*(s2a(i,j,k) + s2a(i,j+1,k))
                v (i,j,k) = vb(i,j) - ziv*zs1a
                vh(i,j,k) = vb(i,j)*zh*(1.-zeta(k)) - ziv2*zs2a
              enddo
            else
              do k=0,nlevp
                v (i,j,k) = vb(i,j)
                vh(i,j,k) = vb(i,j)*zh*(1.-zeta(k))
              enddo
            endif
          endif
        enddo
      enddo

c--------------------------------------------
c        Calculate heati (J/m3/y) over layers
c--------------------------------------------

      do j = 1,ny
        do i = 1,nx
          if (h(i,j).gt.0. .and..not. hw(i,j).gt.hwcut) then
            zdf = sqrt (   (0.5*(dfu(i,j)+dfu(i-1,j)))**2
     *                   + (0.5*(dfv(i,j)+dfv(i,j-1)))**2 )
            do k=1,nlev
              heati(i,j,k) = 2. * (zdf**(powi+1)) * s3a(i,j,k)
            enddo
          endif
        enddo
      enddo

c---------------------------------------------------
c        Calculate w (h*dzeta/dt) on h-grid.
c        (uh,vh and hence w are at layer mid points)
c---------------------------------------------------

      do j=1,ny
        do i=1,nx
          if (h(i,j).gt.0.) then
c           dhdt due to ice_flow+surface_budget+basefrml only
c           (kept in array just for wa below):
            do k=0,nlevp
              zuh = (  uh(i,j,k)*dyu(i,j)-uh(i-1,j,k)*dyu(i-1,j)
     *               + vh(i,j,k)*dxv(i,j)-vh(i,j-1,k)*dxv(i,j-1) )
     *            / darea(i,j)
              if (k.eq.0) dhdt(i,j) = budgsnow(i,j) - budgevap(i,j)
     *                              + budgrain(i,j) - baseperc(i,j)
     *                              + basefrml(i,j) - zuh
              w(i,j,k) = -basefrml(i,j) + dhdt(i,j)*(1.-zeta(k)) + zuh
            enddo
c           surface and base kinematic bcs (assuming rain+melt percolate
c           freely to bottom):
c           w(i,j,0) = budgsnow(i,j) - budgevap(i,j) - budgmelt(i,j)
c           w(i,j,nlevp) = -basefrml(i,j)
          endif
        enddo
      enddo

c        Compute diagnostic vertical velocity wa (on h-grid mid levels)
c        (dz/dt, m/yr). For diagnostic output only.

      do j=1,ny
        jm = max (j-1,1)
        jp = min (j+1,ny)
        do i=1,nx
          if (h(i,j).gt.0.) then
            im = max (i-1,1)
            ip = min (i+1,nx)
            zsx = (hs(ip,j)-hs(im,j)) / (dxu(i,j)+dxu(im,j))
            zhx = (h (ip,j)-h (im,j)) / (dxu(i,j)+dxu(im,j))
            zsy = (hs(i,jp)-hs(i,jm)) / (dyv(i,j)+dyv(i,jm))
            zhy = (h (i,jp)-h (i,jm)) / (dyv(i,j)+dyv(i,jm))
            do k=0,nlevp
              zus = 0.5*(u(i,j,k)+u(i-1,j,k))
              zvs = 0.5*(v(i,j,k)+v(i,j-1,k))
              wa(i,j,k) = -w(i,j,k)
     *                    + dhdt(i,j)*(1.-zeta(k))
     *                    + zus*(zsx-zeta(k)*zhx)
     *                    + zvs*(zsy-zeta(k)*zhy)
            enddo
          endif
        enddo
      enddo

c----------------------------
c        Convert w to zetadot
c----------------------------

      wmax = 0.5/dt
      do j=1,ny
        do i=1,nx
          if (h(i,j).gt.0.) then
            zh = max (h(i,j),1.e-6)
            do k=0,nlevp
              w(i,j,k) = max ( -wmax, min ( wmax, w(i,j,k)/zh ))
            enddo
          endif
        enddo
      enddo

c********************
c****

c---------------------------------------------------------
c        Advect ice tracers.
c        (ice temperature is equivalenced to first tracer)
c---------------------------------------------------------

c        If first call, store (non-uniform) zeta-grid coefficients
c        for upstream parabolic gradients (akpar, bkpar).
c        Nb: w is at layer mid pts, and already on h-grid horizontally,
c        as needed for this parabolic method. (u,v are interpolated
c        horizontally to h-grid mid pts below, and are already on
c        layer mid pts vertically).

      if (firstadv) then
        do k=1,nlev

c         backwards zeta direction (if w +ve):
          if (k.gt.1) then
            dz1 = dzetah(k-1)
            dz2 = dzetah(k-2)
            akpar(0,k) =  (-(dz1+dz2)**2 + dz1**2) / (dz1*dz2*(dz1+dz2))
            akpar(1,k) =  ( (dz1+dz2)**2)          / (dz1*dz2*(dz1+dz2))
            akpar(2,k) =  (- dz1**2)               / (dz1*dz2*(dz1+dz2))
          else
            akpar(0,k) = -1./dzetah(0)
            akpar(1,k) =  1./dzetah(0)
            akpar(2,k) =  0.
          endif

c         forwards zeta direction (if w -ve):
          if (k.lt.nlev) then
            dz1 = dzetah(k)
            dz2 = dzetah(k+1)
            bkpar(0,k) =  (-(dz1+dz2)**2 + dz1**2) / (dz1*dz2*(dz1+dz2))
            bkpar(1,k) =  ( (dz1+dz2)**2)          / (dz1*dz2*(dz1+dz2))
            bkpar(2,k) =  (- dz1**2)               / (dz1*dz2*(dz1+dz2))
          else
            bkpar(0,k) = -1./dzetah(nlev)
            bkpar(1,k) =  1./dzetah(nlev)
            bkpar(2,k) =  0.
          endif

        enddo
        firstadv = .false.
      endif

c===

c------------------
      do n=1,ntrace
c------------------
      call scopy (nx*ny*(nlevp+1), t(1,1,0), 1, tro, 1)

c==================
c====
c     upstream parabolic advection:
      cfllim = 0.2
      do k=1,nlev
        km1 = k-1
        km2 = max (k-2,0)
        kp1 = k+1
        kp2 = min (k+2,nlevp)
        wcfl = cfllim * dzeta(k) / dt
        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.0.) then

              jm1 = max (j-1,1)
              jm2 = max (j-2,1)
              jp1 = min (j+1,ny)
              jp2 = min (j+2,ny)

              im1 = max (i-1,1)
              im2 = max (i-2,1)
              ip1 = min (i+1,nx)
              ip2 = min (i+2,nx)

c             Parabolic interp should only include ice points.
c             Constrain adjacent points to remain within ice.
              if (h(i,jm1).eq.0.) then
                jm1 = j
                jm2 = j
              else if (h(i,jm2).eq.0.) then
                jm2 = jm1
              endif
              if (h(i,jp1).eq.0.) then
                jp1 = j
                jp2 = j
              else if (h(i,jp2).eq.0.) then
                jp2 = jp1
              endif

              if (h(im1,j).eq.0.) then
                im1 = i
                im2 = i
              else if (h(im2,j).eq.0.) then
                im2 = im1
              endif
              if (h(ip1,j).eq.0.) then
                ip1 = i
                ip2 = i
              else if (h(ip2,j).eq.0.) then
                ip2 = ip1
              endif

              ucfl = cfllim * dx(i,j) / dt
              vcfl = cfllim * dy(i,j) / dt

              zu = min (ucfl, max (-ucfl, 0.5*(u(i-1,j,k)+u(i,j,k))))
              zv = min (vcfl, max (-vcfl, 0.5*(v(i,j-1,k)+v(i,j,k))))
              zw = min (wcfl, max (-wcfl, w(i,j,k)))

              t(i,j,k) = tro(i,j,k) + dt *
     *        (  max(zu,0.)/(2.*dx(i,j))
     *           * (-3.*tro(i,j,k) + 4.*tro(im1,j,k) - 1.*tro(im2,j,k))
     *         - min(zu,0.)/(2.*dx(i,j))
     *           * (-3.*tro(i,j,k) + 4.*tro(ip1,j,k) - 1.*tro(ip2,j,k))

     *         + max(zv,0.)/(2.*dy(i,j))
     *           * (-3.*tro(i,j,k) + 4.*tro(i,jm1,k) - 1.*tro(i,jm2,k))
     *         - min(zv,0.)/(2.*dy(i,j))
     *           * (-3.*tro(i,j,k) + 4.*tro(i,jp1,k) - 1.*tro(i,jp2,k))

     *         + max(zw,0.) * (   akpar(0,k)*tro(i,j,k)
     *                          + akpar(1,k)*tro(i,j,km1)
     *                          + akpar(2,k)*tro(i,j,km2) )
     *         - min(zw,0.) * (   bkpar(0,k)*tro(i,j,k)
     *                          + bkpar(1,k)*tro(i,j,kp1)
     *                          + bkpar(2,k)*tro(i,j,kp2) )

     *         )

            endif
          enddo
        enddo
      enddo
c=====
c=====

c----------
      enddo
c----------

c------------------------------------------------------------
c        Step ice, water, sediment and bedrock temperatures
c        due to vertical diffusion, internal and ice-sediment
c        basal frictional heating, small water advection
c------------------------------------------------------------

      call vdif (h, t, hw, tw, sedim, tsed, wsed, tbed,
     *           tsurf, tsurfi,
     *           heati, heatw, heatb, heats, geoflux,
     *           budgrain, budgmelt, baseperc, basefrml, dt)

c---------------------------------------------------------
c        Step ice temperatures due to horizontal diffusion
c---------------------------------------------------------

      !call hdif (t, h, hs, dt)

c*****
c*****

c---------------------------------------------------------
c        Set (diagnostic) surface and bottom tracer values
c---------------------------------------------------------

       do n=ntraca+1,ntrace
         do j=1,ny
           do i=1,nx
             t(i,j,0) = t(i,j,1)
             t(i,j,nlevp) = t(i,j,nlev)
           enddo
         enddo
       enddo

c        Diagnostic check for extreme ice temperatures

      ztmin =  1.e20
      ztmax = -1.e20
      do k=0,nlevp
        do j=1,ny
          do i=1,nx
	    t(i,j,k)=max(t(i,j,k),200.)  !Jer set min temp to 200. K
            if (t(i,j,k).lt.ztmin) then
              itmin = i
              jtmin = j
              ktmin = k
              ztmin = t(i,j,k)
            endif
            if (t(i,j,k).gt.ztmax) then
              itmax = i
              jtmax = j
              ktmax = k
              ztmax = t(i,j,k)
            endif
          enddo
        enddo
      enddo
      ztmin = ztmin - tmelt
      ztmax = ztmax - tmelt
      if (ztmin.lt.-120.) then
         write (6,*) 'Warning (icetherm): extreme min ice temperature:'
         write (6,*) 'i,j,k,lon,lat,tmin: ', itmin, jtmin, ktmin,
     *                alond(itmin,jtmin), alatd(itmin,jtmin), ztmin

      endif
      if (ztmax.gt.5.) then
         write (6,*) 'Warning (icetherm): extreme max ice temperature:'
         write (6,*) 'i,j,k,lon,lat,tmax: ', itmax, jtmax, ktmax,
     *                alond(itmax,jtmax), alatd(itmax,jtmax), ztmax
      endif
      if (ztmin.lt.-120. .or. ztmax.gt.5.) stop

      return
      end

c-----------------------------------------------------------------------

      subroutine satable (arh, s1a, s2a, s3a)

c        Computes various vertical integrals needed for vertical
c        dependence of internal velocities (horiz:sa, vert:s2a) and,
c        internal frictional heating (s3a).
c        Uses simple finite difference approx. for dep on T (T constant
c        within each layer, between interfaces zetah), but within
c        each layer, does integrals w.r.t. zeta exactly.

c        All these integrals are on the h-grid horizontally and the
c        zeta grid vertically. In therm, s*a are interpolated to the
c        u,v-grids for horizontal velocity calcs.
c        Extra k=0 and nlevp indices for s1a are at surface and base.

c        Vertical integrals for u, uh, and heati (involving arh) are
c        done following Ritz et al.(1997), using:
c          s1a (at mid levels, for u,v)
c          s2a (at mid levels, for uh,vh and hence w)
c          s3a (integral over each layer, for heati)

      use comicephys
      use comicegrid

      dimension
     *  arh(nx,ny,nlev),
     *  s1a(nx,ny,0:nlevp),  s2a(nx,ny,0:nlevp),  s3a(nx,ny,nlev)

c        if first call, compute vertical integral zeta deps (saved)

      if (firstsatable) then
        do k=1,nlev
c---
c         for s1a (u,v)
          dzpow1a(k,1) = (zeta(k)**(powi+1) - zetah(k-1)**(powi+1))
     *                   / (powi+1.)
          dzpow1a(k,2) = (zetah(k)**(powi+1) - zeta(k)**(powi+1))
     *                   / (powi+1.)

c         for s2a (uh,vh):
          dzpow2a(k,1) = (  (zeta(k)**(powi+1))*(zeta(k)-zetah(k-1))
     *                    - (zeta(k)**(powi+2)-zetah(k-1)**(powi+2))
     *                        / (powi+2.)
     *                   ) / (powi+1.)
          dzpow2a(k,2) = (- (zeta(k)**(powi+1))*(zetah(k)-zeta(k))
     *                    + (zetah(k)**(powi+2)-zeta(k)**(powi+2))
     *                       / (powi+2.)
     *                   ) / (powi+1.)

c         for s3a (heati):
          dzpow3a(k)   = (zetah(k)**(powi+2) - zetah(k-1)**(powi+2))
     *                   / (powi+2.)
c---
        enddo
        firstsatable = .false.
      endif

c         Do vertical integrations for s1a,s2a,s3a

      do j=1,ny
        do i=1,nx

c         at layer mid points, plus top and bottom (for s1a,s2a):
          s1a(i,j,nlevp) = 0.
          s2a(i,j,nlevp) = 0.
          s1a(i,j,nlev)  = arh(i,j,nlev)*dzpow1a(nlev,2)
          s2a(i,j,nlev)  = arh(i,j,nlev)*dzpow2a(nlev,2)
          s3a(i,j,nlev)  = arh(i,j,nlev)*dzpow3a(nlev)/dzeta(nlev)
          do k=nlev-1,1,-1
            s1a(i,j,k)= s1a(i,j,k+1) + arh(i,j,k+1)* dzpow1a(k+1,1)
     *                               + arh(i,j,k)  * dzpow1a(k,  2)

            s2a(i,j,k)= s2a(i,j,k+1) + s1a(i,j,k+1)*(zeta(k+1)-zetah(k))
     *                               + arh(i,j,k+1)* dzpow2a(k+1,1)
     *                               + s1a(i,j,k  )*(zetah(k)-zeta(k))
     *                               - arh(i,j,k  )* dzpow2a(k,  2)

            s3a(i,j,k)   =             arh(i,j,k)*dzpow3a(k)/dzeta(k)
          enddo
          s1a(i,j,0) = s1a(i,j,1) + arh(i,j,1)*dzpow1a(1,1)
          s2a(i,j,0) = s2a(i,j,1) + s1a(i,j,1)*(zeta(1)-zetah(0))
     *                            + arh(i,j,1)* dzpow2a(1,1)

        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine vdif (h, t, hw, tw, sedim, tsed, wsed, tbed,
     *                 tsurf, tsurfi,
     *                 heati, heatw, heatb, heats, geoflux,
     *                 budgrain, budgmelt, baseperc, basefrml, dt)

c     Sets up and calls time-implicit tridiagonal solver (tridia_i)
c     for icesheet+water+sedim+bedrock vertical heat diffusion,
c     as one vertical system from ice surface to bottom of bedrock.
c     BCs are tsurf[i] at sfc, and geothermal heat flux at bottom.
c     Handles large/small/no ice, water or sediment. If "small"
c     (i.e., thickness <= h[i,w,s]min), treated as thin film
c     with no specific heat capacity, and temps set to interface
c     temperature...but with latent heat, and basefrml of the
c     film thickness tracked.

c     Scans from top to bottom, incrementing tridiagonal layer # (klay)
c     and setting tridiag terms (d[1-3], rhs), and keeping track of
c     correspondence to prognostic physical layers (laywhich, laylev).
c     At each increment, sets terms for any local processes and heat
c     exchange with the previous layer above (not below). For interfaces
c     where there is a pressure-melting-point phase constraint
c     (liquid water or ice present), impose this constraint time
c     explicitly, and if it would be violated, hold interface temp
c     at the melt pt, set tridiag terms for that, and use the
c     (time-explicit, pre-time-step) convergence/divergence of heat
c     flux at the interface to increment basal freeze/melt (basefrml).

c     Use fractional area for small water (fracw=min(hw/hwmin,1)),
c     and apply the interface phase constraint separately for water
c     and non-water fractions (fracw and 1-fracw). This means the
c     melt-pt constraint of liquid water is not discontinuous
c     for hw=0 vs > 0.

c     For small and large water, just use one homogeneous water
c     temperature tw. For large water, exchange coefficent with
c     top and bottom of water is condliq, indep of hw.

c     For small water, apply heatw, the heat tendency due to
c     advection from subr movewater. (For large water, tw itself
c     is stepped in movewater). This captures supercooling
c     (as sub-ice small water flows to higher elevation/thinner ice
c     where pressure-melt point is higher, and so freezes).
c     heatw is applied only to fracw area, whereas basal shear
c     heating (heatb) is applied to whole grid box (sloppy?).
c     Nb: heatw,heatb,heats are J/m2/y, heati converted to this below.

c     After tridiag system is solved, distribute temps back to physical
c     layers. Then scan large ice and water temps for > or < tmelt
c     respectively, percolating/refreezing any resulting meltwater
c     down through ice, and set baseperc.

      use comicephys
      use comicegrid

      dimension
     * h(nx,ny),           t(nx,ny,0:nlevp),
     * hw(nx,ny),          tw(nx,ny),
     * sedim(nx,ny),       tsed(nx,ny,nsed),   wsed(nx,ny,nsed),
     * tbed(nx,ny,nbed),
     * tsurf(nx,ny),       tsurfi(nx,ny),
     * heati(nx,ny,nlev),  heatw(nx,ny),       heatb(nx,ny),
     * heats(nx,ny,nsed),  geoflux(nx,ny),
     * budgrain(nx,ny),    budgmelt(nx,ny),
     * baseperc(nx,ny),    basefrml(nx,ny)

      parameter (himin=1.0,    himinc=0.1,
     *           hsmin=0.2,   hsminc=0.1,
     *                        hbminc=0.1 )

      dimension
     *  d1(nlev+1+nsed+nbed)
     *, d2(nlev+1+nsed+nbed)
     *, d3(nlev+1+nsed+nbed)
     *, rhs(nlev+1+nsed+nbed)
     *, tnew(nlev+1+nsed+nbed)
     *, laywhich(nlev+1+nsed+nbed)
     *, laylev(nlev+1+nsed+nbed)

      hwmin=hwcut
      nvert=nlev+1+nsed+nbed
      call zero (baseperc, nx*ny)
      call zero (basefrml, nx*ny)

c        Overall loops over nx and ny

c==============
      do j=1,ny
      do i=1,nx
c==============

      call zero (d1, nvert)
      call zero (d2, nvert)
      call zero (d3, nvert)
      call zero (rhs, nvert)
      klay = 0

c        Add large ice layers to tridiag arrays

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (h(i,j).gt.himin) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        do k=1,nlev

          if (k.gt.1) then
            zdtpre   = zdt
            zcondpre = zcond
          endif

          zdz  = h(i,j)*dzeta(k)
          zch  = cheaticea + cheaticeb*(t(i,j,k)-tmelt)       ! J/kg/K
          zdt  = dt/(rhoice*zch*zdz)                          ! y m2 K/J
          zcond = condicea*exp(-condiceb*t(i,j,k))
     *          / max(0.5*zdz, himinc)                         ! J/y/m/K
          zheati= heati(i,j,k)*zdz               ! from J/m3/y to J/m2/y

          klay = klay + 1
          laywhich(klay) = 1     ! 1=lg ice, 2=lg water, 3=lg sed, 4=bed
          laylev(klay) = k

          if (k.eq.1) then
            d2(klay) = d2(klay) + zdt*zcond + 1.
            rhs(klay) = rhs(klay) + t(i,j,k) + zdt*zheati
     *                            + zdt*zcond*tsurfi(i,j)
          else
            zc = zcond*zcondpre / (zcond + zcondpre)
            d1(klay)   = d1(klay)   - zdt   *zc
            d2(klay)   = d2(klay)   + zdt   *zc + 1.
            d2(klay-1) = d2(klay-1) + zdtpre*zc
            d3(klay-1) = d3(klay-1) - zdtpre*zc
            rhs(klay) = rhs(klay) + t(i,j,k) + zdt*zheati
          endif

        enddo
c~~~~~~~~~~
      endif
c~~~~~~~~~~

c        Set quantities above major interface
c        (large ice, small ice or air above,
c         large water, large sed or bedrock below,
c         small ice, small water, small sed or none at interface).

      ztme = tmelt - dtmdh*h(i,j)    ! pressure melting pt at ice bottom

      if (h(i,j).gt.himin) then
        zdz   = h(i,j)*dzeta(nlev)
        zch   = cheaticea + cheaticeb*(t(i,j,nlev)-tmelt)
        zdt0  = dt/(rhoice*zch*zdz)

		zc0 = condicea*exp(-condiceb*t(i,j,nlev))
     *         / max(0.5*zdz, himinc)
        zt0    = t(i,j,nlev)
        ztmina = ztme     ! liq fraction (a)
        ztmaxa = ztme     ! liq fraction (a)
        ztminb = -1.e20   ! non-liq fraction (b)
        ztmaxb = ztme     ! non-liq fraction (b)

      else if (h(i,j).gt.0.) then
		zc0 = condicea*exp(-condiceb*t(i,j,nlev/2)) / himinc
        zt0    = tsurfi(i,j)
        ztmina = ztme
        ztmaxa = ztme
        ztminb = -1.e20
        ztmaxb = ztme

      else
		zc0 = condair
        zt0    = tsurf(i,j)
        ztmina = ztme
        ztmaxa =  1.e20
        ztminb = -1.e20
        ztmaxb =  1.e20
      endif

c        Set quantities below major interface

      if (hw(i,j).gt.hwmin) then
        zdz  = hw(i,j)
        zch  = cheatliq
        zdt1 = dt/(rholiq*zch*zdz)
        zc1  = condliq
        zt1  = tw(i,j)
        fracw= 1.
        zheats = 0.
        iwhich = 2
      else if (sedim(i,j).gt.hsmin) then
        zdz  = sedim(i,j)*dzsed(1)
        zch  = cheatsed
        zdt1 = dt/(rhosed*zch*zdz)
        zc1  = condsed / max (0.5*zdz, hsminc)
        zt1  = tsed(i,j,1)
        fracw= min (hw(i,j)/max(hwmin,1.e-6), 1.)
        zheats = heats(i,j,1)
        iwhich = 3
      else
        zch  = cheatbed
        zdz  = dzbed(1)
        zdt1 = dt/(rhobed*zch*zdz)
        zc1  = condbed / max (0.5*zdz, hbminc)
        zt1  = tbed(i,j,1)
        fracw= min (hw(i,j)/max(hwmin,1.e-6), 1.)
        zheats = 0.
        iwhich = 4
      endif

c        Compute interface temperature for each areal fraction
c        (liquid and non-liquid), apply temperature limits
c        (liquid/ice phase), and compute explicit heating for layer
c        above and below (zheat0, zheat1), effective conductivity
c        between these layers (zcond), and freeze/melt at interface
c        (zfrml)

c        Liquid fraction (fracw)

      zti = (zc0*zt0 + zc1*zt1 + heatb(i,j) + heatw(i,j)) / (zc0+zc1)
      if (zti.lt.ztmina .or. zti.gt.ztmaxa) then
        zti       = max (ztmina, min (ztmaxa, zti))
        zfrml_a   = zc0*(zti-zt0) + zc1*(zti-zt1) -heatb(i,j)-heatw(i,j)
        zheat0_a  = 0.
        zheat1_a  = 0.
        zcond0_a  = 0.
        zcond1_a  = 0.
        zcond0e_a = zc0
        zcond1e_a = zc1

c       limit water->ice freezing to what would freeze himin (m)
c       in 50 yrs (not one time step to avoid dt dependence).
c       This crudely prevents overshoot due to large heat loss to
c       atmosphere from initial open water or very thin ice, which
c       would decrease quickly (few months) as ice cover forms/thickens
c       to a few 10's cm and insulates the water. Alternative would be
c       to solve a quadratic for new h (implicit h time step) and
c       deduce average zfrml.
        zfrml_a = min (zfrml_a, 1.1*rhoice*hfus*himin/50.)
      else
        zheat0_a  = zc0*(heatb(i,j)+heatw(i,j)) / (zc0+zc1)
        zheat1_a  = zc1*(heatb(i,j)+heatw(i,j)) / (zc0+zc1)
        zfrml_a   = 0.
        zcond0_a  = zc0*zc1 / (zc0 + zc1)
        zcond1_a  = zc0*zc1 / (zc0 + zc1)
        zcond0e_a = 0.
        zcond1e_a = 0.
      endif
      ztinter_a = zti

c        Non-liquid fraction (1-fracw)

      zti = (zc0*zt0 + zc1*zt1 + heatb(i,j)) / (zc0+zc1)
      if (zti.lt.ztminb .or. zti.gt.ztmaxb) then
        zti = max (ztminb, min (ztmaxb, zti))
        zfrml_b   = zc0*(zti-zt0) + zc1*(zti-zt1) - heatb(i,j)
        zheat0_b  = 0.
        zheat1_b  = 0.
        zcond0_b  = 0.
        zcond1_b  = 0.
        zcond0e_b = zc0
        zcond1e_b = zc1
      else
        zfrml_b   = 0.
        zheat0_b  = zc0*heatb(i,j) / (zc0+zc1)
        zheat1_b  = zc1*heatb(i,j) / (zc0+zc1)
        zcond0_b  = zc0*zc1 / (zc0 + zc1)
        zcond1_b  = zc0*zc1 / (zc0 + zc1)
        zcond0e_b = 0.
        zcond1e_b = 0.
      endif
      ztinter_b = zti

c        Average liquid and non-liquid fractions

      zheat0  = fracw*zheat0_a  + (1.-fracw)*zheat0_b
      zheat1  = fracw*zheat1_a  + (1.-fracw)*zheat1_b + zheats
      zfrml   = fracw*zfrml_a   + (1.-fracw)*zfrml_b
      ztinter = fracw*ztinter_a + (1.-fracw)*ztinter_b
      zcond0  = fracw*zcond0_a  + (1.-fracw)*zcond0_b
      zcond1  = fracw*zcond1_a  + (1.-fracw)*zcond1_b
      zcond0e = fracw*zcond0e_a + (1.-fracw)*zcond0e_b
      zcond1e = fracw*zcond1e_a + (1.-fracw)*zcond1e_b

c        Add first layer below major interface to tridiag arrays,
c        applying conductivity terms to layer above if any

      klay = klay + 1
      laywhich(klay) = iwhich
      laylev(klay) = 1

      rhs(klay) = rhs(klay) + zt1
      if (klay.gt.1) then
        d1(klay)    = d1(klay)    - zdt1* zcond1
        d2(klay)    = d2(klay)    + zdt1*(zcond1+zcond1e) + 1.
        d2(klay-1)  = d2(klay-1)  + zdt0*(zcond0+zcond0e)
        d3(klay-1)  = d3(klay-1)  - zdt0* zcond0
        rhs(klay)   = rhs(klay)   + zdt1*(zheat1+zcond1e*ztinter)
        rhs(klay-1) = rhs(klay-1) + zdt0*(zheat0+zcond0e*ztinter)
      else
        d2(klay)    = d2(klay)    + zdt1*(zcond1+zcond1e) + 1.
        rhs(klay)   = rhs(klay)   + zdt1*(zcond1*zt0 + zcond1e*ztinter
     *                                    + zheat0+zheat1)
      endif

      basefrml(i,j) = basefrml(i,j) + zfrml/(rhoice*hfus)

c        If large water, add next layer down (sed or bed) to tridiag

c+++++++++++++++++++++++++++++++
      if (hw(i,j).gt.hwmin) then
c+++++++++++++++++++++++++++++++

        zdz  = hw(i,j)
        zch  = cheatliq
        zdt0 = dt/(rholiq*zch*zdz)
        zc0  = condliq
        zt0  = tw(i,j)

        if (sedim(i,j).gt.hsmin) then
          zdz  = sedim(i,j)*dzsed(1)
          zch  = cheatsed
          zdt1 = dt/(rhosed*zch*zdz)
          zc1  = condsed / max (0.5*zdz, hsminc)
          zt1  = tsed(i,j,1)
          zheats = heats(i,j,1)
          iwhich = 3
        else
          zdz  = dzbed(1)
          zch  = cheatbed
          zdt1 = dt/(rhobed*zch*zdz)
          zc1  = condbed / max (0.5*zdz, hbminc)
          zt1  = tbed(i,j,1)
          zheats = 0.
          iwhich = 4
        endif

c       pressure melting pt for water (always at top):
        ztme = tmelt - dtmdh*h(i,j)
        zti = (zc0*zt0 + zc1*zt1) / (zc0+zc1)
        if (zti.lt.ztme) then
          zti = ztme
          zfrml  = zc0*(zti-zt0) + zc1*(zti-zt1)
          zheat0 = 0.
          zheat1 = zheats
          zcond0 = 0.
          zcond1 = 0.
          zcond0e= zc0
          zcond1e= zc1
        else
          zfrml  = 0.
          zheat0 = 0.
          zheat1 = zheats
          zcond0 = zc0*zc1 / (zc0 + zc1)
          zcond1 = zc0*zc1 / (zc0 + zc1)
          zcond0e= 0.
          zcond1e= 0.
        endif
        ztinters = zti

        klay = klay + 1
        laywhich(klay) = iwhich
        laylev(klay) = 1

        d1(klay)    = d1(klay)    - zdt1* zcond1
        d2(klay)    = d2(klay)    + zdt1*(zcond1+zcond1e) + 1.
        d2(klay-1)  = d2(klay-1)  + zdt0*(zcond0+zcond0e)
        d3(klay-1)  = d3(klay-1)  - zdt0* zcond0
        rhs(klay)   = rhs(klay)   + zt1 + zdt1*(zheat1+zcond1e*zti)
        rhs(klay-1) = rhs(klay-1)       + zdt0*(zheat0+zcond0e*zti)

        basefrml(i,j) = basefrml(i,j) + zfrml/(rhoice*hfus)

c++++++++++
      endif
c++++++++++

c        If large sediment, add sed layers #s 2-nsed to tridiag

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if (sedim(i,j).gt.hsmin) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        zdz = sedim(i,j)*dzsed(1)
        zch = cheatsed
        zdt = dt/(rhosed*zch*zdz)
        zcond = condsed / max (0.5*zdz, hsminc)

        do k=2,nsed
          zdtpre   = zdt
          zcondpre = zcond

          zdz = sedim(i,j)*dzsed(k)
          zch = cheatsed
          zdt = dt/(rhosed*zch*zdz)
          zcond = condsed / max (0.5*zdz, hsminc)

          zc = zcond*zcondpre / (zcond + zcondpre)

          klay = klay + 1
          laywhich(klay) = 3
          laylev(klay) = k
          d1(klay)   = d1(klay)   - zdt   *zc
          d2(klay)   = d2(klay)   + zdt   *zc + 1.
          d2(klay-1) = d2(klay-1) + zdtpre*zc
          d3(klay-1) = d3(klay-1) - zdtpre*zc
          rhs(klay)  = rhs(klay)  + tsed(i,j,k) + zdt*heats(i,j,k)
        enddo

c          Add bedrock layer #1 to tridiag

        zdz  = sedim(i,j)*dzsed(nsed)
        zch  = cheatsed
        zdt0 = dt/(rhosed*zch*zdz)
        zc0  = condsed / max (0.5*zdz, hsminc)

        zdz  = dzbed(1)
        zch  = cheatbed
        zdt1 = dt/(rhobed*zch*zdz)
        zc1  = condbed / max (0.5*zdz, hbminc)

        zcond  = zc0*zc1 / (zc0 + zc1)

        klay = klay + 1
        laywhich(klay) = 4
        laylev(klay) = 1

        d1(klay)    = d1(klay)    - zdt1*zcond
        d2(klay)    = d2(klay)    + zdt1*zcond + 1.
        if (klay.gt.1) then
          d2(klay-1)  = d2(klay-1)  + zdt0*zcond
          d3(klay-1)  = d3(klay-1)  - zdt0*zcond
        endif
c>>>>>>>>>>
      endif
c>>>>>>>>>>

c         Add bedrock layers # 2-nbed to tridiag

      zdz = dzbed(1)
      zch = cheatbed
      zdt = dt/(rhobed*zch*zdz)
      zcond = condbed / max (0.5*zdz, hbminc)

      do k=2,nbed
        zdtpre   = zdt
        zcondpre = zcond

        zdz = dzbed(k)
        zch = cheatbed
        zdt = dt/(rhobed*zch*zdz)
        zcond = condbed / max (0.5*zdz, hbminc)

        zc = zcond*zcondpre / (zcond + zcondpre)

        klay = klay + 1
        laywhich(klay) = 4
        laylev(klay) = k
        d1(klay)   = d1(klay)   - zdt*zc
        d2(klay)   = d2(klay)   + zdt*zc + 1.
        d2(klay-1) = d2(klay-1) + zdtpre*zc
        d3(klay-1) = d3(klay-1) - zdtpre*zc
        rhs(klay)  = rhs(klay)  + tbed(i,j,k)
      enddo

c        Apply geothermal heat to bottom bedrock layer
c        (klay, zdt are for this layer, from above)

      rhs(klay)  = rhs(klay) + zdt*geoflux(i,j)

c        Solve system of equations

      call tridia_i (d2, d3, d1, tnew, rhs, klay, 0)

c        Distribute new prognostic temperatures back to individual
c        arrays for large ice, water, sed, bed

      do m=1,klay
        if (laywhich(m).eq.1) then
          t(i,j,laylev(m)) = tnew(m)
        else if (laywhich(m).eq.2) then
          tw(i,j) = tnew(m)
        else if (laywhich(m).eq.3) then
          tsed(i,j,laylev(m)) = tnew(m)
        else if (laywhich(m).eq.4) then
          tbed(i,j,laylev(m)) = tnew(m)
        endif
      enddo

c         Constrain (large and small) ice temps to be <= pressure-melt
c         pt, and percolate/refreeze surface melt/rain and any internal
c         melt down from top to bottom (neglecting mass change in
c         a layer if refreezing occurs)

      if (h(i,j).gt.0.) then
c              zheat = rhoice*hfus* (budgrain(i,j) + budgmelt(i,j))  !bug 11/09
	zheat = rhoice*hfus* (budgrain(i,j) + budgmelt(i,j))*dt !fix
        do k=1,nlev
          zdz  = h(i,j)*dzeta(k)
          zch  = cheaticea + cheaticeb*(t(i,j,k)-tmelt)
          ztme = tmelt - dtmdh*h(i,j)*zeta(k)

          t(i,j,k) = t(i,j,k) + zheat/(rhoice*zch*zdz)
          zheat = rhoice*zch*zdz * max(t(i,j,k)-ztme, 0.)
          t(i,j,k) = min (t(i,j,k), ztme)
c                 if (k.eq.nlev) baseperc(i,j) = zheat/(rhoice*hfus)  !bug 11/09
          if (k.eq.nlev) baseperc(i,j) = zheat/(rhoice*hfus*dt) !fix
        enddo
      else
        baseperc(i,j) = budgrain(i,j) + budgmelt(i,j)
      endif

c         Constrain water temp to be >= pressure-melting pt
c         (always at top), and increment basal freezing rate (m/yr)
c         (shouldn't occur for small water, only large)

c     if (hw(i,j).gt.hwmin) then
      if (hw(i,j).gt.0.) then
        ztme = tmelt - dtmdh*h(i,j)
        zheat = rholiq*cheatliq*hw(i,j) * max(ztme-tw(i,j), 0.)
        tw(i,j) = max (tw(i,j), ztme)
c             basefrml(i,j) = basefrml(i,j) + zheat/(rhoice*hfus)    !bug11/09
        basefrml(i,j) = basefrml(i,j) + zheat/(rhoice*hfus*dt) !fix
      endif

c         Limit basal freeze/melt to total ice or liquid in
c         50 yrs (not in one timestep to aviod dt dependence)

      basefrml(i,j) = max (-h(i,j)/50., min (hw(i,j)/50.,basefrml(i,j)))
      if (h(i,j) .le. 0.) basefrml(i,j) = 0.0  !Jer

c         Set basal freeze/melt rate to zero for certain experiments

c        (a) Set diagnostic temps for large ice (at 0 and nlevp)
c        (b) Set interface temps for small ice, water, sed
c        (c) Set all temps to tmelt for niceity for no ice, water, sed
c        nb: interface temps recorded above, before prognostic timestep

      if (h(i,j).gt.himin) then
        t(i,j,0) = tsurfi(i,j)
        t(i,j,nlevp) = ztinter
      else if (h(i,j).le.himin .and. h(i,j).gt.0.) then
        do k=0,nlevp
          t(i,j,k) = tsurfi(i,j)*(1.-zeta(k)) + ztinter*zeta(k)
        enddo
      else if (h(i,j).eq.0.) then
        do k=0,nlevp
          t(i,j,k) = tmelt
        enddo
      endif

      if (hw(i,j).eq.0.) then
        tw(i,j) = tmelt
      else if (hw(i,j).le.hwmin) then
        tw(i,j) = ztinter_a
      endif

c         Adjust large sediment ice/liquid ratio (wsed=liquid content
c         of pores, 0-1). Assume always saturated with liquid and/or
c         ice, with poros=0.4, ice density=rholiq, bulk cheat=1000
c         J/kg/K (indep of wsed)

      ztme = tmelt - dtmdh * (h(i,j)+(rholiq/rhoice)*hw(i,j))

      if (sedim(i,j).gt.hsmin) then
        zcoeft = rhosed*1000.                     ! kg/m3 * J/kg
        zcoefw = 0.4*rholiq*hfus                  ! m3/m3 * kg/m3 * J/kg
        do k=1,nsed
          if (tsed(i,j,k).gt.ztme .and. wsed(i,j,k).lt.1.) then
            zheat = min ( (tsed(i,j,k)-ztme)*zcoeft,
     *                    (1.-wsed(i,j,k))  *zcoefw )
            tsed(i,j,k) = tsed(i,j,k) - zheat/zcoeft
            wsed(i,j,k) = wsed(i,j,k) + zheat/zcoefw
          else if (tsed(i,j,k).lt.ztme .and. wsed(i,j,k).gt.0.)then
            zheat = min ( (ztme-tsed(i,j,k))*zcoeft,
     *                     wsed(i,j,k)      *zcoefw )
            tsed(i,j,k) = tsed(i,j,k) + zheat/zcoeft
            wsed(i,j,k) = wsed(i,j,k) - zheat/zcoefw
          endif
          wsed(i,j,k) = max (0., min (1., wsed(i,j,k)))   ! for roundoff
        enddo

      else if (sedim(i,j).le.hsmin .and. sedim(i,j).gt.0.) then
        if (ztinter.ge.ztme) then
          do k=1,nsed
            tsed(i,j,k) = ztinter
            wsed(i,j,k) = 1.
          enddo
        else
          do k=1,nsed
            tsed(i,j,k) = ztinter
            wsed(i,j,k) = 0.
          enddo
        endif
      else
        do k=1,nsed
          tsed(i,j,k) = tmelt
          wsed(i,j,k) = 0.
        enddo
      endif

c==========
      enddo
      enddo
c==========

      return
      end

c-----------------------------------------------------------------------

      subroutine hdif (t, h, hs, dt)

c        Horizontal diffusion of heat within ice. Step t through dt.
c        Simple time explicit diffusion, allow for non-orthogonal
c        (sloping) zeta coords.

      use comicephys
      use comicegrid

      dimension t(nx,ny,0:nlevp), h(nx,ny), hs(nx,ny)

c     local:
      dimension horflx(nx,ny), horfly(nx,ny), to(nx,ny,0:nlevp)

      call scopy (nx*ny*(nlevp+1), t, 1, to, 1)

c----------------
      do k=1,nlev
c----------------

        do j=1,ny
          do i=1,nx
            ip = min(i+1,nx)
            if (h(i,j).gt.0. and. h(ip,j).gt.0.) then
              zdtdx = (to(ip,j,k)-to(i,j,k)) / dxu(i,j)
              zdtds = 0.5 *
     *             (to(i,j,k+1)-to(i,j,k-1) + to(ip,j,k+1)-to(ip,j,k-1))
     *             / (zeta(k+1)-zeta(k-1))
              zdsdx = ((hs(ip,j)-hs(i,j)) - zeta(k)*(h(ip,j)-h(i,j)))
     *                / dxu(i,j)
              zcond = condicea*exp(-condiceb*0.5*(to(i,j,k)+to(ip,j,k)))
              horflx(i,j) = zcond * (zdtdx + zdtds*zdsdx) * dyu(i,j)
            else
              horflx(i,j) = 0.
            endif

            jp = min(j+1,ny)
            if (h(i,j).gt.0. and. h(i,jp).gt.0.) then
              zdtdy = (to(i,jp,k)-to(i,j,k)) / dyv(i,j)
              zdtds = 0.5 *
     *             (to(i,j,k+1)-to(i,j,k-1) + to(i,jp,k+1)-to(i,jp,k-1))
     *             / (zeta(k+1)-zeta(k-1))
              zdsdy = ((hs(i,jp)-hs(i,j)) - zeta(k)*(h(i,jp)-h(i,j)))
     *                / dyv(i,j)
              zcond = condicea*exp(-condiceb*0.5*(to(i,j,k)+to(i,jp,k)))
              horfly(i,j) = zcond * (zdtdy + zdtds*zdsdy) * dxv(i,j)
            else
              horfly(i,j) = 0.
            endif
          enddo
        enddo

        do j=1,ny
          do i=1,nx
            im = max (i-1,1)
            jm = max (j-1,1)
            zcheat = cheaticea + cheaticeb*(to(i,j,k)-tmelt)
            t(i,j,k) = to(i,j,k)  + dt * (- horflx(i,j) + horflx(im,j)
     *                                    - horfly(i,j) + horfly(i,jm))
     *                                / (rhoice*zcheat*darea(i,j))
          enddo
        enddo

c----------
      enddo
c----------

      return
      end

c-----------------------------------------------------------------------

      subroutine arrhenius (h, t, arh, arhap, tsurf)

c        Sets temperature-dependent flow coeffs arh,arhap
c        (h-grid). Set for non-ice points in case ice advances in
c        icedyn.

      use comicephys
      use comicegrid

c     passed:
      dimension
     *  h(nx,ny),           t(nx,ny,0:nlevp),
     *  arh(nx,ny,nlev),    arhap(nx,ny),
     *  tsurf(nx,ny)

c     local:
      parameter (
c---------------------
c----
     *    crheoli1 = 2.0 e-16,                       ! warm,  Pa-3 a-1
     *    crheoli2 = 1.66e-16,                       ! cold,  Pa-3 a-1
c-----
c-----

     *    qact1    = 9.545e4,                        ! warm, J/mol
     *    qact2    = 7.820e4,                        ! cold, J/mol
     *    tact     = -6.5,                           ! deg K
     *    rgas     = 8.314,                          ! J/mol/K

c777 *    enhancesheet = 3.,
c777 *    enhanceshelf = 1.)
c    *    enhancesheet = 10.,   ! 777777
c    *    enhanceshelf = 0.1)   ! 777777
c     *    enhancesheet = 5.,
c    *    enhanceshelf = 5.)
cc   *    enhancesheet = 2.,
c    *    enhancesheet = 1.,
cc   *    enhanceshelf = 1.)
     *    enhanceshelf = 0.5)
c    *    enhanceshelf = 5.)

      parameter (zq1 = qact1/rgas,  zq2 = qact2/rgas)

      do j = 1,ny
        do i = 1,nx
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          if (h(i,j).gt.0.) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            arhap(i,j) = 0.
            do k=1,nlev

c.....................
c......................
              ztstar = t(i,j,k) + dtmdh*h(i,j)*zeta(k)
              if (ztstar. ge. 263.15) then
                za = 5.47e10 * exp(-13.9e4/(8.314*ztstar))
              else
                za = 1.14e-5 * exp(-6.0 e4/(8.314*ztstar))
              endif
c....................
c.....

              arh(i,j,k) = enhancesheet*za
c             arhap used only for shelf:
              arhap(i,j) = arhap(i,j)
     *                   + ((enhanceshelf*za)**powir)*dzeta(k)

            enddo

c>>>>>>>>>>>>>
          else
c>>>>>>>>>>>>>

            ztstar = min (tmelt, tsurf(i,j))
c.....................
c......................
            if (ztstar. ge. 263.15) then
              za = 5.47e10 * exp(-13.9e4/(8.314*ztstar))
            else
              za = 1.14e-5 * exp(-6.0 e4/(8.314*ztstar))
            endif
c....................
c.....

            do k=1,nlev
              arh(i,j,k) = enhancesheet * za
            enddo
c           arhap used only for shelf:
            arhap(i,j) = (enhanceshelf*za)**powir
c>>>>>>>>>>>>>>
          endif
c>>>>>>>>>>>>>>

        enddo
      enddo

      return
      end