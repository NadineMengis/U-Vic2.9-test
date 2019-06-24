! source file: /net/mare/home1/eby/as/ism/icedyn.F

c-----------------------------------------------------------------------

      subroutine icedyn (h, hs, hb, hw, sedim, t,
     *                   arhap, s2a0, heatb, budgall,
     *                   maskh, dfu, dfv,
     *                   ua, va, ui, vi, ub, vb, uadv, vadv,
     *                   hu, hv, uw, vw, masku, maskv,
     *                   crhu, crhv, fsedu, fsedv,
     *                   powbu, powbv, fracgu, fracgv,
     *                   muind, mvind, thetau, thetav, hgu, hgv,
     *                   fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *                   numh, itera, iterc, totflow, totneg,
     *                   maskwater, indlake, npoilake, nlake, sealev,
     *                   timeice, dt, ifrest)

c     Steps combined ice sheet-shelf dynamics

      use comicephys
      use comicegrid
      use comicesparse

c     passed:
      dimension
     *  h(nx,ny),              hs(nx,ny),
     *  hb(nx,ny),             hw(nx,ny),
     *  sedim(nx,ny),          t(nx,ny,0:nlevp),  arhap(nx,ny),
     *  s2a0(nx,ny),           heatb(nx,ny),
     *  budgall(nx,ny),        maskh(nx,ny),
     *  dfu(0:nxp,0:nyp),      dfv(0:nxp,0:nyp),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  ui(0:nxp,0:nyp),       vi(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),       vb(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),     vadv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  uw(0:nxp,0:nyp),       vw(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),    maskv(0:nxp,0:nyp),
     *  crhu(0:nxp,0:nyp),     crhv(0:nxp,0:nyp),
     *  fsedu(0:nxp,0:nyp),    fsedv(0:nxp,0:nyp),
     *  powbu(0:nxp,0:nyp),    powbv(0:nxp,0:nyp),
     *  fracgu(0:nxp,0:nyp),   fracgv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp),
     *  thetau(nx,ny),         thetav(nx,ny),
     *  hgu(0:nxp,0:nyp),      hgv(0:nxp,0:nyp),
     *  fluxgrdu(0:nxp,0:nyp), fluxgrdv(0:nxp,0:nyp),
     *  fluxschu(0:nxp,0:nyp), fluxschv(0:nxp,0:nyp),
     *  maskwater(nx,ny),      indlake(npoimax,nlakemax),
     *  npoilake(nlakemax)

c     local:

      dimension
     *  h0(nx,ny),             h1(nx,ny),
     *  visc(nx,ny),
     *  g1(nx,ny),             g2(nx,ny),
     *  coef1(nx,ny),
     *  dhsu(0:nxp,0:nyp),     dhsv(0:nxp,0:nyp),
     *  crhu2(0:nxp,0:nyp),    crhv2(0:nxp,0:nyp),
     *  coefbu(0:nxp,0:nyp),   coefbv(0:nxp,0:nyp)

c     for basal heating calculations:
      dimension
     *  heatbu(0:nxp,0:nyp),   heatbv(0:nxp,0:nyp)

c     cross grid:
      dimension
     *  g2c(0:nx,0:ny),        coef2(0:nx,0:ny),
     *  maskc(0:nx,0:ny)

c     for iteration (a,c) convergence criteria:
      dimension
     *  uaprev (0:nxp,0:nyp),  vaprev (0:nxp,0:nyp),     ! loop a
     *  h1prev(nx,ny)                                    ! loop c

c     for Hindmarsh and Payne (1996) convergence, loop c:
c     dimension delh1(nx,ny), delh1prev(nx,ny)

      parameter (delcrita = .01)   ! m/a  ! .01,.03
      parameter (delcritc = .01)   ! m    ! .01,.03

      parameter (niteramax = 3)
c     parameter (niteramax = 1)

      parameter (nitercmax = 3)
c     parameter (nitercmax = 10)

      character(120) g_st
      real totneg_array(nx,ny), h_old(nx,ny), totneg_land
      integer ntrec
      data ntrec /0/
      save ntrec
c-------------------------
c Start of executable code
c-------------------------

      ubmax =  0.75 * dd0 / dt

      itera = 0
      iterc = 0

c     store original h in h0 (and initial h1):
      call scopy (nx*ny, h, 1, h0, 1)
      call scopy (nx*ny, h, 1, h1, 1)

c        Set masks (maskh, masku, maskv). All masks are set
c        outside iteration loops, ignoring changes of non-zero <->
c        zero ice (h) between iterations.

c        Set maskh (0/1 if no ice/ice, h grid)

      call izero (maskh,  nx*ny)
      do j = 1,ny
        do i = 1,nx
          if (h(i,j).gt.0.) then
            maskh(i,j) = 1
          endif
        enddo
      enddo

c        Set maskc (cross grid)

      call izero (maskc, (nx+1)*(ny+1))
      do j = 1,ny-1
        do i = 1,nx-1
          maskc(i,j)=maskh(i,j)*maskh(i+1,j)*maskh(i,j+1)*maskh(i+1,j+1)
        enddo
      enddo

c        Set masku,v (1 to allow flow, 0 for no flow if no ice to
c        either side, or lower ice than non-ice surface) (u,v grids).
c        (nb: outside loop C, so ignores hs changes within iters).

      call izero (masku, (nxp+1)*(nyp+1))
      do j = 1,ny
        do i = 1,nx-1
          if ( maskh(i,j).eq.1 .and. maskh(i+1,j).eq.1 ) then
            masku(i,j) = 1
          else if ( maskh(i,j).eq.1 .and. maskh(i+1,j).eq.0 .and.
     *              hs(i,j).gt.hs(i+1,j) ) then
            masku(i,j) = 1
          else if ( maskh(i,j).eq.0 .and. maskh(i+1,j).eq.1 .and.
     *              hs(i,j).lt.hs(i+1,j) ) then
            masku(i,j) = 1
          endif
        enddo
      enddo

      call izero (maskv, (nxp+1)*(nyp+1))
      do j = 1,ny-1
        do i = 1,nx
          if ( maskh(i,j).eq.1 .and. maskh(i,j+1).eq.1 ) then
            maskv(i,j) = 1
          else if ( maskh(i,j).eq.1 .and. maskh(i,j+1).eq.0 .and.
     *              hs(i,j).gt.hs(i,j+1) ) then
            maskv(i,j) = 1
          else if (maskh(i,j).eq.0 .and. maskh(i,j+1).eq.1 .and.
     *             hs(i,j).lt.hs(i,j+1) ) then
            maskv(i,j) = 1
          endif
        enddo
      enddo

c        Zero diagnostic quantities (some necessary for finite diffs
c        near ice edges)

      call zero (visc, nx*ny)
      call zero (g1,   nx*ny)
      call zero (g2,   nx*ny)
      call zero (coef1,nx*ny)

      call zero (g2c,  (nx+1)*(ny+1))
      call zero (coef2,(nx+1)*(ny+1))

      call zero (hu,   (nxp+1)*(nyp+1))
      call zero (hv,   (nxp+1)*(nyp+1))
      call zero (dhsu, (nxp+1)*(nyp+1))
      call zero (dhsv, (nxp+1)*(nyp+1))
      call zero (uadv, (nxp+1)*(nyp+1))
      call zero (vadv, (nxp+1)*(nyp+1))

      call zero (dfu,  (nxp+1)*(nyp+1))
      call zero (dfv,  (nxp+1)*(nyp+1))

      do j = 1,ny
        do i = 1,nx-1
          if (masku(i,j).eq.0) then
            ua(i,j)  = 0.
            ui(i,j)  = 0.
            ub(i,j)  = 0.
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (maskv(i,j).eq.0) then
            va(i,j)  = 0.
            vi(i,j)  = 0.
            vb(i,j)  = 0.
          endif
        enddo
      enddo

c        If first timestep of run, zero all velocities. Else, keep
c        for next timestep's initial guess (u,v grids, all passed)
c        (some necessary for finite diffs near ice edges)

      if (firsticedyn)  then
        if (ifrest.ne.1) then
          call zero (ub,   (nxp+1)*(nyp+1))
          call zero (vb,   (nxp+1)*(nyp+1))
          call zero (ui,   (nxp+1)*(nyp+1))
          call zero (vi,   (nxp+1)*(nyp+1))
          call zero (ua,   (nxp+1)*(nyp+1))
          call zero (va,   (nxp+1)*(nyp+1))
        endif
        call resetr (thetau, nx*ny, -1.)
        call resetr (thetav, nx*ny, -1.)
        firsticedyn = .false.
      endif

c        Test for no ice. If so, just apply local forcing
c        and skip all velocity calcs

      numh = 0
      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            numh = numh + 1
          endif
        enddo
      enddo
      if (numh.eq.0) then
        totneg = 0.
        totflow = 0.
        do j=1,ny
          do i=1,nx
            h(i,j) = h(i,j) + budgall(i,j)*dt
            if (h(i,j).lt.0.) then
              totneg = totneg - h(i,j)*darea(i,j)
              h(i,j) = 0.
            endif
          enddo
        enddo
        goto 3500
      endif

c============================
c     Top of loop C (~Picard)
      do iterc = 1,nitercmax
c============================

      call scopy (nx*ny, h1, 1, h1prev, 1)

c        Adjust hw and hs (for change in h):

      call adjustpres (maskwater, indlake, npoilake, nlake,
     *                 h, hb, hw, hs, sealev)

c        Set crh[u,v]2 = crh[u,v] for non-water, =1.e50(flag) for water.
c        Do within loop C after adjustpres to capture interative changes
c        in hw (i.e., pts near g.l.s that switch floating <-> grounded).
c        Then use crh[u,v]2 for rest of icedyn. Arbitrary 1.e50 flag
c        value is used just to set coefb[u,v]=coefbwater below (and
c        in alternative criterion for muind >= 1, ok if large). 1.e50 is
c        never encountered in schoofgl (which only involves g.l. u-pts).
c        In icestep, can be encountered but only for mid-strip water-
c        surrounded pts with hu,dfu=0 (muind=0,masku=0), so * by zeros.

      do j = 1,ny
        do i = 1,nx-1
          if (hw(i,j).gt.hwcut .and. hw(i+1,j).gt.hwcut) then
            crhu2(i,j) = 1.e50
          else
            crhu2(i,j) = crhu(i,j)
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (hw(i,j).gt.hwcut .and. hw(i,j+1).gt.hwcut) then
            crhv2(i,j) = 1.e50
          else
            crhv2(i,j) = crhv(i,j)
          endif
        enddo
      enddo

c       Set hu, hv, dhsu, dhsv, dfu, dfv  (u,v grids)
c       Do inside loop C to capture hs changes within iters.

      do j = 1,ny
        do i = 1,nx-1
          if (masku(i,j).eq.1) then
            hu(i,j) = 0.5*(h(i,j) + h(i+1,j))
            dhsu(i,j) = (hs(i+1,j) - hs(i,j)) / dxu(i,j)
            dfu(i,j)  = rhoice*grav*hu(i,j)*dhsu(i,j)
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (maskv(i,j).eq.1) then
            hv(i,j) = 0.5*(h(i,j) + h(i,j+1))
            dhsv(i,j) = (hs(i,j+1) - hs(i,j)) / dyv(i,j)
            dfv(i,j)  = rhoice*grav*hv(i,j)*dhsv(i,j)
          endif
        enddo
      enddo

      call izero (muind, (nxp+1)*(nyp+1))
      call izero (mvind, (nxp+1)*(nyp+1))

c        Find locations of grounding lines (on u,v grids), doing
c        Schoof-based calculations. Sets [u,v]a, m[u,v]ind = -1,
c        and diagnostics at those locations.

      call schoofgl (h, hs, hb, hw, arhap, maskwater, s2a0,
     *               dfu, dfv, hu, hv,
     *               crhu2, crhv2, powbu, powbv,
     *               ua, va, masku, maskv,
     *               fracgu, fracgv, thetau, thetav, hgu, hgv,
     *               fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *               muind, mvind, sealev, timeice, ubmax)

c        Set nuvtot and muind,mvind, indices from 2-D arrays into
c        linear sequence of selected points for tridia or sparse solns
c        below. Don't pick points that have schoof constraints
c        (with m[u,v]ind set to -1 in schoofgl).
c        (Can do outside A loop, depends only on mask[u,v],hw).

      nuvtot = 0
c----------------------
c----------------------
      do j = 1,ny
        do i = 1,nx
          if (masku(i,j).eq.1 .and. muind(i,j).ne.-1) then
c           Either: include if both adjacent h-points have water:
            if (hw(i,j).gt.hwcut. and. hw(min(i+1,nx),j).gt.hwcut) then
c777        Or: include if crhu2 is slippery enough (1.e-8 is "~stream"
c777        value for powb=2, see basecoef comments, 1.e5="taunorm"):
c777        if (crhu2(i,j).gt. 1.e-8*(1.e5**2)/(1.e5**powbu(i,j)) then
              nuvtot = nuvtot + 1
              muind(i,j) = nuvtot
            endif
          endif

          if (maskv(i,j).eq.1 .and. mvind(i,j).ne.-1) then
            if (hw(i,j).gt.hwcut. and. hw(i,min(j+1,ny)).gt.hwcut) then
c777        if (crhv2(i,j).gt. 1.e-8*(1.e5**2)/(1.e5**powbv(i,j)) then
              nuvtot = nuvtot + 1
              mvind(i,j) = nuvtot
            endif
          endif
        enddo
      enddo
c-----
c-----

c===============================
c     Top of loop A (shelf flow)
      do itera = 1,niteramax
c===============================

      call scopy ((nxp+1)*(nyp+1), ua, 1, uaprev, 1)
      call scopy ((nxp+1)*(nyp+1), va, 1, vaprev, 1)

c>>>>>>>>>>>>>>>>>>>>>>>>>>
      if (nuvtot.gt.0) then ! if no shelf-like pts, skip all shelf calcs
c>>>>>>>>>>>>>>>>>>>>>>>>>>

c        Set visc (h grid)

c        Set g1 ((du/dx)**2 + (dv/dy)**2 + (du/dx)*(dv/dy)) (h grid)

      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            zdudx = (ua(i,j)-ua(i-1,j))/dx(i,j)
            zdvdy = (va(i,j)-va(i,j-1))/dy(i,j)
            g1(i,j) = zdudx**2 + zdvdy**2 + zdudx*zdvdy
          endif
        enddo
      enddo

c        Set g2 ([du/dy+dv/dx]**2}/4) (h grid, first g2c cross gd)

      do j = 1,ny-1
        do i = 1,nx-1
          if (maskc(i,j).ne.0) then
            zdudy = (ua(i,j+1)-ua(i,j))/dyc(i,j)
            zdvdx = (va(i+1,j)-va(i,j))/dxc(i,j)
            g2c(i,j) = 0.25*((zdudy + zdvdx)**2)
          endif
        enddo
      enddo

      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            g2(i,j) = 0.25*(g2c(i-1,j-1)+g2c(i,j-1)+g2c(i-1,j)+g2c(i,j))
          endif
        enddo
      enddo

c        Set visc (h grid)

      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            visc(i,j) = 0.5
     *                / (max(g1(i,j)+g2(i,j), gmin) ** powiv)
          endif
        enddo
      enddo

c        Set coef1 (2*visc*h/(A**1/n)) for ice-shelf eqns (h grid)

      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
c           coef1(i,j) = 2.*visc(i,j)*h(i,j)          / arhap(i,j)
            coef1(i,j) = 2.*visc(i,j)*max(h(i,j),.01) / arhap(i,j)
          endif
        enddo
      enddo

c        Set coef2 (visc*h/(A**1/n)) for ice-shelf eqns (cross grid)

      do j = 1,ny-1
        do i = 1,nx-1
          if (maskc(i,j).ne.0) then
            coef2(i,j) = 0.5 * 0.25 * (   coef1(i,j)   + coef1(i+1,j)
     *                                  + coef1(i,j+1) + coef1(i+1,j+1))
          endif
        enddo
      enddo

c        Set linearized basal sliding coeffs coefb[u,v] for shelf eqns
c        = (crhu2**(-1/m)) (ub2**((1-m)/2m)) (u,v grids)

      do j = 1,ny
        do i = 1,nx-1
          if (masku(i,j).eq.1) then
            if (crhu2(i,j).eq.1.e50) then    ! flag for water, set above
              coefbu(i,j) = coefbwater
            else
              zc = crhu2(i,j) ** (-1./powbu(i,j))
              zpow = -0.5*(1. - 1./powbu(i,j))
              zu = ub(i,j)
              zv = 0.25* (vb(i,j) + vb(i+1,j) + vb(i,j-1) + vb(i+1,j-1))
              coefbu(i,j) = zc * ( max(zu**2+zv**2, ubmin**2)**zpow )
              coefbu(i,j) = max (coefbu(i,j), coefbwater)
            endif
          else
            coefbu(i,j) = 0. ! not used
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (maskv(i,j).eq.1) then
            if (crhv2(i,j).eq.1.e50) then    ! flag for water, set above
              coefbv(i,j) = coefbwater
            else
              zc = crhv2(i,j) ** (-1./powbv(i,j))
              zpow = -0.5*(1. - 1./powbv(i,j))
              zu = 0.25* (ub(i,j) + ub(i,j+1) + ub(i-1,j) + ub(i-1,j+1))
              zv = vb(i,j)
              coefbv(i,j) = zc * ( max(zu**2+zv**2, ubmin**2)**zpow )
              coefbv(i,j) = max (coefbv(i,j), coefbwater)
            endif
          else
            coefbv(i,j) = 0. ! not used
          endif
        enddo
      enddo

c        Tridia or Sparse matrix solution of ice-shelf eqns (for ua,va)

      if (ny.eq.1) then
        call dotridia (ua, hu, dhsu, coef1, coefbu, uw, muind, nuvtot)
      else
        call dosparse (ua, va, hu, hv, dhsu, dhsv,
     *                 coef1, coef2, coefbu, coefbv, uw, vw,
     *                 muind, mvind, maskwater, h, hs, hb, hw,
     *                 timeice, itera, iterc)
      endif

      do j=1,ny-1
        do i=1,nx-1
          ua(i,j) = max (-ubmax, min (ubmax, ua(i,j)))
          va(i,j) = max (-ubmax, min (ubmax, va(i,j)))
        enddo
      enddo

c        thetau,v = 1 - fraction of g.l. longitudinal stress buttressed,
c        used for g.l. locations in next calcgl (called from schoofgl),
c        or just diagnostic if not schoofgl.

      call thetacalc (h, maskwater, arhap, muind, mvind, ua, va,
     *                thetau, thetav)

c>>>>>>>>>>
      endif ! nuvtot > 0 test
c>>>>>>>>>>

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c~~~~~
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c~~~~~

c        Test for loop A convergence (on depth-average velocs ua,va)

      delmaxa = 0.
      do j=0,nyp
        do i=0,nxp
          delmaxa = max (delmaxa, abs(ua(i,j)-uaprev(i,j)),
     *                            abs(va(i,j)-vaprev(i,j)))
        enddo
      enddo

c     monitor loop a convergence:
c     do iu = ioterm, iuout1d, iuout1d-ioterm
c       write (iu,'(a,i3,2(a,f10.4))')
c    *    'itera=',itera, '  delmaxa=',delmaxa,'  delcrita=',delcrita
c     enddo

      if (delmaxa.lt.delcrita) goto 1000

c=======================
c     End loop A (itera)
      enddo
      itera = niteramax
 1000 continue
c=======================

 2000 continue

c        Step through one timestep to get h1, with time-implicit
c        internal flow, and implicit h * explicit ub .
c        Also return internal flow quantities ui, vi

cm    call mascalc (h0, ztoth0, ztota0)

      call icestep (h0,  h,   h1,    hw,
     *              hu,  hv,  dhsu,  dhsv,  dfu,  dfv,
     *              s2a0, crhu2, crhv2, powbu, powbv,
     *              ui, vi, ub, vb, ua, va,
     *              masku, maskv, muind, mvind, dt)

c        Explicitly apply surface mass balance to ice thickness
      h_old(:,:)=h1(:,:)
      do j=1,ny
        do i=1,nx
          h1(i,j) = h1(i,j) + budgall(i,j)*dt
        enddo
      enddo

cm    call mascalc (h1, ztoth1a, ztota1a)
      totneg_array(:,:)=0.
      totneg_land=0.
      do j=1,ny
        do i=1,nx
          if (h1(i,j).lt.0. .and. h_old(i,j) .gt. 1.e-2) then
	    totneg_array(i,j)=h1(i,j)
	    totneg_land=totneg_land+h1(i,j)*darea(i,j)
          endif
        enddo
      enddo
!       ntrec=ntrec+1
!       g_st='totneg_array'
!       call jerncdf_ice_timedep(totneg_array,g_st,g_st,ntrec)
!       g_st='h_old'
!       call jerncdf_ice_timedep(h_old,g_st,g_st,ntrec)
!       g_st='budgall'
!       call jerncdf_ice_timedep(budgall,g_st,g_st,ntrec)
c        Correct for negative h1

      totneg = 0.
      do j=1,ny
        do i=1,nx
          if (h1(i,j).lt.0.) then
            totneg = totneg - h1(i,j)*darea(i,j)
            h1(i,j) = 0.
          endif
        enddo
      enddo

c        Impose boundary constraints,etc, on h for various expts

      totflow = 0.
      if (nx.gt.1) then
        do j=1,ny
          totflow = totflow - h1(1,j) *darea(1,j)
          totflow = totflow - h1(nx,j)*darea(nx,j)
          h1(1, j) = 0.
          h1(nx,j) = 0.
        enddo
      endif
      if (ny.gt.1) then
        do i=1,nx
          totflow = totflow - h1(i,1) *darea(i,1)
          totflow = totflow - h1(i,ny)*darea(i,ny)
          h1(i, 1) = 0.
          h1(i,ny) = 0.
        enddo
      endif

cm    call mascalc (h1, ztoth1b, ztota1b)
cm    write (6,'(a,4f20.14)')
cm   *    '   toth1b, 1b-1a, totflow, totneg:',
cm   *     ztoth1b/totarea, (ztoth1b-ztoth1a)/totarea,
cm   *     totflow/totarea, totneg/totarea

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Hindmarsh and Payne(1996) convergence:
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     do j=1,ny
c       do i=1,nx
c         delh1(i,j) = h1(i,j) - h1prev(i,j)
c       enddo
c     enddo
c
c     zcos = 999.
c     zthed= 999.
c     zrat = 999.
c     if (iterc.ge.2) then
c       znorma = 0.
c       znormb = 0.
c       znormc = 0.
c       zdot   = 0.
c       do j=1,ny
c         do i=1,nx
c           if (h1(i,j).gt.0.) then
c             znorma = znorma + delh1prev(i,j)**2
c             znormb = znormb + delh1(i,j)**2
c             znormc = znormc + (delh1(i,j)-delh1prev(i,j))**2
c             zdot   = zdot   + delh1(i,j)*delh1prev(i,j)
c           endif
c         enddo
c       enddo
c       znorma = sqrt(znorma)
c       znormb = sqrt(znormb)
c       znormc = sqrt(znormc)
c       if (znorma.gt.0. .and. znormb.gt.0. .and. znormc.gt.0.) then
c         zcos = zdot/(znorma*znormb)
c         if (abs(zcos).le.1.) then
c           zthe = acos (zcos)
c           zthed= zthe*180./pi
c           if (abs(zthed).lt.30. .or. abs(zthed).gt.180.-30.) then
cc          if (abs(zthed).lt.45. .or. abs(zthed).gt.180.-45.) then
c             zrat = znorma/znormc
c             do j=1,ny
c               do i=1,nx
c                 h1(i,j) = h1prev(i,j) + delh1(i,j)*zrat
c               enddo
c             enddo
c           endif
c         endif
c       endif
c     endif
c
c     call scopy (nx*ny, delh1, 1, delh1prev, 1)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     End of Hindmarsh and Payne convergence
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c        Set uadv,vadv fluxes for diagnostics (iceshow).

      do j = 1,ny
        do i = 1,nx-1
          if (masku(i,j).eq.1) then
c           centered for sheet ua, upstream for shelf ua:
            if (muind(i,j).eq.0) then
              uadv(i,j) = ua(i,j)*hu(i,j)*dyu(i,j)
            else
              uadv(i,j) =
     *            (  max (ua(i,j),0.)*h(i,j) + min (ua(i,j),0.)*h(i+1,j)
     *             ) * dyu(i,j)
            endif
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (maskv(i,j).eq.1) then
c           centered for sheet va, upstream for shelf va:
            if (mvind(i,j).eq.0) then
              vadv(i,j) = va(i,j)*hv(i,j)*dxv(i,j)
            else
              vadv(i,j) =
     *             ( max (va(i,j),0.)*h(i,j) + min (va(i,j),0.)*h(i,j+1)
     *             ) * dxv(i,j)
            endif
          endif
        enddo
      enddo

c        Update h at t+0.5*dt to be used for next iter's viscs

      do j=1,ny
        do i=1,nx
          h(i,j) = 0.5*(h1(i,j) + h0(i,j))
        enddo
      enddo

c        Test for loop C convergence (on ice thickness h1)

      delmaxc = 0.
      do j=1,ny
        do i=1,nx
          delmaxc = max (delmaxc, abs(h1(i,j)-h1prev(i,j)))
        enddo
      enddo

c     monitor loop c convergence:
c     do iu = ioterm, iuout1d, iuout1d-ioterm
c       write (iu,'(a,i3,2(a,f10.4))')
c    *    'iterc=',iterc, '  delmaxc=',delmaxc,'  delcritc=',delcritc
cc      for Hindmarsh and Payne convergence:
cc      write (iu,'(9x,3(a,f20.8))')
cc   *    '  zcos=',zcos, '  zthe=',zthed, '  zrat=',zrat
c     enddo

      if (delmaxc.lt.delcritc) goto 3000

c=======================
c     End loop C (iterc)
      enddo
      iterc = nitercmax
 3000 continue
c=======================

c        Set new h for this timestep (at t+dt)

      call scopy (nx*ny, h1, 1, h, 1)

c        Calculate basal heating (heatb, J/m2/y) due to basal sliding,
c        equal to taub_u*ub + taub_v*vb. First on u-grid (heatbu)
c        and v-grid (heatbv), then transfer to h-grid (heatv).
c        heatb is passed to icetherm.
c
c        If defined sediment, heatb has been calculated in sedflow
c        (using last time step quantities), allowing for tau-split
c        between sediment and bare rock.

c     ub:
      call zero (heatbu, (nxp+1)*(nyp+1))
      do j = 1,ny
        do i = 1,nx-1
          if ( masku(i,j).eq.1 .and.
     *         .not. (hw(i,j).gt.hwcut. and. hw(i+1,j).gt.hwcut)
     *         .and. muind(i,j).ne.-1     ! added 9/10: avoid Schoof pt
     *       ) then
            heatbu(i,j) = -(1.-fsedu(i,j))*dfu(i,j)*ub(i,j)
          endif
        enddo
      enddo

c     vb:
      call zero (heatbv, (nxp+1)*(nyp+1))
      do j = 1,ny-1
        do i = 1,nx
          if ( maskv(i,j).eq.1 .and.
     *         .not. (hw(i,j).gt.hwcut. and. hw(i,j+1).gt.hwcut)
     *         .and. mvind(i,j).ne.-1     ! added 9/10: avoid Schoof pt
     *       ) then
            heatbv(i,j) = -(1.-fsedv(i,j))*dfv(i,j)*vb(i,j)
          endif
        enddo
      enddo

c     transfer to h-grid:
      call zero (heatb, nx*ny)
      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            heatb(i,j) = 0.5*(heatbu(i-1,j) + heatbu(i,j))
     *                 + 0.5*(heatbv(i,j-1) + heatbv(i,j))
c           heatb(i,j) = 0.
          endif
        enddo
      enddo

c        Skip point for no ice (numh=0, top of loop C)

 3500 continue

c          Final adjust of hw and hs (for change in h):

      call adjustpres (maskwater, indlake, npoilake, nlake,
     *                 h, hb, hw, hs, sealev)

      return
      end

c-----------------------------------------------------------------------

      subroutine dotridia (ua, hu, dhsu, coef1, coefbu, uw,muind,nuvtot)

c     Assembles *tri arrays for 1-D tridiagonal solution (called
c     if ny = 1), calls tridia_i, sets new velocities ua

      use comicephys
      use comicegrid

      dimension
     *  ua(0:nxp,0:nyp),       hu(0:nxp,0:nyp),
     *  dhsu(0:nxp,0:nyp),     coef1(nx,ny),       coefbu(0:nxp,0:nyp),
     *  uw(0:nxp,0:nyp),       muind(0:nxp,0:nyp)

      dimension atri(nx), btri(nx), ctri(nx), rtri(nx), vtri(nx)

c        Tridiagonal solution for ua in 1-D problems (in x, with ny=1).
c        Relies on muind running left to right (gaps in ice are ok).

      j = 1

      call zero (atri, nx)
      call zero (btri, nx)
      call zero (ctri, nx)
      call zero (rtri, nx)
      call zero (vtri, nx)

      do i=1,nx-1
        n1 = muind(i,j)

        if (n1.ge.1) then
          rtri(n1) = rhoice*grav*hu(i,j)*dhsu(i,j) - coefbu(i,j)*uw(i,j)

          z1 = 2.*coef1(i+1,j) / (dx(i+1,j)*dxu(i,j))  ! ua(i+1,j)
          z2 = 2.*coef1(i  ,j) / (dx(i,j  )*dxu(i,j))  ! ua(i-1,j)

c         diagonal:
          atri(n1) = - (z1 + z2) - coefbu(i,j) - sidefac*hu(i,j)

c         superdiagonal:
          if (muind(i+1,j).ge.1) then
            btri(n1) = z1
          else
            rtri(n1) = rtri(n1) - z1*ua(i+1,j)
          endif

c         subdiagonal:
          if (muind(i-1,j).ge.1) then
            ctri(n1) = z2
          else
            rtri(n1) = rtri(n1) - z2*ua(i-1,j)
          endif
        endif

      enddo

      call tridia_i (atri, btri, ctri, vtri, rtri, nuvtot, 0)

c        Set new ua

      do i=1,nx
        n1 = muind(i,j)
        if (n1.ge.1) then
          ua(i,j) = vtri(n1)
        endif
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine dosparse (ua, va, hu, hv, dhsu, dhsv,
     *                     coef1, coef2, coefbu, coefbv, uw, vw,
     *                     muind, mvind, maskwater, h, hs, hb, hw,
     *                     timeice, itera, iterc)

c     Solevs elliptic shelf equations for ua,va, using SOR
c     (Successive Over Relaxation),  just for points with
c     m[u,v]ind >= 1.

      use comicephys
      use comicegrid
      use comicesparse

      dimension
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  dhsu(0:nxp,0:nyp),     dhsv(0:nxp,0:nyp),
     *  coef1(nx,ny),          coef2(0:nx,0:ny),
     *  coefbu(0:nxp,0:nyp),   coefbv(0:nxp,0:nyp),
     *  uw(0:nxp,0:nyp),       vw(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp),
     *  maskwater(nx,ny),
     *  h(nx,ny),              hs(nx,ny),
     *  hb(nx,ny),             hw(nx,ny)

c     parameter (nitersurmax =   10, facsor = 0.5) ! a
c     parameter (nitersurmax =   10, facsor = 0.9) ! b
c     parameter (nitersurmax =   10, facsor = 1.2) ! c
c     parameter (nitersurmax =  100, facsor = 1.2)
c     parameter (nitersurmax =  100, facsor = 0.9)
c     parameter (nitersurmax =  500, nitersurmin = 20, facsor = 0.9)
      parameter (nitersurmax =  100, nitersurmin = 20, facsor = 1.2)

c     for diagnostic dump only:
c     dimension
c    *  uaorig(0:nxp,0:nyp),       vaorig(0:nxp,0:nyp)

c     for checkerboard sor:

      if (firstdosparse) then
        do j=0,ny+1
          do i=0,nx+1
            mcheck(i,j) = mod (i+j,2)
          enddo
        enddo
        firstdosparse = .false.
      endif

c     call scopy ((nxp+1)*(nyp+1), ua, 1, uaorig, 1)
c     call scopy ((nxp+1)*(nyp+1), va, 1, vaorig, 1)

c        Calculate mean veloc, used to test skip out of iteration:

      zuav = 0.
      ndav = 0.
      do j=1,ny
        do i=1,nx
          if (muind(i,j).ge.1) then
            zuav = zuav + ua(i,j)**2
            ndav = ndav + 1
          endif
          if (mvind(i,j).ge.1) then
            zuav = zuav + va(i,j)**2
            ndav = ndav + 1
          endif
        enddo
      enddo
      zuav = sqrt (zuav/ max(ndav,1))

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do itersur = 1,nitersurmax
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      zdumax =  0.

      do icheck = 0,1    ! for checkerboard sor updating
      do j=1,ny
        do i=1,nx-1

          if (muind(i,j).ge.1 .and. mcheck(i,j).eq.icheck) then

c              u-neighbors for ua(i,j):

            z1 = 2.*coef1(i+1,j) / (dx(i+1,j )*dxu(i,j))  ! i+1,j; i,j
            z2 = 2.*coef1(i  ,j) / (dx(i,j   )*dxu(i,j))  ! i,j  ; i-1,j
            z3 =    coef2(i  ,j) / (dyc(i,j  )*dyu(i,j))  ! i,j+1; i,j
            z4 =    coef2(i,j-1) / (dyc(i,j-1)*dyu(i,j))  ! i,j  ; i,j-1

c              v-neighbors for ua(i,j):

            za = coef1(i+1,j) / (dy(i+1,j )*dxu(i,j)) ! i+1,j; i+1,j-1
            zb = coef1(i  ,j) / (dy(i,j   )*dxu(i,j)) ! i,j  ; i,j-1
            zc = coef2(i  ,j) / (dxc(i,j  )*dyu(i,j)) ! i+1,j; i,j
            zd = coef2(i,j-1) / (dxc(i,j-1)*dyu(i,j)) ! i,j-1; i+1,j-1

            zforcu =
     *          z1*(ua(i+1,j)        )     - z2*(       -ua(i-1,j))
     *        + z3*(ua(i,j+1)        )     - z4*(       -ua(i,j-1))

     *        + za*(va(i+1,j)-va(i+1,j-1)) - zb*(va(i,j)-va(i,j-1))
     *        + zc*(va(i+1,j)-va(i,j))     - zd*(va(i+1,j-1)-va(i,j-1))

     *        - rhoice*grav*hu(i,j)*dhsu(i,j) + coefbu(i,j)*uw(i,j)

            zclhsu = z1 + z2 + z3 + z4 + coefbu(i,j) + sidefac*hu(i,j)

            zdu = zforcu/zclhsu - ua(i,j)
            ua(i,j) = facsor*zforcu/zclhsu + (1.-facsor)*ua(i,j)
            zdumax = max (zdumax, abs(zdu))

c           if (abs(zdu).gt.zdumax) then
c             idumax = i
c             jdumax = j
c             zdumax = abs(zdu)
c           endif
c           zduav = zduav + zdu**2
c           zuav  = zuav  + (zforcu/zclhsu)**2
c           ndav  = ndav  + 1

          endif

        enddo
      enddo
      enddo

      do icheck = 0,1    ! for checkerboard sor updating
      do j=1,ny-1
        do i=1,nx

          if (mvind(i,j).ge.1 .and. mcheck(i,j).eq.icheck) then

c              v-neighbors for va(i,j):

            z1 = 2.*coef1(i,j+1) / (dy(i,j+1 )*dyv(i,j))  ! va(i,j+1)
            z2 = 2.*coef1(i  ,j) / (dy(i,j   )*dyv(i,j))  ! va(i,j-1)
            z3 =    coef2(i  ,j) / (dxc(i,j  )*dxv(i,j))  ! va(i+1,j)
            z4 =    coef2(i-1,j) / (dxc(i-1,j)*dxv(i,j))  ! va(i-1,j)

c              u-neighbors for va(i,j):

            za = coef1(i,j+1) / (dx(i,j+1 )*dyv(i,j)) ! i,j+1; i-1,j+1
            zb = coef1(i  ,j) / (dx(i,j   )*dyv(i,j)) ! i,j  ; i-1,j
            zc = coef2(i  ,j) / (dyc(i,j  )*dxv(i,j)) ! i,j+1; i,j
            zd = coef2(i-1,j) / (dyc(i-1,j)*dxv(i,j)) ! i-1,j; i-1,j+1

            zforcv =
     *          z1*(va(i,j+1)        )     - z2*(       -va(i,j-1))
     *        + z3*(va(i+1,j)        )     - z4*(       -va(i-1,j))

     *        + za*(ua(i,j+1)-ua(i-1,j+1)) - zb*(ua(i,j)-ua(i-1,j))
     *        + zc*(ua(i,j+1)-ua(i,j))     - zd*(ua(i-1,j+1)-ua(i-1,j))

     *        - rhoice*grav*hv(i,j)*dhsv(i,j) + coefbv(i,j)*vw(i,j)

            zclhsv = z1 + z2 + z3 + z4 + coefbv(i,j) + sidefac*hv(i,j)

            zdv = zforcv/zclhsv - va(i,j)
            va(i,j) = facsor*zforcv/zclhsv + (1.-facsor)*va(i,j)
            zdumax = max (zdumax, abs(zdv))

c           if (abs(zdv).gt.zdumax) then
c             idumax = i
c             jdumax = j
c             zdumax = abs(zdv)
c           endif
c           zduav = zduav + zdv**2
c           zuav  = zuav  + (zforcv/zclhsv)**2
c           ndav  = ndav  + 1

          endif

        enddo
      enddo
      enddo

c     zduav = sqrt(zduav/max(ndav,1))
c     zuav  = sqrt(zuav /max(ndav,1))

c        Print out sor diagnostics

c     if (iterc.eq.1 .and. itera.eq.1 .and. itersur.eq.1 ) then
c       write (6,'(/25("=")/a,i10/25("="))')
c    *     'dosparse: time=', nint(timeice)
c     endif

c     if (itersur.eq.1) then
c       write (6,'(/a,3i8,f12.3/a)')
c    *   'iterc, itera, ndav, zuav=', iterc, itera, ndav, zuav,
c    *   '  itersur      zdumax'
c     endif

c     if (itersur.eq.1 .or. itersur.eq.nitersurmax .or.
c    *   (zdumax.lt..003*zuav .and. itersur.ge.nitersurmin)) then
c       write (6,'(i9, f12.3)') itersur, zdumax
c     endif

c##################
c     diagnostic dump of maps around a given point:
c
c     if (.true.) then
c     if (.false.) then
c
c     io1 = 38 ! 25
c     jo1 = 73 ! 85
c     io2 = io1 + 4
c     jo2 = jo1 + 4
c     iu = 110
c
c     if (iterc.eq.1 .and. itera.eq.1 .and. itersur.eq.1 ) then
c       write (iu,'(/12("=")/a,i8/12("="))') 'time=', nint(timeice)
c     endif
c
c     if (itersur.eq.1 ) then
c
c       write (iu,'(/31("-")/a,2i8)') 'iterc, itera = ', iterc, itera
c
c       write (iu,'(/a,i6)') 'maskwater, h:', itersur
c       do jo=jo2,jo1,-1
c         write (iu,'(5i9, 10x, 5f9.2)')
c    *      (maskwater(io,jo),io=io1,io2), (h(io,jo),io=io1,io2)
c       enddo
c
c       write (iu,'(/a,i6)') 'hb, hs, hw:', itersur
c       do jo=jo2,jo1,-1
c         write (iu,'(5f9.2, 10x, 5f9.2, 10x, 5e11.3)')
c    *      (hb(io,jo),io=io1,io2), (hs(io,jo),io=io1,io2),
c    *      (hw(io,jo),io=io1,io2)
c       enddo
c
c       write (iu,'(/a,i6)') 'muind, mvind:', itersur
c       do jo=jo2,jo1,-1
c         write (iu,'(5i9, 10x, 5i9)')
c    *      (muind(io,jo),io=io1,io2), (mvind(io,jo),io=io1,io2)
c       enddo
c
c       write (iu,'(/a,i6)') 'uaorig, vaorig:', itersur
c       do jo=jo2,jo1,-1
c         write (iu,'(5f9.1, 10x, 5f9.1)')
c    *      (uaorig(io,jo),io=io1,io2), (vaorig(io,jo),io=io1,io2)
c       enddo
c
c     endif
c
c     if (itersur.eq.1 .or. itersur.eq.nitersurmax .or.
c    *   (zdumax.lt..003*zuav .and. itersur.ge.nitersurmin)) then
c     if (.true.) then
c
c       write (iu,'(/a,i6)') 'ua, va:', itersur
c       do jo=jo2,jo1,-1
c         write (iu,'(5f9.1, 10x, 5f9.1)')
c    *      (ua(io,jo),io=io1,io2), (va(io,jo),io=io1,io2)
c       enddo
c
c     endif
c
c     if (zdumax.gt.100000.) then
c       write (6,'(/a)') 'STOP...large zdumax'
c       stop
c     endif
c
c     endif
c##################

c     skip out of sor iteration if max change is small enough:
      if (zdumax.lt..003*zuav .and. itersur.ge.nitersurmin) goto 1000

c>>>>>>>>>>
      enddo ! itersur
c>>>>>>>>>>
 1000 continue

      return
      end

c-----------------------------------------------------------------------

c>>>>>>>>>>>>>>>>>>>>>
c>>>>>>>>>>>>>>>>>>>>>

      subroutine schoofgl (h, hs, hb, hw, arhap, maskwater, s2a0,
     *                     dfu, dfv, hu, hv,
     *                     crhu2, crhv2, powbu, powbv,
     *                     ua, va, masku, maskv,
     *                     fracgu, fracgv, thetau, thetav, hgu, hgv,
     *                     fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *                     muind, mvind, sealev, timeice
     *,			   ubmax)

c     Does Schoof-based calculations for grounding lines (Schoof, JGR,
c     2007, and my 9/27/07 notes).
c     (1) Finds g.l. locations on u grid.
c     (2) Calculates hg, ug, at and upstream of gl (see 9/27/07
c         and 11/10/07 notes) (in subr calcgl).
c     (3) Sets ua,va, velocs at or downstream of g.l., and sets
c         m[u,v]ind = -1 there so ua,va will be
c         (i)  prescribed b.c.s. in dotridia/dosparse, and
c         (ii) treated as shelf velocs in icestep.
c     (4) Also sets diagnostics flux[grd,sch][u,v].

      use comicephys
      use comicegrid

      dimension
     *  h(nx,ny),              hs(nx,ny),
     *  hw(nx,ny),             hb(nx,ny),
     *  arhap(nx,ny),          maskwater(nx,ny),      s2a0(nx,ny),
     *  dfu(0:nxp,0:nyp),      dfv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  crhu2(0:nxp,0:nyp),    crhv2(0:nxp,0:nyp),
     *  powbu(0:nxp,0:nyp),    powbv(0:nxp,0:nyp),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),    maskv(0:nxp,0:nyp),
     *  fracgu(0:nxp,0:nyp),   fracgv(0:nxp,0:nyp),
     *  thetau(nx,ny),         thetav(nx,ny),
     *  hgu(0:nxp,0:nyp),      hgv(0:nxp,0:nyp),
     *  fluxgrdu(0:nxp,0:nyp), fluxgrdv(0:nxp,0:nyp),
     *  fluxschu(0:nxp,0:nyp), fluxschv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp)

      dimension ifglu(0:nxp), ifglv(0:nyp)

      parameter (secpy = 365.*86400.)

      call zero (hgu,      (nxp+1)*(nyp+1))
      call zero (hgv,      (nxp+1)*(nyp+1))
      call zero (fluxgrdu, (nxp+1)*(nyp+1))
      call zero (fluxgrdv, (nxp+1)*(nyp+1))
      call zero (fluxschu, (nxp+1)*(nyp+1))
      call zero (fluxschv, (nxp+1)*(nyp+1))

c        u grid:

c==============
      do j=1,ny
c==============

c          Use ifglu to require at least 3 h-grid floating/water pts
c          (or 2 if at domain edge) next to g.l., otherwise nearby
c          g.l. calcs could impose fluxes twice at same u-pt (in
c          advancing case).

c          If advancing, procedure imposes veloc (ua,va) at 1st u-pt
c          *downstream* of g.l. u-pt. So to test on buttressing,
c          use theta[u,v] from closest floating h cell to g.l.
c          depending on whether u-pts were included in prev
c          elliptical soln or not (muind, mvind).

        call izero (ifglu, nxp+1)
        do i=2,nx-2
          ip3 = min (i+3,nx)
          im2 = max (i-2,1)
          if (masku(i,j).eq.1) then
            if      ( h(i  ,j).gt.0. .and. maskwater(i  ,j).eq.0
     *                               .and. maskwater(i+1,j).ne.0
     *                               .and. maskwater(i+2,j).ne.0
     *                               .and. maskwater(ip3,j).ne.0 ) then
              ifglu(i) = 1
            else if ( h(i+1,j).gt.0. .and. maskwater(i+1,j).eq.0
     *                               .and. maskwater(i  ,j).ne.0
     *                               .and. maskwater(i-1,j).ne.0
     *                               .and. maskwater(im2,j).ne.0 ) then
              ifglu(i) = -1
            endif
          endif
        enddo

c>>>>>>>>>>>>>>>>>>
        do i=2,nx-2
c>>>>>>>>>>>>>>>>>>
          if (iabs(ifglu(i)).eq.1) then

c              Set indices for points around g.l. On h grid,
c              (ia,j) is grounded with ice, (ib,j) is downstream
c              ocean or *open* lake. On u grid, (iau,j) is ~g.l.,
c              (ibu,j) is next downstream pt.

            if (ifglu(i).eq.1) then
              ia  = i
              ib  = i+1
              ic =  i+2
              iau = i
              ibu = i+1
              idir = 1
            else
              ia  = i+1
              ib  = i
              ic  = i-1
              iau = i
              ibu = i-1
              idir = -1
            endif
c           Use closest floating h-box with a valid thetau (from
c           previous iter...see thetacalc):
            if (thetau(ib,j).ne.-1.) then
              ztheta = thetau(ib,j)
            else if (thetau(ic,j).ne.-1.) then
              ztheta = thetau(ic,j)
            else
              ztheta = -1.
            endif

c              Only do if marine ice is not abutting land,
c              and g.l. is  not completely buttressed

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (hs(ia,j).gt.hs(ib,j) .and. ztheta.gt.0.) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c                Set needed Schoof-like values for calcgl:

              pschm = 1./powbu(iau,j)
              zasch = 0.5*(arhap(ia,j)**powi+arhap(ib,j)**powi)
              zbsch = (1./crhu2(iau,j))**(1./powbu(iau,j))
              zpowb = powbu(iau,j)
              zfrac = fracgu(iau,j)
              zdx = 0.5*(dx(ia,j) + dx(ib,j))

c                Calculates hg, ug (g.l. ice thickness, schoof veloc)

              call calcgl ( h(ia,j), hs(ia,j), hb(ia,j),
     *                      h(ib,j), hs(ib,j), hb(ib,j),
     *                      zasch, zbsch, zpowb, zfrac, zdx,
     *                      sealev, ztheta, hg, ug)

c                Calculate Non-schoof "grid" (shallow-ice) flux at g.l.
c                (like expliclit zflux's in icestep)

              zdfv  = 0.25
     *              * ( dfv(ia,j) +dfv(ib,j) +dfv(ia,j-1) +dfv(ib,j-1))
              zdfu  = dfu(iau,j)
              zdf2  = max (zdfu**2 + zdfv**2, .001)

c             internal deformation:
              zs2a  = 2. * (0.5*(s2a0(ia,j) + s2a0(ib,j)))
              zbra   = zs2a * (zdf2**((powi-1.)/2.))
              zfluxg = zbra * zdfu *  hu(iau,j)**2

c             plus basal sliding:
              zcrh   = crhu2(iau,j)
              zbra   = zcrh * (zdf2**((powbu(iau,j)-1.)/2.))
              zfluxg = zfluxg + zbra * zdfu *  hu(iau,j)

c                Impose Schoof velocity, either at g.l. u-pt
c                (if schoof flx > grid flux, g.l. retreat),
c                or at downstream u-pt (if " < ", g.l. advance).
c                Also set muind = -1, to flag schoof-imposed points
c                (used in icestep). Ok to do as we go here, since muind
c                not used above.

              if ( abs(idir*hg*ug) .gt. abs(zfluxg)) then
                iz  = ia
                izu = iau
              else
                iz  = ib
                izu = ibu
              endif
c             account for *upstream* shelf/schoof advection in icestep.
c             limit |ua| for stability, in case h(iz,j) is v. small.
              ua(izu,j) = idir * min (ubmax, (hg*ug)/max(h(iz,j),.01))
              muind(izu,j) = -1

c                Set diagnostic g.l. quantites, returned:

              hgu(iau,j) = hg
              fluxschu(iau,j) = idir*hg*ug
              fluxgrdu(iau,j) = -zfluxg
c~~~~~~~~~~~~~~~~
            endif
c~~~~~~~~~~~~~~~~

          endif
c>>>>>>>>>>>>
        enddo
c>>>>>>>>>>>>
c==========
      enddo
c==========

c        v grid (same comments as above for u grid):

c==============
      do i=1,nx
c==============

        call izero (ifglv, nyp+1)
        do j=2,ny-2
          jp3 = min (j+3,ny)
          jm2 = max (j-2,1)
          if (maskv(i,j).eq.1) then
            if      ( h(i  ,j).gt.0. .and. maskwater(i,j  ).eq.0
     *                               .and. maskwater(i,j+1).ne.0
     *                               .and. maskwater(i,j+2).ne.0
     *                               .and. maskwater(i,jp3).ne.0 ) then
              ifglv(j) = 1
            else if ( h(i,j+1).gt.0. .and. maskwater(i,j+1).eq.0
     *                               .and. maskwater(i,j  ).ne.0
     *                               .and. maskwater(i,j-1).ne.0
     *                               .and. maskwater(i,jm2).ne.0 ) then
              ifglv(j) = -1
            endif
          endif
        enddo

c<<<<<<<<<<<<<<<<<<
        do j=2,ny-2
c<<<<<<<<<<<<<<<<<<
          if (iabs(ifglv(j)).eq.1) then

            if (ifglv(j).eq.1) then
              ja  = j
              jb  = j+1
              jc  = j+2
              jau = j
              jbu = j+1
              idir = 1
            else
              ja  = j+1
              jb  = j
              jc  = j-1
              jau = j
              jbu = j-1
              idir = -1
            endif
            if (thetav(i,jb).ne.-1.) then
              ztheta = thetav(i,jb)
            else if (thetav(i,jc).ne.-1.) then
              ztheta = thetav(i,jc)
            else
              ztheta = -1.
            endif

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (hs(i,ja).gt.hs(i,jb) .and. ztheta.gt.0.) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              pschm = 1./powbv(i,jau)
              zasch = 0.5*(arhap(i,ja)**powi+arhap(i,jb)**powi)
              zbsch = (1./crhv2(i,jau))**(1./powbv(i,jau))
              zpowb = powbv(i,jau)
              zfrac = fracgv(i,jau)
              zdx = 0.5*(dy(i,ja) + dy(i,jb))

              call calcgl ( h(i,ja), hs(i,ja), hb(i,ja),
     *                      h(i,jb), hs(i,jb), hb(i,jb),
     *                      zasch, zbsch, zpowb, zfrac, zdx,
     *                      sealev, ztheta, hg, ug)

              zdfu  = 0.25
     *              * ( dfu(i,ja) +dfu(i,jb) +dfu(i-1,ja) +dfu(i-1,jb))
              zdfv  = dfv(i,jau)
              zdf2  = max (zdfu**2 + zdfv**2, .001)

c             internal deformation:
              zs2a   = 2. * (0.5*(s2a0(i,ja) + s2a0(i,jb)))
              zbra   = zs2a * (zdf2**((powi-1.)/2.))
              zfluxg = zbra * zdfv *  hv(i,jau)**2

c             plus basal sliding:
              zcrh  = crhv2(i,jau)
              zbra  = zcrh * (zdf2**((powbv(i,jau)-1.)/2.))
              zfluxg = zfluxg + zbra * zdfv *  hv(i,jau)

              if ( abs(idir*hg*ug) .gt. abs(zfluxg)) then
                jz  = ja
                jzu = jau
              else
                jz  = jb
                jzu = jbu
              endif
              va(i,jzu) = idir * min (ubmax, (hg*ug)/max(h(i,jz),.01))
              mvind(i,jzu) = -1

              hgv(i,jau) = hg
              fluxschv(i,jau) = idir*hg*ug
              fluxgrdv(i,jau) = -zfluxg
c~~~~~~~~~~~~~~~~
            endif
c~~~~~~~~~~~~~~~~

          endif
c<<<<<<<<<<<<
        enddo
c<<<<<<<<<<<<
c==========
      enddo
c==========

      return
      end

c-----------------------------------------------------------------------

      subroutine calcgl ( h1, hs1, hb1,
     *                    h2, hs2, hb2,
     *                    acoef, bcoef, zpowb, frac, dx,
     *                    sealev, theta, hg, ug)

c     Calculates hg (g.l. depth) and ug (perturb stretching g.l. veloc).
c     h1,hs1,hb1 are upstream, and h2,hs2,hb2 are downstream (all supp),
c     regardless of actual grid direction.

      use comicephys
c     statement function for grid resol/non-resol term, to
c     eliminate explicit stretching perturbation if grid resolves
c     g.l. (see notes 9/27/07).
c     fresol(x) = 1.-exp(-(x**2))
      fresol(x) = 1. ! 777

c        Find upstream  "perturbation" stretching velocity at g.l. (ug)

c     calc g.l. depth (linear interp to u-grid of bed elev, and s.l.)
c     hbg =       0.5*hb1 +  0.5*hb2 ! 777
      hbg = (1.-frac)*hb1 + frac*hb2 ! 777  (a bit better 100 vs 300 ka)
      hg = (sealev-hbg)/rhor
      hsg = hbg + hg

c     1/e-folding length of upstream perturbation (stretching) veloc
      zlen = ( acoef * (0.5*bcoef)**zpowb  * hg**(powi - 2.*zpowb)
     *         * (0.25*rhoip*grav*theta)**(powi-zpowb)
     *       ) ** (1./(zpowb+1.))

c     perturbation (stretching) veloc
      ug = (0.25*rhoip*grav*hg*theta)**(zpowb*(powi+1.)/(zpowb+1.))
     *     * (2.*hg*acoef/bcoef)**(zpowb/(zpowb+1.))
     *     / (0.5*(1.-rhor))**(zpowb/(zpowb+1.)) ! ~Sch eq.16

      ug = ug * fresol(zlen*dx)

      return
      end

c>>>>>
c>>>>>

c-----------------------------------------------------------------------

      subroutine thetacalc (h, maskwater, arhap, muind, mvind, ua, va,
     *                      thetau, thetav)

c     Sets theta[u,v] = 1 - fraction of g.l. longitudinal stress
c     buttressed. These are used in calcgl if schooofgl is defined -
c     if not, just diagnostic.  Nb: theta[u,v] are on h-grid, for
c     E-W and N-S directions respectively, and ice thickness is h
c     (not hgu, hgv as in calcgl), so that if all shelf is free
c     floating, theta=1 exactly (property of elliptic solution).
c
c     Only set if one of the velocity pts on either side of the h-pt
c     was included in previous elliptical soln (m[u,v]ind). Else = -1.

      use comicephys
      use comicegrid

      dimension
     *  h(nx,ny),              maskwater(nx,ny),    arhap(nx,ny),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  thetau(nx,ny),         thetav(nx,ny)

c     parameter (hthetmin = 0.1)
c     parameter (hthetmin = 1.0)
      parameter (hthetmin = 10.)

      do j=1,ny
        do i=1,nx
          if ( maskwater(i,j).ne.0. .and.
     *         (muind(i-1,j).ge.1 .or. muind(i,j).ge.1) .and.
     *         h(i,j).gt.hthetmin ) then
            zdudx = max (0., (ua(i,j)-ua(i-1,j)) / dx(i,j))
            thetau(i,j) = 4. * ( zdudx**(1./powi) / arhap(i,j) )
     *                    / (rhoip*grav*h(i,j))
            thetau(i,j) = min (thetau(i,j), 1.)
          else
            thetau(i,j) = -1.
          endif
        enddo
      enddo

      do j=1,ny
        do i=1,nx
          if ( maskwater(i,j).ne.0. .and.
     *         (mvind(i,j-1).ge.1 .or. mvind(i,j).ge.1)  .and.
     *         h(i,j).gt.hthetmin ) then
            zdvdy = max (0., (va(i,j)-va(i,j-1)) / dy(i,j))
            thetav(i,j) = 4. * ( zdvdy**(1./powi) / arhap(i,j) )
     *                    / (rhoip*grav*h(i,j))
            thetav(i,j) = min (thetav(i,j), 1.)
          else
            thetav(i,j) = -1.
          endif
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine icestep (hold, h, hnew, hw,
     *                    hu,  hv,  dhsu,  dhsv,  dfu,  dfv,
     *                    s2a0, crhu2, crhv2, powbu, powbv,
     *                    ui, vi, ub, vb, ua, va,
     *                    masku, maskv, muind, mvind, dt)

c        Step ice thickness through one timestep, from hold
c        (beginning of timestep) to hnew (end of timestep, returned).
c        Use fully time-implicit Newton-Raphson ADI/tridia for
c        sheet flow. Use explicit ub from icedyn for basal/shelf
c        and deforming sed, but with implicit h in h*ub.

c        Sheet-flow ADI method operates on each horizontal direction
c        in turn, calling subr ice to step the 1-D ice flow equation.
c        Units are meters and years. The 1-D strip arrays *(n1d) are
c        used both for both directions in the domain.

c        For ub advection (sliding/sed def/shelf),
c        use an upstream scheme for h interp (or could use cubic SLT?),
c        since any centered-interp scheme gives spurious buildup
c        at the last non-zero grid box of a flat uniform ice shelf.

      use comicephys
      use comicegrid

c        hold  = initial thickness (m) (h-grid, supplied)
c        h     = mid-timestep ice-sheet thickness (m) (supplied)
c        hnew  = new ice-sheet thickness (m) (h-grid, returned)
c        hw    = liquid water thickness (m) (h-grid, supplied
c        ui = depth-averaged internal-shear u veloc (u-grid, returned)
c        vi = depth-averaged internal-shear v veloc (v-grid, returned)
c        dt = ice-model timestep (yr) (supplied)

      dimension
     *  hold(nx,ny),           h(nx,ny),          hnew(nx,ny),
     *  hw(nx,ny),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  dhsu(0:nxp,0:nyp),     dhsv(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),      dfv(0:nxp,0:nyp),

     *  s2a0(nx,ny),
     *  crhu2(0:nxp,0:nyp),    crhv2(0:nxp,0:nyp),
     *  powbu(0:nxp,0:nyp),    powbv(0:nxp,0:nyp),
     *  ui(0:nxp,0:nyp),       vi(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),       vb(0:nxp,0:nyp),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),    maskv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp)

c     local:
      dimension
     *  pconv(nx,ny)     ! tendency from first ADI stage (m/yr) (h-grid)

c        1-D strip arrays:

       dimension
     *  hold1(nx+ny),   hnew1(nx+ny),   darea1(nx+ny),
     *  pconv1(nx+ny),
     *  e(0:nx+ny),     f(0:nx+ny),     g(0:nx+ny),      flux(0:nx+ny),
     *  ubsp(0:nx+ny),  ubsm(0:nx+ny),  fluxub(0:nx+ny)

      call scopy (nx*ny, hold, 1, hnew, 1)

      call zero (ui,  (nxp+1)*(nyp+1))
      call zero (vi,  (nxp+1)*(nyp+1))
      call zero (ub,   (nxp+1)*(nyp+1))
      call zero (vb,   (nxp+1)*(nyp+1))

      if (powi.eq.3) then
        zpowi = 1.
        zpowim1 = 0.
      else
        zpowi = (powi-1)/2.
        zpowim1 = zpowi - 1.
      endif

c        First ADI pass, doing flow for all E-W strips.
c        Calculate tendency (pconv) due to E-W flow,
c        don't change 2-D ice thickness yet.

      call zero (pconv, nx*ny)

c~~~~~~~~~~~~~~~~~~~~~~
      if (nx.gt.1) then
c~~~~~~~~~~~~~~~~~~~~~~

c>>>>>>>>>>>>>>
      do j=1,ny
c>>>>>>>>>>>>>>

c        Find low-index edge of ice sheet (n1 is adjacent box with h=0).
c        If non-zero ice exists in first box, set n1=1, and b.c. will
c        implicitly be no flow out of left-hand edge of box.

        n1 = -1
        do i = 1,nx
          if (hold(i,j) + dt*max(pconv(i,j),0.).gt.0.) then
            n1 = max (i-1, 1)
            goto 10
          endif
        enddo
   10   if (n1.eq.-1) n1 = 1

c        Find high-index edge of ice sheet(n2 is adjacent box with h=0).
c        If non-zero ice exists in last box, set n2=nx, and b.c. will
c        implicitly be no flow out of right-hand edge of box.

        n2 = -1
        do i = nx,1,-1
          if (hold(i,j) + dt*max(pconv(i,j),0.).gt.0.) then
            n2 = min (i+1, nx)
            goto 20
          endif
        enddo
   20   if (n2.eq.-1) n2 = 1

c*************************
        if (n2.gt.n1) then
c*************************

c          Set 1-D strips (E-W)

        do i=n1,n2
          hold1(i) = hold(i,j)
          pconv1(i)= pconv(i,j)
          darea1(i)= darea(i,j)
        enddo

        do i=n1-1,n2
          ig = max(i,1)
          ip = min(i+1,nx)

c----------------------------------
c         not shelf(floating) nor schoof pts:
          if (muind(i,j).eq.0 .and. masku(i,j).eq.1) then
c----------------------------------

c           use quantities at beginning of timestep (*0, hold)
c           for time-implicit tridiag calcs:
            zdfv   = 0.25
     *           * ( dfv(i,j) +dfv(i+1,j) +dfv(i,j-1) +dfv(i+1,j-1))
            zdfu   = dfu(i,j)
            zdf2   = max (zdfu**2 + zdfv**2, .001)
            zdh    = (hold(ip,j) - hold(ig,j)) / dxu(i,j)

c              Internal deformation:

            zs2a   = 2. * (0.5*(s2a0(ig,j) + s2a0(ip,j)))
            zbra   = zs2a * (zdf2**zpowi)
            zflux  = zbra * zdfu *  hu(i,j)**2
            zdfluxdh  = 2. * zbra * zdfu *  hu(i,j)
            zdfluxda  = ( zbra + zs2a*2.*zpowi*(zdf2**zpowim1)*(zdfu**2)
     *                  ) * (hu(i,j)**2)

            e(i) = ( - zdfluxdh*hu(i,j)
     *               - zdfluxda*rhoice*grav*hu(i,j)*(dhsu(i,j) + zdh)
     *             ) * dyu(i,j)
            f(i) = (zdfluxdh + zdfluxda*rhoice*grav*dhsu(i,j))*dyu(i,j)
            g(i) = ( zdfluxda*rhoice*grav*hu(i,j)/dxu(i,j) ) * dyu(i,j)
            flux(i)  = -zflux * dyu(i,j)
            ui(i,j) = -zflux / max (hu(i,j), 1.e-20) ! diagnostic

c              Plus basal sliding:

            if (powbu(i,j).eq.3.) then
              zpowb = 1.
              zpowbm1 = 0.
            else
              zpowb = (powbu(i,j)-1.)/2.
              zpowbm1 = zpowb - 1.
            endif

            zcrh  = crhu2(i,j)
            zbra  = zcrh * (zdf2**zpowb)
            zflux = zbra * zdfu *  hu(i,j)
            zdfluxdh = zbra * zdfu
            zdfluxda = ( zbra + zcrh*2.*zpowb*(zdf2**zpowbm1)*(zdfu**2)
     *                 ) * hu(i,j)

            e(i) = e(i)
     *           + ( - zdfluxdh*hu(i,j)
     *               - zdfluxda*rhoice*grav*hu(i,j)*(dhsu(i,j) + zdh)
     *             ) * dyu(i,j)
            f(i) = f(i)
     *         + ( zdfluxdh + zdfluxda*rhoice*grav*dhsu(i,j) )*dyu(i,j)
            g(i) = g(i)
     *         + ( zdfluxda*rhoice*grav*hu(i,j) / dxu(i,j) ) * dyu(i,j)
            flux(i) = flux(i) - zflux * dyu(i,j)
            ub(i,j) = -zflux / max (hu(i,j), 1.e-20) ! diagnostic
            ua(i,j) = ui(i,j) + ub(i,j)               ! diagnostic

            fluxub(i) = 0.
            ubsp(i) = 0.
            ubsm(i) = 0.

c-------------
          else
c-------------

c              Contrib of ua, from shelf/schoof, upstrm, implicit wrt h

            e(i) = 0.
            f(i) = 0.
            g(i) = 0.
            flux(i) = 0.
            ub(i,j) = ua(i,j)                         ! diagnostic
            ui(i,j) = 0.                              ! diagnostic

            fluxub(i) = (  max(ua(i,j),0.)*hold(ig,j)
     *                   + min(ua(i,j),0.)*hold(ip,j)) * dyu(i,j)
            ubsp(i) = max(ua(i,j),0.)*dyu(i,j)
            ubsm(i) = min(ua(i,j),0.)*dyu(i,j)

c--------------
          endif
c--------------

        enddo

        e(n1-1) = 0.
        f(n1-1) = 0.
        g(n1-1) = 0.
        flux(n1-1) = 0.
        ubsp  (n1-1) = 0.
        ubsm  (n1-1) = 0.
        fluxub(n1-1) = 0.
        e(n2) = 0.
        f(n2) = 0.
        g(n2) = 0.
        flux(n2) = 0.
        ubsp  (n2) = 0.
        ubsm  (n2) = 0.
        fluxub(n2) = 0.

c          Call iceflow to do E-W flow

        call iceflow (hold1, hnew1,
     *                e, f, g, flux, ubsp, ubsm, fluxub, pconv1,
     *                darea1, n1, n2, nx, dt, j, 0)

c          Store tendency (dh/dt) due to E-W flow only.
c          But if ny=1 (no N-S flow), won't do N-S loop below,
c          so store final result in hnew.

        if (ny.gt.1) then
          do i=n1,n2
            pconv(i,j) = (hnew1(i)-hold1(i)) / dt
          enddo
        else
          do i=n1,n2
            hnew(i,j) = hnew1(i)
          enddo
        endif

c************
        endif   ! n2 > n1
c************
c>>>>>>>>>>
      enddo     ! j loop
c>>>>>>>>>>
c~~~~~~~~~~
      endif     ! nx > 1
c~~~~~~~~~~

c        Second ADI pass, doing flow for all N-S strips.
c        Use tendency (pconv) saved from first ADI pass.
c        Set new ice thickness (hnew).

c~~~~~~~~~~~~~~~~~~~~~~
      if (ny.gt.1) then
c~~~~~~~~~~~~~~~~~~~~~~

c<<<<<<<<<<<<<<
      do i=1,nx
c<<<<<<<<<<<<<<

c        Find low-index edge of ice sheet (n1 is adjacent box with h=0).
c        If non-zero ice exists in first box, set n1=1, and b.c. will
c        implicitly be no flow out of bottom edge of box.

        n1 = -1
        do j = 1,ny
          if (hold(i,j) + dt*max(pconv(i,j),0.).gt.0.) then
            n1 = max (j-1, 1)
            goto 30
          endif
        enddo
   30   if (n1.eq.-1) n1 = 1

c        Find high-index edge of ice sheet(n2 is adjacent box with h=0).
c        If non-zero ice exists in last box, set n2=ny, and b.c. will
c        implicitly be no flow out of top edge of box.

        n2 = -1
        do j = ny,1,-1
          if (hold(i,j) + dt*max(pconv(i,j),0.).gt.0.) then
            n2 = min (j+1, ny)
            goto 40
          endif
        enddo
   40   if (n2.eq.-1) n2 = 1

c*************************
        if (n2.gt.n1) then
c*************************

c          Set 1-D strips (N-S)

        do j=n1,n2
          hold1(j) = hold(i,j)
          pconv1(j)= pconv(i,j)
          darea1(j)= darea(i,j)
        enddo

        do j=n1-1,n2
          jg = max(j,1)
          jp = min(j+1,ny)

c----------------------------------
c         not shelf(floating) nor schoof pts:
          if (mvind(i,j).eq.0 .and. maskv(i,j).eq.1) then
c----------------------------------

c            use quantities at beginning of timestep (*0, hold)
c           for time-implicit tridiag calcs:
            zdfu   = 0.25
     *           * ( dfu(i-1,j) +dfu(i,j) +dfu(i-1,j+1) +dfu(i,j+1))
            zdfv   = dfv(i,j)
            zdf2   = max (zdfu**2 + zdfv**2, .001)
            zdh    = (hold(i,jp) - hold(i,jg)) / dyv(i,j)

c              Internal deformation:

            zs2a  = 2. * (0.5*(s2a0(i,jg) + s2a0(i,jp)))
            zbra  = zs2a * (zdf2**zpowi)
            zflux = zbra * zdfv *  hv(i,j)**2
            zdfluxdh = 2. * zbra * zdfv *  hv(i,j)
            zdfluxda = ( zbra + zs2a*2.*zpowi*(zdf2**zpowim1)*(zdfv**2)
     *                 ) * (hv(i,j)**2)

            e(j) = ( - zdfluxdh*hv(i,j)
     *               - zdfluxda*rhoice*grav*hv(i,j)*(dhsv(i,j) + zdh)
     *             ) * dxv(i,j)
            f(j) = (zdfluxdh + zdfluxda*rhoice*grav*dhsv(i,j))*dxv(i,j)
            g(j) = ( zdfluxda*rhoice*grav*hv(i,j)/dyv(i,j) ) * dxv(i,j)
            flux(j) = -zflux  * dxv(i,j)
            vi(i,j) = -zflux / max (hv(i,j), 1.e-20) ! diagnostic

c              Plus basal sliding:

            if (powbv(i,j).eq.3.) then
              zpowb = 1.
              zpowbm1 = 0.
            else
              zpowb = (powbv(i,j)-1.)/2.
              zpowbm1 = zpowb - 1.
            endif

            zcrh  = crhv2(i,j)
            zbra  = zcrh * (zdf2**zpowb)
            zflux = zbra * zdfv *  hv(i,j)
            zdfluxdh = zbra * zdfv
            zdfluxda = ( zbra + zcrh*2.*zpowb*(zdf2**zpowbm1)*(zdfv**2)
     *                 ) * hv(i,j)

            e(j) = e(j)
     *           + ( - zdfluxdh*hv(i,j)
     *               - zdfluxda*rhoice*grav*hv(i,j)*(dhsv(i,j) + zdh)
     *             ) * dxv(i,j)
            f(j) = f(j)
     *         + ( zdfluxdh + zdfluxda*rhoice*grav*dhsv(i,j) )*dxv(i,j)
            g(j) = g(j)
     *         + ( zdfluxda*rhoice*grav*hv(i,j) / dyv(i,j) ) * dxv(i,j)
            flux(j) = flux(j) - zflux * dxv(i,j)
            vb(i,j) = -zflux / max (hv(i,j), 1.e-20) ! diagnostic
            va(i,j) = vi(i,j) + vb(i,j)               ! diagnostic

            fluxub(j) = 0.
            ubsp(j) = 0.
            ubsm(j) = 0.

c-------------
          else
c-------------

c              Contrib of va (from shelf/schoof, upstrm, implicit wrt h)

            e(j) = 0.
            f(j) = 0.
            g(j) = 0.
            flux(j) = 0.
            vb(i,j) = va(i,j)                         ! diagnostic
            vi(i,j) = 0.                              ! diagnostic

            fluxub(j) = (  max(va(i,j),0.)*hold(i,jg)
     *                   + min(va(i,j),0.)*hold(i,jp)) * dxv(i,j)
            ubsp(j) = max(va(i,j),0.)*dxv(i,j)
            ubsm(j) = min(va(i,j),0.)*dxv(i,j)
c--------------
          endif
c--------------

        enddo

        e(n1-1) = 0.
        f(n1-1) = 0.
        g(n1-1) = 0.
        flux(n1-1) = 0.
        ubsp  (n1-1) = 0.
        ubsm  (n1-1) = 0.
        fluxub(n1-1) = 0.
        e(n2) = 0.
        f(n2) = 0.
        g(n2) = 0.
        flux(n2) = 0.
        ubsp  (n2) = 0.
        ubsm  (n2) = 0.
        fluxub(n2) = 0.

c          Call iceflow to do N-S flow

        call iceflow (hold1, hnew1,
     *                e, f, g, flux, ubsp, ubsm, fluxub, pconv1,
     *                darea1, n1, n2, ny, dt, j, 1)

c          Store new 2-D ice thickness

        do j=n1,n2
          hnew(i,j) = hnew1(j)
        enddo

c************
        endif   ! n2 > n1
c************
c<<<<<<<<<<
      enddo     ! i loop
c<<<<<<<<<<
c~~~~~~~~~~
      endif     ! ny > 1
c~~~~~~~~~~

      return
      end

c-----------------------------------------------------------------------

      subroutine iceflow (h0, h1,
     *                    e, f, g, flux, ubsp, ubsm, fluxub, yconv,
     *                    dar, n1, n2, nxory, dt, j, ifv)

c        Steps ice-sheet 1-D flow equation through one time step.
c        Time-implicit contributions lead to a tridiagonal matrix
c        equation. Units are meters and years.

c        The current direction is called the "x-direction" (index i)
c        here, and the perpendicular direction is the "y-direction",
c        accounted for either by yconv.
c        Only the x-direction has time-implicit terms.
c        The correspondence with real E-W or N-S directions depends
c        on the current ADI stage in the calling program.

c        Ice-sheet thickness are defined on h-grid (h).
c        Other terms (...u) are defined on a staggered u-grid, with
c        ...u(i) at the interface between h-grid boxes i and i+1.

c        h0   = actual ice thickness at start of time step (m)(supplied)
c        h1   = actual ice thickness at end of time step (m) (returned)
c        yconv= convergence rate due to perpendicular flow (m/y) (supp)
c        dar  = grid-box area (h-grid) (m**2) (supplied)
c        dt   = timestep (yr) (supplied)

      use comicephys
      use comicegrid

      dimension
     *  h0(nx+ny), h1(nx+ny),
     *  e(0:nx+ny), f(0:nx+ny), g(0:nx+ny), flux(0:nx+ny),
     *  ubsp(0:nx+ny), ubsm(0:nx+ny), fluxub(0:nx+ny), yconv(nx+ny),
     *  dar(nx+ny)

      dimension
     *  a(nx+ny), b(nx+ny), c(nx+ny), forc(nx+ny)

      parameter (rimp = 1.)         ! time-implicit vs explicit fraction
c     parameter (rimp = 0.5)        ! time-implicit vs explicit fraction

c        Set implicit vs explicit timesteps

      dti =     rimp *dt
      dte = (1.-rimp)*dt

c       Set the 3 diagonals of the tridiagonal matrix, and the rhs

      do i = n1,n2
        im = max(i-1,1)
        ip = min(i+1,nx)

        a(i) = 1.+ dti*(   (-0.5*f(i)   + g(i))
     *                   + ( 0.5*f(i-1) + g(i-1))
     *                   + ( ubsp (i)   - ubsm (i-1))    ! upstream
c    *                   + ( 0.5*(-ubsp(i-1) + ubsp(i))) ! centered
     *                 ) / dar(i)

        b(i) =     dti*(    -0.5*f(i)   - g(i)
     *                              + ubsm (i)           ! upstream
c    *                          + 0.5*ubsp (i)           ! centered
     *                 ) / dar(i)
        c(i) =     dti*(     0.5*f(i-1) - g(i-1)
     *                              - ubsp (i-1)         ! upstream
c    *                          - 0.5*ubsp (i-1)         ! centered
     *                 ) / dar(i)

        forc(i)= h0(i)
     *           + dti*(  e(i)      - e(i-1)      ) / dar(i)
     *           + dt *( -flux(i)   + flux(i-1)   ) / dar(i)
     *           + dte*( -fluxub(i) + fluxub(i-1) ) / dar(i)
     *           + dt * yconv(i)
      enddo

c        Solve tridiagonal system (steps ice thickness from h0 to h1)

      call tridia_i (a(n1), b(n1), c(n1), h1(n1), forc(n1), n2-n1+1, 0)

      return
      end

c-----------------------------------------------------------------------

      subroutine tridia_i (a, b, c, h, forc, n, ktrid)

c        Solves tridiagonal matrix equation.
c        a,b,c are diagonals of the matrix.
c        a is main diagonal (1:n)
c        b is superdiagonal (1:n-1),
c        c is subdiagonal (2:n),
c        forc is right-hand side,
c        h is vector to be solved for.
c        n is the size of the matrix and vectors.
c        If ktrid is non-zero, uses condition that resulting h(i)
c        must be >=  0.

      dimension a(n), b(n), c(n), h(n), forc(n)
      parameter (nmax = 25000)
      dimension r(nmax), s(nmax), y(nmax)

c        Check that local arrays are big enough

      if (n.gt.nmax) then
        write(6,*)'*** Error in tridia_i: Too many points.',
     *            '   nmax=',nmax,'   n=',n
        stop
      endif

c     amin = 1.e-3
      amin = 1.e-6
      do i = 1,n
        if (abs(a(i)).lt.amin) then
          write (6,9000) i,a(i)
 9000     format(' tridia_i warning: a(',i3,')=',f12.8)
          if (a(i).ge.0.) a(i) = amin
          if (a(i).lt.0.) a(i) = -amin
        endif
      enddo

      s(1) = a(1)
      r(1) = b(1)/s(1)
      if (n.gt.2) then
        do i = 2,n-1
          s(i) = a(i)-c(i)*r(i-1)
          r(i) = b(i)/s(i)
        enddo
      endif
      if (n.ge.2) then
        s(n) = a(n)-c(n)*r(n-1)
      endif

      y(1) = forc(1)/s(1)
      if (n.ge.2) then
        do i = 2,n
          y(i) = (forc(i)-c(i)*y(i-1)) / s(i)
        enddo
      endif

      if (ktrid.eq.0) then
        h(n) = y(n)
        if (n.ge.2) then
          do i = n-1,1,-1
            h(i) = y(i)-r(i)*h(i+1)
          enddo
        endif
      else
        h(n) = max (0., y(n))
        if (n.ge.2) then
          do i = n-1,1,-1
            h(i) = max (0., y(i)-r(i)*h(i+1))
          enddo
        endif
      endif

      return
      end
