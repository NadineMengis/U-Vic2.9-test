      subroutine iceshow1d (h, hs, hb, t,
     *                      basefrml, oceanmelt, budgall,
     *                      tsurf, tsurfi, heati, heatb,
     *                      w, wa, maskh,
     *                      sedim, sedimun, tsed, wsed, heats, 
     *                      quarryrate,
     *                      topbed, topbedeq, tbed,
     *                      hw, tw, maskwater,
     *                      u, v, ua, va, ui, vi, ub, vb,
     *                      uadv, vadv, hu, hv, dfu, dfv,
     *                      muind, mvind, thetau, thetav,
     *                      crhu, crhv, fsedu, fsedv, 
     *                      numh, itera, iterc,
     *                      toth0, tota0, toth, tota, totflow, totneg,
     *                      sealev, dtantann, dtantjan, dtseas, rco2,
     *                      timeice, dt, weirun,
     *                      iloop, nloopend, nyearout1d)

c     Ascii 1-D printout - flowlines, within 1-D or 2-D domains.
c     Flowline for 2-D is defined by i1s,j1s,n1s (in comicegrid),
c     set in first call below.

      use comicephys
      use comicegrid

      dimension 
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny), 
     *  t(nx,ny,0:nlevp),
     *  basefrml(nx,ny),
     *  oceanmelt(nx,ny),    budgall(nx,ny),
     *  tsurf(nx,ny),        tsurfi(nx,ny),
     *  heati(nx,ny,nlev),   heatb(nx,ny), 
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp),
     *  maskh(nx,ny),
     *  sedim(nx,ny),        sedimun(nx,ny),
     *  tsed(nx,ny,nsed),    wsed(nx,ny,nsed),   heats(nx,ny,nsed),   
     *  quarryrate(nx,ny),
     *  topbed(nx,ny),       topbedeq(nx,ny),    tbed(nx,ny,nbed),
     *  hw(nx,ny),           tw(nx,ny),          maskwater(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  ui(0:nxp,0:nyp),        vi(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),      vadv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),        hv(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),       dfv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),     mvind(0:nxp,0:nyp), 
     *  thetau(nx,ny),          thetav(nx,ny),
     *  crhu(0:nxp,0:nyp),      crhv(0:nxp,0:nyp),
     *  fsedu(0:nxp,0:nyp),     fsedv(0:nxp,0:nyp)

      logical first1s
      data first1s /.true./
      save first1s

      dimension 
     *  work1d(nx+ny),
     *  work1dv(0:nx+ny+1),   iwork1dv(0:nx+ny+1), 
     *  work2d(0:nxp,0:nyp),  work2da(0:nxp,0:nyp)
      character*5 cwork1d(nx+ny), cwork1da(nx+ny)

#if defined EISA
      dimension hanal(nx)
#endif

#if defined EISVAR20 || defined EISVAR40
      logical firstvar
      data firstvar /.true./  
      save firstvar, hvar0, tvar0
#endif


c       Skip out if not time

      if (nyearout1d.eq.0) return
      if (.not.(mod(abs(timeice)+0.5*dt,max(float(nyearout1d),dt)).lt.dt
#if defined HOMA || defined HOMB || defined HOMC || defined HOMD || defined HOME
     *            .or. iloop.eq.nloopend
#else
     *            .or. iloop.eq.1 .or. iloop.eq.nloopend
#endif
     *          )
     *   ) return 


c        Overall mass budget 

      totbud  = 0.
      totbudi = 0.
      do j=1,ny
        do i=1,nx
          totbud = totbud     + budgall(i,j) * darea(i,j) * dt
          if (maskh(i,j).eq.1)
     *      totbudi = totbudi + budgall(i,j) * darea(i,j) * dt
        enddo
      enddo
      toterr = (toth-toth0) - (totbud + totneg + totflow) 

c        Print general 0-D quantities  

      print*,trim(isname),' Diagnostics:'
      do iprint=1,2
        if (iprint.eq.1) iu = iuout1d
        if (iprint.eq.2) iu = 6
!         write (iu,'(/80("=")/
!      *               a,i10, a,f7.4, a,2i8, a,i5, 3(a,f8.2), a,f8.3)') 
!      *    'time(yr)=',     nint(timeice),
!      *    '   weirun=',    weirun,
!      *    '   iter[a-c]=', itera,iterc,
!      *    '   sealev=',    nint(sealev),
!      *    '   dtantann=',  dtantann,
!      *    '   dtantjan=',  dtantjan,
!      *    '   dtseas  =',  dtseas,
!      *    '   rco2    =',  rco2
        zz = dt*totarea
        write (iu,'(a,5f14.8,f21.15)')
     *   'd(h), bud, budi, flow, neg, err =',
     *    (toth-toth0)/zz, totbud/zz, totbudi/zz, totflow/zz, totneg/zz,
     *    toterr/zz
        call flush (iu)
      enddo

! c        Print 1-D or flowline output 
! 
! c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! #if defined EISANTA || defined NHA || defined CARB
!       do iprint=1,1
! #else
!       do iprint=1,2
! #endif
!         if (iprint.eq.1) iu = iuout1d
!         if (iprint.eq.2) iu = 6
! c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! c=======================================================================
! #if defined EIS1 || defined EIS2 || defined EIS3 || defined EIS4
! c=======================================================================
! 
!         iinc = max (nx/100, 1)
! c       iinc = 1
! 
!         write (iu,'(a/8x,1001i6)') 'ua:', 
!      *                          (nint(xh(i,(ny+1)/2)*1.e-3),i=1,nx,iinc)
!         do j=ny,1,-iinc
!           jkm = nint (yh((nx+1)/2,j)*1.e-3)
!           write (iu,'(i5,3x,1001i6)') jkm, (nint(ua(i,j)),i=1,nx,iinc)
!         enddo
! 
!         write (iu,'(a/8x,1001i6)') 'va:', 
!      *                          (nint(xh(i,(ny+1)/2)*1.e-3),i=1,nx,iinc)
!         do j=ny,1,-iinc
!           jkm = nint (yh((nx+1)/2,j)*1.e-3)
!           write (iu,'(i5,3x,1001i6)') jkm, (nint(va(i,j)),i=1,nx,iinc)
!         enddo
! 
!         write (iu,'(a/8x,1001i6)') 'h:', 
!      *                          (nint(xh(i,(ny+1)/2)*1.e-3),i=1,nx,iinc)
!         do j=ny,1,-iinc
!           jkm = nint (yh((nx+1)/2,j)*1.e-3)
!           write (iu,'(i5,3x,1001i6)') jkm, (nint(h(i,j)),i=1,nx,iinc)
!         enddo
! 
! c=======================================================================
! #elif defined EISA || defined EISB || defined EISC || defined EISD
! c=======================================================================
! 
!         if (ny.eq.1) then
! 
!           ja = (ny+1)/2
!           ia = (nx+1)/2
!           jkm = nint (yh(ia,ja))*1.e-3
! #if defined EISA
! c         analytic 1-D profile:
!           zm = budgall(ia,ja)
!           za = 2.e-16
!           zk = 2.*za*((rhoice*grav)**powi) / (zm*(powi+2.))
!           zkp = zk**(1./powi)
!           zpowa = (powi+1.)/powi
!           zpowb = 2.*(powi+1.)/powi
!           zl = xh(nx,ja)
!           do i=1,nx
!             zx = abs(xh(i,ja))
!             zh = ( (2./zkp) * (zl**zpowa - zx**zpowa) ) ** (1./zpowb)
!             hanal(i) = zh
!           enddo
!           write (iu,'(a/8x,40i6)') 'h,hanal:', 
!      *                              (nint(xh(i,ja)*1.e-3),i=1,nx) 
!           write (iu,'(i5,3x,40i6)') jkm, (nint(h(i,ja)) ,i=1,nx)
!           write (iu,'(i5,3x,40i6)') jkm, (nint(hanal(i)),i=1,nx)
! #else
!           write (iu,'(a/8x,40i6)') 'h:', 
!      *                              (nint(xh(i,ja)*1.e-3),i=1,nx) 
!           write (iu,'(i5,3x,40i6)') jkm, (nint(h(i,ja))  ,i=1,nx)
! c    *                              (nint(yh(ia,j)*1.e-3),j=1,ny) 
! c         write (iu,'(i5,3x,40i6)') jkm, (nint(h(ia,j))  ,j=1,ny)
! #endif
! 
!         else
! 
! c         2-D maps:
!           jc = (ny+1)/2
!           write (iu,'(a/8x,1001i6)') 'h (m):', 
!      *                              (nint(xh(i,jc)*1.e-3),i=1,nx)
!           do j=ny,1,-1
!             jkm = nint (yh((nx+1)/2,j)*1.e-3)
!             write (iu,'(i5,3x,1001i6)') jkm, (nint(h(i,j)),i=1,nx)
!           enddo
! 
!           write (iu,'(a/8x,1001i6)') '10*tb:', 
!      *                              (nint(xh(i,jc)*1.e-3),i=1,nx)
!           do j=ny,1,-1
!             jkm = nint (yh((nx+1)/2,j)*1.e-3)
!             write (iu,'(i5,3x,1001i6)') jkm,
!      *            (nint(10.*(t(i,j,nlevp)-(tmelt-dtmdh*h(i,j)))),i=1,nx)
!           enddo
! 
! c         horizontal transects:
!           j = (ny+1)/2
!           write(iu,*)
! c         write(iu,'(a,2x,1001i6)')   'km   :',
! c    *                                      (nint(xh(i,j)*1.e-3),i=1,nx)
!           write(iu,'(a,2x,1001i6)')   'h    :',
!      *                                      (nint(h(i,j))       ,i=1,nx)
!           write(iu,'(a,2x,1001f6.1)') 'tsurf:',
!      *                                      (tsurf(i,j)-tmelt,   i=1,nx)
!           write(iu,'(a,2x,1001f6.1)') 'tb   :',
!      *                        (t(i,j,nlevp)-(tmelt-dtmdh*h(i,j)),i=1,nx)
! 
!           if (ny.eq.1) then 
!             do i=1,nx
!               work1d(i)=tsurfi(i,j)+.042*h(i,j)/2.1-(tmelt-dtmdh*h(i,j))
!             enddo
!             write(iu,'(a,2x,1001f6.1)') 'tbanl:', (work1d(i),i=1,nx)
!           endif
! 
!           do i=1,nx
!             work1d(i) = 0.
!             do k=1,nlev
!               work1d(i) = work1d(i) + heati(i,j,k)*dzeta(k)*h(i,j)
!             enddo
!           enddo
!           write(iu,'(a,2x,1001i6)') 'heati:',
!      *                        (nint(1.e3*work1d(i)/(31556926)),  i=1,nx)
!           write(iu,'(a,2x,1001i6)') 'heatb:',
!      *                       (nint(1.e3*heatb(i,j)/(31556926)),  i=1,nx)
!           write(iu,'(a,5x,1001i6)')   'uadv :',
!      *                         (nint(uadv(i,j)*1.e-3/dyu(i,j)),i=1,nx-1)
!           write(iu,'(a,5x,1001i6)')   'ua   :',
!      *                         (nint(ua(i,j))                 ,i=1,nx-1)
!           write(iu,'(a,5x,1001i6)')   'ub   :',
!      *                         (nint(ub(i,j))                 ,i=1,nx-1)
! 
! c         vertical profiles:
!           i1 = (nx+1)/2
! c         i2 = (nx+1)/2 + (nx+1)/4
! c         i2 = (nx+1)/2 + (nx+1)/4 + 1
!           i2 = (nx+1)/2 + nint(0.2*(nx+1))
!           j = (ny+1)/2
!           write (iu,'(/2(8x,7x,a,i5,a)/2a/(f8.3,3f8.2,8x,3f8.2))')
!      *                     'divide (', nint(xh(i1,j)*1.e-3), ' km)',
!      *                     'mid pt (', nint(xh(i2,j)*1.e-3), ' km)',
!      *                     '  1-zeta       t       u      wa',
!      *                     '               t       u      wa',
!      *                     ( 1.-zeta(k),
!      *                       ( t(i,j,k)-(tmelt-dtmdh*h(i,j)*zeta(k)), 
!      *                         0.5*(u(i,j,k)+u(i-1,j,k)), 
!      *                         wa(i,j,k), i=i1,i2,i2-i1 ),
!      *                       k=1,nlev )
! 
! c         domain-wide quantities:
!           totamelt = 0.
!           do j=1,ny
!             do i=1,nx
!               if ( h(i,j).gt.0. .and.
!      *             t(i,j,nlevp)-(tmelt-dtmdh*h(i,j)) .gt.-0.01 )
!      *          totamelt = totamelt + darea(i,j)
!             enddo
!           enddo
!           write (iu,'(/(a,f8.3))') 
!      *      'vol  (10^6 km3)=', toth*1.e-15,
!      *      'area (10^6 km2)=', tota*1.e-12, 
!      *      'melt fraction  =', totamelt / max(tota,1.e-3)
! 
! #if defined EISVAR20 || defined EISVAR40
! c         time series of particular quantities: 
!           if (iprint.eq.1)  then
!             if (firstvar) then
!               hvar0 = h(i1,j)
!               tvar0 = t(i1,j,nlevp)-(tmelt-dtmdh*h(i1,j))
!               write (iueistab,'(a)') '    time    dh    dt  uadv  marg'
!               firstvar = .false.
!             endif
!             imarg = 0
!             do i=nx,(nx+1)/2,-1
!               if (h(i,j).gt.0. .and. imarg.eq.0) imarg = i
!             enddo
!             write (iueistab,'(i8,i6,f6.1,i6,i6)')
!      *        nint(timeice),
!      *        nint(h(i1,j) - hvar0),
!      *        t(i1,j,nlevp)-(tmelt-dtmdh*h(i1,j)) - tvar0,
!      *        nint(uadv(i2,j)*1.e-3/dyu(i2,j)),
!      *        imarg
!             call flush (iueistab)
!           endif
! #endif
! 
!         endif
! 
! c=============================================================
! #elif defined EISLINE || defined HOME || defined EISANTA
! c=============================================================
! 
! c          Set 1-D flowline indices within full 1-D or 2-D domain.
! c          i1s,j1s,n1s are in common (comicegrid)
! 
!         if (first1s) then
! #if defined EISANTA
! c         for flowline from Weddell to Ross Seas, ~through Siple Dome
! c                            alona   alata   alonb   alatb
! c         call setflowline ( -110.,  -83.0,   170.,   -70.)
!           call setflowline (  -40.,   -70.,   178.,   -70.)
! #else
!           n1s = 0 
!           if (nx.ge.ny) then
!             do i = 1, nx, max(nx/100,1) 
!               n1s = n1s + 1 
!               i1s(n1s) = i
!               j1s(n1s) = (ny+1)/2
!               d1s(n1s) = xh(i1s(n1s),j1s(n1s))
!             enddo
!           else
!             do j = 1, ny, max(ny/100,1) 
!               n1s = n1s + 1 
!               i1s(m) = (nx+1)/2
!               j1s(m) = j
!               d1s(n1s) = yh(i1s(n1s),j1s(n1s))
!             enddo
!           endif
! #endif
!           first1s = .false.
!         endif
! 
!         write(iu,'(a,1001i5)') 'i1s :',
!      *                (i1s(m),m=1,n1s)
!         write(iu,'(a,1001i5)') 'j1s :',
!      *                (j1s(m),m=1,n1s)
!         write(iu,'(a,1001i5)') 'km  :',
!      *                (nint(d1s(m)*1.e-3),m=1,n1s)
!         write(iu,*)
!         write(iu,'(a,1001i5)') 'h   :',
!      *                (nint(h(i1s(m),j1s(m))),m=1,n1s)
!         write(iu,'(a,1001i5)') 'hs  :',
!      *                (nint(hs(i1s(m),j1s(m))),m=1,n1s)
! c       write(iu,'(a,1001i5)') 'hbi :',
! c    *                (nint(hs(i1s(m),j1s(m))-h(i1s(m),j1s(m))),m=1,n1s)
!         write(iu,'(a,1001i5)') 'hb  :',
!      *                (nint(hb(i1s(m),j1s(m))),m=1,n1s)
!  
!         do m=1,n1s
!           if (sedim(i1s(m),j1s(m)).gt.0.949) then
!             write (cwork1d (m),'(i5)') nint(sedim(i1s(m),j1s(m)))
!             write (cwork1da(m),'(i5)') nint(sedimun(i1s(m),j1s(m)))
!           else
!             write (cwork1d (m),'(f5.1)') sedim(i1s(m),j1s(m))
!             write (cwork1da(m),'(f5.1)') sedimun(i1s(m),j1s(m))
!           endif
!         enddo
!         write(iu,'(a,1001a5)') 'sed :', (cwork1d(m),m=1,n1s)
!         write(iu,'(a,1001a5)') 'sedu:', (cwork1da(m),m=1,n1s)
! 
!         write(iu,'(a,1001i5)') 'topb:',
!      *                    (nint(topbed(i1s(m),j1s(m))),m=1,n1s)
!         write(iu,'(a,1001i5)') 'topq:',
!      *                    (nint(topbedeq(i1s(m),j1s(m))),m=1,n1s)
!         write(iu,'(a,1001i5)') 'hw  :',
!      *                    (nint(hw(i1s(m),j1s(m))),m=1,n1s)
! #if defined LINEE
!        write(iu,'(a,1001i5)') 'hwX :',
!      *                    (nint(100.*hw(i1s(m),j1s(m))),m=1,n1s)
! #endif
!         write(iu,'(a,1001i5)') 'tw  :',
!      *                  (nint(10.*(tw(i1s(m),j1s(m))-tmelt)),m=1,n1s)
!         write(iu,'(a,1001i5)') 'mask:',
!      *                    (maskwater(i1s(m),j1s(m)),m=1,n1s)
! 
!         call iflowline (muind, mvind, iwork1dv, 0)
!         write(iu,'(a,2x,1001i5)') 'muin:', (iwork1dv(m),m=1,n1s-1)
! 
!         write(iu,'(a,1001i5)') 'zflo:',
!      *    ( nint (max ( 0., h(i1s(m),j1s(m)) 
!      *                      - (sealev-hb(i1s(m),j1s(m)))/rhor )),
!      *      m=1,n1s )
! #if defined LINEE
!         write(iu,'(a,1001i5)') 'zflX:',
!      *    ( nint ( 100. * max (0., h(i1s(m),j1s(m))
!      *                            - (sealev-hb(i1s(m),j1s(m)))/rhor)),
!      *      m=1,n1s )
! #endif
! 
! c       like flowline, but for h-grid thetau,v:
!         do m=1,n1s 
!           if (j1s(m).eq.j1s(min(n1s,m+1))) then
!             work1d(m) = thetau(i1s(m),j1s(m))
!           else
!             work1d(m) = thetav(i1s(m),j1s(m))
!           endif
!         enddo
!         write(iu,'(a,1001i5)') 'thet:', 
!      *                                   (nint(work1d(m)*1.e2),m=1,n1s)
! 
!         write(iu,'(a,1001i5)') 'budg:',
!      *                   (nint(1000.*budgall(i1s(m),j1s(m))),m=1,n1s)
!         write(iu,'(a,1001i5)') 'basf:',
!      *                  (nint(1000.*basefrml(i1s(m),j1s(m))),m=1,n1s)
!         write(iu,'(a,1001i5)') 'ocnm:',
!      *                 (nint(1000.*oceanmelt(i1s(m),j1s(m))),m=1,n1s)
!         write(iu,'(a,1001i5)') 'quar:',
!      *                 (nint(1.e6*quarryrate(i1s(m),j1s(m))),m=1,n1s)
! 
!         do m=1,n1s
!           work1d(m) = 0.
!           do k=1,nlev
!             work1d(m) = work1d(m)
!      *            + heati(i1s(m),j1s(m),k)*dzeta(k)*h(i1s(m),j1s(m))
!           enddo
!         enddo
!         write(iu,'(a,1001i5)') 'heai:',
!      *                  (nint(1000.*work1d(m)/31556926), m=1,n1s)
!         write(iu,'(a,1001i5)') 'heab:',
!      *            (nint(1000.*heatb(i1s(m),j1s(m))/31556926),m=1,n1s)
!         do m=1,n1s
!           work1d(m) = 0.
!           do k=1,nsed
!             work1d(m) = work1d(m) + heats(i1s(m),j1s(m),k)
!           enddo
!         enddo
!         write(iu,'(a,1001i5)') 'heas:',
!      *                  (nint(1000.*work1d(m)/31556926), m=1,n1s)
! 
!         call flowline (ua, va, work1dv, 1)
!         write(iu,'(a,2x,1001i5)') 'ua  :',(nint(work1dv(m)), m=1,n1s-1)
! 
!         call flowline (ui, vi, work1dv, 1)
!         write(iu,'(a,2x,1001i5)') 'ui  :',(nint(work1dv(m)), m=1,n1s-1)
! 
!         call flowline (ub, vb, work1dv, 1)
!         write(iu,'(a,2x,1001i5)') 'ub  :',(nint(work1dv(m)), m=1,n1s-1)
! 
!         do j=0,ny
!           do i=0,nx
!             work2d(i,j)  = uadv(i,j)/dyu(i,j)
!             work2da(i,j) = vadv(i,j)/dxu(i,j)
!           enddo
!         enddo
!         call flowline (work2d, work2da, work1dv, 1)
!         write(iu,'(a,2x,1001i5)') 'uadv:',
!      *                              (nint(work1dv(m)*1.e-2), m=1,n1s-1)
! 
!         call flowline (dfu, dfv, work1dv, 1)
!         write(iu,'(a,2x,1001i5)') 'dfu :',
!      *                               (nint(work1dv(m)/100.), m=1,n1s-1)
! 
!         do j=0,ny
!           do i=0,nx
!             if (crhu(i,j).eq.0.) then
!               work2d(i,j)  = 0.
!             else
!               work2d (i,j) = log10(crhu(i,j))
!             endif
!             if (crhv(i,j).eq.0.) then
!               work2da(i,j) = 0.
!             else
!               work2da(i,j) = log10(crhv(i,j))
!             endif
!           enddo
!         enddo
!         call flowline (work2d, work2da, work1dv, 1)
!         write(iu,'(a,2x,1001i5)') 'cbsu:',
!      *                             (nint(10.*work1dv(m)),m=1,n1s-1)
! 
!         call flowline (fsedu, fsedv, work1dv, 1)
!         write(iu,'(a,2x,1001i5)') 'fseu:',
!      *                                (nint(100.*work1dv(m)),m=1,n1s-1)
! 
! 
!         write (iu,'(/a)') 't*10:'
!         do k=0,nlevp
!           write (iu,'(f5.3,1001i5)') 
!      *     1.-zeta(k), (nint(10.*(t(i1s(m),j1s(m),k)-tmelt)),m=1,n1s)
!         enddo
! 
!         write (iu,'(a,1001i5)') 'tbme:', 
!      *       (nint(10.*(tmelt-dtmdh*h(i1s(m),j1s(m))-tmelt)),m=1,n1s)
! 
!         write (iu,'(/a)') 'tsed:'
!         do k=1,nsed
!           write (iu,'(f5.3,1001i5)') 1.-zsed(k), 
!      *              (nint(10.*(tsed(i1s(m),j1s(m),k)-tmelt)),m=1,n1s)
!         enddo
! 
!         write (iu,'(/a)') 'wsed:'
!         do k=1,nsed
!           write (iu,'(f5.3,1001i5)') 
!      *       1.-zsed(k), (nint(100.*(wsed(i1s(m),j1s(m),k))),m=1,n1s)
!         enddo
! 
!         write (iu,'(/a,1001i5)') 'tbed:',
!      *              (nint(10.*(tbed(i1s(m),j1s(m),1)-tmelt)),m=1,n1s)
! 
!         write (iu,'(/a)') 'u:'
!         do k=1,nlev
!           do j=0,ny
!             do i=0,nx
!               work2d(i,j)  = u(i,j,k)
!               work2da(i,j) = v(i,j,k)
!             enddo
!           enddo
!           call flowline (work2d, work2da, work1dv, 1)
!           write (iu,'(f5.3,2x,1001i5)') 1.-zeta(k), 
!      *                                     (nint(work1dv(m)),m=1,n1s-1)
!         enddo
! 
! c=====
! #endif
! c=====
! 
!         call flush (iu)
! c~~~~~~~~~~
!       enddo
! c~~~~~~~~~~
! 
! #if defined LINEE && defined INITEE  
!       if (iloop.eq.nloopend) then
!         iu = 26  
!         write (iu,
!      *    "(6x,'dimension hwedge(nx+1),hswedge(nx+1),hw_wedge(nx+1)')" )
!         write (iu,"( '      data hwedge /')")
!         write (iu,"(('     * ',5(f11.5,',')))") h 
!         write (iu,"( '     *   0. /')")
! 
!         write (iu,"('      data hswedge /')")
!         write (iu,"(('     * ',5(f11.5,',')))") hs
!         write (iu,"( '     *   0. /')")
! 
!         write (iu,"('      data hw_wedge/')")
!         write (iu,"(('     * ',5(f11.5,',')))") hw
!         write (iu,"( '     *   0. /')")
!       endif
! #endif

      return
      end

c-----------------------------------------------------------------------

#if defined EISANTA

      subroutine setflowline (alona, alata, alonb, alatb)

c     For 2D Antarctica runs, sets coords for 1D-flowline in iceshow1d: 
c     n1s, i1s, j1s, d1s (in comicegrid, returned). 
c     a[lon,lat][a,b] are lons,lats of 2 end points (supplied).

      use comicephys
      use comicegrid

c        Convert lon,lat of 2 end pts to i,j polar stereo coords.

      zcos = cos(alata*pi/180.) 
      zr =   ( (1./zcos) - sqrt ((1./zcos)**2 - 1.) ) * zlambda
      zxa = radius * zr * cos (0.5*pi - alona*pi/180.) 
      zya = radius * zr * sin (0.5*pi - alona*pi/180.) 
      ia = (zxa + 0.5*nx*dx0)/dx0 + 1.0001
      ja = (zya + 0.5*ny*dy0)/dy0 + 1.0001

      zcos = cos(alatb*pi/180.) 
      zr =   ( (1./zcos) - sqrt ((1./zcos)**2 - 1.) ) * zlambda
      zxb = radius * zr * cos (0.5*pi - alonb*pi/180.) 
      zyb = radius * zr * sin (0.5*pi - alonb*pi/180.) 
      ib = (zxb + 0.5*nx*dx0)/dx0 + 1.0001
      jb = (zyb + 0.5*ny*dy0)/dy0 + 1.0001

      ia = max (1, min (nx, ia))
      ja = max (1, min (ny, ja))
      ib = max (1, min (nx, ib))
      jb = max (1, min (ny, jb))

c        Set transect coords (i1s,j1s). If transect is more in x 
c        direction, increment i uniformly by 1 and calculate nearest j,
c        or v.v if transect is more in y direction.

      if (iabs(ia-ib).ge.iabs(ja-jb)) then
        n1s = iabs(ib-ia) + 1
        m = 0
        do i = ia,ib,isign(1,ib-ia)
          m = m + 1
          i1s(m) = i
          j1s(m) = ja + nint((m-1.)*(jb-ja)/(n1s-1.))
        enddo
      else
        n1s = iabs(jb-ja) + 1
        m = 0
        do j = ja,jb,isign(1,jb-ja)
          m = m + 1
          j1s(m) = j
          i1s(m) = ia + nint((m-1.)*(ib-ia)/(n1s-1.))
        enddo
      endif

c       Set distance (m) along transect 

      d1s(1) = 0.5*dd0
      do m = 2,n1s
        d1s(m) = d1s(m-1)
     *         + sqrt (  (xh(i1s(m),j1s(m)) - xh(i1s(m-1),j1s(m-1)))**2
     *                 + (yh(i1s(m),j1s(m)) - yh(i1s(m-1),j1s(m-1)))**2)
      enddo

c        Diagnostic output

      iu = 89
      if (iu.ne.0) then
        write (iu,'(a,i6,/2(a,2f8.2/),a)') 
     *    'setflowline:  n1s=', n1s,
     *    'alona, alata = ', alona, alata,
     *    'alonb, alatb = ', alonb, alatb,
     *    '   m   alond   alatd       d1s'

        do m=1,n1s
          write (iu,'(i4,2f8.2, f10.1)')
     *      m, alond(i1s(m),j1s(m)), alatd(i1s(m),j1s(m)), d1s(m)*1.e-3
        enddo
        call flush (iu)
      endif

      return
      end

#endif

c-----------------------------------------------------------------------

      subroutine flowline (uarr,   varr,  work, ifreal)

      use comicegrid

      entry     iflowline (iuarr, ivarr, iwork, ifreal)

c        Assembles 1-D flowline for a vector quantity (u or v-grid),
c        based on h-grid indices i1s,j1s,n1s (in comicegrid).
c        ifreal = 1 for real, else integer.
c        Distinguish between u or v direction by j1s changing or not.  

      dimension 
     *  uarr (0:nxp,0:nyp), varr (0:nxp,0:nyp),  
     *  iuarr(0:nxp,0:nyp), ivarr(0:nxp,0:nyp),  
     *  work(0:nx+ny+1),    iwork(0:nx+ny+1)  

c        First (zero index) point possibly on domain boundary 
c        (erroneous if i1s or j1s decrease away from nx,ny boundary)

      if (ifreal.eq.1) then
c either:
c       if (j1s(1).eq.j1s(2)) then                         ! u direction
c         work(0) = uarr(i1s(1)-1,j1s(1))
c       else
c         work(0) = varr(i1s(1),j1s(1)-1)                  ! v direction
c       endif
c or:
        work(0) = sqrt (   uarr(i1s(1)-1,j1s(1))**2 
     *                   + varr(i1s(1),j1s(1)-1)**2 )
      else
        if (j1s(1).eq.j1s(2)) then 
          iwork(0) = iuarr(i1s(1)-1,j1s(1))                ! u direction
        else
          iwork(0) = ivarr(i1s(1),j1s(1)-1)                ! v direction
        endif
      endif

c       Rest of (interior) points

      if (ifreal.eq.1) then
        do m=1,n1s-1
c either:
c         if (j1s(m).eq.j1s(m+1)) then                     ! u direction
c           work(m)  =  uarr(min(i1s(m),i1s(m+1)), j1s(m)) 
c         else                                             ! v direciton
c           work(m)  =  varr(i1s(m), min(j1s(m),j1s(m+1)))
c         endif
c or:
          work(m) = sqrt (   uarr(min(i1s(m),i1s(m+1)), j1s(m))**2
     *                     + varr(i1s(m), min(j1s(m),j1s(m+1)))**2 )
        enddo
      else
        do m=1,n1s-1
          if (j1s(m).eq.j1s(m+1)) then 
            iwork(m) = iuarr(min(i1s(m),i1s(m+1)), j1s(m)) ! u direction
          else
            iwork(m) = ivarr(i1s(m), min(j1s(m),j1s(m+1))) ! v direction
          endif
        enddo
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine iceshow2d (h, hs, hb, t,  
     *                    budgall, budgsnow, budgrain, budgmelt,
     *                    baseperc, basefrml, oceanmelt, calvrate,
     *                    tsurf, tsurfi, heati, heatb,
     *                    w, wa, sedim,  quarryrate,
     *                    topbed, topbedeq,
     *                    hw, tw, maskwater,
     *                    u, v, ua, va, ub, vb,
     *                    uadv, vadv, hu, hv, dfu, dfv,
     *                    thetau, thetav, hgu, hgv,
     *                    fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *                    crhu, crhv,
     *                    numh, itera, iterc,
     *                    toth0, tota0, toth, tota, totflow, totneg,
     *                    sealev, timeice, dt, weirun,
     *                    iloop, nloopend, nyearout2d)

c     Ascii 2-D printout - x,y maps

      use comicephys
      use comicegrid

      dimension 
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny), 
     *  t(nx,ny,0:nlevp),    
     *  budgall(nx,ny),      budgsnow(nx,ny),    budgrain(nx,ny),
     *  budgmelt(nx,ny),
     *  baseperc(nx,ny),     basefrml(nx,ny),     
     *  oceanmelt(nx,ny),    calvrate(nx,ny),     
     *  tsurf(nx,ny),        tsurfi(nx,ny), 
     *  heati(nx,ny,nlev),   heatb(nx,ny), 
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp),
     *  sedim(nx,ny),        quarryrate(nx,ny),
     *  topbed(nx,ny),       topbedeq(nx,ny), 
     *  hw(nx,ny),           tw(nx,ny),          maskwater(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),      vadv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),        hv(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),       dfv(0:nxp,0:nyp),
     *  thetau(nx,ny),          thetav(nx,ny),
     *  hgu(0:nxp,0:nyp),       hgv(0:nxp,0:nyp),
     *  fluxgrdu(0:nxp,0:nyp),  fluxgrdv(0:nxp,0:nyp),
     *  fluxschu(0:nxp,0:nyp),  fluxschv(0:nxp,0:nyp),
     *  crhu(0:nxp,0:nyp),      crhv(0:nxp,0:nyp)

      dimension work(nx,ny)
      logical firstshow2
      data firstshow2 /.true./
      save firstshow2

c       Skip out if not time

      if (nyearout2d.eq.0) return
      if (.not.(mod(abs(timeice)+0.5*dt,max(float(nyearout2d),dt)).lt.dt
     *            .or. iloop.eq.1 .or. iloop.eq.nloopend
     *          )
     *   ) return 

c       Set maskdisplay (in common), used in printmap

      call scopy_i (nx*ny, maskwater, 1, maskdisplay, 1)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do iprint=1,1
c     do iprint=1,2
       if (iprint.eq.1) iu = iuout2d
       if (iprint.eq.2) iu = 6
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c          Write 2-D ascii maps

        if (firstshow2) then
          call printmap (timeice,alatd,   ' lat (deg S)',    5., h,iu,0)
          call printmap (timeice,alond,   ' lon (deg E)',   10., h,iu,0)
          call printmap (timeice,topbedeq,'equil bdrk (m)',150., h,iu,0)

          call zero (work, nx*ny)
          do m=1,n1s
            work(i1s(m),j1s(m)) = 1.
          enddo
          call printmap (timeice,work,'flowline',1., h,iu,0)

          if (iprint.eq.2) firstshow2 = .false.
        endif
  
        call itor (maskwater, work, nx*ny)
        call printmap (timeice,work,    'maskwater', 1.,  h, iu,0)
        call printmap (timeice,hw,      'hw(m)         ',40.,  h, iu,0)

        do j=1,ny
          do i=1,nx
            work(i,j) = tw(i,j) - tmelt
            if (maskwater(i,j).eq.0) work(i,j) = 0.
          enddo
        enddo
        call printmap (timeice,work,     'tw(C)        ', 0.1,  h, iu,0)

        call printmap (timeice,baseperc, 'baseperc(m/y)', 0.1, h, iu,0)
        call printmap (timeice,basefrml, 'basefrml(m/y)', .01, h, iu,0)
        call printmap (timeice,oceanmelt,'oceanmelt(m/y)',0.5, h, iu,0)
        call printmap (timeice,calvrate, 'calvrate(m/y)', 0.1, h, iu,0)
        call printmap (timeice,budgall,  'budgall (m/y)', 0.1, h, iu,0)
        call printmap (timeice,budgsnow, 'budgsnow(m/y)', 0.1, h, iu,0)
c       call printmap (timeice,budgsnow, 'budgsnow(m/y)', .01, h, iu,0)
        call printmap (timeice,budgrain, 'budgrain(m/y)', 0.1, h, iu,0)
        call printmap (timeice,budgmelt, 'budgmelt(m/y)', 0.1, h, iu,0)
        call printmap (timeice,hb,'bed+sed elev (m)',    150., h, iu,0)
  
        do j=1,ny
          do i=1,nx
            work(i,j) = tsurf(i,j) - tmelt
          enddo
        enddo
        call printmap (timeice,work,'tsurf (C)', 2., h, iu,0)

        do j=1,ny
          do i=1,nx
            work(i,j) = t(i,j,0) - tmelt
          enddo
        enddo
        call printmap (timeice,work,'ice surf t (C)', 2., h, iu,0)
  
        do j=1,ny
          do i=1,nx
            work(i,j) = t(i,j,nlevp) - (tmelt-dtmdh*h(i,j))
c           work(i,j) = t(i,j,nlevp) -  tmelt
            if (h(i,j).eq.0.) work(i,j) = 0. 
          enddo
        enddo
        call printmap (timeice,work, 'basal t (C)',  1., h, iu,0)
        call printmap (timeice,work, 'basal t (C)',  .1, h, iu,0)
#if defined HEINO
        write (iu,*)
        do j=ny,1,-1
          write (iu,'(81i5)') (nint(10.*work(i,j)),i=1,nx)
        enddo
#endif

        do j=1,ny
          do i=1,nx
            work(i,j) = log10(crhu(i,j))
          enddo
        enddo
        call printmap (timeice, work, 'log10(crhu)', 1., h, iu,1)

        do j=1,ny
          do i=1,nx
            work(i,j) = log10(crhv(i,j))
          enddo
        enddo
        call printmap (timeice, work, 'log10(crhv)', 1., h, iu,2)
  
        do j=1,ny
          do i=1,nx
            work(i,j) = ua(i,j)
          enddo
        enddo
        call printmap (timeice,work, 'ua',  10.0, h, iu,1)
        call printmap (timeice,work, 'ua',  100.0,h, iu,1)
  
        do j=1,ny
          do i=1,nx
            work(i,j) = va(i,j)
          enddo
        enddo
        call printmap (timeice,work, 'va',  10.0, h, iu,2)
        call printmap (timeice,work, 'va',  100.0,h, iu,2)

c---------------------
#if defined SCHOOFGL
c---------------------
        do j=1,ny
          do i=1,nx
            work(i,j) = fluxgrdu(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'fluxgrdu', 1.e3, h, iu,-1)

        do j=1,ny
          do i=1,nx
            work(i,j) = fluxgrdv(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'fluxgrdv', 1.e3, h, iu,-2)

        do j=1,ny
          do i=1,nx
            work(i,j) = fluxschu(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'fluxschu', 1.e3, h, iu,-1)

        do j=1,ny
          do i=1,nx
            work(i,j) = fluxschv(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'fluxschv', 1.e3, h, iu,-2)

        do j=1,ny
          do i=1,nx
            if (abs(fluxschu(i,j)).ne.0.) then
              work(i,j) = min (1.e3, max (-1.e3, 
     *                                    fluxgrdu(i,j)/fluxschu(i,j) ))
            else
              work(i,j) = 0.
            endif
          enddo
        enddo
        call printmap (timeice, work, 'fluxgrdu/schu', 1., h, iu,-1)

        do j=1,ny
          do i=1,nx
            if (abs(fluxschv(i,j)).ne.0.) then
              work(i,j) = min (1.e3, max (-1.e3, 
     *                                    fluxgrdv(i,j)/fluxschv(i,j) ))
            else
              work(i,j) = 0.
            endif
          enddo
        enddo
        call printmap (timeice, work, 'fluxgrdv/schv', .1, h, iu,-2)
c-----
#endif
c-----

        do j=1,ny
          do i=1,nx
            work(i,j) = hgu(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'hgu (m)', 40., h, iu,1)

        do j=1,ny
          do i=1,nx
            work(i,j) = hgv(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'hgv (m)', 40., h, iu,2)

        do j=1,ny
          do i=1,nx
            work(i,j) = 1. - thetau(i,j)
          enddo
        enddo
        call printmap (timeice, work, '1-thetau', .05, h, iu,0)

        do j=1,ny
          do i=1,nx
            work(i,j) = 1. - thetav(i,j)
          enddo
        enddo
        call printmap (timeice, work, '1-thetav', .05, h, iu,0)

c#if defined HEINO
c       write (iu,*)
c       do j=ny-1,1,-1
c         write (iu,'(81i5)') (nint(10.*va(i,j)),i=1,nx)
c       enddo
c#endif
  
c-----------------------------------------------------------------------
#if defined HOMA || defined HOMB || defined HOMC || defined HOMD || defined HOMF
c-----------------------------------------------------------------------
        do j=1,ny
          do i=1,nx
            work(i,j) = u(i,j,0)
          enddo
        enddo
        call printmap (timeice,work, 'utop',  5., h, iu,1)
#if defined HOMF
       write (iu,*) 'u(i,j,0):'
       do j=ny-1,1,-1
         write (iu,'(101i5)') (nint(u(i,j,0)),i=1,nx)
       enddo
#endif
  
        do j=1,ny
          do i=1,nx
            work(i,j) = v(i,j,0)
          enddo
        enddo
        call printmap (timeice,work, 'vtop',  0.1, h, iu,2)
  
        do j=1,ny
          do i=1,nx
            work(i,j) = wa(i,j,0)
          enddo
        enddo
        call printmap (timeice,work, 'wtop',  0.5, h, iu,0)
  
        do j=1,ny
          do i=1,nx
            work(i,j) = -dfu(i,j)/1.e3
          enddo
        enddo
        call printmap (timeice,work, '-dfu (kPa)',  5., h, iu,1)
  
        do j=1,ny
          do i=1,nx
            work(i,j) = -dfv(i,j)/1.e3
          enddo
        enddo
        call printmap (timeice,work, '-dfv (kPa)',  0.1, h, iu,2)
c-----
#endif
c-----
  
        do j=1,ny
          do i=1,nx
            work(i,j) = ub(i,j)
          enddo
        enddo
        call printmap (timeice,work, 'ub',  10.0, h, iu,1)
  
        do j=1,ny
          do i=1,nx
            work(i,j) = vb(i,j)
          enddo
        enddo
        call printmap (timeice,work, 'vb',  10.0, h, iu,2)

c#if defined HEINO
c       write (iu,*)
c       do j=ny-1,1,-1
c         write (iu,'(81i5)') (nint(10.*vb(i,j)),i=1,nx)
c       enddo
c#endif
  
        do j=1,ny
          do i=1,nx
            work(i,j) = sqrt(  (0.5*(ub(i,j)+ub(i-1,j)))**2 
     *                       + (0.5*(vb(i,j)+vb(i,j-1)))**2 )
          enddo
        enddo
        call printmap (timeice,work, '|uvb|',  0.1, h, iu,0)
        call printmap (timeice,work, '|uvb|',  100., h, iu,0)
  
        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.0.) then
              zua = sqrt(  (0.5*(ua(i,j)+ua(i-1,j)))**2 
     *                   + (0.5*(va(i,j)+va(i,j-1)))**2 )
              zub = sqrt(  (0.5*(ub(i,j)+ub(i-1,j)))**2 
     *                   + (0.5*(vb(i,j)+vb(i,j-1)))**2 )
              work(i,j) = min (1.e6, zub / max (zua,1.e-20))
            else
              work(i,j) = 0.
            endif
          enddo
        enddo
        call printmap (timeice,work, '|uvb/uva|',  0.01, h, iu,0)
        call printmap (timeice,work, '|uvb/uva|',  0.1,  h, iu,0)
  
        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.0.) then
              work(i,j) = sqrt(  (0.5*(ua(i,j)+ua(i-1,j)))**2 
     *                         + (0.5*(va(i,j)+va(i,j-1)))**2 )
            else
              work(i,j) = 0.
            endif
          enddo
        enddo
        call printmap (timeice,work, '|uvabs|',  10., h, iu,0)
  
        do j=1,ny
          do i=1,nx
            zua = 0.5*(ua(i,j)+ua(i-1,j))
            zva = 0.5*(va(i,j)+va(i,j-1))
            work(i,j) = 1.e-6*h(i,j)*sqrt(zua**2+zva**2)
          enddo
        enddo
        call printmap (timeice,work, '|uvhabs|', 0.01, h, iu,0)
    
        call zero (work,nx*ny)
        do k=1,nlev
          do j=1,ny
            do i=1,nx
              work(i,j) = work(i,j) + heati(i,j,k)*dzeta(k)*h(i,j)
     *                                / (rhoice*hfus)
            enddo
          enddo
        enddo
        call printmap (timeice,work,'heati (m/y)',0.001,h,iu,0)
    
        call zero (work,nx*ny)
        do j=1,ny
         do i=1,nx
            work(i,j) = heatb(i,j) / (rhoice*hfus)
          enddo
        enddo
        call printmap (timeice,work,'heatb (m/y)',0.001,h,iu,0)

        do j=1,ny
          do i=1,nx
c           if (maskwater(i,j).ne.1) then
              work(i,j) = h(i,j) - (rholiq/rhoice)*(sealev-hb(i,j))
c           else
c             work(i,j) = 0.
c           endif
          enddo
        enddo
        call printmap (timeice,work,'hflot (m)',10.,h,iu,0)

        call printmap (timeice,hs,'hs (m)', 150., h, iu,0)
        call printmap (timeice,h,  'ice thickness (m)', 150., h, iu,0)
        call printmap (timeice,h,  'ice thickness_10 (m)', 10., h, iu,0)
#if defined HEINO
        write (iu,*)
        do j=ny,1,-1
          write (iu,'(81i5)') (nint(h(i,j)),i=1,nx)
        enddo
#endif

        call flush (iu)

c~~~~~~~~~~
      enddo
c~~~~~~~~~~
      return
      end

c-----------------------------------------------------------------------

      subroutine printmap (time, arrin, cdesc, scale, h, iu, igrid)

c        Displays array arr to unit iu in "a1" format.
c        igrid = 0 for h-grid variables, 1 for u-grid, 2 for v-grid.

      use comicegrid

      dimension arrin(nx,ny), h(nx,ny)
      character*(*) cdesc

      dimension arr(nx,ny)
      character*1 carr(nx), cbrr*2000
      character*80 chem
      logical first
      data first /.true./
      save iac, jac

      character crev*4, csul*4, coff*3
      parameter (crev = char(27)//'[7m',           ! reverse video
     *           csul = char(27)//'[4m',           ! underline
     *           coff = char(27)//'[m' )           ! reset all
c     character cbra*1
c     integer*1 icbra
c     equivalence (icbra,cbra)
c     data icbra /O'154000000000000000000'/        ! "ctrl["
c     save crev, csul, coff

      if (first) then
c       crev = cbra // '[7m'
c       csul = cbra // '[4m'
c       coff = cbra // '[m'
        iac = ichar('A') - 1
        jac = ichar('a') - 1
      endif

      call scopy (nx*ny, arrin, 1, arr, 1)

#if defined HOMA || defined HOMB || defined HOMC || defined HOMD
c     periodic boundary conditions over 2:nx-1:
c     do j=1,ny
c       arr(1,j) = arr(nx-1,j)
c       arr(nx,j) = arr(2,j)
c     enddo
#endif
#if defined HOMA || defined HOMC
c     periodic boundary conditions over 2:ny-1:
c     do i=1,nx
c       arr(i,1) = arr(i,ny-1)
c       arr(i,nx) = arr(i,2)
c     enddo
#endif

      if (nx.le.1000) then
        nhem = 1
      else
        nhem = 2
      endif

c-------------------
      do ihem=1,nhem
c-------------------

        if (nhem.eq.1) then
          chem = ' '
          ia = 1
          ib = nx
        else if (ihem.eq.1) then
          chem = 'western: '
          ia = nx/2 + 1
          ib = nx
        else
          chem = 'eastern: '
          ia = 1
          ib = nx/2
        endif

        if (.not. first) write (iu,*) 
        write(iu,100)
     *    nint(time), chem(1:lenchr(chem)), cdesc, scale
#if defined LONLAT
        if (iabs(igrid).eq.0 .or. iabs(igrid).eq.2) then
          write(iu,101) 
     *      (nint(  alond(i,(ny+1)/2)), i=5,ib,5)
        else
          write(iu,101) 
     *      (nint(0.5*(alond(i,(ny+1)/2)+alond(min(i+1,nx),(ny+1)/2))),
     *      i=5,ib,5)
        emdif
#else
        if (iabs(igrid).eq.0 .or. iabs(igrid).eq.2) then
          write(iu,101)
     *      (nint(.001*xh(i,(ny+1)/2)), i=5,ib,5)
        else
          write(iu,101) 
     *      (nint(.001*0.5* (xh(i,(ny+1)/2)+xh(min(i+1,nx),(ny+1/2)))),
     *       i=5,ib,5)
        endif
#endif
  100   format(
     *    'time:',i10,'  ',a,'  ',a,' / ',f11.4,'  (A-Z pos, a-z neg):')
  101   format(/6x,2x,200i5)

c>>>>>>>>>>>>>>>>>>>>>
        do j=ny,1,-1
c>>>>>>>>>>>>>>>>>>>>>

          do i=ia,ib
            if (abs(arr(i,j)).le.1.e-8) then
              carr(i) = '.'
            else
              n = nint(arr(i,j)/scale)
              if (n.eq.0) then
                if (arr(i,j).eq.0.) then
                  carr(i) = '.'
                else if (arr(i,j).gt.0.) then
                  carr(i) = '+'
                else
                  carr(i) = '-'
                endif
              else if (n.ge.1..and. n.le.26) then
                carr(i) = char(iac + n)
              else if (n.le.-1..and. n.ge.-26) then
                carr(i) = char(jac - n)
              else if (n.gt.26) then
                carr(i) = '*'
              else if (n.lt.-26) then
                carr(i) = '#'
              endif
            endif
          enddo

          if (.false.) then

c              Show ocean (not lakes) as blanks

            lenb = 0
            do i=ia,ib 
              lenb = lenb + 1
              if (maskdisplay(i,j).eq.1) then    ! maskdisplay in common
                cbrr(lenb:lenb) = ' '
              else
                cbrr(lenb:lenb) = carr(i)
              endif
            enddo
 
          else

c              For "2-D" calving, reverse-video for ocean/lake edge pts
c              (using h-grid maskdisplay, in common). igrid (passed) is
c              0 for h-grid arrays, 1 for u-grid, 2 for v-grid.
c              Also underline where ice (h > 0, ifi)...not done now.

            istatprev = 0
            lenb = 0
            do i = ia,ib
              ifi = 0
c             if (h(i,j).gt.0.) ifi = 1          ! doesn't work with lwy
              ifo = 0
              im = max (i-1,1)
              ip = min (i+1,nx)
              jm = max (j-1,1)
              jp = min (j+1,ny)
c             if (maskdisplay(i,j).ne.0) ifo = 1     ! for any water pts

c             if (igrid.eq.0) then
                if(
     *           ( maskdisplay(i,j).eq.0 .and. 
     *             (maskdisplay(im,j).ne.0.or.maskdisplay(ip,j).ne.0.or.
     *              maskdisplay(i,jm).ne.0.or.maskdisplay(i,jp).ne.0) 
     *           ) 
     *            ) ifo = 1
c             else if (igrid.eq.1) then
c               if ( (maskdisplay(i,j).eq.0.and.maskdisplay(ip,j).ne.0.)
c    *               .or.
c    *               (maskdisplay(i,j).ne.0.and.maskdisplay(ip,j).eq.0.)
c    *             ) ifo = 1
c             else if (igrid.eq.2) then
c               if ( (maskdisplay(i,j).eq.0.and.maskdisplay(i,jp).ne.0.)
c    *               .or.
c    *               (maskdisplay(i,j).ne.0.and.maskdisplay(i,jp).eq.0.)
c    *             ) ifo = 1
c             endif

              if (igrid.lt.0) ifo = 0  ! for g.l. related quantities ?

              istat = ifi + 2*ifo

              if (istat.ne.istatprev) then
                cbrr(lenb+1:lenb+3) = coff
                lenb = lenb + 3
                if (ifi.eq.1) then
                  cbrr(lenb+1:lenb+4) = csul
                  lenb = lenb + 4
                endif
                if (ifo.eq.1) then
                  cbrr(lenb+1:lenb+4) = crev
                  lenb = lenb + 4
                endif
              endif

              cbrr(lenb+1:lenb+1) = carr(i)
              lenb = lenb + 1

              if (i.eq.ib .and. istat.ne.0) then
                cbrr(lenb+1:lenb+3) = coff
                lenb = lenb + 3
              endif

              istatprev = istat
            enddo

          endif

#if defined LONLAT
          if (iabs(igrid).eq.0 .or. iabs(igrid).eq.1) then
            write(iu,'(f6.1,2x,a)') alatd((nx+1)/2,j), cbrr(1:lenb)
          else 
            write(iu,'(f6.1,2x,a)')
     *         0.5*( alatd((nx+1)/2,j) + alatd((nx+1)/2,min(j+1,ny)) ),
     *         cbrr(1:lenb)
          endif
#else
          if (iabs(igrid).eq.0 .or. iabs(igrid).eq.1) then
            write(iu,'(i6,2x,a)') nint(.001*yh((nx+1)/2,j)),cbrr(1:lenb)
          else
            write(iu,'(i6,2x,a)') 
     *       nint(.001*0.5*(yh((nx+1)/2,j) + yh((nx+1)/2,min(j+1,ny)))),
     *       cbrr(1:lenb)
          endif
#endif

c>>>>>>>>>>>>
        enddo
c>>>>>>>>>>>>

c----------
      enddo
c----------

      first = .false.

      return
      end

c-----------------------------------------------------------------------

      subroutine iceplot1d (h, hs, hb, t,  
     *                      budgall, basefrml, 
     *                      oceanmelt, calvrate,
     *                      tsurf, tsurfi, heati, heatb,
     *                      w, wa, maskh,
     *                      sedim, sedimun, tsed, wsed,heats,quarryrate,
     *                      topbed, topbedeq, tbed,
     *                      hw, tw,
     *                      u, v, ua, va, ui, vi, ub, vb, 
     *                      uadv, vadv, hu, hv, 
     *                      masku, maskv, crhu, crhv,
     *                      fracgu, fracgv, thetau, thetav,
     *                      sealev, dtantann, dtantjan, dtseas, rco2,
     *                      timeice, dt, weirun,
     *                      iloop, nloopend, nyearplot1d)
  
c     Writes ascii output file for 1-D flowline plotting, within 
c     1-D or 2-D domains. Flowline is defined by i1s,j1s,n1s
c     in comicegrid, set in iceshow1d.

      use comicephys
      use comicegrid

      dimension 
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny), 
     *  t(nx,ny,0:nlevp),    
     *  budgall(nx,ny),      basefrml(nx,ny),    
     *  oceanmelt(nx,ny),    calvrate(nx,ny),
     *  tsurf(nx,ny),        tsurfi(nx,ny),
     *  heati(nx,ny,nlev),   heatb(nx,ny), 
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp),
     *  maskh(nx,ny),
     *  sedim(nx,ny),        sedimun(nx,ny), 
     *  tsed(nx,ny,nsed),    wsed(nx,ny,nsed),
     *  heats(nx,ny,nsed),   quarryrate(nx,ny),
     *  topbed(nx,ny),       topbedeq(nx,ny),    tbed(nx,ny,nbed),
     *  hw(nx,ny),           tw(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  ui(0:nxp,0:nyp),        vi(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),      vadv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),        hv(0:nxp,0:nyp),
     *  crhu(0:nxp,0:nyp),      crhv(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),     maskv(0:nxp,0:nyp), 
     *  fracgu(0:nxp,0:nyp),    fracgv(0:nxp,0:nyp),
     *  thetau(nx,ny),          thetav(nx,ny)

      dimension work(0:nxp)
      character*80 cform, cformi

c     number of cross-sec vertical points, and bottom and top extents
#if defined EISLINE && defined LINEB
      parameter (nz=200, zbot = -1000., ztop=3000.)
#elif defined EISLINE && defined LINED
      parameter (nz=200, zbot = -2500., ztop=3000.) 
#elif defined EISLINE && defined LINEE
      parameter (nz=200, zbot = -600., ztop=100.) 
#elif defined EISLINE && defined LINEF
      parameter (nz=200, zbot = -1000., ztop=4000.)
#elif defined EISLINE && defined LINEG
      parameter (nz=200, zbot = -1000., ztop=4000.)
#elif defined EISLINE && defined TRANSECTA
      parameter (nz=200, zbot = -1000., ztop=4000.)
#else
c     parameter (nz=200, zbot = -1000., ztop=3000.) 
      parameter (nz=200, zbot = -2500., ztop=3000.) 
#endif

      dimension 
     *  work1d(nx+ny),
     *  work1dv(0:nx+ny+1),   iwork1dv(0:nx+ny+1),
     *  work2d(0:nxp,0:nyp),  work2da(0:nxp,0:nyp)

c        Skip out if not time

      if (nyearplot1d.eq.0) return
      if ( .not. 
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyearplot1d),dt)).lt.dt
     *     .or. iloop.eq.1 .or. iloop.eq.nloopend
     *     )
     *   ) return 

      iu = iuplot1d
c     cform  = '(a/1001f12.5)'
c     cformi = '(a/1001i12  )'
      cform  = '(a/(300f12.5))'
      cformi = '(a/(300i12  ))'

      if (iloop.eq.1) then
        nspw = max (1,nint(nyearplot1d/dt))
        nwrite = nloopend / nspw
        if (nspw.gt.1) nwrite = nwrite + 1
#if defined EISLINE && defined LINEE
        nwrite =  4000./nyearplot1d + 1    ! bb,cc
#endif
        write (iu,'(a,i8)') 'nwriteplot=',nwrite
      endif

#if defined EISLINE && defined LINEE
      if (nint(timeice).lt.6000)  return ! bb,cc
#endif

c        Write 0-D quantities

      igl = 0
      fgl = 0.
      xgl  = 0.
      do m=1,n1s-1
        if (igl.eq.0 .and. fracgu(i1s(m),j1s(m)).lt.1.) then 
          igl = i1s(m)
          fgl = fracgu(i1s(m),j1s(m))
          xgl = (1.-fgl)*xh(i1s(m),j1s(m))+fgl*xh(i1s(m+1),j1s(m+1))
        endif
      enddo

      write (iu,'(80("="))')
      write (iu,'( a,i9/a,f9.4/4(a,i9/),11(a,f9.3/),a,i9)')
     *  'time     =',  nint(timeice),
     *  'weirun   =',  weirun,
     *  'nx       =',  n1s,
     *  'nlev     =',  nlev,
     *  'nbed     =',  nbed,
     *  'nz       =',  nz,
     *  'xlef     =',  (d1s(1)   - 0.5*dx0)*1.e-3,
     *  'xrit     =',  (d1s(n1s) + 0.5*dx0)*1.e-3,
     *  'zbot     =',  zbot,
     *  'ztop     =',  ztop,
     *  'sealev   =',  sealev,
     *  'dtantann =',  dtantann,
     *  'dtantjan =',  dtantjan,
     *  'dtseas   =',  dtseas,     ! new..change ginice...
     *  'rco2     =',  rco2,       ! new..change ginice...
     *  'xgl      =',  xgl*1.e-3,
     *  'fgl      =',  fgl,
     *  'igl      =',  igl

c        Write 1-D quantities (h-grid)

c     write(iu,cform) 'xh       :', (xh(i1s(m),j1s(m)),      m=1,n1s)
      write(iu,cform) 'd1s      :', (d1s(m)*1.e-3,           m=1,n1s)
      write(iu,cform) 'h        :', (h(i1s(m),j1s(m)),       m=1,n1s)
      write(iu,cform) 'hs       :', (hs(i1s(m),j1s(m)),      m=1,n1s)
      write(iu,cform) 'hb       :', (hb(i1s(m),j1s(m)),      m=1,n1s)
      write(iu,cform) 'sedim    :', (sedim(i1s(m),j1s(m)),   m=1,n1s)
      write(iu,cform) 'sedimun  :', (sedimun(i1s(m),j1s(m)), m=1,n1s)
      do m=1,n1s
        work(m) = quarryrate(i1s(m),j1s(m))*1.e6
      enddo
      write(iu,cform) 'quarryrat:', (work(m),                m=1,n1s)
      write(iu,cform) 'hw       :', (hw(i1s(m),j1s(m)),      m=1,n1s)
      write(iu,cform) 'topbedeq :', (topbedeq(i1s(m),j1s(m)),m=1,n1s)
      write(iu,cform) 'thetau   :', (thetau(i1s(m),j1s(m)),  m=1,n1s)
      write(iu,cform) 'budgall  :', (budgall(i1s(m),j1s(m)), m=1,n1s)
      write(iu,cform) 'basefrml :', (basefrml(i1s(m),j1s(m)),m=1,n1s)
      write(iu,cform) 'oceanmelt:',(oceanmelt(i1s(m),j1s(m)),m=1,n1s)
      write(iu,cform) 'tw       :', (tw(i1s(m),j1s(m))-tmelt,m=1,n1s)
      do m=1,n1s
        work(m) = 0.
        do k=1,nlev
          work(m) = work(m) 
     *            + heati(i1s(m),j1s(m),k)*dzeta(k)*h(i1s(m),j1s(m))
        enddo
      enddo
      write(iu,cform) 'heati    :', (work(m)/31556926,       m=1,n1s)
      write(iu,cform) 'heatb    :', 
     *                        (heatb(i1s(m),j1s(m))/31556926,m=1,n1s)
      do m=1,n1s
        work(m) = t(i1s(m),j1s(m),nlevp) 
     *          - (tmelt-dtmdh*h(i1s(m),j1s(m)))
      enddo
      write (iu,cform) 'tbme    :', (work(m),                m=1,n1s)

c        Write 1-D quantities (u-grid)

      call flowline (ua, va, work1dv, 1)
      write(iu,cform) 'ua       :', (work1dv(m),             m=0,n1s)
      call flowline (ub, vb, work1dv, 1)
      write(iu,cform) 'ub       :', (work1dv(m),             m=0,n1s)

      do j=0,ny
        do i=0,nx
          if (crhu(i,j).eq.0.) then
            work2d(i,j)  = 0.
          else
            work2d (i,j) = log10(crhu(i,j))
          endif
          if (crhv(i,j).eq.0.) then
            work2da(i,j) = 0.
          else
            work2da(i,j) = log10(crhv(i,j))
          endif
        enddo
      enddo
      call flowline (work2d, work2da, work1dv, 1)
      write(iu,cform) 'crhu     :', (work1dv(m),             m=0,n1s)

c        Write 2-D quantities (h-grid)

      call crossec ( t,          tmelt, h, hs, hb, topbed, sedim, 
     *              't',   iu, 1, 1,  0,nlevp,nz, zbot, ztop)

      call crossec ( tsed,       tmelt, h, hs, hb, topbed, sedim,
     *              'tsed',iu, 2, 1,   1,nsed,nz, zbot, ztop)

      call crossec ( tbed,       tmelt, h, hs, hb, topbed, sedim,
     *              'tbed',iu, 3, 1,   1,nbed,nz, zbot, ztop)

c        Write 2-D quantities (u-grid)

      call crossecuv (u,   v,    0.  , h, hs, hb, topbed, sedim, 
     *                'u',  iu, 1, 2, 0, nlevp, nz, zbot, ztop)

      call flush (iu)

      return
      end


c-----------------------------------------------------------------------

      subroutine crossec (arr,       shift, h, hs, hb, topbed, sedim,
     *                    cname, iu, iflag, igrid, n1,n2,nz, zbot, ztop)

      use comicephys
      use comicegrid

      entry      crossecuv(arru, arrv, shift, h, hs, hb, topbed, sedim,
     *                    cname, iu, iflag, igrid, n1,n2,nz, zbot, ztop)

c     Called by iceplot1d, for flowline profiles.
c
c     Writes a 2-d ascii cross-section (x vs z, at j, for plotting) 
c     of array arr or arru, converting from vertical within-ice,
c     within-sed, or within-bed "sigma" coord to physical z coord. 
c     Accumulate only those points in x,z space in the cross sec
c     (ice, sed or bedrock), and write only those to iu for gf.
c
c     iflag = 1 for ice sheet , 2 for sediment, 3 for bedrock
c     igrid = 1 for h-grid (arr), 2 for u-grid (arru, interp here to h)

c     Uses flowline definition i1s,j1s,n1s in comicegrid.

      dimension arr(nx,ny,n1:n2), 
     *          arru(0:nxp,0:nyp,n1:n2), arrv(0:nxp,0:nyp,n1:n2),
     *          h(nx,ny), hs(nx,ny), hb(nx,ny),
     *          topbed(nx,ny), sedim(nx,ny)
      character*(*) cname

      dimension cross(nx*200),  indcross(nx*200,2)

      parameter (nfine=10000)
      dimension indfine   (nfine,2), weifine(nfine), 
     *          indfinesed(nfine,2), weifinesed(nfine), 
     *          indfinebed(nfine,2), weifinebed(nfine) 
      save indfine,    weifine, 
     *     indfinesed, weifinesed,
     *     indfinebed, weifinebed  

      logical first
      save first
      data first /.true./

c        Set regularly spaced fine resolution index into zeta/zsed/zbed
c        ([ind,wei]fine for ice, [ind,wei]finesed for sedim, 
c         [ind,wei]finebed for bedrock) 

      if (first) then

c       for ice sheet:
        do m=1,nfine
          zfine = (m-.5)/nfine
          do k=0,nlev
            if (zeta(k).le.zfine .and. zeta(k+1).ge.zfine) then
              indfine(m,1) = k 
              indfine(m,2) = k+1 
              weifine(m) = (zeta(k+1)-zfine) / (zeta(k+1)-zeta(k))
              goto 10
            endif
          enddo
          write (6,*) 'Error (crossec): setting indfine, weifine'
          stop
   10     continue
        enddo

c       for sedim:
        do m=1,nfine
          zfine = (m-.5)/nfine
          if (zfine.le.zsed(1)) then
            indfinesed(m,1) = 1 
            indfinesed(m,2) = 1 
            weifinesed(m) = 1.
          else if (zfine.ge.zsed(nsed)) then
            indfinesed(m,1) = nsed 
            indfinesed(m,2) = nsed 
            weifinesed(m) = 1.
          else
            do k=1,nsed-1
              if (zsed(k).le.zfine .and. zsed(k+1).ge.zfine) then
                indfinesed(m,1) = k 
                indfinesed(m,2) = k+1
                weifinesed(m) = (zsed(k+1)-zfine) / (zsed(k+1)-zsed(k))
                goto 20
              endif
            enddo
            write(6,*) 'Error (crossec): setting indfinesed, weifinesed'
            stop
   20       continue
          endif
        enddo

c       for bedrock:
        do m=1,nfine
          zfine = (m-.5)/nfine
          if (zfine.le.zbed(1)) then
            indfinebed(m,1) = 1 
            indfinebed(m,2) = 1 
            weifinebed(m) = 1.
          else if (zfine.ge.zbed(nbed)) then
            indfinebed(m,1) = nbed 
            indfinebed(m,2) = nbed 
            weifinebed(m) = 1.
          else
            do k=1,nbed-1
              if (zbed(k).le.zfine .and. zbed(k+1).ge.zfine) then
                indfinebed(m,1) = k 
                indfinebed(m,2) = k+1
                weifinebed(m) = (zbed(k+1)-zfine) / (zbed(k+1)-zbed(k))
                goto 30
              endif
            enddo
            write(6,*) 'Error (crossec): setting indfinebed, weifinebed'
            stop
   30       continue
          endif
        enddo

        first = .false.
      endif

c        Fill cross and indcross,  array, vertically interpolating 
c        wrt zeta or zbed. npoi is counter for points in the cross-sec.

      npoi = 0
c================
      do m=1,n1s
c================
        i = i1s(m)
        j = j1s(m)
        if ( (h(i,j).gt.0.     .and. iflag.eq.1) .or.
     *       (sedim(i,j).gt.0. .and. iflag.eq.2) .or.
     *       (                       iflag.eq.3) ) then

          if (iflag.eq.1) then
            za = hs(i,j) - h(i,j)
            zb = hs(i,j)
          else if (iflag.eq.2) then
            za = hb(i,j) - sedim(i,j)
            zb = hb(i,j)
          else if (iflag.eq.3) then
            za = topbed(i,j) - bedthick
            zb = topbed(i,j)
          endif

c         if (iflag.eq.1) then
c           dz = (ztop-zbot)/nz
c         else
            dz = 0.
c         endif

          do k=1,nz
            z = zbot + ((k-.5)/nz)*(ztop-zbot)
            if (z.ge.za .and. z.le.zb) then
              zsig = (zb-z)/(zb-za)

c             use regularly spaced fine res zeta-index (for speed)
              ifine = max (1, min (nfine, nint(zsig*nfine + 0.5) ))
              if (iflag.eq.1) then
                ka  = indfine(ifine,1)
                kb  = indfine(ifine,2)
                zwa = weifine(ifine)
              else if (iflag.eq.2) then
                ka  = indfinesed(ifine,1)
                kb  = indfinesed(ifine,2)
                zwa = weifinesed(ifine)
              else if (iflag.eq.3) then
                ka  = indfinebed(ifine,1)
                kb  = indfinebed(ifine,2)
                zwa = weifinebed(ifine)
              endif
              npoi = npoi + 1
              if (npoi.gt.(nx*200)) then
                write (6,*) 'Error (crossec): too many points'
                stop
              endif
              if (igrid.eq.1) then
                cross(npoi) = zwa*arr(i,j,ka) + (1.-zwa)*arr(i,j,kb)
              else if (igrid.eq.2) then
c either:
c               if (j1s(m).eq.j1s(min(m+1,n1s))) then
c                 cross(npoi)=     zwa *0.5*(arru(i,j,ka)+arru(i-1,j,ka))
c    *                       + (1.-zwa)*0.5*(arru(i,j,kb)+arru(i-1,j,kb))
c    *                       - shift
c               else
c                 cross(npoi)=     zwa *0.5*(arrv(i,j,ka)+arrv(i,j-1,ka))
c    *                       + (1.-zwa)*0.5*(arrv(i,j,kb)+arrv(i,j-1,kb))
c    *                       - shift
c               endif
c or:
                zu =       zwa *0.5*(arru(i,j,ka)+arru(i-1,j,ka))
     *               + (1.-zwa)*0.5*(arru(i,j,kb)+arru(i-1,j,kb)) -shift
                zv =       zwa *0.5*(arrv(i,j,ka)+arrv(i,j-1,ka))
     *               + (1.-zwa)*0.5*(arrv(i,j,kb)+arrv(i,j-1,kb)) -shift
                cross(npoi) = sqrt(zu**2 + zv**2)
              endif
              indcross(npoi,1) = m
              indcross(npoi,2) = k
            endif
          enddo


        endif
c==========
      enddo
c==========

c        Write cross section to output file

      write (iu,"(a,':')") cname
      write (iu,'(a,i9)')
     *  'npoi     =',  npoi

      if (npoi.gt.0) then
        do m=1,npoi
          write (iu,"(a,'(',i4,',',i4,') = ',f12.5)")
     *      cname, indcross(m,1), indcross(m,2), cross(m)
        enddo
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine icetab1d (h, hs, hb, t,  
     *                     budgall, basefrml, 
     *                     oceanmelt, calvrate, 
     *                     tsurf, tsurfi, heati, heatb,
     *                     w, wa, maskh,
     *                     sedim, tsed, wsed, heats, 
     *                     topbed, topbedeq, tbed,
     *                     hw, tw, maskwater, 
     *                     u, v, ua, va, ui, vi, ub, vb, 
     *                     uadv, vadv, hu, hv, 
     *                     masku, maskv, muind, mvind, crhu, crhv,
     *                     fracgu, fracgv, thetau, thetav,
     *                     fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *                     sealev, dtantann, dtantjan, dtseas, rco2,
     *                     timeice, dt, weirun,
     *                     iloop, nloopend, nyeartab)
 
c     Writes ascii tabular file, time series of particular quantities,
c     for 1-D runs

      use comicephys
      use comicegrid

      dimension 
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny), 
     *  t(nx,ny,0:nlevp),   
     *  budgall(nx,ny),      basefrml(nx,ny),    
     *  oceanmelt(nx,ny),    calvrate(nx,ny),
     *  tsurf(nx,ny),        tsurfi(nx,ny),
     *  heati(nx,ny,nlev),   heatb(nx,ny), 
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp),
     *  maskh(nx,ny),
     *  sedim(nx,ny),        tsed(nx,ny,nsed),   wsed(nx,ny,nsed),
     *  heats(nx,ny,nsed), 
     *  topbed(nx,ny),       topbedeq(nx,ny),    tbed(nx,ny,nbed),
     *  hw(nx,ny),           tw(nx,ny),          maskwater(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  ui(0:nxp,0:nyp),        vi(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),      vadv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),        hv(0:nxp,0:nyp),
     *  crhu(0:nxp,0:nyp),      crhv(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),     maskv(0:nxp,0:nyp), 
     *  muind(0:nxp,0:nyp),     mvind(0:nxp,0:nyp), 
     *  fracgu(0:nxp,0:nyp),    fracgv(0:nxp,0:nyp), 
     *  thetau(nx,ny),          thetav(nx,ny),
     *  fluxgrdu(0:nxp,0:nyp),  fluxgrdv(0:nxp,0:nyp),
     *  fluxschu(0:nxp,0:nyp),  fluxschv(0:nxp,0:nyp)

      logical firsthead
      data firsthead /.true./

      if (nyeartab.eq.0) return
      if ( .not. 
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyeartab),dt)).lt.dt
     *     .or. iloop.eq.1 .or. iloop.eq.nloopend .or. nyeartab.eq.-1
     *     )
     *   ) return 

      j = (ny+1)/2
      iu = iutab

#if defined EISLINE && defined LINEE
      x500 = 24.e3   ! b,c,d,e
#else
      x500 = 500.e3
#endif
      do i=1,nx-1
        if (x500.le.0.5*(xh(i,j)+xh(i+1,j))) then
          i500 = i
          goto 50
        endif
      enddo
      i500 = nx
   50 continue
      h500  = h(i500,j)
      hs500 = hs(i500,j)

c        Total ice volume, non-ocean ice volume (average thicknesses)

      tota    = 0.
      toth    = 0.
      tothnf  = 0.
      totsed  = 0.
      do i=1,nx
        tota = tota +        darea(i,j)
        toth = toth + h(i,j)*darea(i,j)
        if (maskwater(i,j).ne.1) tothnf = tothnf + h(i,j)*darea(i,j)
        totsed = totsed + sedim(i,j)*darea(i,j)
      enddo
      toth   = toth  /tota
      tothnf = tothnf/tota
      totsed = totsed  /tota

c        Position of rightmost grounding line (u-grid) (search for 
c        g.l. from right to left), and thicknesses and elevs each side

      igl = 0
      fgl = 0.
      xgl  = 0.
      hgla = 0.
      hglb = 0.
      hglu = 0.
      zgflot = 0.
      zfluxgrd = 0.
      zfluxsch = 0.
      ztheta = 0.
      do i=nx-1,1,-1
        if (fracgu(i,j).gt.0.) then 
          igl = i
          fgl = fracgu(i,j)
          ztheta = thetau(min(i+2,nx),j)
          xgl = (1.-fgl)*xh(i,j)+fgl*xh(i+1,j)
          hgla = h(i,j)
          hglb = h(i+1,j)
c         hglu = (1.-fgl)*h(i,j)+fgl*h(i+1,j)                       !old
c         hglu = (sealev - (     0.5*hb(i,j)+0.5*hb(i+1,j))) / rhor !888
          hglu = (sealev - ((1.-fgl)*hb(i,j)+fgl*hb(i+1,j))) / rhor !888
          zgflot = h(i,j)-(rholiq/rhoice)*(sealev-hb(i,j))
          zfluxgrd = fluxgrdu(i,j)
          zfluxsch = fluxschu(i,j)
          go to 100
        endif
      enddo
  100 continue

      igsc=0
      do i=1,nx
        if (muind(i,j).eq.-1) then
          igsc = i
          goto 110
        endif
      enddo
  110 continue

c        Total surface budget between ice divide (at max hs) and 
c        rightmost grounding line (bgl)

      bgl = 0.
      isummit = (nx+1)/2
      if (igl.gt.0) then
        hsummit = -1.e20
        do i=1,nx
          if (h(i,j).gt.0. .and. hs(i,j).gt.hsummit) then  
            hsummit = hs(i,j)
            isummit = i
          endif
        enddo
        do i = min(isummit,igl), max(isummit,igl)
          zfrac = 1.
          if (i.eq.isummit .and. i.gt.1) zfrac = 0.5
          bgl = bgl + zfrac*dx(i,j)*budgall(i,j)
        enddo
      endif

c        Write header line(s) first call this run

      if (iloop.eq.1) then
        nspw = max (1,nint(nyeartab/dt))
        nwrite = nloopend / nspw
        if (nspw.gt.1) nwrite = nwrite + 1
        write (iu,'(a,i8,a,i5)') 'nwritetab=',nwrite,'  nx=',nx
        write (iu,'(3a)')
     *    '      time   sl dtanta dtantj dtseas   rco2   hsum',
     *    ' toth  igl igsc    xgl    fgl   hglu   hgla   hglb',
     *    '   zgflot  fluxgrd  fluxsch      bgl    theta'
      endif

      write (iu,'(f10.1, i5, 3f7.2, f7.3, f7.1,
     *            3i5, f7.1, f7.3, 3f7.1,
     *            f9.1, 3i9, f9.5)')
     *  timeice, 
     *  nint(sealev), dtantann, dtantjan, dtseas, rco2, h(isummit,j), 
     *  nint(toth), igl, igsc, xgl*1.e-3, fgl, hglu, hgla, hglb, 
     *  zgflot, nint(zfluxgrd), nint(zfluxsch), nint(bgl), ztheta

      call flush (iu)

c     if (.true.) then
c     if (.false.) then
c      iu = 80
c      write (iu,'(/79("=")/a,f12.1)') 'time=',timeice
c      write(iu,'(a,5x,2i15)')  'i   =',(i,      i=igl,igl+1)  
c      write(iu,'(a,5x,2f15.2)')'xh  =',(1.e-3*xh(i,j), i=igl,igl+1)  
c      write(iu,'(a,5x, f15.2)')'xgl =',1.e-3*xgl
c      write(iu,'(a,5x,2f15.2)')'hs  =',(hs(i,j),i=igl,igl+1)  
c      write(iu,'(a,5x,2f15.2)')'h   =',(h(i,j) ,i=igl,igl+1)  
c      write(iu,'(a,5x,2f15.2)')'hb  =',(hb(i,j),i=igl,igl+1)  
c      write(iu,'(a,5x,2f15.8)')'hw  =',(hw(i,j) ,i=igl,igl+1)  
c      write(iu,'(a,5x,2i15)')'mask=',(maskwater(i,j) ,i=igl,igl+1)  
c      write(iu,'(a,5x,2f15.2)')'zflo=',
c    *           ((h(i,j)-(rholiq/rhoice)*(sealev-hb(i,j))),i=igl,igl+1)
c
c      write(iu,'(a,3f15.2)') 'hu  =',(hu(i,j),i=igl-1,igl+1)  
c      write(iu,'(a,3f15.2)') 'uadv=',(uadv(i,j)/dyu(i,j),i=igl-1,igl+1)
c      write(iu,'(a,3e15.4)') 'crhu=',(crhu(i,j),i=igl-1,igl+1)
c      call flush (iu)
c     endif

      return
      end

c-----------------------------------------------------------------------

      subroutine icetab2d (h, hw, maskwater,
     *                     sealev, dtantann, dtantjan, dtseas, rco2,
     *                     ecc, obl, prec, facice, facorb, facco2,
     *                     timeice, dt, weirun,
     *                     iloop, nloopend, nyeartab, ifrest)
 
c     Writes ascii tabular file, time series of particular quantities,
c     for 2-D runs

      use comicephys
      use comicegrid

      dimension 
     *  h(nx,ny), hw(nx,ny), maskwater(nx,ny)

      if (nyeartab.eq.0) return
      if ( .not. 
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyeartab),dt)).lt.dt
     *     .or. iloop.eq.1 .or. iloop.eq.nloopend
     *     )
     *   ) return 

      iu = iutab

c        Write header line(s) first call, if not a restart

      if (iloop.eq.1) then
        nspw = max (1,nint(nyeartab/dt))
        nwrite = nloopend / nspw
        if (nspw.gt.1) nwrite = nwrite + 1
        write (iu,'(a,i8)') 'nwritetab=',nwrite
        write (iu,'(5a)')
     *    '      time weirun sealev dtanta dtantj dtseas   rco2',
     *    '     ecc     obl    prec facice facorb facco2',
     *    '   toti(km3)  totig(km3)  totif(km3)',
     *    '   tota(km2)  totag(km2)  totaf(km2)',
     *    '    h(m)'
      endif

      toti   = 0.
      tota   = 0.
      totif = 0.
      totaf = 0.
      do j=1,ny
        do i=1,nx
          if (h(i,j).ne.0.) then 
            toti = toti + darea(i,j)*h(i,j)
            tota = tota + darea(i,j)
            if (maskwater(i,j).eq.1) then
              totif = totif + darea(i,j)*h(i,j)
              totaf = totaf + darea(i,j)
            endif
          endif
        enddo
      enddo

      write (iu,
     *       '(i10, f7.3, f7.1, 3f7.2, f7.3, 
     *         f8.5,f8.3,f8.3,  f7.3,f7.1,f7.3, 
     *         6i12, 
     *         f8.1)')
     *  nint(timeice), weirun, sealev, dtantann,dtantjan,dtseas, rco2,
     *  ecc, obl*180./pi, prec*180./pi, facice, facorb, facco2,
     *  nint(toti*1.e-9), nint ((toti-totif)*1.e-9), nint(totif*1.e-9),
     *  nint(tota*1.e-6), nint ((tota-totaf)*1.e-6), nint(totaf*1.e-6),
     *  toti/max(tota,.001)

      call flush (iu)

      return
      end

c-----------------------------------------------------------------------

      subroutine writehis (h, hs, hb, t, 
     *                     budgall, basefrml, 
     *                     oceanmelt, calvrate, 
     *                     tsurf, tsurfi, heati, heatb,
     *                     w, wa, sedim, tsed, wsed, heats, 
     *                     topbed, topbedeq, tbed,
     *                     hw, tw, maskwater,
     *                     u, v, ua, va, ub, vb, uadv, vadv,
     *                     timeice, dt, nwrit,
     *                     iloop, nloopend, nyearhis, hislist, nhislist,
     *                     ifrest)

c     Write to history file

      use comicephys
      use comicegrid

      dimension 
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny), 
     *  t(nx,ny,0:nlevp),    
     *  budgall(nx,ny),      basefrml(nx,ny),    
     *  oceanmelt(nx,ny),    calvrate(nx,ny),
     *  tsurf(nx,ny),        tsurfi(nx,ny),
     *  heati(nx,ny,nlev),   heatb(nx,ny), 
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp),
     *  sedim(nx,ny),        tsed(nx,ny,nsed),   wsed(nx,ny,nsed),
     *  heats(nx,ny,nsed),
     *  topbed(nx,ny),       topbedeq(nx,ny),    tbed(nx,ny,nbed),
     *  hw(nx,ny),           tw(nx,ny),          maskwater(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),      vadv(0:nxp,0:nyp)

      integer hislist(1000)

      dimension work(nx,ny)

c===================
#if defined NETCDF
c===================
#include <netcdf.inc>
      integer*4
     * lenattr, nx4, ny4, nlev4, ival4, jval4,
     * varid, rcode, ncerr, ncid, 
     * xid, yid, levid, timid, minmaxid,
     * minmax(2), ndim, idim(4), start(4), count(4), 
     * n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,
c      following are from netcdf.inc...need to convert to int*4,
c      and use *4 variables as arguments in netcdf calls.
     * ncglobal4, ncunlim4,
     * nf_float4, nf_short4, nf_clobber4, nf_write4, nf_noerr4
      parameter (ncglobal4   = ncglobal,   ncunlim4  = ncunlim,
     *           nf_float4   = nf_float,   nf_short4 = nf_short,
     *           nf_clobber4 = nf_clobber, nf_write4 = nf_write,
     *           nf_noerr4   = nf_noerr)
      save varid, rcode, ncid,
     *     xid, yid, levid, timid, minmaxid,
     *     minmax, idim, start, count,
     *     n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13
      data minmax /1,2/
      data n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13
     *     /1, 2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12, 13/
      character cruntitle*80, cvar*16, cunits*16
      dimension ztmp(1000) !Jer: to allow large output files
# if !defined O_ism
c     parameter (chist = 'history.nc')
      parameter (chist = 'fort.92.nc')
# endif     
c====
#else
c====
# if !defined O_ism 
c     parameter (chist = 'history')
      parameter (chist = 'fort.92')
# endif
c=====
#endif
c=====
#if defined O_ism
      character cyr*80
#endif

c        hislist (namelist, if entered) has precedence over nyearhis

      if (nhislist.ne.0) then
        do m=1,nhislist
          if (abs(timeice-hislist(m)).lt.0.5*dt) goto 100
        enddo
        return
  100   continue
      else if (nyearhis.gt.0) then
        if (
     *    .not.( mod(abs(timeice)+0.5*dt,max(float(nyearhis),dt)).lt.dt
     *           .or. iloop.eq.1 .or. iloop.eq.nloopend
     *          )
     *     ) return 
#if defined HEINO
        if (nint(timeice).lt.150000 .and. nyearhis.lt.500) return
#endif
      else
        return
      endif

c-------------------------
      if (nwrit.gt.0) then
c-------------------------
#if defined NETCDF 
        rcode = nf_open(chist, nf_write4, ncid)
c       next write needed to avoid ifort/netcdf bug (?)
        write (6,*) 'Opening existing Netcdf history file'
#else
#if defined SGI || (defined SUN && defined F77)
        write (6,*) 'Opening existing Netcdf history file'
        open (iuhis, file=chist, status='old',form='unformatted',
     *        access='append', iostat=ios, err=8000)
#else
        open (iuhis, file=chist, status='old',form='unformatted',
     *        position='append', iostat=ios, err=8000)
#endif
#endif
c---------
      else
c---------
#if defined NETCDF 
#if defined O_ism
        !set name of file, to be renewed every time
	!a new UVic runstep starts
        write (cyr, '(i10)') nint(abs(timeice))
        do i=1,lenchr(cyr)
           if (cyr(1:1).eq.' ') cyr = cyr(2:)
        enddo
        if (timeice.lt.0.) cyr = cyr(1:lenchr(cyr)) // 'm'
        chist = trim(isname)// trim(cyr)//'.nc'	
#endif

        rcode = nf_create (chist, nf_clobber4, ncid)
        write (6,*) 'Opening new Netcdf history file'
#else
        open (iuhis, file=chist, status='unknown')
        close(iuhis, status='delete')
        open (iuhis, file=chist, status='new', form='unformatted')
#endif

c...................
#if defined NETCDF 
c...................

c          define title, etc

        cruntitle = 'Ice Sheet Test'
        lenattr = lenchr(cruntitle)
        rcode = nf_put_att_text (ncid, ncglobal4, 'title',
     *                             lenattr, cruntitle)

        rcode = nf_put_att_text (ncid, ncglobal4, 'horiz_grid',
     *                           n7, 'unknown')

        rcode = nf_put_att_text (ncid, ncglobal4, 'projection',
     *                           n6, 'latlon')

c          x dimension

        nx4 = nx
        rcode = nf_def_dim (ncid, 'x', nx4, xid)
        rcode = nf_def_var (ncid, 'x', nf_float4, n1, xid, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n10, 'x distance')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n2, 'km')
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
        do i=1,nx
          ztmp(i) = xh(i,ny/2)*.001
        enddo
        rcode = nf_put_vara_double (ncid, varid, n1, nx4, ztmp)
        rcode = nf_redef (ncid)

c          y dimension

        ny4 = ny
        rcode = nf_def_dim (ncid, 'y', ny4, yid)
        rcode = nf_def_var (ncid, 'y', nf_float4, n1, yid, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n10, 'y distance')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n2, 'km')
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
        do j=1,ny
          ztmp(j) = yh(nx/2,j)*.001
        enddo
        rcode = nf_put_vara_double (ncid, varid, n1, ny4, ztmp)
        rcode = nf_redef (ncid)

c          z dimension

        nlev4 = nlev
        rcode = nf_def_dim (ncid, 'z', nlev4, levid)
        rcode = nf_def_var (ncid, 'z', nf_float4, n1, levid, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n8, 'z height')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n5, 'sigma')
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
        rcode = nf_put_vara_double (ncid, varid, n1, nlev4, zeta)
        rcode = nf_redef (ncid)

c          Time dimension 

        rcode = nf_def_dim (ncid, 'time', ncunlim4, timid)
        rcode = nf_def_var (ncid, 'time', nf_float4, n1, timid, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n4, 'time')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n5, 'years')
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n5, 'f12.3')

c         Variable dimensions and attributes

c       For 2D, no time:

        ndim    = 2
        idim(1) = xid
        idim(2) = yid

        cvar = 'alatd'
        ival4 = lenchr(cvar)
        cunits = 'deg N'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'alond'
        ival4 = lenchr(cvar)
        cunits = 'deg E'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

c       For 2D + time:

        ndim    = 3
        idim(1) = xid
        idim(2) = yid
        idim(3) = timid

        cvar = 'h'
        ival4 = lenchr(cvar)
        cunits = 'm'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'hb'
        ival4 = lenchr(cvar)
        cunits = 'm'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'budgall'
        ival4 = lenchr(cvar)
        cunits = 'm/y'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'sedim'
        ival4 = lenchr(cvar)
        cunits = 'm'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'maskwater'
        ival4 = lenchr(cvar)
        cunits = 'flag'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'thomol'
        ival4 = lenchr(cvar)
        cunits = 'C'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'utop'
        ival4 = lenchr(cvar)
        cunits = 'm/y'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'vtop'
        ival4 = lenchr(cvar)
        cunits = 'm/y'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'ubot'
        ival4 = lenchr(cvar)
        cunits = 'm/y'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'vbot'
        ival4 = lenchr(cvar)
        cunits = 'm/y'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')	
        cvar = 'mext'
        ival4 = lenchr(cvar)
        cunits = 'm^2'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')
	
        cvar = 'mdur'
        ival4 = lenchr(cvar)
        cunits = 'days'
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cname)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')	
c       For 3D + time: (none currently)
c
c       ndim    = 4
c       idim(1) = xid
c       idim(2) = yid
c       idim(3) = levid
c       idim(4) = timid

c         Leave define mode

        rcode = nf_enddef (ncid)
c.....
#endif
c.....
c----------
      endif
c----------

#if defined NETCDF 

c        Add current time to time dimension

      nwrit = nwrit + 1
      ival4 = nwrit 
      rcode = nf_inq_varid (ncid, 'time', varid)
      rcode = nf_put_var1_double (ncid, varid, ival4, timeice)

c        Write fields

c     For 2D, no time or + time:

      start(1) = 1
      count(1) = nx
      start(2) = 1
      count(2) = ny
      start(3) = nwrit
      count(3) = 1

      if (nwrit.eq.1) then
        rcode = nf_inq_varid (ncid, 'alatd', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, alatd)

        rcode = nf_inq_varid (ncid, 'alond', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, alond)
      endif

      rcode = nf_inq_varid (ncid, 'h', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, h)

      rcode = nf_inq_varid (ncid, 'hb', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, hb)

      rcode = nf_inq_varid (ncid, 'budgall', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, budgall)

      rcode = nf_inq_varid (ncid, 'sedim', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, sedim)

      call itor (maskwater, work, nx*ny)
      rcode = nf_inq_varid (ncid, 'maskwater', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = t(i,j,nlevp) + dtmdh*h(i,j) - tmelt
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'thomol', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = u(i,j,0)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'utop', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = v(i,j,0)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'vtop', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = ub(i,j)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'ubot', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = vb(i,j)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'vbot', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)     

      do j=1,ny
        do i=1,nx
          work(i,j) = mext(i,j)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'mext', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)
      
      do j=1,ny
        do i=1,nx
          work(i,j) = mdur(i,j)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'mdur', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)        
c     For 3D + time, additional: (none currently)
c
c     start(3) = 1
c     count(3) = nlev
c     start(4) = nwrit
c     count(4) = 1

c        Return to define mode

      rcode = nf_redef (ncid)

      rcode = nf_close (ncid)                ! to flush, re-opened above

#else

      nwrit = nwrit + 1

      if (nwrit.eq.1) then
        call writehis2 
     *    ('alatd',  'deg N', alatd, 1,nx,1,ny, timeice,iuhis,1.)
        call writehis2 
     *    ('alond',  'deg E', alond, 1,nx,1,ny, timeice,iuhis,1.)
      endif

      call writehis2 
     *  ('h',  'm', h, 1,nx,1,ny, timeice,iuhis,1.)
      call writehis2 
     *  ('hb', 'm', hb, 1,nx,1,ny, timeice,iuhis,1.)
      call writehis2 
     *  ('budgall', 'cm/y', budgall, 1,nx,1,ny,timeice,iuhis,100.)
      call writehis2 
     *  ('sedim', 'm', sedim, 1,nx,1,ny,timeice,iuhis,1.)
      call itor (maskwater, work, nx*ny)
      call writehis2 
     *  ('maskwate', ' ', work, 1,nx,1,ny,timeice,iuhis,1.)

      do j=1,ny
        do i=1,nx
          work(i,j) = t(i,j,nlevp) + dtmdh*h(i,j) - tmelt
        enddo
      enddo
      call writehis2 
     *  ('thomol', 'C', work, 1,nx,1,ny,timeice,iuhis,1.)

      call writehis2 
     *  ('utop', 'm/y', u(0,0,0), 0,nxp,0,nyp,timeice,iuhis,1.)
      call writehis2 
     *  ('vtop', 'm/y', v(0,0,0), 0,nxp,0,nyp,timeice,iuhis,1.)
      call writehis2 
     *  ('ubot', 'm/y', ub, 0,nxp,0,nyp,timeice,iuhis,1.)
      call writehis2 
     *  ('vbot', 'm/y', vb, 0,nxp,0,nyp,timeice,iuhis,1.)

      close (iuhis)                          ! to flush, re-opened above
#endif

      return

c        Fatal error for history file opening

 8000 write (ioterm, '(/3a)') 
     *  '*** Error in opening existing history file ', chist,
     *  ' in append mode'
      stop

      end

c-----------------------------------------------------------------------

      subroutine writehis2 (cfieldin, cunitsin, arr, ix1, ix2, iy1, iy2,
     *                      time, iu, scale)

c Writes 2-D ice-grid field (arr) to binary history file iu.
c Field is multiplied by scale.
c Write header including nx,ny,iglob,jglob, and just write arr
c for ice-grid region (not global).

      use comicegrid

      character*(*) cfieldin, cunitsin
      dimension arr(ix1:ix2,iy1:iy2)

      character*8 cfield, cunits
      dimension worka(nx,ny)
      parameter (vershis = 2.0)

      cfield = cfieldin
      cunits = cunitsin

      do j=1,ny
        do i=1,nx
          worka(i,j) = arr(i,j)*scale
        enddo
      enddo

      write (iu) vershis, cfield, cunits
      write (iu) nx, ny
      write (iu) time
      write (iu) worka

      return
      end

c-----------------------------------------------------------------------


      subroutine mascalc (arr, totv, tota)

c     Calculates total (domain) damount of variable arr (totv),
c     and total area with non-zero arr (tota)  

      use comicephys
      use comicegrid

      dimension arr(nx,ny)

      totv  = 0.
      tota  = 0.
      do j=1,ny
        do i=1,nx
          if (arr(i,j).ne.0.) then 
            totv = totv + darea(i,j)*arr(i,j)
            tota = tota + darea(i,j)
          endif
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

#if defined HEINO

      subroutine iceshowheino (h, hs, hb, t, heatb, u, v, ub, vb, 
     *                         timeice, dt, iloop, nloopend)

c     Writes output files for HEINO ISMIP experiment

      use comicephys
      use comicegrid

      dimension
     *  h(nx,ny), hs(nx,ny), hb(nx,ny),
     *  t(nx,ny,0:nlevp), heatb(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp)

      parameter (npoi = 7)
      dimension xpoi(npoi), ixpoi(npoi),
     *          hpoi(npoi), tpoi(npoi), bpoi(npoi)  
      save xpoi, ixpoi
      data xpoi /3900., 3800., 3700., 3500, 3200., 2900., 2600./

      data iuheivol, iuheiamel, iuheihsed, iuheitsed, iuheiased
     *    /     150,       151,       152,       153,       154 /
      data iuheihpoi, iuheitpoi, iuheibpoi
     *    /      160,       170,       180 /

      logical first
      save first
      data first /.true./

      parameter (nyearheino = 1)

      if (nyearheino.eq.0) return

c        First timestep of whole run (so not restarts), set grid indices
c        for indiv pts 1-7, and write header line to time series files

      if (iloop.eq.1) then
        j = (ny+1)/2
        do ipoi=1,npoi
          do i=1,nx
            if ( abs (.001*(xh(i,j)-xh(1,j)) - xpoi(ipoi)) .lt. .001 )
     *        ixpoi(ipoi) = i
          enddo
c         write (ioterm,*) ipoi,ixpoi(ipoi)
        enddo

        nspw = max (1,nint(nyearheino/dt))
        nwrite = nloopend / nspw 
        if (nspw.gt.1) nwrite = nwrite + 1
        write (iu,'(a,i8)') 'nwriteplot=',nwrite
 
        if (ny.eq.1) then
          if1d = 1
        else
          if1d = 0
        endif
        write (iuheivol, "('nt=',i8,'  if1d=',i2)") nwrite,if1d
        write (iuheiamel,"('nt=',i8,'  if1d=',i2)") nwrite,if1d
        write (iuheihsed,"('nt=',i8,'  if1d=',i2)") nwrite,if1d
        write (iuheitsed,"('nt=',i8,'  if1d=',i2)") nwrite,if1d
        write (iuheiased,"('nt=',i8,'  if1d=',i2)") nwrite,if1d
        do ipoi=1,npoi
          write (iuheihpoi+ipoi-1,"('nt=',i8,'  if1d=',i2)") nwrite,if1d
          write (iuheitpoi+ipoi-1,"('nt=',i8,'  if1d=',i2)") nwrite,if1d
          write (iuheibpoi+ipoi-1,"('nt=',i8,'  if1d=',i2)") nwrite,if1d
        enddo
        first = .false.
      endif

      if (.not.(mod(abs(timeice)+.5*dt,max(float(nyearheino),dt)).lt.dt)
     *   ) return

c        Areal averages

      zvol = 0.
      zmmel = 0.
      zhsed = 0.
      ztsed = 0.
      zmsed = 0.
      zased = 0.
      do j=1,ny
        do i=1,nx
          zdar = dx(i,j)*dy(i,j)
          zthomol = t(i,j,nlevp) + dtmdh*h(i,j)

c         all ice sheet:
          if (h(i,j).gt.0.) then
            zvol = zvol + h(i,j)*zdar
            if (zthomol.ge.tmelt-.5*tramp) zmmel = zmmel + zdar
          endif

c         sediment area: 
          if (maskheino(i,j).eq.2) then
            zhsed = zhsed + h(i,j) *zdar 
            ztsed = ztsed + zthomol*zdar
            if (zthomol.ge.tmelt-.5*tramp)  zmsed = zmsed + zdar
            zased = zased +  zdar
          endif
        enddo
      enddo
      zhsed = zhsed/zased
      ztsed = ztsed/zased

c       Indiv points

      if (nx.eq.81 .and. ny.eq.81) then
        j = (ny+1)/2
        do ipoi=1,npoi
          i = ixpoi(ipoi)  
          hpoi(ipoi) = h(i,j)
          tpoi(ipoi) = t(i,j,nlevp) + dtmdh*h(i,j)
          bpoi(ipoi) = heatb(i,j)/31556926              ! J/m2/y to W/m2
        enddo
      endif

c        Write time series

      write (iuheivol , '(2e14.6)') timeice, zvol *1.e-15 
      write (iuheiamel, '(2e14.6)') timeice, zmmel*1.e-12 

      write (iuheihsed, '(2e14.6)') timeice, zhsed*1.e-3
      write (iuheitsed, '(2e14.6)') timeice, ztsed
      write (iuheiased, '(2e14.6)') timeice, zmsed*1.e-12

      do ipoi=1,npoi
        write (iuheihpoi+ipoi-1, '(2e14.6)') timeice, hpoi(ipoi)
        write (iuheitpoi+ipoi-1, '(2e14.6)') timeice, tpoi(ipoi)
        write (iuheibpoi+ipoi-1, '(2e14.6)') timeice, bpoi(ipoi)
      enddo

      call flush (iuheivol)
      call flush (iuheiamel)
      call flush (iuheihsed)
      call flush (iuheitsed)
      call flush (iuheiased)
      do ipoi=1,npoi
        call flush (iuheihpoi+ipoi-1)
        call flush (iuheitpoi+ipoi-1)
        call flush (iuheibpoi+ipoi-1)
      enddo

      return
      end

#endif

c-----------------------------------------------------------------------

#if defined HOM

      subroutine iceshowhom (h, hs, hb, u, v, ub, vb, wa, dfu, dfv,
     *                       timeice, dt, iloop, nloopend)

c     Writes output files for HOM ISMIP experiment

      use comicephys
      use comicegrid

      dimension
     *  h(nx,ny), hs(nx,ny), hb(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),     wa(nx,ny,0:nlevp),
     *  dfu(0:nxp,0:nyp),       dfv(0:nxp,0:nyp)

      character*80 cfile
      data iuhom /191/

      if (.not.(iloop.eq.nloopend)) return

      idomainkm = nint ((xh(nx-1,1) - xh(1,1))*1.e-3)

#if defined HOMC 

      write (cfile,'(a,i3.3,a)') 'dpo1c', idomainkm, '.txt'
      open (iuhom, file=cfile, status='unknown')
      close (iuhom, status='delete')
      open (iuhom, file=cfile, status='new')

      zutopmax = -1.e20
      zvtopmax = -1.e20
      zwtopmax = -1.e20
      zubotmax = -1.e20
      zvbotmax = -1.e20
      zdfumax = -1.e20
      zdfvmax = -1.e20
      zutopmin = 1.e20
      zvtopmin = 1.e20
      zwtopmin = 1.e20
      zubotmin = 1.e20
      zvbotmin = 1.e20
      zdfumin = 1.e20
      zdfvmin = 1.e20

      do j=1,ny-1
        jm = j-1
        if (jm.eq.0) jm = ny-2
        jwa = j
        if (jwa.eq.1) jwa = ny-1
        do i=1,nx-1
          im = i-1
          if (im.eq.0) im = nx-2
          iwa = i
          if (iwa.eq.1) iwa = nx-1

          zx = float(i-1)/float(nx-2)
          zy = float(j-1)/float(ny-2)

          zutop = 0.5* (u(i,j,0) + u(im,j,0))
          zvtop = 0.5* (v(i,j,0) + v(i,jm,0))
          zwtop = wa(iwa,jwa,0)
          zubot = 0.5* (ub(i,j) + ub(im,j))
          zvbot = 0.5* (vb(i,j) + vb(i,jm))
          zdfu  = -0.5* (dfu(i,j) + dfu(im,j)) / 1.e3   ! kPa
          zdfv  = -0.5* (dfv(i,j) + dfv(i,jm)) / 1.e3   ! kPa
          zdp = 0.
          write (iuhom,'(2f7.4,8e14.6)')
     *      zx, zy, zutop, zvtop, zwtop, zubot, zvbot, zdfu, zdfv, zdp

          zutopmax = max (zutopmax, zutop)
          zvtopmax = max (zvtopmax, zvtop)
          zwtopmax = max (zwtopmax, zwtop)
          zubotmax = max (zubotmax, zubot)
          zvbotmax = max (zvbotmax, zvbot)
          zdfumax = max (zdfumax, zdfu)
          zdfvmax = max (zdfvmax, zdfv)
          zutopmin = min (zutopmin, zutop)
          zvtopmin = min (zvtopmin, zvtop)
          zwtopmin = min (zwtopmin, zwtop)
          zubotmin = min (zubotmin, zubot)
          zvbotmin = min (zvbotmin, zvbot)
          zdfumin = min (zdfumin, zdfu)
          zdfvmin = min (zdfvmin, zdfv)
        enddo
      enddo
      close (iuhom)

       write (199,'(a)') 'HOMC min,max:'
       write (199,'(14x,7e14.6)')
     *  zutopmin,zvtopmin,zwtopmin,zubotmin,zvbotmin,zdfumin,zdfvmin 
       write (199,'(14x,7e14.6)')
     *  zutopmax,zvtopmax,zwtopmax,zubotmax,zvbotmax,zdfumax,zdfvmax 

#endif

#if defined HOMD 

      write (cfile,'(a,i3.3,a)') 'dpo1d', idomainkm, '.txt'
      open (iuhom, file=cfile, status='unknown')
      close (iuhom, status='delete')
      open (iuhom, file=cfile, status='new')

      zutopmax = -1.e20
      zwtopmax = -1.e20
      zubotmax = -1.e20
      zdfumax = -1.e20
      zutopmin = 1.e20
      zwtopmin = 1.e20
      zubotmin = 1.e20
      zdfumin = 1.e20

      j = (ny+1)/2
      do i=1,nx-1
        im = i-1
        if (im.eq.0) im = nx-2
        iwa = i
        if (iwa.eq.1) iwa = nx-1

        zx = float(i-1)/float(nx-2)
        zutop = 0.5* (u(i,j,0) + u(im,j,0))
        zwtop = wa(iwa,j,0)
        zubot = 0.5* (ub(i,j) + ub(im,j))
        zdfu  = -0.5* (dfu(i,j) + dfu(im,j)) / 1.e3   ! kPa
        zdp = 0.
        write (iuhom,'(f7.4,5e14.6)')
     *    zx, zutop, zwtop, zubot, zdfu, zdp

        zutopmax = max (zutopmax, zutop)
        zwtopmax = max (zwtopmax, zwtop)
        zubotmax = max (zubotmax, zubot)
        zdfumax = max (zdfumax, zdfu)
        zutopmin = min (zutopmin, zutop)
        zwtopmin = min (zwtopmin, zwtop)
        zubotmin = min (zubotmin, zubot)
        zdfumin = min (zdfumin, zdfu)
      enddo
      close (iuhom)

      write (199,'(a)') 'HOMD min,max:'
      write (199,'(14x,4e14.6)') zutopmin, zwtopmin, zubotmin, zdfumin
      write (199,'(14x,4e14.6)') zutopmax, zwtopmax, zubotmax, zdfumax

#endif

      return
      end

#endif
