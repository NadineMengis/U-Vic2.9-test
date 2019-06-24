! source file: /net/mare/home1/eby/as/ism/icebed.F
      subroutine bedrock (h, hb, hw, topbed, topbedeq, sedim, sedimeq,
     *                    equiload, dt)

c        Steps bedrock elevation topbed through one timestep, then
c        adjusts hb. Computes deflection of elastic
c        lithosphere +fluid mantle (deflect), then computes
c        time-dependent relaxation of topbed to that profile
c        either locally (kasth=0) or for asthenospheric channel flow.

      use comicephys
      use comicegrid

c        h      = ice-sheet thickness (m) (modified)
c        hb     = bedrock+sedim elevation (m) (modified)
c        topbed = bedrock elevation, not including sedim (m) (modified)
c        topbedeq = equilibrium bedrock elevation (m) (supplied)
c        sedim = sediment thickness (m) (supplied)
c        sedimeq = equil. sed. thickness corresp. to topbedeq (m) (supp)
c        equiload = equil. load on bedrock (m) (supp)
c        dt = bedrock-model timestep (yr) (supplied)

      dimension
     *  h(nx,ny),         hb(nx,ny),        hw(nx,ny),
     *  topbed(nx,ny),    topbedeq(nx,ny),
     *  sedim(nx,ny),     sedimeq(nx,ny),
     *  equiload(nx,ny)

c     kasth=0 for local asthenospheric relax, 1 for channel flow:
c     taulocal = relaxation e-folding time (for kasth=0):
c     asthflow = diffusive bedrock flow coeff (m2/y) (for kasth=1):
      parameter (kasth = 0, taulocal = 3000.)   ! 5000.
c     parameter (kasth = 1, asthflow = 1.e8)

      dimension zmass(nx,ny), deflect(nx,ny), za(nx,ny)

c>>>>>>>>>>>>>>>>>>>

c>>>>>>>>>>>>>>>>>>>

c        For elastic lithospheric flexure (2D: Brotchie and Silvester,
c        1969,JGR,74,22,5240, and 1D: Turcotte and Schubert,1982,p.125)

      parameter (cruststiff = 1.e25) ! Huybrechts, cf. 3.76e23 B+S (N m)

c     for parallelization:
      parameter (nseg=4)
      dimension deflectseg(nx,ny,nseg)

      save  crustlen, crustcoef

c      Kelvin Function kei (0th order, from IMSL, for 2D lithosphere):
c      akeic is stored coarse-resol (dkeic), linearly interpolated to
c      finer-resol akei (dkei) in first call:
       parameter (xkei=15., dkeic = .05,  nkeic=xkei/dkeic + .001,
     *                      dkei  = .005, nkei =xkei/dkei  + .001)
       dimension akeic(0:nkeic), akei(0:nkei)
       save akei
       data akeic /
     * -0.7854,-0.7828,-0.7769,-0.7684,-0.7581,-0.7463,-0.7331,-0.7189,
     * -0.7038,-0.6880,-0.6716,-0.6547,-0.6374,-0.6199,-0.6022,-0.5843,
     * -0.5664,-0.5484,-0.5305,-0.5127,-0.4950,-0.4775,-0.4601,-0.4430,
     * -0.4262,-0.4096,-0.3933,-0.3773,-0.3617,-0.3464,-0.3314,-0.3168,
     * -0.3026,-0.2887,-0.2752,-0.2621,-0.2494,-0.2371,-0.2251,-0.2136,
     * -0.2024,-0.1916,-0.1812,-0.1711,-0.1614,-0.1521,-0.1431,-0.1345,
     * -0.1262,-0.1183,-0.1107,-0.1034,-0.0964,-0.0898,-0.0834,-0.0774,
     * -0.0716,-0.0661,-0.0608,-0.0558,-0.0511,-0.0466,-0.0424,-0.0384,
     * -0.0346,-0.0310,-0.0276,-0.0244,-0.0214,-0.0186,-0.0160,-0.0135,
     * -0.0112,-0.0091,-0.0071,-0.0052,-0.0035,-0.0019,-0.0004, 0.0009,
     *  0.0022, 0.0033, 0.0044, 0.0053, 0.0062, 0.0070, 0.0077, 0.0083,
     *  0.0088, 0.0093, 0.0097, 0.0101, 0.0104, 0.0106, 0.0108, 0.0110,
     *  0.0111, 0.0112, 0.0112, 0.0112, 0.0112, 0.0111, 0.0111, 0.0109,
     *  0.0108, 0.0107, 0.0105, 0.0103, 0.0101, 0.0099, 0.0097, 0.0095,
     *  0.0093, 0.0090, 0.0088, 0.0085, 0.0083, 0.0080, 0.0077, 0.0075,
     *  0.0072, 0.0070, 0.0067, 0.0064, 0.0062, 0.0059, 0.0057, 0.0054,
     *  0.0052, 0.0050, 0.0047, 0.0045, 0.0043, 0.0041, 0.0038, 0.0036,
     *  0.0034, 0.0032, 0.0031, 0.0029, 0.0027, 0.0025, 0.0024, 0.0022,
     *  0.0021, 0.0019, 0.0018, 0.0016, 0.0015, 0.0014, 0.0013, 0.0012,
     *  0.0010, 0.0009, 0.0008, 0.0008, 0.0007, 0.0006, 0.0005, 0.0004,
     *  0.0004, 0.0003, 0.0002, 0.0002, 0.0001, 0.0001, 0.0000, 0.0000,
     *  0.0000,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0002,
     * -0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,
     * -0.0003,-0.0004,-0.0004,-0.0004,-0.0004,-0.0004,-0.0004,-0.0004,
     * -0.0004,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,
     * -0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,
     * -0.0002,-0.0002,-0.0002,-0.0002,-0.0002,-0.0002,-0.0002,-0.0002,
     * -0.0002,-0.0002,-0.0002,-0.0002,-0.0001,-0.0001,-0.0001,-0.0001,
     * -0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,
     * -0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000 /

      parameter (nbox=800)  ! 100

c>>>>>

c>>>>>

c>>>>>>>>>>>>>>>>>>>

c>>>>>>>>>>>>>>>>>>>
c        If first call, set flexural length scale (crustlen) and
c        coefficient for central deflection (crustcoef).
c        Also interpolate akei from akeic for 2D lithospheric calcs,
c        and set 3D cartesian (x,y,z) coordinates (*cart, saved)
c        for computation of great-circle distances between points.

      if (firstbedrock) then

        if (nx.gt.1 .and. ny.gt.1) then
          crustlen = ( cruststiff / (rhobed*grav) ) ** 0.25
          crustcoef = crustlen**2 / (2.*pi*cruststiff)
        else
          crustlen = ( 4.*cruststiff / (rhobed*grav) ) ** 0.25
          crustcoef = crustlen**3 / (8.*cruststiff)
        endif

        do ix=0,nkei
          zx = ix*dkei
          ia = max (0, int(zx/dkeic))
          ib = min (nkeic, ia+1)
          zwei = (zx - ia*dkeic) / dkeic
          akei(ix) = (1.-zwei)*akeic(ia) + zwei*akeic(ib)
        enddo

c       set local effect array (abox)
c       (8*crustlen gets forebulge...see akeic above)
        zdx = 0.5 * ( dx((nx+1)/2,(ny+1)/2) + dy((nx+1)/2,(ny+1)/2) )
        ibox = nint ((8.*crustlen)/zdx) + 1
c       write(6,*) 'bedrock: ibox = ',ibox
        if (ibox.gt.nbox) then
          write(6,*)
     *      '*** bedrock: ibox too large. ibox=',ibox, '  nbox=',nbox
          stop
        endif
        do ia = -ibox,ibox
          do ja = -ibox,ibox
            zx = sqrt(float(ia**2 + ja**2)) * zdx / crustlen
            ix = min (nint(zx/dkei), nkei)
            abox(ia,ja) = akei(ix)
          enddo
        enddo

        firstbedrock = .false.

      endif
c>>>>>
c>>>>>

      do j=1,ny
        do i=1,nx
          zmass(i,j) =   rhoice*h(i,j)     + rholiq*hw(i,j)
     *                 + rhosed*sedim(i,j) - equiload(i,j)
        enddo
      enddo

c>>>>>>>>>>>>>>>>>>>>>
c>>>>

c        Compute elastic lithospheric and equilibrium asthenospheric
c        displacement (deflect, negative downwards)

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (nx.gt.1 .and. ny.gt.1) then
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c          2D: Brotchie anbd Sylvester

c          Loop over multiple latitude segments, saving deflections
c          due to each in deflectseg, for parallelization.

        call zero (deflectseg, nx*ny*nseg)

c=====================
        do iseg=1,nseg
c=====================
          ja = (iseg-1)*(ny/nseg) + 1
          jb =  iseg   *(ny/nseg)
          if (iseg.eq.nseg) jb = ny
c~~~~~~~~~~~~~~~~~~~~
          do j=ja,jb
            do i=1,nx
c~~~~~~~~~~~~~~~~~~~~

              zload = zmass(i,j) * grav

c---
c             For diagnostic testing of deflect with dbx:
c             zcen = sqrt (   (2.*(j-1.)/(ny-.9999) - 1.)**2
c    *                      + (2.*(i-1.)/(nx-.9999) - 1.)**2 )
c             if (nx.eq.1) zcen = abs (2.*(j-1.)/(ny-1) - 1.)
c             zload = 0.
c             if (zcen.lt.0.5) zload = rhoice*1000.*grav
cc            if (j.gt.ny/2) zload = rhoice*1000.*grav
c--
c                 Add this point-load's deflection to deflect(i,j)

              if (abs(zload).gt.1000.) then    ! > ~10 cm of ice equiv

                zdefcen = zload*darea(i,j)*crustcoef  !nb: times darea

c               polar stereographic:
                jja = max (1,  j-ibox)
                jjb = min (ny, j+ibox)
                iia = max (1,  i-ibox)
                iib = min (nx, i+ibox)
                do jj=jja,jjb
                  do ii=iia,iib
                    deflectseg(ii,jj,iseg) = deflectseg(ii,jj,iseg)
     *                                     + zdefcen*abox(ii-i,jj-j)
                  enddo
                enddo

              endif
c~~~~~~~~~~~~~~~~
            enddo
          enddo
c~~~~~~~~~~~~~~~~
c============
        enddo
c============

c          Sum up deflections from each segment (for parallelization)

        call zero (deflect, nx*ny)
        do iseg=1,nseg
          do j=1,ny
            do i=1,nx
              deflect(i,j) = deflect(i,j) + deflectseg(i,j,iseg)
            enddo
          enddo
        enddo

c%%%%%%%%%
      else
c%%%%%%%%%

c          1D: Turcotte and Schubert:

        call zero (deflect, nx*ny)

        do j=1,ny
          do i=1,nx
            zload = zmass(i,j) * grav

c              Add this point-load's deflection to deflect(i,j)

            if (abs(zload).gt.1000.) then      ! > ~10 cm of ice equiv
              zdefcen = zload*darea(i,j)*crustcoef   ! nb: times darea
              zdx = 0.5 * (dx(i,j) + dy(i,j))
              zdefcen = zload*zdx*crustcoef
              do ii=1,nx
                do jj=1,ny
                  zx =  sqrt (   (xh(i,j)-xh(ii,jj))**2
     *                         + (yh(i,j)-yh(ii,jj))**2 )
     *                  / crustlen
                  deflect(ii,jj) = deflect(ii,jj)
     *                           - zdefcen*exp(-zx)*(cos(zx)+sin(zx))
                enddo
              enddo
            endif

          enddo
        enddo

c%%%%%%%%%%
      endif
c%%%%%%%%%%

c>>>>>
c>>>>>

c        Local bedrock relaxation towards equilibrium

c-------------------------
      if (kasth.eq.0) then
c-------------------------

        do j=1,ny
          do i=1,nx
            za(i,j) = topbed(i,j) - topbedeq(i,j) - deflect(i,j)
            topbed(i,j) = topbed(i,j) - (dt/taulocal)*za(i,j)
          enddo
        enddo

c---------
      else
c---------

c          Diffusive (thin-channel) asthenospheric flow

        do j=1,ny
          do i=1,nx
            za(i,j) = topbed(i,j) - topbedeq(i,j) - deflect(i,j)
          enddo
        enddo

        do j=1,ny
          jm = max(j-1,1)
          jp = min(j+1,ny)
          do i=1,nx
            im = max(i-1,1)
            ip = min(i+1,nx)
              topbed(i,j) = topbed(i,j)
     *                 + dt * asthflow
     *                 * (    ((za(ip,j)-za(i,j))/dxu(i,j) ) * dyu(i,j)
     *                      + ((za(im,j)-za(i,j))/dxu(im,j)) * dyu(im,j)
     *                      + ((za(i,jp)-za(i,j))/dyu(i,j) ) * dxu(i,j)
     *                      + ((za(i,jm)-za(i,j))/dyu(i,jm)) * dxu(i,jm)
     *                   )  / darea(i,j)

          enddo
        enddo

c          Boundary conditions (equilibrium at boundaries)

        if (nx.eq.1) then
          topbed(1,1)  = topbedeq(1,1)  + deflect(1,1)
          topbed(1,ny) = topbedeq(1,ny) + deflect(1,ny)
        else if (ny.eq.1) then
          topbed(1,1)  = topbedeq(1,1)  + deflect(1,1)
          topbed(nx,1) = topbedeq(nx,1) + deflect(nx,1)
        else
          do j=1,ny
            if (j.eq.1 .or. j.eq.ny) then
              iskip = 1
            else
              iskip = nx-1
            endif
            do i=1,nx,iskip
              topbed(i,j) = topbedeq(i,j) + deflect(i,j)
            enddo
          enddo
        endif

c----------
      endif
c----------

c        Reset bed+sed elevation hb

      do j=1,ny
        do i=1,nx
          hb(i,j) = topbed(i,j) + sedim(i,j)
        enddo
      enddo

      return
      end
