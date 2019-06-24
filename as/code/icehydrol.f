! source file: /net/mare/home1/eby/as/ism/icehydrol.F
      subroutine initwater (h, hb, hw, sealev)

c     Sets initial water depths of ocean and inland lakes (hw).
c     No sub-ice lakes or thin film initially (hw=0).
c
c     Generally performs each search by searching one grid box
c     from the current set at each iteration. Maintains list(1:nlist)
c     of indices (i,j, packed) of points added in current iteration,
c     and use only those for the next iteration.

c     Ocean and rivers can propagate N-S, E-W and along diagonals,
c     and lakes too.

      use comicephys
      use comicegrid

      dimension
     *  h(nx,ny), hb(nx,ny), hw(nx,ny)

c     local:
      parameter (nwatermax = 1000)
      dimension
     *  icefloat(nx,ny), maskh2o(nx,ny), sill(nx,ny),indsill(nwatermax),
     *  hsa(nx,ny),                      ! hs as if no floating lake ice
     *  list(nx*ny), listprev(nx*ny)

      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook = 8)   ! 4 for E-W/N-S, 8 for diagonals too

c#######################################################################

c#######################################################################

c        Full calculations, considering ocean access, flotation, lakes

c        Perform search for oceans (submerged points that are in
c        contact with any submerged points at domain edges)

      call izero (icefloat, nx*ny)
      call izero (maskh2o, nx*ny)

c----------------------
      do iter = 1,nx*ny
c----------------------

        nlist = 0

        if (iter.eq.1) then

c         initial ocean set is all submerged pts around domain edges:
          do j=1,ny
            iskip = max(nx-1,1)
            if (j.eq.1.or.j.eq.ny) iskip = 1
            do i=1,nx,iskip
              if (rhoice*h(i,j).le.rholiq*(sealev-hb(i,j))) then
                maskh2o(i,j) = 1
                icefloat(i,j) = 1
                nlist = nlist + 1
                list(nlist) = i + (j-1)*nx
              endif
            enddo
          enddo

        else

c         normal iteration: look at adjacent points for new ocean pts:
          do n=1,nlistprev
            j = (listprev(n)-1) / nx + 1
            i =  listprev(n) - (j-1)*nx
            do look=1,nlook
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (maskh2o(ii,jj).eq.0) then
                if (rhoice*h(ii,jj).le.rholiq*(sealev-hb(ii,jj))) then
                  maskh2o(ii,jj) = 1
                  icefloat(ii,jj) = 1
                  nlist = nlist + 1
                  list(nlist) = ii + (jj-1)*nx
                endif
              endif
            enddo
          enddo

        endif

c          Skip out if no new ocean points found

        if (nlist.eq.0) goto 1000

c          Copy current list to previous list

        nlistprev = nlist
        do n=1,nlist
          listprev(n) = list(n)
        enddo

c----------
      enddo
c----------
 1000 continue

c        Set initial nwater, sill, indsill just for ocean

      nwater = 0
      indsill(1) = 0
      do j=1,ny
        do i=1,nx
          if (maskh2o(i,j).eq.1) then
            sill(i,j) = sealev
            if (nwater.eq.0) then
              nwater = 1                         ! 1 is ocean_embayments
              indsill(1) = i + (j-1)*nx
            endif
          else
c           for non-ocean, just in case skip out to 5000 below
            sill(i,j) = h(i,j) + hb(i,j)
          endif
        enddo
      enddo

c     skip out with just ocean, no inland lakes

      goto 5000 ! 777

c        Perform search for lakes, by eliminating non-lakes (land pts
c        that have a downhill path to a coast)

c        Right now maskh2o=1 (ocean) or 0(not).
c        Initialize non-ocean points (maskh2o=0) to *possible*
c        lake pts (maskh2o=-1), temporarily. Also set  hsa as
c        if no floating lake ice.

      do j=1,ny
        do i=1,nx
          if (maskh2o(i,j).eq.1) then
            hsa(i,j) = sealev + (1.-rhoice/rholiq)*h(i,j)
          else
            maskh2o(i,j) = -1
            hsa(i,j) = hb(i,j) + h(i,j)
          endif
        enddo
      enddo

c======================
      do iter = 1,nx*ny
c======================

        nlist = 0

        if (iter.eq.1) then

c         Set initial non-lake set (maskh2o=0) to all coastal points,
c         and all land pts at edges of domain (we don't know about
c         slopes outside):
          do j=1,ny
            do i=1,nx
              if (maskh2o(i,j).eq.-1) then
                iflag = 0
                if (i.eq.1.or.i.eq.nx.or.j.eq.1.or.j.eq.ny) then
                  iflag = 1
                else
                  do look=1,nlook
                    ii = max (1, min (nx, i + ish(look)))
                    jj = max (1, min (ny, j + jsh(look)))
                    if (maskh2o(ii,jj).eq.1) iflag = 1
                  enddo
                endif
                if (iflag.eq.1) then
                  maskh2o(i,j) = 0
                  nlist = nlist + 1
                  list(nlist) = i + (j-1)*nx
                endif
              endif
            enddo
          enddo

        else

c         normal iteration: look uphill for new non-lake points:
          do n=1,nlistprev
            j = (listprev(n)-1) / nx + 1
            i =  listprev(n) - (j-1)*nx
            do look=1,nlook
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if ( hsa(ii,jj).ge.hsa(i,j) .and.
     *             maskh2o(ii,jj).eq.-1 ) then
                maskh2o(ii,jj) = 0
                nlist = nlist + 1
                list(nlist) = ii + (jj-1)*nx
              endif
            enddo
          enddo

        endif

c          Skip out if no new non-lake points found

        if (nlist.eq.0) goto 2000

c          Copy current list to previous list

        nlistprev = nlist
        do n=1,nlist
          listprev(n) = list(n)
        enddo

c==========
      enddo
c==========
 2000 continue

c       Propagate sill-depth constraints through lakes.
c       sill already set above for ocean (maskh2o=1).

      do j=1,ny
        do i=1,nx
          if (maskh2o(i,j).eq.-1) then
            sill(i,j) = 1.e20
          else if (maskh2o(i,j).eq.0) then
            sill(i,j) = 0.
          endif
        enddo
      enddo

      do iter = 1,nx*ny
        nchange = 0
        do j=1,ny
          do i=1,nx
            if (maskh2o(i,j).eq.-1) then
              sillold = sill(i,j)
              do look = 1,nlook
                ii = max (1, min (nx, i + ish(look)))
                jj = max (1, min (ny, j + jsh(look)))
                if (maskh2o(ii,jj).eq.0) then
c                 land pt (cannot be ocean,all coastal pts are non-lake:
                  sill(i,j) = min (sill(i,j), hsa(ii,jj))
                else
c                 lake point:
                  sill(i,j) = min (sill(i,j), sill(ii,jj))
                endif
              enddo
c             remove "emergent land" from lake status:
              if (sill(i,j).le.hsa(i,j)) then
                maskh2o(i,j) = 0
                icefloat(i,j) = 0
                sill(i,j) = 0.
              endif
              if (sillold.ne.sill(i,j)) then
                nchange = nchange + 1
              endif
            endif
          enddo
        enddo
        if (nchange.eq.0) goto 3000
        if (iter.gt.nx*ny-10) then
          write (ioterm,*) '*** Warning (initwater): iter=',iter
        endif
      enddo
      write (ioterm,*)'*** Error (initwater): exceeded iteration: ',iter
      stop
 3000 continue

c        Convert -1 maskh2o pts (currently 0=land, 1=ocean, -1=inland
c        lakes) to labels (2,3,..) that identify contiguous inland lakes

      do j=1,ny
        do i=1,nx

          if (maskh2o(i,j).eq.-1) then
            nwater = nwater + 1
            if (nwater.gt.nwatermax) then
              write(ioterm,*)'*** Error (initwater): too many lakes'
              stop
            endif
            maskh2o(i,j) = nwater
            indsill(nwater) = 0

            nlistprev = 1
            listprev(1) = i + (j-1)*nx

            do iter = 1,nx*ny

              nlist = 0
              do n=1,nlistprev
                ja = (listprev(n)-1) / nx + 1
                ia =  listprev(n) - (ja-1)*nx
                do look = 1,nlook
                  ii = max (1, min (nx, ia + ish(look)))
                  jj = max (1, min (ny, ja + jsh(look)))
                  if ( maskh2o(ii,jj).gt.0) then
                    if (maskh2o(ii,jj).ne.maskh2o(i,j)) then
                      write (ioterm,'(a,6i4)')
     *                  '*** Error 1 (initwater) ii,jj,i,j,maskh2o=',
     *                  ii, jj, i, j, maskh2o(ii,jj), maskh2o(i,j)
                      stop
                    endif
                  else if (maskh2o(ii,jj).eq.-1) then
                    if (sill(ii,jj).ne.sill(i,j)) then
                      write (ioterm,'(a,6i4/a,2f22.15)')
     *                  '*** Error 2 (initwater) ii,jj,i,j,maskh2o=',
     *                  ii, jj, i, j, maskh2o(ii,jj), maskh2o(i,j),
     *                  '    sills=',sill(ii,jj), sill(i,j)
                      stop
                    endif
                    maskh2o(ii,jj) = maskh2o(i,j)
                    nlist = nlist + 1
                    list(nlist) = ii + (jj-1)*nx

                  else if (maskh2o(ii,jj).eq.0) then
c                   set location of first-found lake sill point (land)
                    if (hsa(ii,jj).eq.sill(i,j) .and.
     *                  indsill(nwater).eq.0) then
                      indsill(nwater) = ii + (jj-1)*nx
                    endif

                  endif
                enddo
              enddo

c                Skip out if no new lake points found

              if (nlist.eq.0) goto 4000

c                Copy current list to previous list

              nlistprev = nlist
              do n=1,nlist
                listprev(n) = list(n)
              enddo

            enddo
 4000       continue

          endif

        enddo
      enddo

c        At this point, inland lakes rest on top of either land or ice.
c        [The only floating ice is in ocean embayments, bordered by
c        coastal land with no ice or non-floating ice.]
c        Now look to see if inland lakes float adjacent (non-ocean) ice.
c        This expands inland lakes and creates floating ice on lakes.
c        Assume that lake hydrostatic pressure penetrates horizontally
c        only along lake-land interface; for lakes with ice bottoms,
c        cannot penetrate solid ice (ie, not into adjacent cells
c        with ice...assume continuous ice even if the adjacent ice
c        base is higher than the current cell's ice top).

c~~~~~~~~~~~~~~~~~~~~~~
      do iter = 1,nx*ny
c~~~~~~~~~~~~~~~~~~~~~~

        nfloata = 0
        nfloatb = 0

        do j=1,ny
          do i=1,nx
            if ( maskh2o(i,j).ne.1 .and.
     *           h(i,j).gt.0. .and. icefloat(i,j).eq.0 ) then

              do look=1,nlook
                ii = max (1, min (nx, i + ish(look)))
                jj = max (1, min (ny, j + jsh(look)))

c               Two ways to float ice:

c           (a) local lake over ice, transmitted to base by adjacent
c               ice-free land (which must be higher than local sill):
                if (maskh2o(i,j).ge.2 .and.
     *              maskh2o(ii,jj).eq.0 .and. h(ii,jj).eq.0.) then
                  if (rhoice*h(i,j).le.rholiq*(sill(i,j)-hb(i,j))) then
                    icefloat(i,j) = 1
                    nfloata = nfloata + 1
                  endif

c           (b) land ice, adjacent to lake with land bottom
                else if (maskh2o(i,j).eq.0 .and.
     *                   maskh2o(ii,jj).ge.2 .and. h(ii,jj).eq.0.)then
                  if (rhoice*h(i,j).le.rholiq*(sill(ii,jj)-hb(i,j)))then
                    sill(i,j) = sill(ii,jj)
                    maskh2o(i,j) = maskh2o(ii,jj)
                    icefloat(i,j) = 1
                    nfloatb = nfloatb + 1
                  endif
                endif
              enddo

            endif
          enddo
        enddo

c          Skip out if no new floating-ice points found

        if (nfloata + nfloatb.eq.0) goto 5000

c~~~~~~~~~~
      enddo
c~~~~~~~~~~

 5000 continue

c       Set water thickness hw

      do j=1,ny
        do i=1,nx
          if (maskh2o(i,j).eq.0) then
            hw(i,j) = 0.
          else
            hw(i,j) = sill(i,j) - (rhoice/rholiq)*h(i,j) - hb(i,j)
            if (hw(i,j).lt.0.) then
              write (ioterm,*) '*** Error (initwater): -ve hw, i,j=',i,j
              stop
            endif
          endif
        enddo
      enddo

      return

c####
c#####

      end

c-----------------------------------------------------------------------

      subroutine findwater (maskwater,
     *                      indlake, npoilake, nlake,
     *                      h, hb, hw, sealev, timeice)

c        Locates oceans from scratch, resets maskwater and hw
c        for ocean and adjacent land. Resets lakes
c        (maskwater,indlake,npoilake,nlake)
c        depending on existing hw.

c     maskwater = 0=grounded ice or ice-free land,
c                 1=ocean,
c                 +/- 2,3,4,...=lake number (+ open, - subice) (ret)
c                 Here, lakes all +ve, set to + or - in adjustpres
c     indlake   = (packed) i,j of each lake's points (returned)
c     npoilake  = number of points in each lake      (returned)
c     nlake     = # of lakes (including ocean)        (returned)
c     h         = ice thickness                       (supplied)
c     hb        = bed+sed elevation                   (supplied)
c     hw        = water thickness           (modified for ocean)
c     sealev    = sea level                           (supplied)
c
      use comicephys
      use comicegrid

      dimension
     *  maskwater(nx,ny), indlake(npoimax,nlakemax),
     *  npoilake(nlakemax),
     *  h(nx,ny), hb(nx,ny), hw(nx,ny)

c     local:
      dimension
     *  list(nx*ny),      listprev(nx*ny)

      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook = 8)   ! 4 for E-W/N-S, 8 for diagonals too

      call izero (maskwater, nx*ny)

c        Perform search for oceans (submerged points that are in
c        contact with any submerged points at domain edges),
c        setting maskwater to 1 and setting hw for ocean pts,
c        and also setting hw to 0 for non-floating ocean-adjacent pts.

c----------------------
      do iter = 1,nx*ny
c----------------------

        nlist = 0

        if (iter.eq.1) then

          do j=1,ny
            iskip = nx-1
            if (j.eq.1.or.j.eq.ny) iskip = 1
            do i=1,nx,iskip
              if (sealev-hb(i,j)-rhor*h(i,j).gt.hwcut) then
                maskwater(i,j) = 1
                nlist = nlist + 1
                list(nlist) = i + (j-1)*nx
                hw(i,j) = sealev - hb(i,j) - rhor*h(i,j)
              else
                hw(i,j) = 0.
              endif
            enddo
          enddo

        else

c         normal iteration, look at adjacent points for new ocean pts:
          do n=1,nlistprev
            j = (listprev(n)-1) / nx + 1
            i =  listprev(n) - (j-1)*nx
            do look=1,nlook
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (maskwater(ii,jj).eq.0) then
                if (sealev-hb(ii,jj)-rhor*h(ii,jj).gt.hwcut) then
                  maskwater(ii,jj) = 1
                  nlist = nlist + 1
                  list(nlist) = ii + (jj-1)*nx
                  hw(ii,jj) = sealev - hb(ii,jj) - rhor*h(ii,jj)
                else
                  hw(ii,jj) = 0.
                endif
              endif
            enddo
          enddo

        endif

c          Skip out if no new ocean points found

        if (nlist.eq.0) goto 1000

c          Copy current list to previous list

        nlistprev = nlist
        do n=1,nlist
          listprev(n) = list(n)
        enddo

c----------
      enddo
c----------
 1000 continue

c        Search for inland lakes, based on water thickness hw vs hwcut.
c        When find a new, not-already-found lake point (outer loop),
c        increment nlake, then iteratively search for adjacent lake
c        points, forming list indlake and setting maskwater,
c        and set npoilake when done.

      nlake = 1 ! ocean
      call izero (indlake,  npoimax*nlakemax)
      call izero (npoilake, nlakemax)

c     don't allow any non-ocean water (open lakes or sub-ice)
      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).ne.1) then
c            convert liquid water beneath ice to ice (or not..why?)
c            if (h(i,j).gt.hwcut) h(i,j) = h(i,j) + hw(i,j)/rhor
             hw(i,j) = 0.
          endif
        enddo
      enddo
      return

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do ja=1,ny
        do ia=1,nx
          if (maskwater(ia,ja).eq.0 .and. hw(ia,ja).gt.hwcut) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            nlake = nlake + 1
            if (nlake.gt.nlakemax) then
               write (ioterm, *)
     *           'Error (findwater): max # lakes exceeded, nlake=',nlake
               stop
             endif

c----------------------------
            do iter = 1,nx*ny
c----------------------------

              nlist = 0

              if (iter.eq.1) then

                nlist = 1
                list(1) = ia + (ja-1)*nx
                maskwater(ia,ja) = nlake
                npoi = 1
                indlake(npoi,nlake) = list(1)

              else

c               normal iteration, look at adjacent pts for new lake pts
                do n=1,nlistprev
                  j = (listprev(n)-1) / nx + 1
                  i =  listprev(n) - (j-1)*nx
                  do look=1,nlook
                    ii = max (1, min (nx, i + ish(look)))
                    jj = max (1, min (ny, j + jsh(look)))
                    if (maskwater(ii,jj).eq.0 .and. hw(ii,jj).gt.hwcut)
     *                then
                      maskwater(ii,jj) = nlake
                      nlist = nlist + 1
                      list(nlist) = ii + (jj-1)*nx
                      npoi = npoi + 1
                      if (npoi.gt.npoimax) then
                        write (ioterm,'(/a/a,i6,a,i6)')
     *                    'Error (findwater): max lake size exceeded:',
     *                    '   npoi=',npoi,'  nlake=',nlake
                        stop
                      endif
                      indlake(npoi,nlake) = list(nlist)
                    else if (maskwater(ii,jj).ne.nlake) then
                      write (ioterm,'(/a/a,i6/a,i6/a,2i6/a,2i6)')
     *                  'Error (findwater): lake iteration error:',
     *                  '    maskwater=',maskwater(ii,jj),
     *                  '    nlake    =',nlake,
     *                  '    ii,jj    =',ii,jj,
     *                  '    ia,ja    =',ia,ja
                      stop
                    endif
                  enddo
                enddo

              endif

c                Skip out if no new open-lake points found

              if (nlist.eq.0) goto 2000

c                Copy current list to previous list

              nlistprev = nlist
              do n=1,nlist
                listprev(n) = list(n)
              enddo

c----------------
            enddo
c----------------
 2000       continue
            npoilake(nlake) = npoi

c>>>>>>>>>>>>>>
          endif
        enddo
      enddo
c>>>>>>>>>>>>>>

c        Finally, set all maskwater's within each lake to negative if
c        lake is closed (else maskwater remains positive, open lake)

      do ilake = 2,nlake

        ifopen = 0
        do n=1,npoilake(ilake)
          j = (indlake(n,ilake)-1) / nx + 1
          i =  indlake(n,ilake) - (j-1)*nx
          if (h(i,j).eq.0.) ifopen = 1
        enddo

        if (ifopen.eq.0) then
          do n=1,npoilake(ilake)
            j = (indlake(n,ilake)-1) / nx + 1
            i =  indlake(n,ilake) - (j-1)*nx
            maskwater(i,j) = -maskwater(i,j)
          enddo
        endif

      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine adjustpres (maskwater, indlake, npoilake, nlake,
     *                       h, hb, hw, hs, sealev)

c     Adjusts water thicknesses hw within each separate lake or ocean
c     to equalize hydrostatic water pressure on a commom reference
c     level. Then resets surface elevation hs.
c     Ignore consequent changes to maskwater, indlake, etc (i.e., if
c     change in hw would make ice ground or float, or if hw would
c     overflow sill of an open lake...wait till next call to findwater
c     and/or movewater).
c
c     maskwater = 0,1,+/-2... (grounded ice/ice-free land, ocean, lake#)
c     indlake   = (packed) i,j of each lake's points
c     npoilake  = number of points in each lake
c     nlake     = # of lakes (including ocean)
c     h         = ice thickness
c     hb        = bed+sed elevation
c     hw        = water thickness (modified)
c     hs        = surface elevation (modified)
c     sealev    = sea level

      use comicephys
      use comicegrid

      dimension
     *  maskwater(nx,ny), indlake(npoimax,nlakemax), npoilake(nlakemax),
     *  h(nx,ny), hb(nx,ny), hw(nx,ny), hs(nx,ny)

      dimension pres(npoimax)

c        Set hw and hs for land (no lake) and ocean

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.0) then
            hs(i,j) = hb(i,j) + hw(i,j) + h(i,j)
          else if (maskwater(i,j).eq.1) then
            hw(i,j) = max (0., sealev - hb(i,j) - rhor*h(i,j))
            hs(i,j) = hb(i,j) + hw(i,j) + h(i,j)
          endif
        enddo
      enddo

c        Loop for each lake (open or sub-ice, but not ocean)

c-----------------------
      do ilake = 2,nlake
c-----------------------

c          Calculate pres = pressure/(rholiq*grav) at "zero" elevation
c          for each lake point, and its areal average totp

        totp = 0.
        tota = 0.
        do n=1,npoilake(ilake)
          j = (indlake(n,ilake)-1) / nx + 1
          i =  indlake(n,ilake) - (j-1)*nx
          pres(n) = rhor*h(i,j) + (hb(i,j)+hw(i,j))
          totp = totp + pres(n)*darea(i,j)
          tota = tota +         darea(i,j)
        enddo
        totp = totp/tota

c          Subtract pres-totp from each point's water thickness,
c          to equalize new pres's for all points. This conserves
c          total water in lake, unless...if new hw would be negative
c          (newly grounded ice), just set hw to zero, ignore water
c          non-conservation, and wait till next call to findwater
c          to correct lake indices.
c          This workd both for closed lakes (ice rides above
c          hw, and open lakes (ice floats in hw).

        do n=1,npoilake(ilake)
          j = (indlake(n,ilake)-1) / nx + 1
          i =  indlake(n,ilake) - (j-1)*nx
          hw(i,j) = max (0., hw(i,j) - (pres(n)-totp))
          hs(i,j) = hb(i,j) + hw(i,j) + h(i,j)
        enddo

c----------
      enddo
c----------

      return
      end

c-----------------------------------------------------------------------

      subroutine basecoef (h, hb, hw, t, tsurfi, sedpres, sealev,
     *                     crhu, crhv, powbu, powbv,
     *                     fracgu, fracgv, hgu, hgv)

c     Sets basal flow terms crhu, crhv (u,v grids). These depend on
c     basal conditions mostly on h grid (eg, floating or not), so
c     calculations here deal with interpolation to u or v grid.
c     This will be expanded when add explicit basal hydrology
c     (probably on h-grid), and effect on basal flow.

      use comicephys
      use comicegrid

      dimension
     *  h(nx,ny),             hb(nx,ny),        hw(nx,ny),
     *  t(nx,ny,0:nlevp),     tsurfi(nx,ny),    sedpres(nx,ny),
     *  crhu(0:nxp,0:nyp),    crhv(0:nxp,0:nyp),
     *  powbu(0:nxp,0:nyp),   powbv(0:nxp,0:nyp),
     *  fracgu(0:nxp,0:nyp),  fracgv(0:nxp,0:nyp),
     *  hgu(0:nxp,0:nyp),     hgv(0:nxp,0:nyp)

c for powb=2: zcrh=1.e-20 ~no sliding, 1.e-12 ~some, 1.e-8 ~stream.

c        Set basal sliding coefficients crhu (u grid, iloopuv=1), and
c        crhv (v grid, iloopuv=2), looping over u and v grids to avoid
c        duplicate code. Do for all points regardless of  water...
c        values will be re-set for water within icedyn iterations.
c        Also crh[u,v] will be adjusted in sedflow for deforming sed.

      zpowb = powb  ! powb is a (possibly int) parameter in comicephys
      call resetr (powbu, (nxp+1)*(nyp+1), zpowb)
      call resetr (powbv, (nxp+1)*(nyp+1), zpowb)

      call resetr (crhu, (nxp+1)*(nyp+1), crheolbno)
      call resetr (crhv, (nxp+1)*(nyp+1), crheolbno)

      call resetr (fracgu, (nxp+1)*(nyp+1), 1.)
      call resetr (fracgv, (nxp+1)*(nyp+1), 1.)

c>>>>>>>>>>>>>>>>>>>>>
      do iloopuv = 1,2
c>>>>>>>>>>>>>>>>>>>>>

      if (iloopuv.eq.1) then
        ish = 1
        jsh = 0
        nxtmp = nx-1
        nytmp = ny
      else
        ish = 0
        jsh = 1
        nxtmp = nx
        nytmp = ny-1
      endif

c        Set basic sliding value for various experiments

      do j=1,nytmp
        do i=1,nxtmp
          ia = i + ish
          ja = j + jsh

c.....................................................
c......................

          zcrh = zcrh0 !Jer: ice-sheet-specific zcrh0 read in from namelist file

          zlat = 0.5*(alatd(i,j) + alatd(ia,ja))                 ! 666
          zlon = 0.5*(alond(i,j) + alond(ia,ja))                 ! 666
          if (abs(alond(i,j)-alond(ia,ja)).gt.30.) zlon = 180.   ! 666

c         slippery bed under WAIS ice stream region:             ! 666
c         if ( zlat.gt. -86. .and. zlat.lt. -78. .and.           ! 666
c    *         zlon.gt.-170. .and. zlon.lt.-125. )               ! 666
c    *      zcrh = 1.e-8                                         ! 666

c         slippery bed where modern rebounded topog is below sl:
          if (sedpres(i,j).gt.0. .and. sedpres(ia,ja).gt.0.)     ! 666
c    *      zcrh = 1.e-10                                        ! 666
c    *      zcrh = 1.e-8                                         ! 666
c    *      zcrh = 1.e-7                                         ! 666
     *      zcrh = 1.e-6                                         ! 666

c         intermediate value over Pine Island region:
          if ( zlat.gt. -80. .and.                               ! 666
     *         zlon.gt.-120. .and. zlon.lt.-90. )                ! 666
c    *      zcrh = 1.e-8                                         ! 666
     *      zcrh = 1.e-10                                        ! 666

c         bedrock value over most of EAIS including Transantarctics:
          if ( zlon.lt. -170. .or.  zlon.gt.0.)                  ! 666
     *      zcrh = min (3.e-7, zcrh)                             ! 666
c    *      zcrh = min (1.e-7, zcrh)                             ! 666
c    *      zcrh = 1.e-8                                         ! 666
c    *      zcrh = 1.e-10                                        ! 666

c....................................
c.....

c              If at grounding line, calculate areal fraction of this
c              u or v-grid cell that's grounded (fracg[u,v]). i.e.,
c              estimated sub-grid position of grounding line. Calculate
c              height-above-flotation for the two surrounding h-grid
c              points (one floating with h.a.f. < 0, one grounded with
c              h.a.f. > 0), and linearly interpolated zero-location
c              is the estimated g.l position. [This is equivalent to
c              linearly interpolating ice surface, ice base, and bed
c              elevations between the 2 h-points, and seeing where
c              flotation occurs]. fracg[u,v] is used in icedyn,
c              multiplying the basal coefficient coefbs[u,v].

c              Also calc g.l. ice thickness hg[u,v] (using fracg[u,v]).
c              Currently fracg[u,v] (and hg[u,v]) are diagnostic only
c              if schoofgl not defined. For schoofgl, fracg[u,v] are
c              used in calcgl to calculate hg[u,v] for schoof g.l. flux
c              (repeating same calc as here, but within icedyn C-loop).

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if (     (hw(i,j).gt.hwcut .and. hw(ia,ja).le.hwcut)
     *        .or. (hw(i,j).le.hwcut .and. hw(ia,ja).gt.hwcut) ) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            zpot  = hb(i,j)   + hw(i,j)   + rhor*h(i,j)
            zpota = hb(ia,ja) + hw(ia,ja) + rhor*h(ia,ja)
            if (hw(i,j).gt.hwcut) then
              zaf = zpota - zpot     ! ht above flot (liq m), grounded
              zbf = -hw(i,j)         ! ht above flot (liq m), floating
            else if (hw(ia,ja).gt.hwcut) then
              zaf = zpot - zpota
              zbf = -hw(ia,ja)
            endif
            zfgl = zaf / max (zaf-zbf, 0.1)    ! for safety...should
            zfgl = max (0., min (1., zfgl))  ! have zaf>0, zbf<0
c           zfgl = 0.5 ! 777

c           or as in Pattyn et al (2006) (doesn't work if h=0):
c           if (hw(i,j).gt.hwcut) then
c             zf = (sealev - hb(i,j))   / (rhor*max (h(i,j),.1))
c             zg = (sealev - hb(ia,ja)) / (rhor*max (h(ia,ja),.1))
c           else
c             zg = (sealev - hb(i,j))   / (rhor*max (h(i,j),.1))
c             zf = (sealev - hb(ia,ja)) / (rhor*max (h(ia,ja),.1))
c           endif
c           zfgl = (1.-zg) / max (zf-zg, .01)

c           calc g.l. depth (lin.interp bed elev to u,v grid, vs. s.l.,
c           so only works for marine g.l.s, not inland lakes):
c           zhbg =       0.5*hb(i,j) +  0.5*hb(ia,ja)
c           slightly better in 100/300ka tests(?):
            if (hw(i,j).gt.hwcut) then
              zhbg = (1.-zfgl)*hb(ia,ja) + zfgl*hb(i,j)
            else
              zhbg = (1.-zfgl)*hb(i,j)   + zfgl*hb(ia,ja)
            endif
            zhg  = (sealev-zhbg)/rhor
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else if (hw(i,j).gt.hwcut .and. hw(ia,ja).gt.hwcut) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            zfgl = 0.
            zhg  = 0.
c~~~~~~~~~~~~~
          else
c~~~~~~~~~~~~~
            zfgl = 1.
            zhg  = 0.
c~~~~~~~~~~~~~~
          endif
c~~~~~~~~~~~~~~

          if (iloopuv.eq.1) then
            fracgu(i,j) = zfgl
            hgu(i,j)    = zhg
          else
            fracgv(i,j) = zfgl
            hgv(i,j)    = zhg
          endif

c++++++++++++++++++++++
c++++++++++++++++++++++

c         no basal sliding if base is below pressure melting point:

c#####
c old:
c#####
c         zt = 0.
c         nzt = 0
c         if (h(i,j).gt.0. .and. hw(i,j).le.hwcut) then
c           zt = zt + t(i,j,nlevp)-(tmelt-dtmdh*h(i,j))
c           nzt = nzt + 1
c         endif
c         if (h(ia,ja).gt.0. .and. hw(ia,ja).le.hwcut) then
c           zt = zt + t(ia,ja,nlevp)-(tmelt-dtmdh*h(ia,ja))
c           nzt = nzt + 1
c         endif
c         if (nzt.ge.1) then
c           zt = zt/nzt
c         else
c           if (hw(i,j).le.hwcut .and. hw(ia,ja).le.hwcut) then
c             zt = 0.5*(tsurfi(i,j) + tsurfi(ia,ja)) - tmelt
c           else
c             zt = 0.
c           endif
c         endif
c#####
c new:
c#####
c         more influence of adjacent floating base:

          zt = 0.
          nzt = 0
          if (h(i,j).gt.0.) then
            zt = zt + t(i,j,nlevp)-(tmelt-dtmdh*h(i,j))
            nzt = nzt + 1
          else if (hw(i,j).gt.hwcut) then
            zt = zt + 0.
            nzt = nzt + 1
          endif

          if (h(ia,ja).gt.0.) then
            zt = zt + t(ia,ja,nlevp)-(tmelt-dtmdh*h(ia,ja))
            nzt = nzt + 1
          else if (hw(ia,ja).gt.hwcut) then
            zt = zt + 0.
            nzt = nzt + 1
          endif

          if (nzt.ge.1) then
            zt = zt/nzt
          else
            zt = 0.5*(tsurfi(i,j) + tsurfi(ia,ja)) - tmelt
          endif
c#####

          zm = min (1., max (0., (zt + tramp)/tramp))

c         influence of fracg[u,v]:
          zm = zm*zfgl + 1.*(1.-zfgl)

c         no freezing over sediment areas:
c         if (sedpres(i,j).gt.0. .and. sedpres(ia,ja).gt.0.)
c    *     zm = 1

c         no freezing for WAIS Ross ice stream region:             ! 666
c         zlat = 0.5*(alatd(i,j) + alatd(ia,ja))                   ! 666
c         zlon = 0.5*(alond(i,j) + alond(ia,ja))                   ! 666
c         if ( zlat.gt. -86. .and. zlat.lt. -78. .and.             ! 666
c    *         zlon.gt.-170. .and. zlon.lt.-125. ) zm = 1.         ! 666

c         change form of zm vs crh depeendence:
c         zm = zm ** 0.2
c         zcrh = 10.** ((1.-zm)*log10(crheolbno) + zm*log10(zcrh))
c         zcrh =        (1.-zm)*     (crheolbno) + zm*     (zcrh)

          zcrh = 10.** ((1.-zm)*log10(crheolbno) + zm*log10(zcrh))
c+++++
c+++++

          if (iloopuv.eq.1) then
            crhu(i,j) = zcrh
          else
            crhv(i,j) = zcrh
          endif

        enddo  ! i=1,nxtmp
      enddo    ! j=1,nytmp

c>>>>>>>>>>
      enddo   !iloopuv loop
c>>>>>>>>>>

      return
      end

c-----------------------------------------------------------------------

c**********************
c*****
