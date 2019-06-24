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
#if defined EISLINE || defined EISANTA || defined NHA || defined CARB
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
#else
c####

c       Simple setting of hw (ignoring access to ocean issues)  

      do j=1,ny
        do i=1,nx
          hw(i,j) = max (0., sealev - (rhoice/rholiq)*h(i,j) - hb(i,j))
        enddo
      enddo

      return

c#####
#endif
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

#if defined EISLINE && (defined LINEB || defined LINEE || defined TRANSECTA)

c     simple setting of hw (flotation by ocean at sealev everywhere,
c     ignoring access issues):
      call izero (maskwater, nx*ny)
      call zero (hw, nx*ny)
      do j=1,ny
        do i=nx,1,-1
          if (rhoice*h(i,j).le.rholiq*(sealev-hb(i,j))) then
            maskwater(i,j) = 1
            hw(i,j) = sealev - hb(i,j) - (rhoice/rholiq)*h(i,j)
          endif
        enddo
      enddo
      nlake = 1 ! ocean
      call izero (indlake,    npoimax*nlakemax)
      call izero (npoilake,   nlakemax)
      return

#elif defined EISLINE && defined LINEF

c     simple setting of ocean (scanning r-to-l to first grounded pt):
      call izero (maskwater, nx*ny)
      call zero (hw, nx*ny)
      if (nx.gt.ny) then
        do j=1,ny
          do i=nx,1,-1
            if (rhoice*h(i,j).le.rholiq*(sealev-hb(i,j))) then
              maskwater(i,j) = 1
              hw(i,j) = sealev - hb(i,j) - (rhoice/rholiq)*h(i,j)
            else
              goto 50
            endif
          enddo
   50     continue
        enddo
      else
        do i=1,nx
          do j=ny,1,-1
            if (rhoice*h(i,j).le.rholiq*(sealev-hb(i,j))) then
              maskwater(i,j) = 1
              hw(i,j) = sealev - hb(i,j) - (rhoice/rholiq)*h(i,j)
            else
              goto 52
            endif
          enddo
   52     continue
        enddo
      endif
      nlake = 1 ! ocean
      call izero (indlake,    npoimax*nlakemax)
      call izero (npoilake,   nlakemax)
      return

#endif

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

#if defined NOWATER 
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
#endif

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
#if defined EISA || defined EISB || defined EISC       
c.....................................................
          zcrh = crheolbno
c...................
#elif defined EISD
c...................
          zcrh = 1.0e-3    ! (m/a) Pa-1 
c......................
#elif defined EISLINE
c......................

#if defined LINEB 
          zr = xh(i,j)
          if (zr.lt.300.e3) then
c           zcrh = crheolbno
            zcrh = 1.e-8
c           zcrh = 1.e-6
c           zcrh = 1.e-4
          else
c           zcrh = 1.e-7
c           zcrh = 1.e-6
c           zcrh = 1.e-6
c           zcrh = 3.e-6
c           zcrh = 1.e-5
c           zcrh = 3.e-5
            zcrh = 1.e-4
c           if (hb(i,j).gt.-200. .and. zr.gt.700.e3) zcrh = 1.e6 ! inf
c           zcrh = 3.e-4
          endif
#elif defined LINEC
          zcrh = 1.e-10
#elif defined LINED
          zcrh = 1.e-10
#elif defined LINEE
c         zcrh = 1.e-8
          zcrh = 2.e-5    ! abc
c         if (hb(i,j).gt.hb(1,j)) then
c           zcrh = 0.5*zcrh                      ! wedge 2x as sticky
c         endif
#elif defined LINEF
          if (nx.gt.ny) then
            zr = xh(i,j)
          else
            zr = yh(i,j)
          endif
c         if (zr.lt.100.e3) then
c           zcrh = 3.e-9
cc          zcrh = 3.e-10
c         else
c           zcrh = 3.e-9
c         endif
c         zcrh = 1.e-10
c         zcrh = 1.e-8
c         zcrh = 3.e-8
          zcrh = 3.e-9  ! std
c         zcrh = 1.e-9  ! Fa
c         zcrh = 1.e-8  ! Fb (and Fb_oldaslip wit budgsnow=0.1)
c         zcrh = 3.e-8  ! Fc
c         zcrh = (86400.*365.) / (7.624e6**powbu(i,j))  ! powb=3?,schoof
#elif defined LINEG
          zcrh = 1.e-10
#elif defined TRANSECTA
c         zcrh = 3.e-9
c         zcrh = 1.e-8
c         zcrh = 2.e-8 
c         zcrh = 3.e-8 
c         zcrh = 6.e-8 
          zcrh = 1.e-7 
c         if (zr.gt.1900.e3 .and. zr.lt.2200.e3) then  ! weddell+ross
c         if (zr.lt.200.e3) then                       ! ross only
cc          zcrh = 2.e-7 ! too big ice
c           zcrh = 3.e-7 ! u
cc          zcrh = 6.e-7 ! s,t
cc          zcrh = 1.e-6 ! r
c         else
c           zcrh = 1.e-5 ! s,t,u
c         endif
#endif

c....................
#elif defined HEINO
c....................

          if ( (maskheino(i,j).eq.2 .or. maskheino(ia,ja).eq.2) 
     *       ) then
c           soft sediment:
#if defined HEINO_S1
            zcrh = 100. / (rhoice*grav)
#elif defined HEINO_S2
            zcrh = 200. / (rhoice*grav)
#elif defined HEINO_S3
            zcrh = 1000. / (rhoice*grav)
#else
            zcrh = 500. / (rhoice*grav)
#endif
cx          zload = rhoice * grav * max (0.5*(h(i,j)+h(ia,ja)),10.)!cx
cx          zcrh = 1.e8 / (rhoice*grav*(zload**2))                 !cx
            if (iloopuv.eq.1) then
              powbu(i,j) = 1. ! cx 3
            else
              powbv(i,j) = 1. ! cx 3
            endif
          else 
c           hard rock (used for ocean too):
            zload = rhoice * grav * max (0.5*(h(i,j)+h(ia,ja)), 10.)
            zcrh = 1.e5 / (rhoice*grav*(zload**2))
            if (iloopuv.eq.1) then
              powbu(i,j) = 3.
            else
              powbv(i,j) = 3.
            endif
          endif

c.......................................................................
#elif defined HOMA || defined HOMB || defined HOMC||defined HOMD
c.......................................................................

          zomega = 2.*pi/(xh(nx-1,1)-xh(1,1))
          if (iloopuv.eq.1) then
            zx = (i-0.5)*dx0
            zy = (j-1.0)*dy0
          else
            zx = (i-1.0)*dx0
            zy = (j-0.5)*dy0
          endif
#if defined HOMA || defined HOMB
          zcrh = crheolbno    ! nb: powb=1
#elif defined HOMC
          zcrh = 1. / (1000. + 1000.*sin(zomega*zx)*sin(zomega*zy))
#elif defined HOMD
          zcrh = 1. / (1000. + 1000.*sin(zomega*zx))
#endif

c...................
#elif defined HOME
c...................

          zcrh = crheolbno    ! nb: powb=1

c...................
#elif defined HOMF
c...................

c         nb: powb=1:
c         cslipf = 0.                     ! chom f noslip
          cslipf = 1.                     ! chom f slip
c         blows with cslipf=1:
c         zcrh = (1.-cslipf) * crheolbno + cslipf * 2.140373e-7*1000.
c         ok:
          zcrh = (1.-cslipf) * crheolbno + cslipf * 1.000000e-7*1000.

c......................
#elif defined EISANTA
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
#elif defined NHA || defined CARB
c....................................
          zcrh = 1.0e-10

c....
#else
c....
          zcrh = 1.0e-12 
          zcrh = 1.0e-8
c.....
#endif
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
#if ! defined NOBASET 
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
#endif 
c+++++

# if defined HOMESLIP
          if (isliparol(i).eq.1) zcrh = 0.   ! free slip, in common
#endif
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
#if ! defined NOMOVEW
c**********************

      subroutine movewater (maskwater, h, hb, hw, tw, 
     *                      heatw, baseperc, basefrml,
     *                      uw, vw, dt)

c Advects liquid water thickness (hw) due to horizontal flux between 
c adjacent cells, and applies local source/sink baseperc-basefrml. 
c (For now, no P-E on open water surfaces).
c
c Also advects water temperature (tw) for large water, or sets 
c heating advection tendency (heatw) for small water,
c used in vdif (where small-water tw is set).
c
c Thickness advection: (adapted from channelcompute in hydrology model)
c ~~~~~~~~~~~~~~~~~~~~
c Flow can be any of 8 quadrants. Also, applies negative volume 
c correction, tracked by enegsub, and computes diagnostic totosub, 
c flow from land to ocean. Flow equation is:  
c
c   d(hw)/dt = d/dx [ C(hw) hw^powdepth d(hb+hw+rhor*h]/dx^powslope ]
c
c Flow can occur out of a box in all directions, to avoid flip-flops
c at river junctions when just allowing only 1 (steepest) direction.
c Do time implicitly just with respect to hw term in the water
c potential slope.   
c
c 2 methods (if defines):
c channel_tridia:
c Implicit with contributions from all surrounding points. Tridiagonal
c solution, only for 1-D flowline. (6 out of the 8 directions in 
c iq,jq,etc, are wasted).
c
c channel_sparse:
c As channel_tridia but using sparse solver linbcg (Numerical Recipes)
c and row-indexed sparse storage method.
c
c Time implicit terms use linearized contribution from the slope term 
c in the flow equation. The exponent (slope**powslope) is set to 1 for 
c small slopes (< slopecrit) to avoid infinities at zero slope for 
c powslope < 1).
c
c Use nuptot and mupind to do this, index from full arrays into 
c linear sequence of selected points. 

c Temperature advection: 
c ~~~~~~~~~~~~~~~~~~~~~~
c Done at the end, time explicitly w.r.t. tw, in volume-flux form 
c using stored fluxes of hw across each cell face from hw calcs.

      use comicephys
      use comicegrid
      use comicesparse

      dimension
     *  maskwater(nx,ny), h(nx,ny),  hb(nx,ny), hw(nx,ny), tw(nx,ny),
     *  heatw(nx,ny), baseperc(nx,ny), basefrml(nx,ny),
     *  uw(0:nxp,0:nyp), vw(0:nxp,0:nyp) 

c Local arrays:
      dimension 
     *  hwprev(nx,ny), potw(nx,ny), zdt(nx,ny), 
     *  zfstora(nx,ny,8), zfstorb(nx,ny,8), work(nx,ny),
c       diagnostic: 
     *  streamflux(nx,ny),  ! stream flux out of cell (m3/s)
     *  istreamflux(nx,ny)  ! quadrant for streamflux (1-8)

#if defined CHANNEL_TRIDIA
      dimension atri(nx), btri(nx), ctri(nx), vech(nx), rhs(nx)
#elif defined CHANNEL_SPARSE
      dimension vech(nuvmax)          ! elspa,ijspa,etc, in comicesparse
      dimension rhs(nx*ny)
#endif
      dimension mupind(nx,ny), rhsexp(nx,ny)

      dimension iq(8), jq(8), kop(8), dists(8), dista(8)
      save iq, jq, kop
c directions to adjacent cells, in common in hydgrid.h:
c     k:         1   2   3   4   5   6   7   8
c anticlockwise from E:
c     data iq /  1,  1,  0, -1, -1, -1,  0,  1/
c     data jq /  0,  1,  1,  1,  0, -1, -1, -1/
c     data kop/  5,  6,  7,  8,  1,  2,  3,  4/
c sparse order (same as in sparse matrix in channeclcompute,
c lowest to highest i + (j-1)*nx):
      data iq / -1,  0,  1, -1,  1, -1,  0,  1/
      data jq / -1, -1, -1,  0,  0,  1,  1,  1/
      data kop/  8,  7,  6,  5,  4,  3,  2,  1/

      parameter (hwmin = hwcut)
      parameter (slopecrit = 1.e-10)
      logical first
      data first /.true./

c     statement function for flow coefficient:
      flowcoeff (d, dmin) = 1.  ! need to adjust as fn of d       777   

      call scopy (nx*ny, hw, 1, hwprev, 1)

c        Set nuptot and index mupind 

      nuptot = 0
      do j=1,ny
        do i=1,nx
#if defined CHANNEL_TRIDIA 
          if ( .true. ) then
#elif defined CHANNEL_SPARSE
          if ( .true. ) then  ! perhaps elim hw=0 & deep ocn pts  777   
#endif
            nuptot = nuptot + 1
            mupind(i,j) = nuptot
          else
            mupind(i,j) = 0
          endif
        enddo
      enddo
c     if (first) then
c       write (ioterm,*)
c    *     'movewater: nuptot=',nuptot, ' nx*ny=',nx*ny
c       first = .false.
c     endif

      call zero   (zfstora, nx*ny*8)
      call zero   (zfstorb, nx*ny*8)
      call resetr (zdt,     nx*ny, 1.)
      call zero   (rhsexp,  nx*ny)

#if defined CHANNEL_TRIDIA
      call zero (rhs,  nx)
      call resetr (atri, nx, 1.)
      call zero (btri, nx)
      call zero (ctri, nx)
#elif defined CHANNEL_SPARSE
      call zero (rhs, nuptot)
      call izero (ijspa, nspamax)
      call zero  (elspa, nspamax)
      do n1=1,nuptot
        elspa(n1) = 1.
      enddo
#endif

      do j=1,ny
        do i=1,nx
          zdt(i,j) = dt / darea(i,j)
          zlocal = dt * (baseperc(i,j) - basefrml(i,j)) !negl P-E 777   
     *             * (rhoice/rholiq)
          rhsexp(i,j) = rhsexp(i,j) + zlocal
          n1 = mupind(i,j)
          if (n1.ge.1) then
            rhs(n1) = rhs(n1) + zlocal
          endif
          potw(i,j) = hb(i,j) + hw(i,j) + rhor*h(i,j)
        enddo
      enddo

#if defined omp
!$OMP PARALLEL DO
!$OMP* shared  (iq,jq,dx,dy,hw,potw,
!$OMP*          zfstora,zfstorb)
!$OMP* private (i,j,i2,j2,k,dists,dista,zslope, zsabs,
!$OMP*          zdep,zwid,zcros,
!$OMP*          zcoef,zf,powslope,powconst,powdepth)
!$OMP* schedule(dynamic,2)
#endif

c        Compute zfstora, time-explicit coefficient multiplied
c        by potw(i2,j2)-potw(i,j)), for (i,j) and all its 8 
c        neighbors. Also zfstorb, time-implicit term that will be
c        multiplied by hw(i2,j2)-hw(i,j). To save cpu, only compute 
c        for 4 directions at each i,j, which gets all pairs in domain.
c        In subsequent do loops, all 8 directions used, so in next do 
c        loop below, zfstor[a,b](i,j,[5-8]) are copied from ...[1-4].

      do j=1,ny
        do i=1,nx

c            Look at 4 quadrants (gets all connections over domain
c            just once each)...save in zfstora,b

c~~~~~~~~~~~~~~~~~
          do k=1,4
c~~~~~~~~~~~~~~~~~
            i2 = i + iq(k)
            j2 = j + jq(k)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if (i2.ge.1.and.i2.le.nx.and.j2.ge.1.and.j2.le.ny) then
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              if (j.eq.j2) then
                dists(k) = 0.5*(dx(i,j)+dx(i2,j2))
                dista(k) = 0.5*(dy(i,j)+dy(i2,j2))
              else if (i.eq.i2) then
                dists(k) = 0.5*(dy(i,j)+dy(i2,j2))
                dista(k) = 0.5*(dx(i,j)+dx(i2,j2))
              else 
                dists(k) = sqrt (   (0.5*(dx(i,j)+dx(i2,j2)))**2 
     *                            + (0.5*(dy(i,j)+dy(i2,j2)))**2 )
                dista(k) = 0.25*(dx(i,j)+dx(i2,j2)+dy(i,j)+dy(i2,j2))
              endif

              zslope = potw(i2,j2)-potw(i,j) / dists(k)
              zsabs  = abs(zslope)

              if (zslope.gt.0.) then              ! use upwind thickenss
                zdep = hw(i2,j2) 
              else
                zdep = hw(i,j) 
              endif
              zwid = dista(k)
              zcros = zdep * zwid

c Either Manning:           ! use for ice-free hydrol, rivers ?   777   
c             if (h(i,j).eq.0. .and. h(i2,j2).eq.0.) then
c               zcoef = 2.5/.06    ! 2.5/.15 
c               powdepth = 2./3.
c               if (zsabs.gt.slopecrit) then
c                 powslope = 0.5
c                 powconst = 1.
c               else
c                 powslope = 1.
c                 powconst = slopecrit**(powslope - 0.5)
c               endif
c             else
c Or sub-ice, lakes:
                zcoef = flowcoeff (zdep, hwmin)
                powdepth = 0.
                powslope = 1.
                powconst = 1.
c             endif
              
              zf = ( zcoef * zcros * zdep**powdepth )
     *             * zsabs**(powslope-1.) / (powconst*dists(k))
              zfstora(i,j,k) = zf
              zfstorb(i,j,k) = zf*powslope

c++++++++++++++++
            endif
c++++++++++++++++

c~~~~~~~~~~~~~~
          enddo
c~~~~~~~~~~~~~~

        enddo
      enddo

c        Set arrays and rhs's needed to solve system of eqns for 
c        each method, using zfstor[a,b] from above. 
c        rhs, *tri, *spa are for subset of pts solved
c        time implictly (index n1,n2), and rhsexp(i,j) are for 
c        other points, solved time explicitly.
c        Can't parallelize for sparse method, since off-diag terms
c        in ijspa,elspa are assembled using incremented index kof.   

      kof = nuptot + 1
      do j=1,ny
        do i=1,nx
          kofprev = kof
          n1 = mupind(i,j)
          do k=1,8
            i2 = i + iq(k)
            j2 = j + jq(k)
            if ( (i2.ge.1.and.i2.le.nx.and.j2.ge.1.and.j2.le.ny) ) then
              n2 = mupind(i2,j2)

c             fill in other 4 directions of zfstora,b:
              if (k.gt.4) then
                zfstora(i,j,k) = zfstora(i2,j2,kop(k))
                zfstorb(i,j,k) = zfstorb(i2,j2,kop(k))
              endif

              rhsexp(i,j) = rhsexp(i,j)  + zdt(i,j)  * zfstora(i,j,k)
     *                                   * (potw(i2,j2)-potw(i,j))

              if (n1.ge.1) then
#if defined CHANNEL_TRIDIA
                rhs(n1)  = rhs(n1)       + zdt(i,j)  * zfstora(i,j,k)
     *                                   * (potw(i2,j2)-potw(i,j))
                atri(n1) = atri(n1) + zdt(i,j) * zfstorb(i,j,k)
c               all except k=2 and k=7 are not used (n2=0 for others): 
                if (i2.eq.i+1) then       ! superdiagonal
                  btri(n2) = btri(n2) - zdt(i,j) * zfstorb(i,j,k)
                else if (i2.eq.i-1) then  ! subdiagonal
                  ctri(n2) = ctri(n2) - zdt(i,j) * zfstorb(i,j,k)
                endif
#elif defined CHANNEL_SPARSE
                rhs(n1)  = rhs(n1)     + zdt(i,j) * zfstora(i,j,k)
     *                                 * (potw(i2,j2)-potw(i,j))
                elspa(n1) = elspa(n1)  + zdt(i,j) * zfstorb(i,j,k)
                ijspa(n1) = kofprev + 1 
                if (n2.ge.1) then
                  kof = kof + 1
                  elspa(kof) = elspa(kof) - zdt(i,j)*zfstorb(i,j,k)
                  ijspa(kof) = n2
                endif
#endif
              endif

            endif
          enddo
        enddo
      enddo

#if defined CHANNEL_SPARSE
      ijspa(nuptot+1) = kof + 1
#endif

c        Solve linear system arr*vech = rhs for vech (=hw-hwprev)

      call zero (vech, nuptot)

#if defined CHANNEL_TRIDIA

      call tridia_i (atri, btri, ctri, vech, rhs, nuptot, 0, 1.e-3)
c     call tridia_i (atri, btri, ctri, vech, rhs, nuptot, 1, 1.e-3) !777

#elif defined CHANNEL_SPARSE

      itol = 1
      tol = 1.e-11
      itmax = 10000   ! 1000
      call linbcg (nuptot, rhs, vech, itol, tol, itmax, itlin, errlin)
      if (abs(errlin).gt.2.*tol) write (ioterm,*)
     *  '*** Warning linbcg(hw): itlin=',itlin,' errlin=',errlin

#endif

c        Set new water thickness hw

      do j=1,ny
        do i=1,nx
          n1 = mupind(i,j)
          if (n1.ge.1) then
            hw(i,j) = hwprev(i,j) + vech(n1)
          else
            hw(i,j) = hwprev(i,j) + rhsexp(i,j)
          endif
        enddo
      enddo

c        Adjust time explicit pts for neighboring time implicit changes
c        (conserves water)

      do j=1,ny
        do i=1,nx
          n1 = mupind(i,j)
          if (n1.eq.0) then
            do k=1,8
              i2 = i + iq(k)
              j2 = j + jq(k)
              if (i2.ge.1.and.i2.le.nx.and.j2.ge.1.and.j2.le.ny) then
                n2 = mupind(i2,j2)
                if (n2.ge.1) then
                  hw(i,j)= hw(i,j) + zdt(i,j)*zfstorb(i,j,k)
     *                               *(hw(i2,j2)-hwprev(i2,j2))
                endif
              endif
            enddo
          endif
        enddo
      enddo

c        Crude fix for negative water thickness

      enegsub = 0.
      do j=1,ny
        do i=1,nx
          if (hw(i,j).lt.0.) then
c           write (ioterm,*) 'neg vol i,j=',i,j
            enegsub = enegsub + (-hw(i,j)) * darea(i,j) / dt
            hw(i,j) = 0.
          endif
        enddo
      enddo

c        Set totosub (diagnostic, total stream flux from land to ocean,
c        m3/y), and streamflux,istreamflux (diagnostic, total stream 
c        flux (m3/y) and average quadrant *out* of each grid cell).
c
c        Also compute heat advection out of the cell (zhadv), 
c        and step temperature tw due to advection, or set tendency
c        heatw for small water.

      totosub = 0.
      call zero (streamflux, nx*ny)
      call zero (work, nx*ny)
      call zero (uw, (nxp+1)*(nyp+1))
      call zero (vw, (nxp+1)*(nyp+1))

      do j=1,ny
        do i=1,nx
          n1 = mupind(i,j)

          zhadv = 0.
          do k=1,8
            i2 = i + iq(k)
            j2 = j + jq(k)
            if ( (i2.ge.1.and.i2.le.nx.and.j2.ge.1.and.j2.le.ny) ) then
              n2 = mupind(i2,j2)

c                zs is h(upstream)*width*veloc out of each face (m3/y)

              zs = zfstora(i,j,k) * (potw(i,j)-potw(i2,j2))

c             correction for time implicit:
              if (n1.ge.1) then
                zs = zs + zfstorb(i,j,k)*vech(n1)
                if (n2.ge.1) zs = zs - zfstorb(i,j,k)*vech(n2)
              endif

              if (zs.gt.0.) then
                streamflux(i,j) = streamflux(i,j) + zs
                work(i,j) = work(i,j) + k*zs
              else 
                streamflux(i2,j2) = streamflux(i2,j2) - zs
                work(i2,j2) = work(i2,j2) - kop(k)*zs
              endif

c             land/lake to ocean:
              if (maskwater(i,j).ne.1.and.maskwater(i2,j2).eq.1) then
                totosub = totosub + zs
              endif

c                zhadv is T(upstream)*h*width*veloc out of cell (K m3/y)

              zhadv = zhadv - max(zs,0.)*tw(i,j) - min(zs,0.)*tw(i2,j2) 

c                Deduce water velocities uw,vw on u and v-grids 
c                (just set velocs at the rh and top faces of each 
c                h-grid cell, this will get all u and v-grid points).
c                Do only if surrounding hw's are large (used in icedyn).
c                Could try to factor in diagonal velocities ?     777   

              if (i2.eq.i+1 .and. j2.eq.j) then
                if (hw(i,j).gt.hwmin .and. hw(i2,j2).gt.hwmin) then
                  zdist =  0.5*(dx(i,j)+dx(i2,j2))
                  if (zs.gt.0.) then ! upstream
                    zdep = hwprev(i,j)
                  else
                    zdep = hwprev(i2,j2)
                  endif
                  uw(i,j) = zs / ( max(zdep,hwmin) * zdist )
                endif

              else if (i2.eq.i .and. j2.eq.j+1) then
                if (hw(i,j).gt.hwmin .and. hw(i2,j2).gt.hwmin) then
                  zdist =  0.5*(dy(i,j)+dy(i2,j2))
                  if (zs.gt.0.) then ! upstream
                    zdep = hwprev(i,j)
                  else
                    zdep = hwprev(i2,j2)
                  endif
                  vw(i,j) = zs / ( max(zdep,hwmin) * zdist )
                endif
              endif

            endif
          enddo

c            Set new water temperature (large) or heat tendency (small)

          if (hw(i,j).gt.hwmin) then
            ztme = tmelt - dtmdh*h(i,j)
            ztnew =  ( tw(i,j)*hwprev(i,j) + zhadv*zdt(i,j) 
     *                 + ztme*(baseperc(i,j)-basefrml(i,j))*dt 
     *                       *(rhoice/rholiq)
     *               ) / hw(i,j)
            tw(i,j) = ztnew
            heatw(i,j) = 0.
          else if (hw(i,j).gt.1.e-3) then
c           small and no water tw will be set in vdif
            heatw(i,j) = rholiq*cheatliq*zhadv*zdt(i,j) / dt    ! J/m2/y
          else
            heatw(i,j) = 0.
          endif

        enddo
      enddo

      do j=1,ny
        do i=1,nx
          if (streamflux(i,j).gt.0.) then
            istreamflux(i,j) = nint(work(i,j) / streamflux(i,j))
          else
            istreamflux(i,j) = 0
          endif
        enddo
      enddo

      return
      end

c*****
#endif
c*****
