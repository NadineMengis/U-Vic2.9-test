      subroutine bareland(isn)

      use comicegrid
      use comicephys, only: rhoice, rhobed
      use comicesheet, only: hs, h, hb
      use ice_sheet_storage, only: strtis,endis      
      use subgrid 
#if defined O_icefree_subgrid_topog              
      implicit none

      include "size.h"
      include "param.h"          

      !------------------------------------------------------------
      !Use initial state of ice sheet model to determine ice free
      !topography, assuming isostatic equilibrium once ice is removed.
      !Once this topography is determined, create subgrid-binned
      !elevations for the ice-free topography.  Output results to 
      !netCDF files.

      integer i,j,n,nn,ntrec,isn, count,jstrt,jend
      real rb(nx,ny),diff(nx,ny)
      real sg(imt,jmt,50), avgelev(imt,jmt,50)
      real areatot(imt,jmt),sum(imt,jmt)
      real binsum(imt,jmt), gae(imt,jmt)
      real midelev
      real min(imt,jmt), max(imt,jmt), relief(imt,jmt), ed(imt,jmt)
      real jump, elevtest
      logical found
      character(120) g_st
      
      save     

      !If we're running through first ice sheet, zero...      
      if (isn==strtis) then
        !relief bin array
	ed(:,:) = 0. 
        !subgrid elevation bin array
	sg(:,:,:) = 0.
	!average elevation bin array
	avgelev(:,:,:) = 0.
	!relief array
	relief(:,:) = 0.
	!initial guesses for relief search
	min(:,:)=1.e10
        max(:,:)=0.
      endif
      
      !Loop over ice sheet grid, and calculate rebounded ice-free
      !topography at equilibrium.
      do i=1,nx
        do j=1,ny
	  rb(i,j) = hb(i,j)+h(i,j)*rhoice/rhobed
	enddo
      enddo
      !Limit rebounded topography to .ge. 1.
      where (rb .lt. 0.) rb=1.      

      !get total vertical relief in each climate model grid cell by
      !updating min and max elevations of all global cells that
      !are 'parents' of ice sheet grid cells.  Recall: iceiind(nx,ny)
      !and icejind(nx,ny) are arrays on the ice sheet grid that contain
      !the indexes of the 'parent' ice sheet grid.
      do i=1,nx
        do j=1,ny
	  if (rb(i,j) .lt. min(iceiind(i,j),icejind(i,j))) then
	    min(iceiind(i,j),icejind(i,j)) = rb(i,j)
	  endif
	  if (rb(i,j) .gt. max(iceiind(i,j),icejind(i,j))) then
	    max(iceiind(i,j),icejind(i,j)) = rb(i,j)
	  endif	  
	enddo
      enddo
      
      !get total relief in grid cell.  Do hemisphere-at-a-time.
      !(as this routine is called once per ice sheet instance).
      !Bin reliefs into 1000-m intervals, in order to determine 
      !subgrid elevation bin delta z.
      if (isname=='Antarctica') then
        jstrt=2
        jend=51
      elseif (isname=='Greenland') then
        jstrt=52
        jend=101
      endif       
       
      do i=2,imtm1
        do j=jstrt,jend
          relief(i,j)=max(i,j)-min(i,j)   
          if(relief(i,j).lt.-1.e9) relief(i,j)=0. 
          if(relief(i,j) .ge. 0. .and. relief(i,j) .lt. 1000) ed(i,j)=1.
          if(relief(i,j) .ge. 1000. .and. relief(i,j) .lt. 2000) 
     &	    ed(i,j)=2.
          if(relief(i,j) .ge. 2000. .and. relief(i,j) .lt. 3000) 
     &	    ed(i,j)=3.
          if(relief(i,j) .ge. 3000. .and. relief(i,j) .lt. 4000) 
     &	    ed(i,j)=4.
          if(relief(i,j) .ge. 4000. .and. relief(i,j) .lt. 5000) 
     &	    ed(i,j)=5.
        enddo
      enddo
      
      !Determine which elevations are present in each climate grid
      !cell, and their areal fractions, based on the dz determined
      !above.
      do i=1,nx
        do j=1,ny
	  n=1
	  !set dz to the distance determined in previous loop.
          jump=ed(iceiind(i,j),icejind(i,j))
          found = .false.
          do while (.not. found)
            elevtest=(n-1)*100.
	    
	    !if the rebounded elevation is .le. the test elevation
	    !value add the area of this ice sheet grid cell to the
	    !'parent' climate model subgrid area bin, accumulated 
	    !climate model total area array.  Also add area-weighted
	    !elevation to average elevation bin.
            if (rb(i,j).le.elevtest) then
              sg(iceiind(i,j),icejind(i,j),n) = 
     &        sg(iceiind(i,j),icejind(i,j),n)+darea(i,j)  

              areatot(iceiind(i,j),icejind(i,j)) = 
     &        areatot(iceiind(i,j),icejind(i,j))+darea(i,j)

              avgelev(iceiind(i,j),icejind(i,j),n) = 
     &        avgelev(iceiind(i,j),icejind(i,j),n) +
     &        rb(i,j)*darea(i,j)
              found = .true.
            else
              n=n+jump
            endif	
          enddo
        enddo
      enddo
      
      !if finished running final ice sheet instance, calculate
      !final average elevation for each elevation bin, 
      !fractional area taken up by each 
      !elevation bin, total fractional area coverage (diagnostic
      !which should equal 1 where total ice sheet model coverage 
      !occurs), and total number of bins in a climate model grid
      !cell.
      if (isn == endis) then   
        sum(:,:) = 0. 
	binsum(:,:) = 0.
	gae(:,:) = 0.  
        do i=2,imtm1
          do j=2,jmtm1
	    do n=1,50
	      if (sg(i,j,n) .gt. 0.) then
		avgelev(i,j,n) = avgelev(i,j,n)/sg(i,j,n)
	        sg(i,j,n) = sg(i,j,n)/areatot(i,j)
	        sum(i,j) = sum(i,j) + sg(i,j,n)
		binsum(i,j) = binsum(i,j)+1	
	      endif
	    enddo
	  enddo
        enddo
	
	!Calculate final climate model grid cell average 
	!elevation.
        do i=2,imtm1
          do j=2,jmtm1	
	    do n=1,50
	      gae(i,j) = gae(i,j)+avgelev(i,j,n)*sg(i,j,n)
	    enddo
	  enddo
	enddo	

      endif
      
      !********DIAGNOSTICS*********!
      !Choose which you want to output!

!        g_st = 'relief.nc'
!        call jerncdf_snapshot(relief,g_st)  
        g_st = 'elevinterval'
        call jerncdf_snapshot(ed,g_st)           
            
! 
!       if (isn==strtis) then
!         diff(:,:) = rb(:,:)-hb(:,:)
!         g_st = 'diff.nc'
!         call jerncdf_ice_snapshot(diff,g_st)     
!         g_st = 'hs.nc'
!         call jerncdf_ice_snapshot(hs,g_st)
!         g_st = 'hb.nc'
!         call jerncdf_ice_snapshot(hb,g_st)            
!         g_st = 'rb.nc'
!         call jerncdf_ice_snapshot(rb,g_st)
!         g_st = 'h.nc'
!         call jerncdf_ice_snapshot(h,g_st)	
!       endif
     
       if (isn .eq. endis) then 
!          g_st = 'sum.nc'
!          call jerncdf_snapshot(sum,g_st) 
!          g_st = 'relief.nc'
!          call jerncdf_snapshot(relief,g_st)           
!          ntrec=0
! 	 do n=1,50
! 	 ntrec=ntrec+1
! 	 g_st='elevfrac'
! 	 call jerncdf_timedep(sg(:,:,n),g_st,ntrec)
! 	 g_st='avgelev'
! 	 call jerncdf_timedep(avgelev(:,:,n),g_st,ntrec)       
! 	 enddo        
! 	 g_st = 'numofbins'
! 	 call jerncdf_snapshot(binsum,g_st)
	 g_st='gae'
	 call jerncdf_snapshot(gae,g_st) 
       endif   
      
      !********DIAGNOSTICS END*****! 
#endif      
      return
      end
      
      
