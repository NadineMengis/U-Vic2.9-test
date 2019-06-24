#if !defined O_ism_standalone
      subroutine setism

!=======================================================================
!     High level subroutine to initialize all ice sheets.  
!=======================================================================      

      use ice_sheet_storage, only: strtis,endis
           
      implicit none

      include "size.h"
      include "ism.h" 
      include "switch.h"
      include "atm.h"
      include "ice.h"
      
      character(120) g_st
      
      integer n,ntrec,ioun
      logical runstep,runais,rungis
      namelist /ismin/ runais,rungis

      runais=.false.
      rungis=.false.
      strtis=0.
      endis=0.
      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, ismin, end=101)
101   continue      
      if (runais .and. rungis) then
        strtis=1
        endis=2
      elseif (runais .and. .not. rungis) then
        strtis=1
        endis=1
      elseif (.not. runais .and. rungis) then
        strtis=2
        endis=2
      else
        stop 'Error: neither GIS/AIS specified in control.in.'
      endif
#if defined O_use_bias_SAT_correction && defined O_data_sat
      stop 'Error: both O_use_bias_SAT_correction, O_data_sat defined'
#endif
      isb(:)=0     
            
      do n=strtis,endis
	!initialize ice sheet derived data type arrays
	call set_ice_sheet_DDT(n) 
	call set_pointers(n)	
	!call ice sheet model initialization
	call ism(n)         			
	!initialize ice drainage basins
	call seticebasin(n)	
	!get initial ice basin ice volumes	
	call getbasinvol(n)   
      enddo
      !re-initialize the elevation bins, even if read in from a restart (which has already-updated bins).
      !Note, this was originally only called for initial conditions, but now done redundantly for both
      !initial and restart-initialized runs.
      do n=strtis,endis 	
        call set_pointers(n)  
        call elevbin(n) 	
      enddo
      do n=strtis,endis   	
        call set_pointers(n)  
	call elevsort(n)
      enddo

      call recalc_umsk
	
      !ensure nmbal is zerod... this is somewhat redundant if the restart was
      !written immediately after an ice step, which it should be if ice is called
      !at end of the year
      nmbal(:,:,:)=0. 

      !ensure that no remnant sea ice exists under ice shelves
      where (tmsk .le. 0.) aiceocn(:,:,1)=0.
      where (tmsk .le. 0.) aiceocn(:,:,2)=0.   
      where (tmsk .le. 0.) hice(:,:,1)=0.
      where (tmsk .le. 0.) hice(:,:,2)=0.	       

      !call fluxsort.  Sort of redundant: recalculates 
      !diceheat,isflxm,isflxh.  But needed for global sums 
      !as it also calculates isvol,disvol (and diagnostics 
      !nisvol, sisvol.  So really only diagnostic here.
      runstep=.false.
      call fluxsort(runstep)	
      
      call rivinit
#if defined O_icefree_subgrid_topog
      do n=strtis,endis       
        call set_pointers(n) 
        call bareland(n)	      
      enddo	 
#endif
#if defined O_use_bias_SAT_correction
      do n=strtis,endis       
        call set_pointers(n) 
        call getbias(n) 	    
      enddo    
#endif
      return
      end

      subroutine runism      

!=======================================================================      
!     High level subroutine to run all ice sheets
!=======================================================================      
      use ice_sheet_storage, only: strtis,endis
      
      implicit none     
   
      include "size.h"
      include "ism.h"        
      integer n	
      logical runstep       

      do n=strtis,endis	
        !point to ice sheet
	call set_pointers(n)
	!obtain downscaled climate fields
	call downscale(n)  !pre-dynamics total drainage basin volume set here
      enddo
      do n=strtis,endis
        !point to ice sheet
	call set_pointers(n)
	!run 1 or more ice sheet steps
 	call ism(n)	
      enddo

      !zero basin moisture/heat banks here because residual snowmelt gets added to them in 
      !subroutine elevsort.
      dicevol(:)=0.
      diceheat(:)=0.
      
      do n=strtis,endis
        !point to ice sheet
        call set_pointers(n)
	!accumulate elevation bins, allocate subgrid arrays	       
        call elevbin(n) 
      enddo	
      do n=strtis,endis 
        !point to ice sheet    
        call set_pointers(n)
	!fill allocated subgrid arrays. resdistribute snow cover  
        call elevsort(n)	       
      enddo

      !get updated basin volume
      do n=strtis,endis  
        !set pointers 
	call set_pointers(n)
	!calculate ice volumes in each ice drainage basin        
	call getbasinvol(n)		
      enddo

      runstep=.true.
      !determine change in size of each ice drainage basin and resulting moisture fluxes to climate model
      call fluxsort(runstep)
        
      call recalc_umsk
      call rivinit      
#if defined O_mtlm
      call recalc_mtlm_points	
#endif     
      call recalc_land_albedo	   

      nmbal(:,:,:)=0. 
      ismtice(:,:)=0.
#if defined O_use_bias_SAT_correction
      do n=strtis,endis       
	call set_pointers(n) 
	call getbias(n) 	    
      enddo    
#endif      
      return 
      end
#endif
