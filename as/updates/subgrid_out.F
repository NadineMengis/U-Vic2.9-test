#if defined O_subgrid_out
      subroutine subgrid_out
      
      use subgrid
      
      implicit none
      
      include "size.h" 
      include "ice.h"
      include "temp.h"
         
      real default_out(imt,jmt)
	   
      character(120) fname,vname
      integer ntrec
      data ntrec /0/
      save ntrec      
      
      default_out(:,:)=0.
      ntrec=ntrec+1
      
      fname='seasonal_ism_diagnostics'
      
      vname='albedo'
      call jerncdf_sgi2ice_timedep
     &(albedoiE,non_grid_surf_alb,fname,vname,ntrec)
      
      vname='hsnoliE'
      call jerncdf_sgi2ice_timedep
     &(hsnowiE(:,2),hsnol(:,:,2),fname,vname,ntrec)
     
      vname='dswr_out'
      call jerncdf_sgi2ice_timedep
     &(dswr_out,default_out,fname,vname,ntrec)    
     
      vname='dswr_out'
      call jerncdf_sgi2ice_timedep
     &(dswr_out,default_out,fname,vname,ntrec) 
     
      vname='usens_out'
      call jerncdf_sgi2ice_timedep
     &(usens_out,default_out,fname,vname,ntrec)          
     
      vname='ulwr_out'
      call jerncdf_sgi2ice_timedep
     &(ulwr_out,default_out,fname,vname,ntrec)
     
      vname='mextiE'
      call jerncdf_sgi2ice_timedep
     &(mextiE,default_out,fname,vname,ntrec) 
     
      vname='albedoiE'
      call jerncdf_sgi2ice_timedep
     &(albedoiE,default_out,fname,vname,ntrec)          

      vname='mbalmelt'
      call jerncdf_sgi2ice_timedep
     &(mbalmelt,default_out,fname,vname,ntrec) 

      vname='mbalaccum'
      call jerncdf_sgi2ice_timedep
     &(mbalaccum,default_out,fname,vname,ntrec)
     
      vname='biaiE'
      call jerncdf_sgi2ice_timedep
     &(biaiE,default_out,fname,vname,ntrec)      

    !  vname='hreftotsav'
    !  call jerncdf_sgi2ice_timedep
    ! &(hreftotsav,default_out,fname,vname,ntrec)   
     
      vname='hsnowiE'
      call jerncdf_sgi2ice_timedep
     &(hsnowiE(:,2),default_out,fname,vname,ntrec)      
      vname='ticeiE'
      call jerncdf_sgi2ice_timedep
     &(ticeiE,default_out,fname,vname,ntrec)
     
      vname='ithkiE'
      call jerncdf_sgi2ice_timedep
     &(ithkiE,default_out,fname,vname,ntrec)
     
      vname='satiE'
      call jerncdf_sgi2ice_timedep
     &(satiE,default_out,fname,vname,ntrec)      

      vname='pdd_fac'
      call jerncdf_sgi2ice_timedep
     &(pdd_fac,default_out,fname,vname,ntrec)      
      
      return
      end
#endif
