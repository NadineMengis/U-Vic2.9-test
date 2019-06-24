      module subgrid
     
      integer numiE       !length of 1D ice subgrid elevation array
      real :: dzi = 100.  !vertical resolution of subgrid ice bins
      real :: rtag        !ice sheet/climate consistency restart tag
      real :: himin = 1. !min ice thickness before included in subgridded ice (m)
      
      !subgrid elevation/fraction arrays
      real, allocatable :: eleviE(:)   !(cm): elevation of bin
      real, allocatable :: fraciE(:)   !% bin coverage of climate cell 
      real, allocatable :: iariE(:)      
      integer, allocatable :: isniE(:) !indexes each subgrid bin to ice sheet instance
      
      !subgrid climate arrays
      real, allocatable :: hsnowiE(:,:)
     &	                  ,ticeiE(:)
     &                    ,ismticeiE(:)
     &			  ,ithkiE(:)
     &      		  ,mbalmelt(:,:)
     &		          ,mbalaccum(:,:)
     
      !diagnostic/temporary arrays
      real, allocatable :: precipiE(:)
     &		          ,satiE(:)
     &		          ,albedoiE(:)
     &		          ,mextiE(:)
     &		          ,mduriE(:)

      !temporary remapping arrays
      integer, allocatable :: ibin(:,:,:)       
      real, allocatable :: fraciEsav(:)
     &		          ,hsnowiEsav(:)
     &		          ,fraclEsav(:)
     &		          ,hsnowEsav(:)
     &		          ,iariEsav(:)			  
     &		          ,ticeiEsav(:)
     
#if defined O_use_bias_SAT_correction
      real, allocatable :: biaiE(:,:)
     &		          ,totariE(:)
     &			  ,sbciE(:)      
#endif
#if defined O_refreeze
      real, allocatable ::    meltprev(:)
      real, allocatable ::    hnewsno(:)
     &			  ,   htemp(:)
     &			  ,   counter(:)
     &			  ,   hreftot(:)
     
      real, allocatable ::    meltprevsav(:)
      real, allocatable ::    hnewsnosav(:)
     &			  ,   htempsav(:)
     &			  ,   countersav(:)
     &			  ,   hreftotsav(:)     
#endif        
#if defined O_topography_shading
      real, allocatable ::    snv(:,:) !average components of normal vector to ice sheet surf slope
     &			  ,   zen_h(:,:)  !horizontal components of zenith angle, on ice grids
#endif
#if defined O_subgrid_out
      real, allocatable ::    dswr_out(:)
     &		          ,   ultnt_out(:)
     &	         	  ,   usens_out(:)
     &   		  ,   ulwr_out(:)
     &   		  ,   pdd_fac(:)     
#endif      


      end module
      
      
      
