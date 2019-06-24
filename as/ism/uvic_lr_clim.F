      module get_lisiecki     
      real lr_wt(2115)
      real lr_iso(2115)
      real lr_age(2115)     
      end module
      
      subroutine set_uvic_lr_clim
      
      use comicegrid
      use get_lisiecki
      
      implicit none
      
      real mod_val,lgm_val,diff
      integer unit, n
      
      !load get_lisiecki weight curve      
      open (unit=2,file='data/lr_age.txt')
      read (2,*,end=10) lr_age
 10   continue
      close(2)
      lr_age(:)=-lr_age(:)    
      open (unit=3,file='data/lr_iso.txt')
      read (3,*,end=20) lr_iso
 20   continue
      close(3)      

      mod_val=lr_iso(2115) !present-day
      lgm_val=lr_iso(2097) !18000 kabp
      diff=lgm_val-mod_val
      !get weight array: 0=lgm, 1=modern
      
      do n=1,2115
        lr_wt(n)=(lgm_val-lr_iso(n))/diff
      enddo
      return 
      end
      
      subroutine uvic_lr_clim(zclim,sealev)

      use comicesheet, only: timeice
     &, budgsnow, budgrain, budgevap, budgmelt 
     &, tsurf, tsurfi
      use get_lisiecki  
          
      implicit none
      
      integer n
      real t,w,diff,zclim,sealev
      
      include "size.h"
      include "ism.h"
      
      !Find Liesicki-Raymo times that bracket current ice time
      if (timeice/1000. .lt. lr_age(1)) then
        print*, 'Whoops, timeice is .lt. first lg_age.'
	stop
      elseif (timeice/1000. .gt. lr_age(size(lr_age))) then 
        print*, 'Whoops, timeice is .gt. last lg_age.'
	stop      
      else
        !get bracketing LR times
	n=1
        do while (lr_age(n).lt.timeice/1000.)
          n=n+1
        enddo 	
	!get weight between times.  Weight is 0 at timeice=lr_age(n).  Weight is 1 at timeice=lr_age(n+1) 
	diff=lr_age(n)-lr_age(n-1)
	w=(timeice/1000.-lr_age(n-1))/diff	
	if (w.lt.0.) stop 'Error: w .lt. 0'
	!get zlcim value for oceanmelt from interpolated lr_wt values
	zclim=lr_wt(n-1)*(1.-w)+lr_wt(n)*w	
	!limit to between 0 and 1 (0=lgm, 1=mod)
	zclim=max(0.,(min(1.,zclim)))
	!Get interpolated sealevel between -125m and 0m
	sealev=-125*(1.-zclim)	
      endif   
      zclim_out=zclim
      sealev_out=sealev
      return 
      end
