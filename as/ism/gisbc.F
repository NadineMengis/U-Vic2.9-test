      subroutine gisbc
      
!=======================================================================
!     Accumulation boundary conditions for ice sheet model      
!=======================================================================       
      
      use subgrid 
      
      implicit none
      
      include "size.h"
      include "ism.h"
      include "param.h"
      include "atm.h"
      include "ice.h"
      include "csbc.h"
      include "tmngr.h"
      include "cembm.h"
#if defined O_mon_islapse       
      real monlapse(12)
      data monlapse 
     &/7.9e-5,8.9e-5,7.9e-5,7.3e-5
     &,5.9e-5,4.7e-5,4.6e-5,5.7e-5
     &,6.9e-5,7.3e-5,6.5e-5,7.6e-5/
      integer mon    
#endif   
      integer i,j,n
      integer ntrec
      data ntrec /0/
      save ntrec
      character(120) g_st
      
      integer monprev,monnext
      real modtime,wtnext,wtprev

#if defined O_ism_pdd      
      integer m
      real temp(imt,jmt),prec(imt,jmt),month_accum
      save m,month_accum,temp,prec      
      logical firstcall
      data firstcall /.true./
      save firstcall

      if (firstcall) then

        !Initialize monthly mean accumulator values	 
        m = 1
        month_accum = 0.0
        year_accum(:)=0.0
        temp(:,:)   = 0.0
        prec(:,:)   = 0.0
	
	firstcall = .false.
      endif   

      !Accumulate monthly averages    
      do j=1,jmt
        do i=1,imt
      
        !accumulate P
        prec(i,j)=prec(i,j)+precip(i,j)
      
        !accumulate T
        temp(i,j)=temp(i,j)+at(i,j,2,isat)-elev(i,j)*rlapse

        enddo
      enddo

      !increment the intramonthly counter
      month_accum=month_accum+1.
      
      !print*,'sum(prec)=',sum(prec)
      !print*,'sum(temp)=',sum(temp) 
      !print*,'month_accum=',month_accum
      !print*,'eomon=',eomon   
      
      if (eomon) then
        do j=1,jmt
          do i=1,imt
	    gprecmon(i,j,m) = gprecmon(i,j,m) + (prec(i,j)
     &                       /month_accum*10/secday) !g/cm^2/sec => mm/day    
            gtempmon(i,j,m) = gtempmon(i,j,m) + (temp(i,j)        
     &                       /month_accum  )         !C 
	  enddo
	enddo
        
	call embmbc(gtempmon(:,:,m))
        call embmbc(gprecmon(:,:,m))
         
	!rezero accumulating arrays, intramonthly counter
	month_accum = 0.0
	prec(:,:) = 0.0
	temp(:,:) = 0.0	
	
	!up accumulation count for that month by one (occurs once a year)	
	year_accum(m) = year_accum(m)+1

	!move to the next index of the 12 month array	
	m = m + 1

	!if 12 months filled, go back to first index of 12 month array
	if (m.eq.13) m = 1
	
      endif
#endif

      do j=2,jmtm1
        do i=2,imtm1
          !Get subgrid ice temperature boundary condition	
	  do n=strtiE(i,j),endiE(i,j)
            !Ice temperature
            ismticeiE(n) = ismticeiE(n) + ticeiE(n)
	  enddo
	  ismtice(i,j) = ismtice(i,j) + ticel(i,j)
        enddo
      enddo

      !counter
      iac=iac+1

#if defined O_mon_islapse      
      !Determine monthly GIS lapse rate.
      !Assumes restarts are written Dec. 31.
      !Lapse rates given by Fausto et al. (2009) as observed GIS averages.
      mon=min(12,max(1,ceiling(dayoyr/365.*12.)))
      !set the lapse rate to the lowest of the default lapse and the monthly lapse.  
      !This is meant to decrease the summertime lapse rate over GIS.
      rlapseGIS=min(rlapse,monlapse(mon))
#endif
      ntrec=ntrec+1

      return
      end
