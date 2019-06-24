      subroutine downscale(isn)

 !=======================================================================
 !     Downscale SUBGRIDDED UVic ESCM arrays to ice sheet grid using pre-generated indexing. 
 !=======================================================================  
      use comicegrid, only: nx,ny,is2sg,isname,iceiind,icejind
     &, darea, mdur, mext
      use comicesheet, only: tsurf, tsurfi, h, hw, icedrainage
     &, budgsnow, budgrain, budgevap, budgmelt      
      use subgrid
      use ice_sheet_storage, only: strtis,endis
 
      implicit none

      include "size.h"
      include "atm.h"
      include "ism.h"
      include "cembm.h"
      include "ice.h"
      include "grdvar.h"
      include "param.h"
      include "levind.h"

      integer i,j,n,isn,ii,jj,jstrt,jend
      real fnis,tmp,basinbudg,s2i,totar,nmbalvol,accum,melt
      real accumsum, meltsum
      real temp(nx,ny),maski(nx,ny)
      real nsgsmb,sgsmb
      real icebudgsum,climbudgsum,dmsk(imt,jmt),icesum
      character(120) g_st
      integer ntrec
      data ntrec /0/
      save ntrec      
      logical check
      
      s2i=rhosno/rhoice
 
      !Loop through ice sheet grid, and obtain a mbal value from 'parent'
      !subgridded elevation bin.
      budgsnow(:,:)=0.
      budgmelt(:,:)=0.
      tsurf(:,:)=273.15
      tsurfi(:,:)=273.15
      temp(:,:)=0.
      mext(:,:)=0.
      icebudgsum=0.
      accum=0.
      melt=0.
      accumsum=0.
      meltsum=0.
      do i=1,nx
        do j=1,ny	  
          !if ice exists, obtain surface boundary conditions using is2sg indexing array.
          if (h(i,j) .gt. himin) then    
            melt=mbalmelt(is2sg(i,j),2)*s2i/100. !melt in units of m ice
	    accum=mbalaccum(is2sg(i,j),2)*s2i/100. !ditto for accum
	    
	    !Set net SMB in budgsnow.  Note: possible here for melt to exceed available ice thickness.  
	    !The negative SMB in this case is reduced within subroutine icedyn.  However, since melt already
	    !provided a moisture flux to ocean (and a heat flux out of atmosphere) during subgridded EBM operation, 
	    !these are compensated for by extracting/adding the equivalent moisture/heat flux from the ocean 
	    !the following year (in practice, a negative calving).
	    budgsnow(i,j)=melt+accum
            tsurfi(i,j)=ismticeiE(is2sg(i,j))/max(iac,1.)+273.15
            tsurf(i,j)=tsurfi(i,j)

            !generate some diagnostics
            mext(i,j)=mextiE(is2sg(i,j))
	    mdur(i,j)=mduriE(is2sg(i,j))
	    meltsum=meltsum+melt*darea(i,j)
	    accumsum=accumsum+accum*darea(i,j)	    	    
	    
	    !generate some debugging diagnostics
	    icebudgsum=icebudgsum+budgsnow(i,j)*darea(i,j)
          endif  	  
        enddo
      enddo
      
      if (isname=='Antarctica') then
        jstrt=2
        jend=51
      elseif (isname=='Greenland') then
        jstrt=52
        jend=101
      endif         
      check=.false.             
      do i=2,imtm1
        do j=jstrt,jend
        totar=0.
          !Get volume of ice in nmbal array that needs to be passed to an ice array
        nmbalvol=nmbal(i,j,2)*dxt(i)*cst(j)*dyt(j)/1.e6*s2i !m3 ice
	  if (nmbalvol .gt. 0.) then 
	    maski(:,:)=0.
	    do ii=1,nx
              do jj=1,ny
	        !if ice sheet grid cell is referenced to parent global cell
	        if (iceiind(ii,jj)==i .and. icejind(ii,jj)==j) then
                  if (h(ii,jj).le.himin .and. hw(ii,jj).le.0.) then
                    !Get total ice free dry land area of ice children cells
		    totar=totar+darea(ii,jj)
		    !set maski in at ii,jj to 1, indicating that this is a point that gets some of the 
		    !non-subgridded mass balance volume nmbalvol
		    maski(ii,jj)=1.
                  endif
                endif
              enddo
            enddo
	    
	    if(totar.gt.0) then
	      accum = nmbalvol/totar !nmbalvol in m^3 ice !totar in units of m^2
	    else
	      stop 'Error: totar=0, nmbalvol>0'
	    endif
	    icesum=0.	    
            do ii=1,nx
              do jj=1,ny
	        if (maski(ii,jj).gt.0.5) then
		      if (budgsnow(ii,jj).ne.0.) then
		        stop 'Error: budgsnow .ne. zero in non-sg accum'
		    else
	            budgsnow(ii,jj) = accum
		    tsurf(ii,jj)=ticel(i,j)/max(iac,1.)+273.15
                  endif
		  icesum=icesum+budgsnow(ii,jj)*darea(ii,jj)
                endif
              enddo
            enddo
	    if (abs(icesum-nmbalvol).gt.1.) then        
	      print*, 'WARNING'
	      print*, 'i=',i,'j=',j,': diff in volume here'
	      print*, 'climsum=',nmbalvol
	      print*, 'icesum=', icesum
	      check=.true.
	    endif
	    icebudgsum=icebudgsum+icesum
	  endif 
        enddo
      enddo
      if (check) stop 'Stopped in scale.F.'
      !Set surface temperatures everywhere over ice to a maximum of 273.15K and set surface ice temperature to surface temperature.
      do i=1,nx
        do j=1,ny
	  tsurf(i,j)=min(tsurf(i,j),273.15)
	  tsurfi(i,j)=tsurf(i,j)
      	  budgrain(i,j)=0.
          budgevap(i,j)=0.
	  budgmelt(i,j)=0.
	enddo
      enddo

      !accumulate total climate-side mass of ice that should be passed to ice sheet.
      !get mass balance first from non-subgridded climate cells
      dmsk(:,:)=0.
      nsgsmb=0.
      sgsmb=0.
      if (isname=='Antarctica') then
        dmsk(:,1:50)=1.
      elseif (isname=='Greenland') then
        dmsk(:,51:jmt)=1.        
      endif      
      call areatot(nmbal(:,:,2),dmsk,nsgsmb)
      nsgsmb=nsgsmb/1000000.*s2i
      !then get mass balance from subgridded cells
      do n=1,numiE
        if (isniE(n)==isn) then
	  sgsmb=sgsmb+((mbalaccum(n,2)+mbalmelt(n,2))/100.*iariE(n)*s2i)
	endif
      enddo
      !accumulate in m3 ice
      climbudgsum=nsgsmb+sgsmb
      
      !Get volume of mass balance going into ice sheet (as seen in budgsnow - non SG) and compare against nmbal array.  Compare in units of m3.
      !Difference should be negligible.
!       print*, 'climbudgsum=',climbudgsum,'m3 ice'
!       print*, 'icebudgsum=', icebudgsum,'m3 ice'
!       print*, 'difference=',climbudgsum-icebudgsum,'m3 ice'
      !update previous ice sheet volume with mass balance from climate model
      !and calculate total mass in each sheet basin before dynamics and internal 
      !mass balance calculations.     
      do n=1,nisb
        if (isb(n) == isn) then
          temp(:,:) = 0.
          !mask out all other ice basins
          where(icedrainage.eq.n)temp=budgsnow !m ice
	  basinbudg=0.
	  !calculate basin budget.
          call mascalc(temp,basinbudg,tmp)
	  !calculate pre-ice-step total volume by basin: old basin volume + input mass balance budget
	  !print*, 'basinbudg,n=',basinbudg,n
          predynvol(n)=basinvol(n)+basinbudg
        endif
      enddo
      !Gather hemispheric-specific diagnostics output
      temp=0.
      where(mext.gt.0.)temp=1.
      if (isname=='Antarctica') then      
	sisaccum=accumsum
	sismelt=meltsum
        sissmb=accumsum+meltsum
	sismext=sum(temp(:,:)*darea(:,:))
      elseif (isname=='Greenland') then
	nisaccum=accumsum
	nismelt=meltsum
        nissmb=accumsum+meltsum
	nismext=sum(temp(:,:)*darea(:,:))	
      endif
      
      !Output some icesheet-specific diagnostics
      if (isn==strtis) then
        ntrec=ntrec+1
      endif
      
#if defined O_ism_get_standalone_sbc
      !get surface surface mass balance and surface ice temperatures for ice sheets
      !print to file, and stop program, at first run through the routine.  
      !This assumes the model is being run from a restart that is an acceptable endmember for a standalone glacial/interglacial run
      g_st=trim(isname)//'budgsnow_bc'
      call jerncdf_ice_timedep(budgsnow,g_st,g_st,ntrec)
      g_st=trim(isname)//'tsurfi_bc'
      call jerncdf_ice_timedep(tsurfi,g_st,g_st,ntrec)
      !if (isn==endis) stop
#endif
      
      return
      end
