      subroutine recalc_land_albedo

!======================================================================= 
      !recalculate land albedo (sbc(:,:,ilca)) and Dalton number
!======================================================================= 

      implicit none
      
      include "size.h"
      include "param.h"
      include "atm.h"
      include "csbc.h"
      include "veg.h"
  
      integer i,j,ii,jj,n  
      real z0
      do j=2,jmtm1
   	do i=2,imtm1	
# if defined O_mtlm
   	  if (land_map(i,j) .eq. 0 .and.
     &	      tmsk(i,j) .lt. 1.) then
# else
!          if (isurf(i,j) .gt. 0 .and. tmsk(i,j) .lt. 1.) 

          if (tmsk(i,j) .lt. 1.) then
# endif
            !find nearest full land cell, if partial land
	    if (tmsk(i,j) .gt. 0.) then
              call findclimate(tmsk,i,j,3,0.,ii,jj)	    
	    else
	      ii = i
	      jj = j
	    endif
# if defined O_crop_data
            sbc(i,j,iscal) = 1. - (crops(ii,jj,2)*veg_alb(icrops)
     &    		   + veg_alb(iveg(ii,jj))*(1.-crops(ii,jj,2)))
            z0 = veg_rl(iveg(ii,jj))*(1.-crops(ii,jj,2))
     &         + veg_rl(icrops)*crops(ii,jj,2)
!           surf_dalt(i,j) = 0.16/(log(10./z0)*log(10./(z0*exp(-2.))))
            veg_dalt(i,j) = 0.16/(log(10./z0)*log(73.89/z0))
          endif
        enddo
      enddo
# else
            sbc(i,j,iscal) = 1. - veg_alb(iveg(ii,jj))  
          endif
        enddo
      enddo

      do n=1,nveg
        z0 = veg_rl(n)
        veg_dalt(n) = 0.16/(log(10./z0)*log(73.89/z0))
      enddo
# endif	
      
      return
      end

!##################################################################################################################

      subroutine findclimate(arrayin,i,j,maxrad,crit,ii,jj)
      
!=======================================================================      
!     find indices of nearest cell to (i,j) that satisfies criteria crit.
!     cell at (i,j) is not evaluated.
!     search routines goes outward for maxrad cells until first cell
!     satisfies crit.  Study direction of loop for the specific spiral
!     path taken.  Indices of this cell (ii,jj) returned.  If search is
!     unsuccessful, program will stop (i.e. this search expects success).
!     Inputs:
!     -arrayin:		array to be evaluated (imt,jmt)
!     -i,j:		points to evaluate around
!     -maxrad:		maximum radius (units of grid cells)to search
!     -crit:		criteria to be satisfied
!     -ii,jj:		point within maxrad that satisfies crit
!======================================================================= 
 
      implicit none
      
      include "size.h"
      
      real arrayin(imt,jmt),crit
      integer i,j,rad,ii,jj,maxrad
      
      do rad = 1,maxrad
        do ii=i-rad, i+rad
	  do jj = j-rad, j+rad
	    if (arrayin(ii,jj) .eq. crit .and.
     &                          ii .gt. 1 .and.
     &                          ii .lt. imt .and.
     &				jj .gt. 1 .and.
     &				jj .lt. jmt) return
	  enddo
	enddo
      enddo
      print*,'i,j=',i,j
      stop '=> Error: not satisfied: subroutine findclimate'
      
      return
      end      

!##################################################################################################################

      subroutine seticebasin(isn)

!=======================================================================
!     Obtain Antarctic/Greenland drainage basins on ice grid
!=======================================================================

      use comicegrid
      use comicesheet, only: h,hw,icedrainage 
      
      implicit none

      character(120) :: fname, new_file_name

      logical exists, flag

      integer i, ios, iou, j, k, m, ib(10), ic(10)
      integer basin,isn

      include "size.h"
      include "param.h"
      include "atm.h"
      include "riv.h"
      include "ism.h"

      character(120) g_st
      real temp(imt,jmt)
      real tmpij(imt,jmt)

!     read the netcdf river file to obtain ice drainage basins.
      fname = new_file_name ("rivers.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Error => ", trim(fname), " does not exist."
        stop 'Stopped in icebasin'
      else
        ib(:) = 1
        ic(:) = imt
        ic(2) = jmt
        call openfile (fname, iou)
        call getvara ('rivers', iou, imt*jmt, ib, ic, tmpij
     &, 1., 0.)
        temp = tmpij
      endif
      call embmbc(temp)
      !Set initial drmap to river drainage map (including drainage over oceans).  This is for use in snow redistribution scheme as
      !a map for where to route stranded snow that is assumed melted.
      drmap(:,:)=temp(:,:)
      !Expand Antarctic basins northwards
      do j=13,51
        drmap(:,j)=drmap(:,12)
      enddo
      !Explicitly set all N Hemisphere to Greenland drainage basin
      drmap(:,52:jmt)=33
      
      !set very south pole basins manually
      do i=1,imt
       temp(i,1) = temp(i,3)
       temp(i,2) = temp(i,3)
      enddo

      if (isname == 'Antarctica') then
        !Downscale global river drainage basin map to ice sheet grid.
	!Recall, this map has been modified to include several Antarctic
	!drainage basins.
        icedrainage(:,:) = 0.0
        call downscale_to_ice_grid(temp,icedrainage,2)
      elseif (isname == 'Greenland') then
        !Send all runoff from Greenland to drainage basin 33, which is the
	!one over Greenland (look at land_map.nc)
	icedrainage(:,:) = 33
      endif    
      
      !Set integer array isb to ice sheet number to which basin belongs to, 
      !at indexes that correspond to overlying river-basin numbers.  
      !Recall, these numbers are ultimately set by read of ASCII rivermap file 
      !and so are somewhat arbitrary.

      do basin=1,nisb!search through all possible river basin numbers
        do i=1,nx
          do j=1,ny
            if (nint(icedrainage(i,j)).eq.basin) then 
              isb(basin)=isn
              goto 10 
            endif
          enddo
        enddo
   10   continue  
      enddo
      
      return
      end
      
!##################################################################################################################
      
      subroutine getbasinvol(isn)

!======================================================================= 
      !Gather initial ice sheet basin ice volumes and heat and moisture fluxes
!=======================================================================       

      use comicegrid
      use comicesheet, only: h, icedrainage
      implicit none

      include "size.h"
      include "ism.h"

      real temp(nx,ny)      
      integer n,isn
      real tmp

      do n=1,nisb
        if (isb(n) == isn) then	
          !calculate volumes over individual drainage areas, add to global bins
          !mask out all other ice basins
          !calculate updated basin volume after ice step.
          temp(:,:) = 0.	  
	  where(icedrainage.eq.n)temp=h
          call mascalc(temp,basinvol(n),tmp)
        endif
      enddo

      return
      end

!##################################################################################################################      

      subroutine fluxsort(runstep)
      
      use comicegrid, only: ice_accelerator
      use ice_sheet_storage, only: strtis,endis   
      
      implicit none
      
      include "size.h"
      include "ism.h"
      include "pconst.h"
      include "cembm.h"
      include "calendar.h" 
           
      integer n, nn
      logical runstep, exist
      isvol=0.
      sisvol=0.
      nisvol=0.
      disvol=0.      

      if (.not. runstep) then
        !Do tasks that are required only for the first step
        do n=1,nisb
        !zero all fluxes from non-ice-sheet grid basins (which could be imported from a restart which had more ice
	!sheet instances running than the current simulation).
	  exist=.false.
 	  do nn=strtis,endis
 	    if (isb(n) == nn) then   
 	      exist=.true.
 	    endif
 	  enddo
 	  if (.not. exist) then
	    dicevol(n)=0.
	  endif
	enddo
      
      else
        !Do tasks that are required only during a runstep
        isflxm(:)=0.
        isflxh(:)=0.	
        do n=1,nisb
          !if accelerated, climate sees one years worth of ice going in.
	  !So predynvol is one year's budgsnow, plus existing ice volume.
	  !But the ice can grow by more than the input that the climate sees.
	  !So if ice grows, in the absence of sufficiently increased discharge,
	  !basinvol (i.e. updated volume) will be greater than the original 
	  !volume, due to the multiple inputs of SMB.  Hence the negative flux
	  !of moisture.
	
	  !if running, get difference between total ice volume before ice dynamics, and volume after.  
	  !Add this difference to dicevol, which was already being accumulated from 
	  dicevol(n) = dicevol(n) + predynvol(n)-basinvol(n)
        if (dicevol(n).lt.0.) then
	!Surface mass balance code around ablation zone generates maximum potentially possible negative SMB around ablation zones, as elevation bins don't know about individual ice grid thicknesses (only average thickness).  Thus, in some cases subgridded nmbal loss is more than available ice. 
        !To calculate the flux of moisture to the ocean from calving, the routines first add the mass balance field to the ice sheet volume, to get a 'pre-dynamics' volume.  This volume incorporates the subgrid-derived SMB field, which in some cases has generated more melt than is available at particular grid cells.  It can therefore be more negative than the physically consistent SMB field.
        !Next, the total ice volume after the dynamics code is calculated.  However, within the ice dynamics code, ice loss is limited to available ice (to be physically consistent).  If an SMB 'overshoot' cell occurred in the first sum, when the two are compared the ice sheet will have appeared to have gained some artificial volume.
        !As the ice-ocean flux is determined by the pre- and post- dynamic volumes, any overshoot in the absence of significant calving will result in a flux of moisture from ocean, to compensate for artificially large SMB loss (i.e. overshot ablation results in removal of ocean water). 
	!To visualize where this occurs, uncomment the totneg_land output in icedyn.F.
	!A similar issue occurs if ice is accelerated, since the ice can grow by multiples of the input that the climate sends down, if the dV/dt is positive and acceleration is turned on.
 	  print*, 'Note: net negative SMB overshoot in basin #',n
	  print*, 'Ocean flux compensation occurring.'
	endif

        diceheat(n)= -dicevol(n)*1.e6*rhoice*flice*1.e-7
        isflxm(n)= dicevol(n)*1.e6*rhoice/yrlen*secday
     &  	      /ice_accelerator
        isflxh(n)= diceheat(n)/yrlen*secday
     &  	      /ice_accelerator  	      

        enddo
      
      endif
      
      !accumulate sums for global sums, regardless of whether this is initial call, or runstep call
      do n=1,nisb
      isvol=isvol+basinvol(n) 
	if (isb(n) == 1) then
	  sisvol = sisvol + basinvol(n)
	elseif (isb(n) == 2) then 
	  nisvol = nisvol + basinvol(n)
	endif	
	disvol=disvol+dicevol(n)
      enddo   

      return
      end

!##################################################################################################################
      
      subroutine elevbin(isn)
      
!======================================================================= 
      !Generate ice sheet arrays with reference to parent global cell
      !indexes iceiind, icejind
      
      !Allocate global arrays with reference to cell-specific 
      !starting points and number of elevation bins (i.e. count) 
      !in subgrid arrays
      
      !Allocate 1D subgrid arrays       
!=======================================================================             

      use comicegrid, only: 
     &nx,ny,darea,isname,iceiind,icejind,is2sg,ismrflag
      use comicesheet, only: hs,h,hw
      use ice_sheet_storage, only: strtis,endis
      use subgrid
       
      implicit none
      
      include "size.h"        !imt,jmt
      include "param.h"       !imtm1,jmtm1
      include "ism.h"         !ism variables   
      include "atm.h"         !tmsk, subgrid information arrays 
      include "switch.h"      !init
      
      integer i,j,n,ii,jj,nn,kmax,isn,count,jstrt,jend
      integer nsnew,nenew,nsold,neold    
      logical found,init_call
      integer nindex(nx,ny),iind(imt,jmt),jind(imt,jmt)
      real eit,icesum,tmp,temp(nx,ny)
 
      character(120) g_st
!-----------------------------------------------------------------------
      !Do non-ice-sheet-specific tasks
!-----------------------------------------------------------------------

      logical firstcall
      data firstcall /.true./
      save firstcall
       
      init_call=.false.	   
      if (init .and. firstcall) init_call=.true.	      
      kmax=50
      if (isn == strtis) then           
        !Zero goodies.  	
	!Initialize some temporary arrays to be used below
        allocate(ibin(imt,jmt,kmax))	
        ibin(:,:,:) = 0.
        numiE=0
        !As long as not first timestep of initial-condition run, set arrays to incoming strtiE and endiE.
	if (.not. init_call) then
           strtiEsav(:,:)=strtiE(:,:)
           endiEsav(:,:)=endiE(:,:)
 	endif
	!Zero index arrays in prep for filling.		
	strtiE(:,:)=1
	endiE(:,:)=0			
      endif  

      !Downscale imt and jmt arrays to ice grid so each
      !land-based ice cell has a land-or-ice shelf (imt,jmt) associated with it.
      do i=1,imt
        do j=1,jmt
	  iind(i,j) = i
	  jind(i,j) = j
	enddo
      enddo		     
     		     
      call downscale_to_ice_grid(real(iind),temp,2)
      iceiind(:,:)=nint(temp(:,:))
      call downscale_to_ice_grid(real(jind),temp,2)
      icejind(:,:)=nint(temp(:,:))
   
      do i=1,nx
        do j=1,ny
	  if (h(i,j) .gt. himin) then 
            found = .false.
	    n=1
	    !set minimum ice elevation to 0m
	    tmp=max(0.,hs(i,j)) 
            do while (.not. found)
              eit=(n-1)*dzi	    
	      !if ice surface elevation is .le. the test elevation
	      !value add the area of this ice sheet grid cell to the
	      !'parent' climate model subgrid area bin, accumulated 
	      !climate model total area array.  Also add area-weighted
	      !elevation to average elevation bin.
              if (tmp.le.eit) then
                ibin(iceiind(i,j),icejind(i,j),n) = 1
		nindex(i,j)=n
		found = .true.		  
              else
                n=n+1
		if (n .gt. kmax) 
     &		stop 'Error: kmax exceeded: elevbin.'
              endif	
            enddo
	  endif
        enddo
      enddo
      
      if (isn==1) then !AIS
        jstrt=2
        jend=51
      elseif (isn==2) then !GIS
        jstrt=52
        jend=101
      endif
      !Determine length (numiE) of ice sheet subgrid array,
      !assign a subgrid array number to ice sheet cells,
      !and determine starting and ending indices for global grid
      !subgrid elevation bins.
      is2sg(:,:)=0.
      do i=2,imtm1
        do j=jstrt,jend
          !First, save old indices	  	
          found=.false.
          do n=1,kmax
            if (ibin(i,j,n) .ne. 0.) then
	      !bump numiE up by one if ice elevation is present.
              numiE=numiE+1
	      !find any ice sheets cells that correspond to this
	      !subgrid bin, and let them know index (put in is2sg).   
              do ii=1,nx
        	do jj=1,ny
        	  if (nindex(ii,jj) .eq. n
     &  	    .and.iceiind(ii,jj) .eq. i
     &  	    .and. icejind(ii,jj) .eq. j) then
        	    is2sg(ii,jj)=numiE
        	  endif
        	enddo
              enddo
	      !assign strtiE if first ice elevation cell found for
	      !climate grid cell.	
              if (.not. found) then	  
        	strtiE(i,j)=numiE
        	found=.true.
              endif
	      !keep bumping endiE up as new cells are found
              endiE(i,j)=numiE  	      
            endif	     
          enddo
	  !if no ice cells found, then set endiE to be less than strtiE
	  !in hopes that this will keep ice loops from being entered at
	  !all.
          if (.not.found) then
            strtiE(i,j)=1
            endiE(i,j)=0
          endif 	    
        enddo
      enddo     
      !write restart here instead of in icectl in order to include updated indexing
      !arrays which were generated above.
      if (ismrflag) then
        !set ice sheet ID tag to a random number which all ice sheet and climate restarts will be tagged with
	if (isn==strtis) rtag=sum(h+hw)+numiE
        call writeres
      endif      
!-----------------------------------------------------------------------
      !Now do final tasks.
!-----------------------------------------------------------------------

      if (isn==endis) then    
        !If climate model starting from initial conditions, then set indices
	!to those just determined
        if (init_call) then
	  strtiEsav(:,:) = strtiE(:,:)
	  endiEsav(:,:) = endiE(:,:)
	endif
        !set firstcall to false, so init_call will also be set to false next
        !swipe through.
	firstcall=.false.
        !(Re)-allocate subgrid arrays over ice sheet that do not need to be
	!redistributed or are not used in the redistribution process.  Initialize
	!them too.
	if (allocated(eleviE)) then
	  deallocate(eleviE)
	endif	
	allocate(eleviE(numiE))	
	eleviE(:)=0.
	
	if (allocated(precipiE)) then
	  deallocate(precipiE)
	endif	
	allocate(precipiE(numiE))	
	precipiE(:)=0.			
	
	if (.not. allocated(ticeiE)) then
	  allocate(ticeiE(numiE))
	endif		
	allocate(ticeiEsav(size(ticeiE)))	
	ticeiEsav=ticeiE	
	deallocate(ticeiE)	
	allocate(ticeiE(numiE))
	ticeiE(:)=0.        
	
	if (allocated(ticeiE)) then
	  deallocate(ticeiE)
	endif	
	allocate(ticeiE(numiE))
 	ticeiE(:) = 0.
	
	if (allocated(ismticeiE)) then
	  deallocate(ismticeiE)
	endif	
	allocate(ismticeiE(numiE))			
	ismticeiE(:)=0.	
		
	if (allocated(ithkiE)) then
	  deallocate(ithkiE)
	endif	
	allocate(ithkiE(numiE))			
	ithkiE(:)=0.	
	 
	if (.not. allocated(fraciE)) then
	  allocate(fraciE(numiE))
	endif		
	allocate(fraciEsav(size(fraciE)))	
	fraciEsav=fraciE	
	deallocate(fraciE)	
	allocate(fraciE(numiE))
	fraciE(:)=0.	
	if (.not. allocated(mbalmelt)) then
	  allocate(mbalmelt(numiE,2))
	endif	
	deallocate(mbalmelt)	
	allocate(mbalmelt(numiE,2))
	mbalmelt(:,:)=0.	
	
	if (.not. allocated(mbalaccum)) then
	  allocate(mbalaccum(numiE,2))
	endif
	deallocate(mbalaccum)	
	allocate(mbalaccum(numiE,2))
	mbalaccum(:,:)=0.		

	if (.not. allocated(hsnowiE)) then
	  allocate(hsnowiE(numiE,2))
	endif	
	allocate(hsnowiEsav(size(hsnowiE(:,2))))	
	hsnowiEsav=hsnowiE(:,2)	
	deallocate(hsnowiE)	
	allocate(hsnowiE(numiE,2))
	hsnowiE(:,:)=0. 
	
	if (.not. allocated(iariE)) then
	  allocate(iariE(numiE))
	endif	
	allocate(iariEsav(size(iariE(:))))	
	iariEsav=iariE(:)	
	deallocate(iariE)	
	allocate(iariE(numiE))
	iariE(:)=0.
	
	if (allocated(albedoiE)) then
	  deallocate(albedoiE)
	endif	
	allocate(albedoiE(numiE))	
	albedoiE(:)=0.		
	
	if (allocated(mextiE)) then
	  deallocate(mextiE)
	endif	
	allocate(mextiE(numiE))	
	mextiE(:)=0.	
	
	if (allocated(mduriE)) then
	  deallocate(mduriE)
	endif	
	allocate(mduriE(numiE))	
	mduriE(:)=0.	

	if (allocated(satiE)) then
	  deallocate(satiE)
	endif	
	allocate(satiE(numiE))	
	satiE(:)=0.
	
	if (allocated(isniE)) then
	  deallocate(isniE)
	endif	
	allocate(isniE(numiE))	
	isniE(:)=0.

	
#if defined O_refreeze
	if (.not. allocated(hnewsno)) then
	  allocate(hnewsno(numiE))
	endif	
	allocate(hnewsnosav(size(hnewsno(:))))	
	hnewsnosav=hnewsno(:)	
	deallocate(hnewsno)	
	allocate(hnewsno(numiE))
	hnewsno(:)=0.
	
	if (.not. allocated(counter)) then
	  allocate(counter(numiE))
	endif	
	allocate(countersav(size(counter(:))))	
	countersav=counter(:)	
	deallocate(counter)	
	allocate(counter(numiE))
	counter(:)=0.	
	
	if (.not. allocated(htemp)) then
	  allocate(htemp(numiE))
	endif	
	allocate(htempsav(size(htemp(:))))	
	htempsav=htemp(:)	
	deallocate(htemp)	
	allocate(htemp(numiE))
	htemp(:)=0.
	
	if (.not. allocated(hreftot)) then
	  allocate(hreftot(numiE))
	endif	
	allocate(hreftotsav(size(hreftot(:))))	
	hreftotsav=hreftot(:)	
	deallocate(hreftot)	
	allocate(hreftot(numiE))
	hreftot(:)=0.	
	
	if (.not. allocated(meltprev)) then
	  allocate(meltprev(numiE))
	endif	
	allocate(meltprevsav(size(meltprev(:))))
	meltprevsav=meltprev(:)	
	deallocate(meltprev)	
	allocate(meltprev(numiE))
	meltprev(:)=0.     
#endif
#if defined O_topography_shading
	if (allocated(snv)) then
	  deallocate(snv)
	endif	
	allocate(snv(numiE,3))	
	snv(:,1:2)=0.   !set initial normal vectors to a horizontal surface 
	snv(:,3)=1.  
	
	if (allocated(zen_h)) then
	  deallocate(zen_h)
	endif	
	allocate(zen_h(numiE,2))
	zen_h(:,:)=0.  				
#endif

#if defined O_subgrid_out
	if (allocated(dswr_out)) then
	  deallocate(dswr_out)
	endif	
	allocate(dswr_out(numiE))	
	dswr_out(:)=0.
	
	if (allocated(ultnt_out)) then
	  deallocate(ultnt_out)
	endif	
	allocate(ultnt_out(numiE))	
	ultnt_out(:)=0.
	
	if (allocated(usens_out)) then
	  deallocate(usens_out)
	endif	
	allocate(usens_out(numiE))	
	usens_out(:)=0.
	
	if (allocated(ulwr_out)) then
	  deallocate(ulwr_out)
	endif	
	allocate(ulwr_out(numiE))	
	ulwr_out(:)=0.	
	
	if (allocated(pdd_fac)) then
	  deallocate(pdd_fac)
	endif	
	allocate(pdd_fac(numiE))	
	pdd_fac(:)=0.		
#endif
		
	deallocate(ibin)

      endif      
      return
      end  

!##################################################################################################################

      subroutine elevsort(isn)
      
!======================================================================= 
      !Fill 1D ice grid fraction and elevation arrays with values
      !Reduce non-ice-sheet subgrid areal fractions by ice sheet overlay
      !Redistribute surface accumulations (e.g. snow) to new elevation bins   
!======================================================================= 

      use comicegrid, only: nx,ny,darea,is2sg,iceiind,icejind,isname
     &,alond,pi
      use comicesheet, only: hs,h,hw,hb,dtimeice
      use ice_sheet_storage, only: strtis,endis    
      use subgrid
      
      implicit none
      
      include "size.h"        !imt,jmt
      include "param.h"       !imtm1,jmtm1
      include "ism.h"         !ism variables   
      include "atm.h"         !tmsk, subgrid information arrays   
      include "levind.h"      !kmt    
      include "ice.h"	      !hsnol 
      include "grdvar.h"      !cst,dyt,dxt
      include "cembm.h"       !rhosno, rhoice, hsno_max
      include "switch.h"      !init

      real totar(imt,jmt), car(imt,jmt)           
      integer i,j,n,isn,nn,nenew,nsnew,neold,nsold
      integer jstrt,jend,ntrec
      character(120) g_st
      real temp(imt,jmt),vol,tmp,fx
      real tempi(nx,ny)
      logical firstcall
      data firstcall /.true./
      data ntrec /0/
      real, allocatable :: tmpf(:),snovol(:),snnew(:),binar(:)
      real, allocatable :: frisav2(:)
      save totar,firstcall,ntrec
      real nnew,nold
#if defined O_topography_shading
      real c(nx,ny,3),zcx(nx,ny),zcy(nx,ny),lon
      real radians_from_noon_sun
#endif      
      do j=1,jmt
        fx = cst(j)*dyt(j)
	do i=1,imt
	  car(i,j)=dxt(i)*fx
	enddo
      enddo
      
      if (isn==strtis)then
	fracl(:,:) = 0.
	fracbl(:,:) = 1.-tmsk(:,:)
	!fracbl(:,:) = 0.
	fraci(:,:) = 0.	
	fracblis(:,:) = 0.	
	fraco(:,:) = 0.
	totar(:,:) = 0.
	icegrid(:,:)=0
	icemsk(:,:)= 0
	cir(:,:)=0.
	aelev(:,:)=0.
	avgie(:,:)=0.
      endif
#if defined O_topography_shading  
      !get surface slope normal vector components
      call slope_calculator(c)    
      !get zenith angle horizontal components, on stereographic grid (assume sun comes from straight N-S direction)
      do i=1,nx
        do j=1,ny
	  lon=alond(i,j)*pi/180.
	  !two cases, one for AIS and on for GIS?  Case below is for GIS so far, not sure if it's generalizable.
          zcx(i,j)=-sin(lon)
          zcy(i,j)=-cos(lon)
	enddo
      enddo
      do i=1,nx
        do j=1,ny
	  n=is2sg(i,j)	
	  if (h(i,j).gt.himin) then
	    !accumulate normal components of surface slope vector
	    snv(n,1)=snv(n,1)+c(i,j,1)*darea(i,j)
	    snv(n,2)=snv(n,2)+c(i,j,2)*darea(i,j)
	    snv(n,3)=snv(n,3)+c(i,j,3)*darea(i,j)
	    !accumulate horizontal components of zenith vector
	    zen_h(n,1)=zen_h(n,1)+zcx(i,j)*darea(i,j)
	    zen_h(n,2)=zen_h(n,2)+zcy(i,j)*darea(i,j)
	  endif      
	enddo
      enddo      
#endif
      do i=1,nx
        do j=1,ny
	  n=is2sg(i,j)	
	  !accumulate bare land area in climate grid cells	  
	  if (hw(i,j) .le. 0. .and. h(i,j) .le. himin) then 	
	    fracblis(iceiind(i,j),icejind(i,j)) = 
     & 	    fracblis(iceiind(i,j),icejind(i,j))+darea(i,j)       
	  endif
	  if (hw(i,j) .gt. 0. .and. h(i,j) .le. himin) then
	    fraco(iceiind(i,j),icejind(i,j)) = 
     &	    fraco(iceiind(i,j),icejind(i,j))+darea(i,j)
          endif  
	  !accumulate ice-covered area in climate grid cells
          if (h(i,j) .gt. himin) then
	    fraci(iceiind(i,j),icejind(i,j)) =          
     &	    fraci(iceiind(i,j),icejind(i,j))+darea(i,j)
	    !build up total area of ice in fraciE array.  This will be divided by total area later
	    !to give fractional area.
            fraciE(n)=fraciE(n)+darea(i,j)
	    !build up elevation*area in subgrid bin
            eleviE(n)=eleviE(n)+(max(0.,hs(i,j)))*darea(i,j)
	    !build up ice thickness*area in subgrid bin
	    ithkiE(n)=ithkiE(n)+(h(i,j)*darea(i,j))
	    !accumulate ice-covered area in subgrid cell
	    iariE(n)=iariE(n)+darea(i,j)	    
	    !set ice sheet identifier array   
	    isniE(n)=isn   
          endif
	  !accumulate total land/ice area in climate grid cells
	  if (hw(i,j) .le. 0. .or. h(i,j) .gt. himin) then
	    fracl(iceiind(i,j),icejind(i,j)) =
     &	    fracl(iceiind(i,j),icejind(i,j))+darea(i,j)	    
          endif
	  !accumulate total ice sheet grid area in climate grid cells.  Might be sketchy using this... because
	  !around edges of ice sheet grid, overlaid climate grids will only be partially covered with ice sheet
	  !grids.
	  totar(iceiind(i,j),icejind(i,j)) =	  
     &	  totar(iceiind(i,j),icejind(i,j))+darea(i,j)
     
          !accumulate area*ice bed elevation value.
	  aelev(iceiind(i,j),icejind(i,j)) =
     &	  aelev(iceiind(i,j),icejind(i,j))+
     &    darea(i,j)*hb(i,j)
	enddo
      enddo
      
      !Get total ice sheet areas for output.
      tempi(:,:)=0.
      if (isn==1) then
        where (h.gt.himin) tempi=1
	sisiar=sum(darea(:,:)*tempi(:,:))
      elseif (isn==2) then
        where (h.gt.himin) tempi=1      
	nisiar=sum(darea(:,:)*tempi(:,:))        
      endif
      !Once all ice sheets have contributed to accumulations, divide accumulations to get fractional coverages.
      if (isn==1) then
        jstrt=2
        jend=51
      elseif (isn==2) then
        jstrt=52
        jend=101
      endif   
      do j=jstrt,jend            
        do i=2,imtm1
	  if (totar(i,j) .gt. 0.) then	
	    !generate ice grid coverage map
	    icegrid(i,j)=1    
	    !get the fraction of ice sheet cells within parent climate cell that are bare ocean.
	    fraco(i,j) = fraco(i,j)/totar(i,j)
	    !get the fraction of ice sheet cells within parent climate cell that are either ice sheet, shelf, or bare land.
	    fracl(i,j) = fracl(i,j)/totar(i,j)
	    !get the fraction of ice sheet cells within parent climate cell that are ice.
	    fraci(i,j) = fraci(i,j)/totar(i,j) 
	    !get the fraction of ice sheet cells within parent climate cell that are bare land
	    fracblis(i,j) = fracblis(i,j)/totar(i,j)    
	    do n=strtiE(i,j),endiE(i,j)
#if defined O_topography_shading 	      
	      !calculate average normal components of surface slope
	      snv(n,1)=snv(n,1)/fraciE(n)
	      snv(n,2)=snv(n,2)/fraciE(n)	
	      snv(n,3)=snv(n,3)/fraciE(n)
	      zen_h(n,1)=zen_h(n,1)/fraciE(n)
	      zen_h(n,2)=zen_h(n,2)/fraciE(n)      
#endif		    
	      !divide accumulated elevation*area field by area to get area-averaged elevation in each bin    
              eleviE(n)=eleviE(n)/fraciE(n)*100.
	      !divide thickness*area field by area to get area-averaged ice thickness in each bin
	      ithkiE(n)=ithkiE(n)/fraciE(n)		      
	      !When done using fraciE as a divisor, calculate fraciE as the fractional coverage of subgrid bin as a % of the total ice sheet cell area.
              fraciE(n)=fraciE(n)/car(i,j)*10000.                    
	      avgie(i,j)=avgie(i,j)+eleviE(n)*(iariE(n)/totar(i,j)) 
	    enddo 
	    !reset non-subgridded elevation to exposed bedrock, where ice sheet model occurs.
	    elev(i,j)=max(0.,aelev(i,j)/totar(i,j))*100.
	    !get average ice sheet elevation in a grid box
	  endif
	  !Bootfuck some shit.  	
	  if (kmt(i,j) .lt. 0.5) then
	    fracbl(i,j)=1.-fraci(i,j)!-fraco(i,j)
	    if (fracbl(i,j) .lt. 1.e-16) fracbl(i,j) = 0.
	  else
	    fracbl(i,j)=0.
	    fracblis(i,j)=0.
	  endif 
	  !cir=total dry bare land area fraction (calculated as total dry bare ice cells/total ice cells),
	  !SCALED BY total ice area/climate cell area ratio.
	  cir(i,j)=fracblis(i,j)*totar(i,j)/car(i,j)*10000. 	
	enddo
      enddo
      
      !recalculate tmsk.
      tmsk(2:imtm1,jstrt:jend)=0.
      do i=2,imtm1
        do j=jstrt,jend
	  if (kmt(i,j) .gt. 0.) tmsk(i,j)=1.-fraci(i,j)
	  !icemsk(i,j)=icegrid(i,j)
	  icemsk(i,j) = ceiling(fracl(i,j))
	enddo
      enddo	      
      if (isn==endis) then
        call embmbc(elev)
	call embmbc(tmsk)
	call embmbc(fraci)
	call embmbc(fracl)
	call embmbc(fracbl)
	call embmbc(icemsk)
	call embmbc(totar)		
	iac=0
#if defined O_ism_pdd     
        year_accum(:) = 0.
        gtempmon(:,:,:) = 0.0
        gprecmon(:,:,:) = 0.0 
#endif	
        !update ice time for climate   
        dtism=dtimeice                    				
        !If first timestep of initial-condition run, set values to those just determined,
	!and set initial snow thickness to maximum
        if(init .and. firstcall) then
	  fraciEsav(:)=fraciE(:)
	  iariEsav(:)=iariE(:)
	  hsnowiEsav(:)=hsno_max
	  where(icemsk.gt.0) hsnol(:,:,2)=hsno_max*fracbl(:,:)
	  nmbal(:,:,:)=0.
	endif	
        firstcall=.false.
	
        !redistribute surface values due to shifting ice sheet elevation bins.
        !cases: ice retreats, leaving bare land; ice advances, taking over bare land elevation area;
        !ice changes elevation, taking over other ice cells
        allocate(frisav2(size(fraciEsav)))
        frisav2(:)=fraciEsav(:)
	temp(:,:) = 0.
        do i=2,imtm1
          do j=2,jmtm1
	    nsnew=strtiE(i,j)
 	    nenew=endiE(i,j)
 	    nsold=strtiEsav(i,j)
 	    neold=endiEsav(i,j) 	       
 	    if (nsnew.le.nenew) then  	    	     
 	      !allocate short working arrays for convenience
	      allocate(tmpf(nenew-nsnew+1))
 	      allocate(snovol(nenew-nsnew+1))
 	      allocate(snnew(nenew-nsnew+1))
	      allocate(binar(nenew-nsnew+1))    
 	      tmpf(:)=fraciE(nsnew:nenew)
 	      snovol(:)=0.     
 	      !step down from highest updated elevation
 	      do n=size(tmpf),1,-1	      
 	         !step down old elevations
 	        do nn=neold,nsold,-1
 	          !set tmp to the minimum of the old ice area, and
 	          !the remaining area that the updated ice bin needs
 	          !to fill.
 	          tmp=(min(fraciEsav(nn),tmpf(n)))
 	          !reduce total ice area of grid cell that needs to be filled
 	          tmpf(n)=tmpf(n)-tmp
 	          !reduce area of old bin by amount stolen by new cell.
 	          fraciEsav(nn)=fraciEsav(nn)-tmp
 	          !accumulate a snow volume equivalent to the amount stolen
 	          snovol(n)=snovol(n)+tmp*hsnowiEsav(nn)
 	        enddo
 	      enddo
	      !divide by fractional area of new subgrid cells to get real snow thickness 
              hsnowiE(nsnew:nenew,2)=snovol(:)/fraciE(nsnew:nenew) 
  	      !Some updated ice area may remain after above iteration, indicating that
  	      !that the ice has expanded over some bare land fraction.  Deal with this case.
	      tmp=sum(tmpf)
 	      snnew(:)=0.
 	      if (tmp.gt.1.e-15) then
	        !print*, 'Ice expanded here:', i,j
 	        if (tmp .gt. 1.) then
 		  print*, 'ice frac .gt. 1.',tmp,i,j
 		endif
 	        !Get the total true volume of snow (cm3) from the bareland that covers this fractional area.
 	        vol=hsnol(i,j,2)*car(i,j)*(min(tmp,1.))
                binar(:)=iariE(nsnew:nenew)*10000. !m2 to cm2
	        do n=1,size(tmpf)
          	  !Get fraction of real volume that will belong to each subgrid cell.  If area of subgrid cell
	          !is zero, it doesn't get any new volume.  Divide volume by subgrid area (binar).
		  snnew(n)=vol* tmpf(n)/tmp /binar(n)	  
		  !hsnowiE(nsnew+n-1,2)=hsnowiE(nsnew+n-1,2)+snnew(n)
	        enddo
		!Add this snow thickness to existing subgrid snowpack
		hsnowiE(nsnew:nenew,2)=hsnowiE(nsnew:nenew,2)+snnew(:)
  	        !Then remove this volume of snow from the bareland bin
  	        hsnol(i,j,2)=hsnol(i,j,2)-vol/car(i,j)   
 	      endif   			
 	      !Deallocate temporary arrays for next pass through.
 	      deallocate(tmpf,snovol,snnew,binar)
	      
	      !sum up new avgsno for global sums.
	      avgsno(i,j)=0.	    
	      do n=strtiE(i,j),endiE(i,j)
	        avgsno(i,j)=avgsno(i,j)+(hsnowiE(n,2)*fraciE(n))     
	      enddo 	      
	      
	      !if old subgrid exists
	      if (neold .gt. 0.) then
	        do n=0,(nenew-nsnew)
	          !calculate current index of updated arra
                  nnew=nsnew+n
		  !calculate index of old array (set to last subgrid array if new subgrid is larger than old subgrid
		  nold=min((nsold+n),neold)
		  !note: if new subgrid is smaller than old subgrid, old (highest-elev) grid values will simply fall off the map.
		  
	          !Redistribute ice temperatures onto new subgrid bins, as well as possible 
	          !(only used as an initial guess for EMBM, so not necessary to get perfect).  		  
		  ticeiE(nnew)=ticeiEsav(nold)
#if defined O_refreeze		  
	          !Redistribute refrozen fields
		  hnewsno(nnew) =hnewsnosav(nold)
		  counter(nnew) =countersav(nold)
		  htemp(nnew)   =htempsav(nold)
		  hreftot(nnew) =hreftotsav(nold)
		  meltprev(nnew)=meltprevsav(nold)
#endif		  
	        enddo
	      else
	        !these are actually redundant, as during allocation all arrays are set to 0 anyways.
	        ticeiE(nsnew:nenew)=0.
#if defined O_refreeze		
	        hnewsno(nsnew:nenew)=0.
		counter(nsnew:nenew)=0.
		htemp(nsnew:nenew)=0.
		hreftot(nsnew:nenew)=0.
		meltprev(nsnew:nenew)=0. !0=false
#endif			
	      endif      	                 
	    endif
	    
	    !Some old ice area may still exist, indicating the ice has shrunk, giving up
	    !some fraction to bare land.  Deal with this case.
            vol=0.
 	    do n=nsold,neold
	      !If old remaining ice fraction remains, then...
 	      if (fraciEsav(n).gt.1.e-15) then
	        !print*, 'Ice retreated here:',i,j
	      	vol=vol
	        !Get volume of snow in old bin in m3 snow
     &      	+iariEsav(n)*hsnowiEsav(n)/100.
            	!scaled by the area that has been removed
     &      		    *fraciEsav(n)/frisav2(n)
	      endif
 	    enddo
	    !Add this volume to the appropriate ice sheet runoff bin. 
	    dicevol(drmap(i,j))=
     &	    dicevol(drmap(i,j))+vol*rhosno/rhoice
	  
	    !recalculate updated avgsno array for global sums.  This is not necessary for avgmbal because it is 
	    !zeroed after every ice sheet step, but done for redundancy.
	    avgsno(i,j)=0.
	    avgmbal(i,j)=0.   
            do n=strtiE(i,j),endiE(i,j)
              avgsno(i,j)=avgsno(i,j)+hsnowiE(n,2)*fraciE(n)
              avgmbal(i,j)=avgmbal(i,j)+
     &	      (mbalmelt(n,2)+mbalaccum(n,2))*fraciE(n)
            enddo	    
 	  enddo
        enddo			
			
      endif          
      if (isn == endis) then
        deallocate(fraciEsav
     &	          ,hsnowiEsav
     &	          ,iariEsav
     &    	  ,ticeiEsav)
	deallocate(frisav2)
#if defined O_refreeze	
        deallocate(hnewsnosav
     &	          ,countersav
     &	          ,htempsav
     &	          ,hreftotsav
     &	          ,meltprevsav)
#endif
	
! 	g_st='fracl'
! 	call jerncdf_snapshot(real(fracl(:,:)),g_st)
!  	g_st='fracbl'
!  	call jerncdf_snapshot(real(fracbl(:,:)),g_st)
!  	g_st='fraci'
!  	call jerncdf_snapshot(real(fraci(:,:)),g_st)
	!stop
! 	g_st='fracblis'
! 	call jerncdf_snapshot(real(fracblis(:,:)),g_st)	
! 	g_st='fraco'
! 	call jerncdf_snapshot(real(fraco(:,:)),g_st)
! 	g_st='totar'
! 	call jerncdf_snapshot(real(totar(:,:)),g_st)
!  	g_st='icegrid'
!  	call jerncdf_snapshot(real(icegrid(:,:)),g_st)
!  	g_st='icemsk'
!  	call jerncdf_snapshot(real(icemsk(:,:)),g_st)
! 	g_st='cir'
! 	call jerncdf_snapshot(real(cir(:,:)),g_st)
! 	g_st='aelev'
! 	call jerncdf_snapshot(real(aelev(:,:)),g_st)
! 	g_st='avgie'
! 	call jerncdf_snapshot(real(avgie(:,:)),g_st)
! 	g_st='tmsk'
!         call jerncdf_snapshot(real(tmsk(:,:)),g_st)
      endif
      return
      end  
