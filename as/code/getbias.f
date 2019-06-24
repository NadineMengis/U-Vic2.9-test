! source file: /net/mare/home1/eby/as/ism/getbias.F
      subroutine getbias(isn)

      use ice_sheet_storage, only: strtis,endis
      use subgrid, only: biaiE,totariE,numiE,sbciE,isniE,himin
      use comicegrid, only: nx,ny,is2sg,darea,isname
      use comicesheet, only: h

      implicit none

      include "size.h"
      include "param.h"
      include "pconst.h"

      character(120) fname,name,new_file_name
      integer i,j,isn,n,mon,iou,ib(10),ic(10)
      real cb(imt,jmt),tmpij(imtm2,jmtm2)
      real iceb(nx,ny), maxbias, minbias
      logical exists,inqvardef

      character(120) g_st
      if (isn==strtis) then
        if (allocated(biaiE)) deallocate (biaiE)
        allocate (biaiE(numiE,12))
        if (allocated(sbciE)) deallocate (sbciE)
        allocate(sbciE(numiE))
        if (allocated(totariE)) deallocate (totariE)
        allocate(totariE(numiE))
        biaiE(:,:)=0.
        totariE(:)=0.
        sbciE(:)=0.
      endif

      !load global monthly SAT bias maps
      fname = new_file_name("uvic_bias.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Error => ", trim(fname), " does not exist."
        stop 'getbias.F'
      endif
      do i=1,nx
        do j=1,ny
          if (h(i,j).gt.himin) then
            n=is2sg(i,j)
            totariE(n)=totariE(n)+darea(i,j)
          endif
        enddo
      enddo
      ib(:) = 1
      ic(:) = 1
      ic(1) = imtm2
      ic(2) = jmtm2
      call openfile (fname, iou)
      name = 'uvic_sat_bias'
      if (inqvardef(name, iou)) then
	do mon=1,12
	  ib(3) = mon
	  cb(:,:)=0.
	  call getvara (name,iou,imtm2*jmtm2,ib,ic,tmpij,c1,c0)
	  cb(2:imtm1,2:jmtm1) = tmpij
	  call embmbc (cb)
	  call downscale_to_ice_grid(cb,iceb,1)
	  do j=1,ny
	    do i=1,nx
	      if (h(i,j).gt.himin) then
	        n=is2sg(i,j)
	        biaiE(n,mon)=biaiE(n,mon)+iceb(i,j)*darea(i,j)
	      endif
	    enddo
	  enddo
        enddo
        !call closefile (iou)
        else
	  print*, 'Error in getbias: file not found:',fname
        stop
      endif

      if (isn==endis) then
        do mon=1,12
          maxbias=0.
	  minbias=0.
          do n=1,numiE
            if (totariE(n).gt.0) biaiE(n,mon)=biaiE(n,mon)/totariE(n)
	    !try to catch potential bias error
	    if (isniE(n) .eq. isn) then
	      maxbias=max(maxbias,biaiE(n,mon))
	      minbias=min(minbias,biaiE(n,mon))
	    endif
          enddo

	  if (maxbias==0. .and. minbias==0) then
            print*, 'WARNING!!!! CHECK BIASMAP HERE!!!!'
	    print*, 'over ice sheet',isn,', all biases'
	    print*, 'were set to zero in getbias.'
	    print*, 'numiE=',numiE
	  endif
        enddo

      endif
      return
      end

      subroutine calcbias

      use subgrid

      implicit none
      include "tmngr.h"

      integer monprev, monnext
      real modtime,wtprev,wtnext

      !generate subgridded biases from monthly bias averages.  Averages taken to be valid at midpoint of month
      modtime=dayoyr/365.*12.
      monnext=min(12,max(1,ceiling(modtime+0.5)))
      monprev=monnext-1
      if (monprev==0) monprev=12
      !weight to previous month
      wtnext=1.-((real(monnext)-0.5)-modtime)
      !weight to next month
      wtprev=1.-wtnext
      sbciE(:)=biaiE(:,monprev)*wtprev + biaiE(:,monnext)*wtnext

      return
      end
