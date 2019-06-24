      subroutine jerncdf_snapshot(field,fname)
	
      !---------------------------------------------------------------
      !Prints imt x jmt array to named netcdf.  Will overwrite.	

      include "size.h"	
      include "coord.h" 	
      include "pconst.h" 	
      
      integer ntrec, it(10), ib(10), ic(10)
      character(*) :: fname
      integer id_xt, id_yt, iou, n
      real field(1:imt,1:jmt)
      real tmpij(2:imt-1,2:jmt-1), tmpi(2:imt-1), tmpj(2:jmt-1)

      call opennew ('./'//trim(fname)//'.nc', iou)
      call redef (iou)
      call defdim ('xt', iou, imt-2, id_xt)
      call defdim ('yt', iou, jmt-2, id_yt)
      it(:) = id_xt
      call defvar ('xt', iou, 1, it, c0, c0, 'X', 'F'
     &, 'longitude of the t grid', 'grid_longitude', 'degrees_east')
      it(:) = id_yt
      call defvar ('yt', iou, 1, it, c0, c0, 'Y', 'F'
     &, 'latitude of the t grid', 'grid_latitude', 'degrees_north')
      it(1) = id_xt
      call defvar (trim(fname), iou, 2, it, -1000000., 1000000.
     &,   '', 'F', 'output'
     &,   '', 'cm')
      call enddef (iou)
      ib(:) = 1
      ic(:) = imt-2
      tmpi(2:imt-1) = xt(2:imt-1)
      call putvara ('xt', iou, imt-2, ib, ic, tmpi, c1, c0)
      ic(:) = jmt-2
      tmpj(2:jmt-1) = yt(2:jmt-1)
      call putvara ('yt', iou, jmt-2, ib, ic, tmpj, c1, c0)
      ib(:) = 1
      ic(:) = 1
      ic(1) = imt-2
      ic(2) = jmt-2
      tmpij(2:imt-1,2:jmt-1) = field(2:imt-1,2:jmt-1)
      call putvara (trim(fname)
     &, iou, imt*jmt, ib, ic, tmpij, c1, 1)
      call closefile (iou) 
      
      print*, 'Output diagnostic written to ', trim(fname),'.nc'
      
      return
      end subroutine


      subroutine jerncdf_timedep(field,fname,ntrec)
	
      !---------------------------------------------------------------
      !Prints sequential imt x jmt arrays to named netcdf	
      
      implicit none      
      
      include "size.h"	
      include "coord.h" 	
      include "pconst.h" 	
      
      integer it(10), ib(10), ic(10)
      character(*) :: fname
      integer id_xt, id_yt, id_time, iou, ntrec
      real field(1:imt,1:jmt)     
      real tmpij(2:imt-1,2:jmt-1), tmpi(2:imt-1), tmpj(2:jmt-1)

      if (ntrec .eq. 1) then

        call opennew ('./'//trim(fname)//'.nc', iou)
        call redef (iou)
        call defdim ('xt', iou, imt-2, id_xt)
        call defdim ('yt', iou, jmt-2, id_yt)
        call defdim ('time', iou, 0, id_time)
	
        it(1) = id_xt
	it(2) = id_yt
	it(3) = id_time
	
        call defvar ('xt', iou, 1, it(1), c0, c0, 'X', 'F'
     &  , 'longitude of the t grid', 'grid_longitude', 'degrees_east')
        call defvar ('yt', iou, 1, it(2), c0, c0, 'Y', 'F'
     &  , 'latitude of the t grid', 'grid_latitude', 'degrees_north')
        call defvar (trim(fname), iou, 3, it, -1000000., 1000000.
     &  ,   '', 'F', 'output', '', '')   
        call enddef (iou)
      
      else
      
        call openfile('./'//trim(fname)//'.nc',iou)
		
      endif
      
      ib(:) = 1
      ic(:) = imt-2
      tmpi(2:imt-1) = xt(2:imt-1)
      call putvara ('xt', iou, imt-2, ib, ic, tmpi, c1, c0)
      ic(:) = jmt-2
      tmpj(2:jmt-1) = yt(2:jmt-1)
      call putvara ('yt', iou, jmt-2, ib, ic, tmpj, c1, c0)

      ib(:) = 1
      ic(:) = 1
      ic(1) = imt-2
      ic(2) = jmt-2
      ib(3) = ntrec

      tmpij(2:imt-1,2:jmt-1) = field(2:imt-1,2:jmt-1)
      call putvara (trim(fname)         
     &, iou, imt*jmt, ib, ic, tmpij, c1, 1)
      !call closefile (iou) 
      
      print*,'ouput written to ',trim(fname),' at slice ',ntrec
      
      return
      end subroutine

      subroutine jerncdf_ice_snapshot(field,fname)

      !---------------------------------------------------------------
      !Prints nx x ny array to named netcdf.  Will overwrite.	      

      use comicegrid, only: nx,ny,darea
      
      implicit none
      
      integer ib(10), ic(10), iou, id_xt, id_yt
      real field(nx,ny)
      character(*) :: fname

      call opennew ('./'//trim(fname)//'.nc', iou)      
      call redef (iou)
      call defdim ('xt', iou, nx, id_xt)
      call defdim ('yt', iou, ny, id_yt)

      call defvar (trim(fname), iou, 2, (/id_xt,id_yt/), -4.e20 
     &,4.e20, ' ', 'F', 'output', ' ' ,'grid_cells')
      call enddef (iou)

      ib(:) = 1
      ic(:) = 1        
      ic(1) = nx
      ic(2) = ny      
      
      call putvara (trim(fname),iou,nx*ny,ib,ic,field,1.,0.)

      call closefile (iou)     
      
      print*, 'Output diagnostic written to ', trim(fname)     
      
      return
      end subroutine

      subroutine jerncdf_ice_timedep(field,fname,vname,ntrec)

      !---------------------------------------------------------------
      !Prints sequential nx x ny arrays to named netcdf.	      

      use comicegrid
      
      integer it(10), ib(10), ic(10)    
      real field(nx,ny)
      integer ntrec
      character(*) :: fname, vname
      logical exists
      inquire (file='./'//trim(fname)//'.nc', exist=exists)
      
      if (.not. exists) then
        call opennew ('./'//trim(fname)//'.nc', iou)   
      endif
      if (ntrec .eq. 1) then   
        call redef (iou)
        call defdim ('xt', iou, nx, id_xt)
        call defdim ('yt', iou, ny, id_yt)
        call defdim ('time', iou, 0, id_time)	
        it(1) = id_xt
	it(2) = id_yt
	it(3) = id_time
        call defvar ('xt', iou, 1, it(1), 0., 0., 'X', 'F'
     &, 'longitude of the t grid', 'longitude', 'degrees_east')
        call defvar ('yt', iou, 1, it(2), 0., 0., 'Y', 'F'
     &, 'latitude of the t grid', 'latitude', 'degrees_north')
        call defvar (trim(vname), iou, 3, it, -4.e6 
     &,4.e6, ' ', 'F', 'output', ' ' ,'grid_cells')
        call enddef (iou)
      else
        call openfile('./'//trim(fname)//'.nc',iou)
      endif 

      ib(:) = 1
      ic(:) = 1
      ic(1) = nx
      ic(2) = ny
      ib(3) = ntrec

      call putvara (trim(vname), iou, nx*ny, ib, ic, field, 1., 0.)

      !call closefile (iou)     
      
      print*, trim(vname),' written to ', trim(fname),' at slice ',ntrec     
      
      return
      end subroutine
      
      
      subroutine jerncdf_sgi_snapshot(field,fname)

      !---------------------------------------------------------------
      !Prints snapshot of 1D subgrid array to named netcdf.	      

      use subgrid
      
      integer it(10), ib(10), ic(10)    
      real field(numiE)
      character(*) :: fname
      call opennew ('./'//trim(fname)//'.nc', iou)     
      call redef (iou)
      call defdim ('n', iou, 0, id_n)  
      it(1) = id_n
      call defvar ('n', iou, 1, it(1), 0., 0., 'X', 'F'
     &, '', '', '')
      call defvar (trim(fname), iou, 1, it, -4.e6 
     &,4.e6, ' ', 'F', 'output', ' ' ,'blah')
      call enddef (iou)
      ib(:) = 1
      ic(:) = 1
      ic(1) = numiE

      call putvara (trim(fname), iou, numiE, ib, ic, field, 1., 0.)

      !call closefile (iou)     
      
      print*, 'Output written to ', trim(fname)  
      
      return
      end subroutine
      
      subroutine jerncdf_sgi2ice_timedep(sgarr,nsgarr,fname,vname,ntrec)
      !for remapping an arbitrary subgrid array to available ice sheets.  Will overwrite subgrid field on an
      !underlying non-subgridded climate array (mapped to ice sheet).  For example, it will plot subgridded albedo where
      !ice occurs, then plot the non-subgrid-generated albedo where no ice exists in ice sheet grid.
      
      use ice_sheet_storage, only: strtis,endis
      use comicesheet, only: h
      use comicegrid, only: nx,ny,is2sg, isname,iceiind,icejind
      use subgrid, only: numiE,himin     
      
      implicit none
      
      include "size.h"
      
      real sgarr(numiE)
      real, allocatable :: icearr(:,:)
      real nsgarr(imt,jmt)
      integer i,j,n,ntrec
      character(*) fname,vname
      
      do n=strtis,endis
        !point at an ice sheet and allocate diagnostic array
	call set_pointers(n)
	allocate (icearr(nx,ny))
	icearr(:,:)=0.
     
        do i=1,nx
          do j=1,ny
	    if (h(i,j).le.himin) then 
	      icearr(i,j)=nsgarr(iceiind(i,j),icejind(i,j))
	    else
	      icearr(i,j)=sgarr(is2sg(i,j)) 
	    endif 
          enddo
        enddo

        call jerncdf_ice_timedep(icearr,fname,vname,ntrec)
	deallocate (icearr)	
      enddo    
      
      return
      end subroutine
