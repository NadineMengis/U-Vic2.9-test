! source file: /net/mare/home1/eby/as/ism/interpolator.F
      subroutine downscale_to_ice_grid(arro,arri,dobilin)

      use comicegrid

      implicit none

      !Sets the global grid cell locations (and corresponding weights)
      !that each ice sheet cell uses to interpolate its global field
      !values from.  Thanks Wikypedia for teaching me all I know.

      ! bottom left corner (SW)  = lonind(1),latind(1),wei(1)
      ! bottom right corner (SE) = lonind(2),latind(1),wei(2)
      ! top left corner (NW)     = lonind(1),latind(2),wei(3)
      ! top right corner (NE)    = lonind(2),latind(2),wei(4)

      !Initial output (saved in local routine)
      !lonind(2) = lons of nearest global grid cell centers ([west east])
      !latind(2) = lats of nearest global grid cell centers ([south north])
      !wei(nx,ny,4) = bilinear weights (see above for corresponding corners)

      !arro: input global array
      !arri: output ice sheet array
      !dobilin: 0=nearest neighbor, 1=bilinear interpolation

      include "size.h"			!contains imt,jmt
      include "param.h"                 !contains imtm1,jmtm1
      include "coord.h"			!contains xt,yt

      !include "atm.h"			!includes at(i,j,2,isat)

      real ice_lon(nx,ny),ice_lat(nx,ny),glob_lon(imt),glob_lat(jmt)
      integer  i,j,ig,jg,dobilin

      real lonval(2),londiff,latval(2),latdiff
      real arro(imt,jmt),arri(nx,ny), weimax
      character(120) g_st

      !If first swipe through, set downscaling arrays
      if(firstdownscale) then

        !ESCM center-of-grid-cell latitudes
        !xt = longitude of central global cell points
        !     (xt(2:101) = 1.8:358.2) (0 = Greenwich)
        !yt = latitude of central cell points
        !     (yt(2:101) = -89.1:89.1) (-90=S Pole, 90=N Pole)

        !Ice sheet center-of-grid-cell latitudes
        !alond = longitude of central ice cell points
        !      0 = Greenwich, + eastwards, - westwards
        !alatd = latitude (-90=S Pole, 90=N Pole)

        !Interpolation for dummies: latitudes/longitudes are not
        !constant along stereographic i,j grids (unlike lat-lon grids)

        ice_lon = alond
        ice_lat = alatd
        glob_lon = xt
        glob_lat = yt

        !convert ice grid longitudes to all positive
        do i=1,nx
          do  j=1,ny
     	    if (alond(i,j) .lt. 0.) then
	      ice_lon(i,j) = mod (ice_lon(i,j) + 720., 360.)
	    endif
	  enddo
        enddo

        !for each ice cell, find surrounding global cells
        do i=1,nx
          do  j=1,ny
	    !search through global longitude array
	    lonloop:do ig=2,imt
	      !first time the global longitude exceeds the ice longitude,
	      if (glob_lon(ig) .gt. ice_lon(i,j)) then
	        !set surrounding global longitudes
                lonind(i,j,1)=ig-1
	        lonind(i,j,2)=ig
	        lonval(1) = glob_lon(ig-1)
	        lonval(2) = glob_lon(ig)
	        londiff=lonval(2)-lonval(1)
	        exit lonloop
	      endif
	    enddo lonloop

	    !search through global latitude array
	    latloop:do jg=3,jmtm1
	      !first time the global latitude exceeds the ice latitude,
	      if(glob_lat(jg) .gt. ice_lat(i,j)) then
	        !set surrounding global latitudes
                latind(i,j,1)=jg-1
	        latind(i,j,2)=jg
	        latval(1) = glob_lat(jg-1)
	        latval(2) = glob_lat(jg)
	        latdiff=latval(2)-latval(1)
	        exit latloop
	      endif
	    enddo latloop

	    !set bilinear interpolation weights

	    weid(i,j,1)=(lonval(2)-ice_lon(i,j))*(latval(2)-ice_lat(i,j))
     &                  /(londiff*latdiff)
	    weid(i,j,2)=(ice_lon(i,j)-lonval(1))*(latval(2)-ice_lat(i,j))
     &                  /(londiff*latdiff)
	    weid(i,j,3)=(lonval(2)-ice_lon(i,j))*(ice_lat(i,j)-latval(1))
     &                  /(londiff*latdiff)
	    weid(i,j,4)=(ice_lon(i,j)-lonval(1))*(ice_lat(i,j)-latval(1))
     &                  /(londiff*latdiff)
          enddo
        enddo
	firstdownscale=.false.
      endif

      do i=1,nx
        do  j=1,ny
          if (dobilin.eq.1) then
	    arri(i,j) = weid(i,j,1)*arro(lonind(i,j,1),latind(i,j,1))
     &  	      + weid(i,j,2)*arro(lonind(i,j,2),latind(i,j,1))
     &		      + weid(i,j,3)*arro(lonind(i,j,1),latind(i,j,2))
     &		      + weid(i,j,4)*arro(lonind(i,j,2),latind(i,j,2))
          else
            weimax=max(weid(i,j,1),weid(i,j,2),weid(i,j,3),weid(i,j,4))
            if (weimax.eq.weid(i,j,1)) then
              arri(i,j) = arro(lonind(i,j,1),latind(i,j,1))
            else if (weimax.eq.weid(i,j,2)) then
              arri(i,j) = arro(lonind(i,j,2),latind(i,j,1))
            else if (weimax.eq.weid(i,j,3)) then
              arri(i,j) = arro(lonind(i,j,1),latind(i,j,2))
            else if (weimax.eq.weid(i,j,4)) then
              arri(i,j) = arro(lonind(i,j,2),latind(i,j,2))
            endif
	  endif
        enddo
      enddo
      return
      end subroutine

!-----------------------------------------------------------------------

      subroutine upscale_to_global_grid(arro,arri)

      !Upscales ice sheet field arro to global grid arri.  The global
      !grid is first subdivided into a finer mesh, of roughly the same
      !resolution as ice sheet grid.
      !Then, the center latitude of each smaller grid cell (the 'global
      !interpolated grid', or 'gig' cell) is converted to an ij
      !coordinate that corresponds to the ice sheet grid reference, using
      !subroutine ll2ij (a snippet of Pollard code which may not apply
      !perfectly to Greenland).  This ij coordinate is used in the bilinear
      !interpolation to obtain the interpolated 'gig' cell value.  Finally,
      !these gig values are area-weighted and summed to give the global
      !cell value.

      use comicegrid, only: nx,ny,dx0,dy0, isname

      implicit none

      include "size.h"			!contains imt,jmt
      include "param.h"                 !contains imtm1,jmtm1
      include "coord.h"			!contains xt,yt
      include "grdvar.h"		!contains dxt,dyt,cst

      integer ngiglon(imt),ngiglat(jmt)
      real
     &      gwidth                      !EW distance of global cell at lat j
     &     ,gheight			!NW distance of global cell at lat j
      real
     &	    arro(nx,ny),arri(imt,jmt)   !input ice array/output global array
     &     ,circum                      !circumfrence of earth
     &     ,gwlon,gslat		        !west/south edge of global cell
     &     ,dlon,dlat                   !EW/NS distance across global cell
     &     ,garea			!area of global cell (m^2)
     &     ,gigcentlon,gigcentlat       !central lon/lat of gig cell
     &     ,gigarea			!area of gig cell (m^2)
     &     ,gigi,gigj			!central i,j coords of gig cell
     &     ,weiu(4)			!weights of [bl br tl tr] corners
     &     ,arrigig			!bilinear interpolated gig cell value
     &     ,pi,radius,deg2rad           !constants

      integer
     &      i,j,ii,jj,jstart,jend       !counters
     &     ,iind(2),jind(2)             !i,j row indices of ice sheet grid

      parameter (pi	= 3.14159265358979)
      parameter (radius = 6.37122e6)

      circum=2.*pi*radius/360.
      deg2rad=pi/180.

      !Set upscaling arrays
      do j=1,jmt
        !get EW distance across middle of a cell at latitude row j
  	gwidth=dxt(1)*cst(j)/100.

        !determine the EW count of gig cells within one global cell
	!by dividing the global grid width by the ice sheet grid width
	ngiglon(j)=max(1,floor(gwidth/dx0))

	!get NW distance across middle of a cell at latitude row j
	gheight=dyt(j)/100.

	!determine the NS count of gig cells within one global cell
	!by dividing the global grid width by the ice sheet grid width
        ngiglat(j)=max(1,floor(gheight/dy0))
      enddo

      !zero output array
      arri(:,:) = 0.0

      if (isname == 'Antarctica') then
         jstart=2
	 jend=50
      elseif (isname == 'Greenland') then
        jstart=51
	jend=jmtm1
      endif
      do j=jstart,jend!for some reason (in ll2ij) does a mirror image thing

	do i=2,imtm1
	  !determine longitude of western edge of large global cell
	  gwlon=tlon(i,j)-0.5*dxtdeg(i)
	  !determine latitude of southern edge of large global cell
	  gslat=tlat(i,j)-0.5*dytdeg(j)
	  !determine area of large global cell
	  garea=cst(j)*dyt(j)*dxt(i)/10000.
	  !determine d[lat/lon] for gig cells
	  dlon=dxtdeg(i)/ngiglon(j)
	  dlat=dytdeg(j)/ngiglat(j)

	  do jj=1,ngiglat(j)

	    !determine central latitude of gig cell
	    gigcentlat=gslat+(jj-0.5)*dlat
	    !determine area of gig cell
	    gigarea=dlon*cos(gigcentlat*deg2rad)*dlat*circum**2

	    do ii=1,ngiglon(j)
	      !determine central longitude of gig cell
	      gigcentlon=gwlon+(ii-0.5)*dlon
	      !convert to -180=>0=>180 longitude grid
	      if (gigcentlon .gt. 180.) then
	  	gigcentlon=gigcentlon-360.
	      endif

	      !convert centre lon/lat to i,j coordinates on stereographic ice grid

	      if (isname == 'Antarctica') then
	        call ll2ijant(gigcentlon,gigcentlat,gigi,gigj)
	      elseif (isname == 'Greenland') then
	        call ll2ijgre(gigcentlon,gigcentlat,gigi,gigj)
	      endif

	      !set surrounding ice point indexes
	      if(gigi .lt. nx
     &	        .and. gigi .gt. 1.
     &    	.and.
     &	  	 gigj .lt. ny
     &    	.and. gigj .gt. 1.) then

	  	iind(1) = real(floor(gigi))
	  	iind(2) = real(ceiling(gigi))

	  	jind(1) = real(floor(gigj))
	  	jind(2) = real(ceiling(gigj))

	  	!set weights (implicitly divided by 1*1)
	  	weiu(1)=(iind(2)-gigi)*(jind(2)-gigj)
	  	weiu(2)=(gigi-iind(1))*(jind(2)-gigj)
	  	weiu(3)=(iind(2)-gigi)*(gigj-jind(1))
	  	weiu(4)=(gigi-iind(1))*(gigj-jind(1))

	  	!determine values of gig cells by bilinear interpolation
	  	arrigig = weiu(1)*arro(iind(1),jind(1))
     &          	+ weiu(2)*arro(iind(2),jind(1))
     &	        	+ weiu(3)*arro(iind(1),jind(2))
     &	        	+ weiu(4)*arro(iind(2),jind(2))

	  	!area-weight value, add to global grid cell value
	  	arri(i,j) = arri(i,j)+arrigig*(gigarea/garea)

	      endif

	    enddo
	  enddo

	enddo
      enddo

      call embmbc(arri)

      return
      end subroutine

!-----------------------------------------------------------------------

      subroutine downscale_corrector(arro,arromsk,arri,arrimsk)

      !Use to adjust the sum of arri to equal the sum of arro,
      !after interpolation.
      !arro      =input global array
      !arri      =input ice array to be corrected
      !arromsk   =mask of pertinent (=1) global array points
      !arrimsk   =mask of pertinent (=1) ice array points to be corrected

      use comicegrid

      implicit none

      include "size.h"           !includes imt,jmt

      real     arro(imt,jmt),arromsk(imt,jmt)
     &        ,arri(nx,ny),arrimsk(nx,ny)
     &        ,arrosum,arrisum,arrisumcor,area,ratio

      integer i,j

      !original total sum
      call areatot(arro,arromsk,arrosum)
      print*, 'arrosum=', arrosum
      !uncorrected total sum
      arrisum=0.
      do i=1,nx
       do j=1,ny
	arrisum=arrisum+(arri(i,j)*darea(i,j)*1.e4*arrimsk(i,j))
       enddo
      enddo
      print*, 'arrisum=', arrisum
      if (arrisum .ne. 0.) then
        ratio=arrosum/arrisum
        do i=1,nx
          do j=1,ny
	    arri(i,j)=arri(i,j)*ratio
	  enddo
        enddo
      else
        arri(:,:) = 0.0
      endif

      !corrected total sum
      arrisumcor=0.
      do i=1,nx
       do j=1,ny
	arrisumcor=arrisumcor+(arri(i,j)*darea(i,j)*1.e4*arrimsk(i,j))
       enddo
      enddo
      print*, 'arrisumcor=', arrisumcor

      if (arrisumcor .ne. 0.0) then
        if(abs(1.-arrisum/arrosum) .gt. .1) then
          print*, 'Warning: uncorrected downscaled field differs
     &	  from original global field sum by', 1.-(arrisum/arrosum)
        endif

        if(abs(1.-arrisumcor/arrosum) .gt. 1.+1.e-14) then
          print*,'=> interpolator.F error: significant nonconservation
     &    during downscaling correction, arrisumcor/arrosum='
     &    ,arrisumcor/arrosum
          stop
        endif
        print*, 'downscale correction accuracy=',
     &   abs(1.-arrisumcor/arrosum)
      endif

      return
      end subroutine

!-----------------------------------------------------------------------

      subroutine upscale_corrector(arro,arromsk,arri,arrimsk)

      !Use to adjust the sum of arri to equal the sum of arro,
      !after interpolation.
      !arro      =input ice array
      !arri      =input global array to be corrected
      !arromsk   =mask of pertinent ice array points
      !arrimsk   =mask of pertinent global array points

      use comicegrid

      implicit none

      include "size.h"           !includes imt,jmt
      include "param.h"          !includes imtm1,jmtm1
      include "pconst.h"         !includes epsln

      real     arro(nx,ny),arromsk(nx,ny)
     &        ,arri(imt,jmt),arrimsk(imt,jmt)
     &        ,arrosum,arrisum,arrisumcor,area
     &	      ,ratio

      integer i,j
      character(120) g_st

      area=dx0*dy0*10000. !m^2 to cm^2, to correspond to areatot calc

      !original total sum
      arrosum=0.
      do i=1,nx
        do j=1,ny
	  arrosum=arrosum+(arro(i,j)*darea(i,j)*1.e4*arromsk(i,j))
        enddo
      enddo

      !uncorrected total sum
      call areatot(arri,arrimsk,arrisum)

      if (arrisum .gt. 0.) then
        ratio=arrosum/arrisum
        do i=2,imtm1
          do j=2,jmtm1
	    arri(i,j)=arri(i,j)*ratio
	  enddo
        enddo
      else
        arri(:,:) = 0.0
      endif

      !corrected total sum
      call areatot(arri,arrimsk,arrisumcor)

      if (arrisumcor .gt. 0.0) then
        if(abs(1.-arrisum/(arrosum)) .gt. .1) then
          print*, 'Error: uncorrected upscaled field differs by',
     &    (abs(1.-(arrisum/(arrosum)))*100.),'%'
          stop
        endif
        if(abs(1.-arrisumcor/(arrosum)) .gt. 1.e-14) then
          print*,'Error: significant nonconservation
     &    during upscaling correction, arrisumcor/arrosum='
     &    ,arrisumcor/(arrosum)
          stop
        endif
        print*, 'upscale correction accuracy=',
     &abs(1.-arrisumcor/arrosum)
      endif

      call embmbc(arri)

      return
      end subroutine

!-----------------------------------------------------------------------
!Utility subroutines
!-----------------------------------------------------------------------

      subroutine global_cell_area(i,j,area)

      implicit none

      include "size.h"
      include "grdvar.h"

      integer i,j
      real    area, fx

      area=0.
      area = dxt(i)*cst(j)*dyt(j)

      return
      end subroutine

!-----------------------------------------------------------------------

      subroutine ll2ijant(lon,lat,steri,sterj)

      use comicegrid

      implicit none

      real     zcos,zr,zxa,zya
      real     lon,lat,steri,sterj

!     Convert lon,lat to i,j polar stereo coords (Pollard code snippet)
      zcos = cos(lat*pi/180.)
      zr = ( (1./zcos) - sqrt ((1./zcos)**2 - 1.) ) * zlambda
      zxa = radius * zr * cos (0.5*pi - lon*pi/180.)
      zya = radius * zr * sin (0.5*pi - lon*pi/180.)
      steri = (zxa + 0.5*nx*dx0)/dx0 + 1.0001
      sterj = (zya + 0.5*ny*dy0)/dy0 + 1.0001

      return
      end subroutine

!-----------------------------------------------------------------------

      subroutine ll2ijgre(lon,lat,steri,sterj)

      use comicegrid

      implicit none

      real     zcos,zr,zxa,zya
      real     lon,lat,steri,sterj
      lon=lon+39.
      zcos = cos(lat*pi/180.)
      zr = ( (1./zcos) - sqrt ((1./zcos)**2 - 1.) ) * zlambda
      !EW distance from point to bottom left corner
      zxa = radius * zr * cos (0.5*pi - lon*pi/180.) + xoffa
      !NS distance from point to bottom left corner
      zya = radius * zr * sin (0.5*pi - (lon)*pi/180.) + yoffa

      steri = (zxa + 0.5*nx*dx0)/dx0 + 1.0001
      sterj = (zya + 0.5*ny*dy0)/dy0 + 1.0001
      !flip array with center latitude as axis
      sterj = ny - sterj
      return
      end subroutine

