      subroutine slope_calculator(c)

#if defined O_topography_shading      
      use comicegrid
      use comicesheet
      
      implicit none
      
      real dhdx(nx,ny),dhdy(nx,ny)
      real a(3),b(3)
      real c(nx,ny,3) !components of vector normal to surface slope
      integer i,j
      character(120) g_st
      
      !calculates normal vector to the slope at each grid point, based on a look-around of 8 neighbouring points.
      dhdx(:,:)=0.
      dhdy(:,:)=0.
      c(:,:,1:2)=0. 
      c(:,:,3)=1. !set initial normal vector to unit count in vertical direction, implying a horizontal surface    
      do i=2,nx-1
        do j=2,ny-1
	  a(:)=0.
	  b(:)=0
	  !set components of slope vectors in x and y directions for cross product calculation
	  !a=<x-direction slope vector>, b=<y-direction slope vector>
	  !1=x component
	  !2=y component
	  !3=z component
	  a(1)=2*dx(i,j)   
	  a(2)=0.      
	  a(3)=hs(i+1,j)-hs(i-1,j)
	  b(1)=0.
	  b(2)=2*dy(i,j)
	  b(3)=hs(i,j+1)-hs(i,j-1)

	  !compute cross product vector
          c(i,j,1) = a(2)*b(3) - a(3)*b(2)
          c(i,j,2) = a(3)*b(1) - a(1)*b(3)
          c(i,j,3) = a(1)*b(2) - a(2)*b(1)
	enddo
      enddo
#endif      
      return
      end
!          ****************************************************************************     
      FUNCTION ANGVEC (V1, V2) RESULT (COSTHETA)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: V1, V2
      DOUBLE PRECISION :: COSTHETA
      COSTHETA = DOT_PRODUCT(V1,V2)/
     &(SQRT(DOT_PRODUCT(V1,V1))*SQRT(DOT_PRODUCT(V2,V2)))
      RETURN
      END FUNCTION ANGVEC
      
      
      
	  
