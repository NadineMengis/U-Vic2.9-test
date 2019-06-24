      subroutine momout
      
      implicit none
      
      include "size.h"
      include "tmngr.h"
      include "switch.h"
      
      character (120) :: fname      
!-----------------------------------------------------------------------
!     save restart
!-----------------------------------------------------------------------

      if (restrt) then
        if (restts) then
          call def_rest (0)
          call def_rest_mom (0, fname)
          call mom_rest_out (fname, 1, imt, 1, jmt)
        endif
        if (eorun) then
          call def_rest (1)
          call def_rest_mom (1, fname)
          call mom_rest_out (fname, 1, imt, 1, jmt)
        endif
      endif      
      
      return
      end
