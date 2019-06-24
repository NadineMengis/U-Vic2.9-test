#if defined O_ism_standalone
      program standalone      

      use ice_sheet_storage, only: strtis,endis
           
      implicit none
      
      integer n, ioun
      logical runais, rungis
      namelist /ismin/ runais, rungis  
              
      runais=.false.
      rungis=.false.
      strtis=0.
      endis=0.
      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, ismin, end=101)
101   continue      
      if (runais .and. rungis) then
        strtis=1
        endis=2
      elseif (runais .and. .not. rungis) then
        strtis=1
        endis=1
      elseif (.not. runais .and. rungis) then
        strtis=2
        endis=2
      else
        stop 'Error: neither GIS/AIS specified in control.in.'
      endif        

      do n=strtis,endis
	call set_ice_sheet_DDT(n)    
	call set_pointers(n)	
	call ism(n)              
      enddo     

      do n=strtis,endis
	call set_pointers(n)      
 	call ism(n)
      enddo

      end program
#endif
