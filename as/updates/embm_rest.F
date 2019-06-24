      subroutine embm_rest_in (fname, ids, ide, jds, jde)

#if defined O_embm
!=======================================================================
!     input routine for atmospheric restarts

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname              = file name
!     ids, ide ...       = start and end index for data domain
!=======================================================================
# if defined O_ism      
      use subgrid
# endif      
      implicit none

      character(*) :: fname
      character(32) :: nstamp
      character(3) :: a3
      character(120) :: var1, var2
      integer iou, ln, n, ntrec, ids, ide, jds, jde, ig
      integer ils, ile, jls, jle, ib(10), ic(10)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      logical inqvardef, exists, use_sealev_data, use_ice_data

      real tmp
      real, allocatable :: tmpi(:), tmpij(:,:)

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "atm.h"
      include "cembm.h"
      include "coord.h"
      include "csbc.h"
      include "grdvar.h"
# if defined O_ice_cpts && defined O_ice
      include "cpts.h"
# endif
      include "ice.h"
# if defined O_ice_evp && defined O_ice
      include "evp.h"
# endif
      include "tmngr.h"
      include "switch.h"
# if defined O_ism
      include "ism.h"
      real temp(imt,jmt)
      integer i,j
# endif
!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)
      ntrec = 1

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global read domain)
!-----------------------------------------------------------------------
      ils = max(ids,1)
      ile = min(ide,imt)
      jls = max(jds,1)
      jle = min(jde,jmt)

      allocate ( tmpij(ils:ile,jls:jle) )

!-----------------------------------------------------------------------
!     read 1d data (t)
!-----------------------------------------------------------------------
      tmp = nats
      call getvars ('nats', iou, ntrec, tmp, c1, c0)
      nats = tmp
      tmp = dayoyr
      call getvars ('dayoyr', iou, ntrec, tmp, c1, c0)
      dayoyr = tmp
      tmp = itt
      call getvars ('itt', iou, ntrec, tmp, c1, c0)
      itt = tmp
      tmp = irstdy
      call getvars ('irstdy', iou, ntrec, tmp, c1, c0)
      irstdy = tmp
      tmp = msrsdy
      call getvars ('msrsdy', iou, ntrec, tmp, c1, c0)
      msrsdy = tmp
      tmp = totaltime
      call getvars ('totaltime', iou, ntrec, tmp, c1, c0)
      totaltime = tmp
      tmp = year0
      call getvars ('year', iou, ntrec, tmp, c1, c0)
      nyear = tmp
      tmp = month0
      call getvars ('month', iou, ntrec, tmp, c1, c0)
      nmonth = tmp
      tmp = day0
      call getvars ('day', iou, ntrec, tmp, c1, c0)
      nday = tmp
      tmp = hour0
      call getvars ('hour', iou, ntrec, tmp, c1, c0)
      nhour = tmp
      tmp = min0
      call getvars ('minute', iou, ntrec, tmp, c1, c0)
      nmin = tmp
      tmp = sec0
      call getvars ('second', iou, ntrec, tmp, c1, c0)
      nsec = tmp
      call mkstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      if (init_time_in) then
        itt = 0
        irstdy = 0
        msrsdy = 0
        relyr = 0.0
        call mkstmp (stamp, year0, month0, day0, hour0, min0, sec0)
      endif
# if defined O_carbon && !defined O_carbon_uncoupled && !defined O_co2ccn_data && !defined O_co2ccn_data_transient
      tmp = co2ccn
      call getvars ('co2ccn', iou, ntrec, tmp, c1, c0)
      co2ccn = tmp
# endif
# if defined O_carbon_14 && defined O_carbon_14_coupled && !defined O_c14ccn_data && !defined O_c14ccn_data_transient
      tmp = c14ccn
      call getvars ('c14ccn', iou, ntrec, tmp, c1, c0)
      c14ccn = tmp
# endif
# if defined O_co2emit_track_co2
      tmp = ntrack_co2
      call getvars ('ntrack_co2', iou, ntrec, tmp, c1, c0)
      if (ntrack_co2 .ne. nint(tmp)) then
        print*, "Warning => ignoring restart values for co2 tracking"
      else
        ntrack_co2 = nint(tmp)
        tmp = itrack_co2
        call getvars ('itrack_co2', iou, ntrec, tmp, c1, c0)
        itrack_co2 = min(nint(tmp), ntrack_co2-1)
        ib(1) = 1
        ic(1) = ntrack_co2
        call getvara ('track_co2', iou, ic(1), ib, ic, track_co2
     &, c1, c0)
      endif
# endif
# if defined O_co2emit_track_sat || defined O_embm_vcs
      tmp = ntrack_sat
      call getvars ('ntrack_sat', iou, ntrec, tmp, c1, c0)
      if (ntrack_sat .ne. nint(tmp)) then
        print*, "Warning => ignoring restart values for sat tracking"
      else
        ntrack_sat = nint(tmp)
        tmp = itrack_sat
        call getvars ('itrack_sat', iou, ntrec, tmp, c1, c0)
        itrack_sat = min(nint(tmp), ntrack_sat-1)
        ib(1) = 1
        ic(1) = ntrack_sat
        call getvara ('track_sat', iou, ic(1), ib, ic, track_sat
     &, c1, c0)
      endif
# endif
      use_ice_data = .true.
# if defined O_ice_data_transient
!     for large changes in ice_yr ignore restart data
      tmp = ice_yr
      use_ice_data = .false.
      if (inqvardef('ice_yr', iou)) then
        call getvars ('ice_yr', iou, ntrec, tmp, c1, c0)
        if (abs(tmp - ice_yr) .lt. 2.) use_ice_data = .true.
      endif
      if (.not. use_ice_data) then 
        print*, ' '
        print*, '==> Warning: restart ice year is different'
        print*, ' '
      endif
# endif
      use_sealev_data = .true.
# if defined O_sealev_data_transient
!     for large changes in sealev_yr ignore restart data
      tmp = sealev_yr
      use_sealev_data = .false.
      if (inqvardef('sealev_yr', iou)) then
        call getvars ('sealev_yr', iou, ntrec, tmp, c1, c0)
        if (abs(tmp - sealev_yr) .lt. 2.) use_sealev_data = .true.
      endif
      if (.not. use_sealev_data) then 
        print*, ' '
        print*, '==> Warning: restart sea level year is different'
        print*, ' '
      endif
# endif
# if defined O_sealev || defined O_sealev_data
      tmp = sealev
      call getvars ('sealev', iou, ntrec, tmp, c1, c0)
      if (.not. inqvardef('sealev', iou)) then
        dsealev = sealev
      else
#  if defined O_sealev_data_transient
        if (use_sealev_data) then
          sealev = tmp
          dsealev = 0.
        endif
#  endif
        if (tmp .ne. sealev) then
          print*, ' '
          print*, '==> Warning: restart sea level is different, ',
     &      'using newly specified sea level'
          print*, ' '
          dsealev = sealev - tmp
        endif
      endif
# endif
# if defined O_ism
!-----------------------------------------------------------------------
!     read 1d data into ice basin heat and water reservoir arrays
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = 200
      allocate(tmpi(200))
      
      tmpi(:)=dicevol(:)
      call getvara ('dicevol', iou, 200, ib, ic, tmpi, c1, c0)
      dicevol(:)=tmpi(:)
 
      tmpi(:)=diceheat(:)
      call getvara ('diceheat', iou, 200, ib, ic, tmpi, c1, c0)
      diceheat(:)=tmpi(:)
      
      tmpi(:)=isflxm(:)
      call getvara ('isflxm', iou, 200, ib, ic, tmpi, c1, c0)
      isflxm(:)=tmpi(:)
      
      tmpi(:)=isflxh(:)
      call getvara ('isflxh', iou, 200, ib, ic, tmpi, c1, c0)
      isflxh(:)=tmpi(:)
      
      deallocate(tmpi)
      
      tmp=dtism
      call getvars ('dtism', iou, 1, tmp, c1, c0)
      dtism=tmp     
      
      tmp=real(numiE)
      call getvars ('numiE', iou, 1, tmp, c1, c0)                 
      numiE=nint(tmp)
      
      tmp=rtag
      call getvars ('rtag', iou, 1, tmp, c1, c0)
      rtag=tmp
      
      tmp=nisaccum
      call getvars ('L_accummassN', iou, 1, tmp, c1, c0)
      nisaccum=tmp

      tmp=sisaccum
      call getvars ('L_accummassS', iou, 1, tmp, c1, c0)
      sisaccum=tmp

      tmp=nismelt            
      call getvars ('L_meltmassN', iou, 1, tmp, c1, c0)
      nismelt=tmp

      tmp=sismelt            
      call getvars ('L_meltmassS', iou, 1, tmp, c1, c0)  
      sismelt=tmp

      tmp=nisiar            
      call getvars ('L_iceareaN', iou, 1, tmp, c1, c0)
      nisiar=tmp

      tmp=sisiar      
      call getvars ('L_iceareaS', iou, 1, tmp, c1, c0)
      sisiar=tmp

      tmp=nissmb      
      call getvars ('L_icesmbN', iou, 1, tmp, c1, c0)
      nissmb=tmp

      tmp=sissmb            
      call getvars ('L_icesmbS', iou, 1, tmp, c1, c0)
      sissmb=tmp

      tmp=nisvol            
      call getvars ('L_icevolN', iou, 1, tmp, c1, c0)
      nisvol=tmp

      tmp=sisvol            
      call getvars ('L_icevolS', iou, 1, tmp, c1, c0)
      sisvol=tmp

      tmp=nismext            
      call getvars ('L_meltextentN', iou, 1, tmp, c1, c0)           
      nismext=tmp

      tmp=sismext            
      call getvars ('L_meltextentS', iou, 1, tmp, c1, c0)     
      sismext=tmp            

      if (numiE .gt. 0) then                     
        allocate(
     &   tmpi(numiE)
     &  ,eleviE(numiE)
     &  ,fraciE(numiE)
     &  ,iariE(numiE)     
     &  ,isniE(numiE)     
     &  ,hsnowiE(numiE,2)      
     &  ,ticeiE(numiE)
     &  ,ismticeiE(numiE)     
     &  ,ithkiE(numiE)     
     &  ,mbalmelt(numiE,2)     
     &  ,mbalaccum(numiE,2) 
     &  ,precipiE(numiE)
     &  ,satiE(numiE)
     &  ,albedoiE(numiE)
     &  ,mextiE(numiE)
     &  ,mduriE(numiE)
# if defined O_refreeze     
     &  ,hnewsno(numiE)
     &  ,counter(numiE)     
     &  ,htemp(numiE)     
     &  ,hreftot(numiE)     
     &  ,meltprev(numiE)          
# endif                  
     &  )
        tmpi(:)=0.
        eleviE(:)=0.
        fraciE(:)=0.
        iariE(:)=0.     
        isniE(:)=0.     
        hsnowiE(:,:)=0.      
        ticeiE(:)=0.
        ismticeiE(:)=0.     
        ithkiE(:)=0.	
        mbalmelt(:,:)=0.     
        mbalaccum(:,:)=0. 
        precipiE(:)=0.
        satiE(:)=0.
        albedoiE(:)=0.
        mextiE(:)=0.
        mduriE(:)=0.
# if defined O_refreeze      
        hnewsno(:)=0.
        counter(:)=0.
        htemp(:)=0.
        hreftot(:)=0.      
        meltprev(:)=0.      		
# endif
        tmpi(:)=0.      
        ib(1) = 1
        ic(1) = numiE
	
	tmpi(:)=eleviE(:)
        call getvara ('eleviE',iou, numiE, ib, ic, tmpi, c1, c0)
        eleviE(:)=tmpi(:) 
        
	tmpi(:)=fraciE(:)
	call getvara ('fraciE',iou, numiE, ib, ic, tmpi, c1, c0)
        fraciE(:)=tmpi(:)
        
	tmpi(:)=iariE(:)
	call getvara ('iariE',iou, numiE, ib, ic, tmpi, c1, c0)
        iariE(:)=tmpi(:)      
        
	tmpi(:)=isniE(:)
	call getvara ('isniE',iou, numiE, ib, ic, tmpi, c1, c0)
        isniE(:)=nint(tmpi(:))       
# if !defined O_icefree_init

	tmpi(:)=hsnowiE(:,1)
	call getvara ('hsnowiE1',iou, numiE, ib, ic, tmpi, c1, c0)
        hsnowiE(:,1)=tmpi(:)
        
	tmpi(:)=hsnowiE(:,2)
	call getvara ('hsnowiE2',iou, numiE, ib, ic, tmpi, c1, c0)
        hsnowiE(:,2)=tmpi(:)
# endif        
        
	tmpi(:)=ticeiE(:)
	call getvara ('ticeiE',iou, numiE, ib, ic, tmpi, c1, c0)
        ticeiE(:)=tmpi(:)          
        
	tmpi(:)=ismticeiE(:)
	call getvara ('ismticeiE',iou, numiE, ib, ic, tmpi, c1, c0)
        ismticeiE(:)=tmpi(:)       
        
	tmpi(:)=ithkiE(:)
	call getvara ('ithkiE',iou, numiE, ib, ic, tmpi, c1, c0)
        ithkiE(:)=tmpi(:)       
        
	tmpi(:)=mbalmelt(:,1)
	call getvara ('mbalmelt1',iou, numiE, ib, ic, tmpi, c1, c0)
        mbalmelt(:,1)=tmpi(:)                 
        
	tmpi(:)=mbalmelt(:,2)
	call getvara ('mbalmelt2',iou, numiE, ib, ic, tmpi, c1, c0)
        mbalmelt(:,2)=tmpi(:)       
        
	tmpi(:)=mbalaccum(:,1)
	call getvara ('mbalaccum1',iou, numiE, ib, ic, tmpi, c1, c0)
        mbalaccum(:,1)=tmpi(:)       
        
	tmpi(:)=mbalaccum(:,2)
	call getvara ('mbalaccum2',iou, numiE, ib, ic, tmpi, c1, c0)
        mbalaccum(:,2)=tmpi(:)
        
	tmpi(:)=precipiE(:)
	call getvara ('precipiE',iou, numiE, ib, ic, tmpi, c1, c0)
        precipiE(:)=tmpi(:) 
        
	tmpi(:)=satiE(:)
	call getvara ('satiE',iou, numiE, ib, ic, tmpi, c1, c0)
        satiE(:)=tmpi(:)       
        
	tmpi(:)=albedoiE(:)
	call getvara ('albedoiE',iou, numiE, ib, ic, tmpi, c1, c0)
        albedoiE(:)=tmpi(:)       
        
	tmpi(:)=mextiE(:)
	call getvara ('mextiE',iou, numiE, ib, ic, tmpi, c1, c0)
        mextiE(:)=tmpi(:)       
        
	tmpi(:)=mduriE(:)
	call getvara ('mduriE',iou, numiE, ib, ic, tmpi, c1, c0)
        mduriE(:)=tmpi(:)

# if defined O_refreeze      
        tmpi(:)=hnewsno(:)
	call getvara ('hnewsno',iou, numiE, ib, ic, tmpi, c1, c0)
        hnewsno(:)=tmpi(:)
        
	tmpi(:)=counter(:)
	call getvara ('counter',iou, numiE, ib, ic, tmpi, c1, c0)
        counter(:)=tmpi(:)      
        
	tmpi(:)=htemp(:)
	call getvara ('htemp',iou, numiE, ib, ic, tmpi, c1, c0)
        htemp(:)=tmpi(:)      
        
	tmpi(:)=hreftot(:)
	call getvara ('hreftot',iou, numiE, ib, ic, tmpi, c1, c0)
        hreftot(:)=tmpi(:)      
        
	tmpi(:)=meltprev(:)
	call getvara ('meltprev',iou, numiE, ib, ic, tmpi, c1, c0)
        meltprev(:)=tmpi(:)
# endif                           
        deallocate(tmpi)
      endif
# endif
!-----------------------------------------------------------------------
!     read 3d data (x,y,t)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = jle-jls+1
      ib(3) = ntrec
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)
      do n=1,nat
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        var1 = 'at1_'//trim(a3)
        var2 = 'at2_'//trim(a3)
        if (trim(mapat(n)) .eq. 'sat') then
          if (inqvardef('slat1', iou)) var1 = 'slat1'
          if (inqvardef('slat2', iou)) var2 = 'slat2'
        elseif (trim(mapat(n)) .eq. 'shum') then
          if (inqvardef('shum1', iou)) var1 = 'shum1'
          if (inqvardef('shum2', iou)) var2 = 'shum2'
        elseif (trim(mapat(n)) .eq. 'co2') then
          if (inqvardef('co21', iou)) var1 = 'co21'
          if (inqvardef('co22', iou)) var2 = 'co22'
        endif
        tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,1,n)
        call getvara(trim(var1), iou, ln, ib, ic, tmpij, c1, c0)
        at(ils:ile,jls:jle,1,n) = tmpij(ils:ile,jls:jle)
        tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,2,n)
        call getvara(trim(var2), iou, ln, ib, ic, tmpij, c1, c0)
        at(ils:ile,jls:jle,2,n) = tmpij(ils:ile,jls:jle)
      enddo
      tmpij(ils:ile,jls:jle) = rh(ils:ile,jls:jle)
      call getvara ('rh', iou, ln, ib, ic, tmpij, c1, c0)
      rh(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = precip(ils:ile,jls:jle)
      call getvara ('precip', iou, ln, ib, ic, tmpij, c1, c0)
      precip(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isst)
      call getvara ('sbc_sst', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isst) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isss)
      call getvara ('sbc_sss', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isss) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ihflx)
      call getvara ('sbc_hflx', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,ihflx) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isflx)
      call getvara ('sbc_sflx', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isflx) = tmpij(ils:ile,jls:jle)
# if defined O_carbon
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdic)
      call getvara ('sbc_ssdic', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issdic) = tmpij(ils:ile,jls:jle)
#  if defined O_carbon_14
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issc14)
      call getvara ('sbc_ssc14', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issc14) = tmpij(ils:ile,jls:jle)
#  endif
# endif
# if defined O_npzd_o2
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isso2)
      call getvara ('sbc_sso2', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isso2) = tmpij(ils:ile,jls:jle)
# endif
# if defined O_npzd_alk
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issalk)
      call getvara ('sbc_ssalk', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issalk) = tmpij(ils:ile,jls:jle)
# endif

# if defined O_npzd_nitrogen
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issno3)
      call getvara ('sbc_ssno3', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,issno3) = tmpij(ils:ile,jls:jle)
# endif

!# if defined O_ism || defined O_landice_data_transient
!      if (use_ice_data) then
!        tmpij(ils:ile,jls:jle) = tmsk(ils:ile,jls:jle)
!         call getvara ('tmsk', iou, ln, ib, ic, tmpij, c1, c0)
!         tmsk(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
!         tmpij(ils:ile,jls:jle) = hicel(ils:ile,jls:jle,2)
!         call getvara ('hicel', iou, ln, ib, ic, tmpij, c1, c0)
!         hicel(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)
!         tmpij(ils:ile,jls:jle) = aicel(ils:ile,jls:jle,2)
!         call getvara ('aicel', iou, ln, ib, ic, tmpij, c1, c0)
!         aicel(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)
!       endif
!# endif

# if defined O_embm_adiff || defined O_embm_awind
      tmpij(ils:ile,jls:jle) = rtbar(ils:ile,jls:jle)
      call getvara ('rtbar', iou, ln, ib, ic, tmpij, c1, c0)
      rtbar(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = atbar(ils:ile,jls:jle)
      call getvara ('atbar', iou, ln, ib, ic, tmpij, c1, c0)
      atbar(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
# endif
# if defined O_embm_awind
      tmpij(ils:ile,jls:jle) = awx(ils:ile,jls:jle)
      call getvara ('awx', iou, ln, ib, ic, tmpij, c1, c0)
      awx(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = awy(ils:ile,jls:jle)
      call getvara ('awy', iou, ln, ib, ic, tmpij, c1, c0)
      awy(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
# endif
      tmpij(ils:ile,jls:jle) = soilm(ils:ile,jls:jle,1)
      call getvara ('soilm1', iou, ln, ib, ic, tmpij, c1, c0)
      soilm(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = soilm(ils:ile,jls:jle,2)
      call getvara ('soilm2', iou, ln, ib, ic, tmpij, c1, c0)
      soilm(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = surf(ils:ile,jls:jle)
      call getvara ('surf', iou, ln, ib, ic, tmpij, c1, c0)
      surf(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
# if defined O_ice
      tmpij(ils:ile,jls:jle) = ticeo(ils:ile,jls:jle)
      call getvara ('ticeo', iou, ln, ib, ic, tmpij, c1, c0)
      ticeo(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = ticel(ils:ile,jls:jle)
      call getvara ('ticel', iou, ln, ib, ic, tmpij, c1, c0)
      ticel(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)      
      
      tmpij(ils:ile,jls:jle) = hice(ils:ile,jls:jle,1)
      call getvara ('hice1', iou, ln, ib, ic, tmpij, c1, c0)
      hice(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = hice(ils:ile,jls:jle,2)
      call getvara ('hice2', iou, ln, ib, ic, tmpij, c1, c0)
      hice(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)

      tmpij(ils:ile,jls:jle) = aiceocn(ils:ile,jls:jle,1)
      call getvara ('aiceocn1', iou, ln, ib, ic, tmpij, c1, c0)
      aiceocn(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = aiceocn(ils:ile,jls:jle,2)
      call getvara ('aiceocn2', iou, ln, ib, ic, tmpij, c1, c0)
      aiceocn(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)
# if !defined O_ism      
      tmpij(ils:ile,jls:jle) = aicelan(ils:ile,jls:jle,1)
      call getvara ('aicelan1', iou, ln, ib, ic, tmpij, c1, c0)
      aicelan(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = aicelan(ils:ile,jls:jle,2)
      call getvara ('aicelan2', iou, ln, ib, ic, tmpij, c1, c0)
      aicelan(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle) 
# endif           
# if !defined O_icefree_init            
      tmpij(ils:ile,jls:jle) = hsnol(ils:ile,jls:jle,1)
      call getvara ('hsnol1', iou, ln, ib, ic, tmpij, c1, c0)
      hsnol(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = hsnol(ils:ile,jls:jle,2)
      call getvara ('hsnol2', iou, ln, ib, ic, tmpij, c1, c0)
      hsnol(ils:ile,jls:jle,2)= tmpij(ils:ile,jls:jle)
# endif      
      tmpij(ils:ile,jls:jle) = hsnoo(ils:ile,jls:jle,1)
      call getvara ('hsnoo1', iou, ln, ib, ic, tmpij, c1, c0)
      hsnoo(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = hsnoo(ils:ile,jls:jle,2)
      call getvara ('hsnoo2', iou, ln, ib, ic, tmpij, c1, c0)
      hsnoo(ils:ile,jls:jle,2)= tmpij(ils:ile,jls:jle) 

#  if defined O_ism      
      tmpij(ils:ile,jls:jle) = nmbal(ils:ile,jls:jle,1)
      call getvara ('nmbal1', iou, ln, ib, ic, tmpij, c1, c0)
      nmbal(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = nmbal(ils:ile,jls:jle,2)
      call getvara ('nmbal2', iou, ln, ib, ic, tmpij, c1, c0)
      nmbal(ils:ile,jls:jle,2)= tmpij(ils:ile,jls:jle) 
      tmpij(ils:ile,jls:jle) = strtiE(ils:ile,jls:jle)
      call getvara ('strtiE', iou, ln, ib, ic, tmpij, c1, c0)
      strtiE(ils:ile,jls:jle)= nint(tmpij(ils:ile,jls:jle))
      tmpij(ils:ile,jls:jle) = endiE(ils:ile,jls:jle)
      call getvara ('endiE', iou, ln, ib, ic, tmpij, c1, c0)
      endiE(ils:ile,jls:jle)= nint(tmpij(ils:ile,jls:jle))        
      tmpij(ils:ile,jls:jle) = fracblis(ils:ile,jls:jle)
      call getvara ('fracblis', iou, ln, ib, ic, tmpij, c1, c0)
      fracblis(ils:ile,jls:jle)= tmpij(ils:ile,jls:jle)       
      tmpij(ils:ile,jls:jle) = fraco(ils:ile,jls:jle)      
      call getvara ('fraco', iou, ln, ib, ic, tmpij, c1, c0)
      fraco(ils:ile,jls:jle)= tmpij(ils:ile,jls:jle)      
      tmpij(ils:ile,jls:jle) = fraci(ils:ile,jls:jle)      
      call getvara ('fraci', iou, ln, ib, ic, tmpij, c1, c0)
      fraci(ils:ile,jls:jle)= tmpij(ils:ile,jls:jle)       
      tmpij(ils:ile,jls:jle) = fracl(ils:ile,jls:jle)      
      call getvara ('fracl', iou, ln, ib, ic, tmpij, c1, c0)
      fracl(ils:ile,jls:jle)= tmpij(ils:ile,jls:jle)        
      tmpij(ils:ile,jls:jle) = fracbl(ils:ile,jls:jle)      
      call getvara ('fracbl', iou, ln, ib, ic, tmpij, c1, c0)
      fracbl(ils:ile,jls:jle)= tmpij(ils:ile,jls:jle)       
      tmpij(ils:ile,jls:jle) = cir(ils:ile,jls:jle)      
      call getvara ('cir', iou, ln, ib, ic, tmpij, c1, c0)
      cir(ils:ile,jls:jle)= tmpij(ils:ile,jls:jle)           
      tmpij(ils:ile,jls:jle) = icemsk(ils:ile,jls:jle)      
      call getvara ('icemsk', iou, ln, ib, ic, tmpij, c1, c0)
      icemsk(ils:ile,jls:jle)= tmpij(ils:ile,jls:jle)       
      tmpij(ils:ile,jls:jle) = avgsno(ils:ile,jls:jle)      
      call getvara ('avgsno', iou, ln, ib, ic, tmpij, c1, c0)
      avgsno(ils:ile,jls:jle)= tmpij(ils:ile,jls:jle)  	 
      tmpij(ils:ile,jls:jle) = avgmbal(ils:ile,jls:jle)      
      call getvara ('avgmbal', iou, ln, ib, ic, tmpij, c1, c0)
      avgmbal(ils:ile,jls:jle)= tmpij(ils:ile,jls:jle)       
      tmpij(ils:ile,jls:jle) = tmsk(ils:ile,jls:jle)      
      call getvara ('tmsk', iou, ln, ib, ic, tmpij, c1, c0)
      tmsk(ils:ile,jls:jle)= tmpij(ils:ile,jls:jle)  	 
      tmpij(ils:ile,jls:jle) = elev(ils:ile,jls:jle)      
      call getvara ('elev', iou, ln, ib, ic, tmpij, c1, c0)
      elev(ils:ile,jls:jle)= tmpij(ils:ile,jls:jle)     
#  endif      
# endif
# if defined O_ice_cpts && defined O_ice
      do n=1,ncat
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        tmpij(ils:ile,jls:jle) = hseff(ils:ile,jls:jle,1,n)
        call getvara ('hseff1_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        hseff(ils:ile,jls:jle,1,n) = tmpij(ils:ile,jls:jle)
        tmpij(ils:ile,jls:jle) = hseff(ils:ile,jls:jle,2,n)
        call getvara ('hseff2_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        hseff(ils:ile,jls:jle,2,n) = tmpij(ils:ile,jls:jle)
        tmpij(ils:ile,jls:jle) = heff(ils:ile,jls:jle,1,n)
        call getvara ('heff1_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        heff(ils:ile,jls:jle,1,n) = tmpij(ils:ile,jls:jle)
        tmpij(ils:ile,jls:jle) = heff(ils:ile,jls:jle,2,n)
        call getvara ('heff2_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        heff(ils:ile,jls:jle,2,n) = tmpij(ils:ile,jls:jle)
        tmpij(ils:ile,jls:jle) = A(ils:ile,jls:jle,1,n)
        call getvara ('A1_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
        A(ils:ile,jls:jle,1,n) = tmpij(ils:ile,jls:jle)
        tmpij(ils:ile,jls:jle) = A(ils:ile,jls:jle,2,n)
        call getvara ('A2_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
        A(ils:ile,jls:jle,2,n) = tmpij(ils:ile,jls:jle)
        tmpij(ils:ile,jls:jle) = ts(ils:ile,jls:jle,1,n)
        call getvara ('ts1_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
        ts(ils:ile,jls:jle,1,n) = tmpij(ils:ile,jls:jle)
        tmpij(ils:ile,jls:jle) = ts(ils:ile,jls:jle,2,n)
        call getvara ('ts2_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
        ts(ils:ile,jls:jle,2,n) = tmpij(ils:ile,jls:jle)
      enddo
      do n=1,ntilay
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        tmpij(ils:ile,jls:jle) = E(ils:ile,jls:jle,1,n)
        call getvara ('E1_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
        E(ils:ile,jls:jle,1,n) = tmpij(ils:ile,jls:jle)
        tmpij(ils:ile,jls:jle) = E(ils:ile,jls:jle,2,n)
        call getvara ('E2_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
        E(ils:ile,jls:jle,2,n) = tmpij(ils:ile,jls:jle)
      enddo
# endif
# if defined O_ice_cpts_roth_press && defined O_ice_cpts && defined O_ice
      tmpij(ils:ile,jls:jle) = strength(ils:ile,jls:jle,1)
      call getvara ('strength1', iou, ln, ib, ic, tmpij, c1, c0)
      strength(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = strength(ils:ile,jls:jle,2)
      call getvara ('strength2', iou, ln, ib, ic, tmpij, c1, c0)
      strength(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)
# endif
# if defined O_ice_evp && defined O_ice
      tmpij(ils:ile,jls:jle) = uice(ils:ile,jls:jle)
      call getvara ('uice', iou, ln, ib, ic, tmpij, c1, c0)
      uice(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = vice(ils:ile,jls:jle)
      call getvara ('vice', iou, ln, ib, ic, tmpij, c1, c0)
      vice(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isu)
      call getvara ('su', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isu) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isv)
      call getvara ('sv', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isv) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,igu)
      call getvara ('gu', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,igu) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,igv)
      call getvara ('gv', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,igv) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig11n(ils:ile,jls:jle)
      call getvara ('sig11n', iou, ln, ib, ic, tmpij, c1, c0)
      sig11n(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig11s(ils:ile,jls:jle)
      call getvara ('sig11s', iou, ln, ib, ic, tmpij, c1, c0)
      sig11s(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig11e(ils:ile,jls:jle)
      call getvara ('sig11e', iou, ln, ib, ic, tmpij, c1, c0)
      sig11e(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig11w(ils:ile,jls:jle)
      call getvara ('sig11w', iou, ln, ib, ic, tmpij, c1, c0)
      sig11w(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig22n(ils:ile,jls:jle)
      call getvara ('sig22n', iou, ln, ib, ic, tmpij, c1, c0)
      sig22n(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig22s(ils:ile,jls:jle)
      call getvara ('sig22s', iou, ln, ib, ic, tmpij, c1, c0)
      sig22s(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig22e(ils:ile,jls:jle)
      call getvara ('sig22e', iou, ln, ib, ic, tmpij, c1, c0)
      sig22e(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig22w(ils:ile,jls:jle)
      call getvara ('sig22w', iou, ln, ib, ic, tmpij, c1, c0)
      sig22w(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig12n(ils:ile,jls:jle)
      call getvara ('sig12n', iou, ln, ib, ic, tmpij, c1, c0)
      sig12n(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig12s(ils:ile,jls:jle)
      call getvara ('sig12s', iou, ln, ib, ic, tmpij, c1, c0)
      sig12s(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig12e(ils:ile,jls:jle)
      call getvara ('sig12e', iou, ln, ib, ic, tmpij, c1, c0)
      sig12e(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sig12w(ils:ile,jls:jle)
      call getvara ('sig12w', iou, ln, ib, ic, tmpij, c1, c0)
      sig12w(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
# endif
# if defined O_sulphate_data_transient
      tmpij(ils:ile,jls:jle) = sulph(ils:ile,jls:jle,2)
      call getvara ('sulph', iou, ln, ib, ic, tmpij, c1, c0)
      sulph(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)
# endif

      call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      print*, '=> Atm restart read from ',trim(fname),' on ', nstamp

      deallocate ( tmpij )

      return
      end

      subroutine embm_rest_def (fname)
!=======================================================================
!     definition routine for atmospheric restarts

!   inputs:
!     fname = file name
!=======================================================================

      use subgrid

      implicit none

      character(*) :: fname
      character(3) :: a3
      
      integer iou, n, igs, ige, ig, jgs, jge, jg, it(10), iu(10)
      integer id_time, id_xt, id_xu, id_yt, id_yu, id_xt_e, id_xu_e
      integer id_yt_e, id_yu_e, id_track

      real c1e20

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "atm.h"
# if defined O_ice_cpts && defined O_ice
      include "cpts.h"
# endif
      include "ice.h"
# if defined O_ice_evp && defined O_ice
      include "evp.h"
# endif
      include "iounit.h"
      include "tmngr.h"
# if defined O_ism
      include "ism.h"
      integer id_nisb, id_sgdim, id_numiE
# endif
      
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

!-----------------------------------------------------------------------
!     set global write domain size
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     start definitions
!-----------------------------------------------------------------------
      call redef (iou)

!-----------------------------------------------------------------------
!     write global attributes
!-----------------------------------------------------------------------
      call putatttext (iou, 'global', 'Conventions', 'CF-1.0')
      call putatttext (iou, 'global', 'experiment_name', expnam)
      call putatttext (iou, 'global', 'run_stamp', runstamp)

!-----------------------------------------------------------------------
!     define dimensions
!-----------------------------------------------------------------------
      call defdim ('time', iou, 0, id_time)
      call defdim ('xt', iou, ig, id_xt)
      call defdim ('yt', iou, jg, id_yt)
      call defdim ('xu', iou, ig, id_xu)
      call defdim ('yu', iou, jg, id_yu)
      call defdim ('xt_edges', iou, ig+1, id_xt_e)
      call defdim ('yt_edges', iou, jg+1, id_yt_e)
      call defdim ('xu_edges', iou, ig+1, id_xu_e)
      call defdim ('yu_edges', iou, jg+1, id_yu_e)
      call defdim ('track', iou, mtrack, id_track)
#if defined O_ism
      call defdim ('nisb', iou, nisb, id_nisb)
      call defdim ('numiE', iou, 1, id_numiE)            
      if (numiE .gt. 0) then             
      call defdim ('sgdim', iou, numiE, id_sgdim)
      endif
#endif
      
!-----------------------------------------------------------------------
!     define 1d data (t)
!-----------------------------------------------------------------------
      it(1) = id_time
      call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &, 'time', 'time', 'year since 0000-01-01 00:00:00')
      call putatttext (iou, 'time', 'calendar', calendar)
      call defvar ('nats', iou, 1, it, c0, c0, ' ', 'D'
     &, 'nats', ' ',' ')
      call defvar ('dayoyr', iou, 1, it, c0, c0, ' ', 'D'
     &, 'dayoyr', ' ',' ')
      call defvar ('itt', iou, 1, it, c0, c0, ' ', 'D'
     &, 'itt', ' ',' ')
      call defvar ('irstdy', iou, 1, it, c0, c0, ' ', 'D'
     &, 'irstdy', ' ',' ')
      call defvar ('msrsdy', iou, 1, it, c0, c0, ' ', 'D'
     &, 'msrsdy', ' ',' ')
      call defvar ('totaltime', iou, 1, it, c0, c0, ' ', 'D'
     &, 'totaltime', ' ',' ')
      call defvar ('year', iou, 1, it, c0, c0, ' ', 'D'
     &, 'year', ' ',' ')
      call defvar ('month', iou, 1, it, c0, c0, ' ', 'D'
     &, 'month', ' ',' ')
      call defvar ('day', iou, 1, it, c0, c0, ' ', 'D'
     &, 'day', ' ',' ')
      call defvar ('hour', iou, 1, it, c0, c0, ' ', 'D'
     &, 'hour', ' ',' ')
      call defvar ('minute', iou, 1, it, c0, c0, ' ', 'D'
     &, 'minute', ' ',' ')
      call defvar ('second', iou, 1, it, c0, c0, ' ', 'D'
     &, 'second', ' ',' ')
      call defvar ('co2ccn', iou, 1, it, c0, c0, ' ', 'D'
     &, 'co2ccn', ' ',' ')
# if defined O_carbon_14
      call defvar ('c14ccn', iou, 1, it, c0, c0, ' ', 'D'
     &, 'c14ccn', ' ',' ')
# endif
# if defined O_co2emit_track_co2
      call defvar ('itrack_co2', iou, 1, it, c0, c0, ' ', 'D'
     &, 'itrack_co2', ' ',' ')
      call defvar ('ntrack_co2', iou, 1, it, c0, c0, ' ', 'D'
     &, 'ntrack_co2', ' ',' ')
      it(1) = id_track
      call defvar ('track_co2', iou, 1, it, c0, c0, ' ', 'D'
     &, 'co2 for tracking average', ' ', ' ')
      it(1) = id_time
# endif
# if defined O_co2emit_track_sat || defined O_embm_vcs
      call defvar ('itrack_sat', iou, 1, it, c0, c0, ' ', 'D'
     &, 'itrack_sat', ' ',' ')
      call defvar ('ntrack_sat', iou, 1, it, c0, c0, ' ', 'D'
     &, 'ntrack_sat', ' ',' ')
      it(1) = id_track
      call defvar ('track_sat', iou, 1, it, c0, c0, ' ', 'D'
     &, 'sat for tracking average', ' ', ' ')
      it(1) = id_time
# endif
# if defined O_ism || defined O_landice_data
      call defvar ('ice_yr', iou, 1, it, c0, c0, ' ', 'D'
     &, 'ice_yr', ' ',' ')
# endif
# if defined O_sealev || defined O_sealev_data
      call defvar ('sealev_yr', iou, 1, it, c0, c0, ' ', 'D'
     &, 'sealev_yr', ' ',' ')
      call defvar ('sealev', iou, 1, it, c0, c0, ' ', 'D'
     &, 'sealev', ' ',' ')
# endif

!-----------------------------------------------------------------------
!     define 1d data (x, y or z)
!-----------------------------------------------------------------------
      it(1) = id_xt
      call defvar ('xt', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude of the t grid', 'grid_longitude', 'degrees_east')
      it(1) = id_yt
      call defvar ('yt', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude of the t grid', 'grid_latitude', 'degrees_north')
      it(1) = id_xu
      call defvar ('xu', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude of the u grid', 'grid_longitude', 'degrees_east')
      it(1) = id_yu
      call defvar ('yu', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude of the u grid', 'grid_latitude', 'degrees_north')
      it(1) = id_xt_e
      call defvar ('xt_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'longitude of t grid edges', ' ', 'degrees')
      it(1) = id_yt_e
      call defvar ('yt_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'latitude of t grid edges', ' ', 'degrees')
      it(1) = id_xu_e
      call defvar ('xu_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'longitude of u grid edges', ' ', 'degrees')
      it(1) = id_yu_e
      call defvar ('yu_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'latitude of u grid edges', ' ', 'degrees')

#if defined O_ism
!-----------------------------------------------------------------------
!     define 1d ice sheet basin data (basinnumber) and ice sheet subgrid data structure
!-----------------------------------------------------------------------
      it(1) = id_nisb
      call defvar ('dicevol', iou, 1, it, -1.e20,1.e20,' ','D'
     &, 'ice drainage volume', 'dicevol', 'm^3')

      call defvar ('diceheat',iou,1,it,-1.e20,1.e20,' ','D'
     &, 'ice drainage heat', 'diceheat', 'J')

      call defvar ('isflxm',iou,1,it,-1.e20,1.e20,' ','D'
     &, 'ice drainage moisture flux', 'isflxm', 'cm^3/s')

      call defvar ('isflxh',iou,1,it,-1.e20,1.e20,' ','D'
     &, 'ice drainage heat flux', 'isflxh', 'J/s')
     
      it(1) = id_numiE
      call defvar ('numiE', iou, 1, it, c0, c0, ' ', 'D'
     &, 'numiE', ' ',' ')
      call defvar ('dtism', iou, 1, it, c0, c0, ' ', 'D'
     &, 'dtism', ' ',' ')  
      call defvar ('rtag', iou, 1, it, c0, c0, ' ', 'D'
     &, 'rtag', ' ',' ')    
      call defvar ('L_accummassN', iou, 1, it, c0, c0, ' ', 'D'
     &, 'L_accummassN', ' ',' ')    
      call defvar ('L_accummassS', iou, 1, it, c0, c0, ' ', 'D'
     &, 'L_accummassS', ' ',' ')
      call defvar ('L_meltmassN', iou, 1, it, c0, c0, ' ', 'D'
     &, 'L_meltmassN', ' ',' ') 
      call defvar ('L_meltmassS', iou, 1, it, c0, c0, ' ', 'D'
     &, 'L_meltmassS', ' ',' ')                   
      call defvar ('L_iceareaN', iou, 1, it, c0, c0, ' ', 'D'
     &, 'L_iceareaN', ' ',' ')     
      call defvar ('L_iceareaS', iou, 1, it, c0, c0, ' ', 'D'
     &, 'L_iceareaS', ' ',' ')
      call defvar ('L_icesmbN', iou, 1, it, c0, c0, ' ', 'D'
     &, 'L_icesmbN', ' ',' ')     
      call defvar ('L_icesmbS', iou, 1, it, c0, c0, ' ', 'D'
     &, 'L_icesmbS', ' ',' ')     
      call defvar ('L_icevolN', iou, 1, it, c0, c0, ' ', 'D'
     &, 'L_icevolN', ' ',' ')     
      call defvar ('L_icevolS', iou, 1, it, c0, c0, ' ', 'D'
     &, 'L_icevolS', ' ',' ')     
      call defvar ('L_meltextentN', iou, 1, it, c0, c0, ' ', 'D'
     &, 'L_meltextentN', ' ',' ')     
      call defvar ('L_meltextentS', iou, 1, it, c0, c0, ' ', 'D'
     &, 'L_meltextentS', ' ',' ')  
      if (numiE .gt. 0.) then
      it(1) = id_sgdim     
      call defvar ('eleviE',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'eleviE', 'eleviE', 'cm') 
      call defvar ('fraciE',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'fraciE', 'fraciE', '%')  
      call defvar ('iariE',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'iariE', 'iariE', 'm^2')           
      call defvar ('isniE',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'isniE', 'isniE', 'n')      
      call defvar ('hsnowiE1',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'hsnowiE1', 'hsnowiE1', 'cm') 
      call defvar ('hsnowiE2',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'hsnowiE2', 'hsnowiE2', 'cm')      
      call defvar ('ticeiE',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'ticeiE', 'ticeiE', 'C') 
      call defvar ('ismticeiE',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'ismticeiE', 'ismticeiE', 'C')      
      call defvar ('ithkiE',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'ithkiE', 'ithkiE', 'm')          
      call defvar ('mbalmelt1',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'mbalmelt1', 'mbalmelt1', 'cm') 
      call defvar ('mbalmelt2',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'mbalmelt2', 'mbalmelt2', 'cm')      
      call defvar ('mbalaccum1',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'mbalaccum1', 'mbalaccum1', 'cm')      
      call defvar ('mbalaccum2',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'mbalaccum2', 'mbalaccum2', 'cm')
       call defvar ('precipiE',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'precipiE', 'precipiE', '')  
       call defvar ('satiE',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'satiE', 'satiE', 'C')      
       call defvar ('albedoiE',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'albedoiE', 'albedoiE', '')      
       call defvar ('mextiE',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'mextiE', 'mextiE', '')      
       call defvar ('mduriE',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'mduriE', 'mduriE', '')
# if defined O_refreeze
       call defvar ('hnewsno',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'hnewsno', 'hnewsno', '') 
       call defvar ('counter',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'counter', 'counter', '')    
       call defvar ('htemp',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'htemp', 'htemp', '')	      
       call defvar ('hreftot',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'hreftot', 'hreftot', '')    
       call defvar ('meltprev',iou,1,it,-1.e30,1.e30,' ','D'
     &, 'meltprev', 'meltprev', '')    
# endif     
           
      endif
#endif   
!-----------------------------------------------------------------------
!     define 3d data (x,y,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      it(3) = id_time
      iu(3) = id_time
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          call defvar ('slat1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &,     'sea level atmospheric temperature at tau', ' ', ' ')
          call defvar ('slat2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &,     'sea level atmospheric temperature at tau+1', ' ', ' ')
        elseif (trim(mapat(n)) .eq. 'shum') then
          call defvar ('shum1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &,     'atmospheric surface specific humidity at tau', ' ', ' ')
          call defvar ('shum2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &,     'atmospheric surface specific humidity at tau+1', ' ', ' ')
        elseif (trim(mapat(n)) .eq. 'co2') then
          call defvar ('co21', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &,     'atmospheric co2 at tau', ' ', ' ')
          call defvar ('co22', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &,     'atmospheric co2 at tau+1', ' ', ' ')
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          call defvar ('at1_'//trim(a3), iou , 3, it, -c1e20, c1e20, ' '
     &,     'D', 'at1_'//trim(a3), ' ', ' ')
          call defvar ('at2_'//trim(a3), iou , 3, it, -c1e20, c1e20, ' '
     &,     'D', 'at2_'//trim(a3), ' ', ' ')
        endif
      enddo
      call defvar ('rh', iou, 3, it,  -c1e20, c1e20, ' ', 'D'
     &, 'rh', ' ', ' ')
      call defvar ('precip', iou, 3, it,  -c1e20, c1e20, ' ', 'D'
     &, 'precip', ' ', ' ')
      call defvar ('sbc_sst', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sst', ' ', ' ')
      call defvar ('sbc_sss', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sss', ' ', ' ')
      call defvar ('sbc_hflx', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sst', ' ', ' ')
      call defvar ('sbc_sflx', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sss', ' ', ' ')
# if defined O_carbon
      call defvar ('sbc_ssdic', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssdic', ' ', ' ')
#  if defined O_carbon_14
      call defvar ('sbc_ssc14', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssc14', ' ', ' ')
#  endif
# endif
# if defined O_npzd_o2
      call defvar ('sbc_sso2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sso2', ' ', ' ')
# endif
# if defined O_npzd_alk
      call defvar ('sbc_ssalk', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssalk', ' ', ' ')
# endif

# if defined O_npzd_nitrogen
      call defvar ('sbc_ssno3', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ssno3', ' ', ' ')
# endif
!# if defined O_ism || defined O_landice_data
!      call defvar ('tmsk', iou, 3, it, -c1e20, c1e20, ' ', 'D'
!     &, 'tmsk', ' ', ' ')
!      call defvar ('hicel', iou, 3, it, -c1e20, c1e20, ' ', 'D'
!     &, 'hicel', ' ', ' ')
!      call defvar ('aicel', iou, 3, it, -c1e20, c1e20, ' ', 'D'
!     &, 'aicel', ' ', ' ')
!# endif
# if defined O_embm_adiff || defined O_embm_awind
      call defvar ('rtbar', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'rtbar', ' ', ' ')
      call defvar ('atbar', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'atbar', ' ', ' ')
# endif
# if defined O_embm_awind
      call defvar ('awx', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'awx', ' ', ' ')
      call defvar ('awy', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'awy', ' ', ' ')
# endif
      call defvar ('soilm1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'soilm1', ' ', ' ')
      call defvar ('soilm2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'soilm2', ' ', ' ')
      call defvar ('surf', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'surf', ' ', ' ')
# if defined O_ice
      call defvar ('ticel', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ticel', ' ', ' ')
      call defvar ('ticeo', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'ticeo', ' ', ' ')     
      call defvar ('hice1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'hice1', ' ', ' ')
      call defvar ('hice2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'hice2', ' ', ' ')

      call defvar ('aiceocn1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'aiceocn1', ' ', ' ')
      call defvar ('aiceocn2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'aiceocn2', ' ', ' ')
# if !defined O_ism     
      call defvar ('aicelan1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'aicelan1', ' ', ' ')
      call defvar ('aicelan2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'aicelan2', ' ', ' ')     
# endif
      call defvar ('hsnol1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'hsnol1', ' ', ' ')
      call defvar ('hsnol2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'hsnol2', ' ', ' ')
      call defvar ('hsnoo1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'hsnoo1', ' ', ' ')
      call defvar ('hsnoo2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'hsnoo2', ' ', ' ') 
#  if defined O_ism
      call defvar ('nmbal1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'nmbal1', ' ', ' ')
      call defvar ('nmbal2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'nmbal2', ' ', ' ')
      call defvar ('strtiE', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'strtiE', ' ', ' ') 
      call defvar ('endiE', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'endiE', ' ', ' ')    

      call defvar ('fracblis', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'fracblis', ' ', ' ')         
      call defvar ('fraco', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'fraco', ' ', ' ')      
      call defvar ('fraci', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'fraci', ' ', ' ')         
      call defvar ('fracl', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'fracl', ' ', ' ')  
      call defvar ('fracbl', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'fracbl', ' ', ' ')         
      call defvar ('cir', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'cir', ' ', ' ')  
      call defvar ('icemsk', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'icemsk', ' ', ' ')         
      call defvar ('avgsno', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'avgsno', ' ', ' ')  
      call defvar ('avgmbal', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'avgmbal', ' ', ' ')         
       call defvar ('tmsk', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'tmsk', ' ', ' ')  
      call defvar ('elev', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'elev', ' ', ' ')              
#  endif     
               
# endif
# if defined O_ice_cpts && defined O_ice
      do n=1,ncat
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        call defvar ('hseff1_'//trim(a3), iou, 3, it, -c1e20, c1e20
     &,   ' ', 'D', 'hseff1_'//trim(a3), ' ', ' ')
        call defvar ('hseff2_'//trim(a3), iou, 3, it, -c1e20, c1e20
     &,   ' ', 'D', 'hseff2_'//trim(a3), ' ', ' ')
        call defvar ('heff1_'//trim(a3), iou, 3, it, -c1e20, c1e20
     &,   ' ', 'D', 'heff1_'//trim(a3), ' ', ' ')
        call defvar ('heff2_'//trim(a3), iou, 3, it, -c1e20, c1e20
     &,   ' ', 'D', 'heff2_'//trim(a3), ' ', ' ')
        call defvar ('A1_'//trim(a3), iou, 3, it, -c1e20, c1e20
     &,   ' ', 'D', 'A1_'//trim(a3), ' ', ' ')
        call defvar ('A2_'//trim(a3), iou, 3, it, -c1e20, c1e20
     &,   ' ', 'D', 'A2_'//trim(a3), ' ', ' ')
        call defvar ('ts1_'//trim(a3), iou, 3, it, -c1e20, c1e20
     &,   ' ', 'D', 'ts1_'//trim(a3), ' ', ' ')
        call defvar ('ts2_'//trim(a3), iou, 3, it, -c1e20, c1e20
     &,   ' ', 'D', 'ts2_'//trim(a3), ' ', ' ')
      enddo
      do n=1,ntilay
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        call defvar ('E1_'//trim(a3), iou, 3, it, -c1e20, c1e20
     &,   ' ', 'D', 'E1_'//trim(a3), ' ', ' ')
        call defvar ('E2_'//trim(a3), iou, 3, it, -c1e20, c1e20
     &,   ' ', 'D', 'E2_'//trim(a3), ' ', ' ')
      enddo
# endif
# if defined O_ice_cpts_roth_press && defined O_ice_cpts && defined O_ice
      call defvar ('strength1', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'strength1', ' ', ' ')
      call defvar ('strength2', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'strength2', ' ', ' ')
# endif
# if defined O_ice_evp && defined O_ice
      call defvar ('uice', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'uice', ' ', ' ')
      call defvar ('vice', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'vice', ' ', ' ')
      call defvar ('su', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'su', ' ', ' ')
      call defvar ('sv', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'sv', ' ', ' ')
      call defvar ('gu', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'gu', ' ', ' ')
      call defvar ('gv', iou, 3, iu, -c1e20, c1e20, ' ', 'D'
     &, 'gv', ' ', ' ')
      call defvar ('sig11n', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig11n', ' ', ' ')
      call defvar ('sig11s', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig11s', ' ', ' ')
      call defvar ('sig11e', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig11e', ' ', ' ')
      call defvar ('sig11w', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig11w', ' ', ' ')
      call defvar ('sig22n', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig22n', ' ', ' ')
      call defvar ('sig22s', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig22s', ' ', ' ')
      call defvar ('sig22e', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig22e', ' ', ' ')
      call defvar ('sig22w', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig22w', ' ', ' ')
      call defvar ('sig12n', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig12n', ' ', ' ')
      call defvar ('sig12s', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig12s', ' ', ' ')
      call defvar ('sig12e', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig12e', ' ', ' ')
      call defvar ('sig12w', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sig12w', ' ', ' ')
# endif
# if defined O_sulphate_data_transient
      call defvar ('sulph', iou, 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sulph', ' ', ' ')
# endif
      call enddef (iou)

      return
      end

      subroutine embm_rest_out (fname, ids, ide, jds, jde)
!=======================================================================
!     output routine for atmospheric restarts

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname              = file name
!     ids, ide ...       = start and end index for data domain
!=======================================================================

      use subgrid

      implicit none

      character(*) :: fname
      character(3) :: a3
      character(32) :: nstamp
      character(120) :: var1, var2

      integer i, iou, j, ln, n, ntrec, ids, ide, jds, jde, igs, ige, ig
      integer jgs, jge, jg, ils, ile, jls, jle, ib(10), ic(10)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      logical inqvardef, exists

      real time, tmp, c1e20
      real, allocatable :: tmpij(:,:)
      real, allocatable :: tmpi(:), tmpj(:)
      real, allocatable :: tmpie(:), tmpje(:)

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "atm.h"
      include "cembm.h"
      include "coord.h"
      include "csbc.h"
      include "grdvar.h"
# if defined O_ice_cpts && defined O_ice
      include "cpts.h"
# endif
      include "ice.h"
# if defined O_ice_evp && defined O_ice
      include "evp.h"
# endif
      include "tmngr.h"
      include "iounit.h"
      include "switch.h"
# if defined O_ism
      include "ism.h"
      character(120) g_st
# endif

      real xt_e(imt+1), xu_e(imt+1), yt_e(jmt+1), yu_e(jmt+1)

      c1e20 = 1.e20
      nstamp = stamp

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)
      ntrec = 1

!-----------------------------------------------------------------------
!     set global write domain size
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global write domain)
!-----------------------------------------------------------------------
      ils = max(ids,igs)
      ile = min(ide,ige)
      jls = max(jds,jgs)
      jle = min(jde,jge)

      allocate ( tmpij(ils:ile,jls:jle) )

!-----------------------------------------------------------------------
!     write 1d data (t)
!-----------------------------------------------------------------------
      if (init_time_out) then
        tmp = 0.
        call putvars ('time', iou, ntrec, tmp, c1, c0)
        tmp = 0.
        call putvars ('itt', iou, ntrec, tmp, c1, c0)
        tmp = 0.
        call putvars ('irstdy', iou, ntrec, tmp, c1, c0)
        tmp = 0.
        call putvars ('msrsdy', iou, ntrec, tmp, c1, c0)
        call mkstmp (nstamp, year0, month0, day0, hour0, min0, sec0)
      else
        tmp = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call putvars ('time', iou, ntrec, tmp, c1, c0)
        tmp = itt
        call putvars ('itt', iou, ntrec, tmp, c1, c0)
        tmp = iday(imodeltime)
        call putvars ('irstdy', iou, ntrec, tmp, c1, c0)
        tmp = msday(imodeltime)
        call putvars ('msrsdy', iou, ntrec, tmp, c1, c0)
      endif
      tmp = nats
      call putvars ('nats', iou, ntrec, tmp, c1, c0)
      tmp = dayoyr
      call putvars ('dayoyr', iou, ntrec, tmp, c1, c0)
      call putvars ('totaltime', iou, ntrec, totaltime, c1, c0)
      call rdstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      tmp = nyear
      call putvars ('year', iou, ntrec, tmp, c1, c0)
      tmp = nmonth
      call putvars ('month', iou, ntrec, tmp, c1, c0)
      tmp = nday
      call putvars ('day', iou, ntrec, tmp, c1, c0)
      tmp = nhour
      call putvars ('hour', iou, ntrec, tmp, c1, c0)
      tmp = nmin
      call putvars ('minute', iou, ntrec, tmp, c1, c0)
      tmp = nsec
      call putvars ('second', iou, ntrec, tmp, c1, c0)
      tmp = co2ccn
      call putvars ('co2ccn', iou, ntrec, tmp, c1, c0)
# if defined O_carbon_14
      tmp = c14ccn
      call putvars ('c14ccn', iou, ntrec, tmp, c1, c0)
# endif
# if defined O_co2emit_track_co2
      tmp = itrack_co2
      call putvars ('itrack_co2', iou, ntrec, tmp, c1, c0)
      tmp = ntrack_co2
      call putvars ('ntrack_co2', iou, ntrec, tmp, c1, c0)     
      ib(1) = 1
      ic(1) = mtrack
      call putvara ('track_co2', iou, ig, ib, ic, track_co2, c1, c0)
# endif
# if defined O_co2emit_track_sat || defined O_embm_vcs
      tmp = itrack_sat
      call putvars ('itrack_sat', iou, ntrec, tmp, c1, c0)
      tmp = ntrack_sat
      call putvars ('ntrack_sat', iou, ntrec, tmp, c1, c0)
      ib(1) = 1
      ic(1) = mtrack
      call putvara ('track_sat', iou, ig, ib, ic, track_sat, c1, c0)
# endif
# if defined O_ism || defined O_landice_data
      tmp = ice_yr
      call putvars ('ice_yr', iou, ntrec, tmp, c1, c0)
# endif
# if defined O_sealev || defined O_sealev_data
      tmp = sealev_yr
      call putvars ('sealev_yr', iou, ntrec, tmp, c1, c0)
      tmp = sealev
      call putvars ('sealev', iou, ntrec, tmp, c1, c0)
# endif
# if defined O_ism
!-----------------------------------------------------------------------
!     write 1d data into ice basin heat and water reservoir arrays and ice sheet subgrid array
!-----------------------------------------------------------------------
      allocate (tmpi(nisb))
      ib(1) = 1
      ic(1) = nisb
      
      tmpi(:) = dicevol(:)
      call putvara ('dicevol', iou, nisb, ib, ic, tmpi,c1,c0)
      
      tmpi(:) = diceheat(:)
      call putvara ('diceheat', iou, nisb, ib, ic, tmpi,c1,c0)
      
      tmpi(:) = isflxm(:)
      call putvara ('isflxm', iou, nisb, ib, ic, tmpi,c1,c0)
      
      tmpi(:) = isflxh(:)
      call putvara ('isflxh', iou, nisb, ib, ic, tmpi,c1,c0)
      deallocate(tmpi)
      
      tmp=dtism
      call putvars ('dtism', iou, ntrec, tmp, c1, c0)      
      
      tmp=numiE
      call putvars ('numiE', iou, ntrec, tmp, c1, c0)
      
      tmp=rtag
      call putvars ('rtag', iou, ntrec, tmp, c1, c0) 
      
      tmp=nisaccum
      call putvars ('L_accummassN', iou, ntrec, tmp, c1, c0)      
      
      tmp=sisaccum
      call putvars ('L_accummassS', iou, ntrec, tmp, c1, c0)       
      
      tmp=nismelt
      call putvars ('L_meltmassN', iou, ntrec, tmp, c1, c0)       
      
      tmp=sismelt
      call putvars ('L_meltmassS', iou, ntrec, tmp, c1, c0)       
      
      tmp=nisiar
      call putvars ('L_iceareaN', iou, ntrec, tmp, c1, c0)      
      
      tmp=sisiar
      call putvars ('L_iceareaS', iou, ntrec, tmp, c1, c0)       
      
      tmp=nissmb
      call putvars ('L_icesmbN', iou, ntrec, tmp, c1, c0)       
      
      tmp=sissmb
      call putvars ('L_icesmbS', iou, ntrec, tmp, c1, c0)
      
      tmp=nisvol
      call putvars ('L_icevolN', iou, ntrec, tmp, c1, c0)      
      
      tmp=sisvol
      call putvars ('L_icevolS', iou, ntrec, tmp, c1, c0)      
      
      tmp=nismext      
      call putvars ('L_meltextentN', iou, ntrec, tmp, c1, c0)      
      
      tmp=sismext
      call putvars ('L_meltextentS', iou, ntrec, tmp, c1, c0)              
      
      if (numiE.gt.0.) then
        allocate (tmpi(numiE))
        ic(1) = numiE         
        
	tmpi(:) = eleviE(:)
        call putvara ('eleviE', iou, numiE, ib, ic, tmpi,c1,c0)
        
	tmpi(:) = fraciE(:)
        call putvara ('fraciE', iou, numiE, ib, ic, tmpi,c1,c0)  
        
	tmpi(:) = iariE(:)
        call putvara ('iariE', iou, numiE, ib, ic, tmpi,c1,c0)           
        
	tmpi(:) = real(isniE(:))   
	call putvara ('isniE', iou, numiE, ib, ic, tmpi,c1,c0)      
        
	tmpi(:) = hsnowiE(:,1)
        call putvara ('hsnowiE1', iou, numiE, ib, ic, tmpi,c1,c0)      
        
	tmpi(:) = hsnowiE(:,2)
        call putvara ('hsnowiE2', iou, numiE, ib, ic, tmpi,c1,c0)
        
	tmpi(:) = ticeiE(:)
        call putvara ('ticeiE', iou, numiE, ib, ic, tmpi,c1,c0)      
        
	tmpi(:) = ismticeiE(:)
        call putvara ('ismticeiE', iou, numiE, ib, ic, tmpi,c1,c0)      
        
	tmpi(:) = ithkiE(:)
        call putvara ('ithkiE', iou, numiE, ib, ic, tmpi,c1,c0)        
        
	tmpi(:) = mbalmelt(:,1)
        call putvara ('mbalmelt1', iou, numiE, ib, ic, tmpi,c1,c0)
        
	tmpi(:) = mbalmelt(:,2)
        call putvara ('mbalmelt2', iou, numiE, ib, ic, tmpi,c1,c0)      
        
	tmpi(:) = mbalaccum(:,1)
        call putvara ('mbalaccum1', iou, numiE, ib, ic, tmpi,c1,c0)      
        
	tmpi(:) = mbalaccum(:,2)
        call putvara ('mbalaccum2', iou, numiE, ib, ic, tmpi,c1,c0)  

        tmpi(:) = precipiE(:)        
	call putvara ('precipiE', iou, numiE, ib, ic, tmpi,c1,c0)      
        
	tmpi(:) = satiE(:)
        call putvara ('satiE', iou, numiE, ib, ic, tmpi,c1,c0)     
        
	tmpi(:) = albedoiE(:)
        call putvara ('albedoiE', iou, numiE, ib, ic, tmpi,c1,c0)      
        
	tmpi(:) = mextiE(:)
        call putvara ('mextiE', iou, numiE, ib, ic, tmpi,c1,c0)     
        
	tmpi(:) = mduriE(:)
        call putvara ('mduriE', iou, numiE, ib, ic, tmpi,c1,c0)
      
#  if defined O_refreeze
        tmpi(:) = hnewsno(:)
        call putvara ('hnewsno', iou, numiE, ib, ic, tmpi,c1,c0)
        
	tmpi(:) = counter(:)
        call putvara ('counter', iou, numiE, ib, ic, tmpi,c1,c0)      
        
	tmpi(:) = htemp(:)
        call putvara ('htemp', iou, numiE, ib, ic, tmpi,c1,c0)      
        
	tmpi(:) = hreftot(:)
        call putvara ('hreftot', iou, numiE, ib, ic, tmpi,c1,c0)      
        
	tmpi(:) = meltprev(:)
        call putvara ('meltprev', iou, numiE, ib, ic, tmpi,c1,c0)            
#  endif      

      deallocate(tmpi)      
      endif      
# endif
!-----------------------------------------------------------------------
!     write 1d data (x or y)
!-----------------------------------------------------------------------
      allocate ( tmpi(igs:ige) )
      allocate ( tmpj(jgs:jge) )
      allocate ( tmpie(igs:ige+1) )
      allocate ( tmpje(jgs:jge+1) )

      ib(1) = 1
      ic(1) = ig
      tmpi(igs:ige) = xt(igs:ige)
      call putvara ('xt', iou, ig, ib, ic, tmpi, c1, c0)
      tmpi(igs:ige) = xu(igs:ige)
      call putvara ('xu', iou, ig, ib, ic, tmpi, c1, c0)

      ic(1) = jg
      tmpj(jgs:jge) = yt(jgs:jge)
      call putvara ('yt', iou, jg, ib, ic, tmpj, c1, c0)
      tmpj(jgs:jge) = yu(jgs:jge)
      call putvara ('yu', iou, jg, ib, ic, tmpj, c1, c0)

      ic(1) = ig + 1
      call edge_maker (1, xt_e, xt, dxt, xu, dxu, imt)
      tmpie(igs:ige+1) = xt_e(igs:ige+1)
      call putvara ('xt_edges', iou, ig+1, ib, ic, tmpie, c1, c0)
      call edge_maker (2, xu_e, xt, dxt, xu, dxu, imt)
      tmpie(igs:ige+1) = xu_e(igs:ige+1)
      call putvara ('xu_edges', iou, ig+1, ib, ic, tmpie, c1, c0)

      ic(1) = jg + 1
      call edge_maker (1, yt_e, yt, dyt, yu, dyu, jmt)
      tmpje(jgs:jge+1) = yt_e(jgs:jge+1)
      call putvara ('yt_edges', iou, jg+1, ib, ic, tmpje, c1, c0)
      call edge_maker (2, yu_e, yt, dyt, yu, dyu, jmt)
      tmpje(jgs:jge+1) = yu_e(jgs:jge+1)
      call putvara ('yu_edges', iou, jg+1, ib, ic, tmpje, c1, c0)

      deallocate ( tmpi )
      deallocate ( tmpj )
      deallocate ( tmpie )
      deallocate ( tmpje )

!-----------------------------------------------------------------------
!     write 3d data (x,y,t)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = jle-jls+1
      ib(3) = ntrec
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          var1 = 'slat1'
          var2 = 'slat2'
        elseif (trim(mapat(n)) .eq. 'shum') then
          var1 = 'shum1'
          var2 = 'shum2'
        elseif (trim(mapat(n)) .eq. 'co2') then
          var1 = 'co21'
          var2 = 'co22'
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          var1 = 'at1_'//trim(a3)
          var2 = 'at2_'//trim(a3)
        endif
        tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,1,n)
        call putvara(trim(var1), iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,2,n)
        call putvara(trim(var2), iou, ln, ib, ic, tmpij, c1, c0)
      enddo
      tmpij(ils:ile,jls:jle) = rh(ils:ile,jls:jle)
      call putvara ('rh', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = precip(ils:ile,jls:jle)
      call putvara ('precip', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isst)
      call putvara ('sbc_sst', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isss)
      call putvara ('sbc_sss', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ihflx)
      call putvara ('sbc_hflx', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isflx)
      call putvara ('sbc_sflx', iou, ln, ib, ic, tmpij, c1, c0)
# if defined O_carbon
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issdic)
      call putvara ('sbc_ssdic', iou, ln, ib, ic, tmpij, c1, c0)
#  if defined O_carbon_14
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issc14)
      call putvara ('sbc_ssc14', iou, ln, ib, ic, tmpij, c1, c0)
#  endif
# endif
# if defined O_npzd_o2
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isso2)
      call putvara ('sbc_sso2', iou, ln, ib, ic, tmpij, c1, c0)
# endif
# if defined O_npzd_alk
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issalk)
      call putvara ('sbc_ssalk', iou, ln, ib, ic, tmpij, c1, c0)
# endif

# if defined O_npzd_nitrogen
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,issno3)
      call putvara ('sbc_ssno3', iou, ln, ib, ic, tmpij, c1, c0)
# endif
!# if defined O_ism || defined O_landice_data
!      tmpij(ils:ile,jls:jle) =  tmsk(ils:ile,jls:jle)
!      call putvara ('tmsk', iou, ln, ib, ic, tmpij, c1, c0)
!      tmpij(ils:ile,jls:jle) = hicel(ils:ile,jls:jle,2)
!      call putvara ('hicel', iou, ln, ib, ic, tmpij, c1, c0)
!      tmpij(ils:ile,jls:jle) = aicel(ils:ile,jls:jle,2)
!      call putvara ('aicel', iou, ln, ib, ic, tmpij, c1, c0)
!# endif

# if defined O_embm_adiff || defined O_embm_awind
      tmpij(ils:ile,jls:jle) = rtbar(ils:ile,jls:jle)
      call putvara ('rtbar', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = atbar(ils:ile,jls:jle)
      call putvara ('atbar', iou, ln, ib, ic, tmpij, c1, c0)
# endif
# if defined O_embm_awind
      tmpij(ils:ile,jls:jle) = awx(ils:ile,jls:jle)
      call putvara ('awx', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = awy(ils:ile,jls:jle)
      call putvara ('awy', iou, ln, ib, ic, tmpij, c1, c0)
# endif
      tmpij(ils:ile,jls:jle) = soilm(ils:ile,jls:jle,1)
      call putvara ('soilm1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = soilm(ils:ile,jls:jle,2)
      call putvara ('soilm2', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = surf(ils:ile,jls:jle)
      call putvara ('surf', iou, ln, ib, ic, tmpij, c1, c0)
# if defined O_ice
           
      tmpij(ils:ile,jls:jle) = ticel(ils:ile,jls:jle)
      call putvara ('ticel', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = ticeo(ils:ile,jls:jle)
      call putvara ('ticeo', iou, ln, ib, ic, tmpij, c1, c0)      
      tmpij(ils:ile,jls:jle) = hice(ils:ile,jls:jle,1)      
      call putvara ('hice1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hice(ils:ile,jls:jle,2)
      call putvara ('hice2', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = aiceocn(ils:ile,jls:jle,1)
      call putvara ('aiceocn1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = aiceocn(ils:ile,jls:jle,2)
      call putvara ('aiceocn2', iou, ln, ib, ic, tmpij, c1, c0)
# if !defined O_ism
      tmpij(ils:ile,jls:jle) = aicelan(ils:ile,jls:jle,1)
      call putvara ('aicelan1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = aicelan(ils:ile,jls:jle,2)
      call putvara ('aicelan2', iou, ln, ib, ic, tmpij, c1, c0)      
# endif
      tmpij(ils:ile,jls:jle) = hsnol(ils:ile,jls:jle,1)
      call putvara ('hsnol1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hsnol(ils:ile,jls:jle,2)
      call putvara ('hsnol2', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hsnoo(ils:ile,jls:jle,1)
      call putvara ('hsnoo1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hsnoo(ils:ile,jls:jle,2)
      call putvara ('hsnoo2', iou, ln, ib, ic, tmpij, c1, c0)

#  if defined O_ism      
      tmpij(ils:ile,jls:jle) = nmbal(ils:ile,jls:jle,1)
      call putvara ('nmbal1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = nmbal(ils:ile,jls:jle,2)
      call putvara ('nmbal2', iou, ln, ib, ic, tmpij, c1, c0) 
      tmpij(ils:ile,jls:jle) = real(strtiE(ils:ile,jls:jle))
      call putvara ('strtiE', iou, ln, ib, ic, tmpij, c1, c0) 
      tmpij(ils:ile,jls:jle) = real(endiE(ils:ile,jls:jle))
      call putvara ('endiE', iou, ln, ib, ic, tmpij, c1, c0)
      
      tmpij(ils:ile,jls:jle) = real(fracblis(ils:ile,jls:jle))
      call putvara ('fracblis', iou, ln, ib, ic, tmpij, c1, c0)      
      tmpij(ils:ile,jls:jle) = real(fraco(ils:ile,jls:jle))
      call putvara ('fraco', iou, ln, ib, ic, tmpij, c1, c0)       
       tmpij(ils:ile,jls:jle) = real(fraci(ils:ile,jls:jle))
      call putvara ('fraci', iou, ln, ib, ic, tmpij, c1, c0)      
      tmpij(ils:ile,jls:jle) = real(fracl(ils:ile,jls:jle))
      call putvara ('fracl', iou, ln, ib, ic, tmpij, c1, c0)     
      tmpij(ils:ile,jls:jle) = real(fracbl(ils:ile,jls:jle))
      call putvara ('fracbl', iou, ln, ib, ic, tmpij, c1, c0)      
      tmpij(ils:ile,jls:jle) = real(cir(ils:ile,jls:jle))
      call putvara ('cir', iou, ln, ib, ic, tmpij, c1, c0)       
       tmpij(ils:ile,jls:jle) = real(icemsk(ils:ile,jls:jle))
      call putvara ('icemsk', iou, ln, ib, ic, tmpij, c1, c0)      
      tmpij(ils:ile,jls:jle) = real(avgsno(ils:ile,jls:jle))
      call putvara ('avgsno', iou, ln, ib, ic, tmpij, c1, c0)       
      tmpij(ils:ile,jls:jle) = real(avgmbal(ils:ile,jls:jle))
      call putvara ('avgmbal', iou, ln, ib, ic, tmpij, c1, c0)      
       tmpij(ils:ile,jls:jle) = real(tmsk(ils:ile,jls:jle))
       call putvara ('tmsk', iou, ln, ib, ic, tmpij, c1, c0)       
       tmpij(ils:ile,jls:jle) = real(elev(ils:ile,jls:jle))
      call putvara ('elev', iou, ln, ib, ic, tmpij, c1, c0)             
#  endif            
# endif
# if defined O_ice_cpts && defined O_ice
      do n=1,ncat
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        tmpij(ils:ile,jls:jle) = hseff(ils:ile,jls:jle,1,n)
        call putvara ('hseff1_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        tmpij(ils:ile,jls:jle) = hseff(ils:ile,jls:jle,2,n)
        call putvara ('hseff2_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        tmpij(ils:ile,jls:jle) = heff(ils:ile,jls:jle,1,n)
        call putvara ('heff1_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        tmpij(ils:ile,jls:jle) = heff(ils:ile,jls:jle,2,n)
        call putvara ('heff2_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        tmpij(ils:ile,jls:jle) = A(ils:ile,jls:jle,1,n)
        call putvara ('A1_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = A(ils:ile,jls:jle,2,n)
        call putvara ('A2_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = ts(ils:ile,jls:jle,1,n)
        call putvara ('ts1_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = ts(ils:ile,jls:jle,2,n)
        call putvara ('ts2_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
      enddo
      do n=1,ntilay
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        tmpij(ils:ile,jls:jle) = E(ils:ile,jls:jle,1,n)
        call putvara ('E1_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = E(ils:ile,jls:jle,2,n)
        call putvara ('E2_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)
      enddo
# endif
# if defined O_ice_cpts_roth_press && defined O_ice_cpts && defined O_ice
      tmpij(ils:ile,jls:jle) = strength(ils:ile,jls:jle,1)
      call putvara ('strength1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = strength(ils:ile,jls:jle,2)
      call putvara ('strength2', iou, ln, ib, ic, tmpij, c1, c0)
# endif
# if defined O_ice_evp && defined O_ice
      tmpij(ils:ile,jls:jle) = uice(ils:ile,jls:jle)
      call putvara ('uice', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = vice(ils:ile,jls:jle)
      call putvara ('vice', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isu)
      call putvara ('su', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isv)
      call putvara ('sv', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,igu)
      call putvara ('gu', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,igv)
      call putvara ('gv', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig11n(ils:ile,jls:jle)
      call putvara ('sig11n', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig11s(ils:ile,jls:jle)
      call putvara ('sig11s', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig11e(ils:ile,jls:jle)
      call putvara ('sig11e', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig11w(ils:ile,jls:jle)
      call putvara ('sig11w', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig22n(ils:ile,jls:jle)
      call putvara ('sig22n', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig22s(ils:ile,jls:jle)
      call putvara ('sig22s', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig22e(ils:ile,jls:jle)
      call putvara ('sig22e', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig22w(ils:ile,jls:jle)
      call putvara ('sig22w', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig12n(ils:ile,jls:jle)
      call putvara ('sig12n', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig12s(ils:ile,jls:jle)
      call putvara ('sig12s', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig12e(ils:ile,jls:jle)
      call putvara ('sig12e', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sig12w(ils:ile,jls:jle)
      call putvara ('sig12w', iou, ln, ib, ic, tmpij, c1, c0)
# endif
# if defined O_sulphate_data_transient
      tmpij(ils:ile,jls:jle) = sulph(ils:ile,jls:jle,2)
      call putvara ('sulph', iou, ln, ib, ic, tmpij, c1, c0)
# endif

      call rdstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      print*, '=> Atm restart written to ',trim(fname),' on ', nstamp

      deallocate ( tmpij )
#endif

      return
      end
