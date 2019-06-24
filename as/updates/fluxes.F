      subroutine fluxes (is, ie, js, je)

#if defined O_embm
!=======================================================================
!     calculate energy and moisture fluxes

!     Note: evaporation and precipitation are in g cm-2 s-1
!           and humidities are in g g-1

!     for Thompson and Warren outgoing radiation (see: Thompson S.J.,
!     and S.G. Warren 'parameterization of outgoing ...'J. Atmos. Sci.,
!     39, 2667-2680, 1982.
!=======================================================================

      implicit none

      integer i, ie, iem1, imax, is, isp1, iter
      integer j, je, jem1, jmax, js, jsp1, maxit, n
      
      logical track

      real b00, b10, b20, b01, b11, b21, b02, b12, b22, b03, b13, b23
      real delta, df, dt, dultnt, dulwr, dusens, emax, f, fb, ff, fg
      real fh, fl, fm, qair, qlnd, rhrh, sr, scrit, ssh, tair, teff
      real telev, tlnd, tlold, tol, tol2, ultnt, ulwr, usens, wspd
      real vcs, avg_sat, C2K
      
      real tmp,mski,rmski,tsl

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
      include "atm.h"
      include "csbc.h"

# if defined O_ice_cpts && defined O_ice
      include "cpts.h"
# endif
      include "ice.h"
# if defined O_embm_vcs
#  include "tmngr.h"
# endif
      include "veg.h"
      
      character(120) :: g_st
      integer ntrec
      data ntrec /0/
      save ntrec
      real terr  
     
      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

!-----------------------------------------------------------------------
!     set appropriate constants
!-----------------------------------------------------------------------
      fb = 0.94*rhoatm*cpatm
      maxit = 10
      tol = 0.01
      emax = 0.0
      imax = 0
      jmax = 0
      scrit = 0.75*soilmax
      ff = rhoatm*vlocn
      C2K = 273.15

!     Thomson and Warren constants

      b00 = 2.3829382e2
      b10 = -3.47968e1
      b20 = 1.02790e1
      b01 = 2.60065
      b11 = -1.62064
      b21 = 6.34856e-1
      b02 = 4.40272e-3
      b12 = -2.26092e-2
      b22 = 1.12265e-2
      b03 = -2.05237e-5
      b13 = -9.67e-5
      b23 = 5.62925e-5

      tol2 = tol*2.0

# if defined O_embm_vcs
      vcs = 0.0      
      track = .true.      
      avg_sat = 0.
      do n=1,ntrack_sat
        if (track_sat(n) .ge. 1.e10) track = .false.
        avg_sat = avg_sat + track_sat(n)
      enddo
      if (ntrack_sat .gt. 0) avg_sat = avg_sat/ntrack_sat
      if (track .and. relyr + year0 .ge. vcsyri) 
     &  vcs = (avg_sat - vcsref)*vcsfac
     
# endif
      
      evapl(:,:) = 0.
      evapo(:,:) = 0.
        
      do j=jsp1,jem1
        do i=isp1,iem1
          dnswro(i,j) = 0.0
          dnswrl(i,j) = 0.0
	  outlwr(i,j) = 0.0	  
	  upltntl(i,j) = 0.0
	  upltnto(i,j) = 0.0	  
	  upsensl(i,j) = 0.0
	  upsenso(i,j) = 0.0	  
	  uplwrl(i,j) = 0.0
	  uplwro(i,j) = 0.0	  
          evapl(i,j) = 0.
          evapo(i,j) = 0.
	  
!-----------------------------------------------------------------------
!         set the incoming short wave
!-----------------------------------------------------------------------
          dnswro(i,j) = 0.0
          dnswrl(i,j) = 0.0

          if (tmsk(i,j) .gt. 0.) then
# if defined O_sulphate_data_transient
            dnswro(i,j) = solins(i,j)*sbc(i,j,iaca)*pass
     &                   *(sbc(i,j,iscao) - sulph(i,j,2))
# else
            dnswro(i,j) = solins(i,j)*sbc(i,j,iaca)*pass*sbc(i,j,iscao)
# endif
          endif
          if (tmsk(i,j) .lt. 1.) then
# if defined O_sulphate_data_transient
            dnswrl(i,j) = solins(i,j)*sbc(i,j,iaca)*pass
     &                   *(sbc(i,j,iscal) - sulph(i,j,2))
# else
            dnswrl(i,j) = solins(i,j)*sbc(i,j,iaca)*pass*sbc(i,j,iscal)
# endif
          endif	  
!-----------------------------------------------------------------------
!         set wind speed and effective elevated air temperature
!-----------------------------------------------------------------------
          wspd = sbc(i,j,iws)

# if defined O_data_sat
        !where ice grid land exists, apply observed temps.  Else, use model temps.
	!Do this selectively via icemsk.
        mski =real(icemsk(i,j))
	rmski=1.-real(icemsk(i,j))
       !lapse observed temperature to sea level then bring back up with Mike's 'effective temp for humidity' lapse curve to form 
        tsl = rlapse*elev(i,j)
	telev = elev(i,j)
        teff =(sbc(i,j,i1800obs)
     &	     +tsl*(1.-rf1*exp(max(-1.,-telev/rf2))))*mski
     &       +(at(i,j,2,isat)      
     &	     -telev*rlapse*rf1*exp(max(-1.,-telev/rf2)))*rmski
# else
          telev = elev(i,j)
#  if defined O_landice_data
     &          + hicel(i,j,2)
#  endif
#  if defined O_sealev || defined O_sealev_data
     &          + elev_sealev(i,j)
#  endif
          teff = at(i,j,2,isat)
#  if defined O_use_bias_SAT_correction && defined O_use_bias_teff_correction
     &                - sbc(i,j,ibia) 
#   if !defined O_global_bias 
     &  	      * real(icemsk(i,j)) 
#   endif                 
#  endif	  
     &         - telev*rlapse*rf1*exp(max(-1.,-telev/rf2))
# endif
!-----------------------------------------------------------------------
!         calculate outgoing longwave radiation
!-----------------------------------------------------------------------
          rhrh = rh(i,j)*rh(i,j)
          outlwr(i,j) = 1.0e3*(b00 + b10*rh(i,j) + b20*rhrh
     &                + (b01 + b11*rh(i,j) + b21*rhrh)*teff
     &                + (b02 + b12*rh(i,j) + b22*rhrh)*teff**2
     &                + (b03 + b13*rh(i,j) + b23*rhrh)*teff**3)
# if defined O_carbon && defined O_carbon_co2_2d
     &                - co2for*alog(at(i,j,2,ico2)/280.0)
# else
     &                - anthro
# endif
# if defined O_aggfor_data_transient
     &                - aggfor + aggfor_os
# endif
# if defined O_embm_vcs
     &                - vcs
# endif

# if defined O_data_sat
          !where ice grid land exists, apply observed temps.  Else, use model temps.
	  !Do this selectively via icemsk.
          tair = sbc(i,j,i1800obs)*mski
     &	  +(at(i,j,2,isat)-elev(i,j)*rlapse)*rmski   
# else
          tair = at(i,j,2,isat) - elev(i,j)*rlapse
#  if defined O_landice_data
     &           - hicel(i,j,2)*rlapse
#  endif
#  if defined O_sealev || defined O_sealev_data
     &           - elev_sealev(i,j)*rlapse
#  endif
#  if defined O_use_bias_SAT_correction	  
     &                - sbc(i,j,ibia)
#   if !defined O_global_bias 
      !if bias only wanted over ice sheets (i.e. above flag) mask above bias by ice sheet model land/ice.
     &                *real(icemsk(i,j))
#   endif                      
#  endif
# endif	

          if (tmsk(i,j) .gt. 0.0) then


!-----------------------------------------------------------------------
!           calculations only for ocean points
!-----------------------------------------------------------------------
            dt = sbc(i,j,isst) - tair
            fg = dalt_o*wspd

!-----------------------------------------------------------------------
!           calculate evaporation or sublimation (ensure it is positive)
!-----------------------------------------------------------------------
            ssh = cssh*exp(17.67*sbc(i,j,isst)/(sbc(i,j,isst) + 243.5))
            evapo(i,j) = max(c0, rhoatm*fg*(ssh - at(i,j,2,ishum)))
            upltnto(i,j) = vlocn*evapo(i,j)

!-----------------------------------------------------------------------
!           calculate upward sensible heat flux
!-----------------------------------------------------------------------
            upsenso(i,j) = fb*fg*(dt)
!
!-----------------------------------------------------------------------
!           calculate upward longwave re-radiation
!-----------------------------------------------------------------------
            uplwro(i,j) = esocn*(sbc(i,j,isst) + C2K)**4
     &                 - esatm*(tair + C2K)**4
          endif
# if defined O_mtlm

          if (tmsk(i,j) .lt. 1. .and. land_map(i,j) .ne. 0) then
!----------------------------------------------------------------------
!           set fluxes over land from the land model
!---------------------------------------------------------------------
            upltntl(i,j) = 0.0
            evapl(i,j) = sbc(i,j,ievap)
            upsensl(i,j) = sbc(i,j,isens)
            uplwrl(i,j) = sbc(i,j,ilwr)
	  endif
# endif
# if defined O_mtlm
          if (tmsk(i,j) .lt. 1..and. land_map(i,j) .eq. 0 ) then
# else	  
	  if (tmsk(i,j) .lt. 1.) then 
# endif	  
	  
!-----------------------------------------------------------------------
!            calculations only for land points

!           find land temperature by balancing the surface heat budget
!             dwsr = ultnt + usens + ulwr
!           using Newton's method:
!             t(i+1) = t(i) - f(t(i))/df(t(i))
!           where:
!             f(t(i)) = dwsr - ultnt - usens - ulwr
!             -df(t(i)) = dultnt - dusens - dulwr
!-----------------------------------------------------------------------
            tlnd = surf(i,j)
            tlold = tlnd
            fm = esatm*(tair + C2K)**4
            fg = rhoatm
!           calculate stomatal resistance
# if defined O_crop_data
            sr = (1.-crops(i,j,2))*veg_rs(iveg(i,j))
     &         + crops(i,j,2)*veg_rs(icrops)
            dalt_v = veg_dalt(i,j)
# else
            sr = veg_rs(iveg(i,j))
            dalt_v = veg_dalt(iveg(i,j))
# endif
!           add in aerodynamic resistance
            sr = sr + 1.0/(dalt_v*wspd + epsln)
!           set beta parameter for calculating actual evaporation
            fh = min(max(c0+epsln, (soilm(i,j,lf)/soilmax)**(0.25)),c1)
!           set coefficients for latent heat (fl) and evaporation (fg)
            fl = fh*ff/(sr)
            fg = fh*fg/(sr)
            dusens = fb*dalt_v*wspd

!-----------------------------------------------------------------------
!           start loop for all land grid points
!-----------------------------------------------------------------------
            qair = rh(i,j)*cssh*exp(17.67*tair/(tair + 243.5))
            iter = 0
            delta = tol2
            do while (abs(delta) .gt. tol .and. iter .le. maxit)
              iter = iter + 1
              qlnd = cssh*exp(17.67*tlnd/(tlnd + 243.5))
              if (qlnd .gt. qair) then
                ultnt = fl*(qlnd - qair)
                dultnt = fl*qlnd*17.67*243.5/(tlnd + 243.5)**2
              else
                ultnt = 0.0
                dultnt = 0.0
              endif
              usens = dusens*(tlnd - tair)
              ulwr = eslnd*(tlnd + C2K)**4 - fm
              dulwr = 4.0*eslnd*(tlnd + C2K)**3
              f = dnswrl(i,j) - ultnt - usens - ulwr
              df = dultnt + dusens + dulwr
              delta = f/df
              tlnd = tlnd + delta
            enddo
            if (iter .gt. maxit) then
!             if not converged, set to last converged temperature
              if (abs(delta) .gt. emax) then
                emax = abs(delta)
                imax = i
                jmax = j
                tlnd = tlold
		terr = tlold  !Jer
              endif
            endif

!-----------------------------------------------------------------------
!           calculate fluxes on land
!-----------------------------------------------------------------------
            surf(i,j) = tlnd
            qlnd = cssh*exp(17.67*tlnd/(tlnd + 243.5))
            evapl(i,j) = max(c0, fg*(qlnd - qair))
	    tmp = soilm(i,j,lf)/dts/(1.-tmsk(i,j))
            evapl(i,j) = max(c0, min(tmp, evapl(i,j)))
# if defined O_ice
            fluxl(i,j,ishum) = evapl(i,j)!*(1.0 - aicelan(i,j,2))
     &				         !*(1.-tmsk(i,j))	    
# else
            fluxl(i,j,ishum) = evapl(i,j)!*(1.-tmsk(i,j))

# endif        
	    
	    upltntl(i,j) = vlocn*evapl(i,j)
            upsensl(i,j) = dusens*(tlnd - tair)
            uplwrl(i,j) = eslnd*(tlnd + C2K)**4 - fm

!           ensure fluxes are balanced since land can't absorb error
            upsensl(i,j) = dnswrl(i,j) - upltntl(i,j) - uplwrl(i,j)
	  endif
        enddo
      enddo
    
      if (emax .gt. 0.0) then
      
      
        write (stdout,*)
     &  '==> Warning: land surface temperature not converging: '
     &, 'emax, i, j, soilm:', emax, imax, jmax, soilm(imax,jmax,2)
      endif    
      
      ntrec=ntrec+1     
      
      return
      end

      subroutine precipitate (is, ie, js, je)

!=======================================================================
!     calculate precipitation explicitly and update humidity
!=======================================================================
      use subgrid

      implicit none

      integer i, ie, iem1, is, isp1, j, je, jem1, js, jsp1, k, n, negq

      real fb, fc, qmax, rate, tair, teff, telev, soiltemp, pson, psot
      real ssh, tmp

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "scalar.h"
      include "cembm.h"
      include "atm.h"
      include "switch.h"
      include "csbc.h"
      include "ism.h"      
# if defined O_ice_cpts && defined O_ice
      include "cpts.h"
# endif
      include "ice.h"
# if defined O_mtlm
      include "mtlm.h"
      include "levind.h"
      real hs(imt,jmt)
# endif
      real fbl,fis,fnis
      data negq /0/
      save negq
      real precipis,nmbal_snois,psnois, mbalis, hsnois,rhis    
      character(120) g_st
      real mski, rmski, tsl,temp(imt,jmt)
      integer ntrec
      data ntrec /0/
      save ntrec       
      temp(:,:)=0.
      if (eoyear) negq = 0

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

!-----------------------------------------------------------------------
!     set appropriate constants
!-----------------------------------------------------------------------
      fb = rhoatm*shq/dts
      fc = dts/rhosno
!     maximum relative humidity after rain
# if defined O_mtlm
      call unloadland (POINTS, LYING_SNOW, imt, jmt, land_map, hs)
!     convert from kg/m2 to cm
      hs(:,:) = hs(:,:)*0.1/rhosno
# endif

      psnoo(:,:) = 0.0
      psnol(:,:) = 0.0
      nmbal_sno(:,:) = 0.0  
      avgsno(:,:) = 0.
      avgmbal(:,:) = 0.
      precip(:,:) = 0.
      precipltnt(:,:) = 0.
      rh(:,:) = 0.
      k = 0
      
      do j=jsp1,jem1
        do i=isp1,iem1
	
	  precipis=0.
          nmbal_snois=0.
          psnois=0.
          mbalis=0.
          hsnois=0.
	  rhis=0.

	  !fnis: fraction of grid cell not covered by ice sheet (can be ocean or land)
	  fnis=1.-fraci(i,j)	  
	  !fbl: fraction of grid cell that is non-ice-sheet land.
	  fbl=fracbl(i,j)
!-----------------------------------------------------------------------
!         check if specific humidity is greater than rhmax of saturation
!-----------------------------------------------------------------------

# if defined O_data_sat
        !where ice grid land exists, apply observed temps.  Else, use model temps.
	!Do this selectively via icemsk.
          mski =real(icemsk(i,j))
	  rmski=1.-real(icemsk(i,j))
       !lapse observed temperature to sea level then bring back up with Mike's 'effective temp for humidity' lapse curve to form teff
          tsl = rlapse*elev(i,j)
	  telev = elev(i,j)
          teff =(sbc(i,j,i1800obs)
     &	     +tsl*(1.-rf1*exp(max(-1.,-telev/rf2))))*mski
     &       +(at(i,j,2,isat)      
     &	     -telev*rlapse*rf1*exp(max(-1.,-telev/rf2)))*rmski
# else
          telev = elev(i,j)
#  if defined O_landice_data
     &          + hicel(i,j,2)
#  endif
#  if defined O_sealev || defined O_sealev_data
     &          + elev_sealev(i,j)
#  endif
          teff = at(i,j,2,isat)
#  if defined O_use_bias_SAT_correction && defined O_use_bias_teff_correction	  
     &                - sbc(i,j,ibia)
#   if !defined O_global_bias 
     &  	      * real(icemsk(i,j)) 
#   endif          
#  endif	  
     &         - telev*rlapse*rf1*exp(max(-1.,-telev/rf2))
# endif

          ssh = cssh*exp(17.67*teff/(teff + 243.5))
          qmax = rhmax*ssh
          if (at(i,j,2,ishum) .gt. qmax) then
            tmp = fb*(at(i,j,2,ishum) - qmax)
            precip(i,j) = precip(i,j) + tmp
            at(i,j,2,ishum) = at(i,j,2,ishum) - tmp/fb
	    
     &	    *fnis
	    
            rh(i,j) = rhmax
          else
            rh(i,j) = max(rh(i,j), at(i,j,2,ishum)/(ssh + epsln))
            rh(i,j) = max(c0, min(c1, rh(i,j)))
          endif

!-----------------------------------------------------------------------
!         calculate snowfall (hsno at tau was set in the ice model)
!-----------------------------------------------------------------------
# if defined O_data_sat
          !where ice grid land exists, apply observed temps.  Else, use model temps.
	  !Do this selectively via icemsk.
          tair = sbc(i,j,i1800obs)*mski
     &	  +(at(i,j,2,isat)-elev(i,j)*rlapse)*rmski   
# else
!         tair may be adjusted by a snowfall offset temperature tsno          
	  tair = at(i,j,2,isat) - elev(i,j)*rlapse - tsno
#  if defined O_landice_data
     &           - hicel(i,j,2)*rlapse
#  endif
#  if defined O_sealev || defined O_sealev_data
     &           - elev_sealev(i,j)*rlapse
#  endif
#  if defined O_use_bias_SAT_correction	  
     &                - sbc(i,j,ibia)
#   if !defined O_global_bias 
      !if bias only wanted over ice sheets (i.e. above flag) mask above bias by ice sheet model land/ice.
     &                *real(icemsk(i,j))
#   endif                      
#  endif
# endif	

# if defined O_ice_cpts && defined O_ice

!code chunk removed for clarity
      stop 'O_ice_cpts not compatible with ice sheet or subgrid'.

# elif defined O_ice
          !calculate snowfall over ocean
          if (tmsk(i,j) .gt. 0.) then
            if (tair .le. c0 .and. hsnoo(i,j,2) .lt. hsno_max)
     &        psnoo(i,j)=min((hsno_max-hsnoo(i,j,2))/fc,precip(i,j))
     &        *tmsk(i,j)
!           only allow snow where there is sea ice
            psnoo(i,j) = psnoo(i,j)*aiceocn(i,j,2)
            hsnoo(i,j,2) = hsnoo(i,j,2) + fc*psnoo(i,j)
            if (addflxa) fluxo(i,j,ishum) = fluxo(i,j,ishum)
     &                         -         dts*psnoo(i,j)
          endif

          !calculate snowfall over land	
          if (tmsk(i,j) .lt. 1.) then
#  if defined O_mtlm
            hs(i,j) = hs(i,j) + hsnol(i,j,2)
            if (tair .le. c0 .and. hs(i,j) .lt. hsno_max)
     &        psnol(i,j) = min((hsno_max - hs(i,j))/fc, precip(i,j))
#  else
             if (tair .le. c0 .and. hsnol(i,j,2) .lt. hsno_max)
     &        psnol(i,j)=min((hsno_max-hsnol(i,j,2))/fc,precip(i,j))
#  endif
	  
	    if (icemsk(i,j) .gt. 0) then

              if (tair .le. c0) then
#if !defined O_constant_ice
	        nmbal_sno(i,j) =(precip(i,j)-psnol(i,j))
     &		*min(cir(i,j),fbl)
#endif
	        
		nmbal(i,j,2) = nmbal(i,j,2) + fc*nmbal_sno(i,j)		
		hsnol(i,j,2) = hsnol(i,j,2) + fc*psnol(i,j)*fbl
		psnol(i,j) = psnol(i,j)*fbl
	      endif
	      !NOTE::::  This code commented out in this version: when one ice sheet switched off, but model restarted from a 2-sheet model, shelf regions are not captured by mtlm.  But above code does not provide psnol for non-mtlm, non-icemsk, points.  So, as a fix, commented out read-in of tmsk in embm_rest, which resets all tmsk to kmt-derived values for land and ocean, which gets rid of any non-mtlm that is not associated with cie sheet grids.
            !non-icemsk, but still non-mtlm, grid
 	    elseif (land_map(i,j).eq.0 .and. icemsk(i,j).eq.0) then
              if (tair .le. c0) then
 	        hsnol(i,j,2) = hsnol(i,j,2) + fc*psnol(i,j)*fbl	
		psnol(i,j) = psnol(i,j)*fbl	
 		temp(i,j)=1.
              endif
	    endif
#  if !defined O_mtlm
           if (icemsk(i,j) .eq. 0) then	      
	     hsnol(i,j,2) = hsnol(i,j,2) + fc*psnol(i,j)
           endif
#  endif
 	   
	  endif	  

          !calculate snow accumulation over ice sheets using subgrid scheme
	  call fluxis(i,j,precipis,nmbal_snois,psnois,rhis)	  
	  !get total precipitation field for passing to latent heat calculation in solve.  Scale the full
	  !precip flux by the non-ice-sheet area and add on the ice-sheet-derived precipis flux.
	  precipltnt(i,j)=precip(i,j)*fnis+precipis
	  !add accumulation to psnol and nmbal fields to be added to fluxl correction below, and for 
	  !latent heat correction in solve.	  
	  psnol(i,j)=psnol(i,j)+psnois
	  nmbal_sno(i,j)=nmbal_sno(i,j)+nmbal_snois
	  !remove ice sheet precip from atmospheric moisture bin
          at(i,j,2,ishum) = at(i,j,2,ishum) + rhis
# endif

!-----------------------------------------------------------------------
!         update soilm and allocate surplus soil moisture to runoff
!-----------------------------------------------------------------------
# if defined O_mtlm
	  if (tmsk(i,j) .lt. 1. .and. land_map(i,j) .eq. 0) then
# else
          if (tmsk(i,j) .lt. 1.) then
# endif
	    fluxl(i,j,ishum) = fluxl(i,j,ishum) 
     &	    	             - precip(i,j)*fbl
     &	    	             - precipis    
     &			     + psnol(i,j)    
     &			     + nmbal_sno(i,j)   
            soiltemp = soilm(i,j,2)
            soilm(i,j,2) = soilm(i,j,lf) - dts*fluxl(i,j,ishum)
            soilm(i,j,2) = max(c0, soilm(i,j,2))
            soilm(i,j,1) = soiltemp
            if (soilm(i,j,2) .gt. soilmax) then
              runoff(i,j) = (soilm(i,j,2) - soilmax)/dts
              soilm(i,j,2) = soilmax
            else
              runoff(i,j) = c0
            endif
          endif	    
        enddo
      enddo
       	
      call embmbc (psnol)
      call embmbc (psnoo)
      call embmbc (nmbal_sno)   
# if defined O_ice_cpts && defined O_ice
      call embmbc (hseff(1,1,idx,1))
# else
      call embmbc (hsnoo(1,1,2))
      call embmbc (hsnol(1,1,2))
      call embmbc (nmbal(1,1,2))      
# endif
       
      return
      end

      subroutine co2forc
!=======================================================================
!     calculate global average CO2 forcing
!=======================================================================
# if !defined O_carbon_co2_2d

      implicit none

      real yr

      include "cembm.h"

!-----------------------------------------------------------------------
!     relative forcing from 280 ppmv (anthro)
!-----------------------------------------------------------------------
      anthro = co2for*alog(co2ccn/280.0)
# endif
#endif

      return
      end
