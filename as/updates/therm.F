      subroutine therm  (is, ie, js, je)

#if defined O_ice && defined O_embm
!=======================================================================
!     thermodynamic ice model

!     Note: if run with embm this routine must be called after "fluxes"
!           and before "solve"

!     calculates ice and open water growth rates based on the surface
!     energy budget. see Parkinson and Washington, JGR, Vol.84, C1,
!     311-337, 1979 and Hibler, JPO, Vol.9, 815-846, 1979

!     heat and fresh water fluxes between the ocean and atmosphere
!     are adjusted depending on ice growth or melt. ice thickness is
!     changed by the amount of growth or melt and sublimation
!=======================================================================

      implicit none

      integer i, ie, iem1, imax, index, iter, is, isp1, j, je, jem1
      integer jmax, js, jsp1, maxit

      real ai, aicelan3, aiceocn3, al, amin, ao, as, as_crops, ca  
      real delta, df, dh
      real dha, dhflxi, dhflxs, dhi, hice3, dho, dhs
      real hsnoo3, hsnol3, dhss, nmbal3 
      real hsno_aft_sub
      real dhstot, dhtot, dswr, dt, dultnt, dulwr, dusens, emax, f, fa
      real fas, fb, fbot, fcond, fd, fds, fe, ff, ffs, fh, fl, fls, fm
      real fn, fptf, fpts, ftopi, ftopo, ho, hsextra, qair,qice, sla
      real sub, tair, tcdh, ti, tiold, tol, tol2, ultnt, ulwr, usens
      real wspd, zintfc, C2K
      real snwthk,icear
      real fbl,fsc      
      real evapis,dsmis,dswris,ultntis,usensis,ulwris
      real sca

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "cembm.h"
      include "atm.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"
      include "coord.h"
      include "grdvar.h"
      include "veg.h"
# if defined O_mtlm
      include "mtlm.h"
# endif
      include "ism.h"
      include "temp.h"
      
      real under_calb_lan,under_calb_ocn,under_calb
      real snow_coalb, wt
      character(120) g_st
      real mski,rmski, temp(imt,jmt)
      integer ntrec
      data ntrec /0/
      save ntrec
      
      non_grid_surf_alb(:,:)=0.
      temp(:,:)=0.
      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

      fa = dts/(rhoice*flice)
      fb = 0.94*rhoatm*cpatm
      fd = rhoatm/rhoice
      fe = rhoatm*slice
      ff = rhoice*flice
      fh = 21.8746*265.5
      fas = dts/(rhosno*flice)
      fds = rhoatm/rhosno
      ffs = rhosno*flice
      ho = 1.0
      amin = 0.15
      maxit = 10
      tol = 0.01
      emax = 0.0
      imax = 0
      jmax = 0
      fptf = 0.0
# if defined O_ice_evp
      index = 1
# else
      index = lf
# endif
      sla = zw(1)*secday/dampice/2.389e-8
      tol2 = tol*2.0
      C2K = 273.15 
       
      do j=jsp1,jem1
        do i=isp1,iem1
          hsnol3 = hsnol(i,j,index)
          hsnoo3 = hsnoo(i,j,index)	  
          hice3 = hice(i,j,index)
          aicelan3 = 0.
	  aiceocn3 = 0.
	  nmbal3 = nmbal(i,j,index)	  
	  hsno_aft_sub = 0.
	  evapis = 0.
	  dsmis = 0.
          dswris = 0.
	  ultntis = 0.
	  usensis = 0.
	  ulwris = 0.
	  
!-----------------------------------------------------------------------
!         set the incoming shortwave over snow and ice
!-----------------------------------------------------------------------
!         if snow is less than 25 cm linearly reduce to underlying 
!         albedo.  Over land, if land-ice exists underneath, use ice
!         as underlying albedo.  If no ice exists underneath, use prescribed
!         iveg albedo.
	  
	  under_calb_lan=0.
	  under_calb_ocn=0.
	  if (tmsk(i,j) .lt. 1.) then
	  !where land exists, blend ice albedo with underlying bare
	  !land albedo
	  
	  under_calb_lan=1.-veg_alb(iveg(i,j))
	  endif
	  if (tmsk(i,j) .gt. 0.) then
	  !where ocean exists, use ice albedo (for sea ice calcs)
            under_calb_ocn=sice_calb
	  endif

          !blend underlying albedos over partial cells
	  under_calb=under_calb_lan*(1.-tmsk(i,j))
     &              +under_calb_ocn*tmsk(i,j)	 

	  if (icemsk(i,j) .gt. 0.) then
	    !set weight to vary between 0 and 1 if between ltlim and utlim
	    wt=min(1.,max(0.,(ticel(i,j)-ltlim)/(utlim-ltlim)))
	    !set snow_coalb to drift linearly from cold to warm snow albedos between ltlim and utlim temperatures
	    snow_coalb=csno_calb*(1.-wt)+wsno_calb*wt
	  else
	    snow_coalb=csno_calb
	  endif
	  !merge underlying albedos with snow albedo if less than
	  !25 cm of snow.	  
          snwthk = max(hsnol(i,j,2),hsnoo(i,j,2))	  
	  as = max(0.,min(snwthk/25., c1))
	  
          ca = under_calb*(c1 - as) + snow_coalb*as
	  non_grid_surf_alb(i,j)=1.-ca
# if defined O_sulphate_data_transient
          dswr = solins(i,j)*sbc(i,j,iaca)*pass*(ca - sulph(i,j,2))
# else
          dswr = solins(i,j)*sbc(i,j,iaca)*pass*ca
# endif
          wspd = sbc(i,j,iws)

# if defined O_ice_evp
          aiceocn(i,j,index) = min(c1, aiceocn(i,j,index))
# endif

# if  defined O_mtlm
          if (tmsk(i,j) .lt. 1. .and. land_map(i,j) .eq. 0) then
# else
          if (tmsk(i,j) .lt. 1.) then
# endif

!-----------------------------------------------------------------------
!           partial/full land points
!-----------------------------------------------------------------------

            ai = aicelan(i,j,2)	    
# if defined O_data_sat
            !where ice grid land exists, apply observed temps.  Else, use model temps.
	    !Do this selectively via icemsk.
            mski =real(icemsk(i,j))
	    rmski=1.-real(icemsk(i,j))
            tair = sbc(i,j,i1800obs)*mski
     &	    +(at(i,j,2,isat)-elev(i,j)*rlapse)*rmski
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
            aicelan3 = 0.0
	    temp(i,j)=under_calb
            as = hsnol(i,j,2)/(100.0*veg_smd(iveg(i,j)))

# if defined O_crop_data
            as_crops = hsnol(i,j,2)/(100.0*veg_smd(icrops))
            as_crops = min(c1, max(c0, as_crops))
            as = (1.-crops(i,j,2))*as + crops(i,j,2)*as_crops
# endif
!           limit snow coverage between 0 and 1
            aicelan3 = min(c1, max(c0, as))
# if defined O_landice_data
            aicelan3 = max(aicelan3, aicel(i,j,2))
# endif

            al = 1.0 - ai
	    !fbl: fraction of grid cell covered by bare non-ice-or-snow-covered land
	    fbl=al*fracbl(i,j)
	    !fsc: fraction of grid cell covered by non-ice-sheet snow cover
	    fsc=ai*fracbl(i,j)
!-----------------------------------------------------------------------
!           start loop for grid points with snow or ice
!-----------------------------------------------------------------------
            if (ai .gt. c0) then

!-----------------------------------------------------------------------
!             find snow temperature by balancing the surface heat budget
!               dwsr = ultnt + usens + ulwr
!             using Newton's method:
!               t(i+1) = t(i) - f(t(i))/df(t(i))
!             where:
!               f(t(i)) = dwsr - ultnt - usens - ulwr
!               -df(t(i)) = dultnt - dusens - dulwr
!-----------------------------------------------------------------------

              ti = ticel(i,j)
              tiold = ticel(i,j)
              fm = esatm*(tair + C2K)**4
              fls = fe*dalt_i*wspd
              dusens = fb*dalt_i*wspd
              qair = rh(i,j)*cssh*exp(17.67*tair/(tair + 243.5))
              iter = 0
              delta = tol2
	      
              do while (abs(delta) .gt. tol .and. iter .le. maxit)
                iter = iter + 1
                dt = ti - tair
                qice = cssh*exp(21.8746*ti/(ti + 265.5))
                if (qice .gt. qair) then
                  ultnt = fls*(qice - qair)
                  dultnt = fls*qice*21.8746*265.5/(ti + 265.5)**2
                else
                  ultnt = 0.0
                  dultnt = 0.0
                endif
                usens = dusens*dt
                ulwr = esice*(ti + C2K)**4 - fm
                dulwr = 4.0*esice*(ti + C2K)**3
                f = dswr - ultnt - usens - ulwr
                df = dultnt + dusens + dulwr
                delta = f/df
                ti = ti + delta
              enddo
              if (iter .gt. maxit) then
!               if not converged, set to last converged temperature
                if (abs(delta) .gt. emax .and. ti .lt. fptf) then
                  emax = abs(delta)
                  imax = i
                  jmax = j
                  ti = tiold
                endif
              endif

!-----------------------------------------------------------------------
!             set maximum tice to freezing and calculate fluxes
!-----------------------------------------------------------------------
              ti = min(ti, fptf)
              dt = ti - tair
              qice = cssh*exp(21.8746*ti/(ti + 265.5))
              sub = max(c0, fds*dalt_i*wspd*(qice - qair))
              usens = dusens*dt
              ulwr = esice*(ti + C2K)**4 - fm

!             ensure that snow sublimated does not exceed effective bareland hsno thickness
              dha = -dts*sub*fsc

	      if (-dha .gt. hsnol3) then
	        dha = -hsnol3
		!scale sub back up to full flux
		sub = -dha/dts/fsc
	      endif

              ultnt = rhosno*slice*sub
              sub = sub*rhosno

              ftopi = dswr - ulwr - usens - ultnt
     
!-----------------------------------------------------------------------
!             calculate total change in snow volume on land
!             allocate snowmelt to flux for input to bucket
!-----------------------------------------------------------------------
              dhs = 0.0
              if (ticel(i,j) .ge. fptf .and. ftopi .gt. 0.0)
     &          dhs = -fas*ftopi*fsc	      
              !determine snow left after sublimation
	      hsno_aft_sub = hsnol3 + dha	      
              
	      if (-dhs .gt. hsno_aft_sub) then
                dhs = -hsno_aft_sub	      
	      endif
	      
              !snow reduced by scaled melt plus sublimation
	      hsnol3 = hsnol3 + dhs + dha
	      
              !land flux increased by scaled melt
	      !REDUCE evaporative flux from land bucket by bareland fraction
              fluxl(i,j,ishum) = fluxl(i,j,ishum)*fbl
	      !THEN add on fluxes from snow melt.
              fluxl(i,j,ishum) = fluxl(i,j,ishum)	      
     & 	                       + dhs*rhosno/dts	      
              ticel(i,j) = ti
              evapl(i,j) = evapl(i,j)*fbl + sub*fsc
     
!             add ice/snow covered area fluxes to land fluxes
              dnswrl(i,j) = dnswrl(i,j)*fbl + dswr*fsc
              upltntl(i,j) = upltntl(i,j)*fbl + ultnt*fsc
              upsensl(i,j) = upsensl(i,j)*fbl + usens*fsc
              uplwrl(i,j) = uplwrl(i,j)*fbl + ulwr*fsc 

!             ensure fluxes are balanced since land can't absorb error
              upsensl(i,j) = dnswrl(i,j) - upltntl(i,j) - uplwrl(i,j)
     &                   + dhs*ffs/dts
	      	      
!-----------------------------------------------------------------------
!     if no snow on land
!-----------------------------------------------------------------------

            else
              ticel(i,j) = 0.0              
	      fluxl(i,j,ishum) = fluxl(i,j,ishum)*fbl
              evapl(i,j) = evapl(i,j)*fbl	      	      
              dnswrl(i,j) = dnswrl(i,j)*fbl
              upltntl(i,j) = upltntl(i,j)*fbl 
              upsensl(i,j) = upsensl(i,j)*fbl 
              uplwrl(i,j) = uplwrl(i,j)*fbl  
            endif

  	    call thermis(i,j,evapis,dsmis
     &                   ,dswris,ultntis,usensis,ulwris)

!    	    !Add mass accumulations and fluxes from ice sheet.	    
  	    fluxl(i,j,ishum) = fluxl(i,j,ishum)
     &      + dsmis*rhosno/dts         
             evapl(i,j) = evapl(i,j)
     &      + evapis   
     
            !Add energy fluxes from ice sheet.  
	    dnswrl(i,j) = dnswrl(i,j) + dswris 
            upltntl(i,j) = upltntl(i,j) + ultntis 
            upsensl(i,j) = upsensl(i,j) + usensis 
            uplwrl(i,j) = uplwrl(i,j) + ulwris

          endif
# if !defined O_ice_cpts	  
          if (tmsk(i,j) .gt. 0.0) then
	  
!-----------------------------------------------------------------------
!           ocean points
!-----------------------------------------------------------------------

            ai = aiceocn(i,j,2)

            ao = 1.0 - ai
            tair = at(i,j,2,isat)
#  if defined O_sealev_data || defined O_sealev_data
     &           - elev_sealev(i,j)*rlapse
#  endif
# if defined O_use_bias_SAT_correction && defined O_global_bias	
  
     &                - sbc(i,j,ibia)                 
# endif	
!-----------------------------------------------------------------------
!           calculate fluxes to and from the ocean (without ice)
!-----------------------------------------------------------------------
            ftopo = (dnswro(i,j)-uplwro(i,j)-upsenso(i,j)-upltnto(i,j))
     &             ! *tmsk(i,j)
            fbot = (sla*(frzpt(i,j) - sbc(i,j,isst)))
     &             ! *tmsk(i,j) 
!           calculate growth of ice in open water areas
            dho = fa*(fbot - ftopo)

            if (ai .ne. 0.0) then

!-----------------------------------------------------------------------
!             find ice temperature by balancing the surface heat budget:
!               tcdh*(ti - fpts) = dswr - ultnt - usens - ulwr
!             using Newton's method:
!               t(i+1) = t(i) - f(t(i))/df(t(i))
!             where:
!               f(t(i)) = dswr - ultnt - usens - ulwr - tcdh*(ti - fpts)
!               -df(t(i)) = dultnt + dusens + dulwr + tcdh
!-----------------------------------------------------------------------
              tcdh = condice/(hice(i,j,2) + 6.5*hsnoo(i,j,2))
              ti = ticeo(i,j)
              tiold = ticeo(i,j)
              fpts = frzpt(i,j)
              fm = esatm*(tair + C2K)**4
              fn = 4.0*esice

              fl = fe*dalt_i*wspd
              dusens = fb*dalt_i*wspd
              qair = at(i,j,2,ishum)
              iter = 0
              delta = tol2
              do while (abs(delta) .gt. tol .and. iter .le. maxit)
                iter = iter + 1
                dt = ti - tair
                qice = cssh*exp(21.8746*ti/(ti + 265.5))
                if (qice .gt. qair) then
                  ultnt = fl*(qice - qair)
                  dultnt = fl*qice*fh/(ti + 265.5)**2
                else
                  ultnt = 0.0
                  dultnt = 0.0
                endif
                usens = dusens*dt
                ulwr = esice*(ti + C2K)**4 - fm
                dulwr = fn*(ti + C2K)**3
                f = dswr - ultnt - usens - ulwr - tcdh*(ti - fpts)
                df = dultnt + dusens + dulwr + tcdh
                delta = f/df
                ti = ti + delta
              enddo
              if (iter .gt. maxit) then
!               if not converged, set to last converged temperature
                if (abs(delta) .gt. emax .and. ti .lt. fptf) then
                  emax = abs(delta)
                  imax = i
                  jmax = j
                  ti = tiold
                endif
              endif

!-----------------------------------------------------------------------
!             set maximum tice to freezing and calculate fluxes
!-----------------------------------------------------------------------
              ti = min(ti, fptf)
              dt = ti - tair
              qice = cssh*exp(21.8746*ti/(ti + 265.5))
              sub = max(c0, dalt_i*wspd*(qice - qair))
	      ultnt = fe*sub	      
              fcond = tcdh*(ti - fpts)!*tmsk(i,j)
              if (hsnoo(i,j,index) .gt. 0.0) then
                sub = fds*sub
                dha = -dts*sub
                sub = sub*ai*rhosno
              else
                sub = fd*sub
                dha = -dts*sub
                sub = sub*ai*rhoice
              endif
!	      dha = dha*tmsk(i,j)
!	      sub = sub*tmsk(i,j)     

	      usens = dusens*dt                
              ulwr  = esice*(ti + C2K)**4 - fm 

!-----------------------------------------------------------------------
!             add ice covered area fluxes to ocean area fluxes
!-----------------------------------------------------------------------
              ticeo(i,j) = ti
              dnswro(i,j) = dnswro(i,j)*ao + dswr*ai
              upltnto(i,j) = upltnto(i,j)*ao + ultnt*ai
              upsenso(i,j) = upsenso(i,j)*ao + usens*ai
              uplwro(i,j) = uplwro(i,j)*ao + ulwr*ai
              ftopi = (dswr - ulwr - usens - ultnt)!*tmsk(i,j)

!-----------------------------------------------------------------------
!             calculate change in ice volume due to sublimation of
!             ice (dha). adjust evaporation to the atmosphere to
!             account for sublimation from ice. subtract this
!             adjustment from the ocean freshwater flux
!-----------------------------------------------------------------------
	      if (addflxa) fluxo(i,j,ishum) = fluxo(i,j,ishum) 
     &	      + dts*sub*tmsk(i,j)
              evapo(i,j) = evapo(i,j)*ao + sub 
            else
              ticeo(i,j) = sbc(i,j,isst)
              ftopi = 0.0
              fcond = 0.0
              dha = 0.0
            endif

!-----------------------------------------------------------------------
!           calculate total change in ice volume (dh)
!-----------------------------------------------------------------------
            dha = dha*ai*tmsk(i,j)!scale sub. amount by tmsk
            dhflxs = 0.0
            dhs = 0.0
            if (hsnoo(i,j,index) .le. 0.0) then
!             total growth of ice from the ocean
              dhi = ai*fa*(fbot - ftopi) + ao*dho
              dhi = dhi*tmsk(i,j)
!             total growth (loss + sublimation limited to total amount)
              dh = max(-hice(i,j,index), dhi + dha)
!             adjust ocean fluxes for ice growth or melt + sublimation
              dhflxi = dh - dha
            else
!             total growth of ice from the ocean
              !flu balance from conductino through snow and from ocean
              dhi = ai*fa*(fbot - fcond) + ao*dho
              dhi = dhi*tmsk(i,j)
!             loss of snow due to melt
              if (ticeo(i,j) .ge. fptf) dhs = ai*fas*(fcond - ftopi)
              dhs = dhs*tmsk(i,j)
!             total loss of snow including sublimation
              dhs = dhs + dha
!             check if snow loss greater than total
              if (-dhs .gt. hsnoo(i,j,index)) then
!               take extra melt from ice
                dhi = dhi + rhosno/rhoice*(dhs + hsnoo(i,j,index))
!               remove all snow
                dhs = -hsnoo(i,j,index)
              endif
!             adjust ocean fluxes for snow melt and sublimation
              dhflxs = dhs - dha
              
!             ice loss limited to total
              dh = max(-hice(i,j,index), dhi)
!             adjust ocean fluxes for ice growth or melt
              dhflxi = dh
            endif

!-----------------------------------------------------------------------
!           calculate new area and thickness from thermodynamics
!-----------------------------------------------------------------------
!           use minimum area (amin) of open water
            ai = max(amin, aiceocn(i,j,index))
            aiceocn3 = aiceocn(i,j,index) + ((1.0 - ai)*max(c0, dho)/ho
     &         + 0.5*min(c0, dhi)*ai/(hice(i,j,index) + epsln))
!           update ice and snow thickness
            hice3 = hice(i,j,index) + dh
            hsnoo3 = hsnoo(i,j,index) + dhs

!           lower ice area where thickness is < 1 cm
            aiceocn3 = min(aiceocn3, hice3)
!           max ice area where thickness is > 10 m
            aiceocn3 = max(aiceocn3, hice3*0.001)
            aiceocn3 = max(c0, min(c1, aiceocn3))
            if (aiceocn3 .eq. 0.0) then
              dhflxs = dhflxs - hsnoo3
              hsnoo3 = 0.0
            endif
	   
!           check if the weight of the snow pushes the ice/snow
!           interface below the waterline (if so, change snow to ice)
            zintfc = hice3 - (rhosno*hsnoo3 + rhoice*hice3)/rhoocn
            if (zintfc .lt. 0.0) then
              dhss = rhoice/rhosno*zintfc
              if (-dhss .gt. hsnoo3) then
                write(*,*) '==> Warning: dhss is too large: ',dhss
                dhss = -hsnoo3
              endif
              hice3 = hice3 - rhosno/rhoice*dhss
              hsnoo3 = hsnoo3 + dhss
            endif
            hsnoo3 = max(hsnoo3, c0)
	    
!-----------------------------------------------------------------------
!           adjust fluxes to the ocean due to ice melt or growth
!-----------------------------------------------------------------------
            if (addflxa) then
             !fluxo is an EFFECTIVE accumulator, so these additions of
             !thickness should be after they have been scaled down by
             !tmsk
              fluxo(i,j,isat) = fluxo(i,j,isat) + 
     &                         (ff*dhflxi + ffs*dhflxs)
              fluxo(i,j,ishum) = fluxo(i,j,ishum) - 
     &                         (rhoice*dhflxi + rhosno*dhflxs)

#  if defined O_convect_brine
              ao = 1.0 - aiceocn(i,j,2)
              dho = ao*fa*(fbot - ftopo)
              if (dho .gt. 0.) then
                fluxo(i,j,ishum) = fluxo(i,j,ishum) + rhoice*dho
                cbf(i,j,0) = cbf(i,j,0) - rhoice*dho
                cba(i,j,0) = cba(i,j,0) + ao*dts
              endif
              ai = aiceocn(i,j,2)
              dhi = ai*fa*(fbot - ftopi)
              if (dhi .gt. 0.) then
                fluxo(i,j,ishum) = fluxo(i,j,ishum) + rhoice*dhi
                cbf(i,j,1) = cbf(i,j,1) - rhoice*dhi
                cba(i,j,1) = cba(i,j,1) + ai*dts
              endif
#   endif
            endif
#  if defined O_plume_brine 
            subflux(i,j,ishum) = subflux(i,j,ishum) - rhoice*dhflxi
#  endif
# endif

          endif

!-----------------------------------------------------------------------
!         shuffle time levels
!-----------------------------------------------------------------------
          hice(i,j,1) = hice(i,j,2)
          hice(i,j,2) = hice3
	  hsnol(i,j,1) = hsnol(i,j,2)
          hsnol(i,j,2) = hsnol3
	  nmbal(i,j,1) = nmbal(i,j,2)
          nmbal(i,j,2) = nmbal3 	
	  aiceocn(i,j,1) = aiceocn(i,j,2)
          aiceocn(i,j,2) = aiceocn3          	  
	  aicelan(i,j,1) = aicelan(i,j,2)
          aicelan(i,j,2) = aicelan3
	          
	  hsnoo(i,j,1) = hsnoo(i,j,2)
          hsnoo(i,j,2) = hsnoo3
        enddo
      enddo                    
      
      if (emax .gt. 0.0) write (stdout,*)
     &  '==> Warning: ice temperature not converging: emax, i, j:'
     &, emax, imax, jmax
     
!-----------------------------------------------------------------------
!     set boundary conditions
!-----------------------------------------------------------------------
      call embmbc (hice(1,1,2))
      call embmbc (aicelan(1,1,2))
      call embmbc (aiceocn(1,1,2))      
      call embmbc (hsnoo(1,1,2))
      call embmbc (hsnol(1,1,2))
      call embmbc (nmbal(1,1,2))   
      call embmbc (ticel)     

      ntrec=ntrec+1 
#endif     
      return
      end
