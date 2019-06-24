! source file: /net/mare/home1/eby/as/updates/thermis.F
      subroutine thermis(i,j,evapis,dsmis
     &                   ,dswris,ultntis,usensis,ulwris)

      use subgrid

      implicit none

      integer i, imax, index, iter, j, jmax, maxit, n

      real as, asl, ca,  delta, df, dha, dhs
      real dhtot, dswr, dt, dultnt, dulwr, dusens, emax
      real f, fa, fas, fb, fds, fe, ffs
      real fls, fm, fptf, ftopi, qair
      real qice, sub, tair, tinew, tiold, tol, tol2, ultnt
      real ulwr, usens, wspd, C2K, mbal3, nmbal_sub
      real hsno_asub, nmbal_melt, hsno3
      real evapis,dsmis
      real dswris,ultntis,usensis,ulwris
      real mbalmelt3, sca
      real snow_coalb,lapse_rate,wt
      real hr,pr,shci,lhci,rhoc,rhopc,refreeze

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "cembm.h"
      include "atm.h"
      include "ice.h"
      include "coord.h"
      include "grdvar.h"
      include "ism.h"
      include "veg.h"
      include "calendar.h"

      fa = dts/(rhoice*flice)
      fb = 0.94*rhoatm*cpatm
      fe = rhoatm*slice
      fas = dts/(rhosno*flice)
      fds = rhoatm/rhosno
      ffs = rhosno*flice
      maxit = 10
      tol = 0.01
      emax = 0.0
      imax = 0
      jmax = 0
      fptf = 0.0
      shci=2.054 !Jg^1K^1
      lhci=334.  !Jg^1
      rhoc=rhosno !g/cm^3
      rhopc=0.830 !g/cm^3

!# if defined 1
!      index = 1
!# else
      index = lf
!# endif
      tol2 = tol*2.0
      C2K = 273.15

      dsmis = 0.
      evapis = 0.
      dswris = 0.
      ultntis = 0.
      usensis = 0.
      ulwris = 0.

      lapse_rate=rlapse
      do n=strtiE(i,j),endiE(i,j)
        hsno3 = hsnowiE(n,index)
	mbalmelt3 = mbalmelt(n,index)
        wspd = sbc(i,j,iws)
        nmbal_melt = 0.
        nmbal_sub  = 0.
        hsno_asub = 0.
	dhs=0.
	dha=0.
	sub=0.
!-----------------------------------------------------------------------
!           land ice points
!-----------------------------------------------------------------------
!           calculate area of ice and snow
!           if snow is less than 10 cm, linearly reduce area (drifting)
	as = min(hsnowiE(n,2)*0.1, c1)
!           also limit snow area by surface masking depth.
        asl = hsnowiE(n,2)/(100.0*veg_smd(7)*as + epsln)
        asl = min(c1, max(c0, asl))
        tair = at(i,j,2,isat) - eleviE(n)*lapse_rate
     &                - sbciE(n)
	satiE(n)=tair
!-----------------------------------------------------------------------
!           set the incoming shortwave over snow and ice
!-----------------------------------------------------------------------
        as = as*asl
	!set weight to vary between 0 and 1 if between ltlim and utlim
	wt=min(1.,max(0.,(ticeiE(n)-ltlim)/(utlim-ltlim)))
	!set snow_coalb to drift linearly from cold to warm snow albedos between ltlim and utlim temperatures
	snow_coalb=csno_calb*(1.-wt)+wsno_calb*wt
	ca = (gice_calb)*(1.-as) + (snow_coalb)*as
        ca = min(c1, max(c0, ca))
	albedoiE(n)=1.-ca
        dswr = solins(i,j)*sbc(i,j,iaca)*pass*ca
!-----------------------------------------------------------------------
!             find snow temperature by balancing the surface heat budget
!               dwsr = ultnt + usens + ulwr
!             using Newton's method:
!               t(i+1) = t(i) - f(t(i))/df(t(i))
!             where:
!               f(t(i)) = dwsr - ultnt - usens - ulwr
!               -df(t(i)) = dultnt - dusens - dulwr
!-----------------------------------------------------------------------
        tinew = ticeiE(n)
        tiold = ticeiE(n)
        fm = esatm*(tair + C2K)**4
        fls = fe*dalt_i*wspd
        dusens = fb*dalt_i*wspd
        qair = rh(i,j)*cssh*exp(17.67*tair/(tair + 243.5))
        iter = 0
        delta = tol2
        do while (abs(delta) .gt. tol .and. iter .le. maxit)
          iter = iter + 1
          dt = tinew - tair
          qice = cssh*exp(21.8746*tinew/(tinew + 265.5))
          if (qice .gt. qair) then
            ultnt = fls*(qice - qair)
            dultnt = fls*qice*21.8746*265.5/(tinew + 265.5)**2
          else
            ultnt = 0.0
            dultnt = 0.0
          endif
          usens = dusens*dt
          ulwr = esice*(tinew + C2K)**4 - fm
          dulwr = 4.0*esice*(tinew + C2K)**3
          f = dswr - ultnt - usens - ulwr
          df = dultnt + dusens + dulwr
          delta = f/df
          tinew = tinew + delta
        enddo
        if (iter .gt. maxit) then
!               if not converged, set to last converged temperature
          if (abs(delta) .gt. emax .and. tinew .lt. fptf) then
            emax = abs(delta)
            imax = i
            jmax = j
            tinew = tiold
          endif
        endif

!-----------------------------------------------------------------------
!             set maximum tice to freezing and calculate fluxes
!-----------------------------------------------------------------------
        tinew = min(tinew, fptf)
        dt = tinew - tair
        qice = cssh*exp(21.8746*tinew/(tinew + 265.5))
        sub = max(c0, fds*dalt_i*wspd*(qice - qair))
        usens = dusens*dt
        ulwr = esice*(tinew + C2K)**4 - fm

!             ensure that snow sublimated does not exceed hsnol
        dha = -dts*sub
	if (-dha .gt. hsno3) then
          !nmbal = diff between snowpack and total potential sublimation
	  nmbal_sub = max(-ithkiE(n)*100.*rhoice/rhosno,
     &			   (dha + hsno3))
	  !set snow thickness change due to all available snow
	  dha = -hsno3

	  !sublimation flux to atm = sublimated snow plus sublimated ice
          sub = -dha/dts - nmbal_sub/dts
	endif

        ultnt = rhosno*slice*sub
        sub = sub*rhosno

        ftopi = dswr - ulwr - usens - ultnt

!-----------------------------------------------------------------------
!             calculate total change in snow volume on land
!             allocate snowmelt to flux for input to bucket
!-----------------------------------------------------------------------
        dhs = 0.0
        dhtot = 0.0
	if (ticeiE(n) .ge. fptf .and. ftopi .gt. 0.0) then
          dhs = -fas*ftopi
          !accumulate ice sheet melt diagnostics
	  mduriE(n)=mduriE(n)+dts
	  mextiE(n)=1
      !calculate a meltwater reservoir that must be refilled before overflow (i.e. runoff) can occur.
      !Reservoir represents amount of refreezing (specified as a snow thickness) and pore-filling
      !that is needed to bring the column of snow that has accumulated since the last melt event up to the
      !freezing point.
      !After Pfeffer_et_al_1991: RUNOFF REFREEZI1NG AND SEA LEVEL CHANGE, JGR.
      !ISSUE: what if transient melt, then re-freeze, then another melt?   Then volume of snow to be refilled
      !will be very small.

          !If melt event occurs after a previous timestep that didnt get melt: set new refreezing reservoir.
	  !This assumes that last melt made an impermeable boundary over which any new snow (hnewsno) has accumulated.
          if (meltprev(n)<0.5) then
            !get average snow temperature change to get to 0C.  Take average snow temperature since
	    !last melt event.
	    if (counter(n) .gt. 0) then
              htemp(n)=0.-htemp(n)/counter(n)
	    else
	      htemp(n)=0.
	    endif
	    !water required to bring snow column to freezing (in units of snow thickness)
	    hr=hnewsno(n)*htemp(n)*shci/lhci
	    !water required to fill porespace (in units of snow thickness)
	    pr=hnewsno(n)*(rhopc - rhoc)/rhoc
	    !total refreezing reservoir
	    hreftot(n)= hr + pr
! 	    print*, '******FIRST MELT EVENT******'
! 	    print*, 'at:',i,j
! 	    print*, 'htemp(n)=',htemp(n)
! 	    print*, 'hr=',hr
! 	    print*, 'pr=',pr
! 	    print*, 'hreftot(n)=',hreftot(n)
! 	    print*, 'eleviE(n)=',eleviE(n)
	    !rezero accumulators
	    htemp(n)=0.
            counter(n)=0.
            hnewsno(n)=0.

	  endif
! 	  print*, '********REFREEZING EVENT******'
! 	  print*, 'at:',i,j
! 	  print*, 'dhs before refreeze=',dhs
	  !set refreeze thickness to dhs or hrestot, if it is smaller
	  refreeze=min(hreftot(n),-dhs)
	  !reduce the refreezing reservoir by refreeze
	  hreftot(n)=hreftot(n)-refreeze
	  !and reduce the ice loss as well
	  dhs=dhs+refreeze
	  !let next step know that this one had melting involved
          meltprev(n)=1.
! 	 print*, 'refreeze=',refreeze
! 	 print*, 'hreftot(n)=',hreftot(n)
! 	 print*, 'dhs after refreeze=',dhs

	else
          !if not melting conditions, then:
	  !continue accumulating temperature and counter
          htemp(n)=htemp(n)+ticeiE(n)
          counter(n)=counter(n) + 1
	  !inform next swipe through that this swipe didn't get any melting.
	  meltprev(n)=0.
	endif
        !determine snow left after sublimation
	hsno_asub = hsno3 + dha
	if (-dhs .gt. hsno_asub) then
	  !nmbal = diff between snowpack and total potential melt
	  nmbal_melt = max(-ithkiE(n)*100.*rhoice/rhosno,
     &			   (dhs + hsno_asub))
	  !set snow thickness change due to all available snow
	  dhs = -hsno_asub
	endif
	!Sum change in snow thickness and accumulate effective thickness change
        hsno3 = hsno3 + dhs + dha
	!add snow/ice melt to soil moisture
        dsmis = dsmis+(nmbal_melt+dhs)*fraciE(n)
	!add sublimation to evaporation
        evapis = evapis + sub*fraciE(n)
  	!reduce icesheet thickness by any ice sublimation/melt
	ithkiE(n) = ithkiE(n) +
     &  	(nmbal_sub + nmbal_melt)*rhosno/rhoice/100.
        ithkiE(n) = max(0.0,ithkiE(n))

	hsnowiE(n,1) = hsnowiE(n,2)
	hsnowiE(n,2) = hsno3

	mbalmelt3 = mbalmelt3 + nmbal_sub + nmbal_melt
	mbalmelt(n,1) = mbalmelt(n,2)
	mbalmelt(n,2) = mbalmelt3

        usens = dswr - ultnt - ulwr
     &  	     + dhs*ffs/dts + nmbal_melt*ffs/dts
!-----------------------------------------------------------------------
!             add ice/snow covered area fluxes to land fluxes
!-----------------------------------------------------------------------
        ticeiE(n) = tinew
        dswris = dswris + dswr*fraciE(n)
        ultntis = ultntis + ultnt*fraciE(n)
        usensis = usensis + usens*fraciE(n)
        ulwris = ulwris + ulwr*fraciE(n)
      enddo
      return
      end
