! source file: /net/mare/home1/eby/as/updates/fluxis.F
      subroutine fluxis(i,j,precipis,nmbal_snois,psnois,rhis)

      use subgrid

      implicit none

      include "size.h"
      include "atm.h"
      include "cembm.h"
      include "param.h"
      include "csbc.h"
      include "pconst.h"
      include "ism.h"

      real fb, fc, qmax, rate, tair, teff, telev, soiltemp
      real ssh, tmp, tsl
      integer i,j,n, index
      real precipi, psnoi, mbalaccum3
      real precipis,nmbal_snois,psnois,rhis,tsnois

      index = lf

      fb = rhoatm*shq/dts
      fc = dts/rhosno

      precipis=0.
      nmbal_snois=0.
      psnois=0.
      rhis  =0.
      tsnois=-2.

      do n=strtiE(i,j),endiE(i,j)
        precipi = 0.
	psnoi = 0.
	tmp = 0.
	mbalaccum3 = mbalaccum(n,index)
        telev = eleviE(n)
        teff = at(i,j,2,isat)
     &         - telev*rlapse*rf1*exp(max(-1.,-telev/rf2))

        ssh = cssh*exp(17.67*teff/(teff + 243.5))
        qmax = rhmax*ssh
        if (at(i,j,2,ishum) .gt. qmax) then
          tmp = fb*(at(i,j,2,ishum) - qmax)
          precipi = tmp
	  precipiE(n) = tmp
          rhis = rhis - tmp/fb*fraciE(n)
        endif

	tair = at(i,j,2,isat) - tsnois - telev*rlapse
     &                - sbciE(n)
	if (tair .le. c0) then
          psnoi=min((hsno_max-hsnowiE(n,2))/fc,precipi)
	  psnois = psnois + psnoi*fraciE(n)
	  mbalaccum3 = mbalaccum3 + fc*(precipi - psnoi)
	  nmbal_snois = nmbal_snois + (precipi - psnoi)*fraciE(n)
	  hsnowiE(n,2) = hsnowiE(n,2) + fc*psnoi
          !accumulated refreezing snowpack reservoir
          hnewsno(n)=min(hsno_max,hnewsno(n)+fc*precipi)
	endif

	precipis = precipis + precipi*fraciE(n)

	mbalaccum(n,1) = mbalaccum(n,2)
	mbalaccum(n,2) = mbalaccum3
	avgmbal(i,j) =
     &	avgmbal(i,j)+(mbalmelt(n,2)+mbalaccum(n,2))*fraciE(n)
	avgsno(i,j) = avgsno(i,j)+hsnowiE(n,2)*fraciE(n)
      enddo

      return
      end
