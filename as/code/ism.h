! source file: /net/mare/home1/eby/as/ism/ism.h
      !UVic ESCM-derived ice sheet boundary condition arrays (imt x jmt)
      real
     & nmbal            !net mass balance accumulated by EMBM (effective cm snow)
     &,nmbal_sno	!mass balance due to snowfall (cm snow)
     &,ismtice          !ice surface temperature

      !drainage-basin arrays
      real
     & isflxm           !flux of ice sheet-derived water to ocean
     &,isflxh           !equivalent flux of heat to ocean

      !ice volume
      real
     & dicevol         !difference in pre- and post- ice sheet basin volumes (m^3)
     &,disvol          !total change in ice sheet volume (m^3)
     &,diceheat	       !latent heat equivalent of dicevol (J)
     &,basinvol        !volume of ice in drainage basins (m^3)
     &,predynvol       !volume of ice before budget addition (m^3)
     &,isvol           !global ice sheet volume calculated on ice sheet grid (m^3)
     &,nisvol          !NH ice sheet volume calculated on ice sheet grid (m^3)
     &,sisvol          !SH ice sheet volume calculated on ice sheet grid (m^3)
     &,nissmb          !NH net surface mass balance (m^3/year)
     &,sissmb          !SH net surface masss balance (m^3/year)
     &,nisiar          !NH ice area (m^2)
     &,sisiar          !SH ice area (m^2)
     &,nismext         !NH melt extent (m^2)
     &,sismext         !SH melt extent (m^2)
     &,nisaccum        !NH accumulation (m^3 ice)
     &,sisaccum        !SH accumulation (m^3 ice)
     &,nismelt         !NH melt (m^3 ice)
     &,sismelt         !SH melt (m^3 ice)
      real
     &iac              !counter recording number of additions sbc conditions
     &,dtism           !ice sheet timestep (years)

      real avgie       !climate-grid-average ice elevation

      integer
     & nisb            !ice sheet basin index array
     &,drmap
     &,isb             !cross-reference array from ASCII to ice sheet basins

      character(120) :: cmd
      integer ier,ishell

      parameter(nisb=200)

      logical iceinit
      common /ism_couple_l/
     & iceinit

      common /ism_couple_r/
     & nmbal(imt,jmt,2),nmbal_sno(imt,jmt),ismtice(imt,jmt)
     &,isflxm(nisb),isflxh(nisb),predynvol(nisb)
     &,disvol,dicevol(nisb),diceheat(nisb),isvol,basinvol(nisb)
     &,nisvol,sisvol,nissmb,sissmb,nisiar,sisiar,nismext,sismext
     &,nisaccum,sisaccum,nismelt,sismelt
     &,iac
     &,dtism
     &,avgie(imt,jmt)

      common /ism_couple_i/
     & drmap(imt,jmt),isb(nisb)

      real zclim_out,sealev_out
      common /ism_couple_r/ zclim_out,sealev_out

