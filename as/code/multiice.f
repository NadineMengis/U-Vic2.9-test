! source file: /net/mare/home1/eby/as/ism/multiice.F
      subroutine set_ice_sheet_DDT(n)

!=======================================================================
!     Set all constants and allocate storage arrays for each ice sheet
!=======================================================================

      use ice_sheet_storage

      implicit none

      integer n,ioun

      real aisifrest,aisnyearres,aisnyearhis,aisaccel,aisyear0
     &    ,aiszclim,aissealev
     &    ,aisenhancesheet,aiszcrh,aisgeoflux
      real gisifrest,gisnyearres,gisnyearhis,gisaccel,gisyear0
     &    ,giszclim,gissealev
     &    ,gisenhancesheet,giszcrh,gisgeoflux
      namelist /ais/ aisifrest,aisnyearres,aisnyearhis
     &,aisaccel,aisyear0,aiszclim,aissealev
     &,aisenhancesheet,aiszcrh,aisgeoflux
      namelist /gis/ gisifrest,gisnyearres,gisnyearhis
     &,gisaccel,gisyear0,giszclim,gissealev
     &,gisenhancesheet,giszcrh,gisgeoflux

      !Explicitly set constants for each ice sheet
      if(n.eq.1) then
        is(n)%isname = 'Antarctica'

	!!!originally in namelist input
        is(n)%nyearstart=0
        is(n)%nyearres=1
        is(n)%ifrest=1
        is(n)%dtimeice=0.1
        is(n)%ice_accelerator=1
        is(n)%isyear0=0.
        is(n)%nyearhis=1
        is(n)%nyearout1d=0
        is(n)%nyearout2d=0
        is(n)%nyeartab=0
        is(n)%nyearplot1d=0
	is(n)%zclim=1.
	is(n)%sealev=0.

        !!!originally in comicegrid.h
        is(n)%pi     = 3.14159265358979
        is(n)%radius = 6.37122e6

        is(n)%dx0=10.e3
        is(n)%dy0=is(n)%dx0
        is(n)%dd0=is(n)%dx0
        is(n)%nx=282
        is(n)%ny=282
	is(n)%stdparallel=-71.
	is(n)%xoffa=0.
	is(n)%yoffa=0.

        is(n)%nxp=is(n)%nx+1
        is(n)%nyp=is(n)%ny+1
        is(n)%nlev=10
        is(n)%nlevp=is(n)%nlev+1
        is(n)%nsed = 1
        is(n)%nsedp =is(n)%nsed+1
        is(n)%nbed = 1
        is(n)%nbedp =is(n)%nbed+1

        is(n)%ntraca=1
        is(n)%ntracb=0
        is(n)%ntrace=is(n)%ntraca+is(n)%ntracb

        is(n)%interm=5
        is(n)%ioterm=6
        is(n)%iunamel=10
        is(n)%iuokend=11
        is(n)%iuoutfine=18
        is(n)%iuout2d=19
        is(n)%iuout1d=20
        is(n)%iuplot1d=21
        is(n)%iutab=22
        is(n)%iusedbud=23
        is(n)%iusedtrack=24
        is(n)%iures=30
        is(n)%iunest=31
        is(n)%iueis=40
        is(n)%iueistab=41
        is(n)%iubed=45
        is(n)%iubedout=0
        is(n)%iunh =47
        is(n)%iubas=48
        is(n)%iutsin=49
        is(n)%iuhis=55
        is(n)%iuto =60
        is(n)%iupr=61
        is(n)%iuta=62
        is(n)%iuclim=63

        is(n)%gmin = (0.01/is(n)%dd0)**2
        is(n)%ubmin  = 0.30            ! 0.01

        is(n)%nuvsmall = 6

        is(n)%ntr = 1
        is(n)%dtr = 1.

        is(n)%npoimax = max(1,(is(n)%nx*is(n)%ny)/5)
        is(n)%nlakemax = 50

        is(n)%ninmx=is(n)%nx+2
        is(n)%ninmy=is(n)%ny+2
	is(n)%ismrflag=.false.
        !!!

        !!!Formerly in comicephys.h
        !Basic quantities

        is(n)%tmelt = 273.15
        is(n)%hfus  = 0.335 e6

        is(n)%rhoice = 910.
        is(n)%rhobed = 3370.
        is(n)%rhosed = 2390. ! kg/m3

        is(n)%rholiq = 1028.
        is(n)%grav   = 9.81

        is(n)%rhor = is(n)%rhoice/is(n)%rholiq
        is(n)%rhoip = (1.-is(n)%rhoice/is(n)%rholiq)*is(n)%rhoice

        is(n)%powi = 3

        is(n)%powiv = (is(n)%powi-1.)/(2.*is(n)%powi)
        is(n)%powir = (1./is(n)%powi)

        is(n)%powb = 2
        is(n)%pows = is(n)%powb - 1

        is(n)%coefbwater = .001  ! Pa/(m/a)

        is(n)%hwcut = 0.2   ! can't be zero in movewater

        is(n)%tramp = .5      ! deg K (sliding, in basecoef)
        is(n)%tramps= .5       ! deg K (def sed, in sedflow)

        is(n)%geoflux_unif = .0546 * 31556926
        is(n)%geoflux_wais = .070  * 31556926

        is(n)%dtmdh = 8.66e-4

        is(n)%condicea = 2.1*31556926                     ! J/a/m/K
        is(n)%condiceb = 0.                             ! J/a/m/K
        is(n)%cheaticea= 2009.                           ! J/kg/K
        is(n)%cheaticeb= 0.                               ! J/kg/K^2

        is(n)%condsed  = 3.3*31556926                    ! J/a/m/K
        is(n)%cheatsed = 1000.                            ! J/kg/K

        is(n)%condbed  = 3.3*31556926                     ! J/a/m/K
        is(n)%cheatbed = 1000.                       ! J/kg/K

        is(n)%condliq = 70.*31556926                     ! J/a/m2/K
        is(n)%cheatliq = 4218.                            ! J/kg/K

        is(n)%condair = 10.*31556926                      ! J/a/m2/K

        is(n)%rlapse = .0080           ! lapse rate correction (K/m)

        is(n)%sidefac=0.

        is(n)%crheolbno = 1.e-20
        !!!

        !!!Formerly in comicesparse.h
        is(n)%nuvmax=2*is(n)%nx*is(n)%ny
        is(n)%nspamax=is(n)%nuvmax*9 + 2
        !!!

        !!!formerly saved logical data statements
        is(n)%firstclim = .true.
        is(n)%firsttherm = .true.
        is(n)%firstocn = .true.
        is(n)%firstcall = .true.

	is(n)%firstbedrock = .true.
	is(n)%firsticedyn = .true.
	is(n)%firstdosparse = .true.
	is(n)%firstsedflow = .true.
	is(n)%firstsedbudg = .true.
	is(n)%firstadv = .true.
	is(n)%firstsatable = .true.
	is(n)%firstdownscale=.true.
	!!!
        !read available namelist variables
        call getunit (ioun, 'control.in', 'f s r')
        read  (ioun, ais, end=101)
101     continue
        is(n)%nyearres=aisnyearres
        is(n)%ifrest=aisifrest

        is(n)%nyearhis=aisnyearhis
        is(n)%ice_accelerator=aisaccel
	is(n)%isyear0=aisyear0
	is(n)%sealev=aissealev
	is(n)%zclim=aiszclim

	is(n)%enhancesheet=aisenhancesheet
	is(n)%zcrh0=aiszcrh
	is(n)%geoflux_unif=aisgeoflux * 31556926

        call relunit (ioun)
      endif
      if (n.eq.2) then

        is(n)%isname = 'Greenland'
	!!!originally in namelist input
        is(n)%nyearstart=0
        is(n)%nyearres=1
        is(n)%ifrest=1
        is(n)%dtimeice=0.0625
        is(n)%ice_accelerator=1
	is(n)%isyear0=0.
        is(n)%nyearhis=1
        is(n)%nyearout1d=0
        is(n)%nyearout2d=0
        is(n)%nyeartab=0
        is(n)%nyearplot1d=0
	is(n)%zclim=1.
	is(n)%sealev=0.

        !!!originally in comicegrid.h
        is(n)%pi     = 3.14159265358979
        is(n)%radius = 6.37122e6

        is(n)%dx0=10.e3
        is(n)%dy0=is(n)%dx0
        is(n)%dd0 =is(n)%dx0
        is(n)%nx=76*2-1
        is(n)%ny=140*2-1
	is(n)%stdparallel=71.

        is(n)%nxp=is(n)%nx+1
        is(n)%nyp=is(n)%ny+1
        is(n)%nlev=10
        is(n)%nlevp=is(n)%nlev+1
        is(n)%nsed = 1
        is(n)%nsedp =is(n)%nsed+1
        is(n)%nbed = 1
        is(n)%nbedp =is(n)%nbed+1

        is(n)%ntraca=1
        is(n)%ntracb=0
        is(n)%ntrace=is(n)%ntraca+is(n)%ntracb

        is(n)%interm=5
        is(n)%ioterm=6
        is(n)%iunamel=10
        is(n)%iuokend=11
        is(n)%iuoutfine=18
        is(n)%iuout2d=19
        is(n)%iuout1d=20
        is(n)%iuplot1d=21
        is(n)%iutab=22
        is(n)%iusedbud=23
        is(n)%iusedtrack=24
        is(n)%iures=30
        is(n)%iunest=31
        is(n)%iueis=40
        is(n)%iueistab=41
        is(n)%iubed=45
        is(n)%iubedout=0
        is(n)%iunh =47
        is(n)%iubas=48
        is(n)%iutsin=49
        is(n)%iuhis=55
        is(n)%iuto =60
        is(n)%iupr=61
        is(n)%iuta=62
        is(n)%iuclim=63

        is(n)%gmin = (0.01/is(n)%dd0)**2
        is(n)%ubmin  = 0.30            ! 0.01

        is(n)%nuvsmall = 6

        is(n)%ntr = 1
        is(n)%dtr = 1.

        is(n)%npoimax = max(1,(is(n)%nx*is(n)%ny)/5)
        is(n)%nlakemax = 50

        is(n)%ninmx=is(n)%nx+2
        is(n)%ninmy=is(n)%ny+2
	is(n)%ismrflag=.false.
        !!!

        !!!Formerly in comicephys.h
        !Basic quantities

        is(n)%tmelt = 273.15
        is(n)%hfus  = 0.335 e6

        is(n)%rhoice = 910.
        is(n)%rhobed = 3370.
        is(n)%rhosed = 2390. ! kg/m3

        is(n)%rholiq = 1028.
        is(n)%grav   = 9.81

        is(n)%rhor = is(n)%rhoice/is(n)%rholiq
        is(n)%rhoip = (1.-is(n)%rhoice/is(n)%rholiq)*is(n)%rhoice

        is(n)%powi = 3

        is(n)%powiv = (is(n)%powi-1.)/(2.*is(n)%powi)
        is(n)%powir = (1./is(n)%powi)

        is(n)%powb = 2
        is(n)%pows = is(n)%powb - 1

        is(n)%coefbwater = .001  ! Pa/(m/a)

        is(n)%hwcut = 0.2   ! can't be zero in movewater

        is(n)%tramp = .5      ! deg K (sliding, in basecoef)
        is(n)%tramps= .5       ! deg K (def sed, in sedflow)

        is(n)%geoflux_unif = .0546 * 31556926
        is(n)%geoflux_wais = .070  * 31556926

        is(n)%dtmdh = 8.66e-4

        is(n)%condicea = 2.1*31556926                     ! J/a/m/K
        is(n)%condiceb = 0.                             ! J/a/m/K
        is(n)%cheaticea= 2009.                           ! J/kg/K
        is(n)%cheaticeb= 0.                               ! J/kg/K^2

        is(n)%condsed  = 3.3*31556926                    ! J/a/m/K
        is(n)%cheatsed = 1000.                            ! J/kg/K

        is(n)%condbed  = 3.3*31556926                     ! J/a/m/K
        is(n)%cheatbed = 1000.                       ! J/kg/K

        is(n)%condliq = 70.*31556926                     ! J/a/m2/K
        is(n)%cheatliq = 4218.                            ! J/kg/K

        is(n)%condair = 10.*31556926                      ! J/a/m2/K

        is(n)%rlapse = .0080           ! lapse rate correction (K/m)

        is(n)%sidefac=0.

        is(n)%crheolbno = 1.e-20
        !!!

        !!!Formerly in comicesparse.h
        is(n)%nuvmax=2*is(n)%nx*is(n)%ny
        is(n)%nspamax=is(n)%nuvmax*9 + 2
        !!!

        !!!formerly saved logical data statements
        is(n)%firstclim = .true.
        is(n)%firsttherm = .true.
        is(n)%firstocn = .true.
        is(n)%firstcall = .true.

	is(n)%firstbedrock = .true.
	is(n)%firsticedyn = .true.
	is(n)%firstdosparse = .true.
	is(n)%firstsedflow = .true.
	is(n)%firstsedbudg = .true.
	is(n)%firstadv = .true.
	is(n)%firstsatable = .true.
	is(n)%firstdownscale=.true.
	!!!

	!read available namelist variables
        call getunit (ioun, 'control.in', 'f s r')
        read  (ioun, gis, end=102)
102     continue
        is(n)%nyearres=gisnyearres
        is(n)%ifrest=gisifrest

        is(n)%nyearhis=gisnyearhis
	is(n)%ice_accelerator=gisaccel
	is(n)%isyear0=gisyear0
	is(n)%sealev=gissealev
	is(n)%zclim=giszclim
	is(n)%enhancesheet=gisenhancesheet
	is(n)%zcrh0=giszcrh
	is(n)%geoflux_unif=gisgeoflux * 31556926

        call relunit (ioun)
      endif
      !allocate DDT storage arrays.  Note: compared against deallocate routine: all match (Dec 4/09)

      !!!Formerly in comicegrid.h
      !h grid variables
      allocate(
     &is(n)%darea(is(n)%nx,is(n)%ny),
     &is(n)%xh(is(n)%nx,is(n)%ny),
     &is(n)%yh(is(n)%nx,is(n)%ny),
     &is(n)%dx(is(n)%nx,is(n)%ny),
     &is(n)%dy(is(n)%nx,is(n)%ny),
     &is(n)%alond(is(n)%nx,is(n)%ny),
     &is(n)%alatd(is(n)%nx,is(n)%ny),
     &is(n)%dxu(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%dyu(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%dxv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%dyv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%dxc(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%dyc(0:is(n)%nxp,0:is(n)%nyp)
     &)

      !v grid variables
      allocate(
     &is(n)%dzeta(is(n)%nlev),
     &is(n)%zeta(0:is(n)%nlevp),
     &is(n)%dzetah(0:is(n)%nlev),
     &is(n)%zetah(0:is(n)%nlev),
     &is(n)%zsed(0:is(n)%nsedp),
     &is(n)%zsedm(0:is(n)%nsed),
     &is(n)%dzsed(is(n)%nsed),
     &is(n)%zbed(0:is(n)%nbedp),
     &is(n)%zbedm(0:is(n)%nbed),
     &is(n)%dzbed(is(n)%nbed)
     &)

      !other
      allocate(
     &is(n)%maskdisplay(is(n)%nx,is(n)%ny),
     &is(n)%i1s(is(n)%nx+is(n)%ny),
     &is(n)%j1s(is(n)%nx+is(n)%ny),
     &is(n)%d1s(is(n)%nx+is(n)%ny)
     &)

      allocate(
     &is(n)%sedtrack(is(n)%nx,is(n)%ny,0:is(n)%ntr),
     &is(n)%itrtop(is(n)%nx,is(n)%ny)
     &)

      allocate(
     &is(n)%xhin(is(n)%ninmx,is(n)%ninmy),
     &is(n)%yhin(is(n)%ninmx,is(n)%ninmy)
     &)
      !!!

      !!!Formerly in comicesparse.h
      allocate(
     &is(n)%elspa(is(n)%nspamax),
     &is(n)%ijspa(is(n)%nspamax),
     &is(n)%mcheck(0:is(n)%nxp,0:is(n)%nyp)
     &)
      !!!

      !!!Formerly declared in icectl.F
      allocate(
     &is(n)%h(is(n)%nx,is(n)%ny),
     &is(n)%hs(is(n)%nx,is(n)%ny),
     &is(n)%hb(is(n)%nx,is(n)%ny),
     &is(n)%t(is(n)%nx,is(n)%ny,0:is(n)%nlevp),
     &is(n)%arhap(is(n)%nx,is(n)%ny),
     &is(n)%s2a0(is(n)%nx,is(n)%ny),
     &is(n)%heati(is(n)%nx,is(n)%ny,is(n)%nlev),
     &is(n)%heatb(is(n)%nx,is(n)%ny),
     &is(n)%budgsnow(is(n)%nx,is(n)%ny),
     &is(n)%budgrain(is(n)%nx,is(n)%ny),
     &is(n)%budgevap(is(n)%nx,is(n)%ny),
     &is(n)%budgmelt(is(n)%nx,is(n)%ny),
     &is(n)%baseperc(is(n)%nx,is(n)%ny),
     &is(n)%basefrml(is(n)%nx,is(n)%ny),
     &is(n)%oceanmelt(is(n)%nx,is(n)%ny),
     &is(n)%calvrate(is(n)%nx,is(n)%ny),
     &is(n)%budgall(is(n)%nx,is(n)%ny),
     &is(n)%tsurf(is(n)%nx,is(n)%ny),
     &is(n)%tsurfi(is(n)%nx,is(n)%ny),
     &is(n)%tsurfo(is(n)%nx,is(n)%ny),
     &is(n)%w(is(n)%nx,is(n)%ny,0:is(n)%nlevp),
     &is(n)%wa(is(n)%nx,is(n)%ny,0:is(n)%nlevp),
     &is(n)%maskh(is(n)%nx,is(n)%ny),
     &is(n)%geoflux(is(n)%nx,is(n)%ny),

     &is(n)%sedim(is(n)%nx,is(n)%ny),
     &is(n)%sedimeq(is(n)%nx,is(n)%ny),
     &is(n)%sedimun(is(n)%nx,is(n)%ny),
     &is(n)%sedimold(is(n)%nx,is(n)%ny),
     &is(n)%disturb(is(n)%nx,is(n)%ny),
     &is(n)%tsed(is(n)%nx,is(n)%ny,is(n)%nsed),
     &is(n)%wsed(is(n)%nx,is(n)%ny,is(n)%nsed),
     &is(n)%heats(is(n)%nx,is(n)%ny,is(n)%nsed),
     &is(n)%quarryrate(is(n)%nx,is(n)%ny),

     &is(n)%topbed(is(n)%nx,is(n)%ny),
     &is(n)%topbedeq(is(n)%nx,is(n)%ny),
     &is(n)%tbed(is(n)%nx,is(n)%ny,is(n)%nbed),
     &is(n)%equiload(is(n)%nx,is(n)%ny),

     &is(n)%hw(is(n)%nx,is(n)%ny),
     &is(n)%tw(is(n)%nx,is(n)%ny),
     &is(n)%heatw(is(n)%nx,is(n)%ny),
     &is(n)%maskwater(is(n)%nx,is(n)%ny),
     &is(n)%maskpres(is(n)%nx,is(n)%ny),
     &is(n)%sedpres(is(n)%nx,is(n)%ny),
     &is(n)%indlake(is(n)%npoimax,is(n)%nlakemax),
     &is(n)%npoilake(is(n)%nlakemax)
     &)

       allocate(
     &is(n)%icedrainage(is(n)%nx,is(n)%ny),
     &is(n)%temp(is(n)%nx,is(n)%ny)
     &)

       allocate(
     &is(n)%u(0:is(n)%nxp,0:is(n)%nyp,0:is(n)%nlevp),
     &is(n)%v(0:is(n)%nxp,0:is(n)%nyp,0:is(n)%nlevp),
     &is(n)%ua(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%va(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%ui(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%vi(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%ub(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%vb(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%uadv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%vadv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%hu(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%hv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%masku(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%maskv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%crhu(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%crhv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%powbu(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%powbv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%fracgu(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%fracgv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%muind(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%mvind(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%thetau(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%thetav(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%hgu(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%hgv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%fluxgrdu(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%fluxgrdv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%fluxschu(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%fluxschv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%dfu(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%dfv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%fsedu(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%fsedv(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%uw(0:is(n)%nxp,0:is(n)%nyp),
     &is(n)%vw(0:is(n)%nxp,0:is(n)%nyp)
     &)

      !storage of all arrays saved originally
      !in save statements in lower subroutines
      allocate(
     &is(n)%indhx(4,is(n)%nx,is(n)%ny),
     &is(n)%indhy(4,is(n)%nx,is(n)%ny),
     &is(n)%weih(4,is(n)%nx,is(n)%ny),
     &is(n)%arcocn(is(n)%nx,is(n)%ny),
     &is(n)%distocn(is(n)%nx,is(n)%ny),
     &is(n)%distgl(is(n)%nx,is(n)%ny),
     &is(n)%akpar(0:2,is(n)%nlev),
     &is(n)%bkpar(0:2,is(n)%nlev),
     &is(n)%dzpow1a(is(n)%nlev,2),
     &is(n)%dzpow2a(is(n)%nlev,2),
     &is(n)%dzpow3a(is(n)%nlev)
     &)

      allocate(
     &is(n)%lonind(is(n)%nx,is(n)%ny,2),
     &is(n)%latind(is(n)%nx,is(n)%ny,2),
     &is(n)%weid(is(n)%nx,is(n)%ny,4)
     &)

      allocate(
     & is(n)%iceiind(is(n)%nx,is(n)%ny),
     & is(n)%icejind(is(n)%nx,is(n)%ny)
     &)

      allocate(
     & is(n)%is2sg(is(n)%nx,is(n)%ny)
     &)

      allocate(
     & is(n)%albedo(is(n)%nx,is(n)%ny),
     & is(n)%mext(is(n)%nx,is(n)%ny),
     & is(n)%mdur(is(n)%nx,is(n)%ny)
     &)

      allocate(
     &is(n)%abox(-800:800,-800:800)
     &,is(n)%hislist(1000)
     &)

      return
      end

      subroutine deallocate_ice_pointers

!=======================================================================
!     Deallocate storage arrays at end of model run
!=======================================================================

      use ice_sheet_storage

      implicit none

      integer n

      do n=strtis,endis
        deallocate(
     &  is(n)%darea,
     &  is(n)%xh,
     &  is(n)%yh,
     &  is(n)%dx,
     &  is(n)%dy,
     &  is(n)%alond,
     &  is(n)%alatd,
     &  is(n)%dxu,
     &  is(n)%dyu,
     &  is(n)%dxv,
     &  is(n)%dyv,
     &  is(n)%dxc,
     &  is(n)%dyc
     &  )

      !v grid variables
        deallocate(
     &  is(n)%dzeta,
     &  is(n)%zeta,
     &  is(n)%dzetah,
     &  is(n)%zetah,
     &  is(n)%zsed,
     &  is(n)%zsedm,
     &  is(n)%dzsed,
     &  is(n)%zbed,
     &  is(n)%zbedm,
     &  is(n)%dzbed
     &  )

      !other
        deallocate(
     &  is(n)%maskdisplay,
     &  is(n)%i1s,
     &  is(n)%j1s,
     &  is(n)%d1s
     &  )

        deallocate(
     &  is(n)%sedtrack,
     &  is(n)%itrtop
     &  )

        deallocate(
     &  is(n)%xhin,
     &  is(n)%yhin
     &  )
      !!!

      !!!Formerly in comicesparse.h
        deallocate(
     &  is(n)%elspa,
     &  is(n)%ijspa,
     &  is(n)%mcheck
     &  )
      !!!

      !!!Formerly declared in icectl.F
        deallocate(
     &  is(n)%h,
     &  is(n)%hs,
     &  is(n)%hb,
     &  is(n)%t,
     &  is(n)%arhap,
     &  is(n)%s2a0,
     &  is(n)%heati,
     &  is(n)%heatb,
     &  is(n)%budgsnow,
     &  is(n)%budgrain,
     &  is(n)%budgevap,
     &  is(n)%budgmelt,
     &  is(n)%baseperc,
     &  is(n)%basefrml,
     &  is(n)%oceanmelt,
     &  is(n)%calvrate,
     &  is(n)%budgall,
     &  is(n)%tsurf,
     &  is(n)%tsurfi,
     &  is(n)%tsurfo,
     &  is(n)%w,
     &  is(n)%wa,
     &  is(n)%maskh,
     &  is(n)%geoflux,

     &  is(n)%sedim,
     &  is(n)%sedimeq,
     &  is(n)%sedimun,
     &  is(n)%sedimold,
     &  is(n)%disturb,
     &  is(n)%tsed,
     &  is(n)%wsed,
     &  is(n)%heats,
     &  is(n)%quarryrate,

     &  is(n)%topbed,
     &  is(n)%topbedeq,
     &  is(n)%tbed,
     &  is(n)%equiload,

     &  is(n)%hw,
     &  is(n)%tw,
     &  is(n)%heatw,
     &  is(n)%maskwater,
     &  is(n)%maskpres,
     &  is(n)%sedpres,
     &  is(n)%indlake,
     &  is(n)%npoilake
     &  )

         deallocate(
     &  is(n)%icedrainage,
     &  is(n)%temp
     &  )

         deallocate(
     &  is(n)%u,
     &  is(n)%v,
     &  is(n)%ua,
     &  is(n)%va,
     &  is(n)%ui,
     &  is(n)%vi,
     &  is(n)%ub,
     &  is(n)%vb,
     &  is(n)%uadv,
     &  is(n)%vadv,
     &  is(n)%hu,
     &  is(n)%hv,
     &  is(n)%masku,
     &  is(n)%maskv,
     &  is(n)%crhu,
     &  is(n)%crhv,
     &  is(n)%powbu,
     &  is(n)%powbv,
     &  is(n)%fracgu,
     &  is(n)%fracgv,
     &  is(n)%muind,
     &  is(n)%mvind,
     &  is(n)%thetau,
     &  is(n)%thetav,
     &  is(n)%hgu,
     &  is(n)%hgv,
     &  is(n)%fluxgrdu,
     &  is(n)%fluxgrdv,
     &  is(n)%fluxschu,
     &  is(n)%fluxschv,
     &  is(n)%dfu,
     &  is(n)%dfv,
     &  is(n)%fsedu,
     &  is(n)%fsedv,
     &  is(n)%uw,
     &  is(n)%vw
     &  )

      !storage of all arrays saved originally
      !in save statements in lower subroutines
        deallocate(
     &  is(n)%indhx,
     &  is(n)%indhy,
     &  is(n)%weih,
     &  is(n)%arcocn,
     &  is(n)%distocn,
     &  is(n)%distgl,
     &  is(n)%akpar,
     &  is(n)%bkpar,
     &  is(n)%dzpow1a,
     &  is(n)%dzpow2a,
     &  is(n)%dzpow3a
     &  )

        deallocate(
     &  is(n)%lonind,
     &  is(n)%latind,
     &  is(n)%weid
     &  )

        deallocate(
     &  is(n)%iceiind,
     &  is(n)%icejind)

        deallocate(
     &  is(n)%is2sg)

        deallocate(
     &  is(n)%albedo,
     &  is(n)%mext,
     &  is(n)%mdur
     &  )

        deallocate(
     &  is(n)%abox
     &  ,is(n)%hislist
     &  )
      enddo

      return
      end

      subroutine set_pointers(n)

!=======================================================================
!     Swing pointers to point at ice sheet n
!=======================================================================

      use ice_sheet_storage
      use comicegrid
      use comicephys
      use comicesparse
      use comicesheet

      implicit none

      integer n

      isname=>is(n)%isname

      nx=>is(n)%nx
      ny=>is(n)%ny
      nxp=>is(n)%nxp
      nyp=>is(n)%nyp
      nlev=>is(n)%nlev
      nlevp=>is(n)%nlevp
      nsed=>is(n)%nsed
      nsedp=>is(n)%nsedp
      nbed=>is(n)%nbed
      nbedp=>is(n)%nbedp

      dx0=>is(n)%dx0
      dy0=>is(n)%dy0
      dd0=>is(n)%dd0
      pi=>is(n)%pi
      radius=>is(n)%radius
      totarea=>is(n)%totarea
      stdparallel=>is(n)%stdparallel
      xoffa=>is(n)%xoffa
      yoffa=>is(n)%yoffa
      zlambda=>is(n)%zlambda
      bedthick=>is(n)%bedthick

      darea=>is(n)%darea
      xh=>is(n)%xh
      yh=>is(n)%yh
      dx=>is(n)%dx
      dy=>is(n)%dy
      alond=>is(n)%alond
      alatd=>is(n)%alatd
      dxu=>is(n)%dxu
      dyu=>is(n)%dyu
      dxv=>is(n)%dxv
      dyv=>is(n)%dyv
      dxc=>is(n)%dxc
      dyc=>is(n)%dyc

      dzeta=>is(n)%dzeta
      zeta=>is(n)%zeta
      dzetah=>is(n)%dzetah
      zetah=>is(n)%zetah
      zsed=>is(n)%zsed
      zsedm=>is(n)%zsedm
      dzsed=>is(n)%dzsed
      zbed=>is(n)%zbed
      zbedm=>is(n)%zbedm
      dzbed=>is(n)%dzbed

      ntraca=>is(n)%ntraca
      ntracb=>is(n)%ntracb
      ntrace=>is(n)%ntrace
      interm=>is(n)%interm
      ioterm=>is(n)%ioterm
      iunamel=>is(n)%iunamel
      iuokend=>is(n)%iuokend
      iuoutfine=>is(n)%iuoutfine
      iuout2d=>is(n)%iuout2d
      iuout1d=>is(n)%iuout1d
      iuplot1d=>is(n)%iuplot1d
      iutab=>is(n)%iutab
      iusedbud=>is(n)%iusedbud
      iusedtrack=>is(n)%iusedtrack
      iures=>is(n)%iures
      iunest=>is(n)%iunest
      iueis=>is(n)%iueis
      iueistab=>is(n)%iueistab
      iubed=>is(n)%iubed
      iubedout=>is(n)%iubedout
      iunh=>is(n)%iunh
      iubas=>is(n)%iubas
      iutsin=>is(n)%iutsin
      iuhis=>is(n)%iuhis
      iuto=>is(n)%iuto
      iupr=>is(n)%iupr
      iuta=>is(n)%iuta
      iuclim=>is(n)%iuclim
      nuvsmall=>is(n)%nuvsmall
      n1s=>is(n)%n1s
      ntr=>is(n)%ntr

      gmin=>is(n)%gmin
      ubmin=>is(n)%ubmin

      maskdisplay=>is(n)%maskdisplay
      i1s=>is(n)%i1s
      j1s=>is(n)%j1s

      d1s=>is(n)% d1s

      dtr=>is(n)%dtr
      timebot=>is(n)%timebot

      sedtrack=>is(n)%sedtrack

      itrtop=>is(n)%itrtop

      npoimax=>is(n)%npoimax
      nlakemax=>is(n)%nlakemax
      ninmx=>is(n)%ninmx
      ninmy=>is(n)%ninmy
      ismrflag=>is(n)%ismrflag

      xhin=>is(n)%xhin
      yhin=>is(n)%yhin

      nxin=>is(n)%nxin
      nyin=>is(n)%nyin

      tmelt=>is(n)%tmelt
      hfus=>is(n)%hfus
      rhoice=>is(n)%rhoice
      rhobed=>is(n)%rhobed
      rhosed=>is(n)%rhosed
      rholiq=>is(n)%rholiq
      grav=>is(n)%grav
      rhor=>is(n)%rhor
      rhoip=>is(n)%rhoip
      powir=>is(n)%powir
      powiv=>is(n)%powiv
      powi=>is(n)%powi
      powb=>is(n)%powb
      pows=>is(n)%pows

      coefbwater=>is(n)%coefbwater
      hwcut=>is(n)%hwcut
      tramp=>is(n)%tramp
      tramps=>is(n)%tramps
      geoflux_unif=>is(n)%geoflux_unif
      geoflux_wais=>is(n)%geoflux_wais
      dtmdh=>is(n)%dtmdh
      condicea=>is(n)%condicea
      condiceb=>is(n)%condiceb
      cheaticea=>is(n)%cheaticea
      cheaticeb=>is(n)%cheaticeb
      condsed=>is(n)%condsed
      cheatsed=>is(n)%cheatsed
      condbed=>is(n)%condbed
      cheatbed=>is(n)%cheatbed
      condliq=>is(n)%condliq
      cheatliq=>is(n)%cheatliq
      condair=>is(n)%condair
      rlapse=>is(n)%rlapse
      totquar=>is(n)%totquar
      totpelag=>is(n)%totpelag
      totdump=>is(n)%totdump
      totslump=>is(n)%totslump
      totquara=>is(n)%totquara
      totpelaga=>is(n)%totpelaga
      totdumpa=>is(n)%totdumpa
      totslumpa=>is(n)%totslumpa
      totsed=>is(n)%totsed
      totsedprev=>is(n)%totsedprev
      timesedprev=>is(n)%timesedprev
      sidefac=>is(n)%sidefac
      crheolbno=>is(n)%crheolbno

      nuvmax=>is(n)%nuvmax
      nspamax=>is(n)%nspamax

      elspa=>is(n)%elspa
      ijspa=>is(n)%ijspa
      mcheck=>is(n)%mcheck

      h=>is(n)%h
      hs=>is(n)%hs
      hb=>is(n)%hb
      t=>is(n)%t
      arhap=>is(n)%arhap
      s2a0=>is(n)%s2a0
      heati=>is(n)%heati
      heatb=>is(n)%heatb
      budgsnow=>is(n)%budgsnow
      budgrain=>is(n)%budgrain
      budgevap=>is(n)%budgevap
      budgmelt=>is(n)%budgmelt
      baseperc=>is(n)%baseperc
      basefrml=>is(n)%basefrml
      oceanmelt=>is(n)%oceanmelt
      calvrate=>is(n)%calvrate
      budgall=>is(n)%budgall
      tsurf=>is(n)%tsurf
      tsurfi=>is(n)%tsurfi
      tsurfo=>is(n)%tsurfo
      w=>is(n)%w
      wa=>is(n)%wa
      maskh=>is(n)%maskh
      geoflux=>is(n)%geoflux

      sedim=>is(n)%sedim
      sedimeq=>is(n)%sedimeq
      sedimun=>is(n)%sedimun
      sedimold=>is(n)%sedimold
      disturb=>is(n)%disturb
      tsed=>is(n)%tsed
      wsed=>is(n)%wsed
      heats=>is(n)%heats
      quarryrate=>is(n)%quarryrate

      topbed=>is(n)%topbed
      topbedeq=>is(n)%topbedeq
      tbed=>is(n)%tbed
      equiload=>is(n)%equiload

      hw=>is(n)%hw
      tw=>is(n)%tw
      heatw=>is(n)%heatw
      sedpres=>is(n)%sedpres

      maskwater=>is(n)%maskwater
      maskpres=>is(n)%maskpres
      indlake=>is(n)%indlake
      npoilake=>is(n)%npoilake

      icedrainage=>is(n)%icedrainage

      u=>is(n)%u
      v=>is(n)%v
      ua=>is(n)%ua
      va=>is(n)%va
      ui=>is(n)%ui
      vi=>is(n)%vi
      ub=>is(n)%ub
      vb=>is(n)%vb
      uadv=>is(n)%uadv
      vadv=>is(n)%vadv
      hu=>is(n)%hu
      hv=>is(n)%hv
      crhu=>is(n)%crhu
      crhv=>is(n)%crhv
      powbu=>is(n)%powbu
      powbv=>is(n)%powbv
      fracgu=>is(n)%fracgu
      fracgv=>is(n)%fracgv
      thetau=>is(n)%thetau
      thetav=>is(n)%thetav
      hgu=>is(n)%hgu
      hgv=>is(n)%hgv
      fluxgrdu=>is(n)%fluxgrdu
      fluxgrdv=>is(n)%fluxgrdv
      fluxschu=>is(n)%fluxschu
      fluxschv=>is(n)%fluxschv
      dfu=>is(n)%dfu
      dfv=>is(n)%dfv
      fsedu=>is(n)%fsedu
      fsedv=>is(n)%fsedv
      uw=>is(n)%uw
      vw=>is(n)%vw

      masku=>is(n)%masku
      maskv=>is(n)%maskv
      muind=>is(n)%muind
      mvind=>is(n)%mvind

      firstclim=>is(n)%firstclim
      firsttherm=>is(n)%firsttherm
      firstocn=>is(n)%firstocn

      hislist=>is(n)%hislist

      firstcall=>is(n)%firstcall

      nyearend=>is(n)%nyearend
      nyearendin=>is(n)%nyearendin
      nloopstart=>is(n)%nloopstart
      nloopend=>is(n)%nloopend
      nloopendall=>is(n)%nloopendall
      nyearendall=>is(n)%nyearendall

      nhislist=>is(n)%nhislist

      dtimetherm=>is(n)%dtimetherm
      dtimehyd=>is(n)%dtimehyd
      dtimebed=>is(n)%dtimebed
      dtimesed=>is(n)%dtimesed
      dtimeocn=>is(n)%dtimeocn
      nyearsedbud=>is(n)%nyearsedbud

      xglu=>is(n)%xglu
      nwrit=>is(n)%nwrit
      timeice=>is(n)%timeice
      timeicein=>is(n)%timeicein
      rco2inter=>is(n)%rco2inter
      nlake=>is(n)%nlake
      temp=>is(n)%temp

      indhx=>is(n)%indhx
      indhy=>is(n)%indhy
      weih=>is(n)%weih
      arcocn=>is(n)%arcocn
      distocn=>is(n)%distocn
      distgl=>is(n)%distgl
      akpar=>is(n)%akpar
      bkpar=>is(n)%bkpar
      dzpow1a=>is(n)%dzpow1a
      dzpow2a=>is(n)%dzpow2a
      dzpow3a=> is(n)%dzpow3a

      firstbedrock=>is(n)%firstbedrock
      firsticedyn=>is(n)%firsticedyn
      firstdosparse=>is(n)%firstdosparse
      firstsedflow=>is(n)%firstsedflow
      firstsedbudg=>is(n)%firstsedbudg
      firstadv=> is(n)%firstadv
      firstsatable=>is(n)%firstsatable
      firstdownscale=>is(n)%firstdownscale

      abox=>is(n)%abox
      ibox=>is(n)%ibox

      chist=>is(n)%chist

      lonind=>is(n)%lonind
      latind=>is(n)%latind
      weid=>is(n)%weid

      iceiind=>is(n)%iceiind
      icejind=>is(n)%icejind

      is2sg=>is(n)%is2sg

      albedo=>is(n)%albedo
      mext=>is(n)%mext
      mdur=>is(n)%mdur

      facice=>is(n)%facice
      facorb=>is(n)%facorb
      facco2=>is(n)%facco2
      weirun=>is(n)%weirun
      toth=>is(n)%toth
      tota=>is(n)%tota
      totflow=>is(n)%totflow
      totneg=>is(n)%totneg

      numh=>is(n)%numh
      itera=>is(n)%itera
      iterc=>is(n)%iterc

      nyearstart=>is(n)%nyearstart
      nyearres=>is(n)%nyearres
      ifrest=>is(n)%ifrest
      dtimeice=>is(n)%dtimeice
      ice_accelerator=>is(n)%ice_accelerator
      isyear0=>is(n)%isyear0
      nyearhis=>is(n)%nyearhis
      nyearout1d=>is(n)%nyearout1d
      nyearout2d=>is(n)%nyearout2d
      nyeartab=>is(n)%nyeartab
      nyearplot1d=>is(n)%nyearplot1d

      zclim=>is(n)%zclim
      sealev=>is(n)%sealev
      enhancesheet=>is(n)%enhancesheet
      zcrh0=>is(n)%zcrh0

      return
      end
