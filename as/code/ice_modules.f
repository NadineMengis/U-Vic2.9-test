! source file: /net/mare/home1/eby/as/ism/ice_modules.F
      !This file contains modules that are used in
      !multi-ice sheet simulations.  Modules are
      !similar to *.h/common block constructions
      !but are more versatile as they allow for
      !allocatable arrays to be shared between
      !subroutines.

      module ice_sheet_storage

      !-------------------------------------------
      !This module contains the DDT structure in
      !which all ice sheet variables are stored.
      !During an ice sheet timestep, pointers
      !(contained in following modules) are swung
      !to point at the components of these structures,
      !one ice sheet at a time.
      !This DDT structure is called 'isvars'.
      !A higher-level DDT structure, called 'is',
      !stores multiple 'isvars' structures (one
      !per ice sheet).
      !Note: checked for match to passing
      !modules: same # of variables (Dec4/09)
      !-------------------------------------------

      save

      !Set number of ice sheets to run here.

      integer strtis,endis

      !define ALL non-local ice sheet variables in
      !DDT 'isvars'.  Define all arrays as pointers,
      !so they can be allocatable.  Currently, all
      !scalars are defined as actual variables.  In
      !the future, they could be converted to pointers,
      !but this may be redundant.

      type :: isvars

      character(120) :: isname

      integer
     &nx,
     &ny,
     &nxp,
     &nyp,
     &nlev,
     &nlevp,
     &nsed,
     &nsedp,
     &nbed,
     &nbedp

      real
     &dx0,
     &dy0,
     &dd0,
     &pi,
     &radius,
     &totarea,
     &stdparallel,
     &xoffa,
     &yoffa,
     &zlambda,
     &bedthick

      real, pointer ::
     &darea(:,:),
     &xh(:,:),
     &yh(:,:),
     &dx(:,:),
     &dy(:,:),
     &alond(:,:),
     &alatd(:,:),
     &dxu(:,:),
     &dyu(:,:),
     &dxv(:,:),
     &dyv(:,:),
     &dxc(:,:),
     &dyc(:,:)

      real, pointer ::
     &dzeta(:),
     &zeta(:),
     &dzetah(:),
     &zetah(:),
     &zsed(:),
     &zsedm(:),
     &dzsed(:),
     &zbed(:),
     &zbedm(:),
     &dzbed(:)

      integer
     &ntraca,
     &ntracb,
     &ntrace,
     &interm,
     &ioterm,
     &iunamel,
     &iuokend,
     &iuoutfine,
     &iuout2d,
     &iuout1d,
     &iuplot1d,
     &iutab,
     &iusedbud,
     &iusedtrack,
     &iures,
     &iunest,
     &iueis,
     &iueistab,
     &iubed,
     &iubedout,
     &iunh,
     &iubas,
     &iutsin,
     &iuhis,
     &iuto,
     &iupr,
     &iuta,
     &iuclim,
     &nuvsmall,
     &n1s,
     &ntr

      real
     &gmin,
     &ubmin

      integer, pointer ::
     &maskdisplay(:,:),
     &i1s(:),
     &j1s(:)

      real, pointer ::
     &d1s(:)

      real
     &dtr,
     &timebot

      real, pointer ::
     &sedtrack(:,:,:)
      integer, pointer ::
     &itrtop(:,:)

      integer
     &npoimax,
     &nlakemax,
     &ninmx,
     &ninmy

      real, pointer ::
     &xhin(:,:),
     &yhin(:,:)

      integer
     &nxin,
     &nyin

      real
     &tmelt,
     &hfus,
     &rhoice,
     &rhobed,
     &rhosed,
     &rholiq,
     &grav,
     &rhor,
     &rhoip,
     &powir,
     &powiv

      integer
     &powi,
     &powb,
     &pows

      real
     &coefbwater,
     &hwcut,
     &tramp,
     &tramps,
     &geoflux_unif,
     &geoflux_wais,
     &dtmdh,
     &condicea,
     &condiceb,
     &cheaticea,
     &cheaticeb,
     &condsed,
     &cheatsed,
     &condbed,
     &cheatbed,
     &condliq,
     &cheatliq,
     &condair,
     &rlapse,
     &totquar,
     &totpelag,
     &totdump,
     &totslump,
     &totquara,
     &totpelaga,
     &totdumpa,
     &totslumpa,
     &totsed,
     &totsedprev,
     &timesedprev,
     &sidefac,
     &crheolbno

      integer
     &nuvmax,
     &nspamax
      real, pointer ::
     &elspa(:)
      integer, pointer ::
     &ijspa(:)
      integer, pointer ::
     &mcheck(:,:)

      real, pointer ::
     &h(:,:),
     &hs(:,:),
     &hb(:,:),
     &t(:,:,:),
     &arhap(:,:),
     &s2a0(:,:),
     &heati(:,:,:),
     &heatb(:,:),
     &budgsnow(:,:),
     &budgrain(:,:),
     &budgevap(:,:),
     &budgmelt(:,:),
     &baseperc(:,:),
     &basefrml(:,:),
     &oceanmelt(:,:),
     &calvrate(:,:),
     &budgall(:,:),
     &tsurf(:,:),
     &tsurfi(:,:),
     &tsurfo(:,:),
     &w(:,:,:),
     &wa(:,:,:),
     &geoflux(:,:),

     &sedim(:,:),
     &sedimeq(:,:),
     &sedimun(:,:),
     &sedimold(:,:),
     &disturb(:,:),
     &tsed(:,:,:),
     &wsed(:,:,:),
     &heats(:,:,:),
     &quarryrate(:,:),

     &topbed(:,:),
     &topbedeq(:,:),
     &tbed(:,:,:),
     &equiload(:,:),

     &hw(:,:),
     &tw(:,:),
     &heatw(:,:),
     &sedpres(:,:)

      integer, pointer ::
     &maskh(:,:),
     &maskwater(:,:),
     &maskpres(:,:),
     &indlake(:,:),
     &npoilake(:)

      real, pointer ::
     &icedrainage(:,:)

      real, pointer ::
     &u(:,:,:),
     &v(:,:,:),
     &ua(:,:),
     &va(:,:),
     &ui(:,:),
     &vi(:,:),
     &ub(:,:),
     &vb(:,:),
     &uadv(:,:),
     &vadv(:,:),
     &hu(:,:),
     &hv(:,:),
     &crhu(:,:),
     &crhv(:,:),
     &powbu(:,:),
     &powbv(:,:),
     &fracgu(:,:),
     &fracgv(:,:),
     &thetau(:,:),
     &thetav(:,:),
     &hgu(:,:),
     &hgv(:,:),
     &fluxgrdu(:,:),
     &fluxgrdv(:,:),
     &fluxschu(:,:),
     &fluxschv(:,:),
     &dfu(:,:),
     &dfv(:,:),
     &fsedu(:,:),
     &fsedv(:,:),
     &uw(:,:),
     &vw(:,:)

      integer, pointer ::
     &masku(:,:),
     &maskv(:,:),
     &muind(:,:),
     &mvind(:,:)

      logical ::
     &firstclim,
     &firsttherm,
     &firstocn

      integer, pointer ::
     &hislist(:)

      logical
     &firstcall

      integer
     &nyearstart,
     &nyearres,
     &ifrest,
     &ice_accelerator,
     &nyearhis,
     &nyearout1d,
     &nyearout2d,
     &nyeartab,
     &nyearplot1d

      real
     &isyear0

      integer
     &nyearend,
     &nyearendin,
     &nloopstart,
     &nloopend,
     &nloopendall,
     &nyearendall

      integer
     &nhislist

      real
     &dtimetherm,
     &dtimehyd,
     &dtimebed,
     &dtimesed,
     &dtimeocn

      integer nyearsedbud

      real xglu
      integer nwrit
      real timeice
      real dtimeice
      real timeicein
      real sealev
      real rco2inter
      integer nlake
      real, pointer :: temp(:,:)

      integer, pointer ::
     &indhx(:,:,:),
     &indhy(:,:,:)

      real, pointer ::
     &weih(:,:,:)

      real, pointer ::
     &arcocn(:,:),
     &distocn(:,:),
     &distgl(:,:)

      real, pointer ::
     &akpar(:,:),
     &bkpar(:,:)

      real, pointer ::
     &dzpow1a(:,:),
     &dzpow2a(:,:),
     &dzpow3a(:)

      logical
     &firstbedrock,
     &firsticedyn,
     &firstdosparse,
     &firstsedflow,
     &firstsedbudg,
     &firstadv,
     &firstsatable,
     &firstdownscale

      real, pointer ::
     &abox(:,:)
      integer
     &ibox

      character(20) ::
     &chist

      integer, pointer ::
     &lonind(:,:,:),
     &latind(:,:,:)

      real, pointer ::
     &weid(:,:,:)

      integer, pointer ::
     &iceiind(:,:),
     &icejind(:,:)

      integer, pointer ::
     &is2sg(:,:)

      real ::
     &facice,
     &facorb,
     &facco2,
     &weirun,
     &toth,
     &tota,
     &totflow,
     &totneg

      integer ::
     &numh,
     &itera,
     &iterc

      real, pointer ::
     &albedo(:,:),
     &mext(:,:),
     &mdur(:,:)

      logical ::
     &ismrflag

      real ::
     &rtagism

      real ::
     &zclim

      real ::

     &enhancesheet,
     &zcrh0

      end type isvars

      !Define the upper-level DDT 'is', which contains 1 to 2  isvars structures.
      !Note: I tried to make 'is' allocatable but the computer then complained it
      !didn't have enough memory.  I suspect it checks to see if it has enough for
      !a plausible # of instances..?

      type(isvars), dimension(2), target :: is

      end module ice_sheet_storage

      !-------------------------------------------
      !The following modules are seen by the ice
      !sheet subroutines during an ice sheet time-
      !step.  Before a timestep all pointers contained
      !in the following modules are
      !swung to point at the correct ice sheet
      !storage arrays (contained in the above module).
      !In this way, as the ice sheet code is run
      !for each ice sheet, the variables that define
      !the ice sheet are used, and updated as
      !necessary.
      !These modules were adapted from .h files
      !originally constructed by David Pollard (hence
      !the naming convention).
      !They have been adapted and expanded to contain
      !every 'non-local' variable (i.e. every variable
      !that needs to be saved between calls to the
      !ice sheet model and that is ice-sheet-specific).
      !Note that every variable in module
      !ice_sheet_storage must have an equivalent
      !in one of the following modules.
      !-------------------------------------------

      module comicegrid

      character(120), pointer ::
     &isname

      integer, pointer ::
     &nx,
     &ny,
     &nxp,
     &nyp,
     &nlev,
     &nlevp,
     &nsed,
     &nsedp,
     &nbed,
     &nbedp

      real, pointer ::
     &dx0,
     &dy0,
     &dd0,
     &pi,
     &radius,
     &totarea,
     &stdparallel,
     &xoffa,
     &yoffa,
     &zlambda,
     &bedthick

      real, pointer ::
     &darea(:,:),
     &xh(:,:),
     &yh(:,:),
     &dx(:,:),
     &dy(:,:),
     &alond(:,:),
     &alatd(:,:),
     &dxu(:,:),
     &dyu(:,:),
     &dxv(:,:),
     &dyv(:,:),
     &dxc(:,:),
     &dyc(:,:)

      real, pointer ::
     &dzeta(:),
     &zeta(:),
     &dzetah(:),
     &zetah(:),
     &zsed(:),
     &zsedm(:),
     &dzsed(:),
     &zbed(:),
     &zbedm(:),
     &dzbed(:)

      integer, pointer ::
     &ntraca,
     &ntracb,
     &ntrace,
     &interm,
     &ioterm,
     &iunamel,
     &iuokend,
     &iuoutfine,
     &iuout2d,
     &iuout1d,
     &iuplot1d,
     &iutab,
     &iusedbud,
     &iusedtrack,
     &iures,
     &iunest,
     &iueis,
     &iueistab,
     &iubed,
     &iubedout,
     &iunh,
     &iubas,
     &iutsin,
     &iuhis,
     &iuto,
     &iupr,
     &iuta,
     &iuclim,
     &nuvsmall,
     &n1s,
     &ntr

      real, pointer ::
     &gmin,
     &ubmin

      integer, pointer ::
     &maskdisplay(:,:),
     &i1s(:),
     &j1s(:)

      real, pointer ::
     &d1s(:)

      real, pointer ::
     &dtr,
     &timebot

      real, pointer ::
     &sedtrack(:,:,:)
      integer, pointer ::
     &itrtop(:,:)

      integer, pointer ::
     &npoimax,
     &nlakemax,
     &ninmx,
     &ninmy

      real, pointer ::
     &xhin(:,:),
     &yhin(:,:)

      integer, pointer ::
     &nxin,
     &nyin

      integer, pointer ::
     &indhx(:,:,:),
     &indhy(:,:,:)

      real, pointer ::
     &weih(:,:,:)
      real, pointer ::
     &arcocn(:,:),
     &distocn(:,:),
     &distgl(:,:)

      real, pointer ::
     &akpar(:,:),
     &bkpar(:,:)

      real, pointer ::
     &dzpow1a(:,:),
     &dzpow2a(:,:),
     &dzpow3a(:)

      logical, pointer ::
     &firstbedrock,
     &firsticedyn,
     &firstdosparse,
     &firstsedflow,
     &firstsedbudg,
     &firstadv,
     &firstsatable,
     &firstdownscale

      integer, pointer::
     &ice_accelerator

      real, pointer::
     &isyear0

      real, pointer ::
     &abox(:,:)
      integer, pointer ::
     &ibox

      character(20), pointer ::
     &chist

      integer, pointer ::
     &lonind(:,:,:),
     &latind(:,:,:)

      real, pointer ::
     &weid(:,:,:)

      integer, pointer ::
     &iceiind(:,:),
     &icejind(:,:)

      integer, pointer ::
     &is2sg(:,:)

      real, pointer ::
     &albedo(:,:),
     &mext(:,:),
     &mdur(:,:)

      logical, pointer ::
     &ismrflag

      real, pointer ::
     &rtagism

      end module

      !-------------------------------------------

      module comicephys
      real, pointer ::
     &tmelt,
     &hfus,
     &rhoice,
     &rhobed,
     &rhosed,
     &rholiq,
     &grav,
     &rhor,
     &rhoip

      integer, pointer ::
     &powi,
     &powb,
     &pows
      real, pointer ::
     &powir,
     &powiv

      real, pointer ::
     &coefbwater,
     &hwcut,
     &tramp,
     &tramps,
     &geoflux_unif,
     &geoflux_wais,
     &dtmdh,
     &condicea,
     &condiceb,
     &cheaticea,
     &cheaticeb,
     &condsed,
     &cheatsed,
     &condbed,
     &cheatbed,
     &condliq,
     &cheatliq,
     &condair,
     &rlapse,
     &totquar,
     &totpelag,
     &totdump,
     &totslump,
     &totquara,
     &totpelaga,
     &totdumpa,
     &totslumpa,
     &totsed,
     &totsedprev,
     &timesedprev,
     &sidefac,
     &crheolbno

      real, pointer ::

     &enhancesheet,
     &zcrh0

      end module

      !-------------------------------------------

      module comicesparse
      integer, pointer ::
     &nuvmax,
     &nspamax
      real, pointer ::
     &elspa(:)
      integer, pointer ::
     &ijspa(:)
      integer, pointer ::
     &mcheck(:,:)
      end module

      !-------------------------------------------

      module comicebed
      end module

      !-------------------------------------------

      module comicesheet
      real, pointer ::
     &h(:,:),
     &hs(:,:),
     &hb(:,:),
     &t(:,:,:),
     &arhap(:,:),
     &s2a0(:,:),
     &heati(:,:,:),
     &heatb(:,:),
     &budgsnow(:,:),
     &budgrain(:,:),
     &budgevap(:,:),
     &budgmelt(:,:),
     &baseperc(:,:),
     &basefrml(:,:),
     &oceanmelt(:,:),
     &calvrate(:,:),
     &budgall(:,:),
     &tsurf(:,:),
     &tsurfi(:,:),
     &tsurfo(:,:),
     &w(:,:,:),
     &wa(:,:,:),
     &geoflux(:,:),

     &sedim(:,:),
     &sedimeq(:,:),
     &sedimun(:,:),
     &sedimold(:,:),
     &disturb(:,:),
     &tsed(:,:,:),
     &wsed(:,:,:),
     &heats(:,:,:),
     &quarryrate(:,:),

     &topbed(:,:),
     &topbedeq(:,:),
     &tbed(:,:,:),
     &equiload(:,:),

     &hw(:,:),
     &tw(:,:),
     &heatw(:,:),
     &sedpres(:,:)

      integer, pointer ::
     &maskh(:,:),
     &maskwater(:,:),
     &maskpres(:,:),
     &indlake(:,:),
     &npoilake(:)

      real, pointer ::
     &icedrainage(:,:)

      real, pointer ::
     &u(:,:,:),
     &v(:,:,:),
     &ua(:,:),
     &va(:,:),
     &ui(:,:),
     &vi(:,:),
     &ub(:,:),
     &vb(:,:),
     &uadv(:,:),
     &vadv(:,:),
     &hu(:,:),
     &hv(:,:),
     &crhu(:,:),
     &crhv(:,:),
     &powbu(:,:),
     &powbv(:,:),
     &fracgu(:,:),
     &fracgv(:,:),
     &thetau(:,:),
     &thetav(:,:),
     &hgu(:,:),
     &hgv(:,:),
     &fluxgrdu(:,:),
     &fluxgrdv(:,:),
     &fluxschu(:,:),
     &fluxschv(:,:),
     &dfu(:,:),
     &dfv(:,:),
     &fsedu(:,:),
     &fsedv(:,:),
     &uw(:,:),
     &vw(:,:)

      integer, pointer ::
     &masku(:,:),
     &maskv(:,:),
     &muind(:,:),
     &mvind(:,:)

      logical, pointer ::
     &firstclim,
     &firsttherm,
     &firstocn

      integer, pointer ::
     &hislist(:)

      logical, pointer ::
     &firstcall

      integer, pointer ::
     &nyearstart,
     &nyearres,
     &ifrest,
     &nyearhis,
     &nyearout1d,
     &nyearout2d,
     &nyeartab,
     &nyearplot1d

      integer, pointer ::
     &nyearend,
     &nyearendin,
     &nloopstart,
     &nloopend,
     &nloopendall,
     &nyearendall

      integer, pointer ::
     &nhislist
      real, pointer ::
     &dtimetherm,
     &dtimehyd,
     &dtimebed,
     &dtimesed,
     &dtimeocn

      integer, pointer ::
     &nyearsedbud

      real, pointer ::
     &xglu
      integer, pointer ::
     &nwrit
      real, pointer ::
     &timeice,
     &timeicein
      real, pointer ::
     &dtimeice
      real, pointer ::
     &sealev
      real, pointer ::
     &rco2inter
      integer, pointer ::
     &nlake
      real, pointer ::
     &temp(:,:)
      real, pointer ::
     &facice,
     &facorb,
     &facco2
      real, pointer ::
     &weirun
      real, pointer ::
     &toth,
     &tota

      integer, pointer ::
     &numh,
     &itera,
     &iterc
      real, pointer ::
     &totflow,
     &totneg

      real, pointer ::
     &zclim

      end module

