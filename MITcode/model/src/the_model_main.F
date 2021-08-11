C $Header: /u/gcmpack/MITgcm/model/src/the_model_main.F,v 1.122 2015/01/20 15:00:55 dgoldberg Exp $
C $Name:  $

CBOI
C
C !TITLE: MITGCM KERNEL CODE SYNOPSIS
C !AUTHORS: mitgcm developers ( support@mitgcm.org )
C !AFFILIATION: Massachussetts Institute of Technology
C !DATE:
C !INTRODUCTION: Kernel dynamical routines
C This document summarises MITgcm code under the model/ subdirectory.
C The code under model/ ( src/ and inc/ ) contains most of
C the driver routines for the baseline forms of the kernel equations in the
C MITgcm algorithm. Numerical code for much of the baseline forms of
C these equations is also under the model/ directory. Other numerical code
C used for the kernel equations is contained in packages in the pkg/
C directory tree.
C Code for auxiliary equations and alternate discretizations of the kernel
C equations and algorithm can also be found in the pkg/ directory tree.
C
C \subsection{Getting Help and Reporting Errors and Problems}
C If you have questions please subscribe and e-mail support@mitgcm.org.
C We also welcome reports of errors and inconsistencies in the code or
C in the accompanying documentation. Please feel free to send these
C to support@mitgcm.org. For further information and to review
C problems reported to support@mitgcm.org please visit http://mitgcm.org.
C
C \subsection{MITgcm Kernel Code Calling Sequence}
C \bv
C
C Invocation from WRAPPER level...
C
C  |
C  |-THE_MODEL_MAIN :: Primary driver for the MITgcm algorithm
C    |              :: Called from WRAPPER level numerical
C    |              :: code invocation routine. On entry
C    |              :: to THE_MODEL_MAIN separate thread and
C    |              :: separate processes will have been established.
C    |              :: Each thread and process will have a unique ID
C    |              :: but as yet it will not be associated with a
C    |              :: specific region in decomposed discrete space.
C    |
C    |-INITIALISE_FIXED :: Set fixed model arrays such as topography,
C    | |                :: grid, solver matrices etc..
C    | |
C    | |-INI_PARMS :: Routine to set kernel model parameters.
C    | |           :: Kernel parameters are read from file "data"
C    | |           :: in directory in which code executes.
C    | |
C    | |-PACKAGES_BOOT      :: Start up the optional package environment.
C    | |                    :: Runtime selection of active packages.
C    | |-PACKAGES_READPARMS :: read each package input parameter file
C    | | |- ${PKG}_READPARMS
C    | |
C    | |-SET_PARMS :: Finalise model parameter setting (if fct of pkg usage)
C    | |
C    | |-INI_MODEL_IO   :: Initialise Input/Output setting
C    | |  |-MNC_INIT    :: Initialise MITgcm NetCDF interface (MNC)(see pkg/mnc)
C    | |  |-MNC_CW_INIT :: Initialise MNC grid and variable types  (see pkg/mnc)
C    | |  |-MON_INIT    :: Initialises monitor package ( see pkg/monitor )
C    | |
C    | |-INI_GRID       :: Control grid array (vert. and horiz.) initialisation.
C    | | |              :: Grid arrays are held and described in GRID.h.
C    | | |-LOAD_GRID_SPACING    :: Load grid spacing (vector) from files
C    | | |-INI_VERTICAL_GRID    :: Set up vertical grid and coordinate
C    | | |-INI_CARTESIAN_GRID   :: Cartesian horiz. grid initialisation
C    | | |                      :: (calculate grid from kernel parameters).
C    | | |-INI_SPHERICAL_POLAR_GRID :: Spherical polar horiz. grid setting
C    | | |                          :: (calculate grid from kernel parameters).
C    | | |-INI_CURVILINEAR_GRID :: General orthogonal, structured horiz. grid
C    | | |                      :: initialisation; input from raw grid files
C    | | |                      :: (LONC.bin, LATC.bin, DXF.bin, ... ) or per
C    | | |                      :: face file: horizGridFile(.faceXXX.bin)
C    | | |-INI_CYLINDER_GRID    :: Cylindrical horiz. grid setting
C    | |
C    | |-LOAD_REF_FILES   :: Read-in reference vertical profiles (T,S,Rho)
C    | |-INI_EOS          :: Initialise Equation Of State (EOS) coefficients
C    | |-SET_REF_STATE    :: Set reference pressure/geopotential, reference
C    | |                  :: stratification (for implicit IGW), vertical
C    | |                  :: velocity scaling factor and anelastic ref. density
C    | |-SET_GRID_FACTORS :: Set grid factors (fct of k) for deep-atmosphere
C    | |
C    | |-INI_DEPTHS    :: Read (from "bathyFile") or set bathymetry/orography.
C    | |-INI_MASKS_ETC :: Derive horizontal and vertical cell fractions and
C    | |               :: land masking for solid-fluid boundaries.
C    | |
C    | |-PACKAGES_INIT_FIXED  :: do all packages fixed-initialisation setting
C    | | |- ${PKG}_INIT_FIXED
C    | |
C    | |-INI_GLOBAL_DOMAIN :: Initialise domain related (global) quantities.
C    | |-INI_LINEAR_PHISURF :: Set ref. surface Bo_surf
C    | |
C    | |-INI_CORI          :: Set coriolis term. zero, f-plane, beta-plane,
C    | |                   :: sphere options are coded.
C    | |-INI_CG2D          :: 2D conjugate grad solver initialisation.
C    | |-INI_CG3D          :: 3D conjugate grad solver initialisation.
C    | |
C    | |-CONFIG_SUMMARY    :: Provide synopsis of kernel setup. Includes
C    | |                   :: annotated table of kernel parameter settings.
C    | |
C    | |-PACKAGES_CHECK    :: call each package configuration checking S/R
C    | | |- ${PKG}_CHECK
C    | |
C    | |-CONFIG_CHECK      :: Check config and parameter consistency.
C    | |
C    | |-WRITE_GRID        :: write grid fields to output files
C    | |-CPL_EXCH_CONFIGS  :: exchange config with coupler-interface
C    |
C    |-CTRL_UNPACK         :: Control vector support package. see pkg/ctrl
C    |-COST_DEPENDENT_INIT :: ( see pkg/cost )
C    |
C    |-ADTHE_MAIN_LOOP :: Derivative evaluating form of main time stepping loop
C    !                 :: Automatically generated by TAMC/TAF.
C    |
C    |-THE_MAIN_LOOP   :: Main timestepping loop routine.
C    | |
C    | |-INITIALISE_VARIA :: Set the initial conditions for time evolving fields
C    | | |
C #ifdef ALLOW_AUTODIFF
C    | | |-INI_DEPTHS         \
C    | | |-CTRL_DEPTH_INI      \
C    | | |-UPDATE_MASKS_ETC     } ALLOW_DEPTH_CONTROL case
C    | | |-UPDATE_CG2D         /
C #endif
C    | | |-INI_NLFS_VARS :: Initialise all Non-Lin Free-Surf arrays (SURFACE.h)
C    | | |-INI_DYNVARS   :: Initialise to zero all DYNVARS.h arrays
C    | | |-INI_NH_VARS   :: Initialise to zero all NH_VARS.h arrays
C    | | |
C    | | |-INI_FIELDS    :: Control initialising model fields to non-zero
C    | | | |-INI_VEL     :: Initialize 3D flow field.
C    | | | |-INI_THETA   :: Set model initial temperature field.
C    | | | |-INI_SALT    :: Set model initial salinity field.
C    | | | |-INI_PSURF   :: Set model initial free-surface height/pressure.
C    | | | |-READ_PICKUP :: Read in main model pickup files to restart a run.
C    | | |
C    | | |-INI_MIXING   :: Initialise diapycnal diffusivity.
C    | | |
C    | | |-INI_FORCING  :: Set model initial forcing fields, either
C    | | |   |          :: set in-line or from file as shown here:
C    | | |   |-READ_FLD_XY_RS(zonalWindFile)
C    | | |   |-READ_FLD_XY_RS(meridWindFile)
C    | | |   |-READ_FLD_XY_RS(surfQnetFile)
C    | | |   |-READ_FLD_XY_RS(EmPmRfile)
C    | | |   |-READ_FLD_XY_RS(thetaClimFile)
C    | | |   |-READ_FLD_XY_RS(saltClimFile)
C    | | |   |-READ_FLD_XY_RS(surfQswFile)
C    | | |
C    | | |-AUTODIFF_INIT_VARIA :: (see pkg/autodiff )
C    | | |
C    | | |-PACKAGES_INIT_VARIABLES :: Does initialisation of time evolving
C    | | | | ${PKG}_INIT_VARIA     :: package data.
C    | | |
C    | | |-COST_INIT_VARIA     :: ( see pkg/cost )
C    | | |-CONVECTIVE_ADJUSTMENT_INI :: Apply conv. adjustment to initial state
C    | | |
C    | | |-CALC_R_STAR    :: Calculate the new level thickness factor (r* coord)
C    | | |-UPDATE_R_STAR  :: Update the level thickness fraction (r* coord).
C    | | |-UPDATE_SIGMA   :: Update the level thickness fraction (sigma-coord).
C    | | |-CALC_SURF_DR   :: Calculate the new surface level thickness.
C    | | |-UPDATE_SURF_DR :: Update the surface-level thickness fraction.
C    | | |
C    | | |-UPDATE_CG2D    :: Update 2D conjugate grad. for Free-Surf.
C    | | |
C    | | |-INTEGR_CONTINUITY :: Integrate the continuity Equation
C    | | | |-INTEGRATE_FOR_W :: Integrate for vertical velocity
C    | | | |-OBCS_APPLY_W    :: Open boundary package (see pkg/obcs).
C    | | | |-UPDATE_ETAH     :: Update Surface height/pressure
C    | | |
C    | | |-CALC_R_STAR    :: Calculate the new level thickness factor (r* coord)
C    | | |-CALC_SURF_DR   :: Calculate the new surface level thickness.
C    | | |
C    | | |-STATE_SUMMARY    :: Summarise model prognostic variables.
C    | | |
C    | | |-MONITOR          :: Monitor state (see pkg/monitor)
C    | | |
C    | | |-DO_STATEVARS_TAVE :: Time averaging package ( see pkg/timeave ).
C    | | |  |-TIMEAVE_STATVARS :: Accumulate main model state variables
C    | | |  |-PTRACERS_TIMEAVE :: Accumulate passive tracers variables
C    | | |
C    | | |-DO_THE_MODEL_IO  :: Controlling routine for IO
C    | | | |-WRITE_STATE    ::  Write model state variables.
C    | | | |-TIMEAVE_STATV_WRITE :: Write Time averaged output (see pkg/timeave)
C    | | | |-FIZHI_WRITE_STATE :: Write Fizhi pkg output (see pkg/fizhi)
C    | | | |-AIM_WRITE_TAVE    :: Write AIM  pkg output (see pkg/aim_v23)
C    | | | |-LAND_OUTPUT       :: Write Land pkg output (see pkg/land)
C    | | | |-OBCS_OUTPUT       :: Write OBCS pkg output (see pkg/obcs)
C    | | | |-GMREDI_OUTPUT     :: Write GM-Redi pkg output (see pkg/gmredi)
C    | | | |-KPP_OUTPUT        :: Write KPP  pkg output (see pkg/kpp)
C    | | | |-PP81_OUTPUT       :: Write PP81 pkg output (see pkg/pp81)
C    | | | |-KL10_OUTPUT       :: Write KL10 pkg output (see pkg/kl10)
C    | | | |-MY82_OUTPUT       :: Write MY82 pkg output (see pkg/my82)
C    | | | |-OPPS_OUTPUT       :: Write OPPS pkg output (see pkg/opps)
C    | | | |-GGL90_OUTPUT      :: Write GGL90 pkg output (see pkg/ggl90)
C    | | | |-SBO_CALC          :: Compute SBO diagnostics (see pkg/sbo)
C    | | | |-SBO_OUTPUT        :: Write SBO  pkg output   (see pkg/sbo)
C    | | | |-SEAICE_OUTPUT     :: Write SeaIce pkg output (see pkg/seaice)
C    | | | |-SHELFICE_OUTPUT   :: Write ShelfIce pkg output (see pkg/shelfice)
C    | | | |-BULKF_OUTPUT      :: Write Bulk-Force output (see pkg/bulK_force)
C    | | | |-THSICE_OUTPUT     :: Write ThSIce pkg output (see pkg/thsice)
C    | | | |-PTRACERS_OUTPUT   :: Write pTracers pkg output (see pkg/ptracers)
C    | | | |-MATRIX_OUTPUT     :: Write Matrix pkg output (see pkg/matrix)
C    | | | |-GCHEM_OUTPUT      :: Write Geochemistry pkg output (see pkg/gchem)
C    | | | |-CPL_OUTPUT        :: Write Coupler-Interface output (see
C    | | | |                   :: pkg/atm_compon_interf, pkg/ocn_compon_interf)
C    | | | |-LAYERS_CALC       :: Calculate layers diagnostics (see pkg/layers)
C    | | | |-LAYERS_OUTPUT     :: Write Layers pkg output (see pkg/layers)
C    | | | |-DIAGNOSTICS_WRITE :: Write pkg/diagnostics output
C    | | |
C====|>| ****************************
C====|>| BEGIN MAIN TIMESTEPPING LOOP
C====|>| ****************************
C    | |-COST_AVERAGESFIELDS :: time-averaged Cost function terms (see pkg/cost)
C    | |-PROFILES_INLOOP     :: ( see pkg/profiles )
C    | /
C    | |-MAIN_DO_LOOP    :: Open-AD case: Main timestepping loop routine
C    | \                    otherwise: just call FORWARD_STEP
C    | |
C/\  | |-FORWARD_STEP        :: Step forward a time-step ( AT LAST !!! )
C/\  | | |
C/\  | | |-AUTODIFF_INADMODE_UNSET :: Set/reset some adjoint flags
C/\  | | |-RESET_NLFS_VARS   :: Reset some Non-Lin Free-Surf vars (Adjoint)
C/\  | | |-UPDATE_R_STAR     :: Reset r-star factor variables     (Adjoint)
C/\  | | |-UPDATE_SURF_DR    :: Reset NLFS surface thickness vars (Adjoint)
C/\  | | |
C/\  | | |-PTRACERS_SWITCH_ONOFF    :: Set/reset pTracers time-stepping switch
C/\  | | |-DIAGNOSTICS_SWITCH_ONOFF :: Activate/de-activate diagnostics
C/\  | | |-DO_STATEVARS_DIAGS ( 0 ) :: fill-up state variable diagnostics
C/\  | | |
C/\  | | |-NEST_CHILD_SETMEMO :: Nesting interface
C/\  | | |-NEST_PARENT_IO_1   :: Nesting interface
C/\  | | |
C/\  | | |-LOAD_FIELDS_DRIVER :: Control loading of input fields from files
C/\  | | |
C/\  | | |-BULKF_FORCING      :: Calculate surface forcing (see pkg/bulk_force)
C/\  | | |-CHEAPAML           :: Cheap AML driver ( see pkg/cheapaml )
C/\  | | |-CTRL_MAP_FORCING   :: Control vector support package. (see pkg/ctrl)
C/\  | | |-DUMMY_IN_STEPPING  :: Autodiff package ( pkg/autodiff ).
C/\  | | |
C/\  | | |-CPL_EXPORT_MY_DATA :: Send coupling fields to coupler
C/\  | | |-CPL_IMPORT_EXTERNAL_DATA :: Receive coupling fields from coupler
C/\  | | |
C/\  | | |-OASIS_PUT     :: Oasis coupler interface
C/\  | | |-OASIS_GET     :: Oasis coupler interface
C/\  | | |
C/\  | | |-EBM_DRIVER    :: Calculate EBM type atmospheric forcing (see pkg/ebm)
C/\  | | |
C/\  | | |-DO_ATMOSPHERIC_PHYS :: Atmospheric physics computation
C/\  | | | |
C/\  | | | |-UPDATE_OCEAN_EXPORTS     :: ( see pkg/fizhi )
C/\  | | | |-UPDATE_EARTH_EXPORTS     :: ( see pkg/fizhi )
C/\  | | | |-UPDATE_CHEMISTRY_EXPORTS :: ( see pkg/fizhi )
C/\  | | | |-FIZHI_WRAPPER            :: ( see pkg/fizhi )
C/\  | | | |-STEP_FIZHI_FG            :: ( see pkg/fizhi )
C/\  | | | |-FIZHI_UPDATE_TIME        :: ( see pkg/fizhi )
C/\  | | | |
C/\  | | | |-ATM_PHYS_DRIVER          :: ( see pkg/atm_phys )
C/\  | | | |
C/\  | | | |-AIM_DO_PHYSICS           :: ( see pkg/aim_v23 )
C/\  | | |
C/\  | | |-DO_OCEANIC_PHYS     :: Oceanic (& seaice) physics computation
C/\  | | | |
C/\  | | | |-OBCS_CALC         :: Open boundary. package (see pkg/obcs).
C/\  | | | |
C/\  | | | |-FRAZIL_CALC_RHS   :: Compute FRAZIL tendencies ( see pkg/frazil )
C/\  | | | |-THSICE_MAIN       :: Thermodynamic sea-ice driver (see pkg/thsice)
C/\  | | | |-SEAICE_MODEL      :: Sea-ice model driver (see pkg/seaice )
C/\  | | | |-SEAICE_COST_SENSI   :: Sea-ice cost-function (see pkg/seaice )
C/\  | | | |-SHELFICE_THERMODYNAMICS :: Compute ShelfIce thermo (pkg/shelfice)
C/\  | | | |-ICEFRONT_THERMODYNAMICS :: Compute IceFront thermo (pkg/icefront)
C/\  | | | |
C/\  | | | |-SALT_PLUME_DO_EXCH   :: (see pkg/salt_plume )
C/\  | | | |-FREEZE_SURFACE       :: Prevent SST to fall below TFreeze
C/\  | | | |-OCN_APPLY_IMPORT     :: Apply imported fields from coupler
C/\  | | | |-EXTERNAL_FORCING_SURF:: Compute appropriately dimensioned
C/\  | | | |                      :: surface forcing terms.
C/\  | | | |-FIND_RHO_2D @ p(k)   :: Calculate [rho(T,S,p)-Rho_0] of a slice
C/\  | | | |-FIND_RHO_2D @ p(k-1) :: Calculate [rho(T,S,p)-Rho_0] of a slice
C/\  | | | |-GRAD_SIGMA           :: Calculate isoneutral gradients
C/\  | | | |-CALC_IVDC       :: Set Implicit Vertical Diffusivity for Convection
C/\  | | | |-CALC_OCE_MXLAYER        :: Diagnose Oceanic Mixed Layer depth
C/\  | | | |
C/\  | | | |-SALT_PLUME_CALC_DEPTH   :: (see pkg/salt_plume )
C/\  | | | |-SALT_PLUME_VOLFRAC      :: (see pkg/salt_plume )
C/\  | | | |-SALT_PLUME_APPLY (Temp) :: (see pkg/salt_plume )
C/\  | | | |-SALT_PLUME_APPLY (Salt) :: (see pkg/salt_plume )
C/\  | | | |-SALT_PLUME_FORCING_SURF :: (see pkg/salt_plume )
C/\  | | | |-KPP_CALC           :: Compute KPP  vertical mixing ( see pkg/kpp )
C/\  | | | |-PP81_CALC          :: Compute PP81 vertical mixing ( see pkg/pp81 )
C/\  | | | |-KL10_CALC          :: Compute KL10 vertical mixing ( see pkg/kl10 )
C/\  | | | |-MY82_CALC          :: Compute MY82 vertical mixing ( see pkg/kl10 )
C/\  | | | |-GGL90_CALC         :: Compute GGL90 vertical mixing (see pkg/ggl10)
C/\  | | | |-GMREDI_CALC_TENSOR :: Compute GM-Redi tensor ( see pkg/gmredi )
C/\  | | | |-DWNSLP_CALC_FLOW   :: Compute Down-Slope flow  (see pkg/down_slope)
C/\  | | | |-BBL_CALC_RHS       :: Compute BBL tendencies ( see pkg/bbl )
C/\  | | | |-MYPACKAGE_CALC_RHS :: Compute mypackage tendencies (pkg/mypackage)
C/\  | | | |
C/\  | | | |-GMREDI_DO_EXCH     :: ( see pkg/gmredi )
C/\  | | | |-KPP_DO_EXCH        :: ( see pkg/kpp )
C/\  | | | |-DIAGS_RHO_G        :: Compute some density related diagnostics
C/\  | | | |-DIAGS_OCEANIC_SURF_FLUX :: Diagnose oceanic surface fluxes
C/\  | | | |-SALT_PLUME_DIAGNOSTICS_FILL :: (see pkg/salt_plume )
C/\  | | | |-ECCO_PHYS          :: ( see pkg/ecco )
C/\  | | |
C/\  | | |-STREAMICE_TIMESTEP   :: ( see pkg/streamice )
C/\  | | |
C/\  | | |-GCHEM_CALC_TENDENCY  :: geochemistry driver routine (see pkg/gchem)
C/\  | | |
C/\  | | |-LONGSTEP_AVERAGE        :: Averaging state vars ( see pkg/longstep )
C/\  | | |-LONGSTEP_THERMODYNAMICS :: Step forward tracers ( see pkg/longstep )
C/\  | | |
C/\  | | |-THERMODYNAMICS       :: theta, salt + tracer equations driver.
C/\  | | | |                         (synchronous time-stepping case)
C/\  | | | |-CALC_WSURF_TR          :: Compute T & S Linear-Free-Surf correction
C/\  | | | |-PTRACERS_CALC_WSURF_TR :: Compute Tracers Linear-Free-Surf correct.
C/\  | | | |
C/\  | | | |-GMREDI_RESIDUAL_FLOW :: Get the flow field used to advect tracers
C/\  | | | |
C/\  | | | |-TEMP_INTEGRATE       :: Step forward Prognostic Eq for Temperature.
C/\  | | | | |
C/\  | | | | |-ADAMS_BASHFORTH3   :: Extrapolate tracer forward in time (AB-3)
C/\  | | | | |-ADAMS_BASHFORTH2   :: Extrapolate tracer forward in time (AB-2)
C/\  | | | | |-CALC_3D_DIFFUSIVITY :: set vertical diffusivity
C/\  | | | | |
C/\  | | | | |-GAD_SOM_ADVECT     :: Second Order Moment (SOM) advection
C/\  | | | | |-GAD_ADVECTION      :: Generalised advection driver (multi-dim
C/\  | | | | |                         advection case) (see pkg/gad).
C/\  | | | | |-CALC_ADV_FLOW      :: set 3-D flow field to advect tracer
C/\  | | | | |-APPLY_FORCING_T    :: Problem specific forcing for temperature.
C/\  | | | | |-GAD_CALC_RHS       :: Calculate Advection-Diffusion tendency terms
C/\  | | | | |
C/\  | | | | |-ADAMS_BASHFORTH3   :: Extrapolate tendency forward in time (AB-3)
C/\  | | | | |-ADAMS_BASHFORTH2   :: Extrapolate tendency forward in time (AB-2)
C/\  | | | | |-FREESURF_RESCALE_G :: Re-scale Gt for free-surface height.
C/\  | | | | |-DWNSLP_APPLY       :: Add pkg/down_slope tendency
C/\  | | | | |
C/\  | | | | |-TIMESTEP_TRACER    :: Step tracer field forward in time
C/\  | | | | |
C/\  | | | | |-GAD_IMPLICIT_R     :: Solve vertical implicit Advect-Diffus. eqn.
C/\  | | | | |-IMPLDIFF           :: Solve vertical implicit diffusion equation.
C/\  | | | | |-CYCLE_AB_TRACER    :: Cycle time-stepping arrays for tracer field
C/\  | | | | |-CYCLE_TRACER       :: Cycle time-stepping arrays for tracer field
C/\  | | | |
C/\  | | | |-SALT_INTEGRATE       :: Step forward Prognostic Eq for Salinity.
C/\  | | | | |                       same sequence of calls as in TEMP_INTEGRATE
C/\  | | | |
C/\  | | | |-PTRACERS_INTEGRATE   :: Integrate other tracer(s) (see pkg/ptracers).
C/\  | | | | |                       same sequence of calls as in TEMP_INTEGRATE
C/\  | | | | |-OBCS_APPLY_PTRACER :: Open boundary package for pTracers
C/\  | | | |
C/\  | | | |-OBCS_APPLY_TS        :: Open boundary package (see pkg/obcs ).
C/\  | | |
C/\  | | |-LONGSTEP_AVERAGE        :: Averaging state vars ( see pkg/longstep )
C/\  | | |-LONGSTEP_THERMODYNAMICS :: Step forward tracers ( see pkg/longstep )
C/\  | | |
C/\  | | |-DO_STAGGER_FIELDS_EXCHANGES :: Update overlap regions of arrays
C/\  | | |                                 Theta & Salt (implicit IGW case)
C/\  | | |
C/\  | | |-DYNAMICS       :: Momentum equations driver.
C/\  | | | |
C/\  | | | |-CALC_GRAD_PHI_SURF :: Calculate the gradient of the surface
C/\  | | | |                       Potential anomaly.
C/\  | | | |-CALC_VISCOSITY   :: Calculate net vertical viscosity
C/\  | | | |-MOM_CALC_3D_STRAIN :: Calculates the strain tensor of 3D flow field
C/\  | | | |-OBCS_COPY_UV_N   :: for Stevens bndary Conditions (see pkg/obcs)
C/\  | | | |
C/\  | | | |-CALC_PHI_HYD     :: Integrate the hydrostatic relation.
C/\  | | | |-MOM_FLUXFORM     :: Flux Form momentum eqn. (pkg/mom_fluxform)
C/\  | | | |-MOM_VECINV       :: Vector Invariant momentum eqn (pkg/mom_vecinv)
C/\  | | | |-MOM_CALC_SMAG_3D :: Calculate Smagorinsky 3D (harmonic) viscosities
C/\  | | | |-MOM_UV_SMAG_3D   :: Calculate U,V mom. tendency due to Smag 3D Visc
C/\  | | | |-TIMESTEP         :: Step horizontal momentum fields forward in time
C/\  | | | |
C/\  | | | |-MOM_U_IMPLICIT_R :: Solve implicitly vertical Adv-Diffus equation.
C/\  | | | |-IMPLDIFF         :: Solve vertical implicit diffusion equation.
C/\  | | | |-OBCS_SAVE_UV_N   :: for Stevens bndary Conditions (see pkg/obcs)
C/\  | | | |-OBCS_APPLY_UV    :: Apply Open bndary Conditions to provisional U,V
C/\  | | | |-IMPLDIFF         :: (CD-Scheme) Solve vertical impl. diffus. eqn
C/\  | | | |
C/\  | | | |-CALC_GW          :: Vert. momentum tendency terms (Non-Hydrostatic)
C/\  | | | | |-MOM_W_SMAG_3D  :: Calculate W mom. tendency due to Smag 3D Visc
C/\  | | | |-TIMESTEP_WVEL    :: Step vert mom forward in time (Non-Hydrostatic)
C/\  | | |
C/\  | | |-MNC_UPDATE_TIME    :: Update MNC time record (see pkg/mnc)
C/\  | | |
C/\  | | |-UPDATE_R_STAR  :: Update the level thickness fraction (r* coord).
C/\  | | |-UPDATE_SIGMA   :: Update the level thickness fraction (sigma-coord).
C/\  | | |-UPDATE_R_STAR  :: Update the level thickness fraction.
C/\  | | |-UPDATE_SURF_DR :: Update the surface-level thickness fraction.
C/\  | | |-UPDATE_CG2D    :: Update 2D conjugate grad. for Free-Surf.
C/\  | | |
C/\  | | |-SHAP_FILT_APPLY_UV  :: Apply Shapiro Filter to provisional velocity
C/\  | | |-ZONAL_FILT_APPLY_UV :: Apply  Zonal  Filter to provisional velocity
C/\  | | |
C/\  | | |-SOLVE_FOR_PRESSURE  :: Find surface pressure.
C/\  | | | |-CALC_DIV_GHAT     :: Form the RHS of the surface pressure eqn.
C/\  | | | |-CG2D              :: Two-dim pre-con. conjugate-gradient.
C/\  | | | |-PRE_CG3D          :: Finish to set the RHS of the 3-D pressure eqn.
C/\  | | | |-CG3D              :: Three-dim pre-con. conjugate-gradient solver.
C/\  | | | |-POST_CG3D         :: finalise solution of NH and Free-Surf pressure
C/\  | | |
C/\  | | |-MOMENTUM_CORRECTION_STEP :: Finalise momentum stepping
C/\  | | | |-CALC_GRAD_PHI_SURF  :: Return DDx and DDy of surface pressure
C/\  | | | |-CORRECTION_STEP     :: Pressure correction to momentum
C/\  | | | |-OBCS_APPLY_UV       :: Open boundary package (see pkg/obcs).
C/\  | | | |-SHAP_FILT_APPLY_UV  :: Apply Shapiro Filter to latest velocity
C/\  | | | |-ZONAL_FILT_APPLY_UV :: Apply  Zonal  Filter to latest velocity
C/\  | | |
C/\  | | |-INTEGR_CONTINUITY   :: Integrate continuity equation (see above)
C/\  | | |
C/\  | | |-CALC_R_STAR    :: Calculate the new level thickness factor (r* coord)
C/\  | | |-CALC_SURF_DR   :: Calculate the new surface level thickness.
C/\  | | |
C/\  | | |-DO_STAGGER_FIELDS_EXCHANGES :: Update overlap regions of arrays
C/\  | | |                             uVel,vVel & wVel (stagger-time-step case)
C/\  | | |
C/\  | | |-DO_STATEVARS_DIAGS ( 1 ) :: fill-up state variable diagnostics
C/\  | | |
C/\  | | |-THERMODYNAMICS       :: theta, salt + tracer Eq. driver (see above).
C/\  | | |                         (staggered time-stepping case)
C/\  | | |
C/\  | | |-TRACERS_CORRECTION_STEP :: Finalise tracer stepping:
C/\  | | | |                       ::  apply filter, conv.adjustment
C/\  | | | |-TRACERS_IIGW_CORRECTION :: apply Implicit IGW adjustment to T & S
C/\  | | | |-SHAP_FILT_APPLY_TS        :: Apply Shapiro Filter to latest T & S
C/\  | | | |-ZONAL_FILT_APPLY_TS       :: Apply  Zonal  Filter to latest T & S
C/\  | | | |-PTRACERS_ZONAL_FILT_APPLY :: Apply  Zonal Filter to pTracers
C/\  | | | |-SALT_FILL                 :: Fill up negative Salt
C/\  | | | |-OPPS_INTERFACE            :: ( see pkg/opps )
C/\  | | | |-CONVECTIVE_ADJUSTMENT     :: Apply convective adjustment
C/\  | | | |-MATRIX_STORE_TENDENCY_IMP :: ( see pkg/matrix )
C/\  | | |
C/\  | | |-LONGSTEP_AVERAGE        :: Averaging state vars ( see pkg/longstep )
C/\  | | |-LONGSTEP_THERMODYNAMICS :: Step forward tracers ( see pkg/longstep )
C/\  | | |
C/\  | | |-GCHEM_FORCING_SEP :: Tracer forcing for gchem pkg (if tracer
C/\  | | |                   :: dependent tendencies calculated separately)
C/\  | | |
C/\  | | |-DO_FIELDS_BLOCKING_EXCHANGES :: Sync up overlap regions.
C/\  | | |
C/\  | | |-DO_STATEVARS_DIAGS ( 2 ) :: fill-up state variable diagnostics
C/\  | | |
C/\  | | |-GRIDALT_UPDATE    :: ( see pkg/gridalt )
C/\  | | |-STEP_FIZHI_CORR   :: ( see pkg/fizhi )
C/\  | | |
C/\  | | |-FLT_MAIN          :: Step forward Floats (see pkg/flt)
C/\  | | |
C/\  | | |-DO_STATEVARS_TAVE :: Time averaging package (see above)
C/\  | | |
C/\  | | |-NEST_PARENT_IO_2  :: Nesting interface
C/\  | | |-NEST_CHILD_TRANSP :: Nesting interface
C/\  | | |
C/\  | | |-MONITOR          :: Monitor package (pkg/monitor).
C/\  | | |
C/\  | | |-COST_TILE        :: ( see pkg/cost )
C/\  | | |
C/\  | | |-DO_THE_MODEL_IO  :: Controlling routine for IO (see above)
C/\  | | |
C/\  | | |-PTRACERS_RESET   :: Re-initialize PTRACERS ( see pkg/ptracers )
C/\  | | |
C/\  | | |-DO_WRITE_PICKUP  :: Controlling routine for writing files to restart
C/\  | | | |-PACKAGES_WRITE_PICKUP :: Write pickup files for each package
C/\  | | | | |                     ::  which needs it to restart
C/\  | | | | |-GAD_WRITE_PICKUP       :: Write Generic AdvDiff pickups for SOM
C/\  | | | | |                        :: advection scheme (pkg/generic_advdiff)
C/\  | | | | |-CD_CODE_WRITE_PICKUP   :: Write CD-code pickups (see pkg/cd_code)
C/\  | | | | |-OBCS_WRITE_PICKUP      :: Write OBCS pickups    (see pkg/obcs)
C/\  | | | | |-GGL90_WRITE_PICKUP     :: Write GGL90 pickups   (see pkg/ggl90)
C/\  | | | | |-BBL_WRITE_PICKUP       :: Write BBL pickups     (see pkg/bbl)
C/\  | | | | |-CHEAPAML_WRITE_PICKUP  :: Write CheapAML pickups (pkg/cheapaml)
C/\  | | | | |-FLT_WRITE_PICKUP       :: Write Floats pickups  (see pkg/flt)
C/\  | | | | |-PTRACERS_WRITE_PICKUP  :: Write pTracers pickups (pkg/ptracers)
C/\  | | | | |-GCHEM_WRITE_PICKUP     :: Write Geo-Chem pickups (see pkg/gchem)
C/\  | | | | |-SEAICE_WRITE_PICKUP    :: Write SeaIce pickups  (see pkg/seaice)
C/\  | | | | |-STREAMICE_WRITE_PICKUP :: Write StreamIce pickups (pkg/streamice)
C/\  | | | | |-SHELFICE_WRITE_PICKUP :: Write ShelfIce pickups (pkg/shelfice)
C/\  | | | | |-THSICE_WRITE_PICKUP    :: Write ThSIce pickups  (see pkg/thsice)
C/\  | | | | |-LAND_WRITE_PICKUP      :: Write Land   pickups  (see pkg/land)
C/\  | | | | |-ATM_PHYS_WRITE_PICKUP  :: Write Atm-Phys pickups (pkg/atm_phys)
C/\  | | | | |-FIZHI_WRITE_PICKUP     :: Write Fizhi pickups   (see pkg/fizhi)
C/\  | | | | |-FIZHI_WRITE_VEGTILES   :: Write Fizhi VegTiles  (see pkg/fizhi)
C/\  | | | | |-FIZHI_WRITE_DATETIME   :: Write Fizhi DateTime  (see pkg/fizhi)
C/\  | | | | |-CPL_WRITE_PICKUP       :: Write Coupling-Interface pickups
C/\  | | | | |-MYPACKAGE_WRITE_PICKUP :: Write pkg/mypackage pickups
C/\  | | | |
C/\  | | | |-WRITE_PICKUP          :: Write main model pickup files.
C/\  | | |
C/\  | | |-AUTODIFF_INADMODE_SET   :: Set/reset some adjoint flags
C    | |
C<===|=| **************************
C<===|=| END MAIN TIMESTEPPING LOOP
C<===|=| **************************
C    | |
C    | |-COST_AVERAGESFIELDS :: Time-averaged Cost function terms (see pkg/cost)
C    | |-PROFILES_INLOOP     :: ( see pkg/profiles )
C    | |-COST_FINAL          :: Cost function package (see pkg/cost)
C    |
C    |-CTRL_PACK       :: Control vector support package (see pkg/ctrl)
C    |
C    |-GRDCHK_MAIN     :: Gradient check package (see pkg/grdchk)
C    |
C    |-TIMER_PRINTALL  :: Computational timing summary
C    |
C    |-COMM_STATS      :: Summarise inter-proc and inter-thread communication
C    |                 :: events.
C \ev
C
CEOI

#include "PACKAGES_CONFIG.h"      Created by convert_cpp_cmd2defines with the following command line arguments: 定义了哪些包
      例子：
      #undef  ALLOW_STREAMICE
      #undef  ALLOW_THSICE
      #undef  ALLOW_ZONAL_FILT
      /*   */
      /*  Enabled packages:  */
      #define ALLOW_CAL
      #define ALLOW_DEBUG
      #define ALLOW_DIAGNOSTICS
      #define ALLOW_DIC
      #define ALLOW_EXCH2
      #define ALLOW_EXF
      例子：

#include "CPP_OPTIONS.h"  模型的主CPP选项文件：控制在model/src代码中编译哪些可选特性。
#include "AD_CONFIG.h"    Created by convert_cpp_cmd2defines with the following command line arguments:  
                          Forward version -bAD_CONFIG_H -UALLOW_ADJOINT_RUN -UALLOW_TANGENTLINEAR_RUN
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"    用于AutoDiff（AutoDiff）包的CPP options文件：控制在此包代码中编译哪些可选特性
#endif
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"        CPP options file for Control（ctrl）包：控制在此包代码中编译哪些可选功能
#endif

CBOP
C     !ROUTINE: THE_MODEL_MAIN

C     !INTERFACE:
      SUBROUTINE THE_MODEL_MAIN(myThid)        通过传入myThid确定哪一个线程并进行计算

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE THE_MODEL_MAIN
C     | o Master controlling routine for model using the MITgcm     使用MITgcm UV并行包装器的模型主控程序
C     |   UV parallel wrapper.
C     *==========================================================*
C     | THE_MODEL_MAIN is invoked by the MITgcm UV parallel
C     | wrapper with a single integer argument "myThid". This
C     | variable identifies the thread number of an instance of
C     | THE_MODEL_MAIN. Each instance of THE_MODEL_MAIN works
C     | on a particular region of the models domain and
C     | synchronises with other instances as necessary. The
C     | routine has to "understand" the MITgcm parallel
C     | environment and the numerical algorithm. Editing this
C     | routine is best done with some knowledge of both aspects.
用一个整数参数“myThid”调用MITgcm UV并行包装器的THE MODEL MAIN。此变量标识THE MODEL MAIN实例的线程号。
模型的每个实例主要工作在模型域的特定区域，并根据需要与其他实例同步。例程必须“理解”MITgcm并行环境和数值算法。
编辑这个例程最好是在两个方面都有所了解的情况下完成的。
C     | Notes
C     | =====
C     | C*P* comments indicating place holders for which code is
C     |      presently being developed.
标明代码正在发展
C     *==========================================================*
C     \ev


C     !CALLING SEQUENCE:
C     THE_MODEL_MAIN()
C       |
C       |
C       |--INITIALISE_FIXED     初始化固定---设置模型配置（固定矩阵大小等）--例如地形，水文，时间步长，网格等
C       |   o Set model configuration (fixed arrays)
C       |     Topography, hydrography, timestep, grid, etc..
C       |
C       |--CTRL_UNPACK      o Derivative mode. Unpack control vector.  导数/微分模式。打开控制向量。   
C       |
C       |--ADTHE_MAIN_LOOP  o Main timestepping loop for combined  用于集成了向前和向后模式积分的主时间步进循环  
C       |                     prognostic and reverse mode integration.
C       |
C       |--THE_MAIN_LOOP    o Main timestepping loop for pure prognostic  纯向前积分的主时间步进循环
C       |                     integration.
C       |
C       |--CTRL_PACK        o Derivative mode. Unpack control vector.   导数/微分模式。关闭控制向量。
C       |
C       |--GRDCHK_MAIN      o Gradient check control routine.   梯度检查控制程序
C       |
C       |--TIMER_PRINTALL   o Print out timing statistics.   打印计时统计数据
C       |
C       |--COMM_STATS       o Print out communication statistics.   打印通信统计数据。

C     !USES:
      IMPLICIT NONE

C     == Global variables === 全局变量
#include "SIZE.h"      定义模型网格大小，并行tiles数量等
#include "EEPARAMS.h"  “执行环境”的参数。这些都由特定的数值模型和执行环境支持例程使用。
                       线程数、读取“模型”参数文件的单元号、错误计数器、err信息、父子嵌套接口等
#include "PARAMS.h"    定义模型“参数”的头文件。模型标准输入文件中的值存储在这里保存的变量中。描述参数的注释也可以在这里找到。
                       例如data里面的参数的定义，data.pkg里边的use的定义，各种文件的文件名char等
#include "DYNVARS.h"   动态的模型变量，在计算中循环的
                       %   _RL  etaN  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)    
                       %   _RL  uVel (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
                       %   _RL  theta(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
#ifdef ALLOW_AUTODIFF   是一个pkg，涉及到自动微分，AD，在手册第七章，为了配置与AD相关的设置，需要启用以下软件包：pkg/autodiff • pkg/ctrl • pkg/cost • pkg/grdchk
# include "tamc.h"   此标头配置关于3级检查点的时间步进循环的拆分
#endif
#ifdef ALLOW_CTRL
# include "ctrl.h"
# include "optim.h"
#endif

与处理成本函数（大J，因变量，大J（u）=大J（M（u）））有关的所有方面(参数设置、初始化、累加、最终评估）由包pkg/cost控制。
与自变量处理相关的方面由包pkg/ctrl控制，将在下一节中讨论（这是第七节））


C     !INPUT/OUTPUT PARAMETERS:                输入输出参数
C     == Routine arguments ==
C     myThid :: Thread number for this instance of the routine.   此实例的线程号
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     Note: Under the multi-threaded model myIter and myTime are local
C           variables passed around as routine arguments.
C           Although this is fiddly it saves the need to impose
C           additional synchronisation points when they are updated.
在多线程模型下，myIter和myTime是作为例程参数传递的局部变量。虽然这是微小难弄的的，但当它们被更新时节省了需要附加的同步点。
C     myTime :: Time counter for this thread      此线程的时间计数器
C     myIter :: Iteration counter for this thread 此线程的迭代计数器
      INTEGER myIter
      _RL     myTime
      LOGICAL exst
      LOGICAL lastdiva
CEOP

C--   set default:
      exst     = .TRUE.
      lastdiva = .TRUE.

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('THE_MODEL_MAIN',myThid)
#endif

#if defined(USE_PAPI) || defined(USE_PCL_FLOPS_SFP) || defined(USE_PCL_FLOPS) || defined(USE_PCL)
      CALL TIMER_CONTROL('','INIT','THE_MODEL_MAIN',myThid)
#endif
C--   This timer encompasses the whole code   这个计时器包含了整个代码----这里开始的，在后面那个是结束的
      CALL TIMER_START('ALL                    [THE_MODEL_MAIN]',myThid)
7
#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_CALL('INITIALISE_FIXED',myThid)
#endif
C--   Set model configuration (fixed arrays)  设置模型配置 fixed arrays
      CALL TIMER_START('INITIALISE_FIXED       [THE_MODEL_MAIN]',myThid)
      CALL INITIALISE_FIXED( myThid )
      CALL TIMER_STOP ('INITIALISE_FIXED       [THE_MODEL_MAIN]',myThid)

      myTime = startTime
      myIter = nIter0

#if ( defined (ALLOW_ADMTLM) )                      最外层判断 如果定义了ALLOW_ADMTLM    TIMER_START和TIMER_STOP关注一下

      STOP 'should never get here; ADMTLM_DSVD calls ADMTLM_DRIVER'

#elif ( defined (ALLOW_AUTODIFF))      ---实际的最外层判断-如果定义了ALLOW_AUTODIFF

# ifndef EXCLUDE_CTRL_PACK        如果没有定义排除 ctrl pack
      IF ( useCTRL ) THEN         如果用ctrl
         inquire( file='costfinal', exist=exst )
         IF ( .NOT. exst ) THEN---------------------------考虑是否调用ctrl——unpack
            IF ( (optimcycle.NE.0 .OR. .NOT.doinitxx)
     &           .AND. doMainUnpack ) THEN
               CALL TIMER_START('CTRL_UNPACK   [THE_MODEL_MAIN]',myThid)
               CALL CTRL_UNPACK( .TRUE. , myThid )        调用ctrl——unpack
               CALL TIMER_STOP ('CTRL_UNPACK   [THE_MODEL_MAIN]',myThid)
            ENDIF
         ENDIF
       ENDIF
# endif /* EXCLUDE_CTRL_PACK */

# ifdef ALLOW_COST              如果定义了ALLOW COST
      CALL COST_DEPENDENT_INIT ( myThid )  调用cost init初始化
# endif

# if ( defined (ALLOW_TANGENTLINEAR_RUN) )        ---如果定义了ALLOW_TANGENTLINEAR_RUN  （TLE）就跑G_THE_MAIN_LOOP

#  ifdef ALLOW_DEBUG           如果定义了debug，就检查G_THE_MAIN_LOOP
      IF (debugMode) CALL DEBUG_CALL('G_THE_MAIN_LOOP',myThid)
#  endif
      CALL TIMER_START('G_THE_MAIN_LOOP           [TANGENT RUN]',myThid)
      CALL G_THE_MAIN_LOOP ( myTime, myIter, myThid )                        直接跑 G_THE_MAIN_LOOP
      CALL TIMER_STOP ('G_THE_MAIN_LOOP           [TANGENT RUN]',myThid)

# elif ( defined (ALLOW_ADJOINT_RUN) || \        --else if 若是定义了ALLOW_ADJOINT_RUN 和 ALLOW_ECCO_OPTIMIZATION 
         defined (ALLOW_ECCO_OPTIMIZATION) )

#  ifdef ALLOW_DIVIDED_ADJOINT                   -还定义了ALLOW_DIVIDED_ADJOINT 就跑MDTHE_MAIN_LOOP
C-- The following assumes the TAF option '-pure'
      inquire( file='costfinal', exist=exst )
      IF ( .NOT. exst) THEN
#   ifdef ALLOW_DEBUG     如果定义了debug，就检查MDTHE_MAIN_LOOP
         IF (debugMode) CALL DEBUG_CALL('MDTHE_MAIN_LOOP',myThid)
#   endif
         CALL TIMER_START('MDTHE_MAIN_LOOP            [MD RUN]', myThid)
         CALL MDTHE_MAIN_LOOP ( myTime, myIter, myThid )                     跑   MDTHE_MAIN_LOOP( myTime, myIter, myThid )
         CALL TIMER_STOP ('MDTHE_MAIN_LOOP            [MD RUN]', myThid)
         CALL COST_FINAL_STORE ( myThid, lastdiva )          把最后的cost储存
      ELSE
#   ifdef ALLOW_DEBUG
         IF (debugMode) CALL DEBUG_CALL('ADTHE_MAIN_LOOP',myThid)
#   endif
         CALL TIMER_START('ADTHE_MAIN_LOOP       [ADJOINT RUN]', myThid)
         CALL ADTHE_MAIN_LOOP (  myThid )                                     跑     ADTHE_MAIN_LOOP(  myThid )  
         CALL TIMER_STOP ('ADTHE_MAIN_LOOP       [ADJOINT RUN]', myThid)
         CALL COST_FINAL_RESTORE ( myThid, lastdiva )
      ENDIF

#  else /* ALLOW_DIVIDED_ADJOINT undef */       -如果没有定义ALLOW_DIVIDED_ADJOINT   就跑ADTHE_MAIN_LOOP
#   ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_CALL('ADTHE_MAIN_LOOP',myThid)
#   endif
      CALL TIMER_START('ADTHE_MAIN_LOOP          [ADJOINT RUN]', myThid)
      CALL ADTHE_MAIN_LOOP ( myThid )
      CALL TIMER_STOP ('ADTHE_MAIN_LOOP          [ADJOINT RUN]', myThid)
#  endif /* ALLOW_DIVIDED_ADJOINT */

# else /* forward run only within AD setting */    这里是虽然定义了AUTODIFF但是没有定义别的什么 就只在AD设置里往前运行 就跑THE_MAIN_LOOP

#  ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_CALL('THE_MAIN_LOOP',myThid)
#  endif
C--   Call time stepping loop of full model                              调用全模型的时间步进回路
      CALL TIMER_START('THE_MAIN_LOOP          [THE_MODEL_MAIN]',myThid)
      CALL THE_MAIN_LOOP( myTime, myIter, myThid )                             THE_MAIN_LOOP( myTime, myIter, myThid )  
      CALL TIMER_STOP ('THE_MAIN_LOOP          [THE_MODEL_MAIN]',myThid)

# endif /* forward run only within AD setting */

# ifndef EXCLUDE_CTRL_PACK          大循环已经结束了，如果最后没有定义排除CTRL PACK ，那么就CTRL PACK
      IF ( useCTRL ) THEN
      IF ( lastdiva .AND. doMainPack ) THEN
         CALL TIMER_START('CTRL_PACK           [THE_MODEL_MAIN]',myThid)
         CALL CTRL_PACK( .FALSE. , myThid )
         CALL TIMER_STOP ('CTRL_PACK           [THE_MODEL_MAIN]',myThid)
         IF ( ( optimcycle.EQ.0 .OR. (.NOT. doMainUnpack) )
     &        .AND. myIter.EQ.nIter0 ) THEN
            CALL TIMER_START('CTRL_PACK     [THE_MODEL_MAIN]',myThid)
            CALL CTRL_PACK( .TRUE. , myThid )
            CALL TIMER_STOP ('CTRL_PACK     [THE_MODEL_MAIN]',myThid)
         ENDIF
      ENDIF
      ENDIF
# endif /* EXCLUDE_CTRL_PACK */

# ifdef ALLOW_GRDCHK               如果定义了梯度检查，那么就梯度检查
      IF ( useGrdchk .AND. lastdiva ) THEN
         CALL TIMER_START('GRDCHK_MAIN         [THE_MODEL_MAIN]',myThid)
         CALL GRDCHK_MAIN( myThid )
         CALL TIMER_STOP ('GRDCHK_MAIN         [THE_MODEL_MAIN]',myThid)
      ENDIF
# endif

#else /* ALL AD-related undef */     没有任何AD相关   如果没有定义ALLOW_ADMTLM（最外层）  就直接跑THE_MAIN_LOOP  直接是简单向前模型

# ifdef ALLOW_DEBUG             如果定义了debug，就检查THE_MAIN_LOOP    DEBUG_CALL
      IF (debugMode) CALL DEBUG_CALL('THE_MAIN_LOOP',myThid)
# endif
C--   Call time stepping loop of full model
      CALL TIMER_START('THE_MAIN_LOOP          [THE_MODEL_MAIN]',myThid)
      CALL THE_MAIN_LOOP( myTime, myIter, myThid )
      CALL TIMER_STOP ('THE_MAIN_LOOP          [THE_MODEL_MAIN]',myThid)

#endif /* ALLOW_TANGENTLINEAR_RUN ALLOW_ADJOINT_RUN ALLOW_ADMTLM */ 全部结束 ALLOW_TANGENTLINEAR_RUN ALLOW_ADJOINT_RUN ALLOW_ADMTLM

#ifdef ALLOW_STREAMICE        如果定义了STREAMICE 且use了，就调用它
      IF (useStreamIce) THEN
        CALL STREAMICE_FINALIZE_PETSC
      ENDIF
#endif

#ifdef ALLOW_MNC              如果定义了MNC，且use了，就调用它
      IF (useMNC) THEN
C       Close all open NetCDF files
        _BEGIN_MASTER( myThid )
        CALL MNC_FILE_CLOSE_ALL( myThid )
        _END_MASTER( myThid )
      ENDIF
#endif

C--   This timer encompasses the whole code      这个计时器包含了整个代码----这里结束的，在前面那个是开始的 
      CALL TIMER_STOP ('ALL                    [THE_MODEL_MAIN]',myThid)

C--   Write timer statistics            写入计时器统计信息      
      IF ( myThid .EQ. 1 ) THEN
       CALL TIMER_PRINTALL( myThid )    计时器_打印所有（myThid）
       CALL COMM_STATS                  通信统计
      ENDIF

C--   Check threads synchronization :   检查线程同步
      CALL BAR_CHECK( 9, myThid )

#ifdef ALLOW_DEBUG                      如果定义了debug，就最后DEBUG_LEAVE它，结束 
      IF (debugMode) CALL DEBUG_LEAVE('THE_MODEL_MAIN',myThid)
#endif

      RETURN
      END
