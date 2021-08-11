C $Header: /u/gcmpack/MITgcm/pkg/exf/EXF_OPTIONS.h,v 1.41 2017/04/14 23:14:48 jmc Exp $
C $Name:  $

CBOP
C !ROUTINE: EXF_OPTIONS.h
C !INTERFACE:
C #include "EXF_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for EXternal Forcing (EXF) package:
C | Control which optional features to compile in this package code.
控制要在此包代码中编译的可选功能
C *==================================================================*
CEOP

#ifndef EXF_OPTIONS_H
#define EXF_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_EXF
#ifdef ECCO_CPPOPTIONS_H

C-- When multi-package option-file ECCO_CPPOPTIONS.h is used (directly included
C    in CPP_OPTIONS.h), this option file is left empty since all options that
C   are specific to this package are assumed to be set in ECCO_CPPOPTIONS.h

#else /* ndef ECCO_CPPOPTIONS_H */

C-- Package-specific Options & Macros go here

C   --------------------
C   pkg/exf CPP options:
C   (see also table below on how to combine options)

C   > ( EXF_VERBOSE ) < replaced with run-time integer parameter "exf_debugLev"
C
C   >>> ALLOW_ATM_WIND <<<
C       If defined, set default value of run-time param. "useAtmWind" to True.
C       If useAtmWind=True, read-in and use wind vector (uwind/vwind)
C       to compute surface wind stress.
        如果定义了这个，设置一下run-time的默认值，useAtmWind就是true了，
        如果这个是true，那么就读入并且使用风矢量uwind和vwind来计算表面风应力
C
C   >>> ALLOW_ATM_TEMP <<<
C       This is the main EXF option controlling air-sea buoyancy fluxes:
C      If undefined, net heat flux (Qnet) and net fresh water flux
C       (EmP or EmPmR) are set according to hfluxfile & sfluxfile setting.
C      If defined, net heat flux and net fresh water flux are computed
C       from sum of various components (radiative SW,LW + turbulent heat
C       fluxes SH,LH ; Evap, Precip and optionally RunOff) thus ignoring
C       hfluxfile & sfluxfile.
C      In addition, it allows to read-in from files atmospheric temperature
C       and specific humidity, net radiative fluxes, and precip.
C       Also enable to read-in Evap (if EXF_READ_EVAP is defined) or
C       turbulent heat fluxes (if ALLOW_READ_TURBFLUXES is defined).
C
C   >>> ALLOW_DOWNWARD_RADIATION <<<
C       If defined, downward long-wave and short-wave radiation
C       can be read-in form files to compute net lwflux and swflux.
C
C   >>> ALLOW_ZENITHANGLE <<<
C       If defined, ocean albedo反照率 varies with the 天顶角zenith angle, and
C       incoming fluxes at the top of the atmosphere are computed
C
C   >>> ALLOW_BULKFORMULAE <<<
C       Allows the use of bulk formulae in order to estimate
C       turbulent fluxes (Sensible,Latent,Evap) at the ocean surface.
C       允许使用体积公式来估算海表面的湍流通量（感热、潜热、蒸发）
C   >>> EXF_CALC_ATMRHO
C       Calculate the local air density as function of temp, humidity
C       and pressure
C       根据温度、湿度和压力计算局部空气密度
C   >>> EXF_READ_EVAP <<<
C       If defined, evaporation field is read-in from file;
C     Note: if ALLOW_BULKFORMULAE is defined, evap that is computed from
C       atmospheric state will be replaced by read-in evap but computed
C       latent heat flux will be kept.
C       如果定义了，会从文件中读入蒸发，但是若又定义了体估算，从大气状态计算的蒸发会被替代
        计算的潜热通量却将被保留
C   >>> ALLOW_READ_TURBFLUXES <<<
C       If defined, turbulent heat fluxes (sensible and latent) can be read-in
C       from files (but overwritten if ALLOW_BULKFORMULAE is defined).
C       如果定义了，感热和潜热会从文件中读入，但是如果定义了体计算，还是会被覆盖掉
C   >>> ALLOW_RUNOFF <<<
C       If defined, river and glacier runoff can be read-in from files.
C       如果已定义，则可以从文件中读入河流和冰川径流。
C   >>> ALLOW_SALTFLX <<<
C       If defined, upward salt flux can be read-in from files.
C       如果定义了，可以从文件中读入向上的盐通量。
C   >>> ALLOW_RUNOFTEMP <<<
C       If defined, river and glacier runoff temperature
C       can be read-in from files.
C
C   >>> ATMOSPHERIC_LOADING <<<
C       If defined, atmospheric pressure can be read-in from files.
C   WARNING: this flag is set (define/undef) in CPP_OPTIONS.h
C            and cannot be changed here (in EXF_OPTIONS.h)
C       如果已定义，则可以从文件中读入大气压力。
        这个flag在CPP_OPTIONS.h中定义，并且不能在这里改变其值
C   >>> EXF_SEAICE_FRACTION <<<
C       If defined, seaice fraction can be read-in from files (areaMaskFile)
C       如果已定义，则可以从文件（areaMaskFile）读入seaice部分
C   >>> ALLOW_CLIMSST_RELAXATION <<<
C       Allow the relaxation to a monthly climatology of sea surface
C       temperature, e.g. the Reynolds climatology.
C       允许放松到每月的海面温度气候学，例如雷诺气候学。
C   >>> ALLOW_CLIMSSS_RELAXATION <<<
C       Allow the relaxation to a monthly climatology of sea surface
C       salinity, e.g. the Levitus climatology.
C
C   >>> USE_EXF_INTERPOLATION <<<
C       Allows to provide input field on arbitrary Lat-Lon input grid
C       (as specified in EXF_NML_04) and to interpolate to model grid.
C     Note: default is to interpolate unless {FLD}_interpMethod is set to 0
C       允许在任意横向输入网格上提供输入字段并插值到模型网格
C   ====================================================================
C
C    The following CPP options:
C       ALLOW_ATM_WIND / useAtmWind (useWind)
C       ALLOW_ATM_TEMP               (TEMP)
C       ALLOW_DOWNWARD_RADIATION     (DOWN)
C       ALLOW_BULKFORMULAE           (BULK)
C       EXF_READ_EVAP                (EVAP)
C       ALLOW_READ_TURBFLUXES        (TURB)
C
C    permit all ocean-model forcing configurations listed in the 2 tables below.
C    The first configuration (A1,B1) is the flux-forced, ocean model.
C    Configurations A2,B3 and A2,B4 use pkg/exf open-water bulk formulae
C    to compute, from atmospheric variables, the missing surface fluxes.
C    The forcing fields in the rightmost column are defined in EXF_FIELDS.h
C    (ocean-model surface forcing field are defined in model/inc/FFIELDS.h)
C允许下表中列出的所有海洋模型强迫配置。第一种结构（A1，B1）是通量强迫的海洋模型。配置A2、B3和A2、B4使用pkg/exf开放水域体积公式，从大气变量计算缺失的地表通量。最右边一栏中的强迫场在EXF\u fields.h中定义（海洋模型表面强迫场在model/inc/FFIELDS.h中定义）
C    (A) Surface momentum flux: [model: fu,fv ; exf: ustress,vstress]
C        表面动量通量             模型中是fu和fv，exf中是ustress和vstress
C    # |useWind|        actions
C   ---|-------|-------------------------------------------------------------
C   (1)| False | Read-in ustress,vstress (if needed in B, compute wind-speed)
C      |       | 读入ustress，vstress（如果B中需要，计算风速）
C   (2)| True  | Read-in uwind,vwind ; compute wind stress ustress,vstress.
                 读入uwind和vwind，用于计算风应力
C   ---|-------|-------------------------------------------------------------
C
C    (B) Surface buoyancy flux:
C        [ net heat flux: Qnet (exf: hflux), net short-wave: Qsw (exf: swflux)
C          fresh-water flux: EmPmR (exf: sflux) and saltFlux (exf: saltflx) ]
C          净热通量    Qnet        hflux      净短波     Qsw        swflux
           淡水通量    EmPmR       sflux      盐通量     saltFlux   saltflx
C    # |TEMP |DOWN |BULK |EVAP |TURB |            actions
C   ---|-----|-----|-----|-----|-----|-------------------------------------
C   (1)|  -  |  -  |  -  |  -  |  -  | Read-in hflux, swflux and sflux.
C      |     |     |     |     |     |
C   (2)|  -  | def |  -  |  -  |  -  | Read-in hflux, swdown and sflux.
C      |     |     |     |     |     | Compute swflux.
C      |     |     |     |     |     |
C   (3)| def | def | def |  -  |  -  | Read-in atemp, aqh, swdown, lwdown,
C      |     |     |     |     |     |  precip, and runoff.
C      |     |     |     |     |     | Compute hflux, swflux and sflux.
C      |     |     |     |     |     |
C   (4)| def |  -  | def |  -  |  -  | Read-in atemp, aqh, swflux, lwflux,
C      |     |     |     |     |     |  precip, and runoff.
C      |     |     |     |     |     | Compute hflux and sflux.
C      |     |     |     |     |     |
C   (5)| def | def |  -  | def | def | Read-in hs, hl, swdown, lwdown,
C      |     |     |     |     |     |  evap, precip and runoff.
C      |     |     |     |     |     | Compute hflux, swflux and sflux.
C      |     |     |     |     |     |
C   (6)| def |  -  |  -  | def | def | Read-in hs, hl, swflux, lwflux,
C      |     |     |     |     |     |  evap, precip and runoff.
C      |     |     |     |     |     | Compute  hflux and sflux.
C
C   =======================================================================

C-  Bulk formulae related flags.
#define ALLOW_ATM_TEMP
#define ALLOW_ATM_WIND
#define ALLOW_DOWNWARD_RADIATION
#ifdef ALLOW_ATM_TEMP
C Note: To use ALLOW_BULKFORMULAE or EXF_READ_EVAP, needs #define ALLOW_ATM_TEMP
# define ALLOW_BULKFORMULAE
# undef  ALLOW_BULK_LARGEYEAGER04
# undef  EXF_READ_EVAP
# ifndef ALLOW_BULKFORMULAE
C  Note: To use ALLOW_READ_TURBFLUXES, ALLOW_ATM_TEMP needs to
C        be defined but ALLOW_BULKFORMULAE needs to be undef
#  define ALLOW_READ_TURBFLUXES
# endif
#endif /* ALLOW_ATM_TEMP */

C-  Other forcing fields
#define ALLOW_RUNOFF
#undef  ALLOW_RUNOFTEMP
#define ALLOW_SALTFLX

#if (defined (ALLOW_BULKFORMULAE) && defined (ATMOSPHERIC_LOADING))
C Note: To use EXF_CALC_ATMRHO, both ALLOW_BULKFORMULAE
C       and ATMOSPHERIC_LOADING need to be defined
# undef EXF_CALC_ATMRHO
#endif

C-  Zenith Angle/Albedo related flags.
#ifdef ALLOW_DOWNWARD_RADIATION
# undef ALLOW_ZENITHANGLE
#endif

C-  Use 海洋发射率ocean_emissivity*lwdown in lwFlux. This flag should be defined
C   unless to reproduce old results (obtained with inconsistent old code)
#ifdef ALLOW_DOWNWARD_RADIATION
# define EXF_LWDOWN_WITH_EMISSIVITY
#endif

C-  Relaxation to monthly climatologies.
#define ALLOW_CLIMSST_RELAXATION
#define ALLOW_CLIMSSS_RELAXATION

C-  Allows to read-in seaice fraction from files (areaMaskFile)
#undef EXF_SEAICE_FRACTION

C-  Use spatial interpolation to interpolate
C   forcing files from input grid to model grid.
#undef USE_EXF_INTERPOLATION
C   for interpolated vector fields, rotate towards model-grid axis
C   using old rotation formulae (instead of grid-angles)
#undef EXF_USE_OLD_VEC_ROTATION
C   for interpolation around N & S pole, use the old formulation
C   (no pole symmetry, single vector-comp interp, reset to 0 zonal-comp @ N.pole)
#undef EXF_USE_OLD_INTERP_POLE

#define EXF_INTERP_USE_DYNALLOC
#if ( defined USE_EXF_INTERPOLATION && defined EXF_INTERP_USE_DYNALLOC && defined USING_THREADS )
# define EXF_IREAD_USE_GLOBAL_POINTER
#endif

#endif /* ndef ECCO_CPPOPTIONS_H */
#endif /* ALLOW_EXF */
#endif /* EXF_OPTIONS_H */
