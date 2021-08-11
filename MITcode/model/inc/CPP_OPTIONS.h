C $Header: /u/gcmpack/MITgcm/model/inc/CPP_OPTIONS.h,v 1.54 2016/11/28 22:55:00 jmc Exp $
C $Name:  $

#ifndef CPP_OPTIONS_H
#define CPP_OPTIONS_H

CBOP
C !ROUTINE: CPP_OPTIONS.h
C !INTERFACE:
C #include "CPP_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | main CPP options file for the model:         模型的主CPP选项文件：控制在model/src代码中编译哪些可选特性。
C | Control which optional features to compile in model/src code.
C *==================================================================*
CEOP

C CPP flags controlling particular source code features   控制特定源代码功能的CPP标志

C o Shortwave heating as extra term in external_forcing.F   短波加热作为外部强迫的附加项
C Note: this should be a run-time option   这应该是一个运行时选项
#undef SHORTWAVE_HEATING

C o Include/exclude Geothermal Heat Flux at the bottom of the ocean  包括/不包括海底地热通量
#undef ALLOW_GEOTHERMAL_FLUX

C o Include/exclude phi_hyd calculation code  包括/排除功率因数液压计算代码
#define INCLUDE_PHIHYD_CALCULATION_CODE

C o Include/exclude call to S/R CONVECT  包括/排除对S/R Converct的调用
#define INCLUDE_CONVECT_CALL

C o Include/exclude call to S/R CALC_DIFFUSIVITY  包括/排除对S/R计算扩散率的调用
#define INCLUDE_CALC_DIFFUSIVITY_CALL

C o Allow full 3D specification of vertical diffusivity  允许垂直扩散率的全三维规格
#undef ALLOW_3D_DIFFKR

C o Allow latitudinally varying BryanLewis79 vertical diffusivity  允许沿纬度变化的BryanLewis79垂直扩散率
#undef ALLOW_BL79_LAT_VARY

C o Include/exclude Implicit vertical advection code   包含/排除隐式垂直平流代码
#define INCLUDE_IMPLVERTADV_CODE

C o Include/exclude combined Surf.Pressure and Drag Implicit solver code  包含/排除组合冲浪。压力和拖动隐式求解器代码
#undef ALLOW_SOLVE4_PS_AND_DRAG

C o Include/exclude AdamsBashforth-3rd-Order code  包含/排除AdamsBashforth-3阶代码
#undef ALLOW_ADAMSBASHFORTH_3

C o Include/exclude nonHydrostatic code  包括/排除非静水压代码
#undef ALLOW_NONHYDROSTATIC

C o Allow to account for heating due to friction (and momentum dissipation)  考虑摩擦（和动量耗散）产生的热量
#undef ALLOW_FRICTION_HEATING

C o Allow mass source or sink of Fluid in the interior  允许内部流体的质量源或质量汇
C   (3-D generalisation of oceanic real-fresh water flux)  海洋真实淡水通量的三维概化
#undef ALLOW_ADDFLUID

C o Include pressure loading code  包括压力加载代码
#define ATMOSPHERIC_LOADING

C o exclude/allow external forcing-fields load  排除/允许外部强制场加载
C   this allows to read & do simple linear time interpolation of oceanic 
C   forcing fields, if no specific pkg (e.g., EXF) is used to compute them.
如果不使用特定的pkg（例如EXF）来计算海洋强迫场，则可以读取并进行简单的线性时间插值。
#undef EXCLUDE_FFIELDS_LOAD

C o Include/exclude balancing surface forcing fluxes code  包括/排除平衡表面强制通量代码
#undef ALLOW_BALANCE_FLUXES

C o Include/exclude balancing surface forcing relaxation code  包括/排除平衡面强制松弛代码
#undef ALLOW_BALANCE_RELAX

C o Include/exclude GM-like eddy stress in momentum code  在动量代码中包括/排除类GM涡流应力
#undef ALLOW_EDDYPSI

C o Use "Exact Convervation" of fluid in Free-Surface formulation    在自由面公式中使用流体的“精确对流”，使d/dt（eta）完全等于-Div.输运
C   so that d/dt(eta) is exactly equal to - Div.Transport
#define EXACT_CONSERV

C o Allow the use of Non-Linear Free-Surface formulation
C   this implies that surface thickness (hFactors) vary with time  允许使用非线性自由表面公式这意味着表面厚度（hFactors）随时间而变化
#undef NONLIN_FRSURF

C o Include/exclude code for single reduction Conjugate-Gradient solver  单归约共轭梯度求解器的包含/排除代码
#define ALLOW_SRCG

C o Choices for implicit solver routines solve_*diagonal.F  隐式解算器例程求解*diagonal.F的选择
C   The following has low memory footprint, but not suitable for AD  以下具有较低的内存占用，但不适用于AD
#undef SOLVE_DIAGONAL_LOWMEMORY
C   The following one suitable for AD but does not vectorize  以下一个适用于AD，但不向量化
#undef SOLVE_DIAGONAL_KINNER

C o ALLOW isotropic scaling of harmonic and bi-harmonic terms when
C   using an locally isotropic spherical grid with (dlambda) x (dphi*cos(phi))
C *only for use on a lat-lon grid*
C   Setting this flag here affects both momentum and tracer equation unless
C   it is set/unset again in other header fields (e.g., GAD_OPTIONS.h).
C   The definition of the flag is commented to avoid interference with
C   such other header files.
C   The preferred method is specifying a value for viscAhGrid or viscA4Grid
C   in data which is then automatically scaled by the grid size;
C   the old method of specifying viscAh/viscA4 and this flag is provided
C   for completeness only (and for use with the adjoint).
当使用具有（dlambda）x（dphi*cos（phi））*的局部各向同性球形网格时，允许谐波和双谐波项的各向同性缩放*仅用于lat-lon网格*
在此处设置此标志会影响动量和示踪方程，除非在其他标题字段（例如，GAD\u OPTIONS.h）中再次设置/取消设置。
标记的定义被注释以避免与其他头文件的干扰。首选方法是在数据中为viscAhGrid或viscA4Grid指定一个值，然后根据网格大小自动缩放该值；
指定viscAh/viscA4和该标志的旧方法仅用于完整性（并用于伴随）。
C#define ISOTROPIC_COS_SCALING

C o This flag selects the form of COSINE(lat) scaling of bi-harmonic term.
C *only for use on a lat-lon grid*
C   Has no effect if ISOTROPIC_COS_SCALING is undefined.
C   Has no effect on vector invariant momentum equations.
C   Setting this flag here affects both momentum and tracer equation unless
C   it is set/unset again in other header fields (e.g., GAD_OPTIONS.h).
C   The definition of the flag is commented to avoid interference with
C   such other header files.
此标志选择双调和项的余弦（lat）缩放形式*如果各向同性COS缩放未定义，则仅用于lat-lon网格*没有效果。
对矢量不变动量方程没有影响。在这里设置这个标志会影响动量和示踪方程，
除非它在其他标题字段（例如GAD\u OPTIONS.h）中再次设置/取消设置。标记的定义被注释以避免与其他头文件的干扰。
C#define COSINEMETH_III

C o Use "OLD" UV discretisation near boundaries (*not* recommended)  在边界附近使用“旧”UV离散化（*不推荐）
C   Note - only works with pkg/mom_fluxform and "no_slip_sides=.FALSE."
C          because the old code did not have no-slip BCs
#undef OLD_ADV_BCS

C o Use LONG.bin, LATG.bin, etc., initialization for ini_curviliear_grid.F  使用LONG.bin、LATG.bin等初始化ini\u curvilier\u grid.F
C   Default is to use "new" grid files (OLD_GRID_IO undef) but OLD_GRID_IO
C   is still useful with, e.g., single-domain curvilinear configurations.
#undef OLD_GRID_IO

C o Use old EXTERNAL_FORCING_U,V,T,S subroutines (for backward compatibility)  使用旧的外部子例程（为了向后兼容）
#undef USE_OLD_EXTERNAL_FORCING

C o Execution environment support options     执行环境支持选项
#include "CPP_EEOPTIONS.h"

C o Include/exclude single header file containing multiple packages options
C   (AUTODIFF, COST, CTRL, ECCO, EXF ...) instead of the standard way where
C   each of the above pkg get its own options from its specific option file.
C   Although this method, inherited from ECCO setup, has been traditionally
C   used for all adjoint built, work is in progress to allow to use the
C   standard method also for adjoint built.
包含/排除包含多个包选项（AUTODIFF、COST、CTRL、ECCO、EXF…）的单个头文件，而不是上述每个包从其特定选项文件中获取自己的选项的标准方式。
虽然这种方法继承自ECCO设置，传统上用于所有伴随生成，但允许对伴随生成也使用标准方法的工作仍在进行中。
c#ifdef PACKAGES_CONFIG_H
c# include "ECCO_CPPOPTIONS.h"
c#endif

#endif /* CPP_OPTIONS_H */
