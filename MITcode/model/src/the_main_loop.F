C $Header: /u/gcmpack/MITgcm/model/src/the_main_loop.F,v 1.181 2017/03/20 19:38:16 mmazloff Exp $
C $Name:  $

#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#endif
#ifdef ALLOW_SEAICE
# include "SEAICE_OPTIONS.h"
#endif
#ifdef ALLOW_GENERIC_ADVDIFF
# include "GAD_OPTIONS.h"
#endif
#ifdef ALLOW_GMREDI
# include "GMREDI_OPTIONS.h"
#endif
#ifdef ALLOW_STREAMICE
# include "STREAMICE_OPTIONS.h"
#endif
#ifdef ALLOW_GGL90
# include "GGL90_OPTIONS.h"
#endif
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
#endif
#ifdef ALLOW_COST
# include "COST_OPTIONS.h"
#endif
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_ECCO
# include "ECCO_OPTIONS.h"
#endif
#ifdef ALLOW_PROFILES
# include "PROFILES_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: THE_MAIN_LOOP
C     !INTERFACE:
      SUBROUTINE THE_MAIN_LOOP( myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *================================================================*
C     | SUBROUTINE the_main_loop
C     | o Run the ocean model and evaluate the specified cost function.
C     *================================================================*
C     |
C     | THE_MAIN_LOOP is the toplevel routine for the Tangent Linear and    the_main_loop是TAMC的顶级例程
C     | Adjoint Model Compiler (TAMC).
C     | For this purpose the initialization       为了这个目的，模型的初始化被分成了俩个部分
C     | of the model was split into two parts. Those parameters that do     不依赖于具体的某个模型运行的那些参数在INITIALISE FIXED里初始化
C     | not depend on a specific model run are set in INITIALISE_FIXED,     
C     | whereas those that do depend on the specific realization are        依赖于某个具体实现的在INITIALISE VARIAL里初始化
C     | initialized in INITIALISE_VARIA.
C     | This routine is to be used in conjuction with the MITgcmuv
C     | checkpoint 37.
C     *================================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables ==    全局变量
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

c**************************************
#ifdef ALLOW_AUTODIFF      如果定义了AUTODIFF 
# ifndef ALLOW_OPENAD  若定义了openad

c These includes are needed for  这些头文件是ad-checkpoing需要的，提供场用于储存
c AD-checkpointing.
c They provide the fields to be stored.

#  include "AUTODIFF_MYFIELDS.h"
#  include "GRID.h"
#  include "DYNVARS.h"
#  include "SURFACE.h"
#  include "FFIELDS.h"
#  include "EOS.h"
#  include "AUTODIFF.h"

#  ifdef ALLOW_GENERIC_ADVDIFF
#   include "GAD.h"
#   include "GAD_SOM_VARS.h"
#  endif
#  ifdef ALLOW_MOM_FLUXFORM
#   include "MOM_FLUXFORM.h"
#  endif
#  ifdef ALLOW_CD_CODE
#   include "CD_CODE_VARS.h"
#  endif
#  ifdef ALLOW_PTRACERS
#   include "PTRACERS_SIZE.h"
#   include "PTRACERS_FIELDS.h"
#   include "PTRACERS_START.h"
#  endif
#  ifdef ALLOW_GCHEM
#   include "GCHEM_FIELDS.h"
#  endif
#  ifdef ALLOW_CFC
#   include "CFC.h"
#  endif
#  ifdef ALLOW_DIC
#   include "DIC_VARS.h"
#   include "DIC_LOAD.h"
#   include "DIC_ATMOS.h"
#   include "DIC_CTRL.h"
#   include "DIC_COST.h"
#  endif
#  ifdef ALLOW_BLING
#   include "BLING_VARS.h"
#   include "BLING_LOAD.h"
#  endif
#  ifdef ALLOW_OBCS
#   include "OBCS_PARAMS.h"
#   include "OBCS_FIELDS.h"
#   include "OBCS_SEAICE.h"
#   ifdef ALLOW_PTRACERS
#    include "OBCS_PTRACERS.h"
#   endif
#  endif
#  ifdef ALLOW_EXF
#   include "EXF_FIELDS.h"
#   ifdef ALLOW_BULKFORMULAE
#    include "EXF_CONSTANTS.h"
#   endif
#  endif /* ALLOW_EXF */
#  ifdef ALLOW_SEAICE
#   include "SEAICE_SIZE.h"
#   include "SEAICE.h"
#   include "SEAICE_PARAMS.h"
#   include "SEAICE_COST.h"
#   include "SEAICE_TRACER.h"
#  endif
#  ifdef ALLOW_SALT_PLUME
#   include "SALT_PLUME.h"
#  endif
#  ifdef ALLOW_THSICE
#   include "THSICE_SIZE.h"
#   include "THSICE_VARS.h"
#   include "THSICE_COST.h"
#  endif
#  ifdef ALLOW_SHELFICE
#   include "SHELFICE.h"
#   include "SHELFICE_COST.h"
#  endif
#  ifdef ALLOW_STREAMICE
#   include "STREAMICE.h"
#   include "STREAMICE_ADV.h"
#   include "STREAMICE_BDRY.h"
#   include "STREAMICE_CG.h"
#  endif
#  ifdef ALLOW_EBM
#   include "EBM.h"
#  endif
#  ifdef ALLOW_RBCS
#   include "RBCS_SIZE.h"
#   include "RBCS_FIELDS.h"
#  endif
#  ifdef ALLOW_OFFLINE
#   include "OFFLINE.h"
#  endif
#  ifdef ALLOW_CG2D_NSA
#   include "CG2D.h"
#  endif
#  ifdef ALLOW_DIVIDED_ADJOINT
#   ifdef ALLOW_USE_MPI
#    include "mpif.h"
#   endif
#  endif

#  include "tamc.h"

# ifdef ALLOW_GGL90
#  include "GGL90.h"
# endif
# ifdef ALLOW_PROFILES
#  include "PROFILES_SIZE.h"
#  include "profiles.h"
# endif
# ifdef ALLOW_ECCO
#  include "ecco_cost.h"
# endif

# endif /* undef ALLOW_OPENAD */       结束即没有定义openad

# ifdef ALLOW_CTRL                     在AUTODIFF中，若定义了CTRL，就加入对应头文件
# include "CTRL_SIZE.h"
# include "ctrl.h"
# include "ctrl_dummy.h"
# include "CTRL_GENARR.h"
# include "CTRL_OBCS.h"
# endif
# ifdef ALLOW_COST                     若定义了COST，也加入cost.h
#  include "cost.h"
# endif

#endif /* ALLOW_AUTODIFF */            结束if AUTODIFF
c**************************************

C     !INPUT/OUTPUT PARAMETERS:      输入输出参数
C     == Routine arguments ==
C     note: under the multi-threaded model myIter and                 又是这段话
C           myTime are local variables passed around as routine
C           arguments. Although this is fiddly it saves the need to
C           impose additional synchronisation points when they are
C           updated.
在多线程模型下，myIter和myTime是作为例程参数传递的局部变量。虽然这是微小难弄的的，但当它们被更新时节省了需要附加的同步点。
C     myTime :: time counter for this thread
C     myIter :: iteration counter for this thread
C     myThid :: thread number for this instance of the routine.
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !FUNCTIONS:
C     == Functions ==

C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER iloop
#ifdef STORE_LOADEDREC_TEST
      INTEGER bi,bj
#endif /* STORE_LOADEDREC_TEST */   存储加载测试
CEOP

#ifdef ALLOW_DEBUG          debug      DEBUG_ENTER
      IF (debugMode) CALL DEBUG_ENTER('THE_MAIN_LOOP',myThid)
#endif

#ifdef ALLOW_AUTODIFF_TAMC
c--   Initialize storage for the cost function evaluation.    为成本函数评估初始化存储
CADJ  INIT dummytape = common, 1    这是TAF的命令，没有lisense用不了
c--   Initialize storage for the outermost loop.  初始化最外层循环的存储
# ifndef ALLOW_AUTODIFF_WHTAPEIO
CADJ  INIT tapelev_init   = USER
CADJ  INIT tapelev_ini_bibj_k   = USER
# else
CADJ  INIT tapelev_init   = common, 1
# endif
c
# ifdef ALLOW_TAMC_CHECKPOINTING
#  if (defined (AUTODIFF_2_LEVEL_CHECKPOINT))
CADJ  INIT tapelev2 = USER
#  elif (defined (AUTODIFF_4_LEVEL_CHECKPOINT))
CADJ  INIT tapelev4 = USER
#  else
CADJ  INIT tapelev3 = USER
#  endif
# endif
#endif

#ifdef ALLOW_AUTODIFF
      nIter0 = NINT( (startTime-baseTime)/deltaTClock )
      ikey_dynamics = 1
#endif

#ifdef ALLOW_AUTODIFF_TAMC
# ifdef NONLIN_FRSURF
CADJ STORE hFacC = tapelev_init, key = 1
CADJ STORE recip_hFacC = tapelev_init, key = 1
# endif
#endif

#ifdef ALLOW_OPENAD      如果定义了openad
# ifdef ALLOW_THETA0_CONTROL     如果温度需要控制
c$openad INDEPENDENT(xx_theta_dummy)
# endif
# ifdef ALLOW_SALT0_CONTROL
c$openad INDEPENDENT(xx_salt_dummy)
# endif
# ifdef ALLOW_HFLUX0_CONTROL
c$openad INDEPENDENT(xx_hflux0)
# endif
# ifdef ALLOW_SFLUX0_CONTROL
c$openad INDEPENDENT(xx_sflux0)
# endif
# ifdef ALLOW_TAUU0_CONTROL
c$openad INDEPENDENT(xx_tauu0)
# endif
# ifdef ALLOW_TAUV0_CONTROL
c$openad INDEPENDENT(xx_tauv0)
# endif
# ifdef ALLOW_DIFFKR_CONTROL
c$openad INDEPENDENT(xx_diffkr_dummy)
# endif
# ifdef ALLOW_KAPGM_CONTROL
c$openad INDEPENDENT(xx_kapgm)
# endif
# ifdef ALLOW_TR10_CONTROL
c$openad INDEPENDENT(xx_tr1)
# endif
# ifdef ALLOW_HFLUXM_CONTROL
c$openad INDEPENDENT(xx_hfluxm)
# endif
# ifdef ALLOW_GENARR2D_CONTROL
cphc$openad INDEPENDENT(xx_genarr2d)
c$openad INDEPENDENT(xx_genarr2d_dummy)
# endif
# ifdef ALLOW_GENARR3D_CONTROL
cphc$openad INDEPENDENT(xx_genarr3d)
c$openad INDEPENDENT(xx_genarr3d_dummy)
# endif
# ifdef ALLOW_GENTIM2D_CONTROL
cphc$openad INDEPENDENT(xx_gentim2d)
c$openad INDEPENDENT(xx_gentim2d_dummy)
# endif
#endif /* ALLOW_OPENAD */  结束openad

#ifdef ALLOW_DEBUG      DEBUG_CALL  INITIALISE_VARIA   也即在使用之前debug一下
      IF (debugMode) CALL DEBUG_CALL('INITIALISE_VARIA',myThid)
#endif
C--   Set initial conditions (variable arrays)      设置初始条件（变量数组）
      CALL TIMER_START('INITIALISE_VARIA    [THE_MAIN_LOOP]', myThid)
      CALL INITIALISE_VARIA( myThid )                             初始化这个例程的条件 INITIALISE_VARIA( myThid ) 
      CALL TIMER_STOP ('INITIALISE_VARIA    [THE_MAIN_LOOP]', myThid)
#ifdef ECCO_CTRL_DEPRECATED
#ifdef ALLOW_ECCO
cph: avoid renewed call of initialise_varia in recomputation loop
cph: in some circumstances
CADJ STORE sbar_gen,tbar_gen  = onetape
#endif
#endif

#ifdef ALLOW_SHOWFLOPS
      CALL TIMER_START('SHOWFLOPS_INIT      [THE_MAIN_LOOP]', myThid)
      CALL SHOWFLOPS_INIT( myThid )
      CALL TIMER_STOP ('SHOWFLOPS_INIT      [THE_MAIN_LOOP]', myThid)
#endif

c--   Do the model integration.              进行模型积分
      CALL TIMER_START('MAIN LOOP           [THE_MAIN_LOOP]', myThid)          计时开始

c     >>>>>>>>>>>>>>>>>>>>>>>>>>>   LOOP   <<<<<<<<<<<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>>>>>>>>>  STARTS  <<<<<<<<<<<<<<<<<<<<<<<<<<<<

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef ALLOW_OPENAD                     如果定义了openad，
# ifdef ALLOW_AUTODIFF                     又定义了autodiff
#  ifdef ALLOW_TAMC_CHECKPOINTING            还定义了checkpoint
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      max_lev4=nTimeSteps/(nchklev_1*nchklev_2*nchklev_3)+1
      max_lev3=nTimeSteps/(nchklev_1*nchklev_2)+1
      max_lev2=nTimeSteps/nchklev_1+1

c**************************************
#   ifdef ALLOW_DIVIDED_ADJOINT     DIVIDED_ADJOINT 
CADJ loop = divided
#   endif
c**************************************

#   ifdef AUTODIFF_4_LEVEL_CHECKPOINT
      do ilev_4 = 1,nchklev_4                   这些do循环的结束体都在后面，从
         if(ilev_4.le.max_lev4) then            这里的if判断也是
c**************************************
#ifdef ALLOW_AUTODIFF_WHTAPEIO     允许自动差异
            CALL AUTODIFF_WHTAPEIO_SYNC( 4 , 0, myThid )
#endif
            CALL AUTODIFF_STORE( myThid )
#include "checkpoint_lev4_directives.h"
            CALL AUTODIFF_RESTORE( myThid )
#ifdef ALLOW_AUTODIFF_WHTAPEIO
            CALL AUTODIFF_WHTAPEIO_SYNC( 4 , 1, myThid )
#endif
c**************************************
c--     Initialise storage for the middle loop.   初始化中间循环的存储
CADJ    INIT tapelev3 = USER
#   endif /* AUTODIFF_4_LEVEL_CHECKPOINT */   结束4层的判断

#   ifndef AUTODIFF_2_LEVEL_CHECKPOINT        开始2层的判断
      do ilev_3 = 1,nchklev_3                        do循环
         if(ilev_3.le.max_lev3) then                 if判断
c**************************************
#ifdef ALLOW_AUTODIFF_WHTAPEIO
            CALL AUTODIFF_WHTAPEIO_SYNC( 3 , 0, myThid )
#endif
            CALL AUTODIFF_STORE( myThid )
#include "checkpoint_lev3_directives.h"
            CALL AUTODIFF_RESTORE( myThid )
#ifdef ALLOW_AUTODIFF_WHTAPEIO
            CALL AUTODIFF_WHTAPEIO_SYNC( 3 , 1, myThid )
#endif
c**************************************
c--     Initialise storage for the middle loop.  初始化中间循环的存储
CADJ    INIT tapelev2 = USER
#   endif /* AUTODIFF_2_LEVEL_CHECKPOINT */     结束2层的判断

        do ilev_2 = 1,nchklev_2                 这应该是默认的3层检查点   do循环
         if(ilev_2.le.max_lev2) then                                    if判断
c**************************************
#ifdef ALLOW_AUTODIFF_WHTAPEIO
            CALL AUTODIFF_WHTAPEIO_SYNC( 2 , 0, myThid )   自动微分同步
#endif
            CALL AUTODIFF_STORE( myThid )                  自动微分存储
#include "checkpoint_lev2_directives.h"
            CALL AUTODIFF_RESTORE( myThid )                自动微分恢复
#ifdef ALLOW_AUTODIFF_WHTAPEIO
            CALL AUTODIFF_WHTAPEIO_SYNC( 2 , 1, myThid )   自动微分同步
#endif
c**************************************

# endif /* ALLOW_TAMC_CHECKPOINTING */           结束了ALLOW_AUTODIFF下的检查点层数判断，还是在OPENAD里面

c**************************************
c--
c--       Initialize storage for the innermost loop.               初始化最内部循环的存储
c--       Always check common block sizes for the checkpointing!   始终检查检查点的公共块大小        都结合了CADJ
c--
CADJ INIT comlev1        = COMMON,nchklev_1
CADJ INIT comlev1_bibj   = COMMON,nchklev_1*nsx*nsy*nthreads_chkpt
CADJ INIT comlev1_bibj_k = COMMON,nchklev_1*nsx*nsy*nr*nthreads_chkpt
c--
#ifdef ALLOW_CTRL
CADJ INIT ctrltape       = COMMON,1
#endif
c--
#   ifdef ALLOW_KPP
CADJ INIT comlev1_kpp    = COMMON,nchklev_1*nsx*nsy
CADJ INIT comlev1_kpp_k  = COMMON,nchklev_1*nsx*nsy*nr
#   endif /* ALLOW_KPP */
c--
#   ifdef ALLOW_GMREDI
CADJ INIT comlev1_gmredi_k_gad
CADJ &    = COMMON,nchklev_1*nsx*nsy*nr*nthreads_chkpt*maxpass
#   endif /* ALLOW_GMREDI */
c--
#   ifdef ALLOW_PTRACERS
CADJ INIT comlev1_bibj_ptracers = COMMON,
CADJ &    nchklev_1*nsx*nsy*nthreads_chkpt*PTRACERS_num
CADJ INIT comlev1_bibj_k_ptracers = COMMON,
CADJ &    nchklev_1*nsx*nsy*nthreads_chkpt*PTRACERS_num*nr
#   endif /* ALLOW_PTRACERS */
c--
#   ifndef DISABLE_MULTIDIM_ADVECTION
CADJ INIT comlev1_bibj_gad = COMMON,
CADJ &    nchklev_1*nsx*nsy*nthreads_chkpt*maxpass
CADJ INIT comlev1_bibj_k_gad = COMMON,
CADJ &    nchklev_1*nsx*nsy*nr*nthreads_chkpt*maxpass
CADJ INIT comlev1_bibj_k_gad_pass = COMMON,
CADJ &    nchklev_1*nsx*nsy*nr*nthreads_chkpt*maxpass*maxpass
#   endif /* DISABLE_MULTIDIM_ADVECTION */
c--
#   ifdef ALLOW_MOM_COMMON
#   ifndef AUTODIFF_DISABLE_LEITH
CADJ INIT comlev1_mom_ijk_loop
CADJ &     = COMMON,nchklev_1*
CADJ &       (snx+2*olx)*nsx*(sny+2*oly)*nsy*nr*nthreads_chkpt
#   endif /* AUTODIFF_DISABLE_LEITH */
#   endif /* ALLOW_MOM_COMMON */
c--
#   if (defined (ALLOW_EXF) && defined (ALLOW_BULKFORMULAE))
CADJ INIT comlev1_exf_1
CADJ &     = COMMON,nchklev_1*snx*nsx*sny*nsy*nthreads_chkpt
CADJ INIT comlev1_exf_2
CADJ &     = COMMON,niter_bulk*nchklev_1*snx*nsx*sny*nsy*nthreads_chkpt
#   endif /* ALLOW_BULKFORMULAE */
c--
#   ifdef ALLOW_SEAICE
#    ifdef SEAICE_ALLOW_DYNAMICS
CADJ INIT comlev1_dynsol = COMMON,nchklev_1*MPSEUDOTIMESTEPS
#     ifdef SEAICE_LSR_ADJOINT_ITER
CADJ INIT comlev1_lsr = COMMON,nchklev_1*MPSEUDOTIMESTEPS*SOLV_MAX_FIXED
CADJ INIT comlev1_bibj_lsr = COMMON,
CADJ &  nchklev_1*nsx*nsy*nthreads_chkpt*MPSEUDOTIMESTEPS*SOLV_MAX_FIXED
#     else /* make the common blocks smaller to reduce memory footprint */
CADJ INIT comlev1_lsr = COMMON,nchklev_1*MPSEUDOTIMESTEPS
CADJ INIT comlev1_bibj_lsr = COMMON,
CADJ &  nsx*nsy*nthreads_chkpt*nchklev_1*MPSEUDOTIMESTEPS
#     endif
#    endif
#    ifdef SEAICE_ALLOW_EVP
CADJ INIT comlev1_evp = COMMON,nEVPstepMax*nchklev_1
CADJ INIT comlev1_bibj_evp = COMMON,
CADJ &  nsx*nsy*nthreads_chkpt*nEVPstepMax*nchklev_1
#    endif
CML#   ifdef SEAICE_MULTICATEGORY
CMLCADJ INIT comlev1_multdim
CMLCADJ &    = COMMON,nchklev_1*nsx*nsy*nthreads_chkpt*nitd
CML#   endif
#    ifndef DISABLE_MULTIDIM_ADVECTION
CADJ INIT comlev1_bibj_k_gadice = COMMON,
CADJ &    nchklev_1*nsx*nsy*nthreads_chkpt*maxpass
CADJ INIT comlev1_bibj_k_gadice_pass = COMMON,
CADJ &    nchklev_1*nsx*nsy*nthreads_chkpt*maxpass*maxpass
#    endif /* DISABLE_MULTIDIM_ADVECTION */
# endif /* ALLOW_SEAICE */
c--
#   ifdef ALLOW_THSICE
CADJ INIT comlev1_thsice_1 = COMMON,
CADJ &    nchklev_1*snx*nsx*sny*nsy*nthreads_chkpt
CADJ INIT comlev1_thsice_2 = COMMON,
CADJ &    nchklev_1*snx*nsx*sny*nsy*nlyr*nthreads_chkpt
CADJ INIT comlev1_thsice_3 = COMMON,
CADJ &    nchklev_1*snx*nsx*sny*nsy*MaxTsf*nthreads_chkpt
CADJ INIT comlev1_thsice_4 = COMMON,
CADJ &    nchklev_1*nsx*nsy*maxpass*nthreads_chkpt
CADJ INIT comlev1_thsice_5 = COMMON,
CADJ &    nchklev_1*snx*nsx*sny*nsy*MaxTsf*(niter_bulk+1)*nthreads_chkpt
CADJ INIT comlev1_thsice_s4t = COMMON,
CADJ &    nchklev_1*nsx*nsy*maxtsf*nthreads_chkpt
#   endif /* ALLOW_THSICE */
c--
#   ifdef ALLOW_STREAMICE
CADJ INIT comlev1_stream_nl = COMMON,nchklev_1*streamice_max_nl
CADJ INIT comlev1_stream_front = COMMON,nchklev_1*4
CADJ INIT comlev1_stream_ij
CADJ &     = COMMON,nchklev_1*4*(snx+2)*nsx*(sny+2)*nsy
CADJ INIT comlev1_stream_hybrid
CADJ &     = COMMON,nchklev_1*snx*nsx*sny*nsy*nr*nthreads_chkpt
#   endif
c--
#   ifdef ALLOW_CG2D_NSA
CADJ INIT comlev1_cg2d
CADJ &     = COMMON,nchklev_1*nthreads_chkpt
CADJ INIT comlev1_cg2d_iter
CADJ &     = COMMON,nchklev_1*nthreads_chkpt*numItersMax
#   endif
c--
c**************************************

#ifdef STORE_LOADEDREC_TEST   存储加载测试
          DO bj = myByLo(myThid), myByHi(myThid)
           DO bi = myBxLo(myThid), myBxHi(myThid)
            loadedRec(bi,bj) = 0
           ENDDO
          ENDDO
#endif /* STORE_LOADEDREC_TEST */

#ifdef ALLOW_TAMC_CHECKPOINTING    如果定义了checkpoint，那么循环的步数需要一些调整        这里是离开了AUTODIFF判断，只在OPENAD里面

          do ilev_1 = 1,nchklev_1           这又一个do循环

c--         The if-statement below introduces a some flexibility in the    下面的if语句在选择3元组时引入了一些灵活性
c--         choice of the 3-tupel ( nchklev_1, nchklev_2, nchklev_3 ).

            iloop = (ilev_2 - 1)*nchklev_1 + ilev_1
#    ifndef AUTODIFF_2_LEVEL_CHECKPOINT
     &            + (ilev_3 - 1)*nchklev_2*nchklev_1
#    endif
#    ifdef AUTODIFF_4_LEVEL_CHECKPOINT
     &            + (ilev_4 - 1)*nchklev_3*nchklev_2*nchklev_1
#    endif

            if ( iloop .le. nTimeSteps ) then     又一个if判断

#  else /* ALLOW_TAMC_CHECKPOINTING  undefined */

      DO iloop = 1, nTimeSteps      这里跟下面的一个ENDDO配套，是没定义check的   这里跟下面的主循环iloop是同一个开始，没定义AUTODIFF没定义check，或者没定义openad定义了autodiff

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  endif /* ALLOW_TAMC_CHECKPOINTING */
# endif /* ALLOW_AUTODIFF */
#endif /* undef ALLOW_OPENAD */
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef ALLOW_AUTODIFF
      DO iloop = 1, nTimeSteps 
#endif /* ALLOW_AUTODIFF */

c--     >>> Loop body start <<<              循环体启动

#ifdef ALLOW_AUTODIFF
      nIter0 = NINT( (startTime-baseTime)/deltaTClock )
C--   Reset the model iteration counter and the model time.
      myIter = nIter0 + (iloop-1)
      myTime = startTime + deltaTClock*(iloop-1)
# ifndef ALLOW_OPENAD
      ikey_dynamics = ilev_1
# endif
#endif

#ifdef ALLOW_ECCO
      IF ( useECCO ) then
# ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL('cost_averagesfields',myThid)
# endif
c--   Accumulate time averages of temperature, salinity              温度、盐度的平均累积时间
        CALL TIMER_START('COST_AVERAGESFIELDS [MAIN_DO_LOOP]',myThid)
        CALL COST_AVERAGESFIELDS( myTime, myThid )                   平均成本函数场
        CALL TIMER_STOP ('COST_AVERAGESFIELDS [MAIN_DO_LOOP]',myThid)
      ENDIF
#endif /* ALLOW_ECCO */

#ifdef ALLOW_PROFILES
      IF (usePROFILES) THEN
#ifdef ALLOW_DEBUG
        IF (debugMode) CALL DEBUG_CALL('profiles_inloop',myThid)
#endif
c--     Accumulate in-situ time averages of theta, salt, and SSH.    累积θ、salt和SSH的原位（现场观测）时间平均值。
#ifdef ALLOW_AUTODIFF
C--   Reset the model iteration counter and the model time.   重置模型迭代计数器和模型时间
        myIter = nIter0 + (iloop-1)
        myTime = startTime + deltaTClock*(iloop-1)
#endif
        CALL TIMER_START('PROFILES_INLOOP    [MAIN_DO_LOOP]', myThid)
        CALL PROFILES_INLOOP( myTime, myThid )
        CALL TIMER_STOP ('PROFILES_INLOOP    [MAIN_DO_LOOP]', myThid)
      ENDIF
#endif

      CALL TIMER_START('MAIN_DO_LOOP        [THE_MAIN_LOOP]', myThid)
      CALL MAIN_DO_LOOP( iloop, myTime, myIter, myThid )                又一个下一级的主循环出来了MAIN_DO_LOOP
      CALL TIMER_STOP ('MAIN_DO_LOOP        [THE_MAIN_LOOP]', myThid)

c--     >>> Loop body end <<<   循环体结束
#ifndef ALLOW_OPENAD    如果没有定义openad
# ifdef ALLOW_AUTODIFF    如果定义了autodiff
#   ifdef ALLOW_TAMC_CHECKPOINTING   还定义了checkpoint       结束俩次if和do
            endif
          enddo
          endif
        enddo
#    ifndef AUTODIFF_2_LEVEL_CHECKPOINT     俩层check   结束if，do和if
        endif
      enddo
#    endif
#    ifdef AUTODIFF_4_LEVEL_CHECKPOINT      四层check   结束if，do和if
       endif
      enddo
#    endif
#   else /* ndef ALLOW_TAMC_CHECKPOINTING */       没有定义check的话，就结束DO和一次if
      ENDDO
#   endif /* ALLOW_TAMC_CHECKPOINTING */
# else /* ALLOW_AUTODIFF */   如果也没有定义autodiff的话
       land
                          就直接结束最重要的那次DO循环就行了，导向的是下一级的例程  MAIN_DO_LOOP( iloop, myTime, myIter, myThid )
# endif /* ALLOW_AUTODIFF */
#endif /* undef ALLOW_OPENAD */

c     >>>>>>>>>>>>>>>>>>>>>>>>>>>   LOOP   <<<<<<<<<<<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>>>>>>>>>  STOPS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<   循环结束

#ifdef ALLOW_ECCO 
c--     Accumulate time averages of temperature, salinity                温度、盐度的平均累积时间
      CALL TIMER_START('COST_AVERAGESFIELDS [THE_MAIN_LOOP]',myThid)
      CALL COST_AVERAGESFIELDS( endtime, myThid )
      CALL TIMER_STOP ('COST_AVERAGESFIELDS [THE_MAIN_LOOP]',myThid)
#endif /* ALLOW_ECCO */

#ifdef ALLOW_PROFILES
      IF (usePROFILES) THEN
c--   Accumulate in-situ time averages of temperature, salinity, SSH.   累积温度、盐度、SSH的原位时间平均值
        CALL TIMER_START('PROFILES_INLOOP    [THE_MAIN_LOOP]', myThid)
        CALL PROFILES_INLOOP( endtime, myThid )
        CALL TIMER_STOP ('PROFILES_INLOOP    [THE_MAIN_LOOP]', myThid)
c--   Compute all profiles cost function contributions.                 计算所有profiles的成本函数贡献
        CALL COST_PROFILES( myIter, myTime, myThid )
      ENDIF
#endif

#ifdef ALLOW_ECCO
      IF ( useECCO ) then
c--   Compute all ecco cost function contributions.                     计算所有ecco成本函数贡献
        CALL TIMER_START('ECCO_COST_DRIVER   [THE_MAIN_LOOP]', myThid)
        CALL ECCO_COST_DRIVER ( endTime, myIter, myThid )
        CALL TIMER_STOP ('ECCO_COST_DRIVER   [THE_MAIN_LOOP]', myThid)
      ENDIF
#endif /* ALLOW_ECCO */

#ifdef ALLOW_COST
c--   Sum all cost function contributions.             所有成本函数贡献之和
      CALL TIMER_START('COST_FINAL         [ADJOINT SPIN-DOWN]', myThid)
      CALL COST_FINAL ( myThid )                       COST_FINAL
      CALL TIMER_STOP ('COST_FINAL         [ADJOINT SPIN-DOWN]', myThid)

# ifdef ALLOW_OPENAD
c$openad DEPENDENT(fc)
# endif /* ALLOW_OPENAD */

#endif /* ALLOW_COST */

      _BARRIER
      CALL TIMER_STOP ('MAIN LOOP           [THE_MAIN_LOOP]', myThid)            计时结束

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('THE_MAIN_LOOP',myThid)
#endif

      RETURN
      END
