C $Header: /u/gcmpack/MITgcm/model/inc/PARAMS.h,v 1.286 2017/04/04 23:19:33 jmc Exp $
C $Name:  $
C

CBOP
C     !ROUTINE: PARAMS.h
C     !INTERFACE:
C     #include PARAMS.h

C     !DESCRIPTION:
C     Header file defining model "parameters".  The values from the
C     model standard input file are stored into the variables held
C     here. Notes describing the parameters can also be found here.
定义模型“参数”的头文件。模型标准输入文件中的值存储在这里保存的变量中。描述参数的注释也可以在这里找到。

CEOP

C--   Contants
C     Useful physical values
      Real*8 PI
      PARAMETER ( PI    = 3.14159265358979323844D0   )         pi
      Real*8 deg2rad
      PARAMETER ( deg2rad = 2.D0*PI/360.D0           )         转化角度和弧度

C--   COMMON /PARM_C/ Character valued parameters used by the model.   模型使用的字符值参数
C     buoyancyRelation :: Flag used to indicate which relation to use to   用于指示使用哪个关系来获取浮力的标志
C                         get buoyancy.
C     eosType         :: choose the equation of state:    选择状态方程
C                        LINEAR, POLY3, UNESCO, JMD95Z, JMD95P, MDJWF, IDEALGAS    这些中的一个
C     pickupSuff      :: force to start from pickup files (even if nIter0=0)   强制从pickup文件开始（即使nIter0=0）并读取具有此后缀的拾取文件
C                        and read pickup files with this suffix (max 10 Char.)
C     mdsioLocalDir   :: read-write tiled file from/to this directory name  将tiled 文件从 -该目录- 读入或者写入该目录 而非当前目录
C                        (+ 4 digits Processor-Rank) instead of current dir.
C     adTapeDir       :: read-write checkpointing tape files from/to this  将检查点tape 文件从 -该目录- 读入或者写入该目录 而非当前目录
C                        directory name instead of current dir. Conflicts  与上面那个冲突了，所以只能设置其中一个
C                        mdsioLocalDir, so only one of the two can be set.
C     tRefFile      :: File containing reference Potential Temperat.  tRef (1.D) 包含参考potential temperature的文件，位势温度，tRef（都是1D的）
C     sRefFile      :: File containing reference salinity/spec.humid. sRef (1.D) 包含参考盐度/大气中的一个量的文件
C     rhoRefFile    :: File containing reference density profile rhoRef (1.D)  包含参考密度分布的文件
C     gravityFile   :: File containing gravity vertical profile (1.D)          包含重力垂直分布的文件
C     delRFile      :: File containing vertical grid spacing delR  (1.D array) 包含垂直网格间距delR的文件
C     delRcFile     :: File containing vertical grid spacing delRc (1.D array) 包含垂直网格间距delRc的文件  中心
C     hybSigmFile   :: File containing hybrid-sigma vertical coord. coeff. (2x 1.D) 包含混合西格玛垂直坐标系数的文件 2*1D
C     delXFile      :: File containing X-spacing grid definition (1.D array)  包含X间距栅格定义的文件
C     delYFile      :: File containing Y-spacing grid definition (1.D array)  包含Y间距栅格定义的文件
C     horizGridFile :: File containing horizontal-grid definition   包含水平网格定义的文件（仅当使用曲线网格时）curvilinear_grid
C                        (only when using curvilinear_grid)
C     bathyFile       :: File containing bathymetry. If not defined bathymetry   包含水深测量的文件。如果未定义，则从内联函数获取水深测量。
C                        is taken from inline function.
C     topoFile        :: File containing the topography of the surface (unit=m)   包含表面地形的文件（单位=m）
C                        (mainly used for the atmosphere = ground height).        主要用于大气=地面高度
C     addWwallFile    :: File containing 2-D additional Western  cell-edge wall   包含二维附加西部单元边缘墙的文件
C     addSwallFile    :: File containing 2-D additional Southern cell-edge wall   包含二维附加南部单元边缘墙的文件
C                        (e.g., to add "thin-wall" where it is =1)     在1处添加薄壁
C     hydrogThetaFile :: File containing initial hydrographic data (3-D)  包含位温的初始水文数据的文件 三维的
C                        for potential temperature.
C     hydrogSaltFile  :: File containing initial hydrographic data (3-D)  包含盐度的初始水文数据的文件 三维的
C                        for salinity.
C     diffKrFile      :: File containing 3D specification of vertical diffusivity  包含垂直扩散率三维规范的文件
C     viscAhDfile     :: File containing 3D specification of horizontal viscosity  包含水平粘度三维规范的文件
C     viscAhZfile     :: File containing 3D specification of horizontal viscosity
C     viscA4Dfile     :: File containing 3D specification of horizontal viscosity
C     viscA4Zfile     :: File containing 3D specification of horizontal viscosity
C     zonalWindFile   :: File containing zonal wind data   包含纬向风数据的文件
C     meridWindFile   :: File containing meridional wind data  包含经向风数据的文件
C     thetaClimFile   :: File containing surface theta climataology used   包含松弛项中使用的表面θ气候学的文件    lambda(theta-theta*)
C                        in relaxation term -lambda(theta-theta*)
C     saltClimFile    :: File containing surface salt climataology used    包含松弛项中使用的表面盐度气候学的文件    lambda(theta-theta*)
C                        in relaxation term -lambda(salt-salt*)
C     surfQfile       :: File containing surface heat flux, excluding SW   包含表面热流的文件，不包括SW（旧版本，为向后兼容而保留）
C                        (old version, kept for backward compatibility)
C     surfQnetFile    :: File containing surface net heat flux             包含表面净热通量的文件
C     surfQswFile     :: File containing surface shortwave radiation       包含表面短波辐射的文件
C     EmPmRfile       :: File containing surface fresh water flux          包含地表淡水通量的文件
C           NOTE: for backward compatibility EmPmRfile is specified in
C                 m/s when using external_fields_load.F.  It is converted
C                 to kg/m2/s by multiplying by rhoConstFresh.
注意：为了向后兼容，当使用external_fields_load.F时，EmPmRfile以m/s为单位指定。It is converted to kg/m2/s by multiplying by rhoConstFresh.
C     saltFluxFile    :: File containing surface salt flux    表面盐通量
C     pLoadFile       :: File containing pressure loading     包含压力加载的文件
C     addMassFile     :: File containing source/sink of fluid in the interior    包含内部流体源/汇的文件
C     eddyPsiXFile    :: File containing zonal Eddy streamfunction data     包含纬向涡流函数数据的文件
C     eddyPsiYFile    :: File containing meridional Eddy streamfunction data包含子午涡流函数数据的文件
C     the_run_name    :: string identifying the name of the model "run"  标识模型“run”名称的字符串
      COMMON /PARM_C/
     &                buoyancyRelation, eosType,
     &                pickupSuff, mdsioLocalDir, adTapeDir,
     &                tRefFile, sRefFile, rhoRefFile, gravityFile,
     &                delRFile, delRcFile, hybSigmFile,
     &                delXFile, delYFile, horizGridFile,
     &                bathyFile, topoFile, addWwallFile, addSwallFile,
     &                viscAhDfile, viscAhZfile,
     &                viscA4Dfile, viscA4Zfile,
     &                hydrogThetaFile, hydrogSaltFile, diffKrFile,
     &                zonalWindFile, meridWindFile, thetaClimFile,
     &                saltClimFile,
     &                EmPmRfile, saltFluxFile,
     &                surfQfile, surfQnetFile, surfQswFile,
     &                lambdaThetaFile, lambdaSaltFile,
     &                uVelInitFile, vVelInitFile, pSurfInitFile,
     &                pLoadFile, addMassFile,
     &                eddyPsiXFile, eddyPsiYFile, geothermalFile,
     &                the_run_name
      CHARACTER*(MAX_LEN_FNAM) buoyancyRelation
      CHARACTER*(6)  eosType
      CHARACTER*(10) pickupSuff
      CHARACTER*(MAX_LEN_FNAM) mdsioLocalDir
      CHARACTER*(MAX_LEN_FNAM) adTapeDir
      CHARACTER*(MAX_LEN_FNAM) tRefFile
      CHARACTER*(MAX_LEN_FNAM) sRefFile
      CHARACTER*(MAX_LEN_FNAM) rhoRefFile
      CHARACTER*(MAX_LEN_FNAM) gravityFile
      CHARACTER*(MAX_LEN_FNAM) delRFile
      CHARACTER*(MAX_LEN_FNAM) delRcFile
      CHARACTER*(MAX_LEN_FNAM) hybSigmFile
      CHARACTER*(MAX_LEN_FNAM) delXFile
      CHARACTER*(MAX_LEN_FNAM) delYFile
      CHARACTER*(MAX_LEN_FNAM) horizGridFile
      CHARACTER*(MAX_LEN_FNAM) bathyFile, topoFile
      CHARACTER*(MAX_LEN_FNAM) addWwallFile, addSwallFile
      CHARACTER*(MAX_LEN_FNAM) hydrogThetaFile, hydrogSaltFile
      CHARACTER*(MAX_LEN_FNAM) diffKrFile
      CHARACTER*(MAX_LEN_FNAM) viscAhDfile
      CHARACTER*(MAX_LEN_FNAM) viscAhZfile
      CHARACTER*(MAX_LEN_FNAM) viscA4Dfile
      CHARACTER*(MAX_LEN_FNAM) viscA4Zfile
      CHARACTER*(MAX_LEN_FNAM) zonalWindFile
      CHARACTER*(MAX_LEN_FNAM) meridWindFile
      CHARACTER*(MAX_LEN_FNAM) thetaClimFile
      CHARACTER*(MAX_LEN_FNAM) saltClimFile
      CHARACTER*(MAX_LEN_FNAM) surfQfile
      CHARACTER*(MAX_LEN_FNAM) surfQnetFile
      CHARACTER*(MAX_LEN_FNAM) surfQswFile
      CHARACTER*(MAX_LEN_FNAM) EmPmRfile
      CHARACTER*(MAX_LEN_FNAM) saltFluxFile
      CHARACTER*(MAX_LEN_FNAM) uVelInitFile
      CHARACTER*(MAX_LEN_FNAM) vVelInitFile
      CHARACTER*(MAX_LEN_FNAM) pSurfInitFile
      CHARACTER*(MAX_LEN_FNAM) pLoadFile
      CHARACTER*(MAX_LEN_FNAM) addMassFile
      CHARACTER*(MAX_LEN_FNAM) eddyPsiXFile
      CHARACTER*(MAX_LEN_FNAM) eddyPsiYFile
      CHARACTER*(MAX_LEN_FNAM) geothermalFile
      CHARACTER*(MAX_LEN_FNAM) lambdaThetaFile
      CHARACTER*(MAX_LEN_FNAM) lambdaSaltFile
      CHARACTER*(MAX_LEN_PREC/2) the_run_name

C--   COMMON /PARM_I/ Integer valued parameters used by the model.   模型使用的整数值参数
C     cg2dMaxIters        :: Maximum number of iterations in the     最大迭代次数，在二维梯度解算器中
C                            two-dimensional con. grad solver.
C     cg2dChkResFreq      :: Frequency with which to check residual  检查残差的频率
C                            in con. grad solver.
C     cg2dPreCondFreq     :: Frequency for updating cg2d preconditioner  更新cg2d预处理器的频率(非线性自由表面)
C                            (non-linear free-surf.)
C     cg2dUseMinResSol    :: =0 : use last-iteration/converged solution  使用上一次迭代/收敛解
C                            =1 : use solver minimum-residual solution   使用解算器最小剩余解
C     cg3dMaxIters        :: Maximum number of iterations in the     最大迭代次数，在三维梯度解算器中
C                            three-dimensional con. grad solver.
C     cg3dChkResFreq      :: Frequency with which to check residual  检查残差的频率
C                            in con. grad solver.
C     printResidualFreq   :: Frequency for printing residual in CG iterations  CG迭代中打印残差的频率
C     nIter0              :: Start time-step number of for this run       此运行的开始时间步数
C     nTimeSteps          :: Number of timesteps to execute               要执行的时间步数
C     nTimeSteps_l2       :: Number of inner timesteps to execute per timestep  每个时间步要执行的内部时间步数
C     selectCoriMap       :: select setting of Coriolis parameter map:   选择科里奥利参数的设置
C                           =0 f-Plane (Constant Coriolis, = f0)    f平面近似
C                           =1 Beta-Plane Coriolis (= f0 + beta.y)  β平面近似
C                           =2 Spherical Coriolis (= 2.omega.sin(phi))
C                           =3 Read Coriolis 2-d fields from files. 从文件中读取
C     selectSigmaCoord    :: option related to sigma vertical coordinate   与垂直坐标有关的选项
C     nonlinFreeSurf      :: option related to non-linear free surface     与非线性自由面有关的选项（0线性，1非线性）
C                           =0 Linear free surface ; >0 Non-linear
C     select_rStar        :: option related to r* vertical coordinate      与r*垂直坐标相关的选项（0用r,>0用r*）
C                           =0 (default) use r coord. ; > 0 use r*
C     selectNHfreeSurf    :: option for Non-Hydrostatic (free-)Surface formulation:   非静水（自由）表面配方选项
C                           =0 (default) hydrostatic surf. ; > 0 add NH effects.
C     selectP_inEOS_Zc    :: select which pressure to use in EOS (for z-coords)     选择EOS中要使用的压力（对于z坐标）
C                           =0: simply: -g*rhoConst*z    简单的rgh
C                           =1: use pRef = integral{-g*rho(Tref,Sref,pRef)*dz}
C                           =2: use hydrostatic dynamical pressure     静水动压力
C                           =3: use full (Hyd+NH) dynamical pressure   使用全（Hyd+NH）动态压力
C     selectAddFluid      :: option to add mass source/sink of fluid in the interior   在内部添加流体质量源/汇的选项
C                            (3-D generalisation of oceanic real-fresh water flux)     海洋真实淡水通量的三维概化
C                           =0 off ; =1 add fluid ; =-1 virtual flux (no mass added)   关闭、加流体、虚通量（未加质量）
C     selectImplicitDrag  :: select Implicit treatment of bottom/top drag    选择“下/上拖动的隐式处理”
C                           = 0: fully explicit    完全显式的
C                           = 1: implicit on provisional velocity  隐式-临时速度
C                                (i.e., before grad.Eta increment)  （即）在梯度Eta增量前
C                           = 2: fully implicit (combined with Impl Surf.Press)  完全隐式
C     momForcingOutAB     :: =1: take momentum forcing contribution   从Adams-Bashforth time stepping中取出或者加入动量强迫贡献。
C                            out of (=0: in) Adams-Bashforth time stepping.
C     tracForcingOutAB    :: =1: take tracer (Temp,Salt,pTracers) forcing contribution
C                            out of (=0: in) Adams-Bashforth time stepping. 从Adams-Bashforth time stepping中取出或者加入tracer (Temp,Salt,pTracers)强迫贡献。
C     tempAdvScheme       :: Temp. Horiz.Advection scheme selector  温度水平平流方案选择器
C     tempVertAdvScheme   :: Temp. Vert. Advection scheme selector  温度垂直
C     saltAdvScheme       :: Salt. Horiz.advection scheme selector  盐度水平
C     saltVertAdvScheme   :: Salt. Vert. Advection scheme selector  盐度垂直
C     selectKEscheme      :: Kinetic Energy scheme selector (Vector Inv.)   动能方案选择器（矢量输入）
C     selectVortScheme    :: Scheme selector for Vorticity term (Vector Inv.)涡度项的方案选择器（矢量输入）
C     selectBotDragQuadr  :: quadratic bottom drag discretisation option:   二次底拖曳离散化选项：
C                           =0: average KE from grid center to U & V location  从网格中心到U&V位置的平均KE
C                           =1: use local velocity norm @ U & V location    在U&V位置使用当地速度标准
C                           =2: same with wet-point averaging of other component  与其他组分的湿点平均值相同
C     readBinaryPrec      :: Precision used for reading binary files  用于读取二进制文件的Precision 32 64
C     writeStatePrec      :: Precision used for writing model state.  用于写入模型状态的Precision
C     writeBinaryPrec     :: Precision used for writing binary files  用于写入二进制文件的Precision
C     rwSuffixType        :: controls the format of the mds file suffix.  控制mds文件后缀的格式
C                          =0 (default): use iteration number (myIter, I10.10);  使用迭代数（myIter，I10.10）；
C                          =1: 100*myTime (100th sec); =2: myTime (seconds);     100*myTime（第100秒）=2:我的时间（秒）；
C                          =3: myTime/360 (10th of hr); =4: myTime/3600 (hours). myTime/360（小时的10分之一）=4:myTime/3600（小时）。
C     monitorSelect       :: select group of variables to monitor      选择要监视的变量组
C                            =1 : dynvars ; =2 : + vort ; =3 : + surface    
C-    debugLevel          :: controls printing of algorithm intermediate results    控制算法中间结果和统计的打印
C                            and statistics ; higher -> more writing
C-    plotLevel           :: controls printing of field maps ; higher -> more flds  控制场地图的打印

      COMMON /PARM_I/
     &        cg2dMaxIters, cg2dChkResFreq,
     &        cg2dPreCondFreq, cg2dUseMinResSol,
     &        cg3dMaxIters, cg3dChkResFreq,
     &        printResidualFreq,
     &        nIter0, nTimeSteps, nTimeSteps_l2, nEndIter,
     &        selectCoriMap,
     &        selectSigmaCoord,
     &        nonlinFreeSurf, select_rStar,
     &        selectNHfreeSurf, selectP_inEOS_Zc,
     &        selectAddFluid, selectImplicitDrag,
     &        momForcingOutAB, tracForcingOutAB,
     &        tempAdvScheme, tempVertAdvScheme,
     &        saltAdvScheme, saltVertAdvScheme,
     &        selectKEscheme, selectVortScheme,
     &        selectBotDragQuadr,
     &        readBinaryPrec, writeBinaryPrec, writeStatePrec,
     &        rwSuffixType, monitorSelect, debugLevel, plotLevel
      INTEGER cg2dMaxIters
      INTEGER cg2dChkResFreq
      INTEGER cg2dPreCondFreq
      INTEGER cg2dUseMinResSol
      INTEGER cg3dMaxIters
      INTEGER cg3dChkResFreq
      INTEGER printResidualFreq
      INTEGER nIter0
      INTEGER nTimeSteps
      INTEGER nTimeSteps_l2
      INTEGER nEndIter
      INTEGER selectCoriMap
      INTEGER selectSigmaCoord
      INTEGER nonlinFreeSurf
      INTEGER select_rStar
      INTEGER selectNHfreeSurf
      INTEGER selectP_inEOS_Zc
      INTEGER selectAddFluid
      INTEGER selectImplicitDrag
      INTEGER momForcingOutAB, tracForcingOutAB
      INTEGER tempAdvScheme, tempVertAdvScheme
      INTEGER saltAdvScheme, saltVertAdvScheme
      INTEGER selectKEscheme
      INTEGER selectVortScheme
      INTEGER selectBotDragQuadr
      INTEGER readBinaryPrec
      INTEGER writeStatePrec
      INTEGER writeBinaryPrec
      INTEGER rwSuffixType
      INTEGER monitorSelect
      INTEGER debugLevel
      INTEGER plotLevel

C--   COMMON /PARM_L/ Logical valued parameters used by the model.  模型使用的逻辑值参数。
C- Coordinate + Grid params:   坐标轴和网格参数
C     fluidIsAir       :: Set to indicate that the fluid major constituent  指示流体的主要成分是空气
C                         is air
C     fluidIsWater     :: Set to indicate that the fluid major constituent  指示流体的主要成分是水
C                         is water
C     usingPCoords     :: Set to indicate that we are working in a pressure 设置为指示我们在压力类型坐标（p或p*）中工作。
C                         type coordinate (p or p*).
C     usingZCoords     :: Set to indicate that we are working in a height   设置为指示我们在高度类型坐标（z或z*）中工作
C                         type coordinate (z or z*)
C     usingCartesianGrid :: If TRUE grid generation will be in a cartesian  如果真  网格生成将在笛卡尔坐标系中进行。
C                           coordinate frame.
C     usingSphericalPolarGrid :: If TRUE grid generation will be in a       如果真  网格生成将在球形极坐标系中进行
C                                spherical polar frame.
C     rotateGrid      :: rotate grid coordinates to geographical coordinates 根据Euler角phiEuler、thetaEuler、psiEuler将栅格坐标旋转到地理坐标
C                        according to Euler angles phiEuler, thetaEuler, psiEuler
C     usingCylindricalGrid :: If TRUE grid generation will be Cylindrical   如果真  网格生成将是圆柱形的
C     usingCurvilinearGrid :: If TRUE, use a curvilinear grid (to be provided)  如果为真，则使用曲线网格（待提供）
C     hasWetCSCorners :: domain contains CS-type corners where dynamics is solved  域包含求解动力学的CS型角点
C     deepAtmosphere :: deep model (drop the shallow-atmosphere approximation)  深模式（浅层大气近似）
C     setInterFDr    :: set Interface depth (put cell-Center at the middle)     设置界面深度（将单元格中心置于中间）
C     setCenterDr    :: set cell-Center depth (put Interface at the middle)     设置单元格中心深度（将接口置于中间）
C- Momentum params:  动量方程参数
C     no_slip_sides  :: Impose "no-slip" at lateral boundaries.    在横向边界施加“无滑移”
C     no_slip_bottom :: Impose "no-slip" at bottom boundary.       在底部边界施加“无滑移”
C     bottomVisc_pCell :: account for partial-cell in bottom visc. (no-slip BC)      （无滑移边界条件）
C     useSmag3D      :: Use isotropic 3-D Smagorinsky 使用各向同性三维Smagorinsky
C     useFullLeith   :: Set to true to use full Leith viscosity(may be unstable  设置为true以使用full Leith粘度（在不规则网格上可能不稳定）
C                       on irregular grids)
C     useStrainTensionVisc:: Set to true to use Strain-Tension viscous terms  设置为true以使用应变-张力-粘性项
C     useAreaViscLength :: Set to true to use old scaling for viscous lengths, 如果设置为true，则对粘性长度使用旧的缩放
C                          e.g., L2=Raz.  May be preferable for cube sphere.   e、 g.，L2=拉兹。可能更适合立方体球体。
C     momViscosity  :: Flag which turns momentum friction terms on and off.   打开和关闭动量摩擦项的标志。
C     momAdvection  :: Flag which turns advection of momentum on and off.     打开和关闭动量平流的标志。
C     momForcing    :: Flag which turns external forcing of momentum on       打开和关闭外部动力的标志。
C                      and off.
C     momPressureForcing :: Flag which turns pressure term in momentum equation 打开和关闭动量方程中的压力项的标志。
C                          on and off.
C     metricTerms   :: Flag which turns metric terms on or off.   打开或关闭metric terms的标志。
C     useNHMTerms   :: If TRUE use non-hydrostatic metric terms.  如果为真，则使用非流体静力metric terms。度量系数
C     useCoriolis   :: Flag which turns the coriolis terms on and off.   打开和关闭科里奥利术语的标志
C     use3dCoriolis :: Turns the 3-D coriolis terms (in Omega.cos Phi) on - off   打开-关闭三维科里奥利项（在ωcos Phi中）
C     useCDscheme   :: use CD-scheme to calculate Coriolis terms.  用CD格式计算科里奥利项。
C     vectorInvariantMomentum :: use Vector-Invariant form (mom_vecinv package)  使用向量不变形式（mom_vecinv package）
C                                (default = F = use mom_fluxform package)        (default = F = use mom_fluxform package)
C     useJamartWetPoints :: Use wet-point method for Coriolis (Jamart & Ozer 1986) 用湿点法测定科里奥利
C     useJamartMomAdv :: Use wet-point method for V.I. non-linear term   V.I.非线性项采用湿点法
C     upwindVorticity :: bias interpolation of vorticity in the Coriolis term  科里奥利项中涡度的偏差插值
C     highOrderVorticity :: use 3rd/4th order interp. of vorticity (V.I., advection)   使用3/4阶插值。涡度（V.I.，平流） 
C     useAbsVorticity :: work with f+zeta in Coriolis terms                work with f+zeta in Coriolis terms
C     upwindShear     :: use 1rst order upwind interp. (V.I., vertical advection)   使用一阶逆风对讲机(垂直平流）
C     momStepping    :: Turns momentum equation time-stepping off  关闭动量方程时间步进
C     calc_wVelocity :: Turns vertical velocity calculation off    关闭垂直速度计算
C- Temp. & Salt params:   温盐参数
C     tempStepping   :: Turns temperature equation time-stepping on/off  打开/关闭温度方程时间步进
C     saltStepping   :: Turns salinity equation time-stepping on/off  打开/关闭盐度方程时间步进
C     addFrictionHeating :: account for frictional heating  考虑摩擦加热
C     tempAdvection  :: Flag which turns advection of temperature on and off.  打开和关闭温度平流的标志。
C     tempVertDiff4  :: use vertical bi-harmonic diffusion for temperature   对温度使用垂直双谐波扩散
C     tempIsActiveTr :: Pot.Temp. is a dynamically active tracer   Pot.Temp。是一个动态活动的跟踪器
C     tempForcing    :: Flag which turns external forcing of temperature on/off  打开/关闭外部温度强制的标志
C     saltAdvection  :: Flag which turns advection of salinity on and off.  打开和关闭盐度平流的标志。
C     saltVertDiff4  :: use vertical bi-harmonic diffusion for salinity  对盐度使用垂直双谐波扩散
C     saltIsActiveTr :: Salinity  is a dynamically active tracer   盐度是一种动态活性示踪剂
C     saltForcing    :: Flag which turns external forcing of salinity on/off   打开/关闭外部盐度强制的标志
C     maskIniTemp    :: apply mask to initial Pot.Temp.     apply mask to initial Pot.Temp.
C     maskIniSalt    :: apply mask to initial salinity      apply mask to initial salinity
C     checkIniTemp   :: check for points with identically zero initial Pot.Temp.   检查初始位温相同为零的点
C     checkIniSalt   :: check for points with identically zero initial salinity    检查初始盐度相同为零的点
C- Pressure solver related parameters (PARM02)     压力解算器相关参数
C     useSRCGSolver  :: Set to true to use conjugate gradient
C                       solver with single reduction (only one call of
C                       s/r mpi_allreduce), default is false
如果设置为true，则使用带有单个缩减的共轭梯度解算器 （仅调用一次s/r mpi\u allreduce），默认值为false
C- Time-stepping & free-surface params:            时间步长和自由曲面参数：
C     rigidLid            :: Set to true to use rigid lid       设置为true以使用刚性lid
C     implicitFreeSurface :: Set to true to use implicit free surface   设置为true以使用隐式自由曲面
C     uniformLin_PhiSurf  :: Set to true to use a uniform Bo_surf in the
C                            linear relation Phi_surf = Bo_surf*eta
C     uniformFreeSurfLev  :: TRUE if free-surface level-index is uniform (=1)  如果自由表面水平指数是均匀的，则为真（=1）
C     exactConserv        :: Set to true to conserve exactly the total Volume   设置为true以精确保存总体积
C     linFSConserveTr     :: Set to true to correct source/sink of tracer  设置为true以校正由于线性自由表面而导致的表面跟踪器的源/汇
C                            at the surface due to Linear Free Surface
C     useRealFreshWaterFlux :: if True (=Natural BCS), treats P+R-E flux 如果为真（=自然BCS），则将P+R-E通量视为真正的淡水（=>改变海平面）如果为F，则将P+R-E转化为盐通量（无SL效应）
C                         as a real Fresh Water (=> changes the Sea Level)
C                         if F, converts P+R-E to salt flux (no SL effect)
C     storePhiHyd4Phys :: store hydrostatic potential for use in Physics/EOS  存储用于物理/EOS的静水压势
C                         this requires specific code for restart & exchange   这需要重新启动和交换的特定代码
C     quasiHydrostatic :: Using non-hydrostatic terms in hydrostatic algorithm  在静水力学算法中使用非静水力学项
C     nonHydrostatic   :: Using non-hydrostatic algorithm   使用非静力学算法
C     use3Dsolver      :: set to true to use 3-D pressure solver    设置为true以使用三维压力解算器
C     implicitIntGravWave :: treat Internal Gravity Wave implicitly  隐式处理重力内波
C     staggerTimeStep   :: enable a Stagger time stepping U,V (& W) then T,S   启用交错时间步进U、V和W，然后启用T、S
C     applyExchUV_early :: Apply EXCH to U,V earlier, just before integr_continuity  将EXCH应用于U，V early，正好在integr_continuity之前
C     doResetHFactors   :: Do reset thickness factors @ beginning of each time-step  在每个时间步开始时重置厚度因子
C     implicitDiffusion :: Turns implicit vertical diffusion on  打开隐式垂直扩散
C     implicitViscosity :: Turns implicit vertical viscosity on  打开隐式垂直粘度
C     tempImplVertAdv   :: Turns on implicit vertical advection for Temperature  打开温度的隐式垂直平流
C     saltImplVertAdv   :: Turns on implicit vertical advection for Salinity     打开盐度的隐式垂直平流
C     momImplVertAdv    :: Turns on implicit vertical advection for Momentum     打开动量的隐式垂直平流
C     multiDimAdvection :: Flag that enable multi-dimension advection            启用多维平流的标志
C     useMultiDimAdvec  :: True if multi-dim advection is used at least once     如果至少使用一次多维平流，则为True
C     momDissip_In_AB   :: if False, put Dissipation tendency contribution  如果为假，则将耗散趋势贡献从Adams bashorth时间步中去掉。
C                          out off Adams-Bashforth time stepping.
C     doAB_onGtGs       :: if the Adams-Bashforth time stepping is used, always  如果使用了Adams-Bashforth时间步进，总是在tracer趋势上应用AB（而不是在tracer上）
C                          apply AB on tracer tendencies (rather than on Tracer)
C- Other forcing params -    其他强制参数
C     balanceEmPmR    :: substract global mean of EmPmR at every time step   减去每个时间步EmPmR的全局平均值
C     balanceQnet     :: substract global mean of Qnet at every time step    减去每个时间步Qnet的全局平均值
C     balancePrintMean:: print substracted global means to STDOUT   将减去的全局平均值打印到标准输出
C     doThetaClimRelax :: Set true if relaxation to temperature   如果需要对温度气候进行松弛，则将其设置为true。
C                        climatology is required.
C     doSaltClimRelax  :: Set true if relaxation to salinity   如果需要对盐度气候进行松弛，则将其设置为true。
C                        climatology is required.
C     balanceThetaClimRelax :: substract global mean effect at every time step  减去每个时间步的全局平均效果
C     balanceSaltClimRelax :: substract global mean effect at every time step   减去每个时间步的全局平均效果
C     allowFreezing  :: Allows surface water to freeze and form ice  允许地表水冻结并形成冰周期
C     periodicExternalForcing :: Set true if forcing is time-dependant   如果强制与时间相关，则设置为true
C- I/O parameters -   I/O参数
C     globalFiles    :: Selects between "global" and "tiled" files.
C                       On some platforms with MPI, option globalFiles is either
C                       slow or does not work. Use useSingleCpuIO instead.
在“全局”和“平铺”文件之间进行选择。在一些使用MPI的平台上，选项globalFiles要么速度慢要么不起作用。改用useSingleCpuIO。
C     useSingleCpuIO :: moved to EEPARAMS.h  搬到EEPARAMS.h了
C     pickupStrictlyMatch :: check and stop if pickup-file do not stricly match  如果pickup文件不完全匹配，请检查并停止
C     startFromPickupAB2 :: with AB-3 code, start from an AB-2 pickup  使用AB-3代码，从AB-2pickup 开始
C     usePickupBeforeC54 :: start from old-pickup files, generated with code from
C                           before checkpoint-54a, Jul 06, 2004.  从旧的拾取文件开始，使用代码从checkpoint-54a之前生成，2004年7月6日。
C     pickup_write_mdsio :: use mdsio to write pickups  使用mdsio写入pickup
C     pickup_read_mdsio  :: use mdsio to read  pickups  使用mdsio读取pickup
C     pickup_write_immed :: echo the pickup immediately (for conversion)   echo重复 映现  the pickup immediately   用于转换
C     writePickupAtEnd   :: write pickup at the last timestep   在最后一个时间步写入拾取
C     timeave_mdsio      :: use mdsio for timeave output   将mdsio用于timeave输出
C     snapshot_mdsio     :: use mdsio for "snapshot" (dumpfreq/diagfreq) output   将mdsio用于“快照”（dumpfreq/diagfreq）输出
C     monitor_stdio      :: use stdio for monitor output   使用stdio作为（monitor）监视器输出
C     dumpInitAndLast :: dumps model state to files at Initial (nIter0)      除了dumpFreq iter的倍数之外，在初始（nIter0）最后一次迭代时将模型状态转储到文件中。
C                        & Last iteration, in addition multiple of dumpFreq iter.

      COMMON /PARM_L/
     & fluidIsAir, fluidIsWater,
     & usingPCoords, usingZCoords,
     & usingCartesianGrid, usingSphericalPolarGrid, rotateGrid,
     & usingCylindricalGrid, usingCurvilinearGrid, hasWetCSCorners,
     & deepAtmosphere, setInterFDr, setCenterDr,
     & no_slip_sides, no_slip_bottom, bottomVisc_pCell, useSmag3D,
     & useFullLeith, useStrainTensionVisc, useAreaViscLength,
     & momViscosity, momAdvection, momForcing,
     & momPressureForcing, metricTerms, useNHMTerms,
     & useCoriolis, use3dCoriolis,
     & useCDscheme, vectorInvariantMomentum,
     & useEnergyConservingCoriolis, useJamartWetPoints, useJamartMomAdv,
     & upwindVorticity, highOrderVorticity,
     & useAbsVorticity, upwindShear,
     & momStepping, calc_wVelocity, tempStepping, saltStepping,
     & addFrictionHeating,
     & tempAdvection, tempVertDiff4, tempIsActiveTr, tempForcing,
     & saltAdvection, saltVertDiff4, saltIsActiveTr, saltForcing,
     & maskIniTemp, maskIniSalt, checkIniTemp, checkIniSalt,
     & useSRCGSolver,
     & rigidLid, implicitFreeSurface,
     & uniformLin_PhiSurf, uniformFreeSurfLev,
     & exactConserv, linFSConserveTr, useRealFreshWaterFlux,
     & storePhiHyd4Phys, quasiHydrostatic, nonHydrostatic,
     & use3Dsolver, implicitIntGravWave, staggerTimeStep,
     & applyExchUV_early, doResetHFactors,
     & implicitDiffusion, implicitViscosity,
     & tempImplVertAdv, saltImplVertAdv, momImplVertAdv,
     & multiDimAdvection, useMultiDimAdvec,
     & momDissip_In_AB, doAB_onGtGs,
     & balanceEmPmR, balanceQnet, balancePrintMean,
     & balanceThetaClimRelax, balanceSaltClimRelax,
     & doThetaClimRelax, doSaltClimRelax,
     & allowFreezing,
     & periodicExternalForcing,
     & globalFiles,
     & pickupStrictlyMatch, usePickupBeforeC54, startFromPickupAB2,
     & pickup_read_mdsio, pickup_write_mdsio, pickup_write_immed,
     & writePickupAtEnd,
     & timeave_mdsio, snapshot_mdsio, monitor_stdio,
     & outputTypesInclusive, dumpInitAndLast

      LOGICAL fluidIsAir
      LOGICAL fluidIsWater
      LOGICAL usingPCoords
      LOGICAL usingZCoords
      LOGICAL usingCartesianGrid
      LOGICAL usingSphericalPolarGrid, rotateGrid
      LOGICAL usingCylindricalGrid
      LOGICAL usingCurvilinearGrid, hasWetCSCorners
      LOGICAL deepAtmosphere
      LOGICAL setInterFDr
      LOGICAL setCenterDr

      LOGICAL no_slip_sides
      LOGICAL no_slip_bottom
      LOGICAL bottomVisc_pCell
      LOGICAL useSmag3D
      LOGICAL useFullLeith
      LOGICAL useStrainTensionVisc
      LOGICAL useAreaViscLength
      LOGICAL momViscosity
      LOGICAL momAdvection
      LOGICAL momForcing
      LOGICAL momPressureForcing
      LOGICAL metricTerms
      LOGICAL useNHMTerms

      LOGICAL useCoriolis
      LOGICAL use3dCoriolis
      LOGICAL useCDscheme
      LOGICAL vectorInvariantMomentum
      LOGICAL useEnergyConservingCoriolis
      LOGICAL useJamartWetPoints
      LOGICAL useJamartMomAdv
      LOGICAL upwindVorticity
      LOGICAL highOrderVorticity
      LOGICAL useAbsVorticity
      LOGICAL upwindShear
      LOGICAL momStepping
      LOGICAL calc_wVelocity
      LOGICAL tempStepping
      LOGICAL saltStepping
      LOGICAL addFrictionHeating
      LOGICAL tempAdvection
      LOGICAL tempVertDiff4
      LOGICAL tempIsActiveTr
      LOGICAL tempForcing
      LOGICAL saltAdvection
      LOGICAL saltVertDiff4
      LOGICAL saltIsActiveTr
      LOGICAL saltForcing
      LOGICAL maskIniTemp
      LOGICAL maskIniSalt
      LOGICAL checkIniTemp
      LOGICAL checkIniSalt
      LOGICAL useSRCGSolver
      LOGICAL rigidLid
      LOGICAL implicitFreeSurface
      LOGICAL uniformLin_PhiSurf
      LOGICAL uniformFreeSurfLev
      LOGICAL exactConserv
      LOGICAL linFSConserveTr
      LOGICAL useRealFreshWaterFlux
      LOGICAL storePhiHyd4Phys
      LOGICAL quasiHydrostatic
      LOGICAL nonHydrostatic
      LOGICAL use3Dsolver
      LOGICAL implicitIntGravWave
      LOGICAL staggerTimeStep
      LOGICAL applyExchUV_early
      LOGICAL doResetHFactors
      LOGICAL implicitDiffusion
      LOGICAL implicitViscosity
      LOGICAL tempImplVertAdv
      LOGICAL saltImplVertAdv
      LOGICAL momImplVertAdv
      LOGICAL multiDimAdvection
      LOGICAL useMultiDimAdvec
      LOGICAL momDissip_In_AB
      LOGICAL doAB_onGtGs
      LOGICAL balanceEmPmR
      LOGICAL balanceQnet
      LOGICAL balancePrintMean
      LOGICAL doThetaClimRelax
      LOGICAL doSaltClimRelax
      LOGICAL balanceThetaClimRelax
      LOGICAL balanceSaltClimRelax
      LOGICAL allowFreezing
      LOGICAL periodicExternalForcing
      LOGICAL globalFiles
      LOGICAL pickupStrictlyMatch
      LOGICAL usePickupBeforeC54
      LOGICAL startFromPickupAB2
      LOGICAL pickup_read_mdsio, pickup_write_mdsio
      LOGICAL pickup_write_immed, writePickupAtEnd
      LOGICAL timeave_mdsio, snapshot_mdsio, monitor_stdio
      LOGICAL outputTypesInclusive
      LOGICAL dumpInitAndLast

C--   COMMON /PARM_R/ "Real" valued parameters used by the model.   模型使用的“实值”参数
C     cg2dTargetResidual
C          :: Target residual for cg2d solver; no unit (RHS normalisation)   cg2d解算器的目标残差；没有单位
C     cg2dTargetResWunit
C          :: Target residual for cg2d solver; W unit (No RHS normalisation)  cg2d解算器的目标残差；W 单位（无RHS标准化）
C     cg3dTargetResidual
C               :: Target residual for cg3d solver.      cg3d解算器的目标残差。
C     cg2dpcOffDFac :: Averaging weight for preconditioner off-diagonal.   平均权的预处理器关闭对角线
C     Note. 20th May 1998
我有个奇怪的发现！在模型文件中，我们论证了这里使用的预处理器的形式（见有限体积、不可压缩的Navier-Stokes模型……，Marshall等人）。
代数给出了一个简单的0.5因子来计算ac和aCw的平均值，从而得到一个对称的预调节器。通过使用一个0.51的因子，也就是说，将预条件中的非对角项稍微缩小一点，
我成功地得到了一个测试用例中收敛的迭代次数，从192到134！需要进一步调查！目前，我引入了一个参数cg2dpcOffDFac，默认值为0.51，但可以在运行时设置。
C           I made a weird discovery! In the model paper we argue
C           for the form of the preconditioner used here ( see
C           A Finite-volume, Incompressible Navier-Stokes Model
C           ...., Marshall et. al ). The algebra gives a simple
C           0.5 factor for the averaging of ac and aCw to get a
C           symmettric pre-conditioner. By using a factor of 0.51
C           i.e. scaling the off-diagonal terms in the
C           preconditioner down slightly I managed to get the
C           number of iterations for convergence in a test case to
C           drop form 192 -> 134! Need to investigate this further!
C           For now I have introduced a parameter cg2dpcOffDFac which
C           defaults to 0.51 but can be set at runtime.
C     delR      :: Vertical grid spacing ( units of r ).          垂直网格间距（r单位）。
C     delRc     :: Vertical grid spacing between cell centers (r unit).  单元格中心之间的垂直网格间距（r单位）。
C     delX      :: Separation between cell faces (m) or (deg), depending
C     delY         on input flags. Note: moved to header file SET_GRID.h  单元面之间的间距（m）或（deg），取决于输入标志。注意：移到头文件SET_GRID.h 
C     xgOrigin   :: Origin of the X-axis (Cartesian Grid) / Longitude of Western   X轴原点（笛卡尔网格）/最南面纬度（纬度网格）
C                :: most cell face (Lat-Lon grid) (Note: this is an "inert"  （注意：这是一个“惰性”参数，但它使地理参考变得简单。）
C                :: parameter but it makes geographical references simple.)
C     ygOrigin   :: Origin of the Y-axis (Cartesian Grid) / Latitude of Southern   Y轴原点（笛卡尔网格）/最南面纬度（纬度网格）
C                :: most face (Lat-Lon grid).
C     rSphere    :: Radius of sphere for a spherical polar grid ( m ). 球面极坐标网的球面半径（m）
C     recip_rSphere :: Reciprocal radius of sphere ( m^-1 ).           球面倒数半径（m^-1）
C     radius_fromHorizGrid :: sphere Radius of input horiz. grid (Curvilinear Grid)   输入水平的球体半径。网格（曲线网格）
C     seaLev_Z   :: the reference height of sea-level (usually zero)   海平面的基准高度（通常为零）
C     top_Pres   :: pressure (P-Coords) or reference pressure (Z-Coords) at the top   顶部压力（P坐标）或参考压力（Z坐标）
C     rSigmaBnd  :: vertical position (in r-unit) of r/sigma transition (Hybrid-Sigma)   r/sigma转换（混合sigma）的垂直位置（以r为单位）
C     gravity    :: Acceleration due to constant gravity ( m/s^2 )   恒重力加速度（m/s^2）
C     recip_gravity :: Reciprocal gravity acceleration ( s^2/m )     重力反加速度（s^2/m）
C     gBaro      :: Accel. due to gravity used in barotropic equation ( m/s^2 )   加速。由于正压方程中使用的重力（m/s^2）
C     gravFacC   :: gravity factor (vs surf. gravity) vert. profile at cell-Center  重力系数（vs。重力）垂直。cell中心轮廓
C     gravFacF   :: gravity factor (vs surf. gravity) vert. profile at cell-interF  重力系数（vs。重力）垂直。cell间轮廓
C     rhoNil     :: Reference density for the linear equation of state  线性状态方程的参考密度
C     rhoConst   :: Vertically constant reference density (Boussinesq)  垂直恒定参考密度（Boussinesq）
C     rho1Ref    :: reference vertical profile for density (anelastic)  密度参考垂直剖面（滞弹性）
C     rhoFacC    :: normalized (by rhoConst) reference density at cell-Center  细胞中心标准化（rhoConst）参考密度
C     rhoFacF    :: normalized (by rhoConst) reference density at cell-interFace 单元界面归一化（rhoConst）参考密度
C     rhoConstFresh :: Constant reference density for fresh water (rain)  淡水（雨水）恒定参考密度
C     thetaConst :: Constant reference for potential temperature  位温恒定参考值
C     tRef       :: reference vertical profile for potential temperature  潜在温度的参考垂直剖面
C     sRef       :: reference vertical profile for salinity/specific humidity  盐度/比湿度的参考垂直剖面
C     pRef4EOS   :: reference pressure used in EOS (case selectP_inEOS_Zc=1)   EOS中使用的参考压力（case selectP\u inEOS\u Zc=1）
C     phiRef     :: reference potential (press/rho, geopot) profile (m^2/s^2)  参考电位（press/rho，geopot）剖面（m^2/s^2）
C     dBdrRef    :: vertical gradient of reference buoyancy  [(m/s/r)^2]:  参考浮力的垂直梯度[（m/s/r）^2]：z坐标：=N^2\u参考=Brunt-Vaissala频率[s^-2]p坐标：=-（d.alpha/dp）\u参考[（m^2.s/kg）^2]
C                :: z-coord: = N^2_ref = Brunt-Vaissala frequency [s^-2]
C                :: p-coord: = -(d.alpha/dp)_ref          [(m^2.s/kg)^2]
C     rVel2wUnit :: units conversion factor (Non-Hydrostatic code),  单位换算系数（非静水压代码）
C                :: from r-coordinate vertical velocity to vertical velocity [m/s]. 从r坐标的垂直速度到垂直速度[m/s]。
C                :: z-coord: = 1 ; p-coord: wSpeed [m/s] = rVel [Pa/s] * rVel2wUnit z坐标：=1；p坐标：W速度[m/s]=rVel[Pa/s]*rVel2wUnit
C     wUnit2rVel :: units conversion factor (Non-Hydrostatic code),  单位换算系数（非静水压代码），
C                :: from vertical velocity [m/s] to r-coordinate vertical velocity. 从垂直速度[m/s]到r坐标垂直速度。z坐标：=1；p坐标：rVel[Pa/s]=W速度[m/s]*WUnit2水平
C                :: z-coord: = 1 ; p-coord: rVel [Pa/s] = wSpeed [m/s] * wUnit2rVel
C     mass2rUnit :: units conversion factor (surface forcing),   单位换算系数（表面作用力）
C                :: from mass per unit area [kg/m2] to vertical r-coordinate unit.  从单位面积质量[kg/m2]到垂直r坐标单位。
C                :: z-coord: = 1/rhoConst ( [kg/m2] / rho = [m] ) ;   z-coord: = 1/rhoConst
C                :: p-coord: = gravity    ( [kg/m2] *  g = [Pa] ) ;   p-coord: = gravity   
C     rUnit2mass :: units conversion factor (surface forcing),
C                :: from vertical r-coordinate unit to mass per unit area [kg/m2].  单位换算系数（表面作用力），从垂直r坐标单位到单位面积质量[kg/m2]。
C                :: z-coord: = rhoConst  ( [m] * rho = [kg/m2] ) ;
C                :: p-coord: = 1/gravity ( [Pa] /  g = [kg/m2] ) ;
C     f0         :: Reference coriolis parameter ( 1/s )  参考科里奥利参数
C                   ( Southern edge f for beta plane )
C     beta       :: df/dy ( s^-1.m^-1 )    beta平面的beta
C     fPrime     :: Second Coriolis parameter ( 1/s ), related to Y-component  第二个科里奥利参数（1/s），与旋转的Y分量有关（参考值=2.Omega.Cos（Phi））
C                   of rotation (reference value = 2.Omega.Cos(Phi))
C     omega      :: Angular velocity ( rad/s )   角速度（rad/s）
C     rotationPeriod :: Rotation period (s) (= 2.pi/omega)   旋转周期（=2.pi/omega）
C     viscArNr   :: vertical profile of Eddy viscosity coeff.   垂直混合动量的涡粘系数垂直剖面（单位：r^2/s）
C                   for vertical mixing of momentum ( units of r^2/s )
C     viscAh     :: Eddy viscosity coeff. for mixing of  动量横向混合的涡粘性系数（m^2/s）
C                   momentum laterally ( m^2/s )
C     viscAhW    :: Eddy viscosity coeff. for mixing of vertical 涡流粘滞系数。对于横向垂直动量的混合，静水压模型没有影响，如果未设置（m^2/s），则默认为viscAhD；
C                   momentum laterally, no effect for hydrostatic 如果水平可变，则不使用viscAhD。使用粘度。
C                   model, defaults to viscAhD if unset ( m^2/s )
C                   Not used if variable horiz. viscosity is used.
C     viscA4     :: Biharmonic viscosity coeff. for mixing of   双调和粘滞系数。横向动量混合（m^4/s）
C                   momentum laterally ( m^4/s )
C     viscA4W    :: Biharmonic viscosity coeff. for mixing of vertical 双调和粘滞系数。对于横向垂直动量的混合，静水压模型没有影响，如果未设置（m^2/s），则默认为viscA4D；
C                   momentum laterally, no effect for hydrostatic 如果水平可变，则不使用viscA4D。使用粘度。
C                   model, defaults to viscA4D if unset ( m^2/s ) 
C                   Not used if variable horiz. viscosity is used.
C     viscAhD    :: Eddy viscosity coeff. for mixing of momentum laterally  涡流粘滞系数。横向动量混合（作用于散度部分）（m^2/s）
C                   (act on Divergence part) ( m^2/s )
C     viscAhZ    :: Eddy viscosity coeff. for mixing of momentum laterally  涡流粘滞系数。横向动量混合（作用于涡度部分）（m^2/s）
C                   (act on Vorticity  part) ( m^2/s )
C     viscA4D    :: Biharmonic viscosity coeff. for mixing of momentum laterally  双调和粘滞系数。横向动量混合（作用于发散部分）（m^4/s）
C                   (act on Divergence part) ( m^4/s )
C     viscA4Z    :: Biharmonic viscosity coeff. for mixing of momentum laterally  双调和粘滞系数。横向动量混合（作用于涡度部分）（m^4/s）
C                   (act on Vorticity  part) ( m^4/s )
C     smag3D_coeff :: Isotropic 3-D Smagorinsky coefficient (-)  各向同性三维Smagorinsky系数（-）
C     viscC2leith  :: Leith non-dimensional viscosity factor (grad(vort))   Leith无因次粘滞系数（grad（vort））
C     viscC2leithD :: Modified Leith non-dimensional visc. factor (grad(div))  改进的Leith无量纲visc。系数（grad（div））
C     viscC4leith  :: Leith non-dimensional viscosity factor (grad(vort))   Leith无因次粘滞系数（grad（vort））
C     viscC4leithD :: Modified Leith non-dimensional viscosity factor (grad(div))  修正的Leith无因次粘度因子（grad（div））
C     viscC2smag   :: Smagorinsky non-dimensional viscosity factor (harmonic)  Smagorinsky非维粘度系数（谐波）
C     viscC4smag   :: Smagorinsky non-dimensional viscosity factor (biharmonic)  Smagorinsky无因次粘滞系数（双调和）
C     viscAhMax    :: Maximum eddy viscosity coeff. for mixing of  最大涡粘性系数。横向动量混合（m^2/s）
C                    momentum laterally ( m^2/s )
C     viscAhReMax  :: Maximum gridscale Reynolds number for eddy viscosity  涡流粘度系数的最大网格雷诺数。用于横向混合动量（非尺寸）
C                     coeff. for mixing of momentum laterally (non-dim)
C     viscAhGrid   :: non-dimensional grid-size dependent viscosity   非尺寸网格尺寸相关粘度
C     viscAhGridMax:: maximum and minimum harmonic viscosity coefficients ...  无因次网格尺寸相关的最大和最小谐波粘滞系数。
C     viscAhGridMin::  in terms of non-dimensional grid-size dependent visc.  
C     viscA4Max    :: Maximum biharmonic viscosity coeff. for mixing of  最大双调和粘度系数。横向动量混合（m^4/s）
C                     momentum laterally ( m^4/s )
C     viscA4ReMax  :: Maximum Gridscale Reynolds number for   双调和粘性系数的最大网格尺度雷诺数。横向动量（不暗淡）
C                     biharmonic viscosity coeff. momentum laterally (non-dim)
C     viscA4Grid   :: non-dimensional grid-size dependent bi-harmonic viscosity  无因次网格尺寸相关的双调和粘性
C     viscA4GridMax:: maximum and minimum biharmonic viscosity coefficients ...
C     viscA4GridMin::  in terms of non-dimensional grid-size dependent viscosity  基于无因次网格尺寸相关粘度的最大和最小双调和粘度系数
C     diffKhT   :: Laplacian diffusion coeff. for mixing of  拉普拉斯扩散系数。横向混合热量（m^2/s）
C                 heat laterally ( m^2/s )
C     diffK4T   :: Biharmonic diffusion coeff. for mixing of  双调和扩散系数。横向混合热量（m^4/s）
C                 heat laterally ( m^4/s )
C     diffKrNrT :: vertical profile of Laplacian diffusion coeff.  垂直混合热的拉普拉斯扩散系数的垂直分布（r^2/s单位）
C                 for mixing of heat vertically ( units of r^2/s )
C     diffKr4T  :: vertical profile of Biharmonic diffusion coeff.  双调和扩散系数的垂直分布。垂直混合热量（r^4/s单位）
C                 for mixing of heat vertically ( units of r^4/s )
C     diffKhS  ::  Laplacian diffusion coeff. for mixing of  拉普拉斯扩散系数。用于横向混合盐（m^2/s）
C                 salt laterally ( m^2/s )
C     diffK4S   :: Biharmonic diffusion coeff. for mixing of  双调和扩散系数。用于横向混合盐（m^4/s）
C                 salt laterally ( m^4/s )
C     diffKrNrS :: vertical profile of Laplacian diffusion coeff.  拉普拉斯扩散系数的垂直剖面。用于垂直混合盐（r^2/s单位），
C                 for mixing of salt vertically ( units of r^2/s ),
C     diffKr4S  :: vertical profile of Biharmonic diffusion coeff.  双调和扩散系数的垂直分布。用于垂直混合盐（r^4/s单位）
C                 for mixing of salt vertically ( units of r^4/s )
C     diffKrBL79surf :: T/S surface diffusivity (m^2/s) Bryan and Lewis, 1979  T/S表面扩散率（m^2/S）Bryan和Lewis，1979年
C     diffKrBL79deep :: T/S deep diffusivity (m^2/s) Bryan and Lewis, 1979     T/S深扩散系数（m^2/S）Bryan和Lewis，1979年
C     diffKrBL79scl  :: depth scale for arctan fn (m) Bryan and Lewis, 1979    arctan fn（m）深度标度Bryan and Lewis，1979
C     diffKrBL79Ho   :: depth offset for arctan fn (m) Bryan and Lewis, 1979   arctan fn（m）深度偏移Bryan and Lewis，1979
C     BL79LatVary    :: polarwise of this latitude diffKrBL79 is applied with  这个纬度的极性diffKrBL79逐渐向赤道过渡到diffKrBLEQ
C                       gradual transition to diffKrBLEQ towards Equator
C     diffKrBLEQsurf :: same as diffKrBL79surf but at Equator 与diffKrBL79surf相同，但在赤道
C     diffKrBLEQdeep :: same as diffKrBL79deep but at Equator
C     diffKrBLEQscl  :: same as diffKrBL79scl but at Equator
C     diffKrBLEQHo   :: same as diffKrBL79Ho but at Equator
C     deltaT    :: Default timestep ( s )   默认时间步长
C     deltaTClock  :: Timestep used as model "clock". This determines the
C                    IO frequencies and is used in tagging output. It can
C                    be totally different to the dynamical time. Typically
C                    it will be the deep-water timestep for accelerated runs.
C                    Frequency of checkpointing and dumping of the model state
C                    are referenced to this clock. ( s )
用作模型“时钟”的时间步长。它确定IO频率并用于标记输出。它可以完全不同于动态时间。
通常情况下，这将是加速运行的深水时间步。检查点的频率和模型状态的转储都参照这个时钟秒）
C     deltaTMom    :: Timestep for momemtum equations ( s ) 动量方程的时间步长 
C     dTtracerLev  :: Timestep for tracer equations ( s ), function of level k   示踪剂方程的时间步长，k级函数
C     deltaTFreeSurf :: Timestep for free-surface equation ( s )   自由曲面方程的时间步长
C     freeSurfFac  :: Parameter to turn implicit free surface term on or off  用于打开或关闭隐式自由曲面项的参数
C                     freeSurFac = 1. uses implicit free surface  使用隐式自由曲面
C                     freeSurFac = 0. uses rigid lid              使用刚性盖
C     abEps        :: Adams-Bashforth-2 stabilizing weight     Adams-Bashforth-2稳定配重
C     alph_AB      :: Adams-Bashforth-3 primary factor         Adams-Bashforth-3主要因素
C     beta_AB      :: Adams-Bashforth-3 secondary factor       Adams-Bashforth-3次要因素
C     implicSurfPress :: parameter of the Crank-Nickelson time stepping :  曲柄-尼克尔森时间步进参数：表面压力梯度的隐式部分
C                     Implicit part of Surface Pressure Gradient ( 0-1 )
C     implicDiv2DFlow :: parameter of the Crank-Nickelson time stepping :  Crank-Nickelson时间步参数：正压气流辐散的隐式部分（0-1）
C                     Implicit part of barotropic flow Divergence ( 0-1 )
C     implicitNHPress :: parameter of the Crank-Nickelson time stepping :  曲柄-尼克尔森时间步进参数：非静水压力梯度隐式部分（0-1）
C                     Implicit part of Non-Hydrostatic Pressure Gradient ( 0-1 )  
C     hFacMin      :: Minimum fraction size of a cell (affects hFacC etc...)  单元格的最小分数大小（影响hFacC等…）
C     hFacMinDz    :: Minimum dimensional size of a cell (affects hFacC etc..., m)   单元的最小尺寸（影响hFacC等，m）
C     hFacMinDp    :: Minimum dimensional size of a cell (affects hFacC etc..., Pa)  单元的最小尺寸
C     hFacMinDr    :: Minimum dimensional size of a cell (-> hFacC etc..., r units)  单元的最小尺寸（hFacC etc..., r units）
C     hFacInf      :: Threshold (inf and sup) for fraction size of surface cell  控制消失和产生水平的表面单元分数大小的阈值
C     hFacSup          that control vanishing and creating levels
C     tauCD         :: CD scheme coupling timescale ( s )  CD方案耦合时间刻度
C     rCD           :: CD scheme normalised coupling parameter (= 1 - deltaT/tauCD)  CD方案归一化耦合参数（=1-deltaT/tauCD）
C     epsAB_CD      :: Adams-Bashforth-2 stabilizing weight used in CD scheme  Adams-Bashforth-2稳定配重在CD方案中的应用
C     baseTime      :: model base time (time origin) = time @ iteration zero   模型基准时间（时间原点）=时间@迭代零点
C     startTime     :: Starting time for this integration ( s ).  此集成的开始时间。
C     endTime       :: Ending time for this integration ( s ).    此集成的结束时间。
C     chkPtFreq     :: Frequency of rolling check pointing ( s ). 滚动检查点的频率。
C     pChkPtFreq    :: Frequency of permanent check pointing ( s ). 永久性检查点的频率。
C     dumpFreq      :: Frequency with which model state is written to 将模型状态写入后处理文件的频率。
C                      post-processing files ( s ).
C     diagFreq      :: Frequency with which model writes diagnostic output  模型写入中间量诊断输出的频率。
C                      of intermediate quantities.
C     afFacMom      :: Advection of momentum term tracer parameter  动量项示踪参数平流
C     vfFacMom      :: Momentum viscosity tracer parameter          动量粘度示踪参数
C     pfFacMom      :: Momentum pressure forcing tracer parameter   动量压力强迫示踪参数
C     cfFacMom      :: Coriolis term tracer parameter               科里奥利项示踪参数
C     foFacMom      :: Momentum forcing tracer parameter            动量强迫示踪参数          
C     mtFacMom      :: Metric terms tracer parameter                度量项跟踪参数
C     cosPower      :: Power of cosine of latitude to multiply viscosity  纬度余弦乘以粘度的幂
C     cAdjFreq      :: Frequency of convective adjustment           对流平差频率
C
C     taveFreq      :: Frequency with which time-averaged model state  将时间平均模型状态写入后处理文件的频率
C                      is written to post-processing files ( s ).
C     tave_lastIter :: (for state variable only) fraction of the last time    （仅适用于状态变量）最后一个时间步长（每个taveFreq周期）的
C                      step (of each taveFreq period) put in the time average. 分数放入时间平均值(第一个iter的分数=1-tave（最后一个iter）
C                      (fraction for 1rst iter = 1 - tave_lastIter)
C     tauThetaClimRelax :: Relaxation to climatology time scale ( s ).  松弛到气候学时间尺度
C     tauSaltClimRelax :: Relaxation to climatology time scale ( s ).   松弛到气候学时间尺度
C     latBandClimRelax :: latitude band where Relaxation to Clim. is applied,    latitude band松弛，在这个范围内的松弛一下
C                         i.e. where |yC| <= latBandClimRelax
C     externForcingPeriod :: Is the period of which forcing varies (eg. 1 month)   强迫变化的周期（如1个月）
C     externForcingCycle :: Is the repeat time of the forcing (eg. 1 year)         强迫的重复时间（如1年）
C                          (note: externForcingCycle must be an integer  （注意：externForcingCycle必须是externForcingPeriod的整数倍）
C                           number times externForcingPeriod)
C     convertFW2Salt :: salinity, used to convert Fresh-Water Flux to Salt Flux  盐度，用于将淡水通量转换为盐通量（如果设置为-1，则使用模型表面（局部）值）
C                       (use model surface (local) value if set to -1)
C     temp_EvPrRn :: temperature of Rain & Evap.  雨水和蒸发排放温度。
C     salt_EvPrRn :: salinity of Rain & Evap.  雨水和蒸发的盐分
C     temp_addMass :: temperature of addMass array  添加质量的温度
C     salt_addMass :: salinity of addMass array     添加质量的盐度
C        (notes: a) tracer content of Rain/Evap only used if both  只当NonLin_FrSurf & useRealFreshWater设置时，降水和蒸发的示踪物含量才使用
C                     NonLin_FrSurf & useRealFreshWater are set.
C                b) use model surface (local) value if set to UNSET_RL) 使用模型表面（局部）值（如果设置为未设置）
C     hMixCriteria:: criteria for mixed-layer diagnostic      混合层诊断标准
C     dRhoSmall   :: parameter for mixed-layer diagnostic     混合层诊断参数
C     hMixSmooth  :: Smoothing parameter for mixed-layer diag (default=0=no smoothing)   混合层diag的平滑参数（默认值=0=无平滑）
C     ivdc_kappa  :: implicit vertical diffusivity for convection [m^2/s]    对流的隐式垂直扩散率[m^2/s]
C     sideDragFactor     :: side-drag scaling factor (used only if no_slip_sides)  侧阻力比例因子（仅在无滑动侧时使用）（默认值=2：全阻力=1:半滑（BC）
C                           (default=2: full drag ; =1: gives half-slip BC)
C     bottomDragLinear    :: Linear    bottom-drag coefficient (units of [r]/s)  线性底部阻力系数（单位[r]/s）
C     bottomDragQuadratic :: Quadratic bottom-drag coefficient (units of [r]/m)  二次底阻力系数（单位[r]/m）
C               (if using zcoordinate, units becomes linear: m/s, quadratic: [-])（如果使用zcoordinate，单位变为线性：m/s，二次：[-]）
C     smoothAbsFuncRange :: 1/2 of interval around zero, for which FORTRAN ABS  0附近的1/2间隔，FORTRAN ABS将被一个更平滑的函数代替(affects myabs, mymin, mymax)
C                           is to be replace by a smoother function
C                           (affects myabs, mymin, mymax)
C     nh_Am2        :: scales the non-hydrostatic terms and changes internal scales   缩放非静水压项并更改内部比例
C                      (i.e. allows convection at different Rayleigh numbers)   （即允许不同瑞利数下的对流）
C     tCylIn        :: Temperature of the cylinder inner boundary  圆筒内边界温度
C     tCylOut       :: Temperature of the cylinder outer boundary  圆筒外边界温度
C     phiEuler      :: Euler angle, rotation about original z-axis 欧拉角，绕原始z轴旋转
C     thetaEuler    :: Euler angle, rotation about new x-axis      欧拉角，绕新x轴旋转
C     psiEuler      :: Euler angle, rotation about new z-axis      欧拉角，绕新z轴旋转
      COMMON /PARM_R/ cg2dTargetResidual, cg2dTargetResWunit,
     & cg2dpcOffDFac, cg3dTargetResidual,
     & delR, delRc, xgOrigin, ygOrigin, rSphere, recip_rSphere,
     & radius_fromHorizGrid, seaLev_Z, top_Pres, rSigmaBnd,
     & deltaT, deltaTMom, dTtracerLev, deltaTFreeSurf, deltaTClock,
     & abEps, alph_AB, beta_AB,
     & f0, beta, fPrime, omega, rotationPeriod,
     & viscFacAdj, viscAh, viscAhW, smag3D_coeff,
     & viscAhMax, viscAhGrid, viscAhGridMax, viscAhGridMin,
     & viscC2leith, viscC2leithD,
     & viscC2smag, viscC4smag,
     & viscAhD, viscAhZ, viscA4D, viscA4Z,
     & viscA4, viscA4W, viscA4Max,
     & viscA4Grid, viscA4GridMax, viscA4GridMin,
     & viscAhReMax, viscA4ReMax,
     & viscC4leith, viscC4leithD, viscArNr,
     & diffKhT, diffK4T, diffKrNrT, diffKr4T,
     & diffKhS, diffK4S, diffKrNrS, diffKr4S,
     & diffKrBL79surf, diffKrBL79deep, diffKrBL79scl, diffKrBL79Ho,
     & BL79LatVary,
     & diffKrBLEQsurf, diffKrBLEQdeep, diffKrBLEQscl, diffKrBLEQHo,
     & tauCD, rCD, epsAB_CD,
     & freeSurfFac, implicSurfPress, implicDiv2DFlow, implicitNHPress,
     & hFacMin, hFacMinDz, hFacInf, hFacSup,
     & gravity, recip_gravity, gBaro,
     & gravFacC, recip_gravFacC, gravFacF, recip_gravFacF,
     & rhoNil, rhoConst, recip_rhoConst, rho1Ref,
     & rhoFacC, recip_rhoFacC, rhoFacF, recip_rhoFacF, rhoConstFresh,
     & thetaConst, tRef, sRef, pRef4EOS, phiRef, dBdrRef,
     & rVel2wUnit, wUnit2rVel, mass2rUnit, rUnit2mass,
     & baseTime, startTime, endTime,
     & chkPtFreq, pChkPtFreq, dumpFreq, adjDumpFreq,
     & diagFreq, taveFreq, tave_lastIter, monitorFreq, adjMonitorFreq,
     & afFacMom, vfFacMom, pfFacMom, cfFacMom, foFacMom, mtFacMom,
     & cosPower, cAdjFreq,
     & tauThetaClimRelax, tauSaltClimRelax, latBandClimRelax,
     & externForcingCycle, externForcingPeriod,
     & convertFW2Salt, temp_EvPrRn, salt_EvPrRn,
     & temp_addMass, salt_addMass, hFacMinDr, hFacMinDp,
     & ivdc_kappa, hMixCriteria, dRhoSmall, hMixSmooth,
     & sideDragFactor, bottomDragLinear, bottomDragQuadratic, nh_Am2,
     & smoothAbsFuncRange,
     & tCylIn, tCylOut,
     & phiEuler, thetaEuler, psiEuler

      _RL cg2dTargetResidual
      _RL cg2dTargetResWunit
      _RL cg3dTargetResidual
      _RL cg2dpcOffDFac
      _RL delR(Nr)
      _RL delRc(Nr+1)
      _RL xgOrigin
      _RL ygOrigin
      _RL rSphere
      _RL recip_rSphere
      _RL radius_fromHorizGrid
      _RL seaLev_Z
      _RL top_Pres
      _RL rSigmaBnd
      _RL deltaT
      _RL deltaTClock
      _RL deltaTMom
      _RL dTtracerLev(Nr)
      _RL deltaTFreeSurf
      _RL abEps, alph_AB, beta_AB
      _RL f0
      _RL beta
      _RL fPrime
      _RL omega
      _RL rotationPeriod
      _RL freeSurfFac
      _RL implicSurfPress
      _RL implicDiv2DFlow
      _RL implicitNHPress
      _RL hFacMin
      _RL hFacMinDz
      _RL hFacMinDp
      _RL hFacMinDr
      _RL hFacInf
      _RL hFacSup
      _RL viscArNr(Nr)
      _RL viscFacAdj
      _RL viscAh
      _RL viscAhW
      _RL viscAhD
      _RL viscAhZ
      _RL smag3D_coeff
      _RL viscAhMax
      _RL viscAhReMax
      _RL viscAhGrid, viscAhGridMax, viscAhGridMin
      _RL viscC2leith
      _RL viscC2leithD
      _RL viscC2smag
      _RL viscA4
      _RL viscA4W
      _RL viscA4D
      _RL viscA4Z
      _RL viscA4Max
      _RL viscA4ReMax
      _RL viscA4Grid, viscA4GridMax, viscA4GridMin
      _RL viscC4leith
      _RL viscC4leithD
      _RL viscC4smag
      _RL diffKhT
      _RL diffK4T
      _RL diffKrNrT(Nr)
      _RL diffKr4T(Nr)
      _RL diffKhS
      _RL diffK4S
      _RL diffKrNrS(Nr)
      _RL diffKr4S(Nr)
      _RL diffKrBL79surf
      _RL diffKrBL79deep
      _RL diffKrBL79scl
      _RL diffKrBL79Ho
      _RL BL79LatVary
      _RL diffKrBLEQsurf
      _RL diffKrBLEQdeep
      _RL diffKrBLEQscl
      _RL diffKrBLEQHo
      _RL tauCD, rCD, epsAB_CD
      _RL gravity,       recip_gravity
      _RL gBaro
      _RL gravFacC(Nr),   recip_gravFacC(Nr)
      _RL gravFacF(Nr+1), recip_gravFacF(Nr+1)
      _RL rhoNil
      _RL rhoConst,      recip_rhoConst
      _RL rho1Ref(Nr)
      _RL rhoFacC(Nr),   recip_rhoFacC(Nr)
      _RL rhoFacF(Nr+1), recip_rhoFacF(Nr+1)
      _RL rhoConstFresh
      _RL thetaConst
      _RL tRef(Nr)
      _RL sRef(Nr)
      _RL pRef4EOS(Nr)
      _RL phiRef(2*Nr+1)
      _RL dBdrRef(Nr)
      _RL rVel2wUnit(Nr+1), wUnit2rVel(Nr+1)
      _RL mass2rUnit, rUnit2mass
      _RL baseTime
      _RL startTime
      _RL endTime
      _RL chkPtFreq
      _RL pChkPtFreq
      _RL dumpFreq
      _RL adjDumpFreq
      _RL diagFreq
      _RL taveFreq
      _RL tave_lastIter
      _RL monitorFreq
      _RL adjMonitorFreq
      _RL afFacMom
      _RL vfFacMom
      _RL pfFacMom
      _RL cfFacMom
      _RL foFacMom
      _RL mtFacMom
      _RL cosPower
      _RL cAdjFreq
      _RL tauThetaClimRelax
      _RL tauSaltClimRelax
      _RL latBandClimRelax
      _RL externForcingCycle
      _RL externForcingPeriod
      _RL convertFW2Salt
      _RL temp_EvPrRn
      _RL salt_EvPrRn
      _RL temp_addMass
      _RL salt_addMass
      _RL ivdc_kappa
      _RL hMixCriteria
      _RL dRhoSmall
      _RL hMixSmooth
      _RL sideDragFactor
      _RL bottomDragLinear
      _RL bottomDragQuadratic
      _RL smoothAbsFuncRange
      _RL nh_Am2
      _RL tCylIn, tCylOut
      _RL phiEuler, thetaEuler, psiEuler

C--   COMMON /PARM_A/ Thermodynamics constants ?   热力学常数
      COMMON /PARM_A/ HeatCapacity_Cp    热容量
      _RL HeatCapacity_Cp

C--   COMMON /PARM_ATM/ Atmospheric physical parameters (Ideal Gas EOS, ...)    大气物理参数（理想气体状态方程）
C     celsius2K :: convert centigrade (Celsius) degree to Kelvin  将摄氏度转换成开尔文
C     atm_Po    :: standard reference pressure                    标准参考压力
C     atm_Cp    :: specific heat (Cp) of the (dry) air at constant pressure   恒压下（干）空气的比热
C     atm_Rd    :: gas constant for dry air                       干空气气体常数
C     atm_kappa :: kappa = R/Cp (R: constant of Ideal Gas EOS)    kappa=R/Cp（R：理想气体状态方程常数）
C     atm_Rq    :: water vapour specific volume anomaly relative to dry air   相对于干空气的水汽比容异常
C                  (e.g. typical value = (29/18 -1) 10^-3 with q [g/kg])   （例如，典型值=（29/18-1）10^-3和q[g/kg]）
C     integr_GeoPot :: option to select the way we integrate the geopotential  选择整合位势的方式
C                     (still a subject of discussions ...) （仍在讨论中……）
C     selectFindRoSurf :: select the way surf. ref. pressure (=Ro_surf) is        Ro_surf
C             derived from the orography. Implemented: 0,1 (see INI_P_GROUND)
      COMMON /PARM_ATM/
     &            celsius2K,
     &            atm_Cp, atm_Rd, atm_kappa, atm_Rq, atm_Po,
     &            integr_GeoPot, selectFindRoSurf
      _RL celsius2K
      _RL atm_Po, atm_Cp, atm_Rd, atm_kappa, atm_Rq
      INTEGER integr_GeoPot, selectFindRoSurf

C Logical flags for selecting packages      用于选择包的逻辑标志
      LOGICAL useGAD
      LOGICAL useOBCS
      LOGICAL useSHAP_FILT
      LOGICAL useZONAL_FILT
      LOGICAL useOPPS
      LOGICAL usePP81
      LOGICAL useKL10
      LOGICAL useMY82
      LOGICAL useGGL90
      LOGICAL useKPP
      LOGICAL useGMRedi
      LOGICAL useDOWN_SLOPE
      LOGICAL useBBL
      LOGICAL useCAL
      LOGICAL useEXF
      LOGICAL useBulkForce
      LOGICAL useEBM
      LOGICAL useCheapAML
      LOGICAL useAUTODIFF
      LOGICAL useGrdchk
      LOGICAL useSMOOTH
      LOGICAL usePROFILES
      LOGICAL useECCO
      LOGICAL useCTRL
      LOGICAL useSBO
      LOGICAL useFLT
      LOGICAL usePTRACERS
      LOGICAL useGCHEM
      LOGICAL useRBCS
      LOGICAL useOffLine
      LOGICAL useMATRIX
      LOGICAL useFRAZIL
      LOGICAL useSEAICE
      LOGICAL useSALT_PLUME
      LOGICAL useShelfIce
      LOGICAL useStreamIce
      LOGICAL useICEFRONT
      LOGICAL useThSIce
      LOGICAL useLand
      LOGICAL useATM2d
      LOGICAL useAIM
      LOGICAL useAtm_Phys
      LOGICAL useFizhi
      LOGICAL useGridAlt
      LOGICAL useDiagnostics
      LOGICAL useREGRID
      LOGICAL useLayers
      LOGICAL useMNC
      LOGICAL useRunClock
      LOGICAL useEMBED_FILES
      LOGICAL useMYPACKAGE
      COMMON /PARM_PACKAGES/
     &        useGAD, useOBCS, useSHAP_FILT, useZONAL_FILT,
     &        useOPPS, usePP81, useKL10, useMY82, useGGL90, useKPP,
     &        useGMRedi, useBBL, useDOWN_SLOPE,
     &        useCAL, useEXF, useBulkForce, useEBM, useCheapAML,
     &        useGrdchk, useSMOOTH, usePROFILES, useECCO, useCTRL,
     &        useSBO, useFLT, useAUTODIFF,
     &        usePTRACERS, useGCHEM, useRBCS, useOffLine, useMATRIX,
     &        useFRAZIL, useSEAICE, useSALT_PLUME, useShelfIce,
     &        useStreamIce, useICEFRONT, useThSIce, useLand,
     &        useATM2D, useAIM, useAtm_Phys, useFizhi, useGridAlt,
     &        useDiagnostics, useREGRID, useLayers, useMNC,
     &        useRunClock, useEMBED_FILES,
     &        useMYPACKAGE

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
