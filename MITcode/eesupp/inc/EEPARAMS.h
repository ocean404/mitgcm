C $Header: /u/gcmpack/MITgcm/eesupp/inc/EEPARAMS.h,v 1.36 2014/02/13 04:10:25 jmc Exp $
C $Name:  $
CBOP
C     !ROUTINE: EEPARAMS.h
C     !INTERFACE:
C     include "EEPARAMS.h"
C
C     !DESCRIPTION:
C     *==========================================================*
C     | EEPARAMS.h                                               |
C     *==========================================================*
C     | Parameters for "execution environemnt". These are used   | “执行环境”的参数。这些都由特定的数值模型和执行环境支持例程使用。 
C     | by both the particular numerical model and the execution |
C     | environment support routines.                            |
C     *==========================================================*
CEOP

C     ========  EESIZE.h  ========================================

C     MAX_LEN_MBUF  :: Default message buffer max. size 默认消息缓冲区最大大小  
C     MAX_LEN_FNAM  :: Default file name max. size 文件名最大值
C     MAX_LEN_PREC  :: Default rec len for reading "parameter" files 读取参数文件时最大的长度

      INTEGER MAX_LEN_MBUF
      PARAMETER ( MAX_LEN_MBUF = 512 )
      INTEGER MAX_LEN_FNAM
      PARAMETER ( MAX_LEN_FNAM = 512 )
      INTEGER MAX_LEN_PREC
      PARAMETER ( MAX_LEN_PREC = 200 )

C     MAX_NO_THREADS  :: Maximum number of threads allowed.  允许的最大线程数
CC    MAX_NO_PROCS    :: Maximum number of processes allowed.允许的最大进程数
CC    MAX_NO_BARRIERS :: Maximum number of distinct thread "barriers"不同线程“屏障”的最大数目
      INTEGER MAX_NO_THREADS
      PARAMETER ( MAX_NO_THREADS =  4 )
c     INTEGER MAX_NO_PROCS
c     PARAMETER ( MAX_NO_PROCS   =  70000 )
c     INTEGER MAX_NO_BARRIERS
c     PARAMETER ( MAX_NO_BARRIERS = 1 )

C     Particularly weird and obscure voodoo numbers  特别奇怪和晦涩的巫毒数字
C     lShare :: This wants to be the length in
C               [148]-byte words of the size of
C               the address "window" that is snooped
C               on an SMP bus. By separating elements in
C               the global sum buffer we can avoid generating
C               extraneous invalidate traffic between
C               processors. The length of this window is usually
C               a cache line i.e. small O(64 bytes).
C               The buffer arrays are usually short arrays
C               and are declared REAL ARRA(lShare[148],LBUFF).
C               Setting lShare[148] to 1 is like making these arrays
C               one dimensional.
      INTEGER cacheLineSize
      INTEGER lShare1
      INTEGER lShare4
      INTEGER lShare8
      PARAMETER ( cacheLineSize = 256 )
      PARAMETER ( lShare1 =  cacheLineSize )
      PARAMETER ( lShare4 =  cacheLineSize/4 )
      PARAMETER ( lShare8 =  cacheLineSize/8 )

CC    MAX_VGS  :: Maximum buffer size for Global Vector Sum 全局向量和的最大缓冲区大小
c     INTEGER MAX_VGS
c     PARAMETER ( MAX_VGS = 8192 )

C     ========  EESIZE.h  ========================================

C     Symbolic values 符号值
C     precXXXX :: precision used for I/O   使用的精度
      INTEGER precFloat32
      PARAMETER ( precFloat32 = 32 )
      INTEGER precFloat64
      PARAMETER ( precFloat64 = 64 )

C     Real-type constant for some frequently used simple number (0,1,2,1/2):   一些常用简单数的实型常数
      _RS     zeroRS, oneRS, twoRS, halfRS
      PARAMETER ( zeroRS = 0.0 _d 0 , oneRS  = 1.0 _d 0 )
      PARAMETER ( twoRS  = 2.0 _d 0 , halfRS = 0.5 _d 0 )
      _RL     zeroRL, oneRL, twoRL, halfRL
      PARAMETER ( zeroRL = 0.0 _d 0 , oneRL  = 1.0 _d 0 )
      PARAMETER ( twoRL  = 2.0 _d 0 , halfRL = 0.5 _d 0 )

C     UNSET_xxx :: Used to indicate variables that have not been given a value   用于指示尚未给定值的变量
      Real*8  UNSET_FLOAT8
      PARAMETER ( UNSET_FLOAT8 = 1.234567D5 )         D是双精度型数，E是实型数
      Real*4  UNSET_FLOAT4
      PARAMETER ( UNSET_FLOAT4 = 1.234567E5 )
      _RL     UNSET_RL
      PARAMETER ( UNSET_RL     = 1.234567D5 )
      _RS     UNSET_RS
      PARAMETER ( UNSET_RS     = 1.234567D5 )
      INTEGER UNSET_I
      PARAMETER ( UNSET_I      = 123456789  )

C     debLevX  :: used to decide when to print debug messages  用于决定何时打印调试消息
      INTEGER debLevZero
      INTEGER debLevA, debLevB,  debLevC, debLevD, debLevE
      PARAMETER ( debLevZero=0 )
      PARAMETER ( debLevA=1 )
      PARAMETER ( debLevB=2 )
      PARAMETER ( debLevC=3 )
      PARAMETER ( debLevD=4 )
      PARAMETER ( debLevE=5 )

C     SQUEEZE_RIGHT      :: Flag indicating right blank space removal
C                           from text field.  指示从文本字段中删除右空格的标志。
C     SQUEEZE_LEFT       :: Flag indicating left blank space removal
C                           from text field.  指示从文本字段中删除左空格的标志。
C     SQUEEZE_BOTH       :: Flag indicating left and right blank
C                           space removal from text field.  指示从文本字段中删除左右空格的标志。
C     PRINT_MAP_XY       :: Flag indicating to plot map as XY slices 指示将地图打印为XY切片的标志
C     PRINT_MAP_XZ       :: Flag indicating to plot map as XZ slices 指示将地图打印为XZ切片的标志
C     PRINT_MAP_YZ       :: Flag indicating to plot map as YZ slices 指示将地图打印为YZ切片的标志
C     commentCharacter   :: Variable used in column 1 of parameter
C                           files to indicate comments. 在参数文件的第1列中用于指示注释的变量。
C     INDEX_I            :: Variable used to select an index label
C     INDEX_J               for formatted input parameters.  用于为格式化输入参数选择索引标签的变量。
C     INDEX_K
C     INDEX_NONE
      CHARACTER*(*) SQUEEZE_RIGHT
      PARAMETER ( SQUEEZE_RIGHT = 'R' )
      CHARACTER*(*) SQUEEZE_LEFT
      PARAMETER ( SQUEEZE_LEFT = 'L' )
      CHARACTER*(*) SQUEEZE_BOTH
      PARAMETER ( SQUEEZE_BOTH = 'B' )
      CHARACTER*(*) PRINT_MAP_XY
      PARAMETER ( PRINT_MAP_XY = 'XY' )
      CHARACTER*(*) PRINT_MAP_XZ
      PARAMETER ( PRINT_MAP_XZ = 'XZ' )
      CHARACTER*(*) PRINT_MAP_YZ
      PARAMETER ( PRINT_MAP_YZ = 'YZ' )
      CHARACTER*(*) commentCharacter
      PARAMETER ( commentCharacter = '#' )
      INTEGER INDEX_I
      INTEGER INDEX_J
      INTEGER INDEX_K
      INTEGER INDEX_NONE
      PARAMETER ( INDEX_I    = 1,
     &            INDEX_J    = 2,
     &            INDEX_K    = 3,
     &            INDEX_NONE = 4 )

C     EXCH_IGNORE_CORNERS :: Flag to select ignoring or
C     EXCH_UPDATE_CORNERS    updating of corners during an edge exchange.  在边交换期间选择忽略或更新角点的标志。
      INTEGER EXCH_IGNORE_CORNERS
      INTEGER EXCH_UPDATE_CORNERS
      PARAMETER ( EXCH_IGNORE_CORNERS = 0,
     &            EXCH_UPDATE_CORNERS = 1 )

C     FORWARD_SIMULATION   正向模拟
C     REVERSE_SIMULATION   反向模拟
C     TANGENT_SIMULATION   切向模拟
      INTEGER FORWARD_SIMULATION
      INTEGER REVERSE_SIMULATION
      INTEGER TANGENT_SIMULATION
      PARAMETER ( FORWARD_SIMULATION = 0,
     &            REVERSE_SIMULATION = 1,
     &            TANGENT_SIMULATION = 2 )

C--   COMMON /EEPARAMS_L/ Execution environment public logical variables.  执行环境公共逻辑变量。
C     eeBootError    :: Flags indicating error during multi-processing  多处理期间指示错误的标志
C     eeEndError     :: initialisation and termination.  初始化和终止
C     fatalError     :: Flag used to indicate that the model is ended with an error  用于指示模型以错误结束的标志
C     debugMode      :: controls printing of debug msg (sequence of S/R calls).  控制调试消息（S/R调用序列）的打印。
C     useSingleCpuIO :: When useSingleCpuIO is set, MDS_WRITE_FIELD outputs from
C                       master MPI process only. -- NOTE: read from main parameter
C                       file "data" and not set until call to INI_PARMS. 当设置useSingleCpuIO时，MDS\u WRITE\u字段仅从主MPI进程输出。-注意：从主参数文件“data”读取，直到调用INI参数才设置。
C     useSingleCpuInput :: When useSingleCpuInput is set, EXF_INTERP_READ
C                       reads forcing files from master MPI process only.当设置useSingleCpuInput时，EXF\u INTERP\u READ只从主MPI进程读取强制文件。
C                       -- NOTE: read from main parameter file "data"
C                          and defaults to useSingleCpuInput = useSingleCpuIO
C     printMapIncludesZeros  :: Flag that controls whether character constant
C                               map code ignores exact zero values. 控制字符常量映射代码是否忽略精确零值的标志。
C     useCubedSphereExchange :: use Cubed-Sphere topology domain.   使用立方体球体拓扑域。
C     useCoupler     :: use Coupler for a multi-components set-up.  使用耦合器进行多组件设置
C     useNEST_PARENT :: use Parent Nesting interface (pkg/nest_parent)  使用父嵌套接口
C     useNEST_CHILD  :: use Child  Nesting interface (pkg/nest_child)   使用子嵌套接口
C     useOASIS       :: use OASIS-coupler for a multi-components set-up.使用OASIS耦合器进行多组件设置。
      COMMON /EEPARAMS_L/
c    &  eeBootError, fatalError, eeEndError,
     &  eeBootError, eeEndError, fatalError, debugMode,
     &  useSingleCpuIO, useSingleCpuInput, printMapIncludesZeros,
     &  useCubedSphereExchange, useCoupler,
     &  useNEST_PARENT, useNEST_CHILD, useOASIS,
     &  useSETRLSTK, useSIGREG
      LOGICAL eeBootError
      LOGICAL eeEndError
      LOGICAL fatalError
      LOGICAL debugMode
      LOGICAL useSingleCpuIO
      LOGICAL useSingleCpuInput
      LOGICAL printMapIncludesZeros
      LOGICAL useCubedSphereExchange
      LOGICAL useCoupler
      LOGICAL useNEST_PARENT
      LOGICAL useNEST_CHILD
      LOGICAL useOASIS
      LOGICAL useSETRLSTK
      LOGICAL useSIGREG

C--   COMMON /EPARAMS_I/ Execution environment public integer variables.  执行环境公共整数变量
C     errorMessageUnit    :: Fortran IO unit for error messages  错误信息的Fortran IO unit 
C     standardMessageUnit :: Fortran IO unit for informational messages  信息性消息的Fortran IO unit
C     maxLengthPrt1D :: maximum length for printing (to Std-Msg-Unit) 1-D array  打印（至标准消息单元）1-D阵列的最大长度
C     scrUnit1      :: Scratch file 1 unit number  暂存文件1单元号
C     scrUnit2      :: Scratch file 2 unit number  暂存文件2单元号
C     eeDataUnit    :: Unit # for reading "execution environment" parameter file 单元#用于读取“执行环境”参数文件
C     modelDataUnit :: Unit number for reading "model" parameter file.           读取“模型”参数文件的单元号。
C     numberOfProcs :: Number of processes computing in parallel                 并行计算的进程数
C     pidIO         :: Id of process to use for I/O.                             用于I/O的进程Id
C     myBxLo, myBxHi :: Extents of domain in blocks in X and Y                   每个线程负责的X和Y中块中的域范围。    
C     myByLo, myByHi :: that each threads is responsble for.                     
C     myProcId      :: My own "process" id.                                      我自己的“进程”id
C     myPx          :: My X coord on the proc. grid.                             我在程序网格上的X坐标
C     myPy          :: My Y coord on the proc. grid.                                我在程序网格上的Y坐标
C     myXGlobalLo   :: My bottom-left (south-west) x-index global domain.   我的左下角（西南）的x索引值在全局域中。此处并不是指定该点的x坐标，例如m或度数。
C                      The x-coordinate of this point in for example m or   如果需要的话，模型需要提供一种推断信息的机制。
C                      degrees is *not* specified here. A model needs to
C                      provide a mechanism for deducing that information
C                      if it is needed.
C     myYGlobalLo   :: My bottom-left (south-west) y-index in global domain.同上
C                      The y-coordinate of this point in for example m or
C                      degrees is *not* specified here. A model needs to
C                      provide a mechanism for deducing that information
C                      if it is needed.
C     nThreads      :: No. of threads                                          线程数   
C     nTx, nTy      :: No. of threads in X and in Y                            X和Y的线程数
C                      This assumes a simple cartesian gridding of the threads  这假设了线程的简单笛卡尔网格，这在其他地方不需要，但这使得它更容易
C                      which is not required elsewhere but that makes it easier
C     ioErrorCount  :: IO Error Counter. Set to zero initially and increased   IO错误计数器。最初设置为零，每次发生IO错误时增加1。
C                      by one every time an IO error occurs.
      COMMON /EEPARAMS_I/
     &  errorMessageUnit, standardMessageUnit, maxLengthPrt1D,
     &  scrUnit1, scrUnit2, eeDataUnit, modelDataUnit,
     &  numberOfProcs, pidIO, myProcId,
     &  myPx, myPy, myXGlobalLo, myYGlobalLo, nThreads,
     &  myBxLo, myBxHi, myByLo, myByHi,
     &  nTx, nTy, ioErrorCount
      INTEGER errorMessageUnit
      INTEGER standardMessageUnit
      INTEGER maxLengthPrt1D
      INTEGER scrUnit1
      INTEGER scrUnit2
      INTEGER eeDataUnit
      INTEGER modelDataUnit
      INTEGER ioErrorCount(MAX_NO_THREADS)             整型向量，共有max个整型，目前都是0还没赋值
      INTEGER myBxLo(MAX_NO_THREADS)
      INTEGER myBxHi(MAX_NO_THREADS)
      INTEGER myByLo(MAX_NO_THREADS)
      INTEGER myByHi(MAX_NO_THREADS)
      INTEGER myProcId
      INTEGER myPx
      INTEGER myPy
      INTEGER myXGlobalLo
      INTEGER myYGlobalLo
      INTEGER nThreads
      INTEGER nTx
      INTEGER nTy
      INTEGER numberOfProcs
      INTEGER pidIO

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
