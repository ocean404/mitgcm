C $Header: /u/gcmpack/MITgcm/eesupp/inc/CPP_EEOPTIONS.h,v 1.41 2015/11/04 20:49:37 jmc Exp $
C $Name:  $

CBOP
C     !ROUTINE: CPP_EEOPTIONS.h
C     !INTERFACE:
C     include "CPP_EEOPTIONS.h"
C
C     !DESCRIPTION:
C     *==========================================================*
C     | CPP\_EEOPTIONS.h                                         |
C     *==========================================================*
C     | C preprocessor "execution environment" supporting        |许多选项同时使用编译时和运行时开关实现。这允许选项被完全删除，
C     | flags. Use this file to set flags controlling the        |在运行时成为可选的，或者被永久启用。此约定有助于伴随模型编译器执行的数据依赖性分析。这种数据依赖性分析可能会被运行时开关打乱，因为它无法将这些开关重新编码为在集成期间是固定的。
C     | execution environment in which a model runs - as opposed |使用这些标志的合理方法是在运行时将所有选项设置为可选，但是一旦确定了一个实验配置，就用在编译时设置的适当选项重新生成代码。
C     | to the dynamical problem the model solves.               |
C     | Note: Many options are implemented with both compile time|
C     |       and run-time switches. This allows options to be   |
C     |       removed altogether, made optional at run-time or   |
C     |       to be permanently enabled. This convention helps   |
C     |       with the data-dependence analysis performed by the |
C     |       adjoint model compiler. This data dependency       |
C     |       analysis can be upset by runtime switches that it  |
C     |       is unable to recoginise as being fixed for the     |
C     |       duration of an integration.                        |
C     |       A reasonable way to use these flags is to          |
C     |       set all options as selectable at runtime but then  |
C     |       once an experimental configuration has been        |
C     |       identified, rebuild the code with the appropriate  |
C     |       options set at compile time.                       |
C     *==========================================================*
CEOP

#ifndef _CPP_EEOPTIONS_H_
#define _CPP_EEOPTIONS_H_

C     In general the following convention applies:
C     ALLOW  - indicates an feature will be included but it may 指示将包含某个功能，但它可能有一个运行时标志，以允许打开和关闭该功能
C     CAN      have a run-time flag to allow it to be switched
C              on and off.
C              If ALLOW or CAN directives are "undef'd" this generally  如果ALLOW或CAN指令是“undef'd”，这通常意味着该功能将不可用，即它不会包含在编译代码中，因此没有使用该功能的运行时选项可用。
C              means that the feature will not be available i.e. it
C              will not be included in the compiled code and so no
C              run-time option to use the feature will be available.
C
C     ALWAYS - indicates the choice will be fixed at compile time  指示选项将在编译时修复，因此不存在运行时选项
C              so no run-time option will be present

C=== Macro related options ===
C--   Control storage of floating point operands
C     On many systems it improves performance only to use
C     8-byte precision for time stepped variables.
C     Constant in time terms ( geometric factors etc.. )
C     can use 4-byte precision, reducing memory utilisation and
C     boosting performance because of a smaller working set size.
C     However, on vector CRAY systems this degrades performance.
C     Enable to switch REAL4_IS_SLOW from genmake2 (with LET_RS_BE_REAL4):
#ifdef LET_RS_BE_REAL4
#undef REAL4_IS_SLOW
#else /* LET_RS_BE_REAL4 */
#define REAL4_IS_SLOW
#endif /* LET_RS_BE_REAL4 */

C--   Control use of "double" precision constants.
C     Use D0 where it means REAL*8 but not where it means REAL*16
#define D0 d0

C--   Enable some old macro conventions for backward compatibility
#undef USE_OLD_MACROS_R4R8toRSRL

C=== IO related options ===
C--   Flag used to indicate whether Fortran formatted write
C     and read are threadsafe. On SGI the routines can be thread
C     safe, on Sun it is not possible - if you are unsure then
C     undef this option.
#undef FMTFTN_IO_THREAD_SAFE

C--   Flag used to indicate whether Binary write to Local file (i.e.,  每个tile都有不同文件
C     a different file for each tile) and read are thread-safe.
#undef LOCBIN_IO_THREAD_SAFE

C--   Flag to turn off the writing of error message to ioUnit zero
#undef DISABLE_WRITE_TO_UNIT_ZERO

C--   Alternative formulation of BYTESWAP, faster than
C     compiler flag -byteswapio on the Altix.
#undef FAST_BYTESWAP

C--   Flag defined for eeboot_minimal.F, eeset_parms.F and open_copy_data_file.F仅从进程0写入STDOUT、STDERR和scratch文件
C     to write STDOUT, STDERR and scratch files from process 0 only.
C WARNING: to use only when absolutely confident that the setup is working
C     since any message (error/warning/print) from any proc <> 0 will be lost.
#undef SINGLE_DISK_IO

C=== MPI, EXCH and GLOBAL_SUM related options ===
C--   Flag turns off MPI_SEND ready_to_receive polling in the
C     gather_* subroutines to speed up integrations.
#undef DISABLE_MPI_READY_TO_RECEIVE

C--   Control MPI based parallel processing基于MPI的控制并行处理
CXXX We no longer select the use of MPI via this file (CPP_EEOPTIONS.h)
CXXX To use MPI, use an appropriate genmake2 options file or use
CXXX genmake2 -mpi .
CXXX #undef  ALLOW_USE_MPI

C--   Control use of communication that might overlap computation.
C     Under MPI selects/deselects "non-blocking" sends and receives.
#define ALLOW_ASYNC_COMMUNICATION
#undef  ALLOW_ASYNC_COMMUNICATION
#undef  ALWAYS_USE_ASYNC_COMMUNICATION
C--   Control use of communication that is atomic to computation.
C     Under MPI selects/deselects "blocking" sends and receives.
#define ALLOW_SYNC_COMMUNICATION
#undef  ALWAYS_USE_SYNC_COMMUNICATION

C--   Control XY periodicity in processor to grid mappings
C     Note: Model code does not need to know whether a domain is
C           periodic because it has overlap regions for every box.
C           Model assume that these values have been
C           filled in some way.
#undef  ALWAYS_PREVENT_X_PERIODICITY
#undef  ALWAYS_PREVENT_Y_PERIODICITY
#define CAN_PREVENT_X_PERIODICITY
#define CAN_PREVENT_Y_PERIODICITY

C--   disconnect tiles (no exchange between tiles, just fill-in edges断开分片连接
C     assuming locally periodic subdomain)
#undef DISCONNECTED_TILES

C--   Always cumulate tile local-sum in the same order by applying MPI allreduce
C     to array of tiles ; can get slower with large number of tiles (big set-up)
#define GLOBAL_SUM_ORDER_TILES

C--   Alternative way of doing global sum without MPI allreduce call
C     but instead, explicit MPI send & recv calls. Expected to be slower.
#undef GLOBAL_SUM_SEND_RECV

C--   Alternative way of doing global sum on a single CPU
C     to eliminate tiling-dependent roundoff errors. Note: This is slow.
#undef  CG2D_SINGLECPU_SUM

C=== Other options (to add/remove pieces of code) ===
C--   Flag to turn on checking for errors from all threads and procs
C     (calling S/R STOP_IF_ERROR) before stopping.在停止前打开检查所有线程和进程的错误（如果出错则调用S/R STOP_IF_ERROR）的标志。
#define USE_ERROR_STOP

C--   Control use of communication with other component:
C     allow to import and export from/to Coupler interface. 控制与其他组件通信的使用：**允许从耦合器接口导入和导出**。
#undef COMPONENT_MODULE

#endif /* _CPP_EEOPTIONS_H_ */

#include "CPP_EEMACROS.h"

