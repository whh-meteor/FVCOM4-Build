#===========================================================================
# FVCOM v4.4.10 make.inc - 终极适配 Windows/GitHub Actions
# 解决了路径映射、MPI 模块识别、预处理器兼容性及 fjulian.inc 包含问题
#===========================================================================

#========== 路径自动转换 (关键修正) =========================================
# 使用 cygpath 将 MSYS2 的虚拟路径转换为 gfortran 能识别的 Windows 物理路径
MSYS_ROOT     := $(shell cygpath -m /mingw64)

# 修复 TOPDIR 的物理路径识别
TOPDIR_W      := $(shell cygpath -m $(shell pwd)/..)
TOPDIR         = ..
INSTALLDIR     = $(TOPDIR)/libs/install

LIBPATH        = $(MSYS_ROOT)/lib
INCLUDEPATH    = $(MSYS_ROOT)/include

LIBDIR         = -L$(LIBPATH)
INCDIR         = -I$(INCLUDEPATH)

#---------- 内部库路径补丁 (解决 fjulian.inc 报错关键) ----------
# 放弃使用相对路径 ../libs/julian，直接强制使用物理绝对路径
JULIAN_INC     = -I$(TOPDIR_W)/libs/julian

#===========================================================================
# 预处理器选项 (CPP)
#===========================================================================
CPP            = gcc -E
DEF_FLAGS      = -P -traditional 
COMPILER       = -DGFORTRAN

#--------------------------------------------------------------------------
# 核心依赖库配置 (NETCDF, MPI, METIS)
#--------------------------------------------------------------------------
IOLIBS         = -lnetcdff -lnetcdf
IOINCS         = -I$(INCLUDEPATH)

MPILIB         = -lmsmpi
MPIINCS        = -I$(INCLUDEPATH)

FLAG_411       = -DMETIS_5
PARTINCS       = -I$(INCLUDEPATH)
PARTLIBS       = -lmetis

#==========================================================================
# 模型功能开关 (Flags)
#==========================================================================
FLAG_1  = -DDOUBLE_PRECISION
FLAG_3  = -DWET_DRY
FLAG_4  = -DMULTIPROCESSOR
FLAG_7  = -DDATA_ASSIM
FLAG_8  = -DLIMITED_NO
FLAG_9  = -DSEMI_IMPLICIT
FLAG_10 = -DGCN
FLAG_14 = -DRIVER_FLOAT
FLAG_15 = -DMPDATA
FLAG_28 = -DAIR_PRESSURE
FLAG_44 = -DTVD

#==========================================================================
# 编译器定义
#==========================================================================
FC             = gfortran
CC             = gcc
CXX            = g++
CFLAGS         = -O3

OPT            = -O3 -ffree-line-length-none -fallow-argument-mismatch

CPPFLAGS       = $(DEF_FLAGS) $(COMPILER)
FFLAGS         = $(OPT)
MDEPFLAGS      = --cpp --fext=f90 --file=-
RANLIB         = ranlib
AR             = ar rc

#--------------------------------------------------------------------------
# 汇总预处理标志
#--------------------------------------------------------------------------
CPPARGS = $(CPPFLAGS) $(DEF_FLAGS) $(FLAG_1) $(FLAG_2) \
          $(FLAG_3) $(FLAG_4) $(FLAG_5) $(FLAG_6)  \
          $(FLAG_7) $(FLAG_8) $(FLAG_9) $(FLAG_10) \
          $(FLAG_11) $(FLAG_12) $(FLAG_13) $(FLAG_14) \
          $(FLAG_15) $(FLAG_16) $(FLAG_17) $(FLAG_18) \
          $(FLAG_19) $(FLAG_20) $(FLAG_21) $(FLAG_22) \
          $(FLAG_23) $(FLAG_24) $(FLAG_25) $(FLAG_26) \
          $(FLAG_27) $(FLAG_28) $(FLAG_29) $(FLAG_30) \
          $(FLAG_31) $(FLAG_32) $(FLAG_33) $(FLAG_34) \
          $(FLAG_35) $(FLAG_36) $(FLAG_37) $(FLAG_38) \
          $(FLAG_39) $(FLAG_40) $(FLAG_41) $(FLAG_42) \
          $(FLAG_43) $(FLAG_44) $(FLAG_80) $(FLAG_LAM) $(FLAG_411)\
          $(FLAG_101) $(FLAG_102) $(FLAG_103) $(FLAG_104) $(FLAG_105)\
          $(FLAG_211) $(FLAG_212) $(FLAG_213) $(FLAG_251) $(FLAG_261) $(FLAG_262)\
          $(FLAG_USE_NETCDF4) $(FLAG_USE_COMPRESSION)

#--------------------------------------------------------------------------
# 最终库与包含路径汇总 (包含绝对路径补丁)
#--------------------------------------------------------------------------
INCS = $(JULIAN_INC) $(MPIINCS) $(IOINCS) $(PARTINCS) $(INCDIR) -I.
LIBS = $(LIBDIR) $(PARTLIBS) $(IOLIBS) -ljulian $(MPILIB) $(CLIB)