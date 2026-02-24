#===========================================================================
# FVCOM v4.4.10 make.inc - Windows MSYS2/MinGW 适配版
# 用于桌面软件集成与 GitHub Actions 编译
#===========================================================================

#========== TOPDIR ========================================================
# TOPDIR 是源码所在的根目录。在 GitHub Actions 中，通常会自动识别。
           TOPDIR        = $(shell pwd)/..
           INSTALLDIR    = $(TOPDIR)/libs/install

#===========================================================================
# 预处理器选项 (CPP)
#===========================================================================
           DEF_FLAGS     = -P -traditional 

#===========================================================================
# 环境路径设置 - 适配 MSYS2 Mingw64 标准路径
#===========================================================================
           # 这些变量会被后续的 IOLIBS, PARTLIBS 等引用
           LIBPATH       = /mingw64/lib
           INCLUDEPATH   = /mingw64/include

           LIBDIR        = -L$(LIBPATH)
           INCDIR        = -I$(INCLUDEPATH)

#--------------------------------------------------------------------------
# 标准库 (FVCOM 内部 julian 库)
#--------------------------------------------------------------------------
           DTLIBS        = -ljulian
           DTINCS        = 

#--------------------------------------------------------------------------
# NETCDF 配置 (必须使用 NetCDF-Fortran)
#--------------------------------------------------------------------------
           IOLIBS        = -lnetcdff -lnetcdf
           IOINCS        = -I/mingw64/include

#==========================================================================
# 模型功能开关 (核心配置)
#==========================================================================

# 1. 精度：建议桌面集成使用双精度
           FLAG_1 = -DDOUBLE_PRECISION

# 2. 干湿网格
           FLAG_3 = -DWET_DRY

# 3. 并行计算 (MPI)
           FLAG_4 = -DMULTIPROCESSOR

# 4. METIS 划分库配置 (假设您已按照方案二编译了 METIS_source)
           FLAG_411 = -DMETIS_5
           PARLIB   = 
           PARTINCS = -I$(TOPDIR)/METIS_source/include
           PARTLIBS = -L$(TOPDIR)/METIS_source -lmetis

# 5. 其他物理/数学标志 (保留您的原始选择)
           FLAG_7  = -DDATA_ASSIM
           FLAG_8  = -DLIMITED_NO
           FLAG_9  = -DSEMI_IMPLICIT
           FLAG_10 = -DGCN
           FLAG_14 = -DRIVER_FLOAT
           FLAG_15 = -DMPDATA
           FLAG_28 = -DAIR_PRESSURE
           FLAG_44 = -DTVD

#==========================================================================
# 编译器与平台定义 (针对 Windows/gfortran 深度定制)
#==========================================================================

         # 预处理器路径
         CPP      = /usr/bin/cpp 
         COMPILER = -DGFORTRAN
         
         # 编译器定义
         # 注意：在 MSYS2 下链接 MS-MPI，FC 直接用 gfortran 即可
         FC       = gfortran
         CC       = gcc
         CXX      = g++
         CFLAGS   = -O3

         # 核心编译参数 (关键说明):
         # -O3: 开启最高级优化
         # -ffree-line-length-none: 防止旧代码行过长报错
         # -fallow-argument-mismatch: 核心补丁！解决 FVCOM 4 在现代 gfortran 下的兼容性问题
         OPT      = -O3 -ffree-line-length-none -fallow-argument-mismatch
         
         # 显式链接 Microsoft MPI
         MPILIB   = -lmsmpi
         CLIB     = 

#==========================================================================
# 系统命令汇总
#==========================================================================
         CPPFLAGS  = $(DEF_FLAGS) $(COMPILER)
         FFLAGS    = $(DEBFLGS) $(OPT)
         MDEPFLAGS = --cpp --fext=f90 --file=-
         RANLIB    = ranlib
         AR        = ar rc

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
# 最终库文件链接顺序 (顺序极其重要)
#--------------------------------------------------------------------------
         LIBS = $(LIBDIR) $(CLIB) $(PARLIB) $(PARTLIBS) $(IOLIBS) $(DTLIBS) \
                $(MPILIB) $(GOTMLIB) $(KFLIB) $(BIOLIB) \
                $(OILIB) $(VISITLIB) $(PROJLIBS) $(PETSC_LIB)

         INCS = $(INCDIR) $(IOINCS) $(GOTMINCS) $(BIOINCS) \
                $(VISITINCPATH) $(PROJINCS) $(DTINCS) $(PARTINCS) \
                $(PETSC_FC_INCLUDES)