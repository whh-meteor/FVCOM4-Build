#===========================================================================
# FVCOM v4.4.10 make.inc - Windows MSYS2/MinGW 适配版
# 专门优化用于 GitHub Actions (msys2 运行环境)
#===========================================================================

#========== TOPDIR ========================================================
# 使用相对路径最稳健，src 的上一层就是 FVCOM 根目录
           TOPDIR        = ..
           INSTALLDIR    = $(TOPDIR)/libs/install

#===========================================================================
# 预处理器选项 (CPP)
#===========================================================================
           DEF_FLAGS     = -P -traditional 

#===========================================================================
# 环境路径设置 - 适配 MSYS2 Mingw64 标准路径
#===========================================================================
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
# NETCDF 配置 (针对 MSYS2 环境)
#--------------------------------------------------------------------------
           IOLIBS        = -lnetcdff -lnetcdf
           IOINCS        = -I/mingw64/include

#--------------------------------------------------------------------------
# MPI 配置 (Windows 专属补丁)
# 原文件中可能没有 MPIINCS，这里我们手动定义以确保编译器能找到 mpif.h
#--------------------------------------------------------------------------
           MPIINCS       = -I/mingw64/include
           MPILIB        = -lmsmpi

#==========================================================================
# 模型功能开关 (Flags)
#==========================================================================
           FLAG_1 = -DDOUBLE_PRECISION
           FLAG_3 = -DWET_DRY
           FLAG_4 = -DMULTIPROCESSOR

           # METIS 路径：确保指向 TOPDIR 下的 METIS_source
           FLAG_411 = -DMETIS_5
           PARLIB   = 
           PARTINCS = -I/mingw64/include
           PARTLIBS = -L/mingw64/lib -lmetis

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
         CPP      = cpp 
         COMPILER = -DGFORTRAN
         
         FC       = gfortran
         CC       = gcc
         CXX      = g++
         CFLAGS   = -O3

         # OPT 包含核心补丁 -fallow-argument-mismatch
         OPT      = -O3 -ffree-line-length-none -fallow-argument-mismatch
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
# 最终库文件链接顺序 (调整了 MPI 和 NetCDF 的层级)
#--------------------------------------------------------------------------
         LIBS = $(LIBDIR) $(CLIB) $(PARLIB) $(PARTLIBS) $(IOLIBS) $(DTLIBS) \
                $(MPILIB) $(GOTMLIB) $(KFLIB) $(BIOLIB) \
                $(OILIB) $(VISITLIB) $(PROJLIBS) $(PETSC_LIB)

         # 这里手动加入了 $(MPIINCS) 确保编译器能找到 mpif.h
         INCS = $(INCDIR) $(IOINCS) $(MPIINCS) $(GOTMINCS) $(BIOINCS) \
                $(VISITINCPATH) $(PROJINCS) $(DTINCS) $(PARTINCS) \
                $(PETSC_FC_INCLUDES)