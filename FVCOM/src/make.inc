#===========================================================================
# FVCOM v4.4.10 make.inc - 终极适配 (解压库补丁版)
#===========================================================================

#---------- 1. 物理路径自动转换 ----------
MSYS_ROOT     := $(shell cygpath -m /mingw64)
CURR_DIR      := $(shell pwd)
# 强制定位 src/libs 的物理地址
LIBS_W        := $(shell cygpath -m $(CURR_DIR)/libs)

TOPDIR        = ..
INSTALLDIR    = $(TOPDIR)/libs/install

#---------- 2. 内部库包含路径 (指向解压后的目录) ----------
JULIAN_INC    = -I$(LIBS_W)/julian
METIS_INC     = -I$(LIBS_W)/metis/include
PROJ_INC      = -I$(LIBS_W)/proj/include

#---------- 3. 编译器与预处理器 ----------
CPP           = gcc -E
DEF_FLAGS     = -P -traditional 
COMPILER      = -DGFORTRAN

FC            = gfortran
CC            = gcc
CXX           = g++
OPT           = -O3 -ffree-line-length-none -fallow-argument-mismatch

#---------- 4. 依赖库配置 (NetCDF, MPI, METIS) ----------
# 优先使用 MSYS2 系统提供的 NetCDF 和 MPI
IOLIBS        = -lnetcdff -lnetcdf
MPILIB        = -lmsmpi

#---------- 5. 汇总 INCS 与 LIBS (极其关键) ----------
# 加入 $(JULIAN_INC) 以修复 fjulian.inc 缺失问题
INCS = $(JULIAN_INC) $(METIS_INC) $(PROJ_INC) -I$(MSYS_ROOT)/include -I.
# 链接顺序：必须先链接本地生成的 libjulian.a
LIBS = -L$(LIBS_W)/julian -ljulian -L$(MSYS_ROOT)/lib -lmetis $(IOLIBS) $(MPILIB)

#---------- 6. 预处理标志 ----------
CPPARGS = $(DEF_FLAGS) $(COMPILER) -DDOUBLE_PRECISION -DWET_DRY -DMULTIPROCESSOR \
          -DDATA_ASSIM -DLIMITED_NO -DSEMI_IMPLICIT -DGCN -DRIVER_FLOAT -DMPDATA \
          -DAIR_PRESSURE -DTVD -DMETIS_5