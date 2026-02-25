#===========================================================================
# FVCOM v4.4.10 make.inc - MSYS2 CI 稳定版
#===========================================================================

#---------- 1. 当前目录 ----------
CURR_DIR      := $(shell pwd)
LIBS_DIR      := $(CURR_DIR)/libs

TOPDIR        = ..
INSTALLDIR    = $(TOPDIR)/libs/install

#---------- 2. 内部库路径 ----------
JULIAN_INC    = -I$(LIBS_DIR)/julian
METIS_INC     = -I$(LIBS_DIR)/metis/include
PROJ_INC      = -I$(LIBS_DIR)/proj/include

#---------- 3. 编译器 ----------
CPP           = gcc -E
DEF_FLAGS     = -P -traditional
COMPILER      = -DGFORTRAN

FC            = gfortran
CC            = gcc
CXX           = g++
OPT           = -O3 -ffree-line-length-none -fallow-argument-mismatch

#---------- 4. 使用 pkg-config 获取 NetCDF ----------
NETCDF_CFLAGS := $(shell pkg-config --cflags netcdf-fortran)
NETCDF_LIBS   := $(shell pkg-config --libs netcdf-fortran)

#---------- 5. MPI (如不需要可删除) ----------
# MPILIB        = -lmsmpi

#---------- 6. 汇总 INCS ----------
INCS = $(JULIAN_INC) $(METIS_INC) $(PROJ_INC) $(NETCDF_CFLAGS) -I.

#---------- 7. 汇总 LIBS ----------
# 顺序很重要：先 julian，再 metis，再 netcdf
LIBS = -L$(LIBS_DIR)/julian -ljulian -lmetis $(NETCDF_LIBS) $(MPILIB)

#---------- 8. 预处理标志 ----------
CPPARGS = $(DEF_FLAGS) $(COMPILER) \
          -DDOUBLE_PRECISION -DWET_DRY \
          -DMETIS_5