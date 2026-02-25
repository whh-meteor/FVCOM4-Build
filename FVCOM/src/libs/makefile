#----------- FVCOM/src/libs/Makefile (MSYS2 Simplified Version) --------

SHELL = /bin/sh
include ../make.inc

# Windows/MSYS2 版本只编译 julian
PACKAGES = julian

MYINSTALLDIR = $(INSTALLDIR)

all:
	@echo "===================================="
	@echo "Building FVCOM Libraries (MSYS2 Mode)"
	@echo "===================================="

	# 1. 解压 julian（如果不存在）
	if [ ! -d julian ]; then ./untar.sh julian; fi

	# 2. 编译 Julian
	@echo "Building Julian..."
	cd julian && \
	$(CC) -O3 -c *.c && \
	$(FC) -O3 -c *.for && \
	ar rs libjulian.a *.o

	@echo "Julian build complete."
	@echo ""
	@echo "Metis / Proj / NetCDF will use MSYS2 system libraries."
	@echo "Library build finished."

clean:
	@if [ -d julian ]; then \
		cd julian && rm -f *.o *.a ; \
	fi

clobber:
	rm -rf julian