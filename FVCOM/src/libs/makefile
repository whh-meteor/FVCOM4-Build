#===============================================================================
# FVCOM/src/libs/Makefile - Windows MSYS2 Stable Version
#===============================================================================

SHELL = /bin/sh
RM    = rm -f

include ../make.inc

PACKAGES = julian

all:
	@echo "===================================="
	@echo "Building FVCOM Libraries (MSYS2 Mode)"
	@echo "===================================="

	# 解压 julian（如果不存在）
	if [ ! -d julian ]; then ./untar.sh julian; fi

	@echo "Building Julian..."

	cd julian && \
	$(CC) $(CFLAGS) -c *.c && \
	$(FC) $(FFLAGS) -c *.for && \
	ar rcs libjulian.a *.o

	@echo "Julian build complete."
	@echo ""
	@echo "Metis / NetCDF use MSYS2 system libraries."
	@echo "Library build finished."

clean:
	@if [ -d julian ]; then \
		cd julian && $(RM) *.o *.a ; \
	fi

clobber:
	$(RM) -r julian