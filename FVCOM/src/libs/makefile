#----------- FVCOM/src/libs/Makefile - Fixed for MSYS2 ----------------
SHELL         = /bin/sh
include ../make.inc

# 确保所有需要的包都在列表中
PACKAGES = julian metis proj fproj netcdf

MYINSTALLDIR = $(INSTALLDIR)

all:
	@echo "Starting to build all libraries..."
	# 1. 循环解压
	for item in $(PACKAGES); do \
		if [ ! -d $$item ]; then ./untar.sh $$item; fi; \
	done
	
	# 2. 编译 Julian (FVCOM 核心时间库)
	@echo "Building Julian..."

	# 强制删除 Windows 下会被执行的库文件行
	sed -i '/LIBRARY/d' julian/[Mm]akefile

	cd julian && $(MAKE) CC="$(CC)" FC="$(FC)" AR="ar rs"

	# 3. 编译 Metis
	@echo "Building Metis..."
	cd metis && $(MAKE) CC="$(CC)" AR="ar rs"

	# 4. 编译 Proj & FProj
	@echo "Building Proj..."
	cd proj && ./configure CC="$(CC)" CFLAGS="-O3" --prefix=$(MYINSTALLDIR) --disable-shared
	cd proj && $(MAKE) install
	cd fproj && ./configure CC="$(CC)" FC="$(FC)" --with-proj4=$(MYINSTALLDIR) --prefix=$(MYINSTALLDIR)
	cd fproj && $(MAKE) install

	# 5. 编译 NetCDF (如果系统没安装)
	@echo "Building NetCDF..."
	cd netcdf && ./configure CC="$(CC)" FC="$(FC)" --prefix=$(MYINSTALLDIR)
	cd netcdf && $(MAKE) install

clean:
	for item in $(PACKAGES); do (cd $$item && $(MAKE) clean) || true; done

clobber:
	for item in $(PACKAGES); do (rm -rf $$item) || true; done
	rm -rf $(MYINSTALLDIR)/*