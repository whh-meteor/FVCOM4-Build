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

	# 先规范化行尾并用更鲁棒的正则删除只包含 $(LIBRARY) 或 libjulian.a 的行，然后打印确认
	if [ -f julian/Makefile ]; then \
	  echo "Patching julian/Makefile (normalize CRLF -> LF)"; \
	  dos2unix julian/Makefile >/dev/null 2>&1 || true; \
	  # 删除只包含 $(LIBRARY) 的那一行（考虑空白和 CR） \
	  perl -0777 -i -pe 's/^[ \t]*\$\(\s*LIBRARY\s*\)[ \t\r]*\n//gim' julian/Makefile; \
	  # 保险起见也删除只包含 libjulian.a 的独立行 \
	  perl -0777 -i -pe 's/^[ \t]*libjulian\.a[ \t\r]*\n//gim' julian/Makefile; \
	  echo "=== julian/Makefile (head after patch) ==="; nl -ba julian/Makefile | sed -n '1,80p'; \
	else \
	  echo "julian/Makefile not found!"; \
	fi

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