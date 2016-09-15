FC := gfortran

OS := $(shell uname)

ifeq ($(OS),Darwin)
	LIBBLAS = /System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/libBLAS.dylib
else
	 LIBBLAS = -fopenmp /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a 
endif

SRC = call_pardiso.f90 call_pardiso_z.f90

LIBPARDISO = ${HOME}/Libraries/osx/libpardiso500-MACOS-X86-64.dylib

all: $(SRC)
	 $(FC) -shared -fPIC $(SRC) -o ../lib/PARDISO $(LIBPARDISO)
clean:
	rm ../lib/PARDISO
