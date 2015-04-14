
topdir = PARDISO
libdir = $(topdir)/lib

include Makefile.inc
OS := $(shell uname)

ifeq ($(OS),Darwin)
	LIBBLAS = /System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/libBLAS.dylib
else
	 LIBBLAS = -fopenmp /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a 
endif

#LIB = -L/users/verbof/pardiso -lmetis41_pardiso -lmetis41-P_pardiso -lpardiso -lpils_pardiso  $(LIBBLAS)
#LIB = -L/users/verbof/Libraries/linuxAMD64 -lpardiso500-INTEL1301-X86-64

SRC = call_pardiso.f90 call_pardiso_z.f90
LIB = ${HOME}/Libraries/linuxAMD64/libpardiso500-INTEL1301-X86-64.so

all: $(SRC)
	 $(FL) -shared -fPIC -fcray-pointer $(SRC) -o ../lib/PARDISO  $(LIB)
clean:
	rm ../lib/PARDISO
