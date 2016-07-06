#
#-- 64 bit version ----------------------------------------------------
iFortOptions64=-O3 -xHost -openmp -openmp-report2 -static -static-intel -parallel

MKLPATH=/opt/intel/mkl/lib/intel64
MKL95PATH=/opt/intel/mkl/include/intel64/lp64

libOptions= -L$(MKLPATH) -I$(MKL95PATH) -lmkl_blas95_lp64 -lmkl_lapack95_lp64\
    -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_intel_thread.a $(MKLPATH)/libmkl_core.a -Wl,--end-group -liomp5 -lpthread

LIBDIRS=/usr/local/Slatec
LIBSS=Slatec_AMD64_10p0_LM
#----------------------------------------------------------------------

#-- 32 bit version ----------------------------------------------------
iFortOptions32=-O3 -xHost -openmp -openmp-report2 -static -static-intel -parallel

MKLPATH_32=/opt/intel/mkl/lib/ia32
MKL95PATH_32=/opt/intel/mkl/include/ia32

libOptions_32= -L$(MKLPATH_32) -I$(MKL95PATH_32) -lmkl_blas95 -lmkl_lapack95\
    -Wl,--start-group $(MKLPATH_32)/libmkl_intel.a $(MKLPATH_32)/libmkl_intel_thread.a $(MKLPATH_32)/libmkl_core.a -Wl,--end-group -liomp5 -lpthread

LIBSS_32=Slatec_Intel32_10p1
#----------------------------------------------------------------------


Adal: Mod_Precision.f90 Mod_Init.f90 Mod_Fields.f90 Adal.f90 

	ifort $(iFortOptions64) -o Adal Mod_Precision.f90 Mod_Init.f90 Mod_Fields.f90 Adal.f90 -L$(LIBDIRS) -l$(LIBSS) $(libOptions)

