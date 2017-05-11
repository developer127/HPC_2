# On the Mac, you need to first install Xcode.  Then
# type this command in the Terminal:
#       xcode-select --install
# to ensure that you have the command-line tools installed.
#
# To just compile the library, do:
#       make
# resp.  
#       make VERBOSE=1
#       make CC=gcc OPENMP=1
#       make CC=/opt/local/bin/clang-mp-3.7 OPT=1 OPENMP=1
#
# To run a demo using the library 
#       cd ../Demo ; make
#

LIBRARY = libhpc

ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG -Wpedantic  # -mfpmath=sse #-march=native #
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -fPIC
I = -I../Include
RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)
CP = cp -f

AR_TARGET = $(LIBRARY).a

all: install

HPC = \
 sed_buildS.o \
 sed_dupl.o \
 sed_gaxpy.o \
 sed_print.o \
 sed_gs_constr.o \
 sed_util.o  \
 hpc_cumsum.o \
 hpc_malloc.o \
 mesh_buildRhs.o \
 mesh_load.o  \
 mesh_print.o  \
 mesh_refine.o  \
 mesh_getFixed.o  \
 mesh_util.o  

# sed_compress.o \


# bcrs_compress.o \
# bcrs_gaxpy.o \
# bcrs_print.o \
# bcrs_util.o \
# cs_apat.o \
# cs_dupl.o \
# cs_lapmat_p1_square.o \
# cs_lapmat_p2_square.o \
# cs_lame_2d_p1_square.o \
# gem_compress.o \
# gem_gauss.o \
# gem_gausssol.o \
# gem_entry.o \
# gem_load.o \
# gem_print.o \
# gem_util.o \
# hpc_compress.o \
# hpc_dupl.o \
# hpc_gaxpy.o \
# hpc_util.o 
# hpc_print.o \
# jds_csr2jds.o \
# jds_gaxpy.o \
# jds_print.o \
# jds_util.o  \
# sky_cholesky.o  \
# sky_cholsol.o  \
# sky_load.o  \
# sky_print.o  \
# sky_util.o 

$(HPC): ../Include/hpc.h Makefile

%.o: ../Source/%.c ../Include/hpc.h
	$(CC) $(CF) $(I) -c $<

static: $(AR_TARGET)

$(AR_TARGET): $(HPC)
	$(ARCHIVE)  $@ $^
	- $(RANLIB) $@

# install archive in this directory
install: $(AR_TARGET)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) *.a *.obj *.dll *.dylib *.so *.so.*

