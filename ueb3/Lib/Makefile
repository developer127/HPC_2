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

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -fPIC -g
I = -I../Include
RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)
CP = cp -f

AR_TARGET = $(LIBRARY).a

all: install

HPC = \
 cs_cg.o \
 cs_compress.o \
 cs_entry.o \
 cs_gaxpy.o \
 cs_lapmat_p1_square.o \
 cs_load.o \
 cs_print.o \
 cs_util.o \
 gem_compress.o \
 gem_print.o \
 gem_gaxpy.o \
 gem_gauss.o \
 gem_gausssol.o \
 gem_util.o \
 hpc_cumsum.o \
 hpc_dot.o \
 hpc_malloc.o \
 hpc_veccopy.o \
 hpc_vecinitzero.o \
 hpc_vecprint.o \
 hpc_vecsum.o \
 sky_cholesky.o  \
 sky_cholsol.o  \
 sky_compress.o  \
 sky_gaxpy.o  \
 sky_print.o  \
 sky_util.o 

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

