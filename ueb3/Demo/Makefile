
ifdef OPENMP
  CFLAGS += -fopenmp 
  LDFLAGS += -fopenmp
endif

ifdef OPT
  CFLAGS += -O3 -funroll-loops -DNDEBUG # -Wpedantic
endif

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3
I = -I../Include

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libhpc.a $(LDLIBS)

all: lib hpc_demo

run:
	- ./hpc_demo 4
	- ./hpc_demo 8
	- ./hpc_demo 12
	- ./hpc_demo 16
	- ./hpc_demo 20
	- ./hpc_demo 24
	- ./hpc_demo 28
	- ./hpc_demo 32
	- ./hpc_demo 36
	- ./hpc_demo 40
#	- ./hpc_demo 64

# 32 ->
# Dimension of matrix    = ( 1024, 1024)    ( = 32^2 x 32^2 )
# Number of matrix entries =      4992

lib:
	( cd ../Lib ; $(MAKE) )

hpc_demo: lib hpc_demo.c Makefile
	$(CC) $(CF) $(I) -o hpc_demo hpc_demo.c $(CS)

clean:
	- $(RM) *.o

purge: clean
	- $(RM) -r hpc_demo \
    *.a *.dSYM *.obj *.dll
