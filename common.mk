OPT   += -DBIGSIM        # particle load requires 64 bit int
OPT   += -DLONGIDS       # particle ids are written in 64 bit int
OPT   += -DMAKE_LEAN

#OPT += -DUSE_OMP


#OPT   +=  -DFOF_ONLY     # no subhalos
OPT   += -DGET_GROUPVEL   # group file contains group velocities

## Default halo options for Subfind halos
OPT   += -DSUBFIND

##### AHF halos using in SUSSING mergertree comparison. NOT
## fully implemented yet -- needs two more steps of removing
## the particles in subhalos from the host halo. And the re-ordering
## of the particles based on (pseudo-) Binding Energy rank.
#OPT += -DSUSSING_TREES


#### For Rockstar halos in bgc2 format. Will enable an additional
## step for sorting the particles in a halo based on the distance
## from the center -> mimicking the Binding Energy ordering present
## in Subfind.
#OPT += -DBGC2

OPT += -DWMAP5
# OPT += -DWMAP1
# OPT += -DWMAP3



### Set the compiler -- options are icc/gcc/clang.
CC=gcc
#### Add any compiler specific flags you want
CFLAGS=
#### Add any compiler specific link flags you want
CLINK=

### You should NOT edit below this line
DISTNAME=HINGE
MINOR=0
MAJOR=1

INCLUDE := -I../io -I../utils -I.

### The POSIX_SOURCE flag is required to get the definition of strtok_r
### _BSD_SOURCE is required for strsep
CFLAGS += -Wsign-compare -Wall -Wextra -Wshadow -Wunused -std=gnu11 -g -gdwarf-3 -m64 -fPIC -D_GNU_SOURCE -D_XOPEN_SOURCE=700 -D_XOPEN_SOURCE_EXTENDED -O3 #-Ofast
GSL_CFLAGS := $(shell gsl-config --cflags)
GSL_LIBDIR := $(shell gsl-config --prefix)/lib
GSL_LINK   := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR)


ifneq (USE_OMP,$(findstring USE_OMP,$(OPT)))
  ifneq (clang,$(findstring clang,$(CC)))
     $(warning Recommended compiler for a serial build is clang)
  endif
endif

ifeq (icc,$(findstring icc,$(CC)))
  CFLAGS += -xhost -opt-prefetch #-vec-report6
  ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
		CFLAGS += -qopenmp
		CLINK  += -qopenmp
  endif
else

  ### compiler specific flags for gcc
  ifeq (gcc,$(findstring gcc,$(CC)))
		CFLAGS += -ftree-vectorize -funroll-loops -fprefetch-loop-arrays --param simultaneous-prefetches=8 #-ftree-vectorizer-verbose=6 -fopt-info-vec-missed #-fprofile-use -fprofile-correction #-fprofile-generate
    ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
			CFLAGS += -fopenmp
			CLINK  += -fopenmp
    endif
  endif

  ### compiler specific flags for clang
  ifeq (clang,$(findstring clang,$(CC)))
		CFLAGS += -funroll-loops
    ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
      $(error clang does not support OpenMP - please use gcc/icc for compiling with openmp)
     endif
  endif

  #### common options for gcc and clang
#  CFLAGS  += -march=native
	CFLAGS  += -Wformat=2  -Wpacked  -Wnested-externs -Wpointer-arith  -Wredundant-decls  -Wfloat-equal -Wcast-qual
  CFLAGS  +=  -Wcast-align -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs -Wstrict-prototypes  #-D_POSIX_C_SOURCE=2 -Wpadded -Wconversion
  CLINK += -lm
endif
