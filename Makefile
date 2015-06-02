#=======================================================================
# Makefile for Imcompact3D
#=======================================================================

# Choose pre-processing options
#   -DSHM	   - enable shared-memory implementation
#   -DDOUBLE_PREC  - use double-precision
OPTIONS = -DDOUBLE_PREC

# Choose an FFT engine, available options are:
#   fftw3      - FFTW version 3.x
#   generic    - A general FFT algorithm (no 3rd-party library needed)
FFT= generic

# Paths to FFTW 3
FFTW3_PATH=lib/gcc/fftw3/double-mpi/3.2.2   # full path of FFTW installation if using fftw3 engine above
FFTW3_INCLUDE = -I/opt/fftw/3.3.4.0/sandybridge/include
FFTW3_LIB = -L/opt/fftw/3.3.4.0/sandybridge/lib -lfftw3_mpi -lfftw3f_mpi

# Specify Fortran and C compiler names and flags here
# Normally, use MPI wrappers rather than compilers themselves 
# Supply a Fortran pre-processing flag together with optimisation level flags
# Some examples are given below:

# Intel compiler+Intelmpi
FC = mpiifort
OPTFC = -O3 -real-size 32 -funroll-loops -cpp -zero # -opt-report=2
CC = mpiicc
CFLAGS = -O3 -zero

# Intel compiler+openmpi
#FC = mpif90
#OPTFC = -O3 -real-size 32 -funroll-loops -cpp -zero # -opt-report=2
#CC = mpicc
#CFLAGS = -O3 -zero


# GNU
#FC = mpif90
#OPTFC = -O3 -fdefault-real-8 -funroll-loops -ftree-vectorize -fcray-pointer -cpp
#CC = mpicc
#CFLAGS = -O3

# Cray
#FC = ftn
#OPTFC = -e Fm
#CC = cc
#CFLAGS = 

# PGI
#FC = ftn
#OPTFC = -fast -O3 -Mpreprocess
#CC = cc
#CFLAGS = -O3

# PathScale
#FC = ftn
#OPTFC = -Ofast -cpp
#CC = cc
#CFLAGS = -O3

#-----------------------------------------------------------------------
# Normally no need to change anything below

# include PATH 
ifeq ($(FFT),generic)
  INC=
else ifeq ($(FFT),fftw3)
  INC=$(FFTW3_INCLUDE)
endif

# library path
ifeq ($(FFT),generic)
   LIBFFT=
else ifeq ($(FFT),fftw3)
   LIBFFT=$(FFTW3_LIB)
endif

# List of source files
SRC = decomp_2d.f90 glassman.f90 fft_$(FFT).f90 module_param.f90 io.f90 variables.f90 poisson.f90 schemes.f90 implicit.f90 convdiff.f90 myconjht.f90 user_module.f90 scalar_exp.f90 incompact3d.f90 navier.f90 derive.f90 parameters.f90 tools.f90 visu.f90

#-----------------------------------------------------------------------
# Normally no need to change anything below

ifneq (,$(findstring DSHM,$(OPTIONS)))
SRC := FreeIPC.f90 $(SRC)  
OBJ =	$(SRC:.f90=.o) alloc_shm.o FreeIPC_c.o
else
OBJ =	$(SRC:.f90=.o)
endif	

all: incompact3d

alloc_shm.o: alloc_shm.c
	$(CC) $(CFLAGS) -c $<

FreeIPC_c.o: FreeIPC_c.c
	$(CC) $(CFLAGS) -c $<

incompact3d : $(OBJ)
	$(FC) $(OPTFC) -o $@ $(OBJ) $(LIBFFT) $(N8HPC_LINALG_FFLAGS)
#-mkl
#-L/usr/lib64 -llapack -lblas
%.o : %.f90
	$(FC) $(OPTFC) $(OPTIONS) $(INC) -c $<

.PHONY: clean 
clean:
	rm -f *.o *.mod incompact3d

.PHONY: realclean
realclean: clean
	rm -f *~ \#*\#
