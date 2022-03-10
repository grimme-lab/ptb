PROG = ~/bin/gtb

#--------------------------------------------------------------------------
 OSTYPE=LINUXI
#--------------------------------------------------------------------------

SOURCES = accuracy.f90 la.f90 blas_level1.f90 blas_level2.f90 blas_level3.f90 lapack_eig.f90
$(foreach src, $(SOURCES), $(eval $(src) := $(patsubst %.f90, %.o, $(src))))

MODULES= iso_fortran_env.o bascom.o parcom.o com.o thresholds.o cbascom.o mocom.o \
			$(foreach src, $(SOURCES), $($(src)))

CFILES= symmetry_i.o

OBJS1= 

OBJS2= rdcoord.o readl.o asym.o elem.o readl2.o readline.o prmat.o timing.o echo.o solve.o rdstvint.o \
       ints.o intpack.o pop.o rdmat.o fermi.o reorder.o bas.o occ.o dtrf2.o lin.o energy.o mrec.o \
       cma.o dtrf3.o drsp.o mlpop.o calcrab.o pgtb.o eeq.o ncoord.o prdipole.o itabrow.o \
       sscal.o cbas.o gauss.o printmos.o qsort.o momatch2.o rdtm.o shellocc_ref.o wr_control_atoms.o \
       egtbcorr.o symmetry.o
     
OBJS3= 
      
OBJS4 = 
       
MAINFILES= main.o 

OBJS = $(MODULES) $(CFILES) $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4) $(MAINFILES)
#--------------------------------------------------------------------------

WITH_DFTD4=FALSE

ifeq ($(OSTYPE),LINUXI)
  PREOPTS =
  FC = ifort
  #FC = gfortran
  CC  = gcc
  LINKER = ifort -g -O -static -fopenmp  -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
  LIBS = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
  PREFLAG = -E -P

   FFLAGS = -xHost -O3 -align array64byte
#  FFLAGS = -check all -g -debug all -traceback -fpe0 -xHost
#  FFLAGS = -O1 
#  FFLAGS = -O -axAVX -qopenmp -align array64byte 
  CCFLAGS = -O -DLINUX
endif

ifeq ($(OSTYPE),LINUXL)
  FC = lfc 
  CC = icc 
  LINKER = lfc  
  PREFLAG = -E -P
  CCFLAGS = -O -DLINUX
  LIBS = /home/grimme/source/lapack-3.4.2/liblapack.a /home/grimme/source/lapack-3.4.2/librefblas.a
  FFLAGS = -O --ntrace --tpp --info --prefetch 2
# FFLAGS = -O --trace --tpp --info --chk a,e,s,u --prefetch 2
endif

ifeq ($(WITH_DFTD4),TRUE)
  DFTD4_PATH = $(PWD)/dftd4/install
  LINKER += -I$(DFTD4_PATH)/include -L$(DFTD4_PATH)/lib64
  LIBS := -ldftd4 $(LIBS)
  FFLAGS += -I$(DFTD4_PATH)/include
  PREFLAG += -DWITH_DFTD4
endif


ifeq ($(OSTYPE),SOMETHINGELSE)
  PREOPTS =
  FC = 
  CC = 
  LINKER = 
  LIBS =
  PREFLAG = 
  FFLAGS = 
  CCFLAGS = 
endif


.PHONY: all
.PHONY: clean

all: $(PROG)


#--------------------------------------------------------------------------
# Setting the rules to make the object files:
#--------------------------------------------------------------------------
# implizite Regel zur Erzeugung von *.o aus *.F ausschalten
%.o: %.F
fortran.o : $(CUDA)/src/fortran.c
	gcc $(CCFLAGS) -I$(CUDA)/src/ $(INC) -c  $(CUDA)/src/fortran.c

# aus *.F mache ein *.f
%.f: %.F
	@echo "making $@ from $<"
	$(CC) $(PREFLAG) $(PREOPTS) $< -o $@

# aus *.f mache ein *.o
%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

# aus *.c mache ein *.o
%.o: %.c
	@echo "making $@ from $<"
	$(CC) $(CCFLAGS) -c $< -o $@

# linken
$(PROG): $(OBJS) 
	$(LINKER) $(OBJS) $(LIBS) -o $(PROG)


# dependencies
$(blas_level1.f90): $(accuracy.f90)
$(blas_level2.f90): $(accuracy.f90)
$(blas_level3.f90): $(accuracy.f90)
$(lapack_eig.f90): $(accuracy.f90)
$(la.f90): $(blas_level1.f90) $(blas_level2.f90) $(blas_level3.f90) $(lapack_eig.f90)

#aufraeumen
clean:
	rm -f *.o $(PROG) 
	rm -f *.mod $(PROG) 
	rm -f $(patsubst %.F, %.f, $(wildcard *.F))


