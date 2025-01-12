# --- Variables ---
PROG = ../build/ptb_dev
OSTYPE = LINUXI
WITH_DFTD4 = FALSE

# --- Code-specific variables ---
SOURCES = accuracy.f90 la.f90 blas_level1.f90 blas_level2.f90 blas_level3.f90 lapack_eig.f90
$(foreach src, $(SOURCES), $(eval $(src) := $(patsubst %.f90, %.o, $(src))))

MODULES= bascom.o parcom.o com.o thresholds.o cbascom.o mocom.o dftd4.o aescom.o \
			$(foreach src, $(SOURCES), $($(src)))

CFILES= symmetry_i.o

OBJS1 = purification_settings.f90 timer.o matrix_ops.o metrics.o rdcoord.o readl.o asym.o elem.o readl2.o readline.o prmat.o  timing.o echo.o solve.o rdstvint.o \
       ints.o intpack.o pop.o rdmat.o fermi.o reorder.o bas.o occ.o dtrf2.o lin.o energy.o mrec.o \
       cma.o dtrf3.o drsp.o mlpop.o calcrab.o symmetry.o pgtb.o eeq.o ncoord.o prdipole.o itabrow.o \
       sscal.o cbas.o gauss.o printmos.o qsort.o rdtm.o shellocc_ref.o wr_control_atoms.o \
       egtbcorr.o  guess_qsh.o eptb.o symmetry2.o symtranslib.o spinden.o camm.o json_output.o
     
       
MAINFILES= main.o 

OBJS = $(MODULES) $(CFILES) $(OBJS1) $(MAINFILES)

# --- Compiler-specific  variables ---
PREFLAG = -E -P
CCFLAGS = -O -DLINUX

ifeq ($(OSTYPE),LINUXI)
  FC = ifx
  CC  = icc
  LINKER = $(FC) -g -O -qopenmp  -I$(MKLROOT)/include -qmkl-ilp64=parallel
  LIBS = -lpthread -lm -lstdc++ -qmkl-ilp64=parallel
  FFLAGS = -O3 -align array64byte -i8 -march=core-avx2 -mcmodel=large \
  				-heap-arrays -check bounds -warn interfaces
  CCFLAGS = -O -DLINUX
else ifeq ($(OSTYPE),LINUXL)
  FC = lfc 
  CC = icc 
  LINKER = lfc  
  LIBS = /home/grimme/source/lapack-3.4.2/liblapack.a /home/grimme/source/lapack-3.4.2/librefblas.a
  FFLAGS = -O --ntrace --tpp --info --prefetch 2
else
  $(error Unsupported OSTYPE: $(OSTYPE))
endif

ifeq ($(WITH_DFTD4),TRUE)
  DFTD4_PATH = $(PWD)/dftd4/install
  LINKER += -I$(DFTD4_PATH)/include -L$(DFTD4_PATH)/lib64
  LIBS := -ldftd4 $(LIBS)
  FFLAGS += -I$(DFTD4_PATH)/include
  PREFLAG += -DWITH_DFTD4
endif

# --- Default Target ---
.PHONY: all clean
all: create_build_dir $(PROG)


create_build_dir:
	@mkdir -p ../build


#--------------------------------------------------------------------------
# Setting the rules to make the object files:
#--------------------------------------------------------------------------
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
	rm -rf *.o *.mod *genmod* $(PROG) ../build
