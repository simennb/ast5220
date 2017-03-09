# -*- Makefile -*-

FC        = ifort     # F90 compiler
OPTIM     = -O3       # Optimization flags; set to -g -C for useful debug information

# It should hopefully not be necessary to edit anything below
FITSDIR = /mn/stornext/u3/hke/local/lib      # Directory containing libcfitsio.a
LAPACK  = -L/mn/stornext/u3/hke/local/lib -llapack -lblas
HEALPIX = -L/mn/stornext/u3/hke/local/lib -lhealpix
HEALINC = -I/mn/stornext/u3/hke/local/include
OUTPUT  = cmbspec
FFLAGS  = $(HEALPIX) $(LAPACK)

# List of source files to be compiled
OBJS    = math_tools.o spline_1D_mod.o rk_mod.o bs_mod.o ode_solver.o \
	  params.o time_mod.o cmbspec.o

# Linking stage
cmbspec: $(OBJS)
	$(FC) $(FFLAGS) $(OPTIM) -o $(OUTPUT) $(OBJS) $(LAPACK)

# Dependencies
cmbspec.o       : time_mod.o
time_mod.o      : params.o spline_1D_mod.o

# Compilation of source files
%.o : %.f90
	$(FC) $(OPTIM) $(HEALINC) -c $<

# Clean-up command (write "make clean")
.PHONY: clean
clean:
	rm *.mod *.o *~ cmbspec
