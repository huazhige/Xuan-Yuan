########################################################################################
#                       Makefile for RootFinding code
#                Written by Prof. Dongwook Lee, AMS 209, UCSC
########################################################################################

FC	= gfortran

FFLAGS_DEBUG = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all \
               -fbacktrace

FFLAGS_OPT = -ggdb -O3 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none \
             -Wuninitialized

LIBS = -I../src/ -L../src/

EXE_FILE = advect_diff.exe

OBJS  = advection.o \
	read_initFile_module.o \
	mesh_init.o \
	setup_module.o \
	advect_init.o \
	cfl.o \
	calculate_fluxes.o \
	riemann_solver.o \
	add_flux_div.o \
	bc.o

.PHONY: clean
########################################################################################
#COMPILING AND LINKING USING GENERIC SUFFIX RULE FOR F90

$(EXE_FILE) : $(OBJS)
	@$(FC) $(FFLAGS_OPT) $(OBJS) -o $(EXE_FILE) $(LIBS)
	@echo "code is now linking..."

#LET'S APPLY GENERIC SUFFIX RULE HERE FOR FORTRAN 90

# method 1 using generic suffix rule
#.SUFFIXES:
#.SUFFIXES: .F90 .o
#.F90.o:
#	$(FC) $(FFLAGS_OPT) -c $<

# method 2 using inference rule
%.o: %.F90
	$(FC) $(FFLAGS_DEBUG) -c $<

#######################################################################################
#SOME USEFUL COMMANDS
clean:
	@rm -f *.o *.mod *~ advect_diff.exe *.dat

# debug: clean
debug: FFLAGS_OPT = $(FFLAGS_DEBUG)
debug: $(EXE_FILE)

#######################################################################################
#LET'S DEFINE SOME MODULE DEPENDENCIES!
advection.o: setup_module.o mesh_init.o advect_init.o cfl.o calculate_fluxes.o \
	           add_flux_div.o bc.o

setup_module.o : read_initFile_module.o

mesh_init.o: setup_module.o

advect_init.o: setup_module.o mesh_init.o

cfl.o: setup_module.o

calculate_fluxes.o: setup_module.o riemann_solver.o

riemann_solver.o: setup_module.o mesh_init.o

add_flux_div.o: setup_module.o mesh_init.o

bc.o: setup_module.o


#######################################################################################
