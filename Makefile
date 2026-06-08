# Set compiler, flags and netcdf library according to machine we are using:

FC = mpif90
NETCDF = -I/usr/lib64/gfortran/modules
NETCDFLIB = -L/usr/lib64 -lnetcdff -lnetcdf
FFLAGS = -O3 -pedantic -Wuninitialized -fimplicit-none -Wall -Wextra -funroll-loops --param max-unroll-times=5 #-fcheck=all


TARGET = dumfric

# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

SRCDIR = src
OBJDIR = obj
BINDIR = bin

OBJFILES = main.o time.o grid.o evolve.o mpitools.o params.o init.o vecpot.o ncio.o friction.o outflow.o ohmic.o hyper.o efield.o bfield.o current.o driver.o diagnostics.o ffts.o tests.o
FULLTARGET = $(BINDIR)/$(TARGET)

VPATH = $(SRCDIR):$(OBJDIR)

# Rule to build the fortran files

%.o: %.f90
	@mkdir -p $(BINDIR) $(OBJDIR)
	$(FC) -c $(FFLAGS) $(NETCDF) -J $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.F90
	@mkdir -p $(BINDIR) $(OBJDIR) 
	$(FC) -c $(FFLAGS) $(NETCDF) -J $(OBJDIR) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) -J $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(NETCDFLIB)

.PHONEY: clean
clean:
	@rm -rf *~ $(BINDIR) $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

.PHONEY: tidy
tidy:
	@rm -rf $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log *.nc *.html *.unf
	
# All the dependencies
params.o: params.f90
ffts.o: ffts.f90
mpitools.o: mpitools.f90 params.o
time.o: time.f90 mpitools.o params.o
grid.o: grid.f90 mpitools.o params.o
ncio.o: ncio.f90 mpitools.o params.o grid.o
current.o: current.f90 mpitools.o params.o grid.o ncio.o
outflow.o: outflow.f90 mpitools.o params.o grid.o
ohmic.o: ohmic.f90 mpitools.o params.o grid.o friction.o
hyper.o: hyper.f90 mpitools.o params.o grid.o friction.o
driver.o: driver.f90 mpitools.o grid.o time.o
bfield.o: bfield.f90 mpitools.o params.o grid.o hyper.o ncio.o
efield.o: efield.f90 mpitools.o params.o grid.o friction.o driver.o ohmic.o hyper.o
friction.o: friction.f90 mpitools.o params.o grid.o ncio.o
vecpot.o: vecpot.f90 mpitools.o params.o grid.o ncio.o ffts.o
diagnostics.o: diagnostics.f90 mpitools.o params.o grid.o time.o driver.o bfield.o current.o efield.o vecpot.o hyper.o friction.o outflow.o
evolve.o: evolve.f90 grid.o params.o friction.o outflow.o vecpot.o efield.o bfield.o current.o driver.o ohmic.o hyper.o ncio.o mpitools.o diagnostics.o
init.o: init.f90 params.o mpitools.o grid.o time.o vecpot.o friction.o outflow.o ohmic.o hyper.o driver.o diagnostics.o bfield.o
tests.o: tests.f90 params.o mpitools.o grid.o vecpot.o bfield.o diagnostics.o current.o
main.o: main.f90 evolve.o mpitools.o params.o init.o tests.o
