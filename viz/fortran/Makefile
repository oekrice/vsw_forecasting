# Set FORTRAN90 compiler
FC = mpif90
NETCDF = -I/usr/lib64/gfortran/modules
NETCDFLIB = -L/usr/lib64 -lnetcdff -lnetcdf
FFLAGS = -O3 -pedantic -Wuninitialized -fimplicit-none -Wall -Wextra -funroll-loops --param max-unroll-times=5 -fopenmp #-fcheck=all

# Set build parameters
TARGET = tracer

# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

SRCDIR = src
OBJDIR = obj
BINDIR = bin

OBJFILES = shared.o fluxes.o ffts.o readnc.o trace.o xtrapol.o outeqm.o main.o
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
	@rm -rf $(OBJDIR) *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

	
# All the dependencies
shared.o: shared.f90
ffts.o: ffts.f90
readnc.o: readnc.f90 shared.o ffts.o
fluxes.o: fluxes.f90 readnc.o shared.o
trace.o: trace.f90 shared.o readnc.o
xtrapol.o: xtrapol.f90 shared.o readnc.o ffts.o
outeqm.o: outeqm.f90 shared.o readnc.o ffts.o
main.o: main.f90 trace.o readnc.o shared.o xtrapol.o fluxes.o outeqm.o
