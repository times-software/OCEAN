include ../Makefile.arch

SCRIPTS = QespressoDriver.pl \
          celldm.py converter.py qe_dendip.pl 
FILESF90 = makecoords.x makeatompp.x makeacell.x 
# we also want rhoofg.x, but are having trouble linking the FFTW libraries
#FILESF90 = makecoords.x makeatompp.x makeacell.x nelectron.x rhoofg.x
#FILESF90 = makecoords.x makeatompp.x makeacell.x wfconvert.x nelectron.x rhoofg.x
WFCONVERT_OBJS = wfconvert.o headpar.o grabwf.o wfsubs.o


.SUFFIXES: .f .x .f90

all: makecoords.x makeatompp.x makeacell.x 
# find the FFTW libraries so that you can also build rhoofg.x
#all: makecoords.x makeatompp.x makeacell.x nelectron.x rhoofg.x
#all: makecoords.x makeatompp.x makeacell.x wfconvert.x nelectron.x rhoofg.x

install:
	cp $(SCRIPTS) $(FILESF90) $(INSTDIR)

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(FILESF90); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

clean:
	rm -f *.o *.x

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f

makecoords.x: makecoords.f90
	$(FC) $(OPTIONS)  makecoords.f90 -o makecoords.x

makeatompp.x: makeatompp.f90
	$(FC) $(OPTIONS) makeatompp.f90 -o makeatompp.x

makeacell.x: makeacell.f90
	$(FC) $(OPTIONS) makeacell.f90 -o makeacell.x

rhoofg.x: rho2.f90
	$(FC) $(OPTIONS) rho2.f90 $(FFTWI) $(FFTWL) -o rhoofg.x

headpar.o: headpar.f90
	 $(FC) $(OPTIONS) -c headpar.f90 -o headpar.o

grabwf.o: grabwf.f90
	$(FC) $(OPTIONS) -c grabwf.f90 -o grabwf.o

wfsubs.o: wfsubs.f90
	$(FC) $(OPTIONS) -c wfsubs.f90 -o wfsubs.o
