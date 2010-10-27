SCRIPTS = dendip.pl par_dendip.pl bvecs.pl gvecs2.pl ofermi.pl
FILESF90 = wfconvert.x nelectron.x rhoofg.x
WFCONVERT_OBJS = wfconvert.o headpar.o grabwf.o wfsubs.o


.SUFFIXES: .f .x .f90

all: wfconvert.x nelectron.x rhoofg.x

install:
	cp $(SCRIPTS) $(FILESF90) $(INSTDIR)

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(FILESF90); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f

rhoofg.x: rho2.f90
	$(FC) $(OPTIONS) rho2.f90 $(FFTWI) $(FFTWL) -o rhoofg.x

nelectron.x: nelectron.f90
	$(FC) $(OPTIONS) nelectron.f90 -o nelectron.x

wfconvert.x: $(WFCONVERT_OBJS)
	$(FC) $(OPTIONS) $(WFCONVERT_OBJS) -o wfconvert.x

wfconvert.o: wfconvert2.f90
	$(FC) $(OPTIONS) -c wfconvert2.f90 -o wfconvert.o

headpar.o: headpar.f90
	 $(FC) $(OPTIONS) -c headpar.f90 -o headpar.o

grabwf.o: grabwf.f90
	$(FC) $(OPTIONS) -c grabwf.f90 -o grabwf.o

wfsubs.o: wfsubs.f90
	$(FC) $(OPTIONS) -c wfsubs.f90 -o wfsubs.o
