include ../Makefile.arch

include $(ESPRESSO_DIR)/make.sys


PP = $(ESPRESSO_DIR)/PP
LIBOBJS = $(ESPRESSO_DIR)/flib/ptools.a $(ESPRESSO_DIR)/flib/flib.a $(ESPRESSO_DIR)/clib/clib.a \
					$(ESPRESSO_DIR)/iotk/src/libiotk.a 

SCRIPTS = dendip.pl par_dendip.pl bvecs.pl gvecs2.pl ofermi.pl QespressoDriver.pl \
          celldm.py converter.py qe_dendip.pl 
FILESF90 = makecoords.x makeatompp.x makeacell.x nelectron.x wfconvert.x qeband.x qebocc.x
# we also want rhoofg.x, but are having trouble linking the FFTW libraries
#FILESF90 = makecoords.x makeatompp.x makeacell.x nelectron.x rhoofg.x
#FILESF90 = makecoords.x makeatompp.x makeacell.x wfconvert.x nelectron.x rhoofg.x
WFCONVERT_OBJS = wfconvert.o headpar.o grabwf.o wfsubs.o



all: makecoords.x makeatompp.x makeacell.x nelectron.x wfconvert.x qeband.x qebocc.x
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


makecoords.x: makecoords.f90
	$(FC) $(OPTIONS)  makecoords.f90 -o makecoords.x

makeatompp.x: makeatompp.f90
	$(FC) $(OPTIONS) makeatompp.f90 -o makeatompp.x

makeacell.x: makeacell.f90
	$(FC) $(OPTIONS) makeacell.f90 -o makeacell.x

rhoofg.x: rho2.f90
	$(FC) $(OPTIONS) rho2.f90 $(FFTWI) $(FFTWL) -o rhoofg.x

nelectron.x: nelectron.f90
	$(FC) $(OPTIONS) nelectron.f90 -o nelectron.x

headpar.o: ../PREP/headpar.f90
	 $(FC) $(OPTIONS) -c ../PREP/headpar.f90 -o headpar.o

grabwf.o: ../PREP/grabwf.f90
	$(FC) $(OPTIONS) -c ../PREP/grabwf.f90 -o grabwf.o

wfsubs.o: ../PREP/wfsubs.f90
	$(FC) $(OPTIONS) -c ../PREP/wfsubs.f90 -o wfsubs.o

wfconvert.o: wfconvert.f90
	$(FC) $(OPTIONS) -c wfconvert.f90 -o wfconvert.o

qehead2.o: qehead2.f90
	$(FC) $(OPTIONS) -c qehead2.f90 -o qehead2.o -I$(PP)

qenpw2.o : qenpw2.f90
	$(FC) $(OPTIONS) -c qenpw2.f90 -o qenpw2.o -I$(PP)

qebocc.o: qebocc.f90
	$(FC) $(OPTIONS) -c qebocc.f90 -o qebocc.o -I$(PP)

qeband.o: qeband.f90
	$(FC) $(OPTIONS) -c qeband.f90 -o qeband.o -I$(PP)

wfconvert.x: $(PP)/qexml.o wfconvert.o qehead2.o grabwf.o qenpw2.o wfsubs.o 
	$(LD) $(LDFLAGS) -o $@ wfconvert.o qehead2.o grabwf.o qenpw2.o wfsubs.o $(PP)/qexml.o $(ESPRESSO_DIR)/Modules/libqemod.a -I$(PP) $(LIBOBJS)

qeband.x: $(PP)/qexml.o qeband.o
	$(LD) $(LDFLAGS) -o $@ qeband.o $(PP)/qexml.o $(ESPRESSO_DIR)/Modules/libqemod.a $(LIBOBJS)

qebocc.x: $(PP)/qexml.o qebocc.o
	$(LD) $(LDFLAGS) -o $@ qebocc.o $(PP)/qexml.o $(ESPRESSO_DIR)/Modules/libqemod.a $(LIBOBJS)

