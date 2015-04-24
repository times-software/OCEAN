include ../Makefile.arch
#### Select your QE version here

#### QE 5.1 
#QEXML = $(ESPRESSO_DIR)/Modules/qexml.o
#QEMOD = $(ESPRESSO_DIR)/Modules/
#OPTIONS += -D__QE51
#### QE 5.0 or earlier?
QEXML = $(ESPRESSO_DIR)/PP/src/qexml.o
QEMOD = $(ESPRESSO_DIR)/PP/src





include $(ESPRESSO_DIR)/make.sys
.SUFFIXES :


# QE 5.X
#../flib/ptools.a ../flib/flib.a ../clib/clib.a ../iotk/src/libiotk.a
QEOBJS = $(QEXML) $(ESPRESSO_DIR)/Modules/libqemod.a $(ESPRESSO_DIR)/flib/ptools.a \
         $(ESPRESSO_DIR)/flib/flib.a $(ESPRESSO_DIR)/clib/clib.a $(ESPRESSO_DIR)/iotk/src/libiotk.a


PP = $(ESPRESSO_DIR)/PP/src
LIBOBJS = $(ESPRESSO_DIR)/flib/ptools.a $(ESPRESSO_DIR)/flib/flib.a $(ESPRESSO_DIR)/clib/clib.a \
					$(ESPRESSO_DIR)/iotk/src/libiotk.a 

SCRIPTS = dendip.pl par_dendip.pl bvecs.pl gvecs2.pl ofermi.pl QespressoDriver.pl \
          celldm.py converter.py qe_dendip.pl 
FILESF90 = makecoords.x makeatompp.x makeacell.x nelectron.x qe_wfconvert.x qeband.x qebocc.x
# we also want rhoofg.x, but are having trouble linking the FFTW libraries
#FILESF90 = makecoords.x makeatompp.x makeacell.x nelectron.x rhoofg.x
#FILESF90 = makecoords.x makeatompp.x makeacell.x wfconvert.x nelectron.x rhoofg.x
WFCONVERT_OBJS = wfconvert.o headpar.o grabwf.o wfsubs.o



all: makecoords.x makeatompp.x makeacell.x nelectron.x qe_wfconvert.x qeband.x qebocc.x
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
	$(FC) $(OPTIONS) makecoords.f90 -o makecoords.x

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

grabwf.o: grabwf.f90
	$(FC) $(OPTIONS) -cpp -c grabwf.f90 -o grabwf.o -I$(QEMOD)

wfsubs.o: wfsubs.f90
	$(FC) $(OPTIONS) -c wfsubs.f90 -o wfsubs.o

wfconvert.o: wfconvert.f90
	$(FC) $(OPTIONS) -cpp -c wfconvert.f90 -o wfconvert.o

qehead2.o: qehead2.f90
	$(FC) $(OPTIONS) -c qehead2.f90 -o qehead2.o -I$(QEMOD)

qenpw2.o : qenpw2.f90
	$(FC) $(OPTIONS) -c qenpw2.f90 -o qenpw2.o -I$(QEMOD)

qebocc.o: qebocc.f90
	$(FC) $(OPTIONS) -cpp -c qebocc.f90 -o qebocc.o -I$(QEMOD)

qeband.o: qeband.f90
	$(FC) $(OPTIONS) -cpp -c qeband.f90 -o qeband.o -I$(QEMOD)

qe_wfconvert.x: $(QEOBJS) wfconvert.o qehead2.o grabwf.o qenpw2.o wfsubs.o 
	$(LD) $(LDFLAGS) -o $@ wfconvert.o qehead2.o grabwf.o qenpw2.o wfsubs.o $(QEOBJS) -I$(QEMOD) $(LIBOBJS)

qeband.x: $(QEOBJS) qeband.o
	$(LD) $(LDFLAGS) -o $@ qeband.o $(QEOBJS) $(LIBOBJS)

qebocc.x: $(QEOBJS) qebocc.o
	$(LD) $(LDFLAGS) -o $@ qebocc.o $(QEOBJS) $(LIBOBJS)

