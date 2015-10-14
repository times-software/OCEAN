include ../Makefile.arch

ifndef ESPRESSO_DIR
all:
	@echo "Your build does not have QE support."
install:
	@echo "Your build does not have QE support."
clean:
	rm -f *.o *.x

else

include $(ESPRESSO_DIR)/make.sys
.SUFFIXES :


QEOBJS = $(QEXML) $(ESPRESSO_DIR)/Modules/libqemod.a $(ESPRESSO_DIR)/flib/ptools.a \
         $(ESPRESSO_DIR)/flib/flib.a $(ESPRESSO_DIR)/clib/clib.a $(ESPRESSO_DIR)/iotk/src/libiotk.a


PP = $(ESPRESSO_DIR)/PP/src
LIBOBJS = $(ESPRESSO_DIR)/flib/ptools.a $(ESPRESSO_DIR)/flib/flib.a $(ESPRESSO_DIR)/clib/clib.a \
					$(ESPRESSO_DIR)/iotk/src/libiotk.a 

SCRIPTS = 
FILESF90 = qe_wfconvert.x qeband.x qebocc.x
WFCONVERT_OBJS = wfconvert.o headpar.o grabwf.o wfsubs.o



all: $(FILESF90)

install:
	cp $(SCRIPTS) $(FILESF90) $(INSTDIR)

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(FILESF90); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

clean:
	rm -f *.o *.x


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
endif
