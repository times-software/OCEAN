include ../Makefile.arch

SCRIPTS = dendip.pl par_dendip.pl bvecs.pl gvecs2.pl ofermi.pl OBF_wfn.pl qe_dendip.pl qe_data_file.pl jdftx_dendip.pl
FILESF90 = nelectron.x rhoofg.x wfconvert.x
WFCONVERTGW_OBJS = wfconvert_gw.o headpar.o grabwf.o wfsubs.o
FFT_WRAPPER = ../Modules/FFT_wrapper.o -I../Modules ../Modules/jlmfft.o

ifdef __HAVE_IOTK
	WFCONVERT_OBJS = headpar.o grabwf.o wfsubs.o qehead.o qe_grabwf.o
	INCL_IOTK = -I$(ESPRESSO_DIR)/iotk/src
	LIB_IOTK = $(ESPRESSO_DIR)/iotk/src/libiotk.a
	OPTIONS += -D__HAVE_ESPRESSO
else
	WFCONVERT_OBJS = headpar.o grabwf.o wfsubs.o qe_grabwf_noiotk.o
endif


.SUFFIXES: .f .x .f90

all: wfconvert.x nelectron.x rhoofg.x 
	cd src; make

install:
	cp $(SCRIPTS) $(FILESF90) $(INSTDIR)
	cd src; make install

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(FILESF90); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

clean: 
	rm -f *.o *.x *__genmod.mod *__genmod.f90
	cd src; make clean

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f

rhoofg.x: rho2.f90 
	$(FC) $(OPTIONS) rho2.f90 $(FFT_WRAPPER) $(FFTWL) -o rhoofg.x 

nelectron.x: nelectron.f90
	$(FC) $(OPTIONS) nelectron.f90 -o nelectron.x

wfconvert.x: $(WFCONVERT_OBJS) wfconvert.o
	$(MPIFORT) $(OPTIONS) $(WFCONVERT_OBJS) wfconvert.o -o wfconvert.x $(LIB_IOTK) 

wfconvert_gw.x: $(WFCONVERTGW_OBJS)
	$(FC) $(OPTIONS) $(WFCONVERTGW_OBJS) -o wfconvert_gw.x

wfconvert.o: wfconvert2.f90 $(WFCONVERT_OBJS)
	$(FC) $(OPTIONS) -c wfconvert2.f90 -o wfconvert.o

wfconvert_gw.o: wfconvert_gw.f90
	$(FC) $(OPTIONS) -c wfconvert_gw.f90 -o wfconvert_gw.o

headpar.o: headpar.f90
	 $(FC) $(OPTIONS) -c headpar.f90 -o headpar.o

grabwf.o: grabwf.f90
	$(FC) $(OPTIONS) -c grabwf.f90 -o grabwf.o

wfsubs.o: wfsubs.f90
	$(FC) $(OPTIONS) -c wfsubs.f90 -o wfsubs.o

qehead.o: qehead.f90 
	$(FC) $(OPTIONS) -c qehead.f90 -o qehead.o $(INCL_IOTK)

qe_grabwf.o: qe_grabwf.f90
	$(FC) $(OPTIONS) -c qe_grabwf.f90 -o qe_grabwf.o $(INCL_IOTK) 

qe_grabwf_noiotk.o: qe_grabwf_noiotk.f90
	$(FC) $(OPTIONS) -c qe_grabwf_noiotk.f90 -o qe_grabwf_noiotk.o
