include ../Makefile.arch

SCRIPTS = dendip.pl par_dendip.pl bvecs.pl gvecs2.pl ofermi.pl OBF_wfn.pl qe_dendip.pl
FILESF90 = nelectron.x rhoofg.x wfconvert.x
WFCONVERT_OBJS = wfconvert.o headpar.o grabwf.o wfsubs.o qehead.o qe_grabwf.o
WFCONVERTGW_OBJS = wfconvert_gw.o headpar.o grabwf.o wfsubs.o
FFT_WRAPPER = ../zbridge/zobj/FFT_wrapper.o -I../zbridge/ ../zbridge/zobj/jlmfft.o
INCL_IOTK = -I$(ESPRESSO_DIR)/iotk/src
LIB_IOTK = $(ESPRESSO_DIR)/iotk/src/libiotk.a


.SUFFIXES: .f .x .f90

all: wfconvert.x nelectron.x rhoofg.x 

install:
	cp $(SCRIPTS) $(FILESF90) $(INSTDIR)

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(FILESF90); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

clean: 
	rm -f *.o *.x *__genmod.mod *__genmod.f90

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f

rhoofg.x: rho2.f90 
	$(FC) $(OPTIONS) rho2.f90 $(FFT_WRAPPER) $(FFTWL) -o rhoofg.x 

nelectron.x: nelectron.f90
	$(FC) $(OPTIONS) nelectron.f90 -o nelectron.x

wfconvert.x: $(WFCONVERT_OBJS)
	$(MPIFORT) $(OPTIONS) $(WFCONVERT_OBJS) -o wfconvert.x $(LIB_IOTK)

wfconvert_gw.x: $(WFCONVERTGW_OBJS)
	$(FC) $(OPTIONS) $(WFCONVERTGW_OBJS) -o wfconvert_gw.x

wfconvert.o: wfconvert2.f90
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
