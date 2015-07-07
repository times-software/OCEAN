include ../Makefile.arch

SCRIPTS = dendip.pl par_dendip.pl bvecs.pl gvecs2.pl ofermi.pl OBF_wfn.pl
FILESF90 = nelectron.x rhoofg.x wfconvert.x
WFCONVERT_OBJS = wfconvert.o headpar.o grabwf.o wfsubs.o 
WFCONVERTGW_OBJS = wfconvert_gw.o headpar.o grabwf.o wfsubs.o
FFT_WRAPPER = ../zbridge/zobj/FFT_wrapper.o -I../zbridge/ ../zbridge/zobj/jlmfft.o

.SUFFIXES: .f .x .f90

all: wfconvert.x nelectron.x rhoofg.x 

install:
	cp $(SCRIPTS) $(FILESF90) $(INSTDIR)

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(FILESF90); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

clean: 
	rm -f *.o *.x

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f

rhoofg.x: rho2.f90 
	$(FC) $(OPTIONS) rho2.f90 $(FFT_WRAPPER) -o rhoofg.x 

nelectron.x: nelectron.f90
	$(FC) $(OPTIONS) nelectron.f90 -o nelectron.x

wfconvert.x: $(WFCONVERT_OBJS)
	$(FC) $(OPTIONS) $(WFCONVERT_OBJS) -o wfconvert.x

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
