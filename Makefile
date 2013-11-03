# Master Makefile to compile and install all the programs required to run
# ai2nbse. Please, don't change anything here, add all definitions to the
# Makefile.arch file

include Makefile.arch

SCRIPTS = 
AUX = Version Header

all:
	cd Common;       make "FC = $(FC)" "OPTIONS = $(OPTIONS)"
	cd ABINIT;       make "FC = $(FC)" "OPTIONS = $(OPTIONS)"
	cd ESPRESSO;     ./configure
	cd ESPRESSO;     make pw pp
	cd PREP;         make "FC = $(FC)" "OPTIONS = $(OPTIONS)" "FFTWI = $(FFTWI)" "FFTWL = $(FFTWL)"
	cd QEPREP;       make "FC = $(FC)" "OPTIONS = $(OPTIONS)"
	cd PAW;          make "COMP = $(FC)" "OPTS = $(OPTIONS)"
	cd SCREEN;       make "COMP = $(FC)" "OPTS = $(OPTIONS)"
	cd CNBSE;        make "COMP = $(FC)" "OPTS = $(OPTIONS)"
	cd zbridge;      make "COMP = $(FC)" "OPTS = $(OPTIONS)"

espresso:
	cd ESPRESSO      ./configure
	cd ESPRESSO      make pw pp

prep:
	cd PREP;         make "FC = $(FC)" "OPTIONS = $(OPTIONS)" "FFTWI = $(FFTWI)" "FFTWL = $(FFTWL)"

qeprep:
	cd QEPREP;       make "FC = $(FC)" "OPTIONS = $(OPTIONS)"

bridge:
	cd zbridge;      make "COMP = $(FC)" "OPTS = $(OPTIONS)"

screen:
	cd SCREEN;       make "COMP = $(FC)" "OPTS = $(OPTIONS)"

cnbse:
	cd CNBSE;        make "COMP = $(FC)" "OPTS = $(OPTIONS)"

clean:
	rm -f $(INSTDIR)/*.x $(INSTDIR)/*.script
	rm -r $(INSTDIR)
	mkdir $(INSTDIR)
	cd Common;       make clean
	cd ABINIT;       make clean
	cd ESPRESSO;     make clean
	cd PREP;         make clean
	cd QEPREP;       make clean
	cd SCREEN;       make clean
	cd PAW;          make clean
	cd CNBSE;        make clean 
	cd zbridge;      make clean

install:
	mkdir $(INSTDIR)
	cp $(SCRIPTS) $(AUX) $(INSTDIR)
	cd Common;       make "INSTDIR = $(INSTDIR)" install
	cd ABINIT;       make "INSTDIR = $(INSTDIR)" install
	cd PREP;         make "INSTDIR = $(INSTDIR)" install
	cd QEPREP;       make "INSTDIR = $(INSTDIR)" install
	cd SCREEN;       make "INSTDIR = $(INSTDIR)" install
	cd PAW;          make "INSTDIR = $(INSTDIR)" install
	cd CNBSE;        make "INSTDIR = $(INSTDIR)" install
	cd zbridge;      make "INSTDIR = $(INSTDIR)" install
#	chmod u+x $(INSTDIR)/*.script
	chmod u+x $(INSTDIR)/*.pl
#	ln -sf $(ABINIT) $(INSTDIR)/abinit
#	ln -sf $(CUT3D) $(INSTDIR)/cut3d
#	cp ESPRESSO/PW/pw.x $(INSTDIR)
#	cp ESPRESSO/PP/pp.x $(INSTDIR)
	cp ESPRESSO/PP/qeband.x $(INSTDIR)
	cp ESPRESSO/PP/qebocc.x $(INSTDIR)
	cp ESPRESSO/PP/wfconvert.x $(INSTDIR)

instdev:
	for F in $(SCRIPTS) $(AUX); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	cd Common;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd ABINIT;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd PREP;         make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd QEPREP;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd CNBSE;        make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd SCREEN;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
#	chmod u+x $(INSTDEVDIR)/*.script
	chmod u+x $(INSTDEVDIR)/*.pl
