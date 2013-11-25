# Master Makefile to compile and install all the programs required to run
# ai2nbse. Please, don't change anything here, add all definitions to the
# Makefile.arch file

include Makefile.arch

SCRIPTS = 
AUX = Version Header

all:
	cd Common;       make "FC = $(FC)" "OPTIONS = $(OPTIONS)"
	cd ABINIT;       make "FC = $(FC)" "OPTIONS = $(OPTIONS)"
	cd PREP;       make "FC = $(FC)" "OPTIONS = $(OPTIONS)"
	cd PAW;          make "COMP = $(FC)" "OPTS = $(OPTIONS)"
	cd SCREEN;       make "COMP = $(FC)" "OPTS = $(OPTIONS)"
	cd CNBSE;        make "COMP = $(FC)" "OPTS = $(OPTIONS)"
	cd zbridge;      make "COMP = $(FC)" "OPTS = $(OPTIONS)"

clean:
	cd Common;       make clean
	cd ABINIT;       make clean
	cd PREP;      make clean
	cd SCREEN;       make clean
	cd PAW;   make clean
	cd CNBSE  make clean 
	cd zbridge;  make clean
install:
	mkdir -p $(INSTDIR)
	cp $(SCRIPTS) $(AUX) $(INSTDIR)
	cd Common;       make "INSTDIR = $(INSTDIR)" install
	cd ABINIT;       make "INSTDIR = $(INSTDIR)" install
	cd PREP;       make "INSTDIR = $(INSTDIR)" install
	cd SCREEN;        make "INSTDIR = $(INSTDIR)" install
	cd PAW;          make "INSTDIR = $(INSTDIR)" install
	cd CNBSE;        make "INSTDIR = $(INSTDIR)" install
	cd zbridge;      make "INSTDIR = $(INSTDIR)" install
	chmod u+x $(INSTDIR)/*.pl

instdev:
	for F in $(SCRIPTS) $(AUX); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	cd Common;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd ABINIT;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd PREP;      make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd CNBSE;         make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd SCREEN;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	chmod u+x $(INSTDEVDIR)/*.pl
