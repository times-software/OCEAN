# Master Makefile to compile and install all the programs required to run
# OCEAN. Please, don't change anything here, add all definitions to the
# Makefile.arch file

include Makefile.arch

SCRIPTS = 
AUX = Version Header

all:
	cd Common;       make 
	cd DFT;          make 
	cd PREP;         make 
	cd PAW;          make 
	cd SCREEN;       make 
	cd CNBSE;        make 
	cd zbridge;      make 

clean:
	cd Common;       make clean
	cd DFT;          make clean
	cd PREP;         make clean
	cd SCREEN;       make clean
	cd PAW;          make clean
	cd CNBSE         make clean 
	cd zbridge;      make clean
install:
	mkdir -p $(INSTDIR)
	cp $(SCRIPTS) $(AUX) $(INSTDIR)
	cd Common;       make install
	cd DFT;          make install
	cd PREP;         make install
	cd SCREEN;       make install
	cd PAW;          make install
	cd CNBSE;        make install
	cd zbridge;      make install
	chmod u+x $(INSTDIR)/*.pl

instdev:
	for F in $(SCRIPTS) $(AUX); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	cd Common;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd DFT;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd PREP;      make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd CNBSE;         make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd SCREEN;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	chmod u+x $(INSTDEVDIR)/*.pl
