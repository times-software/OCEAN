# Master Makefile to compile and install all the programs required to run
# OCEAN. Please, don't change anything here, add all definitions to the
# Makefile.arch file

include Makefile.arch

SCRIPTS = 
AUX = Version Header
AUX_EXE = $(ABINIT_EXE) $(CUT3D_EXE) $(PW_EXE) $(PP_EXE)

all: 
	perl append_git.pl
	cd Common;       make 
	cd zbridge;      make 
	cd Modules;      make
	cd DFT;          make 
	cd PREP;         make 
	cd OPF;          make 
	cd SCREEN;       make 
	cd CNBSE;        make 
	cd INTERPOLATE;  make
	cd OCEAN2;			 make
	cd POST;         make

clean:
	cd Common;       make clean
	cd Modules;      make clean
	cd DFT;          make clean
	cd PREP;         make clean
	cd SCREEN;       make clean
	cd OPF;          make clean
	cd CNBSE;        make clean 
	cd zbridge;      make clean
	cd INTERPOLATE;  make clean
	cd OCEAN2;       make clean

install:
	mkdir -p $(INSTDIR)
	cp $(SCRIPTS) $(AUX) $(INSTDIR)
	cp $(AUX_EXE) $(INSTDIR)
	cd Common;       make install
	cd DFT;          make install
	cd PREP;         make install
	cd SCREEN;       make install
	cd OPF;          make install
	cd CNBSE;        make install
	cd zbridge;      make install
	cd INTERPOLATE;  make install
	cd OCEAN2;       make install
	chmod u+x $(INSTDIR)/*.pl
	cd POST;				 make install

instdev:
	for F in $(SCRIPTS) $(AUX); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	cd Common;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd DFT;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd PREP;      make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd CNBSE;         make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	cd SCREEN;       make "INSTDEVDIR = $(INSTDEVDIR)" instdev
	chmod u+x $(INSTDEVDIR)/*.pl

