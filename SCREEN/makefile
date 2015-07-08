include ../Makefile.arch

COMP = $(FC)
OPTS = $(OPTIONS)

all:
	cd zxipp; make "COMP = $(COMP)" "OPTS = $(OPTS)"
	cd zf90;  make "COMP = $(COMP)" "OPTS = $(OPTS)"

install:
	cd zxipp; make "INSTDIR = $(INSTDIR)" install
	cd zf90;  make "INSTDIR = $(INSTDIR)" install
	cp screen.pl OBF_screen.pl OBF_screen_multi.pl core_shift.pl $(INSTDIR)/ 

clean:
	cd zxipp; make clean
	cd zf90; make clean

