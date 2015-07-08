include ../Makefile.arch

COMP = $(FC)
OPTS = $(OPTIONS)

all:
	cd ppmod; make "COMP = $(COMP)" "OPTS = $(OPTS)"
	cd zatom; make "COMP = $(COMP)" "OPTS = $(OPTS)"

install:
	cd ppmod; make "INSTDIR = $(INSTDIR)" install
	cd zatom; make "INSTDIR = $(INSTDIR)" install
	cp paw.pl $(INSTDIR)/ 

clean:
	cd ppmod; make clean
	cd zatom; make clean

