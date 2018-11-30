include ../Makefile.arch

COMP = $(FC)
OPTS = $(OPTIONS) $(LEGACY_REAL)

all:
	cd ppmod; make "COMP = $(COMP)" "OPTS = $(OPTS)"
	cd zatom; make "COMP = $(COMP)" "OPTS = $(OPTS)"

install:
	cd ppmod; make "INSTDIR = $(INSTDIR)" install
	cd zatom; make "INSTDIR = $(INSTDIR)" install
	cp opf.pl $(INSTDIR)/ 
	cp validate_opts.pl $(INSTDIR)/

clean:
	cd ppmod; make clean
	cd zatom; make clean

