include ../Makefile.arch
SCRIPTS = extractPsp.pl inputDatabase.pl pspDatabase.pl

make:

install:
	cp $(SCRIPTS) $(FILES) $(INSTDIR)
	cd json; make install

database:
	cd json; make
