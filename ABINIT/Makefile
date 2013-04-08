FILES   = avec.x 
SCRIPTS = AbinitDriver.pl NewAbinitDriver.pl pp.pl par_ab2.pl


.SUFFIXES: .f .x

all: $(FILES) kgen2.x makecoords.x makeatompp.x makeacell.x kgen_qe.x

clean:
	rm -f $(FILES)
	rm -f kgen.x

install:
	cp $(SCRIPTS) $(FILES) $(AUX) kgen2.x $(INSTDIR)

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(FILES); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(AUX); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	ln -fs $(PWD)/kgen.x $(INSTDEVDIR)

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f

kgen2.x: kgen2.f90
	$(FC) $(OPTIONS) kgen2.f90 -o kgen2.x

makecoords.x: makecoords.f90
	$(FC) $(OPTIONS)   makecoords.f90 -o makecoords.x

makeatompp.x: makeatompp.f90
	 $(FC) $(OPTIONS) makeatompp.f90 -o makeatompp.x

makeacell.x: makeacell.f90
	 $(FC) $(OPTIONS) makeacell.f90 -o makeacell.x

kgen_qe.x: kgen_qe.f90
	$(FC) $(OPTIONS) kgen_qe.f90 -o kgen_qe.x
