FILES   = avec.x 
SCRIPTS = AbinitDriver.pl NewAbinitDriver.pl pp.pl p


.SUFFIXES: .f .x

all: $(FILES) kgen.x

clean:
	rm -f $(FILES)
	rm -f kgen.x

install:
	cp $(SCRIPTS) $(FILES) $(AUX) kgen.x $(INSTDIR)

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(FILES); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(AUX); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	ln -fs $(PWD)/kgen.x $(INSTDEVDIR)

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f

kgen.x: kgen.f90
	$(FC) $(OPTIONS) kgen.f90 -o kgen.x
