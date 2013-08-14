FILES   =
SCRIPTS = ocean.pl par_ocean.pl parse defaults.pl defaults.h
AUX = oparse.h Pquadrature specpnt sphpts

.SUFFIXES: .f .x

all: $(FILES)

clean:
	rm -f $(FILES)

install:
	cp $(SCRIPTS) $(FILES) $(AUX) $(INSTDIR)

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(FILES); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(AUX); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f
