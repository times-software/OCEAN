include ../Makefile.arch
FILES   =
SCRIPTS = ocean.pl par_ocean.pl parse defaults.pl defaults.h parse.pl edges.pl \
          psp_tools/adder.py psp_tools/finder.py psp_tools/interactivefinder.py
AUX = oparse.h Pquadrature specpnt sphpts
LEBDEV = specpnt.5 specpnt.6 specpnt.7 specpnt.8 specpnt.9 specpnt.10 specpnt.15 specpnt.20

.SUFFIXES: .f .x

all: $(FILES)

clean:
	rm -f $(FILES)

install:
	cp $(SCRIPTS) $(FILES) $(AUX) $(LEBDEV) $(INSTDIR)

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(FILES); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(AUX); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(LEBDEV); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f
