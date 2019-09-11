include ../Makefile.arch
FILES   =
SCRIPTS = ocean.pl par_ocean.pl parse defaults.pl defaults.h parse.pl edges.pl structure.pl psp_parser.pl
AUX = oparse.h Pquadrature specpnt sphpts EvenQuadHalf.txt
SPECPNT = specpnt.5 specpnt.6 specpnt.7 specpnt.8 specpnt.9 specpnt.10 \
          specpnt.11 specpnt.12 specpnt.13 specpnt.14 specpnt.15 \
          specpnt.16 specpnt.17 specpnt.18 specpnt.19 specpnt.20 \
          specpnt.21 specpnt.22 specpnt.23 specpnt.24 

.SUFFIXES: .f .x

all: $(FILES)

clean:
	rm -f $(FILES)

install:
	cp $(SCRIPTS) $(FILES) $(AUX) $(SPECPNT) $(INSTDIR)

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(FILES); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(AUX); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(SPECPNT); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f
