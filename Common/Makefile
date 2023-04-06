include ../Makefile.arch
FILES   =
SCRIPTS = ocean.pl par_ocean.pl parse defaults.pl defaults.h parse.pl edges.pl structure.pl psp_parser.pl \
          oparse.json oparse.type.json parseJSON.pl OCEANcompare.pl photon.pl
AUX = oparse.h Pquadrature specpnt sphpts EvenQuadHalf.txt
SPECPNT = specpnt.5 specpnt.6 specpnt.7 specpnt.8 specpnt.9 specpnt.10 \
          specpnt.11 specpnt.12 specpnt.13 specpnt.14 specpnt.15 \
          specpnt.16 specpnt.17 specpnt.18 specpnt.19 specpnt.20 \
          specpnt.21 specpnt.22 specpnt.23 specpnt.24 specpnt.25 \
					specpnt.26 specpnt.27 specpnt.28 specpnt.29 specpnt.30 \
					specpnt.31 specpnt.32 specpnt.33 specpnt.34 specpnt.35 \
          specpnt.36 specpnt.37 specpnt.38 specpnt.39 specpnt.40 

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
