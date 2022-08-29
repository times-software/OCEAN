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
	cp upf2shirley.pl $(INSTDIR)/
	cp validate_opts.pl $(INSTDIR)/
	@if [ -f '$(ONCVPSP_DIR)/src/oncvpsp.x' ] ; \
	then \
     cp $(ONCVPSP_DIR)/src/oncvpsp.x $(INSTDIR) ; \
	fi

clean:
	cd ppmod; make clean
	cd zatom; make clean

