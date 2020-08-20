include ../Makefile.arch
SCRIPTS = extractPsp.pl inputDatabase.pl pspDatabase.pl
AUX = json/ONCVPSP-PBE-PDv0.4-standard.json json/ONCVPSP-PBE-PDv0.4-stringent.json 

make:

install:
	cp $(SCRIPTS) $(FILES) $(AUX) $(INSTDIR)

instdev:
	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
	for F in $(AUX); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

