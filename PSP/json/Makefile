include ../../Makefile.arch

AUX = json/ONCVPSP-PBE-PDv0.4-standard.json json/ONCVPSP-PBE-PDv0.4-stringent.json \
			json/ONCVPSP-PW-PDv0.4-standard.json json/ONCVPSP-PW-PDv0.4-stringent.json \
			json/ONCVPSP-PBESOL-PDv0.4-standard.json json/ONCVPSP-PBESOL-PDv0.4-stringent.json
PSP: ONCVPSP-PW-PDv0.4.upf.json ONCVPSP-PW-PDv0.4.psp8.json \
     ONCVPSP-PBE-PDv0.4.upf.json ONCVPSP-PBE-PDv0.4.psp8.json \
     ONCVPSP-PBESOL-PDv0.4.upf.json ONCVPSP-PBESOL-PDv0.4.psp8.json


all: $(PSP)

install:
	cp *.json $(INSTDIR)

database: $(PSP)

ONCVPSP-PW-PDv0.4.upf.json: ONCVPSP-PW-PDv0.4-standard.json ONCVPSP-PW-PDv0.4-stringent.json
	$(INSTDIR)/pspDatabase.pl $(INSTDIR)/oncvpsp.x ONCVPSP-PW-PDv0.4-standard.json ONCVPSP-PW-PDv0.4-stringent.json

ONCVPSP-PW-PDv0.4.psp8.json: ONCVPSP-PW-PDv0.4-standard.json ONCVPSP-PW-PDv0.4-stringent.json
	$(INSTDIR)/pspDatabase.pl $(INSTDIR)/oncvpsp.x ONCVPSP-PW-PDv0.4-standard.json ONCVPSP-PW-PDv0.4-stringent.json

ONCVPSP-PBE-PDv0.4.upf.json: ONCVPSP-PBE-PDv0.4-standard.json ONCVPSP-PBE-PDv0.4-stringent.json
	$(INSTDIR)/pspDatabase.pl $(INSTDIR)/oncvpsp.x ONCVPSP-PBE-PDv0.4-standard.json ONCVPSP-PBE-PDv0.4-stringent.json

ONCVPSP-PBE-PDv0.4.psp8.json: ONCVPSP-PBE-PDv0.4-standard.json ONCVPSP-PBE-PDv0.4-stringent.json
	$(INSTDIR)/pspDatabase.pl $(INSTDIR)/oncvpsp.x ONCVPSP-PBE-PDv0.4-standard.json ONCVPSP-PBE-PDv0.4-stringent.json

ONCVPSP-PBESOL-PDv0.4.upf.json: ONCVPSP-PBESOL-PDv0.4-standard.json ONCVPSP-PBESOL-PDv0.4-stringent.json
	$(INSTDIR)/pspDatabase.pl $(INSTDIR)/oncvpsp.x ONCVPSP-PBESOL-PDv0.4-standard.json ONCVPSP-PBESOL-PDv0.4-stringent.json

ONCVPSP-PBESOL-PDv0.4.psp8.json: ONCVPSP-PBESOL-PDv0.4-standard.json ONCVPSP-PBESOL-PDv0.4-stringent.json
	$(INSTDIR)/pspDatabase.pl $(INSTDIR)/oncvpsp.x ONCVPSP-PBESOL-PDv0.4-standard.json ONCVPSP-PBESOL-PDv0.4-stringent.json
