include ../../Makefile.arch

AUX = json/ONCVPSP-PBE-PDv0.4-standard.json json/ONCVPSP-PBE-PDv0.4-stringent.json \
			json/ONCVPSP-PW-PDv0.4-standard.json json/ONCVPSP-PW-PDv0.4-stringent.json \
			json/ONCVPSP-PBEsol-PDv0.4-standard.json json/ONCVPSP-PBEsol-PDv0.4-stringent.json \
      json/ONCVPSP-PBE-PDv0.5-stringent.json json/ONCVPSP-PBE-PDv0.5-standard.json 
#      json/ONCVPSP-PBEsol-PDv0.5-stringent.json
PSP  = ONCVPSP-PW-PDv0.4.upf.json ONCVPSP-PW-PDv0.4.psp8.json \
     ONCVPSP-PBE-PDv0.4.upf.json ONCVPSP-PBE-PDv0.4.psp8.json \
     ONCVPSP-PBEsol-PDv0.4.upf.json ONCVPSP-PBEsol-PDv0.4.psp8.json \
     ONCVPSP-PBE-PDv0.5.upf.json ONCVPSP-PBE-PDv0.5.psp8.json 
#     ONCVPSP-PBEsol-PDv0.5.upf.json ONCVPSP-PBEsol-PDv0.5.psp8.json


all: $(PSP)

clean:
	rm $(PSP) 

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

ONCVPSP-PBEsol-PDv0.4.upf.json: ONCVPSP-PBEsol-PDv0.4-standard.json ONCVPSP-PBEsol-PDv0.4-stringent.json
	$(INSTDIR)/pspDatabase.pl $(INSTDIR)/oncvpsp.x ONCVPSP-PBEsol-PDv0.4-standard.json ONCVPSP-PBEsol-PDv0.4-stringent.json

ONCVPSP-PBEsol-PDv0.4.psp8.json: ONCVPSP-PBEsol-PDv0.4-standard.json ONCVPSP-PBEsol-PDv0.4-stringent.json
	$(INSTDIR)/pspDatabase.pl $(INSTDIR)/oncvpsp.x ONCVPSP-PBEsol-PDv0.4-standard.json ONCVPSP-PBEsol-PDv0.4-stringent.json

ONCVPSP-PBE-PDv0.5.upf.json: ONCVPSP-PBE-PDv0.5-standard.json ONCVPSP-PBE-PDv0.5-stringent.json
	$(INSTDIR)/pspDatabase.pl $(INSTDIR)/oncvpsp.x ONCVPSP-PBE-PDv0.5-standard.json ONCVPSP-PBE-PDv0.5-stringent.json

ONCVPSP-PBE-PDv0.5.psp8.json: ONCVPSP-PBE-PDv0.5-standard.json ONCVPSP-PBE-PDv0.5-stringent.json
	$(INSTDIR)/pspDatabase.pl $(INSTDIR)/oncvpsp.x ONCVPSP-PBE-PDv0.5-standard.json ONCVPSP-PBE-PDv0.5-stringent.json

#ONCVPSP-PBEsol-PDv0.5.upf.json: ONCVPSP-PBEsol-PDv0.5-stringent.json
#	$(INSTDIR)/pspDatabase.pl $(INSTDIR)/oncvpsp.x ONCVPSP-PBE-PDv0.5-stringent.json

#ONCVPSP-PBEsol-PDv0.5.psp8.json: ONCVPSP-PBEsol-PDv0.5-stringent.json
#	$(INSTDIR)/pspDatabase.pl $(INSTDIR)/oncvpsp.x ONCVPSP-PBEsol-PDv0.5-stringent.json
