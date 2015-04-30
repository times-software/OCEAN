OPTS = -fast  #-YDEALLOC=ALL
#COMP = /usr/local/usr/absoft32/bin/f90
COMP = pgf90
#
CKSOBJ = zobj/cainkset.o zobj/cks.o \
         zobj/enkread.o \
         zobj/getabb.o zobj/getcel.o zobj/getcoeff.o zobj/getphase.o zobj/getprefs.o zobj/getylm.o zobj/getylmfac.o \
         zobj/nbsecoeffs.o zobj/nbseftread.o zobj/nbseprjnproj.o zobj/nbseprjprep.o \
         zobj/seanitup.o zobj/setlm.o zobj/snatch.o
#
all: zexe/cks.x
#
zexe/cks.x: $(CKSOBJ)
	$(COMP) $(OPTS) $(CKSOBJ) -o zexe/cks.x
#
zobj/cainkset.o: zf90/cainkset.f90
	$(COMP) $(OPTS) -c zf90/cainkset.f90 -o zobj/cainkset.o
zobj/cks.o: zf90/cks.f90 zf90/sphsetnx.h.f90 zf90/sphsetx.h.f90
	$(COMP) $(OPTS) -c zf90/cks.f90 -o zobj/cks.o
zobj/enkread.o: zf90/enkread.f90
	$(COMP) $(OPTS) -c zf90/enkread.f90 -o zobj/enkread.o
zobj/getabb.o: zf90/getabb.f90
	$(COMP) $(OPTS) -c zf90/getabb.f90 -o zobj/getabb.o
zobj/getcel.o: zf90/getcel.f90
	$(COMP) $(OPTS) -c zf90/getcel.f90 -o zobj/getcel.o
zobj/getcoeff.o: zf90/getcoeff.f90
	$(COMP) $(OPTS) -c zf90/getcoeff.f90 -o zobj/getcoeff.o
zobj/getphase.o: zf90/getphase.f90
	$(COMP) $(OPTS) -c zf90/getphase.f90 -o zobj/getphase.o
zobj/getprefs.o: zf90/getprefs.f90
	$(COMP) $(OPTS) -c zf90/getprefs.f90 -o zobj/getprefs.o
zobj/getylm.o: zf90/getylm.f90
	$(COMP) $(OPTS) -c zf90/getylm.f90 -o zobj/getylm.o
zobj/getylmfac.o: zf90/getylmfac.f90
	$(COMP) $(OPTS) -c zf90/getylmfac.f90 -o zobj/getylmfac.o
zobj/nbsecoeffs.o: zf90/nbsecoeffs.f90
	$(COMP) $(OPTS) -c zf90/nbsecoeffs.f90 -o zobj/nbsecoeffs.o
zobj/nbseftread.o: zf90/nbseftread.f90
	$(COMP) $(OPTS) -c zf90/nbseftread.f90 -o zobj/nbseftread.o
zobj/nbseprjnproj.o: zf90/nbseprjnproj.f90
	$(COMP) $(OPTS) -c zf90/nbseprjnproj.f90 -o zobj/nbseprjnproj.o
zobj/nbseprjprep.o: zf90/nbseprjprep.f90
	$(COMP) $(OPTS) -c zf90/nbseprjprep.f90 -o zobj/nbseprjprep.o
zobj/seanitup.o: zf90/seanitup.f90
	$(COMP) $(OPTS) -c zf90/seanitup.f90 -o zobj/seanitup.o
zobj/setlm.o: zf90/setlm.f90
	$(COMP) $(OPTS) -c zf90/setlm.f90 -o zobj/setlm.o
zobj/snatch.o: zf90/snatch.f90
	$(COMP) $(OPTS) -c zf90/snatch.f90 -o zobj/snatch.o

install:
	cp zexe/cks.x $(INSTDIR)/

clean:
	rm -f zobj/*.o
	rm -f zobj/*.x
	rm -f *__genmod.???
