OBJGEN = zobj/limel.o zobj/jimel.o zobj/getylm.o zobj/getprefs.o \
         zobj/nbseprjprep.o zobj/nbseprjnproj.o zobj/fancyvector.o
OBJELS = $(OBJGEN) zobj/elssub.o zobj/confmel.o zobj/sphidriver.o
OBJJTV = $(OBJGEN) zobj/jtvsub.o zobj/cmjtv.o zobj/sdjtv.o zobj/jlmatfetch.o zobj/ckmels.o
#
all: zexe/melcalc.x zexe/dotter.x zexe/meljtv.x
#
zexe/melcalc.x: $(OBJELS)
	$(COMP) $(OPTS) -o zexe/melcalc.x $(OBJELS)
zexe/meljtv.x: $(OBJJTV)
	$(COMP) $(OPTS) -o zexe/meljtv.x $(OBJJTV)
zexe/dotter.x: zf90/dotter.f90
	$(COMP) -o zexe/dotter.x zf90/dotter.f90
#
zobj/confmel.o: zf90/confmel.f90
	$(COMP) -c $(OPTS) -o zobj/confmel.o zf90/confmel.f90 -I../
zobj/cmjtv.o: zf90/cmjtv.f90
	$(COMP) -c $(OPTS) -o zobj/cmjtv.o zf90/cmjtv.f90 -I../
zobj/getprefs.o: zf90/getprefs.f90
	$(COMP) -c $(OPTS) -o zobj/getprefs.o zf90/getprefs.f90
zobj/getylm.o: zf90/getylm.f90
	$(COMP) -c $(OPTS) -o zobj/getylm.o zf90/getylm.f90
zobj/jimel.o: zf90/jimel.f90
	$(COMP) -c $(OPTS) -o zobj/jimel.o zf90/jimel.f90
zobj/elssub.o: zf90/elssub.f90
	$(COMP) -c $(OPTS) -o zobj/elssub.o zf90/elssub.f90
zobj/jlmatfetch.o: zf90/jlmatfetch.f90
	$(COMP) -c $(OPTS) -o zobj/jlmatfetch.o zf90/jlmatfetch.f90
zobj/ckmels.o: zf90/ckmels.f90
	$(COMP) -c $(OPTS) -o zobj/ckmels.o zf90/ckmels.f90
zobj/jtvsub.o: zf90/jtvsub.f90
	$(COMP) -c $(OPTS) -o zobj/jtvsub.o zf90/jtvsub.f90
zobj/limel.o: zf90/limel.f90
	$(COMP) -c $(OPTS) -o zobj/limel.o zf90/limel.f90
zobj/sdjtv.o: zf90/sdjtv.f90
	$(COMP) -c $(OPTS) -o zobj/sdjtv.o zf90/sdjtv.f90
zobj/sphidriver.o: zf90/sphidriver.f90 zf90/sphsetnx.h.f90 zf90/sphsetx.h.f90
	$(COMP) -c $(OPTS) -o zobj/sphidriver.o zf90/sphidriver.f90
#
zobj/nbseprjprep.o: zf90/nbseprjprep.f90
	$(COMP) -c $(OPTS) -o zobj/nbseprjprep.o zf90/nbseprjprep.f90
zobj/nbseprjnproj.o: zf90/nbseprjnproj.f90
	$(COMP) -c $(OPTS) -o zobj/nbseprjnproj.o zf90/nbseprjnproj.f90
#
zobj/fancyvector.o: zf90/fancyvector.f90
	$(COMP) -c $(OPTS) -o zobj/fancyvector.o zf90/fancyvector.f90


install:
	cp zexe/melcalc.x $(INSTDIR)/
	cp zexe/dotter.x $(INSTDIR)/
	cp zexe/meljtv.x $(INSTDIR)/

clean:
	rm -f zobj/*.o
	rm -f zexe/*.x
	rm -f *__genmod.???
	rm -f *.mod
