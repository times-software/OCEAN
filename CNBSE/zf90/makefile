COMP=pgf90
OPTS=-fast

all: zexe/getnval.x zexe/spec_average.x

zexe/getnval.x: zobj/getnval.o zobj/getomega.o
	$(COMP) $(OPTS) zobj/getnval.o zobj/getomega.o -o zexe/getnval.x

zexe/spec_average.x: zobj/spec_average.o
	$(COMP) $(OPTS) zobj/spec_average.o -o zexe/spec_average.x

zobj/getnval.o: getnval.f90
	$(COMP) $(OPTS) -c getnval.f90 -o zobj/getnval.o

zobj/getomega.o: getomega.f90
	$(COMP) $(OPTS) -c getomega.f90 -o zobj/getomega.o

zobj/spec_average.o: spec_average.f90
	$(COMP) $(OPTS) -c spec_average.f90 -o zobj/spec_average.o

install:
	cp zexe/getnval.x $(INSTDIR)/
	cp zexe/spec_average.x $(INSTDIR)/

clean:
	rm -f zobj/*.o
	rm -f zexe/*.x
