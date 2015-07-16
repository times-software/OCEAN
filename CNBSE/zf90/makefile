COMP=pgf90
OPTS=-fast

all: zexe/getnval.x

zexe/getnval.x: zobj/getnval.o zobj/getomega.o
	$(COMP) $(OPTS) zobj/getnval.o zobj/getomega.o -o zexe/getnval.x

zobj/getnval.o: getnval.f90
	$(COMP) $(OPTS) -c getnval.f90 -o zobj/getnval.o

zobj/getomega.o: getomega.f90
	$(COMP) $(OPTS) -c getomega.f90 -o zobj/getomega.o


install:
	cp zexe/getnval.x $(INSTDIR)/

clean:
	rm -f zobj/*.o
	rm -f zexe/*.x
