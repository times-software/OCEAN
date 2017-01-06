all: zexe/fhi2eric.x

zexe/fhi2eric.x: zf90/fhi2eric.f90
	$(COMP) $(OPTS) zf90/fhi2eric.f90 -o zexe/fhi2eric.x 

install:
	cp zexe/fhi2eric.x $(INSTDIR)/

clean:
	rm -f zexe/*.x
