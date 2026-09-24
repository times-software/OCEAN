include ../../Makefile.arch

#OPTIONS += -prof-gen
#OPTIONS += -prof-use 
#OPTIONS += -march=haswell -prec-div -prec-sqrt
#OPTIONS += -check all -g -traceback
#OPTIONS += -qopenmp
OPTIONS += -fcheck=all -g -fbacktrace

OBJS = newthreey.o newgetprefs.o newgetylm.o getomega.o

all: corex.x corescreen.x

coreScreen.o: coreScreen.f90
	$(FC) $(OPTIONS) -c -o coreScreen.o coreScreen.f90 -I../../Modules

corescreen.x: coreScreen.o
	$(FC) $(OPTIONS) $(LDFLAGS) -o corescreen.x coreScreen.o -I../../Modules

coreExchange.o: coreExchange.f90 
	$(FC) $(OPTIONS) -c -o coreExchange.o coreExchange.f90 -I../../Modules 
	
corex.x: coreExchange.o $(OBJS)
	$(FC) $(OPTIONS) $(LDFLAGS) -o corex.x coreExchange.o -I../../Modules $(OBJS)

newthreey.o:
	$(FC) $(OPTIONS) -c -o newthreey.o newthreey.f90

newgetprefs.o: 
	$(FC) $(OPTIONS) -c -o newgetprefs.o newgetprefs.f90

newgetylm.o:
	$(FC) $(OPTIONS) -c -o newgetylm.o newgetylm.f90

getomega.o:
	$(FC) $(OPTIONS) -c -o getomega.o getomega.f90

clean:
	rm coreExchange.o newthreey.o newgetprefs.o newgetylm.o

install:
	cp corex.x $(INSTDIR)/
