include ../Makefile.arch

FOPTS=$(OPTIONS)

all: exciton_plot.x val_exciton_plot.x makewyck.x

val_exciton_plot.x: val_exciton_plot.o periodic.o
	$(FC) $(FOPTS)  -o val_exciton_plot.x val_exciton_plot.o periodic.o $(BLAS)

exciton_plot.x: exciton_plot.o periodic.o
	$(FC)  $(FOPTS) -o exciton_plot.x exciton_plot.o periodic.o $(BLAS)

val_exciton_plot.o: periodic.o periodic.mod OCEAN_val_exciton_plot.f90
	$(FC)  $(FOPTS) -c -o val_exciton_plot.o OCEAN_val_exciton_plot.f90

exciton_plot.o: periodic.o periodic.mod OCEAN_exciton_plot.f90
	$(FC)  $(FOPTS) -c -o exciton_plot.o OCEAN_exciton_plot.f90

makewyck.x: makewyck.f90 periodic.o
	$(FC) $(FOPTS) makewyck.f90 periodic.o -o makewyck.x

periodic.mod: periodic.o

periodic.o: periodic.f90
	$(FC)  $(FOPTS) -c -o periodic.o periodic.f90

install:
	cp exciton_plot.x val_exciton_plot.x makewyck.x $(INSTDIR)/
	cp core-exciton.pl val-exciton.pl $(INSTDIR)/

clean:
	rm *.o *.x
