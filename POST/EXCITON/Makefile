include ../../Makefile.arch

FOPTS=$(OPTIONS) -I../../Modules

all: exciton_plot.x val_exciton_plot.x exciton_band.x

exciton_band.x: OCEAN_exciton_band.f90 poly_interp_mod.o
	$(FC) $(FOPTS) -o exciton_band.x OCEAN_exciton_band.f90 poly_interp_mod.o $(BLAS)

val_exciton_plot.x: val_exciton_plot.o periodic.o
	$(FC) $(FOPTS)  -o val_exciton_plot.x val_exciton_plot.o periodic.o $(BLAS)

exciton_plot.x: exciton_plot.o periodic.o
	$(FC)  $(FOPTS) -o exciton_plot.x exciton_plot.o periodic.o $(BLAS) ../../Modules/ocean_interpolate.o

val_exciton_plot.o: periodic.o periodic.mod OCEAN_val_exciton_plot.f90
	$(FC)  $(FOPTS) -c -o val_exciton_plot.o OCEAN_val_exciton_plot.f90

exciton_plot.o: periodic.o periodic.mod OCEAN_exciton_plot.f90
	$(FC)  $(FOPTS) -c -o exciton_plot.o OCEAN_exciton_plot.f90

periodic.mod: periodic.o

periodic.o: periodic.f90
	$(FC)  $(FOPTS) -c -o periodic.o periodic.f90

poly_interp_mod.o: poly_interp_mod.f90
	$(FC)  $(FOPTS) -c -o poly_interp_mod.o poly_interp_mod.f90
poly_interp_mod.mod: poly_interp_mod.o

install:
	cp exciton_plot.x val_exciton_plot.x exciton_band.x $(INSTDIR)/

clean:
	rm *.o *.x *.mod
