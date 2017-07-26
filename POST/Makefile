include ../Makefile.arch

all: exciton_plot.x

exciton_plot.x: exciton_plot.o periodic.o
	$(FC) -o exciton_plot.x exciton_plot.o periodic.o $(BLAS)

exciton_plot.o: periodic.o periodic.mod OCEAN_exciton_plot.f90
	$(FC) -c -o exciton_plot.o OCEAN_exciton_plot.f90

periodic.mod: periodic.o

periodic.o: periodic.f90
	$(FC) -c -o periodic.o periodic.f90
