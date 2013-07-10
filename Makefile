FC = mpif90
FLAGS = -O0 -ffree-line-length-none -x f95-cpp-input -m64 -DHAVE_FLUSH -DBLAS  -Wall -fcheck=all
LDFLAGS = -m64  -Wall -fcheck=all

FFTWI = -I/home/jtv1/bin/gnu/include
FFTWL = /home/jtv1/bin/gnu/lib/libfftw3.a


BLACS =  -L/home/jtv1/ATLAS/ATLAS_sandy/lib  -lscalapack -llapack -lf77blas -lcblas -latlas


all: ocean.x

OCEANOBJS = OCEAN.o AI_kinds.o OCEAN_mpi.o OCEAN_load_data.o OCEAN_psi.o OCEAN_energies.o OCEAN_haydock.o

ocean.x: $(OCEANOBJS)
	$(FC) $(LDFLAGS) -o ocean.x $(OCEANOBJS) $(BLACS)

OCEAN_mpi.o: OCEAN_mpi.f90 AI_kinds.o
	$(FC) $(FLAGS) -c -o OCEAN_mpi.o OCEAN_mpi.f90 

OCEAN.o: OCEAN.f90 AI_kinds.o
	$(FC) $(FLAGS) -c -o OCEAN.o OCEAN.f90

AI_kinds.o: AI_kinds.f90
	$(FC) $(FLAGS) -c -o AI_kinds.o AI_kinds.f90

OCEAN_load_data.o: OCEAN_load_data.f90
	$(FC) $(FLAGS) -c -o OCEAN_load_data.o OCEAN_load_data.f90

OCEAN_psi.o: OCEAN_psi.f90
	$(FC) $(FLAGS) -c -o OCEAN_psi.o OCEAN_psi.f90

OCEAN_energies.o: OCEAN_energies.f90
	$(FC) $(FLAGS) -c -o OCEAN_energies.o OCEAN_energies.f90

OCEAN_haydock.o: OCEAN_haydock.f90
	$(FC) $(FLAGS) -c -o  OCEAN_haydock.o OCEAN_haydock.f90


clean:
	rm *.o *.mod
