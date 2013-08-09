FC = mpif90
FLAGS = -O3 -ffree-line-length-none -x f95-cpp-input -m64 -DHAVE_FLUSH -DBLAS  -Wall  -g -ffree-line-length-none  -funroll-loops -ftree-vectorizer-verbose=0 -ffast-math -ftree-vectorize -DMPI -fopenmp  -march=corei7 
LDFLAGS = -m64  -Wall  -g -fopenmp -march=corei7

FFTWI = -I/home/jtv1/bin/gnu/include
FFTWL = /home/jtv1/bin/gnu/lib/libfftw3.a


BLACS =  -L/home/jtv1/ATLAS/ATLAS_sandy/lib  -lscalapack -llapack -lf77blas -lcblas -latlas


all: ocean.x

OCEANOBJS = AI_kinds.o OCEAN_mpi.o OCEAN_system.o OCEAN_bloch.o OCEAN_obf.o OCEAN_multiplet.o long_range.o OCEAN_load_data.o OCEAN_psi.o \
            OCEAN_energies.o OCEAN_haydock.o OCEAN.o getabb.o getomega.o gamfcn.o jlmfft.o limel.o jimel.o \
            nbsemkcmel.o intval.o newgetylm.o  newgetprefs.o newthreey.o cainmhsetup.o redtrid.o elsdch.o 

ocean.x: $(OCEANOBJS)
	$(FC) $(LDFLAGS) -o ocean.x $(OCEANOBJS) $(BLACS) $(FFTWL)

OCEAN_mpi.o: OCEAN_mpi.f90 AI_kinds.o
	$(FC) $(FLAGS) -c -o OCEAN_mpi.o OCEAN_mpi.f90 

OCEAN.o: OCEAN.f90 AI_kinds.o
	$(FC) $(FLAGS) -c -o OCEAN.o OCEAN.f90

AI_kinds.o: AI_kinds.f90
	$(FC) $(FLAGS) -c -o AI_kinds.o AI_kinds.f90

OCEAN_system.o: OCEAN_system.f90
	$(FC) $(FLAGS) -c -o  OCEAN_system.o OCEAN_system.f90

OCEAN_load_data.o: OCEAN_load_data.f90 OCEAN_energies.o
	$(FC) $(FLAGS) -c -o OCEAN_load_data.o OCEAN_load_data.f90

OCEAN_psi.o: OCEAN_psi.f90
	$(FC) $(FLAGS) -c -o OCEAN_psi.o OCEAN_psi.f90 $(FFTWI)

OCEAN_energies.o: OCEAN_energies.f90
	$(FC) $(FLAGS) -c -o OCEAN_energies.o OCEAN_energies.f90

OCEAN_haydock.o: OCEAN_haydock.f90 
	$(FC) $(FLAGS) -c -o  OCEAN_haydock.o OCEAN_haydock.f90 $(FFTWI)


getabb.o: getabb.f90
	$(FC) $(FLAGS) -c -o  getabb.o getabb.f90

getomega.o: getomega.f90
	$(FC) $(FLAGS) -c -o getomega.o getomega.f90

gamfcn.o: gamfcn.f90
	$(FC) $(FLAGS) -c -o gamfcn.o gamfcn.f90

long_range.o: long_range.f90
	$(FC) $(FLAGS) -c -o long_range.o long_range.f90 $(FFTWI)

jlmfft.o: jlmfft.f
	 $(FC) $(FLAGS) -c -o jlmfft.o jlmfft.f

OCEAN_multiplet.o: OCEAN_multiplet.f90 OCEAN_psi.o
	$(FC) $(FLAGS) -c -o OCEAN_multiplet.o OCEAN_multiplet.f90


limel.o: limel.f90
	$(FC) $(FLAGS) -c -o limel.o limel.f90

jimel.o: jimel.f90
	$(FC) $(FLAGS) -c -o jimel.o jimel.f90

nbsemkcmel.o: nbsemkcmel.f90
	$(FC) $(FLAGS) -c -o nbsemkcmel.o nbsemkcmel.f90

intval.o: intval.f90
	$(FC) $(FLAGS) -c -o intval.o intval.f90

newgetylm.o: newgetylm.f90
	$(FC) $(FLAGS) -c -o newgetylm.o newgetylm.f90

newgetprefs.o: newgetprefs.f90
	$(FC) $(FLAGS) -c -o newgetprefs.o newgetprefs.f90

newthreey.o: newthreey.f90
	$(FC) $(FLAGS) -c -o newthreey.o newthreey.f90

cainmhsetup.o: cainmhsetup.f90
	$(FC) $(FLAGS) -c -o cainmhsetup.o cainmhsetup.f90

redtrid.o: redtrid.f
	$(FC) $(FLAGS) -c -o redtrid.o redtrid.f

elsdch.o: elsdch.f
	$(FC) $(FLAGS) -c -o elsdch.o elsdch.f

OCEAN_bloch.o: OCEAN_bloch.f90 OCEAN_system.o AI_kinds.o OCEAN_mpi.o
	$(FC) $(FLAGS) -c -o OCEAN_bloch.o OCEAN_bloch.f90 $(FFTWI)

OCEAN_obf.o: OCEAN_obf.f90 OCEAN_system.o AI_kinds.o OCEAN_mpi.o
	$(FC) $(FLAGS) -c -o OCEAN_obf.o OCEAN_obf.f90 $(FFTWI)

clean:
	rm *.o *.mod
