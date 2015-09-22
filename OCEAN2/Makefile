include ../Makefile.arch

FLAGS = $(OPTIONS) -cpp -DHAVE_FLUSH -DBLAS -DMPI -DHAVE_CONTIGUOUS 

#MPIFORT = mpif90
#FLAGS = -O2 -DHAVE_FLUSH -DBLAS -DMPI -DHAVE_CONTIGUOUS -cpp -warn all  -march=corei7 -fopenmp -openmp-report2 -g -traceback
#LDFLAGS = -fopenmp -traceback

#FFTWL = /home/jtv1/bin/intel/lib/libfftw3.a
#FFTWI = -I/home/jtv1/bin/intel/include

# Need to specify scalapack, blacs, blas
BLACS = -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64

CNBSE = ../CNBSE/zcmsnmult/zf90


all: ocean.x

OCEANOBJS = AI_kinds.o OCEAN_mpi.o OCEAN_system.o OCEAN_bloch.o OCEAN_obf.o OCEAN_multiplet.o  \
            long_range.o OCEAN_load_data.o OCEAN_psi.o OCEAN_energies.o OCEAN_haydock.o OCEAN.o \
            getabb.o getomega.o gamfcn.o jlmfft.o jimel.o nbsemkcmel.o intval.o \
            newgetylm.o  newgetprefs.o newthreey.o cainmhsetup.o elsdch.o cinv.o \
            sizereport.o OCEAN_timekeeper.o OCEAN_invdrv.o

ocean.x: $(OCEANOBJS)
	$(MPIFORT) $(LDFLAGS) -o ocean.x $(OCEANOBJS) $(FFTWL) $(BLACS) 

OCEAN_mpi.o: OCEAN_mpi.f90 AI_kinds.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_mpi.o OCEAN_mpi.f90 

OCEAN.o: OCEAN.f90 AI_kinds.o 
	$(MPIFORT) $(FLAGS) -c -o OCEAN.o OCEAN.f90

AI_kinds.o: AI_kinds.f90
	$(MPIFORT) $(FLAGS) -c -o AI_kinds.o AI_kinds.f90

OCEAN_timekeeper.o: OCEAN_timekeeper.f90
	$(MPIFORT) $(FLAGS) -c -o OCEAN_timekeeper.o OCEAN_timekeeper.f90

OCEAN_system.o: OCEAN_system.f90
	$(MPIFORT) $(FLAGS) -c -o  OCEAN_system.o OCEAN_system.f90

OCEAN_load_data.o: OCEAN_load_data.f90 OCEAN_energies.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_load_data.o OCEAN_load_data.f90 

OCEAN_psi.o: OCEAN_psi.f90
	$(MPIFORT) $(FLAGS) -c -o OCEAN_psi.o OCEAN_psi.f90 $(FFTWI)

OCEAN_energies.o: OCEAN_energies.f90
	$(MPIFORT) $(FLAGS) -c -o OCEAN_energies.o OCEAN_energies.f90 $(FFTWI)

OCEAN_haydock.o: OCEAN_haydock.f90 OCEAN_timekeeper.o
	$(MPIFORT) $(FLAGS) -c -o  OCEAN_haydock.o OCEAN_haydock.f90 $(FFTWI)

OCEAN_exact.o: OCEAN_exact.f90
	$(MPIFORT) $(FLAGS) -c -o  OCEAN_exact.o OCEAN_exact.f90

getabb.o: $(CNBSE)/getabb.f90
	$(MPIFORT) $(FLAGS) -c -o  getabb.o $(CNBSE)/getabb.f90

getomega.o: $(CNBSE)/getomega.f90
	$(MPIFORT) $(FLAGS) -c -o getomega.o $(CNBSE)/getomega.f90

gamfcn.o: $(CNBSE)/gamfcn.f90
	$(MPIFORT) $(FLAGS) -c -o gamfcn.o $(CNBSE)/gamfcn.f90

long_range.o: long_range.f90
	$(MPIFORT) $(FLAGS) -c -o long_range.o long_range.f90 $(FFTWI)

jlmfft.o: $(CNBSE)/../zstd/jlmfft.f
	 $(MPIFORT) -c -o jlmfft.o $(CNBSE)/../zstd/jlmfft.f

OCEAN_multiplet.o: OCEAN_multiplet.f90 OCEAN_psi.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_multiplet.o OCEAN_multiplet.f90


limel.o: $(CNBSE)/limel.f90
	$(MPIFORT) $(FLAGS) -c -o limel.o $(CNBSE)/limel.f90

jimel.o: $(CNBSE)/jimel.f90
	$(MPIFORT) $(FLAGS) -c -o jimel.o $(CNBSE)/jimel.f90

nbsemkcmel.o: nbsemkcmel.f90
	$(MPIFORT) $(FLAGS) -c -o nbsemkcmel.o nbsemkcmel.f90

intval.o: $(CNBSE)/intval.f90
	$(MPIFORT) $(FLAGS) -c -o intval.o $(CNBSE)/intval.f90

newgetylm.o: $(CNBSE)/newgetylm.f90
	$(MPIFORT) $(FLAGS) -c -o newgetylm.o $(CNBSE)/newgetylm.f90

newgetprefs.o: $(CNBSE)/newgetprefs.f90
	$(MPIFORT) $(FLAGS) -c -o newgetprefs.o $(CNBSE)/newgetprefs.f90

newthreey.o: $(CNBSE)/newthreey.f90
	$(MPIFORT) $(FLAGS) -c -o newthreey.o $(CNBSE)/newthreey.f90

cainmhsetup.o: $(CNBSE)/cainmhsetup.f90
	$(MPIFORT) $(FLAGS) -c -o cainmhsetup.o $(CNBSE)/cainmhsetup.f90

redtrid.o: redtrid.f
	$(MPIFORT) $(FLAGS) -c -o redtrid.o redtrid.f

elsdch.o: $(CNBSE)/../zstd/elsdch.f
	$(MPIFORT) $(FLAGS) -c -o elsdch.o $(CNBSE)/../zstd/elsdch.f

invdrv.o: $(CNBSE)/invdrv.f90
	$(MPIFORT) $(FLAGS) -c -o invdrv.o $(CNBSE)/invdrv.f90

cinv.o: $(CNBSE)/../zstd/cinv.f
	$(MPIFORT) $(FLAGS) -c -o cinv.o $(CNBSE)/../zstd/cinv.f

sizereport.o: $(CNBSE)/sizereport.f90
	$(MPIFORT) $(FLAGS) -c -o sizereport.o $(CNBSE)/sizereport.f90

OCEAN_bloch.o: OCEAN_bloch.f90 OCEAN_system.o AI_kinds.o OCEAN_mpi.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_bloch.o OCEAN_bloch.f90 $(FFTWI)

OCEAN_obf.o: OCEAN_obf.f90 OCEAN_system.o AI_kinds.o OCEAN_mpi.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_obf.o OCEAN_obf.f90 $(FFTWI)

OCEAN_invdrv.o: OCEAN_invdrv.f90
	$(MPIFORT) $(FLAGS) -c -o OCEAN_invdrv.o OCEAN_invdrv.f90


clean:
	rm *.o *.mod

install:
	cp ocean.x $(INSTDIR)/
