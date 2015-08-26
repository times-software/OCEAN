FC = mpif90
#FLAGS = -O3 -ffree-line-length-none -x f95-cpp-input -m64 -DHAVE_FLUSH -DBLAS  -Wall  -g -ffree-line-length-none  -funroll-loops -ftree-vectorizer-verbose=0 -ffast-math -ftree-vectorize -DMPI -march=corei7 -DHAVE_CONTIGUOUS -fbacktrace  #-fbounds-check -fopenmp
#-fbacktrace -fbounds-check
#LDFLAGS = -m64  -Wall  -g -march=corei7 -fopenmp
FLAGS = -O2 -DHAVE_FLUSH -DBLAS -DMPI -DHAVE_CONTIGUOUS -cpp -march=corei7  -traceback  -warn all #-fopenmp
#-align array64byte
LDFLAGS =  -static-intel -traceback #-fopenmp

#FFTWI = -I/home/jtv1/bin/gnu/include
#FFTWL = /home/jtv1/bin/gnu/lib/libfftw3.a
##FFTWI = -I/users/jtv1/local/fftw-3.3.3/api
##FFTWL =  -L/users/jtv1/local/.libs -lfftw3

##BLACS = -lscalapack -llapack -lf77blas -lcblas -latlas
#BLACS =  -L/home/jtv1/ATLAS/ATLAS_sandy/lib  -lscalapack -llapack -lf77blas -lcblas -latlas
BLACS = -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64
FFTWL = /home/jtv1/bin/intel/lib/libfftw3.a
FFTWI = -I/home/jtv1/bin/intel/include


all: ocean2.x

OCEANOBJS = AI_kinds.o OCEAN_mpi.o OCEAN_system.o OCEAN_bloch.o OCEAN_obf.o OCEAN_multiplet.o  \
            long_range.o OCEAN_load_data.o OCEAN_psi.o OCEAN_energies.o OCEAN_haydock.o OCEAN.o \
            getabb.o getomega.o gamfcn.o jlmfft.o jimel.o nbsemkcmel.o intval.o \
            newgetylm.o  newgetprefs.o newthreey.o cainmhsetup.o elsdch.o cinv.o \
            sizereport.o OCEAN_timekeeper.o OCEAN_invdrv.o \
            OCEAN_val_states.o OCEAN_ladder.o OCEAN_val_energy.o \
            OCEAN_read_tmels.o OCEAN_hyb_louie_levine.o OCEAN_get_rho.o kmapr.o optim.o \
						OCEAN_rixs_holder.o

ocean2.x: $(OCEANOBJS)
	$(FC) $(LDFLAGS) -o ocean2.x $(OCEANOBJS) $(FFTWL) $(BLACS) 

OCEAN_mpi.o: OCEAN_mpi.f90 AI_kinds.o
	$(FC) $(FLAGS) -c -o OCEAN_mpi.o OCEAN_mpi.f90 

OCEAN.o: OCEAN.f90 AI_kinds.o 
	$(FC) $(FLAGS) -c -o OCEAN.o OCEAN.f90

AI_kinds.o: AI_kinds.f90
	$(FC) $(FLAGS) -c -o AI_kinds.o AI_kinds.f90

OCEAN_timekeeper.o: OCEAN_timekeeper.f90
	$(FC) $(FLAGS) -c -o OCEAN_timekeeper.o OCEAN_timekeeper.f90

OCEAN_system.o: OCEAN_system.f90
	$(FC) $(FLAGS) -c -o  OCEAN_system.o OCEAN_system.f90

OCEAN_load_data.o: OCEAN_load_data.f90 OCEAN_energies.o
	$(FC) $(FLAGS) -c -o OCEAN_load_data.o OCEAN_load_data.f90 

OCEAN_psi.o: OCEAN_psi.f90 OCEAN_rixs_holder.o
	$(FC) $(FLAGS) -c -o OCEAN_psi.o OCEAN_psi.f90 $(FFTWI)

OCEAN_energies.o: OCEAN_energies.f90 OCEAN_val_energy.o
	$(FC) $(FLAGS) -c -o OCEAN_energies.o OCEAN_energies.f90 $(FFTWI)

OCEAN_haydock.o: OCEAN_haydock.f90 OCEAN_timekeeper.o
	$(FC) $(FLAGS) -c -o  OCEAN_haydock.o OCEAN_haydock.f90 $(FFTWI)

OCEAN_exact.o: OCEAN_exact.f90
	$(FC) $(FLAGS) -c -o  OCEAN_exact.o OCEAN_exact.f90

getabb.o: getabb.f90
	$(FC) $(FLAGS) -c -o  getabb.o getabb.f90

getomega.o: getomega.f90
	$(FC) $(FLAGS) -c -o getomega.o getomega.f90

gamfcn.o: gamfcn.f90
	$(FC) $(FLAGS) -c -o gamfcn.o gamfcn.f90

long_range.o: long_range.f90
	$(FC) $(FLAGS) -c -o long_range.o long_range.f90 $(FFTWI)

jlmfft.o: jlmfft.f
	 $(FC) -c -o jlmfft.o jlmfft.f

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

invdrv.o: invdrv.f90
	$(FC) $(FLAGS) -c -o invdrv.o invdrv.f90

cinv.o: cinv.f
	$(FC) $(FLAGS) -c -o cinv.o cinv.f

sizereport.o: sizereport.f90
	$(FC) $(FLAGS) -c -o sizereport.o sizereport.f90

OCEAN_bloch.o: OCEAN_bloch.f90 OCEAN_system.o AI_kinds.o OCEAN_mpi.o
	$(FC) $(FLAGS) -c -o OCEAN_bloch.o OCEAN_bloch.f90 $(FFTWI)

OCEAN_obf.o: OCEAN_obf.f90 OCEAN_system.o AI_kinds.o OCEAN_mpi.o
	$(FC) $(FLAGS) -c -o OCEAN_obf.o OCEAN_obf.f90 $(FFTWI)

OCEAN_invdrv.o: OCEAN_invdrv.f90
	$(FC) $(FLAGS) -c -o OCEAN_invdrv.o OCEAN_invdrv.f90

OCEAN_val_states.o: OCEAN_val_states.f90
	$(FC) $(FLAGS) -c -o OCEAN_val_states.o OCEAN_val_states.f90

OCEAN_ladder.o : OCEAN_ladder.f90 OCEAN_hyb_louie_levine.o
	$(FC) $(FLAGS) -c -o OCEAN_ladder.o OCEAN_ladder.f90

OCEAN_bubble.o : OCEAN_bubble.f90
	$(FC) $(FLAGS) -c -o OCEAN_bubble.o OCEAN_bubble.f90

OCEAN_val_energy.o: OCEAN_val_energy.f90
	$(FC) $(FLAGS) -c -o OCEAN_val_energy.o OCEAN_val_energy.f90

OCEAN_read_tmels.o: OCEAN_read_tmels.f90
	$(FC) $(FLAGS) -c -o OCEAN_read_tmels.o OCEAN_read_tmels.f90

OCEAN_hyb_louie_levine.o: OCEAN_hyb_louie_levine.f90
	$(FC) $(FLAGS) -c -o OCEAN_hyb_louie_levine.o OCEAN_hyb_louie_levine.f90

OCEAN_get_rho.o: OCEAN_get_rho.f90
	$(FC) $(FLAGS) -c -o OCEAN_get_rho.o OCEAN_get_rho.f90

kmapr.o: kmapr.f90
	$(FC) $(FLAGS) -c -o kmapr.o kmapr.f90

optim.o: optim.f90
	$(FC) $(FLAGS) -c -o optim.o optim.f90

OCEAN_rixs_holder.o: OCEAN_rixs_holder.f90
	$(FC) $(FLAGS) -c -o OCEAN_rixs_holder.o OCEAN_rixs_holder.f90

clean:
	rm *.o *.mod
