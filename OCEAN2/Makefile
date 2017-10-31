include ../Makefile.arch

# Currently BLAS and MPI are *not* optional
FLAGS = $(OPTIONS) #-Wall -finit-real=nan

CNBSE = ../CNBSE/zcmsnmult/zf90


all: ocean.x o_spect.x rixs_plane.x

OCEANOBJS = AI_kinds.o OCEAN_mpi.o OCEAN_system.o OCEAN_bloch.o OCEAN_obf.o OCEAN_multiplet.o  \
            long_range.o OCEAN_load_data.o OCEAN_psi.o OCEAN_energies.o OCEAN_haydock.o OCEAN.o \
            getabb.o getomega.o gamfcn.o jimel.o nbsemkcmel.o intval.o \
            newgetylm.o  newgetprefs.o newthreey.o cainmhsetup.o elsdch.o cinv.o \
            sizereport.o OCEAN_timekeeper.o OCEAN_invdrv.o OCEAN_val_energy.o \
            OCEAN_read_tmels.o OCEAN_val_states.o OCEAN_bubble.o \
            OCEAN_get_rho.o optim.o OCEAN_hyb_louie_levine.o OCEAN_ladder.o \
						OCEAN_pfy.o OCEAN_rixs_holder.o OCEAN_corewidths.o OCEAN_action.o \
						OCEAN_driver.o OCEAN_gmres.o OCEAN_filenames.o $(FFT_WRAPPER_OBJS) $(CAIN_OBJS)

KINDS = AI_kinds.o

FFT_WRAPPER_OBJS = ../zbridge/zobj/FFT_wrapper.o ../zbridge/zobj/jlmfft.o
FFT_WRAPPER_INCL = -I../zbridge/ 

CAIN_OBJS = ../CNBSE/zcain/zobj/getprefs.o ../CNBSE/zcain/zobj/getylm.o ../CNBSE/zcain/zobj/fancyvector.o

OCEAN_val_energy.o: OCEAN_val_energy.f90 OCEAN_psi.o OCEAN_mpi.o OCEAN_system.o OCEAN_constants.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_val_energy.o OCEAN_val_energy.f90

OCEAN_constants.o: OCEAN_constants.f90 AI_kinds.o 
	$(MPIFORT) $(FLAGS) -c -o OCEAN_constants.o OCEAN_constants.f90

ocean.x: $(OCEANOBJS)
	$(MPIFORT) $(LDFLAGS) -o ocean.x $(OCEANOBJS) $(FFTWL)  $(BLAS) $(OPENMP) $(FFTWL)

o_spect.x: o_spect.f90 AI_kinds.o gamfcn.o
	$(MPIFORT) $(FLAGS) -o o_spect.x o_spect.f90 AI_kinds.o gamfcn.o

rixs_plane.x: rixs_plane.f90 AI_kinds.o
	$(FC) $(FLAGS) -o rixs_plane.x rixs_plane.f90 AI_kinds.o

OCEAN_mpi.o: OCEAN_mpi.f90 AI_kinds.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_mpi.o OCEAN_mpi.f90 

OCEAN.o: OCEAN.f90 AI_kinds.o OCEAN_driver.o OCEAN_system.o OCEAN_timekeeper.o OCEAN_mpi.o OCEAN_psi.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN.o OCEAN.f90

AI_kinds.o: AI_kinds.f90
	$(MPIFORT) $(FLAGS) -c -o AI_kinds.o AI_kinds.f90

OCEAN_timekeeper.o: OCEAN_timekeeper.f90
	$(MPIFORT) $(FLAGS) -c -o OCEAN_timekeeper.o OCEAN_timekeeper.f90

OCEAN_system.o: OCEAN_system.f90 AI_kinds.o OCEAN_mpi.o
	$(MPIFORT) $(FLAGS) -c -o  OCEAN_system.o OCEAN_system.f90

OCEAN_load_data.o: OCEAN_load_data.f90 OCEAN_energies.o OCEAN_val_states.o OCEAN_bubble.o OCEAN_ladder.o  \
									 long_range.o OCEAN_multiplet.o OCEAN_mpi.o OCEAN_psi.o OCEAN_system.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_load_data.o OCEAN_load_data.f90 

OCEAN_psi.o: OCEAN_psi.f90 OCEAN_constants.o OCEAN_timekeeper.o OCEAN_rixs_holder.o OCEAN_system.o OCEAN_mpi.o 
	$(MPIFORT) $(FLAGS) -c -o OCEAN_psi.o OCEAN_psi.f90 $(OPENMP)

OCEAN_energies.o: OCEAN_energies.f90 OCEAN_val_energy.o OCEAN_system.o OCEAN_psi.o OCEAN_constants.o 
	$(MPIFORT) $(FLAGS) -c -o OCEAN_energies.o OCEAN_energies.f90 

OCEAN_haydock.o: OCEAN_haydock.f90 OCEAN_timekeeper.o OCEAN_psi.o OCEAN_pfy.o OCEAN_corewidths.o OCEAN_action.o \
								 OCEAN_filenames.o
	$(MPIFORT) $(FLAGS) -c -o  OCEAN_haydock.o OCEAN_haydock.f90 

OCEAN_action.o: OCEAN_action.f90 OCEAN_timekeeper.o OCEAN_psi.o OCEAN_bubble.o OCEAN_ladder.o OCEAN_mpi.o \
								OCEAN_system.o OCEAN_energies.o OCEAN_multiplet.o long_range.o 
	$(MPIFORT) $(FLAGS) -c -o OCEAN_action.o OCEAN_action.f90

OCEAN_exact.o: OCEAN_exact.f90
	$(MPIFORT) $(FLAGS) -c -o  OCEAN_exact.o OCEAN_exact.f90

getabb.o: $(CNBSE)/getabb.f90
	$(MPIFORT) $(FLAGS) -c -o  getabb.o $(CNBSE)/getabb.f90

getomega.o: $(CNBSE)/getomega.f90
	$(MPIFORT) $(FLAGS) -c -o getomega.o $(CNBSE)/getomega.f90

gamfcn.o: $(CNBSE)/gamfcn.f90
	$(MPIFORT) $(FLAGS) -c -o gamfcn.o $(CNBSE)/gamfcn.f90

long_range.o: long_range.f90 OCEAN_system.o OCEAN_mpi.o OCEAN_obf.o
	$(MPIFORT) $(FLAGS) -c -o long_range.o long_range.f90 $(OPENMP) $(FFT_WRAPPER_INCL)

#jlmfft.o: $(CNBSE)/../zstd/jlmfft.f
#	 $(MPIFORT) -c -o jlmfft.o $(CNBSE)/../zstd/jlmfft.f

OCEAN_multiplet.o: OCEAN_multiplet.f90 OCEAN_psi.o OCEAN_system.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_multiplet.o OCEAN_multiplet.f90 $(OPENMP)


limel.o: $(CNBSE)/limel.f90
	$(MPIFORT) $(FLAGS) -c -o limel.o $(CNBSE)/limel.f90

jimel.o: $(CNBSE)/jimel.f90
	$(MPIFORT) $(FLAGS) -c -o jimel.o $(CNBSE)/jimel.f90

nbsemkcmel.o: nbsemkcmel.f90 $(KINDS)
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
	$(MPIFORT) $(FLAGS) -c -o OCEAN_bloch.o OCEAN_bloch.f90 

OCEAN_obf.o: OCEAN_obf.f90 OCEAN_system.o AI_kinds.o OCEAN_mpi.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_obf.o OCEAN_obf.f90 

OCEAN_invdrv.o: OCEAN_invdrv.f90
	$(MPIFORT) $(FLAGS) -c -o OCEAN_invdrv.o OCEAN_invdrv.f90

OCEAN_read_tmels.o: OCEAN_read_tmels.f90 OCEAN_system.o OCEAN_mpi.o OCEAN_psi.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_read_tmels.o OCEAN_read_tmels.f90

OCEAN_bubble.o: OCEAN_bubble.f90 OCEAN_mpi.o OCEAN_val_states.o OCEAN_system.o AI_kinds.o OCEAN_constants.o OCEAN_psi.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_bubble.o OCEAN_bubble.f90 $(FFT_WRAPPER_INCL) $(OPENMP)

OCEAN_val_states.o: OCEAN_val_states.f90 AI_kinds.o OCEAN_mpi.o OCEAN_system.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_val_states.o OCEAN_val_states.f90

OCEAN_get_rho.o: OCEAN_get_rho.f90 OCEAN_mpi.o AI_kinds.o 
	$(MPIFORT) $(FLAGS) -c -o OCEAN_get_rho.o OCEAN_get_rho.f90 $(FFT_WRAPPER_INCL)

optim.o: optim.f90
	$(MPIFORT) $(FLAGS) -c -o optim.o optim.f90

OCEAN_hyb_louie_levine.o: OCEAN_hyb_louie_levine.f90 OCEAN_constants.o AI_kinds.o OCEAN_system.o OCEAN_mpi.o 
	$(MPIFORT) $(FLAGS) -c -o OCEAN_hyb_louie_levine.o OCEAN_hyb_louie_levine.f90

OCEAN_ladder.o: OCEAN_ladder.f90 OCEAN_psi.o OCEAN_mpi.o OCEAN_val_states.o OCEAN_system.o OCEAN_hyb_louie_levine.o 
	$(MPIFORT) $(FLAGS) -c -o OCEAN_ladder.o OCEAN_ladder.f90 $(FFT_WRAPPER_INCL) $(OPENMP)

OCEAN_pfy.o: OCEAN_pfy.f90 OCEAN_psi.o OCEAN_mpi.o OCEAN_system.o OCEAN_constants.o AI_kinds.o 
	$(MPIFORT) $(FLAGS) -c -o OCEAN_pfy.o OCEAN_pfy.f90

OCEAN_rixs_holder.o: OCEAN_rixs_holder.f90 AI_kinds.o OCEAN_system.o OCEAN_mpi.o OCEAN_filenames.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_rixs_holder.o OCEAN_rixs_holder.f90

OCEAN_corewidths.o: OCEAN_corewidths.f90 AI_kinds.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_corewidths.o OCEAN_corewidths.f90

OCEAN_driver.o: OCEAN_driver.f90 OCEAN_mpi.o OCEAN_system.o OCEAN_psi.o OCEAN_haydock.o OCEAN_gmres.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_driver.o OCEAN_driver.f90

OCEAN_gmres.o: OCEAN_gmres.f90 AI_kinds.o OCEAN_psi.o OCEAN_system.o OCEAN_mpi.o \
							 OCEAN_corewidths.o OCEAN_constants.o OCEAN_action.o OCEAN_filenames.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_gmres.o OCEAN_gmres.f90

OCEAN_filenames.o: OCEAN_filenames.f90 OCEAN_system.o
	$(MPIFORT) $(FLAGS) -c -o OCEAN_filenames.o OCEAN_filenames.f90

clean:
	rm -f *.o *.mod ocean.x o_spect.x rixs_plane.x

install:
	cp ocean.x o_spect.x rixs_plane.x $(INSTDIR)/
