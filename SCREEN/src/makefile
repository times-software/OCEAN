include ../../Makefile.arch

OPTIONS += -warn all -check all -DDEBUG


all: screen_grid.o screen_system.o screen_driver.x

OBJS = screen_driver.o ai_kinds.o ocean_mpi.o screen_grid.o screen_system.o screen_sites.o screen_paral.o \
       screen_energy.o ocean_dft_files.o ocean_legacy_files.o screen_wavefunction.o screen_wvfn_converter.o \
       screen_chi0.o screen_chi_driver.o screen_chi.o screen_centralpotential.o periodic.o schi_sinqr.o \
       schi_direct.o ocean_sphericalharmonics.o

%.o %.mod: %.f90
	$(MPIFORT) $(OPTIONS) -c -o $*.o $<

ocean_mpi.o: ocean_mpi.f90 ai_kinds.mod
ocean_mpi.mod: ocean_mpi.f90 ai_kinds.mod

ai_kinds.o: ai_kinds.f90
ai_kinds.mod: ai_kinds.f90

ocean_constants.o: ocean_constants.f90 ai_kinds.mod
ocean_constants.mod: ocean_constants.f90 ai_kinds.mod

screen_grid.o: screen_grid.f90 ai_kinds.mod ocean_constants.mod ocean_mpi.mod 
screen_grid.mod: screen_grid.f90 ai_kinds.mod ocean_constants.mod ocean_mpi.mod 

screen_system.o: screen_system.f90 ai_kinds.mod ocean_mpi.mod
screen_system.mod: screen_system.f90 ai_kinds.mod ocean_mpi.mod

screen_paral.o: screen_paral.f90 ai_kinds.mod ocean_mpi.mod
screen_paral.mod: screen_paral.f90 ai_kinds.mod ocean_mpi.mod

screen_sites.o: screen_sites.f90 screen_paral.mod screen_grid.mod ocean_mpi.mod ai_kinds.mod screen_wavefunction.mod \
                periodic.mod
screen_sites.mod: screen_sites.f90 screen_paral.mod screen_grid.mod ocean_mpi.mod ai_kinds.mod screen_wavefunction.mod \
                  periodic.mod

screen_energy.o: screen_energy.f90 ai_kinds.mod ocean_mpi.mod screen_system.mod ocean_dft_files.mod
screen_energy.mod: screen_energy.f90 ai_kinds.mod ocean_mpi.mod screen_system.mod ocean_dft_files.mod

ocean_prep.o: ocean_prep.f90 ocean_legacy_files.mod ai_kinds.mod ocean_mpi.mod
ocean_prep.mod: ocean_prep.f90 ocean_legacy_files.mod ai_kinds.mod ocean_mpi.mod

ocean_legacy_files.o: ocean_legacy_files.f90 ai_kinds.mod ocean_mpi.mod 
ocean_legacy_files.mod: ocean_legacy_files.f90 ai_kinds.mod ocean_mpi.mod

ocean_dft_files.o: ocean_dft_files.f90 ai_kinds.mod ocean_mpi.mod ocean_legacy_files.mod
ocean_dft_files.mod: ocean_dft_files.f90 ai_kinds.mod ocean_mpi.mod ocean_legacy_files.mod

screen_wavefunction.o:   screen_wavefunction.f90 ai_kinds.mod screen_system.mod screen_paral.mod
screen_wavefunction.mod: screen_wavefunction.f90 ai_kinds.mod screen_system.mod screen_paral.mod

screen_wvfn_converter.o:   screen_wvfn_converter.f90 ai_kinds.mod ocean_dft_files.mod ocean_legacy_files.mod \
                           screen_sites.mod screen_grid.mod screen_paral.mod
screen_wvfn_converter.mod: screen_wvfn_converter.f90 ai_kinds.mod ocean_dft_files.mod ocean_legacy_files.mod \
                           screen_sites.mod screen_grid.mod screen_paral.mod

screen_chi0.o:   screen_chi0.f90 ai_kinds.mod screen_paral.mod screen_sites.mod screen_grid.mod screen_wavefunction.mod
screen_chi0.mod: screen_chi0.f90 ai_kinds.mod screen_paral.mod screen_sites.mod screen_grid.mod screen_wavefunction.mod

screen_chi.o:   screen_chi.f90 ai_kinds.mod screen_grid.mod ocean_mpi.mod ocean_constants.mod schi_sinqr.mod \
                screen_centralpotential.mod schi_direct.mod
screen_chi.mod: screen_chi.f90 ai_kinds.mod screen_grid.mod ocean_mpi.mod ocean_constants.mod schi_sinqr.mod \
                screen_centralpotential.mod schi_direct.mod

screen_chi_driver.o:   screen_chi_driver.f90 screen_chi0.mod screen_chi.mod screen_centralpotential.mod
screen_chi_driver.mod: screen_chi_driver.f90 screen_chi0.mod screen_chi.mod screen_centralpotential.mod

screen_centralpotential.o:   screen_centralpotential.f90 ai_kinds.mod ocean_mpi.mod
screen_centralpotential.mod: screen_centralpotential.f90 ai_kinds.mod ocean_mpi.mod

schi_sinqr.o:   schi_sinqr.f90 ai_kinds.mod screen_centralpotential.mod
schi_sinqr.mod: schi_sinqr.f90 ai_kinds.mod screen_centralpotential.mod

schi_direct.o:   schi_direct.f90 ai_kinds.mod screen_centralpotential.mod schi_sinqr.mod ocean_sphericalharmonics.mod
schi_direct.mod: schi_direct.f90 ai_kinds.mod screen_centralpotential.mod schi_sinqr.mod ocean_sphericalharmonics.mod

periodic.o:   periodic.f90
periodic.mod: periodic.f90

ocean_sphericalharmonics.o:   ocean_sphericalharmonics.f90 ai_kinds.mod ocean_constants.mod
ocean_sphericalharmonics.mod: ocean_sphericalharmonics.f90 ai_kinds.mod ocean_constants.mod

screen_driver.o: screen_driver.f90 ai_kinds.mod ocean_mpi.mod screen_grid.mod screen_system.mod screen_sites.mod \
                 screen_energy.mod ocean_dft_files.mod screen_wvfn_converter.mod screen_wavefunction.mod \
                 screen_chi_driver.mod screen_centralpotential.mod

screen_driver.x: $(OBJS)
	$(MPIFORT) $(OPTIONS) -o screen_driver.x $(OBJS) $(FFTWL) $(BLAS) $(FFTWL)


clean:
	rm -rf *.mod *.o

