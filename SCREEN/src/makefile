include ../../Makefile.arch

#OPTIONS += -prof-gen
#OPTIONS += -prof-use 
OPTIONS += -march=haswell -fopenmp 
#OPTIONS += -check all -g -traceback


all: screen_driver.x

OBJS = screen_driver.o screen_grid.o screen_system.o screen_sites.o screen_paral.o \
       screen_energy.o screen_wavefunction.o screen_wvfn_converter.o \
       screen_chi0.o screen_chi_driver.o screen_chi.o screen_centralpotential.o schi_sinqr.o \
       schi_direct.o ocean_sphericalharmonics.o ocean_quadrature.o 

MODOBJS = ../../Modules/ocean_mpi.o ../../Modules/ai_kinds.o ../../Modules/ocean_constants.o \
          ../../Modules/ocean_interpolate.o ../../Modules/periodic.o ../../Modules/screen_timekeeper.o \
          ../../Modules/ocean_legacy_files.o ../../Modules/ocean_qe54_files.o \
          ../../Modules/ocean_qe62_files.o ../../Modules/ocean_dft_files.o ../../Modules/screen_opf.o

install:
	cp screen_driver.x $(INSTDIR)/

%.o %.mod: %.f90
	$(MPIFORT) $(OPTIONS) -c -o $*.o $< $(FFTWI) -I../../Modules


screen_grid.o: screen_grid.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_constants.mod ../../Modules/ocean_mpi.mod ocean_quadrature.mod
screen_grid.mod: screen_grid.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_constants.mod ../../Modules/ocean_mpi.mod ocean_quadrature.mod

screen_system.o: screen_system.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod
screen_system.mod: screen_system.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod

screen_paral.o: screen_paral.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod
screen_paral.mod: screen_paral.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod

screen_sites.o: screen_sites.f90 screen_paral.mod screen_grid.mod ../../Modules/ocean_mpi.mod ../../Modules/ai_kinds.mod screen_wavefunction.mod \
                ../../Modules/periodic.mod
screen_sites.mod: screen_sites.f90 screen_paral.mod screen_grid.mod ../../Modules/ocean_mpi.mod ../../Modules/ai_kinds.mod screen_wavefunction.mod \
                  ../../Modules/periodic.mod

screen_energy.o: screen_energy.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod screen_system.mod ../../Modules/ocean_dft_files.mod
screen_energy.mod: screen_energy.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod screen_system.mod ../../Modules/ocean_dft_files.mod

ocean_prep.o: ocean_prep.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod
ocean_prep.mod: ocean_prep.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod

screen_wavefunction.o:   screen_wavefunction.f90 ../../Modules/ai_kinds.mod screen_system.mod screen_paral.mod
screen_wavefunction.mod: screen_wavefunction.f90 ../../Modules/ai_kinds.mod screen_system.mod screen_paral.mod

screen_wvfn_converter.o:   screen_wvfn_converter.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_dft_files.mod  \
                           screen_sites.mod screen_grid.mod screen_paral.mod ../../Modules/screen_opf.mod \
                           ../../Modules/screen_timekeeper.mod ../../Modules/ocean_interpolate.mod
screen_wvfn_converter.mod: screen_wvfn_converter.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_dft_files.mod  \
                           screen_sites.mod screen_grid.mod screen_paral.mod ../../Modules/screen_opf.mod \
                           ../../Modules/screen_timekeeper.mod ../../Modules/ocean_interpolate.mod

screen_chi0.o:   screen_chi0.f90 ../../Modules/ai_kinds.mod screen_paral.mod screen_sites.mod screen_grid.mod \
                 screen_wavefunction.mod ocean_quadrature.mod ../../Modules/screen_timekeeper.mod
screen_chi0.mod: screen_chi0.f90 ../../Modules/ai_kinds.mod screen_paral.mod screen_sites.mod screen_grid.mod \
                 screen_wavefunction.mod ocean_quadrature.mod ../../Modules/screen_timekeeper.mod

screen_chi.o:   screen_chi.f90 ../../Modules/ai_kinds.mod screen_grid.mod ../../Modules/ocean_mpi.mod ../../Modules/ocean_constants.mod schi_sinqr.mod \
                screen_centralpotential.mod schi_direct.mod
screen_chi.mod: screen_chi.f90 ../../Modules/ai_kinds.mod screen_grid.mod ../../Modules/ocean_mpi.mod ../../Modules/ocean_constants.mod schi_sinqr.mod \
                screen_centralpotential.mod schi_direct.mod

screen_chi_driver.o:   screen_chi_driver.f90 screen_chi0.mod screen_chi.mod screen_centralpotential.mod \
                       screen_grid.mod ../../Modules/screen_timekeeper.mod
screen_chi_driver.mod: screen_chi_driver.f90 screen_chi0.mod screen_chi.mod screen_centralpotential.mod \
                       screen_grid.mod ../../Modules/screen_timekeeper.mod

screen_centralpotential.o:   screen_centralpotential.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod screen_system.mod
screen_centralpotential.mod: screen_centralpotential.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod screen_system.mod

schi_sinqr.o:   schi_sinqr.f90 ../../Modules/ai_kinds.mod screen_centralpotential.mod
schi_sinqr.mod: schi_sinqr.f90 ../../Modules/ai_kinds.mod screen_centralpotential.mod

schi_direct.o:   schi_direct.f90 ../../Modules/ai_kinds.mod screen_centralpotential.mod schi_sinqr.mod ocean_sphericalharmonics.mod
schi_direct.mod: schi_direct.f90 ../../Modules/ai_kinds.mod screen_centralpotential.mod schi_sinqr.mod ocean_sphericalharmonics.mod


ocean_quadrature.o:   ocean_quadrature.f90 ../../Modules/ai_kinds.mod
ocean_quadrature.mod: ocean_quadrature.f90 ../../Modules/ai_kinds.mod

ocean_sphericalharmonics.o:   ocean_sphericalharmonics.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_constants.mod
ocean_sphericalharmonics.mod: ocean_sphericalharmonics.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_constants.mod

screen_driver.o: screen_driver.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod screen_grid.mod screen_system.mod screen_sites.mod \
                 screen_energy.mod ../../Modules/ocean_dft_files.mod screen_wvfn_converter.mod screen_wavefunction.mod \
                 screen_chi_driver.mod screen_centralpotential.mod ../../Modules/screen_opf.mod ../../Modules/screen_timekeeper.mod

screen_driver.x: $(OBJS)
	$(MPIFORT) $(OPTIONS) -o screen_driver.x $(OBJS) $(MODOBJS) $(FFTWL) $(BLAS) $(FFTWL)


clean:
	rm -rf *.mod *.o

