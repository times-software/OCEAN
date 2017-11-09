include ../../Makefile.arch

OPTIONS += -stand f03


all: screen_grid.o screen_system.o screen_driver.x

OBJS = screen_driver.o ai_kinds.o ocean_mpi.o screen_grid.o screen_system.o screen_sites.o screen_paral.o \
       screen_energy.o

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

screen_sites.o: screen_sites.f90 screen_paral.mod screen_grid.mod ocean_mpi.mod ai_kinds.mod
screen_sites.mod: screen_sites.f90 screen_paral.mod screen_grid.mod ocean_mpi.mod ai_kinds.mod

screen_energy.o: screen_energy.f90 ai_kinds.mod ocean_mpi.mod screen_system.mod
screen_energy.mod: screen_energy.f90 ai_kinds.mod ocean_mpi.mod screen_system.mod

ocean_prep.o: ocean_prep.f90 ocean_legacy_files.mod ai_kinds.mod ocean_mpi.mod
ocean_prep.mod: ocean_prep.f90 ocean_legacy_files.mod ai_kinds.mod ocean_mpi.mod

ocean_legacy_files.o: ocean_legacy_files.f90 ai_kinds.mod ocean_mpi.mod
ocean_legacy_files.mod: ocean_legacy_files.f90 ai_kinds.mod ocean_mpi.mod

screen_driver.o: screen_driver.f90 ai_kinds.mod ocean_mpi.mod screen_grid.mod screen_system.mod screen_sites.mod \
                 screen_energy.mod

screen_driver.x: $(OBJS)
	$(MPIFORT) $(OPTIONS) -o screen_driver.x $(OBJS)


clean:
	rm -rf *.mod *.o

