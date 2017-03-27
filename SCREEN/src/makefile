include ../../Makefile.arch

OPTIONS += -stand f03

all: screen_grid.o

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

clean:
	rm -rf *.mod *.o

