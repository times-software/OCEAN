include ../Makefile.arch

src = $(wildcard *.f90)
obj = $(src:.f90=.o)
mod = $(src:.f90=.mod)

all: $(obj) $(mod)

%.o %.mod: %.f90
	$(MPIFORT) $(OPTIONS) -c -o $*.o $< $(FFTWI)

ocean_constants.o: ai_kinds.mod
ocean_constants.mod: ai_kinds.mod

ocean_interpolate.o:   ai_kinds.mod
ocean_interpolate.mod: ai_kinds.mod

ocean_mpi.o:   ai_kinds.mod
ocean_mpi.mod: ai_kinds.mod

ocean_legacy_files.o:   ai_kinds.mod ocean_mpi.mod
ocean_legacy_files.mod: ai_kinds.mod ocean_mpi.mod

ocean_qe54_files.o:   ai_kinds.mod ocean_mpi.mod screen_timekeeper.mod
ocean_qe54_files.mod: ai_kinds.mod ocean_mpi.mod screen_timekeeper.mod

ocean_qe62_files.o:   ai_kinds.mod ocean_mpi.mod screen_timekeeper.mod
ocean_qe62_files.mod: ai_kinds.mod ocean_mpi.mod screen_timekeeper.mod

ocean_dft_files.o:   ai_kinds.mod ocean_mpi.mod ocean_legacy_files.mod ocean_qe54_files.mod \
                     ocean_qe62_files.mod
ocean_dft_files.mod: ai_kinds.mod ocean_mpi.mod ocean_legacy_files.mod ocean_qe54_files.mod \
                     ocean_qe62_files.mod

screen_opf.o:   ai_kinds.mod ocean_mpi.mod
screen_opf.mod: ai_kinds.mod ocean_mpi.mod

.PHONY: clean
clean:
	rm -f $(obj) $(mod)
