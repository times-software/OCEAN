include ../../Makefile.arch

OPTIONS +=

all: ocean_prep.x

OBJS = ocean_mpi.o ai_kinds.o ocean_legacy_files.o ocean_qe54_files.o ocean_qe62_files.o \
       screen_timekeeper.o ocean_dft_files.o ocean_prep.o prep_system.o prep_wvfn.o

install:
	cp ocean_prep.x $(INSTDIR)/

%.o %.mod: %.f90
	$(MPIFORT) $(OPTIONS) -c -o $*.o $< $(FFTWI)

ocean_mpi.o: ocean_mpi.f90 ai_kinds.mod
ocean_mpi.mod: ocean_mpi.f90 ai_kinds.mod

ai_kinds.o: ai_kinds.f90
ai_kinds.mod: ai_kinds.f90

ocean_constants.o: ocean_constants.f90 ai_kinds.mod
ocean_constants.mod: ocean_constants.f90 ai_kinds.mod

ocean_legacy_files.o: ocean_legacy_files.f90 ai_kinds.mod ocean_mpi.mod
ocean_legacy_files.mod: ocean_legacy_files.f90 ai_kinds.mod ocean_mpi.mod

ocean_qe54_files.o:   ai_kinds.mod ocean_mpi.mod screen_timekeeper.mod
ocean_qe54_files.mod: ai_kinds.mod ocean_mpi.mod screen_timekeeper.mod

ocean_qe62_files.o:   ai_kinds.mod ocean_mpi.mod screen_timekeeper.mod
ocean_qe62_files.mod: ai_kinds.mod ocean_mpi.mod screen_timekeeper.mod

ocean_dft_files.o: ocean_dft_files.f90 ai_kinds.mod ocean_mpi.mod ocean_legacy_files.mod ocean_qe54_files.mod \
                   ocean_qe62_files.mod
ocean_dft_files.mod: ocean_dft_files.f90 ai_kinds.mod ocean_mpi.mod ocean_legacy_files.mod ocean_qe54_files.mod \
                     ocean_qe62_files.mod

prep_wvfn.o:   prep_wvfn.f90 ocean_dft_files.mod ocean_mpi.mod ai_kinds.mod 
prep_wvfn.mod: prep_wvfn.f90 ocean_dft_files.mod ocean_mpi.mod ai_kinds.mod


ocean_prep.o: ai_kinds.mod ocean_mpi.mod ocean_dft_files.mod prep_system.mod prep_wvfn.mod


ocean_prep.x: $(OBJS)
	$(MPIFORT) $(OPTIONS) -o ocean_prep.x $(OBJS) $(FFTWL) $(BLAS) $(FFTWL)


clean:
	rm -rf *.mod *.o
