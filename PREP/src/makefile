include ../../Makefile.arch

OPTIONS +=

all: ocean_prep.x

OBJS = ocean_prep.o prep_system.o prep_wvfn.o ocean_cks.o  ocean_tmels.o

MODOBJS = ../../Modules/ocean_mpi.o ../../Modules/ai_kinds.o ../../Modules/ocean_constants.o \
          ../../Modules/ocean_interpolate.o ../../Modules/periodic.o ../../Modules/screen_timekeeper.o \
          ../../Modules/ocean_legacy_files.o ../../Modules/ocean_qe54_files.o \
          ../../Modules/ocean_qe62_files.o ../../Modules/ocean_dft_files.o ../../Modules/screen_opf.o \
          ../../Modules/ocean_abi_files.o ../../Modules/ocean_phys.o ../../Modules/FFT_wrapper.o \
          ../../Modules/jlmfft.o

install:
	cp ocean_prep.x $(INSTDIR)/

%.o %.mod: %.f90
	$(MPIFORT) $(OPTIONS) -c -o $*.o $< $(FFTWI) -I../../Modules 


prep_wvfn.o:   prep_wvfn.f90 ../../Modules/ocean_dft_files.mod ../../Modules/ocean_mpi.mod \
               ../../Modules/ai_kinds.mod ../../Modules/screen_opf.mod ocean_cks.mod \
               ../../Modules/screen_timekeeper.mod ocean_tmels.mod
prep_wvfn.mod: prep_wvfn.f90 ../../Modules/ocean_dft_files.mod ../../Modules/ocean_mpi.mod \
               ../../Modules/ai_kinds.mod ../../Modules/screen_opf.mod ocean_cks.mod \
               ../../Modules/screen_timekeeper.mod ocean_tmels.mod ../../Modules/fft_wrapper.mod

ocean_cks.o:   ocean_cks.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_interpolate.mod \
               ../../Modules/screen_opf.mod prep_system.mod
ocean_cks.mod: ocean_cks.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_interpolate.mod \
               ../../Modules/screen_opf.mod prep_system.mod

ocean_tmels.o:   ocean_tmels.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod prep_system.mod
ocean_tmels.mod: ocean_tmels.f90 ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod prep_system.mod

prep_system.o:   prep_system.f90 ../../Modules/periodic.mod ../../Modules/ai_kinds.mod ../../Modules/ocean_phys.mod
prep_system.mod: prep_system.f90 ../../Modules/periodic.mod ../../Modules/ai_kinds.mod ../../Modules/ocean_phys.mod

ocean_prep.o: ../../Modules/ai_kinds.mod ../../Modules/ocean_mpi.mod ../../Modules/ocean_dft_files.mod \
              prep_system.mod prep_wvfn.mod


ocean_prep.x: $(OBJS)
	$(MPIFORT) $(OPTIONS) -o ocean_prep.x $(OBJS) $(MODOBJS) $(FFTWL) $(BLAS) $(FFTWL)


clean:
	rm -rf *.mod *.o
