include ../Makefile.arch

COMP = $(FC)
OPTS = $(OPTIONS) $(LEGACY_REAL)


.SUFFIXES: .f90 .o .x

CONDENSLIST = zobj/condens.o zobj/chkfftreal.o zobj/chkfftrecp.o zobj/jlmfft.o zobj/facpowfind.o \
zobj/gentoreal.o zobj/lxb_getval.o zobj/optim.o zobj/brcapper.o zobj/solidca.o
CONDENS2LIST = zobj/condens2.o zobj/lxb_getval.o zobj/brcapper.o zobj/solidca.o zobj/rho.o
BRDGLIST = zobj/conbbtoug.o zobj/ipack.o zobj/lxb_getval.o
BRLADRLIST = zobj/bridgelad.o zobj/dublinterp.o zobj/lxb_getval.o zobj/dist.o zobj/formx.o zobj/getspacings.o zobj/cross.o
BWAYLIST = zobj/conbbtoux.o zobj/chkfftreal.o zobj/chkfftrecp.o zobj/jlmfft.o zobj/facpowfind.o \
zobj/gentoreal.o zobj/lxb_getval.o zobj/optim.o zobj/formx.o
GWAYLIST = zobj/conugtoux.o zobj/chkfftreal.o zobj/chkfftrecp.o zobj/jlmfft.o zobj/facpowfind.o \
zobj/gentoreal.o zobj/lxb_getval.o zobj/optim.o zobj/formx.o zobj/FFT_wrapper.o
SHIFTLIST = zobj/shifted_conugtoux.o zobj/chkfftreal.o zobj/chkfftrecp.o zobj/jlmfft.o zobj/facpowfind.o \
zobj/gentoreal.o zobj/lxb_getval.o zobj/optim.o zobj/formx.o zobj/FFT_wrapper.o
ORTHLIST = zobj/orthog.o zobj/normalize.o zobj/lxb_getval.o zobj/dptest.o
SETUP2LIST = zobj/setup2.o zobj/enkread.o zobj/lxb_getval.o zobj/melread.o zobj/cross.o zobj/getrecvec.o
LLOBJS = zobj/ll.o zobj/lxb_getval.o zobj/levlou.o
XEQYOBJS = zobj/whom0.o
SPECTOBJS = zobj/spectrum.o zobj/lxb_getval.o
#

all: zexe/setup2.x zexe/orthog.x zexe/conugtoux.x zobj/FFT_wrapper.o zexe/shifted_conugtoux.x zexe/jdftx_density.x
#

zobj/FFT_wrapper.o: zf90/FFT_wrapper.f90
	$(COMP) $(OPTS) -c -cpp zf90/FFT_wrapper.f90 $(FFTWI) -o zobj/FFT_wrapper.o

zexe/bridgegw.x: zobj/bridgegw.o zobj/lxb_getval.o
	$(COMP) $(OPTS) -o zexe/bridgegw.x zobj/bridgegw.o zobj/lxb_getval.o
zexe/bridgelad.x: $(BRLADRLIST)
	$(COMP) $(OPTS) -o zexe/bridgelad.x $(BRLADRLIST)
zexe/conbbtoug.x: $(BRDGLIST)
	$(COMP) $(OPTS) -o zexe/conbbtoug.x $(BRDGLIST)
zexe/conbbtoux.x: $(BWAYLIST)
	$(COMP) $(OPTS) -o zexe/conbbtoux.x $(BWAYLIST)
zexe/conugtoux.x: $(GWAYLIST)
	$(COMP) $(OPTS) -o zexe/conugtoux.x $(GWAYLIST) $(FFTWI) $(FFTWL)
zexe/shifted_conugtoux.x: $(SHIFTLIST)
	$(COMP) $(OPTS) -o zexe/shifted_conugtoux.x $(SHIFTLIST) $(FFTWI) $(FFTWL)
zexe/condens.x: $(CONDENSLIST)
	$(COMP) $(OPTS) -o zexe/condens.x $(CONDENSLIST)
zexe/condens2.x: $(CONDENS2LIST)
	$(COMP) $(OPTS) -o zexe/condens2.x $(CONDENS2LIST)
zexe/ll.x: $(LLOBJS)
	$(COMP) $(OPTS) -o zexe/ll.x $(LLOBJS)
zexe/orthog.x: $(ORTHLIST)
	$(COMP) $(OPTS) -o zexe/orthog.x $(ORTHLIST)
zexe/spect.x: $(SPECTOBJS)
	$(COMP) $(OPTS) -o zexe/spect.x $(SPECTOBJS)
zexe/setup2.x: $(SETUP2LIST)
	$(COMP) $(OPTS) -o zexe/setup2.x $(SETUP2LIST)
zexe/whom0.x: $(XEQYOBJS)
	$(COMP) $(OPTS) -o zexe/whom0.x $(XEQYOBJS)
zexe/jdftx_density.x: zf90/jdftx_density.f90
	$(COMP) $(OPTS) -o zexe/jdftx_density.x zf90/jdftx_density.f90
#

zobj/setup2.o: zf90/setup2.f90
	$(COMP) $(OPTS) -c zf90/setup2.f90 -o zobj/setup2.o
zobj/enkread.o: zf90/enkread.f90
	$(COMP) $(OPTS) -c zf90/enkread.f90 -o zobj/enkread.o
zobj/melread.o: zf90/melread.f90
	$(COMP) $(OPTS) -c zf90/melread.f90 -o zobj/melread.o
zobj/cross.o: zf90/cross.f90
	$(COMP) $(OPTS) -c zf90/cross.f90 -o zobj/cross.o
zobj/getrecvec.o: zf90/getrecvec.f90
	$(COMP) $(OPTS) -c zf90/getrecvec.f90 -o zobj/getrecvec.o
zobj/getspacings.o: zf90/getspacings.f90
	$(COMP) $(OPTS) -c zf90/getspacings.f90 -o zobj/getspacings.o
#
zobj/ll.o: zf77/ll.f
	$(COMP) $(OPTS) -c zf77/ll.f -o zobj/ll.o
zobj/levlou.o: zf77/levlou.f
	$(COMP) $(OPTS) -c zf77/levlou.f -o zobj/levlou.o
zobj/spectrum.o: zf77/spectrum.f
	$(COMP) $(OPTS) -c zf77/spectrum.f -o zobj/spectrum.o
#
zobj/brcapper.o: zf90/brcapper.f90
	$(COMP) $(OPTS) -c zf90/brcapper.f90 -o zobj/brcapper.o
zobj/bridgegw.o: zf90/bridgegw.f90
	$(COMP) $(OPTS) -c zf90/bridgegw.f90 -o zobj/bridgegw.o
zobj/bridgelad.o: zf90/bridgelad.f90
	$(COMP) $(OPTS) -c zf90/bridgelad.f90 -o zobj/bridgelad.o
zobj/chkfftreal.o: zf90/chkfftreal.f90
	$(COMP) $(OPTS) -c zf90/chkfftreal.f90 -o zobj/chkfftreal.o
zobj/chkfftrecp.o: zf90/chkfftrecp.f90
	$(COMP) $(OPTS) -c zf90/chkfftrecp.f90 -o zobj/chkfftrecp.o
zobj/conbbtoug.o: zf90/conbbtoug.f90
	$(COMP) $(OPTS) -c zf90/conbbtoug.f90 -o zobj/conbbtoug.o
zobj/conbbtoux.o: zf90/conbbtoux.f90
	$(COMP) $(OPTS) -c zf90/conbbtoux.f90 -o zobj/conbbtoux.o
zobj/conugtoux.o: zf90/conugtoux.f90 zobj/FFT_wrapper.o
	$(COMP) $(OPTS) -c zf90/conugtoux.f90 -o zobj/conugtoux.o
zobj/shifted_conugtoux.o: zf90/shifted_conugtoux.f90 zobj/FFT_wrapper.o
	$(COMP) $(OPTS) -c zf90/shifted_conugtoux.f90 -o zobj/shifted_conugtoux.o
zobj/condens.o: zf90/condens.f90
	$(COMP) $(OPTS) -c zf90/condens.f90 -o zobj/condens.o
zobj/condens2.o: zf90/condens2.f90
	$(COMP) $(OPTS) -c zf90/condens2.f90 -o zobj/condens2.o
zobj/rho.o: zf90/rho.f90
	$(COMP) $(OPTS) -c zf90/rho.f90 -o zobj/rho.o
zobj/solidca.o: zstd/solidca.f
	$(COMP) $(OPTS) -c zstd/solidca.f -o zobj/solidca.o
zobj/dist.o: zf90/dist.f90
	$(COMP) $(OPTS) -c zf90/dist.f90 -o zobj/dist.o
zobj/dptest.o: zf90/dptest.f90
	$(COMP) $(OPTS) -c zf90/dptest.f90 -o zobj/dptest.o
zobj/dublinterp.o: zf90/dublinterp.f90
	$(COMP) $(OPTS) -c zf90/dublinterp.f90 -o zobj/dublinterp.o
zobj/facpowfind.o: zf90/facpowfind.f90
	$(COMP) $(OPTS) -c zf90/facpowfind.f90 -o zobj/facpowfind.o
zobj/formx.o: zf90/formx.f90
	$(COMP) $(OPTS) -c zf90/formx.f90 -o zobj/formx.o
zobj/gentoreal.o: zf90/gentoreal.f90
	$(COMP) $(OPTS) -c zf90/gentoreal.f90 -o zobj/gentoreal.o
zobj/ipack.o: zf90/ipack.f90
	$(COMP) $(OPTS) -c zf90/ipack.f90 -o zobj/ipack.o
zobj/jlmfft.o: zstd/jlmfft.f
	$(COMP) $(OPTS) -c zstd/jlmfft.f -o zobj/jlmfft.o
zobj/lxb_getval.o: zf90/lxb_getval.f90
	$(COMP) $(OPTS) -c zf90/lxb_getval.f90 -o zobj/lxb_getval.o
zobj/normalize.o: zf90/normalize.f90
	$(COMP) $(OPTS) -c zf90/normalize.f90 -o zobj/normalize.o
zobj/optim.o: zf90/optim.f90
	$(COMP) $(OPTS) -c zf90/optim.f90 -o zobj/optim.o
zobj/orthog.o: zf90/orthog.f90
	$(COMP) $(OPTS) -c zf90/orthog.f90 -o zobj/orthog.o
zobj/whom0.o: zf77/whom0.f
	$(COMP) $(OPTS) -c zf77/whom0.f -o zobj/whom0.o

clean:
	rm -f zexe/setup2.x zexe/condens.x zexe/condens2.x zexe/bridgelad.x \
              zexe/conbbtoux.x zexe/orthog.x zexe/conbbtoug.x \
              zexe/conugtoux.x zexe/bridgegw.x zexe/whom0.x \
              zexe/ll.x zexe/spect.x zobj/bridgegw.o $(CONDENSLIST) \
              $(BRDGLIST) $(BRLADRLIST) $(BWAYLIST) $(GWAYLIST) $(ORTHLIST) \
              $(SETUP2LIST) $(LLOBJS) $(XEQYOBJS) $(SPECTOBJS) $(CONDENS2LIST) zobj/FFT_wrapper.o
	rm -f *__genmod.???
	rm -f *.mod

install:
	cp zexe/*.x $(INSTDIR)

instdev:
	for F in zexe/*.x; do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

