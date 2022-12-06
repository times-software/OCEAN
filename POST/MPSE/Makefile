include ../../Makefile.arch

.SUFFIXES: .f .f90 .o .mod .F
binaries = CellVolume.x calcrs.x kkconv.x selfenergy.x
modsrc= m_constants.f90 m_errormod.f90 m_iofiles.f90 m_padio.f90 m_iomod.f90 m_SelfEnergy.f90 conv.f90
modobj=$(modsrc:.f90=.o)
srcf90= polint.f90 xlogx.f90 trap.f90 bpr1_2.f90 bpr2_2.f90 bpr3_2.f90 \
csigma.f90 csigz.f90 logi.f90 omegaq.f90 quinn.f90 str.f90 terp.f90 kkconv.f90 terpc.f90 rdcmt.f90
srcf77= getomi.f IntGrl.f rdloss.f getdom.f 
obj=$(srcf90:.f90=.o)
objf77=$(srcf77:.f=.o)
FLAGS = $(OPTIONS) -fcheck=all -O0

all: $(modobj) $(obj) $(objf77)
	echo $(modobj) $(obj)
	$(FC) $(FLAGS) $(modobj) $(obj) selfenergy.f90 -o selfenergy.x
	$(FC) $(FLAGS) cellvol.f90 -o CellVolume.x
	$(FC) $(FLAGS) rs.f90 -o calcrs.x
	$(FC) $(FLAGS) eps2exc.f terp.o polint.o getomi.o IntGrl.o rdloss.o rdcmt.o getdom.o -o eps2exc.x
	$(FC) $(FLAGS) $(modobj) $(obj) kkprog.f90  -o kkconv.x

install:
	cp CellVolume.x $(INSTDIR)
	cp calcrs.x $(INSTDIR)
	cp nbse_conv.x $(INSTDIR)
	cp nbseconv.script $(INSTDIR)

clean: 
	rm -f *.o
	rm -f *.mod
	rm -f $(binaries)

.f90.o:
	$(FC) -c $(FLAGS) $(@:.o=.f90)
.f.o:
	$(FC) -c $(FLAGS) $(@:.o=.f)
