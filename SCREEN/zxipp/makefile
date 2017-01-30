#
OBXPL = zobj/xippsin_lm.o zobj/xeqtanx.o zobj/rinvert.o zobj/mkvipt.o zobj/intval.o
OBXPB = zobj/xippsin_box.o zobj/xeqtanx.o zobj/rinvert.o zobj/mkvipt.o zobj/intval.o
OBXP2 = zobj/xippsin.o zobj/xeqtanx.o zobj/rinvert.o zobj/mkvipt.o zobj/intval.o
OBVHO = zobj/vhommod.o zobj/levlou.o zobj/sphj0.o zobj/sphj1.o zobj/getomega.o
OBRSC = zobj/rscombine.o zobj/intval.o
OBXMN = zobj/mkximunu.o zobj/formreytab.o zobj/getprefs.o zobj/ylmeval.o
OBDRX = zobj/dressximunu.o zobj/intval.o zobj/dftder3.o zobj/getc2.o zobj/cacorr.o zobj/cinvert.o
#
TARGETS = zexe/xipps.x zexe/vhommod.x zexe/rscombine.x  zexe/mkximunu.x zexe/dressximunu.x
all: $(TARGETS)
#
#
zexe/xipps_lm.x: $(OBXPL)
	$(COMP) $(OPTS) $(OBXPL) -o zexe/xipps_lm.x
zexe/xipps_box.x: $(OBXPB)
	$(COMP) $(OPTS) $(OBXPB) -o zexe/xipps_box.x -lmkl_intel_lp64  -lmkl_sequential -lmkl_core
zexe/xipps.x: $(OBXP2)
	$(COMP) $(OPTS) $(OBXP2) -o zexe/xipps.x
zexe/rscombine.x: $(OBRSC)
	$(COMP) $(OPTS) $(OBRSC) -o zexe/rscombine.x
zexe/vhommod.x: $(OBVHO)
	$(COMP) $(OPTS) $(OBVHO) -o zexe/vhommod.x
zexe/mkximunu.x: $(OBXMN)
	$(COMP) $(OPTS) $(OBXMN) -o zexe/mkximunu.x
#
zobj/xippsin_box.o: zf90/xipp_box.f90
	$(COMP) $(OPTS) -c zf90/xipp_box.f90 -o zobj/xippsin_box.o
zobj/xippsin_lm.o: zf90/xippsin_lm.f90
	$(COMP) $(OPTS) -c zf90/xippsin_lm.f90 -o zobj/xippsin_lm.o
zobj/getomega.o: zf90/getomega.f90
	$(COMP) $(OPTS) -c zf90/getomega.f90 -o zobj/getomega.o
zobj/intval.o: zf90/intval.f90
	$(COMP) $(OPTS) -c zf90/intval.f90 -o zobj/intval.o
zobj/levlou.o: zf90/levlou.f90
	$(COMP) $(OPTS) -c zf90/levlou.f90 -o zobj/levlou.o
zobj/mkvipt.o: zf90/mkvipt.f90
	$(COMP) $(OPTS) -c zf90/mkvipt.f90 -o zobj/mkvipt.o
zobj/rinvert.o: zf90/rinvert.f90
	$(COMP) $(OPTS) -c zf90/rinvert.f90 -o zobj/rinvert.o
zobj/rscombine.o: zf90/rscombine.f90
	$(COMP) $(OPTS) -c zf90/rscombine.f90 -o zobj/rscombine.o
zobj/sphj0.o: zf90/sphj0.f90
	$(COMP) $(OPTS) -c zf90/sphj0.f90 -o zobj/sphj0.o
zobj/sphj1.o: zf90/sphj1.f90
	$(COMP) $(OPTS) -c zf90/sphj1.f90 -o zobj/sphj1.o
zobj/vhommod.o: zf90/vhommod.f90
	$(COMP) $(OPTS) -c zf90/vhommod.f90 -o zobj/vhommod.o
zobj/xippsin.o: zf90/xippsin.f90
	$(COMP) $(OPTS) -c zf90/xippsin.f90 -o zobj/xippsin.o
zobj/xeqtanx.o: zf90/xeqtanx.f90
	$(COMP) $(OPTS) -c zf90/xeqtanx.f90 -o zobj/xeqtanx.o
zobj/mkximunu.o: zf90/mkximunu.f90
	$(COMP) $(OPTS) -c zf90/mkximunu.f90 -o zobj/mkximunu.o
zobj/formreytab.o: zf90/formreytab.f90
	$(COMP) $(OPTS) -c zf90/formreytab.f90 -o zobj/formreytab.o
zobj/getprefs.o: zf90/getprefs.f90
	$(COMP) $(OPTS) -c zf90/getprefs.f90 -o zobj/getprefs.o
zobj/ylmeval.o: zf90/ylmeval.f90
	$(COMP) $(OPTS) -c zf90/ylmeval.f90 -o zobj/ylmeval.o
zobj/dressximunu.o: zf90/dressximunu.f90
	$(COMP) $(OPTS) -c zf90/dressximunu.f90 -o zobj/dressximunu.o:
zobj/cinvert.o: zf90/cinvert.f90
	$(COMP) $(OPTS) -c zf90/cinvert.f90 -o zobj/cinvert.o
zobj/dftder3.o: zf90/dftder3.f90
	$(COMP) $(OPTS) -c zf90/dftder3.f90 -o zobj/dftder3.o
zobj/cacorr.o: zf90/cacorr.f
	$(COMP) $(OPTS) -c zf90/cacorr.f -o zobj/cacorr.o
zobj/getc2.o: zf90/getc2.f
	$(COMP) $(OPTS) -c zf90/getc2.f -o zobj/getc2.o

install:
	cp $(TARGETS) $(INSTDIR)/

clean:
	rm -f zobj/*.o
	rm -f zexe/*.x
