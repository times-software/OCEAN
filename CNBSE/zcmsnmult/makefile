OPTS = -O2   #-YDEALLOC=ALL
#COMP = /usr/local/usr/absoft32/bin/f90
COMP = ifort
#
CAINM = zobj/cainmultip.o zobj/cainxact.o zobj/getabb.o zobj/haydump.o zobj/fgact.o zobj/ctact.o zobj/soact.o \
        zobj/invdrv.o zobj/rtov.o zobj/vtor.o zobj/sizereport.o \
        zobj/nbsemkcmel.o zobj/newsoact.o zobj/soprep.o zobj/limel.o zobj/jimel.o \
        zobj/newgetylm.o zobj/newgetprefs.o zobj/cainmhsetup.o zobj/newthreey.o zobj/intval.o zobj/redtrid.o \
        zobj/getomega.o zobj/elsdch.o zobj/gamfcn.o zobj/cinv.o zobj/cainchact.o zobj/cainxrvmult.o zobj/loadux.o \
        zobj/cainnqtoxr.o zobj/cainxrtonq.o zobj/jlmfft.o
#
all: zexe/cainmultip.x
#
zexe/cainmultip.x: $(CAINM)
	$(COMP) $(OPTS) $(CAINM) -o zexe/cainmultip.x
#
zobj/cainchact.o: zf90/cainchact.f90
	$(COMP) -c zf90/cainchact.f90 -o zobj/cainchact.o
zobj/cainmhsetup.o: zf90/cainmhsetup.f90 sphsetnx.h.f90 sphsetx.h.f90
	$(COMP) -c zf90/cainmhsetup.f90 -o zobj/cainmhsetup.o
zobj/cainmultip.o: zf90/cainmultip.f90 cainmusetupx.h nbsemusetupnx.h
	$(COMP) $(OPTS) -c zf90/cainmultip.f90 -o zobj/cainmultip.o
zobj/cainnqtoxr.o: zf90/cainnqtoxr.f90
	$(COMP) $(OPTS) -c zf90/cainnqtoxr.f90 -o zobj/cainnqtoxr.o
zobj/cainxact.o: zf90/cainxact.f90
	$(COMP) $(OPTS) -c zf90/cainxact.f90 -o zobj/cainxact.o
zobj/cainxrtonq.o: zf90/cainxrtonq.f90
	$(COMP) $(OPTS) -c zf90/cainxrtonq.f90 -o zobj/cainxrtonq.o
zobj/cainxrvmult.o: zf90/cainxrvmult.f90
	$(COMP) $(OPTS) -c zf90/cainxrvmult.f90 -o zobj/cainxrvmult.o
zobj/cinv.o: zstd/cinv.f
	$(COMP) $(OPTS) -c zstd/cinv.f -o zobj/cinv.o
zobj/ctact.o: zf90/ctact.f90
	$(COMP) $(OPTS) -c zf90/ctact.f90 -o zobj/ctact.o
zobj/elsdch.o: zstd/elsdch.f
	$(COMP) $(OPTS) -c zstd/elsdch.f -o zobj/elsdch.o
zobj/fgact.o: zf90/fgact.f90
	$(COMP) $(OPTS) -c zf90/fgact.f90 -o zobj/fgact.o
zobj/fsassem.o: zf90/fsassem.f90
	$(COMP) $(OPTS) -c zf90/fsassem.f90 -o zobj/fsassem.o
zobj/fshift.o: zf90/fshift.f90
	$(COMP) $(OPTS) -c zf90/fshift.f90 -o zobj/fshift.o
zobj/gamfcn.o: zf90/gamfcn.f90
	$(COMP) $(OPTS) -c zf90/gamfcn.f90 -o zobj/gamfcn.o
zobj/getabb.o: zf90/getabb.f90
	$(COMP) $(OPTS) -c zf90/getabb.f90 -o zobj/getabb.o
zobj/getomega.o: zf90/getomega.f90
	$(COMP) $(OPTS) -c zf90/getomega.f90 -o zobj/getomega.o
zobj/haydump.o: zf90/haydump.f90
	$(COMP) $(OPTS) -c zf90/haydump.f90 -o zobj/haydump.o
zobj/intval.o: zf90/intval.f90
	$(COMP) $(OPTS) -c zf90/intval.f90 -o zobj/intval.o
zobj/invdrv.o: zf90/invdrv.f90
	$(COMP) $(OPTS) -c zf90/invdrv.f90 -o zobj/invdrv.o
zobj/jimel.o: zf90/jimel.f90
	$(COMP) $(OPTS) -c zf90/jimel.f90 -o zobj/jimel.o
zobj/jlmfft.o: zstd/jlmfft.f
	$(COMP) $(OPTS) -c zstd/jlmfft.f -o zobj/jlmfft.o
zobj/limel.o: zf90/limel.f90
	$(COMP) $(OPTS) -c zf90/limel.f90 -o zobj/limel.o
zobj/loadux.o: zf90/loadux.f90
	$(COMP) $(OPTS) -c zf90/loadux.f90 -o zobj/loadux.o
zobj/mktmel.o: zf90/mktmel.f90
	$(COMP) $(OPTS) -c zf90/mktmel.f90 -o zobj/mktmel.o
zobj/nbsemkcmel.o: zf90/nbsemkcmel.f90
	$(COMP) $(OPTS) -c zf90/nbsemkcmel.f90 -o zobj/nbsemkcmel.o
zobj/newgetprefs.o: zf90/newgetprefs.f90
	$(COMP) $(OPTS) -c zf90/newgetprefs.f90 -o zobj/newgetprefs.o
zobj/newgetylm.o: zf90/newgetylm.f90
	$(COMP) $(OPTS) -c zf90/newgetylm.f90 -o zobj/newgetylm.o
zobj/newsoact.o: zf90/newsoact.f90
	$(COMP) $(OPTS) -c zf90/newsoact.f90 -o zobj/newsoact.o
zobj/newthreey.o: zf90/newthreey.f90
	$(COMP) $(OPTS) -c zf90/newthreey.f90 -o zobj/newthreey.o
zobj/redtrid.o: zf77/redtrid.f
	$(COMP) $(OPTS) -c zf77/redtrid.f -o zobj/redtrid.o
zobj/rtov.o: zf90/rtov.f90
	$(COMP) $(OPTS) -c zf90/rtov.f90 -o zobj/rtov.o
zobj/sizereport.o: zf90/sizereport.f90
	$(COMP) $(OPTS) -c zf90/sizereport.f90 -o zobj/sizereport.o
zobj/soact.o: zf90/soact.f90
	$(COMP) $(OPTS) -c zf90/soact.f90 -o zobj/soact.o
zobj/soprep.o: zf90/soprep.f90
	$(COMP) $(OPTS) -c zf90/soprep.f90 -o zobj/soprep.o
zobj/vtor.o: zf90/vtor.f90
	$(COMP) $(OPTS) -c zf90/vtor.f90 -o zobj/vtor.o

install:
	cp zexe/cainmultip.x $(INSTDIR)

clean:
	rm -f zobj/*.o
	rm -f zexe/*.x
