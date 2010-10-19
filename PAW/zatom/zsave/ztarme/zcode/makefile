OPTS = -O
COMP = f90
#
HFOBJS = \
alloy.o amult.o bach.o bintegrate.o config.o \
corepot.o egripper2.o elserf.o ener.o getanu.o \
getexch.o getfg.o getpot.o grandop.o gro.o \
hfca.o hfgi.o hfhead.o hfneck.o hfpa.o \
hfps.o hfrg.o hfvwn.o leadbehv.o melwriter.o \
mkkbfile.o mkkxfile.o mkvaltab.o sdump.o sean1.o \
sean2.o sean4.o sean5.o sig.o spherbes.o \
transp.o vwntst.o weighter.o
#
all: hfk.x 
#
hfk.x: $(HFOBJS)
	$(COMP) $(HFOBJS) -o hfk.x
#
amult.o: amult.f90
	$(COMP) $(OPTS) -c amult.f90
diagq.o: diagq.f90
	$(COMP) $(OPTS) -c diagq.f90
getanu.o: getanu.f90
	$(COMP) $(OPTS) -c getanu.f90
getfg.o: getfg.f90
	$(COMP) $(OPTS) -c getfg.f90
sdump.o: sdump.f90
	$(COMP) $(OPTS) -c sdump.f90
melwriter.o: melwriter.f90
	$(COMP) $(OPTS) -c melwriter.f90
#
.f.o:
	$(COMP) $(OPTS) -c $<
