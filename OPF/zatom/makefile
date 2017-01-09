#
HFOBJS = zobj/hfhead.o zobj/hfhead2.o zobj/hfneck.o \
         zobj/mkkbfile.o zobj/hfpa.o zobj/getpot.o zobj/sig.o \
         \
         zobj/pseudo.o zobj/pseudize.o zobj/ffitx0.o zobj/unscreen.o zobj/mklfc.o \
         \
         zobj/hfrg.o zobj/hfvwn.o zobj/bach.o zobj/gro.o \
         \
         zobj/elsolve.o zobj/radnorm.o zobj/augment.o zobj/elener.o \
         zobj/bcsolve.o \
         \
         zobj/vwntst.o zobj/hfca.o zobj/grandop.o zobj/hfgi.o \
         zobj/elserf.o zobj/mkvaltab.o zobj/corepot.o \
         zobj/recon.o zobj/dipole.o zobj/ssubsh.o zobj/longit.o \
         zobj/config.o zobj/sdump.o \
         zobj/sean1.o zobj/sean2.o zobj/sean4.o \
         zobj/sean5.o zobj/amult.o zobj/transp.o zobj/spherbes.o \
         zobj/leadbehv.o zobj/mkkxfile.o zobj/alloy.o zobj/weighter.o \
         zobj/getkap.o zobj/der1.o zobj/der2.o zobj/besft.o zobj/besjl.o zobj/xidump.o \
         zobj/escanner.o zobj/integ.o zobj/integi.o zobj/intego.o zobj/setxkdk.o \
         zobj/ppdump.o zobj/ppload.o zobj/intval.o zobj/dft.o zobj/setqmm.o zobj/getplead.o \
         zobj/orbmom.o zobj/bintegrate.o zobj/radint.o \
         zobj/freshen.o zobj/rpower.o zobj/mkcorcon.o \
         zobj/getexch.o zobj/screencore.o zobj/orbcont.o \
         zobj/addpot.o zobj/potfigure.o zobj/optript.o zobj/chgocc.o \
         \
         zobj/spartanfip.o zobj/baregrip.o zobj/orthred.o \
         zobj/projdumper.o zobj/escanprep.o zobj/normangnodes.o \
         zobj/reconstruct.o zobj/optradf.o zobj/getfgnew.o zobj/fgcalc.o zobj/getmeznl.o \
         \
         zobj/continuum.o zobj/projaug.o zobj/trck.o
#
#        zobj/fillinpaw.o zobj/egripper2.o 
#        zobj/fipfront.o zobj/pwavesetter.o zobj/newdumper.o 
#        zobj/formfkgknll.o zobj/coullgen.o zobj/formmelnll.o zobj/powlgen.o
#        zobj/coreload.o zobj/pawload.o zobj/geticore.o
#
HPOBJS = zobj/hartpart.o
BTOBJS = zobj/backtran.o
FTOBJS = zobj/elsft.o
KBOBJS = zobj/kbft.o zobj/gro.o zobj/diagq.o zobj/elsdch.o
LLOBJS = zobj/llscreen.o zobj/levlou.o
RTOBJS = zobj/rtabber.o
RTOBJ2 = zobj/rtab2.o
SEOBJS = zobj/kbbessel.o zobj/besjl.o zobj/eispack.o zobj/gro.o
VSOBJS = zobj/vsubtract.o zobj/elserf.o
VTOBJS = zobj/vtransf.o
#
all: zexe/ftrns.x zexe/hfk.x zexe/kbft.x zexe/separate.x zexe/ffiler.x \
     zexe/analyze.x zexe/vsubtract.x zexe/vtransf.x zexe/llscreen.x \
     zexe/backtran.x zexe/rtabber.x zexe/mason.x zexe/rtab2.x zexe/vamend.x zexe/hartpart.x

install:
	cp zexe/hfk.x $(INSTDIR)

clean:
	rm -f zobj/*.o zexe/*.x *__genmod.???

ALLOBJS = $(HFOBJS) $(HPOBJS) $(BTOBJS) $(FTOBJS) $(KBOBJS) $(LLOBJS) $(RTOBJS)$(RTOBJ2) $(SEOBJS) $(VSOBJS) $(VTOBJS)

ALLEXE =  zexe/ftrns.x zexe/hfk.x zexe/kbft.x zexe/separate.x zexe/ffiler.x \
     zexe/analyze.x zexe/vsubtract.x zexe/vtransf.x zexe/llscreen.x \
     zexe/backtran.x zexe/rtabber.x zexe/mason.x zexe/rtab2.x zexe/vamend.x zexe/hartpart.x


$(ALLOBJS): | zobj

zobj:
	@mkdir -p $@

$(ALLEXE): | zexe

zexe:
	@mkdir -p $@
#
zexe/hartpart.x: zobj/hartpart.o
	$(COMP) -o zexe/hartpart.x zobj/hartpart.o
zexe/mason.x: zobj/mason.o
	$(COMP) -o zexe/mason.x zobj/mason.o
zexe/analyze.x: zobj/analyze.o zobj/rinvert.o
	$(COMP) -o zexe/analyze.x zobj/analyze.o zobj/rinvert.o
zexe/backtran.x: $(BTOBJS)
	$(COMP) $(BTOBJS) -o zexe/backtran.x
zexe/ffiler.x: zobj/ffiler.o
	$(COMP) -o zexe/ffiler.x zobj/ffiler.o
zexe/ftrns.x: $(FTOBJS)
	$(COMP) $(FTOBJS) -o zexe/ftrns.x
zexe/hfk.x: $(HFOBJS)
	$(COMP) $(HFOBJS) -o zexe/hfk.x
zexe/kbft.x: $(KBOBJS)
	$(COMP) $(KBOBJS) -o zexe/kbft.x
zexe/llscreen.x: $(LLOBJS)
	$(COMP) $(LLOBJS) -o zexe/llscreen.x
zexe/rtabber.x: $(RTOBJS)
	$(COMP) $(RTOBJS) -o zexe/rtabber.x
zexe/rtab2.x: $(RTOBJ2)
	$(COMP) $(RTOBJ2) -o zexe/rtab2.x
zexe/separate.x: $(SEOBJS)
	$(COMP) $(SEOBJS) -o zexe/separate.x
zexe/vamend.x: zf90/vamend.f90
	$(COMP) -o zexe/vamend.x zf90/vamend.f90
zexe/vsubtract.x: $(VSOBJS)
	$(COMP) $(VSOBJS) -o zexe/vsubtract.x
zexe/vtransf.x: $(VTOBJS)
	$(COMP) $(VTOBJS) -o zexe/vtransf.x
#
zobj/amult.o: zf90/amult.f90
	$(COMP) $(OPTS) -c zf90/amult.f90 -o zobj/amult.o
zobj/besft.o: zf90/besft.f90
	$(COMP) $(OPTS) -c zf90/besft.f90 -o zobj/besft.o
zobj/besjl.o: zf90/besjl.f90
	$(COMP) $(OPTS) -c zf90/besjl.f90 -o zobj/besjl.o
zobj/bintegrate.o: zf90/bintegrate.f90
	$(COMP) $(OPTS) -c zf90/bintegrate.f90 -o zobj/bintegrate.o
zobj/der1.o: zf90/der1.f90
	$(COMP) $(OPTS) -c zf90/der1.f90 -o zobj/der1.o
zobj/der2.o: zf90/der2.f90
	$(COMP) $(OPTS) -c zf90/der2.f90 -o zobj/der2.o
zobj/dft.o: zf90/dft.f90
	$(COMP) $(OPTS) -c zf90/dft.f90 -o zobj/dft.o
zobj/diagq.o: zf90/diagq.f90
	$(COMP) $(OPTS) -c zf90/diagq.f90 -o zobj/diagq.o
zobj/dipole.o: zf90/dipole.f90
	$(COMP) $(OPTS) -c zf90/dipole.f90 -o zobj/dipole.o
zobj/egripper2.o: zf90/egripper2.f90
	$(COMP) $(OPTS) -c zf90/egripper2.f90 -o zobj/egripper2.o
zobj/baregrip.o: zf90/baregrip.f90
	$(COMP) $(OPTS) -c zf90/baregrip.f90 -o zobj/baregrip.o
zobj/getanu.o: zf90/getanu.f90
	$(COMP) $(OPTS) -c zf90/getanu.f90 -o zobj/getanu.o
zobj/getfg.o: zf90/getfg.f90
	$(COMP) $(OPTS) -c zf90/getfg.f90 -o zobj/getfg.o
zobj/getkap.o: zf90/getkap.f90
	$(COMP) $(OPTS) -c zf90/getkap.f90 -o zobj/getkap.o
zobj/kbbessel.o: zf90/kbbessel.f90
	$(COMP) $(OPTS) -c zf90/kbbessel.f90 -o zobj/kbbessel.o
zobj/kbft.o: zf90/kbft.f90
	$(COMP) $(OPTS) -c zf90/kbft.f90 -o zobj/kbft.o
zobj/longit.o: zf90/longit.f90
	$(COMP) $(OPTS) -c zf90/longit.f90 -o zobj/longit.o
zobj/melwriter.o: zf90/melwriter.f90
	$(COMP) $(OPTS) -c zf90/melwriter.f90 -o zobj/melwriter.o
#
# pseudo stuff ...
zobj/pseudo.o: zf90/pseudo.f90
	$(COMP) $(OPTS) -c zf90/pseudo.f90 -o zobj/pseudo.o
zobj/pseudize.o: zf90/pseudize.f90
	$(COMP) $(OPTS) -c zf90/pseudize.f90 -o zobj/pseudize.o
zobj/ffitx0.o: zf90/ffitx0.f90
	$(COMP) $(OPTS) -c zf90/ffitx0.f90 -o zobj/ffitx0.o
zobj/unscreen.o: zf90/unscreen.f90
	$(COMP) $(OPTS) -c zf90/unscreen.f90 -o zobj/unscreen.o
zobj/mklfc.o: zf90/mklfc.f90
	$(COMP) $(OPTS) -c zf90/mklfc.f90 -o zobj/mklfc.o
#
#
zobj/recon.o: zf90/recon.f90
	$(COMP) $(OPTS) -c zf90/recon.f90 -o zobj/recon.o
zobj/sdump.o: zf90/sdump.f90
	$(COMP) $(OPTS) -c zf90/sdump.f90 -o zobj/sdump.o
zobj/ssubsh.o: zf90/ssubsh.f90
	$(COMP) $(OPTS) -c zf90/ssubsh.f90 -o zobj/ssubsh.o
zobj/vtransf.o: zf90/vtransf.f90
	$(COMP) $(OPTS) -c zf90/vtransf.f90 -o zobj/vtransf.o
zobj/xidump.o: zf90/xidump.f90
	$(COMP) $(OPTS) -c zf90/xidump.f90 -o zobj/xidump.o
zobj/hartpart.o: zf90/hartpart.f90
	$(COMP) $(OPTS) -c zf90/hartpart.f90 -o zobj/hartpart.o
#
zobj/bach.o: zstd/bach.f
	$(COMP) $(OPTS) -c zstd/bach.f -o zobj/bach.o
zobj/eispack.o: zstd/eispack.f
	$(COMP) $(OPTS) -c zstd/eispack.f -o zobj/eispack.o
zobj/elsdch.o: zstd/elsdch.f
	$(COMP) $(OPTS) -c zstd/elsdch.f -o zobj/elsdch.o
zobj/elserf.o: zstd/elserf.f
	$(COMP) $(OPTS) -c zstd/elserf.f -o zobj/elserf.o
zobj/rinvert.o: zstd/rinvert.f
	$(COMP) $(OPTS) -c zstd/rinvert.f -o zobj/rinvert.o
zobj/spherbes.o: zstd/spherbes.f
	$(COMP) $(OPTS) -c zstd/spherbes.f -o zobj/spherbes.o
#
zobj/backtran.o: zf77/backtran.f
	$(COMP) $(OPTS) -c zf77/backtran.f -o zobj/backtran.o
zobj/hfca.o: zf77/hfca.f ./alfinv.h
	$(COMP) $(OPTS) -c zf77/hfca.f -o zobj/hfca.o
zobj/hfgi.o: zf77/hfgi.f
	$(COMP) $(OPTS) -c zf77/hfgi.f -o zobj/hfgi.o
zobj/hfhead.o: zf90/hfhead.f90
	$(COMP) $(OPTS) -c zf90/hfhead.f90 -o zobj/hfhead.o
zobj/hfhead2.o: zf77/hfhead2.f
	$(COMP) $(OPTS) -c zf77/hfhead2.f -o zobj/hfhead2.o
zobj/hfneck.o: zf77/hfneck.f ./alfinv.h
	$(COMP) $(OPTS) -c zf77/hfneck.f -o zobj/hfneck.o
zobj/hfpa.o: zf77/hfpa.f ./alfinv.h
	$(COMP) $(OPTS) -c zf77/hfpa.f -o zobj/hfpa.o
zobj/hfrg.o: zf77/hfrg.f
	$(COMP) $(OPTS) -c zf77/hfrg.f -o zobj/hfrg.o
zobj/hfvwn.o: zf77/hfvwn.f
	$(COMP) $(OPTS) -c zf77/hfvwn.f -o zobj/hfvwn.o
zobj/mkkbfile.o: zf77/mkkbfile.f
	$(COMP) $(OPTS) -c zf77/mkkbfile.f -o zobj/mkkbfile.o
zobj/mkkxfile.o: zf77/mkkxfile.f
	$(COMP) $(OPTS) -c zf77/mkkxfile.f -o zobj/mkkxfile.o
zobj/mkvaltab.o: zf77/mkvaltab.f
	$(COMP) $(OPTS) -c zf77/mkvaltab.f -o zobj/mkvaltab.o
zobj/vsubtract.o: zf77/vsubtract.f
	$(COMP) $(OPTS) -c zf77/vsubtract.f -o zobj/vsubtract.o
#
zobj/alloy.o: zf77/alloy.f
	$(COMP) $(OPTS) -c zf77/alloy.f -o zobj/alloy.o
zobj/analyze.o: zf77/analyze.f
	$(COMP) $(OPTS) -c zf77/analyze.f -o zobj/analyze.o
zobj/sean1.o: zf77/sean1.f
	$(COMP) $(OPTS) -c zf77/sean1.f -o zobj/sean1.o
zobj/sean2.o: zf77/sean2.f
	$(COMP) $(OPTS) -c zf77/sean2.f -o zobj/sean2.o
zobj/sean4.o: zf77/sean4.f
	$(COMP) $(OPTS) -c zf77/sean4.f -o zobj/sean4.o
zobj/sean5.o: zf77/sean5.f
	$(COMP) $(OPTS) -c zf77/sean5.f -o zobj/sean5.o
zobj/getpot.o: zf77/getpot.f ./alfinv.h
	$(COMP) $(OPTS) -c zf77/getpot.f -o zobj/getpot.o
zobj/sig.o: zf77/sig.f
	$(COMP) $(OPTS) -c zf77/sig.f -o zobj/sig.o
zobj/gro.o: zf77/gro.f
	$(COMP) $(OPTS) -c zf77/gro.f -o zobj/gro.o
zobj/mason.o: zf77/mason.f
	$(COMP) $(OPTS) -c zf77/mason.f -o zobj/mason.o
zobj/rtabber.o: zf77/rtabber.f
	$(COMP) $(OPTS) -c zf77/rtabber.f -o zobj/rtabber.o
zobj/rtab2.o: zf77/rtab2.f
	$(COMP) $(OPTS) -c zf77/rtab2.f -o zobj/rtab2.o
zobj/elsft.o: zf77/elsft.f
	$(COMP) $(OPTS) -c zf77/elsft.f -o zobj/elsft.o
zobj/levlou.o: zf77/levlou.f
	$(COMP) $(OPTS) -c zf77/levlou.f -o zobj/levlou.o
zobj/llscreen.o: zf77/llscreen.f
	$(COMP) $(OPTS) -c zf77/llscreen.f -o zobj/llscreen.o
zobj/vwntst.o: zf77/vwntst.f
	$(COMP) $(OPTS) -c zf77/vwntst.f -o zobj/vwntst.o
zobj/grandop.o: zf77/grandop.f
	$(COMP) $(OPTS) -c zf77/grandop.f -o zobj/grandop.o
zobj/corepot.o: zf77/corepot.f
	$(COMP) $(OPTS) -c zf77/corepot.f -o zobj/corepot.o
zobj/getexch.o: zf90/getexch.f90
	$(COMP) $(OPTS) -c zf90/getexch.f90 -o zobj/getexch.o
#
zobj/config.o: zf90/config.f90
	$(COMP) $(OPTS) -c zf90/config.f90 -o zobj/config.o
zobj/transp.o: zf90/transp.f90
	$(COMP) $(OPTS) -c zf90/transp.f90 -o zobj/transp.o
zobj/leadbehv.o: zf77/leadbehv.f
	$(COMP) $(OPTS) -c zf77/leadbehv.f -o zobj/leadbehv.o
zobj/elener.o: zf90/elener.f90
	$(COMP) $(OPTS) -c zf90/elener.f90 -o zobj/elener.o
zobj/ffiler.o: zf77/ffiler.f
	$(COMP) $(OPTS) -c zf77/ffiler.f -o zobj/ffiler.o
zobj/weighter.o: zf77/weighter.f
	$(COMP) $(OPTS) -c zf77/weighter.f -o zobj/weighter.o
zobj/escanner.o: zf90/escanner.f90
	$(COMP) $(OPTS) -c zf90/escanner.f90 -o zobj/escanner.o
#
zobj/hfpaold.o: zf77/zoct20/hfpa.f
	$(COMP) $(OPTS) -c zf77/zoct20/hfpa.f -o zobj/hfpaold.o
zobj/integ.o: zf90/integ.f90
	$(COMP) $(OPTS) -c zf90/integ.f90 -o zobj/integ.o
zobj/integi.o: zf90/integi.f90
	$(COMP) $(OPTS) -c zf90/integi.f90 -o zobj/integi.o
zobj/intego.o: zf90/intego.f90
	$(COMP) $(OPTS) -c zf90/intego.f90 -o zobj/intego.o
zobj/setxkdk.o: zf90/setxkdk.f90 ./alfinv.h
	$(COMP) $(OPTS) -c zf90/setxkdk.f90 -o zobj/setxkdk.o
#
zobj/ppload.o: zf90/ppload.f90
	$(COMP) $(OPTS) -c zf90/ppload.f90 -o zobj/ppload.o 
zobj/ppdump.o: zf90/ppdump.f90
	$(COMP) $(OPTS) -c zf90/ppdump.f90 -o zobj/ppdump.o
zobj/mkcorcon.o: zf90/mkcorcon.f90
	$(COMP) $(OPTS) -c zf90/mkcorcon.f90 -o zobj/mkcorcon.o
zobj/intval.o: zstd/intval.f90
	$(COMP) $(OPTS) -c zstd/intval.f90 -o zobj/intval.o
#
zobj/freshen.o: zf90/freshen.f90
	$(COMP) $(OPTS) -c zf90/freshen.f90 -o zobj/freshen.o
zobj/pwavesetter.o: zf90/pwavesetter.f90
	$(COMP) $(OPTS) -c zf90/pwavesetter.f90 -o zobj/pwavesetter.o
zobj/optradf.o: zf90/optradf.f90
	$(COMP) $(OPTS) -c zf90/optradf.f90 -o zobj/optradf.o
zobj/fipfront.o: zf90/fipfront.f90
	$(COMP) $(OPTS) -c zf90/fipfront.f90 -o zobj/fipfront.o
zobj/formfkgknll.o: zf90/formfkgknll.f90
	$(COMP) $(OPTS) -c zf90/formfkgknll.f90 -o zobj/formfkgknll.o
zobj/formmelnll.o: zf90/formmelnll.f90
	$(COMP) $(OPTS) -c zf90/formmelnll.f90 -o zobj/formmelnll.o
zobj/pawload.o: zf90/pawload.f90
	$(COMP) $(OPTS) -c zf90/pawload.f90 -o zobj/pawload.o
zobj/coreload.o: zf90/coreload.f90
	$(COMP) $(OPTS) -c zf90/coreload.f90 -o zobj/coreload.o
zobj/coullgen.o: zf90/coullgen.f90
	$(COMP) $(OPTS) -c zf90/coullgen.f90 -o zobj/coullgen.o
zobj/powlgen.o: zf90/powlgen.f90
	$(COMP) $(OPTS) -c zf90/powlgen.f90 -o zobj/powlgen.o
zobj/normangnodes.o: zf90/normangnodes.f90
	$(COMP) $(OPTS) -c zf90/normangnodes.f90 -o zobj/normangnodes.o
zobj/reconstruct.o: zf90/reconstruct.f90
	$(COMP) $(OPTS) -c zf90/reconstruct.f90 -o zobj/reconstruct.o
zobj/newdumper.o: zf90/newdumper.f90
	$(COMP) $(OPTS) -c zf90/newdumper.f90 -o zobj/newdumper.o
zobj/projdumper.o: zf90/projdumper.f90
	$(COMP) $(OPTS) -c zf90/projdumper.f90 -o zobj/projdumper.o
zobj/escanprep.o: zf90/escanprep.f90
	$(COMP) $(OPTS) -c zf90/escanprep.f90 -o zobj/escanprep.o
#
zobj/fillinpaw.o: zf90/fillinpaw.f90
	$(COMP) $(OPTS) -c zf90/fillinpaw.f90 -o zobj/fillinpaw.o
zobj/spartanfip.o: zf90/spartanfip.f90
	$(COMP) $(OPTS) -c zf90/spartanfip.f90 -o zobj/spartanfip.o
#
zobj/continuum.o: zf90/continuum.f90
	$(COMP) $(OPTS) -c zf90/continuum.f90 -o zobj/continuum.o
zobj/projaug.o: zf90/projaug.f90
	$(COMP) $(OPTS) -c zf90/projaug.f90 -o zobj/projaug.o
zobj/trck.o: zf90/trck.f90
	$(COMP) $(OPTS) -c zf90/trck.f90 -o zobj/trck.o
#
zobj/rpower.o: zf90/rpower.f90
	$(COMP) $(OPTS) -c zf90/rpower.f90 -o zobj/rpower.o
zobj/hpload.o: zf90/hpload.f90
	$(COMP) $(OPTS) -c zf90/hpload.f90 -o zobj/hpload.o
zobj/chgocc.o: zf90/chgocc.f90
	$(COMP) $(OPTS) -c zf90/chgocc.f90 -o zobj/chgocc.o
zobj/optript.o: zf90/optript.f90
	$(COMP) $(OPTS) -c zf90/optript.f90 -o zobj/optript.o
zobj/orthred.o: zf90/orthred.f90
	$(COMP) $(OPTS) -c zf90/orthred.f90 -o zobj/orthred.o
zobj/potfigure.o: zf90/potfigure.f90
	$(COMP) $(OPTS) -c zf90/potfigure.f90 -o zobj/potfigure.o
zobj/radint.o: zf90/radint.f90
	$(COMP) $(OPTS) -c zf90/radint.f90 -o zobj/radint.o
zobj/geticore.o: zf90/geticore.f90
	$(COMP) $(OPTS) -c zf90/geticore.f90 -o zobj/geticore.o
zobj/getfgnew.o: zf90/getfgnew.f90
	$(COMP) $(OPTS) -c zf90/getfgnew.f90 -o zobj/getfgnew.o
zobj/fgcalc.o: zf90/fgcalc.f90
	$(COMP) $(OPTS) -c zf90/fgcalc.f90 -o zobj/fgcalc.o
zobj/getmeznl.o: zf90/getmeznl.f90
	$(COMP) $(OPTS) -c zf90/getmeznl.f90 -o zobj/getmeznl.o
zobj/orbmom.o: zf90/orbmom.f90
	$(COMP) $(OPTS) -c zf90/orbmom.f90 -o zobj/orbmom.o
zobj/screencore.o: zf90/screencore.f90
	$(COMP) $(OPTS) -c zf90/screencore.f90 -o zobj/screencore.o
zobj/orbcont.o: zf90/orbcont.f90
	$(COMP) $(OPTS) -c zf90/orbcont.f90 -o zobj/orbcont.o
zobj/addpot.o: zf90/addpot.f90
	$(COMP) $(OPTS) -c zf90/addpot.f90 -o zobj/addpot.o
zobj/setqmm.o: zf90/setqmm.f90 ./alfinv.h
	$(COMP) $(OPTS) -c zf90/setqmm.f90 -o zobj/setqmm.o
zobj/getplead.o: zf90/getplead.f90 ./alfinv.h
	$(COMP) $(OPTS) -c zf90/getplead.f90 -o zobj/getplead.o
#
zobj/elsolve.o: zf90/elsolve.f90
	$(COMP) $(OPTS) -c zf90/elsolve.f90 -o zobj/elsolve.o
zobj/bcsolve.o: zf90/bcsolve.f90
	$(COMP) $(OPTS) -c zf90/bcsolve.f90 -o zobj/bcsolve.o
zobj/radnorm.o: zf90/radnorm.f90
	$(COMP) $(OPTS) -c zf90/radnorm.f90 -o zobj/radnorm.o
zobj/augment.o: zf90/augment.f90
	$(COMP) $(OPTS) -c zf90/augment.f90 -o zobj/augment.o
