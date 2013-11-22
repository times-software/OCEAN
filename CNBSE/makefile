all:
	cd zcain; make "COMP = $(COMP)" "OPTS = $(OPTS)"
	cd zcks;  make "COMP = $(COMP)" "OPTS = $(OPTS)"
	cd zcmsnmult; make "COMP = $(COMP)" "OPTS = $(OPTS)"
	cd zf90; make "COMP = $(COMP)" "OPTS = $(OPTS)"

install:
	cd zcain; make install "INSTDIR = $(INSTDIR)"
	cd zcks; make install "INSTDIR = $(INSTDIR)" 
	cd zcmsnmult; make install "INSTDIR = $(INSTDIR)"
	cd zf90; make install "INSTDIR = $(INSTDIR)"
	cp cnbse.pl OBF_cnbse.pl $(INSTDIR)/

clean:
	cd zcain; make clean
	cd zcks; make clean
	cd zcmsnmult; make clean
	cd zf90; make clean

