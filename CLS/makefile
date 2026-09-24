include ../Makefile.arch

COMP = $(FC)
OPTS = $(OPTIONS)

all:
	cd src; make 

install:
	cp cls.pl $(INSTDIR)/ 
	cd src; make install

clean:
	cd src; make clean
