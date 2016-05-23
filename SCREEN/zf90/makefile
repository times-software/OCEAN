COMP=pgf90
OPTS=-O3

all: zexe/avg.x zexe/builder.x zexe/pawsetup.x zexe/mkrbfile.x zexe/mkrbfile_mult.x zexe/mpi_avg.x

BUILDER = zobj/builder.o zobj/realu2.o zobj/vcbder2.o zobj/getomega.o zobj/mkvipt.o zobj/rdp.o \
          zobj/intval.o

AVG = zobj/avg.o zobj/snatch.o zobj/getomega.o

MPI_AVG = zobj/mpi_avg.o zobj/snatch.o zobj/getomega.o

zexe/mkrbfile_mult.x: zobj/mkrbfile_mult.o zobj/mkcmesh.o zobj/snatch.o zobj/mkbmesh.o
	$(COMP) $(OPTS) zobj/mkrbfile_mult.o zobj/mkcmesh.o zobj/snatch.o zobj/mkbmesh.o -o zexe/mkrbfile_mult.x

zexe/mkrbfile.x: zobj/mkrbfile.o zobj/mkcmesh.o zobj/snatch.o
	$(COMP) $(OPTS) zobj/mkrbfile.o zobj/mkcmesh.o zobj/snatch.o -o zexe/mkrbfile.x

zexe/pawsetup.x: pawsetup.f90 zobj/periodic.o
	$(COMP) $(OPTS) pawsetup.f90 zobj/periodic.o -o zexe/pawsetup.x

zexe/avg.x: $(AVG)
	$(COMP) $(OPTS) $(AVG) -o zexe/avg.x

zexe/mpi_avg.x: $(MPI_AVG)
	$(MPIFORT) $(OPTS) $(MPI_AVG) -o zexe/mpi_avg.x

zexe/builder.x: $(BUILDER)
	$(COMP) $(OPTS) $(BUILDER) -o zexe/builder.x

zobj/avg.o: avg.f90
	$(COMP) $(OPTS) -c avg.f90 -o zobj/avg.o

zobj/mpi_avg.o: MPI_avg.f90
	$(MPIFORT) $(OPTS) -c MPI_avg.f90 -o zobj/mpi_avg.o

zobj/builder.o: builder.f90 
	$(COMP) $(OPTS) -c builder.f90 -o zobj/builder.o

zobj/getomega.o: getomega.f90
	$(COMP) $(OPTS) -c getomega.f90 -o zobj/getomega.o

zobj/intval.o: intval.f90
	$(COMP) $(OPTS) -c intval.f90 -o zobj/intval.o

zobj/mkcmesh.o: mkcmesh.f90
	$(COMP) $(OPTS) -c mkcmesh.f90 -o zobj/mkcmesh.o

zobj/mkbmesh.o: mkbmesh.f90
	$(COMP) $(OPTS) -c mkbmesh.f90 -o zobj/mkbmesh.o

zobj/mkrbfile.o: mkrbfile.f90
	$(COMP) $(OPTS) -c mkrbfile.f90 -o zobj/mkrbfile.o

zobj/mkrbfile_mult.o: mkrbfile_mult.f90
	$(COMP) $(OPTS) -c mkrbfile_mult.f90 -o zobj/mkrbfile_mult.o

zobj/mkvipt.o: mkvipt.f90
	$(COMP) $(OPTS) -c mkvipt.f90 -o zobj/mkvipt.o

zobj/pawsetup.o: pawsetup.f90 periodic.mod
	$(COMP) $(OPTS) -c pawsetup.f90 -o zobj/pawsetup.o

periodic.mod: zobj/periodic.o
zobj/periodic.o: periodic.f90
	$(COMP) $(OPTS) -c periodic.f90 -o zobj/periodic.o

zobj/rdp.o: rdp.f90
	$(COMP) $(OPTS) -c rdp.f90 -o zobj/rdp.o

zobj/realu2.o: realu2.f90
	$(COMP) $(OPTS) -c realu2.f90 -o zobj/realu2.o

zobj/realu.o: realu.f90
	$(COMP) $(OPTS) -c realu.f90 -o zobj/realu.o

zobj/snatch.o: snatch.f90
	$(COMP) $(OPTS) -c snatch.f90 -o zobj/snatch.o

zobj/vcbder2.o: vcbder2.f90
	$(COMP) $(OPTS) -c vcbder2.f90 -o zobj/vcbder2.o

zobj/vcbder.o: vcbder.f90
	$(COMP) $(OPTS) -c vcbder.f90 -o zobj/vcbder.o

zobj/xxxhed.o: xxxhed.f90
	$(COMP) $(OPTS) -c xxxhed.f90 -o zobj/xxxhed.o

install:
	cp zexe/avg.x $(INSTDIR)/
	cp zexe/mpi_avg.x $(INSTDIR)/
	cp zexe/builder.x $(INSTDIR)/
	cp zexe/pawsetup.x $(INSTDIR)/
	cp zexe/mkrbfile.x $(INSTDIR)/
	cp zexe/mkrbfile_mult.x $(INSTDIR)/

clean:
	rm -f zobj/*.o
	rm -f zexe/*.x
