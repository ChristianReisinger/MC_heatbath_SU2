CXX=g++
CXXFLAGS=-O3
LIBS = -lm

o_files: fields.o io.o ranlux.o ranlxd.o ranlxs.o Wilson_loops.o #smearing_techniques.o <- still SU(3) !!! 

fields.o: fields.cc fields.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

io.o: io.cc io.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

ranlux.o: ranlux.cc ranlux.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

ranlxd.o: ranlxd.h ranlxd.c
	${CXX} ${CXXFLAGS} -c -o ranlxd.o ranlxd.c

ranlxs.o: ranlxs.h ranlxs.c
	${CXX} ${CXXFLAGS} -c -o ranlxs.o ranlxs.c

smearing_techniques.o: smearing_techniques.cc smearing_techniques.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

Wilson_loops.o: Wilson_loops.cc Wilson_loops.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

clean:
	rm -f *~ *.o
