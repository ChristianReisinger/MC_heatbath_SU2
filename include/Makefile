CXX=g++
CXXFLAGS=-O3 -D __SUN_N__=2
LIBS = -lm
INCLUDE = -I"./"
INCLUDE_EXT := ${shell sed -e 's/^/-I"/' -e 's/$$/"/' includes | tr '\n' ' '}


o_files: fields.o io.o ranlux.o ranlxd.o ranlxs.o Wilson_loops.o smearing_techniques.o heatbath.o

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

heatbath.o: heatbath.cc heatbath.hh
	${CXX} ${CXXFLAGS} -c -o $@ $< ${INCLUDE_EXT} ${INCLUDE}

smearing_techniques.o: smearing_techniques.cc smearing_techniques.hh
	${CXX} ${CXXFLAGS} -c -o $@ $< ${INCLUDE_EXT} ${INCLUDE}

Wilson_loops.o: Wilson_loops.cc Wilson_loops.hh
	${CXX} ${CXXFLAGS} -c -o $@ $<

clean:
	rm -f *~ *.o
