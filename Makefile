CXX=g++
CXXFLAGS=-O3
LIBS=-lm

all: MC_heatbath

MC_heatbath: MC_heatbath.o
	${CXX} ${CXXFLAGS} -o $@ MC_heatbath.o includes/fields.o includes/io.o includes/ranlux.o includes/ranlxd.o includes/ranlxs.o includes/Wilson_loops.o -I"includes/" ${LIBS}

MC_heatbath.o: MC_heatbath.cc
	${CXX} ${CXXFLAGS} -c -o $@ $< -I"includes/"

clean:
	rm -f *~ *.o MC_heatbath
