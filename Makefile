SHELL=/bin/bash

CXX=g++
CXXFLAGS=-std=c++11 -O3

BIN_DIR=bin/
BUILD_DIR=build/
SRC_DIR=src/
#INCLUDE_DIR=

#OBJ_NAMES=
BIN_NAMES=MC_heatbath

#OBJS=${OBJ_NAMES:%=${BUILD_DIR}/%.o}
TARGETS=${BIN_NAMES:%=${BIN_DIR}/%}
LIBS=-lm -l:ranlxd.o -l:ranlxs.o

all: ${TARGETS}

${TARGETS}: ${BIN_DIR}/%: ${SRC_DIR}/%.cc #${OBJS} 
	readarray deps < <(depfinder.sh . $<);\
	${CXX} ${CXXFLAGS} -o $@ $^ "$${deps[@]}" ${LIBS}
	
#${OBJS}: ${BUILD_DIR}/%.o: ${SRC_DIR}/%.cc  ${INCLUDE_DIR}/%.hh
#	${CXX} ${CXXFLAGS} -c -o $@ $< -I"${INCLUDE_DIR}" $$(depfinder.sh . $<)

clean:
	rm -f ${BIN_DIR}/* ${BUILD_DIR}/*
