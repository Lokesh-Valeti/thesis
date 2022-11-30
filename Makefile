all: oblivds

CXXFLAGS=-std=c++17 -Wall -ggdb
LDFLAGS=-ggdb
LDLIBS=-lbsd -lboost_system -lboost_context -lboost_thread -lpthread

oblivds: oblivds.o mpcio.o preproc.o
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)

oblivds.o: preproc.hpp mpcio.hpp
mpcio.o: mpcio.hpp
preproc.o: preproc.hpp mpcio.hpp
