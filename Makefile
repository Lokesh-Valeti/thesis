all: oblivds

CXXFLAGS=-std=c++17 -Wall -ggdb
LDFLAGS=-ggdb
LDLIBS=-lbsd -lboost_system -lboost_context -lboost_thread -lpthread

BIN=oblivds
OBJS=oblivds.o mpcio.o preproc.o

$(BIN): $(OBJS)
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)

oblivds.o: preproc.hpp mpcio.hpp
mpcio.o: mpcio.hpp
preproc.o: preproc.hpp mpcio.hpp

clean:
	-rm -f $(BIN) $(OBJS) *.p[01].t*
