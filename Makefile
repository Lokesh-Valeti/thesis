all: oblivds

CXXFLAGS=-std=c++17 -Wall -ggdb
LDFLAGS=-ggdb
LDLIBS=-lbsd -lboost_system -lboost_context -lboost_thread -lpthread

BIN=oblivds
OBJS=oblivds.o mpcio.o preproc.o online.o mpcops.o

$(BIN): $(OBJS)
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)

oblivds.o: preproc.hpp mpcio.hpp types.hpp
mpcio.o: mpcio.hpp types.hpp
preproc.o: preproc.hpp mpcio.hpp types.hpp
online.o: online.hpp mpcops.hpp coroutine.hpp
mpcops.o: mpcops.hpp coroutine.hpp

# Remove the files created by the preprocessing phase
reset:
	-rm -f *.p[01].t*

clean: reset
	-rm -f $(BIN) $(OBJS)
