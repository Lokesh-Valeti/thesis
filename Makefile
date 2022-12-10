all: oblivds

CXXFLAGS=-std=c++17 -Wall -ggdb
LDFLAGS=-ggdb
LDLIBS=-lbsd -lboost_system -lboost_context -lboost_chrono -lboost_thread -lpthread

BIN=oblivds
SRCS=oblivds.cpp mpcio.cpp preproc.cpp online.cpp mpcops.cpp
OBJS=$(SRCS:.cpp=.o)

$(BIN): $(OBJS)
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Remove the files created by the preprocessing phase
reset:
	-rm -f *.p[01].t*

clean: reset
	-rm -f $(BIN) $(OBJS)

depend:
	makedepend -Y -- $(CXXFLAGS) -- $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

oblivds.o: mpcio.hpp types.hpp preproc.hpp online.hpp
mpcio.o: mpcio.hpp types.hpp
preproc.o: types.hpp preproc.hpp mpcio.hpp
online.o: online.hpp mpcio.hpp types.hpp mpcops.hpp coroutine.hpp
mpcops.o: mpcops.hpp types.hpp mpcio.hpp coroutine.hpp
