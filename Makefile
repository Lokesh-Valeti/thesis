all: prac

CXXFLAGS=-march=native -std=c++17 -Wall -Wno-ignored-attributes -ggdb -O3
LDFLAGS=-ggdb
LDLIBS=-lbsd -lboost_system -lboost_context -lboost_chrono -lboost_thread -lpthread

# Enable this to have all communication logged to stdout
# CXXFLAGS += -DVERBOSE_COMMS

BIN=prac
SRCS=prac.cpp mpcio.cpp preproc.cpp online.cpp mpcops.cpp rdpf.cpp \
    cdpf.cpp duoram.cpp cell.cpp
OBJS=$(SRCS:.cpp=.o)
ASMS=$(SRCS:.cpp=.s)

$(BIN): $(OBJS)
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.s: %.cpp
	g++ $(CXXFLAGS) -S -o $@ $^

# Remove the files created by the preprocessing phase
reset:
	-rm -f *.p[012].t*

clean: reset
	-rm -f $(BIN) $(OBJS) $(ASMS)

depend:
	makedepend -Y -- $(CXXFLAGS) -- $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

prac.o: mpcio.hpp types.hpp corotypes.hpp mpcio.tcc preproc.hpp options.hpp
prac.o: online.hpp
mpcio.o: mpcio.hpp types.hpp corotypes.hpp mpcio.tcc rdpf.hpp coroutine.hpp
mpcio.o: bitutils.hpp dpf.hpp prg.hpp aes.hpp rdpf.tcc mpcops.hpp cdpf.hpp
mpcio.o: cdpf.tcc
preproc.o: types.hpp coroutine.hpp corotypes.hpp mpcio.hpp mpcio.tcc
preproc.o: preproc.hpp options.hpp rdpf.hpp bitutils.hpp dpf.hpp prg.hpp
preproc.o: aes.hpp rdpf.tcc mpcops.hpp cdpf.hpp cdpf.tcc
online.o: online.hpp mpcio.hpp types.hpp corotypes.hpp mpcio.tcc options.hpp
online.o: mpcops.hpp coroutine.hpp rdpf.hpp bitutils.hpp dpf.hpp prg.hpp
online.o: aes.hpp rdpf.tcc duoram.hpp duoram.tcc cdpf.hpp cdpf.tcc cell.hpp
mpcops.o: mpcops.hpp types.hpp mpcio.hpp corotypes.hpp mpcio.tcc
mpcops.o: coroutine.hpp bitutils.hpp
rdpf.o: rdpf.hpp mpcio.hpp types.hpp corotypes.hpp mpcio.tcc coroutine.hpp
rdpf.o: bitutils.hpp dpf.hpp prg.hpp aes.hpp rdpf.tcc mpcops.hpp
cdpf.o: bitutils.hpp cdpf.hpp mpcio.hpp types.hpp corotypes.hpp mpcio.tcc
cdpf.o: coroutine.hpp dpf.hpp prg.hpp aes.hpp cdpf.tcc
duoram.o: duoram.hpp types.hpp mpcio.hpp corotypes.hpp mpcio.tcc
duoram.o: coroutine.hpp duoram.tcc mpcops.hpp cdpf.hpp dpf.hpp prg.hpp
duoram.o: bitutils.hpp aes.hpp cdpf.tcc rdpf.hpp rdpf.tcc
cell.o: types.hpp duoram.hpp mpcio.hpp corotypes.hpp mpcio.tcc coroutine.hpp
cell.o: duoram.tcc mpcops.hpp cdpf.hpp dpf.hpp prg.hpp bitutils.hpp aes.hpp
cell.o: cdpf.tcc rdpf.hpp rdpf.tcc cell.hpp options.hpp
