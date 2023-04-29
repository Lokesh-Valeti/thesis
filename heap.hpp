#ifndef __HEAP_HPP__
#define __HEAP_HPP__

#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"
#include "options.hpp"
#include "mpcops.hpp"

class MinHeap {
    private: Duoram < RegAS > oram;

    size_t MAX_SIZE;

    public: size_t num_items = 0;
    MinHeap(int num_players, size_t size) : oram(num_players, size), MAX_SIZE(size){

    };

    
    RegAS extract_min(MPCTIO tio, yield_t & yield, int is_optimized);
    void initialize(MPCTIO tio, yield_t & yield);
    int insert(MPCTIO tio, yield_t & yield, RegAS val);
    int insert_optimized(MPCTIO tio, yield_t & yield, RegAS val);
    int verify_heap_property(MPCTIO tio, yield_t & yield);
    RegXS restore_heap_property(MPCTIO tio, yield_t & yield, RegXS index);
    auto restore_heap_property_optimized(MPCTIO tio, yield_t & yield, RegXS index, size_t depth, size_t layer, typename Duoram<RegAS>::template OblivIndex<RegXS,3> oidx);
    auto restore_heap_property_at_root(MPCTIO tio, yield_t & yield);
    void print_heap(MPCTIO tio, yield_t & yield);
};

void Heap(MPCIO &mpcio,
   const PRACOptions &opts, char **args);
#endif