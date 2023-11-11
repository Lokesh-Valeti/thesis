#ifndef __HEAP_HPP__
#define __HEAP_HPP__

#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"
#include "options.hpp"
#include "mpcops.hpp"

class MinHeap {
    private: 
    Duoram < RegAS > oram;
    size_t MAX_SIZE;
    size_t num_items;

    public: 
    
    MinHeap(int player_num, size_t size) : oram(player_num, size){

    };

    void set_num_items (size_t n) {num_items = n;}
    RegAS extract_min(MPCIO &mpcio, MPCTIO tio, yield_t & yield, int is_optimized);
    void init(MPCTIO tio, yield_t & yield);
    void init(MPCTIO tio, yield_t & yield, size_t which_init);
    int insert(MPCTIO tio, yield_t & yield, RegAS val);
    void insert_optimized(MPCTIO tio, yield_t & yield, RegAS val);
    
    #ifdef HEAP_DEBUG
    void verify_heap_property(MPCTIO tio, yield_t & yield);
    #endif
    
    RegXS restore_heap_property(MPCIO &mpcio, MPCTIO tio, yield_t & yield, RegXS index);
    std::pair<RegXS, RegBS> restore_heap_property_optimized(MPCTIO tio, yield_t & yield, RegXS index, size_t depth, size_t layer, typename Duoram<RegAS>::template OblivIndex<RegXS,3> oidx);
    std::pair<RegXS, RegBS> restore_heap_property_at_explicit_index(MPCTIO tio, yield_t & yield,  size_t index);
    void print_heap(MPCTIO tio, yield_t & yield);
};

void Heap(MPCIO &mpcio,
   const PRACOptions &opts, char **args, int argc);
#endif