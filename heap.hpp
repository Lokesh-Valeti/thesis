#ifndef __HEAP_HPP__
#define __HEAP_HPP__

#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"
 
#include "options.hpp"
#include "mpcops.hpp" 

class HEAP {
  private: 
    Duoram<RegAS> *oram;
    // RegXS root;

     
     size_t MAX_SIZE;

    // std::tuple<RegXS, RegBS> insert(MPCTIO &tio, yield_t &yield, RegXS ptr,
    //     const Node &new_node, Duoram<Node>::Flat &A, int TTL, RegBS isDummy);
    // void insert(MPCTIO &tio, yield_t &yield, const Node &node, Duoram<Node>::Flat &A);


    // int del(MPCTIO &tio, yield_t &yield, RegXS ptr, RegAS del_key,
    //     Duoram<Node>::Flat &A, RegBS F_af, RegBS F_fs, int TTL, 
    //     del_return &ret_struct);

  public:
    size_t num_items = 0;
    HEAP(int num_players, size_t size) {
      this->initialize(num_players, size);
    };

    ~HEAP() {
      if(oram)
        delete oram;
    };

    void  initialize(int num_players, size_t size);
    RegAS extract_min(MPCTIO tio, yield_t &yield);
    int   insert(MPCTIO tio, yield_t &yield, RegAS val);
    int   verify_heap_property(MPCTIO tio, yield_t &yield);
    RegAS restore_heap_property(MPCTIO tio, yield_t &yield, RegAS index);
    RegAS restore_heap_property_at_root(MPCTIO tio, yield_t &yield, size_t index);
};

 void Heap(MPCIO &mpcio,
    const PRACOptions &opts, char **args);
 
 //#include "heap.tcc"
#endif
