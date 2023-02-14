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

     size_t num_items = 0;
     size_t MAX_SIZE;

    // std::tuple<RegXS, RegBS> insert(MPCTIO &tio, yield_t &yield, RegXS ptr,
    //     const Node &new_node, Duoram<Node>::Flat &A, int TTL, RegBS isDummy);
    // void insert(MPCTIO &tio, yield_t &yield, const Node &node, Duoram<Node>::Flat &A);


    // int del(MPCTIO &tio, yield_t &yield, RegXS ptr, RegAS del_key,
    //     Duoram<Node>::Flat &A, RegBS F_af, RegBS F_fs, int TTL, 
    //     del_return &ret_struct);

  public:
    HEAP(int num_players, size_t size) {
      this->initialize(num_players, size);
    };

    ~HEAP() {
      if(oram)
        delete oram;
    };

    void initialize(int num_players, size_t size);
    RegAS extract_min(MPCTIO tio, yield_t &yield, std::vector<RegAS> A);
    int insert(MPCTIO tio, yield_t &yield, std::vector<RegAS>& A, RegAS val);
    int verify_heap_property(MPCTIO tio, yield_t &yield, std::vector<RegAS> A);

    // void insert(MPCTIO &tio, yield_t &yield, Node &node);
    // int del(MPCTIO &tio, yield_t &yield, RegAS del_key); 

    // // Display and correctness check functions
    // void pretty_print(MPCTIO &tio, yield_t &yield);
    // void pretty_print(const std::vector<Node> &R, value_t node,
    //     const std::string &prefix, bool is_left_child, bool is_right_child);
    // void check_bst(MPCTIO &tio, yield_t &yield);
    // std::tuple<bool, address_t> check_bst(const std::vector<Node> &R,
    //     value_t node, value_t min_key, value_t max_key);

};

 void Heap(MPCIO &mpcio,
    const PRACOptions &opts, char **args);
 
 //#include "heap.tcc"
#endif
