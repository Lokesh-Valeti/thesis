#include <functional>

#include "types.hpp"
#include "duoram.hpp"
#include "bst.hpp"

// This file demonstrates how to implement custom ORAM wide cell types.
// Such types can be structures of arbitrary numbers of RegAS and RegXS
// fields.  The example here imagines a node of a binary search tree,
// where you would want the key to be additively shared (so that you can
// easily do comparisons), the pointers field to be XOR shared (so that
// you can easily do bit operations to pack two pointers and maybe some
// tree balancing information into one field) and the value doesn't
// really matter, but XOR shared is usually slightly more efficient.

struct Node {
    RegAS key;
    RegXS pointers;
    RegXS value;

// Field-access macros so we can write A[i].NODE_KEY instead of
// A[i].field(&Node::key)

#define NODE_KEY field(&Node::key)
#define NODE_POINTERS field(&Node::pointers)
#define NODE_VALUE field(&Node::value)

    // For debugging and checking answers
    void dump() const {
        printf("[%016lx %016lx %016lx]", key.share(), pointers.share(),
            value.share());
    }

    // You'll need to be able to create a random element, and do the
    // operations +=, +, -=, - (binary and unary).  Note that for
    // XOR-shared fields, + and - are both really XOR.

    inline void randomize() {
        key.randomize();
        pointers.randomize();
        value.randomize();
    }

    inline Node &operator+=(const Node &rhs) {
        this->key += rhs.key;
        this->pointers += rhs.pointers;
        this->value += rhs.value;
        return *this;
    }

    inline Node operator+(const Node &rhs) const {
        Node res = *this;
        res += rhs;
        return res;
    }

    inline Node &operator-=(const Node &rhs) {
        this->key -= rhs.key;
        this->pointers -= rhs.pointers;
        this->value -= rhs.value;
        return *this;
    }

    inline Node operator-(const Node &rhs) const {
        Node res = *this;
        res -= rhs;
        return res;
    }

    inline Node operator-() const {
        Node res;
        res.key = -this->key;
        res.pointers = -this->pointers;
        res.value = -this->value;
        return res;
    }

    // Multiply each field by the local share of the corresponding field
    // in the argument
    inline Node mulshare(const Node &rhs) const {
        Node res = *this;
        res.key.mulshareeq(rhs.key);
        res.pointers.mulshareeq(rhs.pointers);
        res.value.mulshareeq(rhs.value);
        return res;
    }

    // You need a method to turn a leaf node of a DPF into a share of a
    // unit element of your type.  Typically set each RegAS to
    // dpf.unit_as(leaf) and each RegXS or RegBS to dpf.unit_bs(leaf).
    // Note that RegXS will extend a RegBS of 1 to the all-1s word, not
    // the word with value 1.  This is used for ORAM reads, where the
    // same DPF is used for all the fields.
    inline void unit(const RDPF &dpf, DPFnode leaf) {
        key = dpf.unit_as(leaf);
        pointers = dpf.unit_bs(leaf);
        value = dpf.unit_bs(leaf);
    }

    // Perform an update on each of the fields, using field-specific
    // MemRefs constructed from the Shape shape and the index idx
    template <typename Sh, typename U>
    inline static void update(Sh &shape, yield_t &shyield, U idx,
            const Node &M) {
        run_coroutines(shyield,
            [&shape, &idx, &M] (yield_t &yield) {
                Sh Sh_coro = shape.context(yield);
                Sh_coro[idx].NODE_KEY += M.key;
            },
            [&shape, &idx, &M] (yield_t &yield) {
                Sh Sh_coro = shape.context(yield);
                Sh_coro[idx].NODE_POINTERS += M.pointers;
            },
            [&shape, &idx, &M] (yield_t &yield) {
                Sh Sh_coro = shape.context(yield);
                Sh_coro[idx].NODE_VALUE += M.value;
            });
    }
};

// I/O operations (for sending over the network)

template <typename T>
T& operator>>(T& is, Node &x)
{
    is >> x.key >> x.pointers >> x.value;
    return is;
}

template <typename T>
T& operator<<(T& os, const Node &x)
{
    os << x.key << x.pointers << x.value;
    return os;
}

// This macro will define I/O on tuples of two or three of the cell type

DEFAULT_TUPLE_IO(Node)

std::tuple<RegBS, RegBS> compare_keys(Node n1, Node n2, MPCTIO tio, yield_t &yield) {
  CDPF cdpf = tio.cdpf(yield);
  auto [lt, eq, gt] = cdpf.compare(tio, yield, n2.key - n1.key, tio.aes_ops());
  RegBS lteq = lt^eq;
  return {lteq, gt};
}

RegBS check_ptr_zero(RegXS ptr, MPCTIO tio, yield_t &yield) {
  CDPF cdpf = tio.cdpf(yield);
  RegAS ptr_as;
  mpc_xs_to_as(tio, yield, ptr_as, ptr);
  RegAS zero;
  auto [lt, eq, gt] = cdpf.compare(tio, yield, ptr_as - zero, tio.aes_ops());
  return eq;
}

// Assuming pointer of 64 bits is split as:
// - 32 bits Left ptr
// - 32 bits Right ptr
// < Left, Right>

inline RegXS extractLeftPtr(RegXS pointer){ 
  return ((pointer&(0xFFFFFFFF00000000))>>32); 
}

inline RegXS extractRightPtr(RegXS pointer){
  return (pointer&(0x00000000FFFFFFFF)); 
}

inline void setLeftPtr(RegXS &pointer, RegXS new_ptr){ 
  pointer&=(0x00000000FFFFFFFF);
  pointer+=(new_ptr<<32);
}

inline void setRightPtr(RegXS &pointer, RegXS new_ptr){
  pointer&=(0xFFFFFFFF00000000);
  pointer+=(new_ptr);
}

std::tuple<RegXS, RegBS> insert(MPCTIO &tio, yield_t &yield, RegXS ptr, Node new_node, Duoram<Node>::Flat A, int TTL, RegBS isDummy) {
  if(TTL==0) {
    RegBS zero;
    return {ptr, zero};
  }

  Node cnode = A[ptr];
  //Compare key
  auto [lteq, gt] = compare_keys(cnode, new_node, tio, yield);

  // Depending on [lteq, gt] select the next ptr/index as
  // upper 32 bits of cnode.pointers if lteq
  // lower 32 bits of cnode.pointers if gt 
  RegXS left = extractLeftPtr(cnode.pointers);
  RegXS right = extractRightPtr(cnode.pointers);
  
  RegXS next_ptr;
  mpc_select(tio, yield, next_ptr, gt, left, right, 32);

  CDPF dpf = tio.cdpf(yield);
  size_t &aes_ops = tio.aes_ops();
  //RegAS next_ptr_as;
  //mpc_xs_to_as(tio, yield, next_ptr_as, next_ptr);
  RegBS F_z = dpf.is_zero(tio, yield, next_ptr, aes_ops); 
  RegBS F_i;

  if(tio.player()==0) {
    isDummy^=1;
  }
  mpc_and(tio, yield, F_i, (isDummy), F_z);
  if(tio.player()==0) {
    isDummy^=1;
  }
   
  isDummy^=F_i;
  auto [wptr, direction] = insert(tio, yield, next_ptr, new_node, A, TTL-1, isDummy);
  
  RegXS ret_ptr;
  RegBS ret_direction;
  mpc_select(tio, yield, ret_ptr, F_i, wptr, ptr);
  //ret_direction = direction + F_p(direction - gt)
  mpc_and(tio, yield, ret_direction, F_i, direction^gt);
  ret_direction^=direction;  

  return {ret_ptr, ret_direction};
}


// Insert(root, ptr, key, TTL, isDummy) -> (new_ptr, wptr, wnode, f_p)
void insert(MPCTIO &tio, yield_t &yield, RegXS &root, Node node, Duoram<Node>::Flat A, size_t &num_items) {
  if(num_items==0) {
    Node zero;
    A[0] = zero;
    A[1] = node;
    (root).set(1*tio.player());
    num_items++;
    return;
  }
  else {
    // Insert node into next free slot in the ORAM
    int new_id = 1 + num_items;
    int TTL = num_items++;
    A[new_id] = node;
    RegXS new_addr;
    new_addr.set(new_id * tio.player());
    RegBS isDummy;

    //Do a recursive insert
    auto [wptr, direction] = insert(tio, yield, root, node, A, TTL, isDummy);

    //Complete the insertion by reading wptr and updating its pointers
    RegXS pointers = A[wptr].NODE_POINTERS;
    RegXS left_ptr = extractLeftPtr(pointers);
    RegXS right_ptr = extractRightPtr(pointers);
    RegXS new_right_ptr, new_left_ptr;
    mpc_select(tio, yield, new_right_ptr, direction, right_ptr, new_addr);
    if(tio.player()==0) {
      direction^=1;
    }
    mpc_select(tio, yield, new_left_ptr, direction, left_ptr, new_addr);
    setLeftPtr(pointers, new_left_ptr);
    setRightPtr(pointers, new_right_ptr);
    A[wptr].NODE_POINTERS = pointers;
  }
  
}

void newnode(Node &a) {
  a.key.randomize(8);
  a.pointers.set(0);
  a.value.randomize();
}

// Now we use the node in various ways.  This function is called by
// online.cpp.

void bst(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
{
    nbits_t depth=5;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }

    MPCTIO tio(mpcio, 0, opts.num_threads);
    run_coroutines(tio, [&tio, depth] (yield_t &yield) {
        size_t size = size_t(1)<<depth;
        Duoram<Node> oram(tio.player(), size);
        auto A = oram.flat(tio, yield);

        size_t num_items = 0;
        RegXS root;

        Node c; 
        for(int i = 0; i<30; i++) {
          newnode(c);
          insert(tio, yield, root, c, A, num_items);
        }
        
        if (depth < 10) {
            oram.dump();
            auto R = A.reconstruct();
            if (tio.player() == 0) {
                for(size_t i=0;i<R.size();++i) {
                    printf("\n%04lx ", i);
                    R[i].dump();
                }
                printf("\n");
            }
        }
    });
}
