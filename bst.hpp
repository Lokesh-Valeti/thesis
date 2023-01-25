#ifndef __NODE_HPP__
#define __NODE_HPP__

#include "mpcio.hpp"
#include "options.hpp"
#include "cell.hpp"

struct Node: public Cell {
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

void bst(MPCIO &mpcio,
    const PRACOptions &opts, char **args);

#endif
