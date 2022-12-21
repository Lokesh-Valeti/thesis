#ifndef __DUORAM_HPP__
#define __DUORAM_HPP__

// Implementation of the 3-party protocols described in:
// Adithya Vadapalli, Ryan Henry, Ian Goldberg, "Duoram: A
// Bandwidth-Efficient Distributed ORAM for 2- and 3-Party Computation".

// A Duoram object is like physical memory: it's just a flat address
// space, and you can't access it directly.  Instead, you need to access
// it through a "Shape", such as Flat, Tree, Path, etc.  Shapes can be
// nested, so you can have a Path of a Subtree of a Tree sitting on the
// base Duoram.  Each Shape's parent must remain in scope (references to
// it must remain valid) for the lifetime of the child Shapre.  Each
// shape is bound to a context, which is a thread-specific MPCTIO and a
// coroutine-specific yield_t.  If you launch new threads and/or
// coroutines, you'll need to make a copy of the current Shape for your
// new context, and call set_context() on it.  Be sure not to call
// set_content() on a Shape shared with other threads or coroutines.

// This is templated, because you can have a Duoram of additively shared
// (RegAS) or XOR shared (RegXS) elements, or std::arrays of those to
// get "wide" memory cells.

// The initial implementation is focused on additive shares.

template <typename T>
class Duoram {
    // The computational parties have three vectors: the share of the
    // database itself, the party's own blinding factors for its
    // database share, and the _other_ computational party's blinded
    // database share (its database share plus its blind).

    // The player number (0 and 1 for the computational parties and 2
    // for the server) and the size of the Duoram
    int player;
    size_t oram_size;

    // The server has two vectors: a copy of each computational party's
    // blind.  The database vector will remain empty.

    std::vector<T> database;         // computational parties only

    std::vector<T> blind;            // computational parties use this name
    std::vector<T> &p0_blind;        // server uses this name

    std::vector<T> peer_blinded_db;  // computational parties
    std::vector<T> &p1_blind;        // server

public:
    // The type of this Duoram
    using type = T;

    // Pass the player number and desired size
    Duoram(int player, size_t size);

    // Get the size
    inline size_t size() { return oram_size; }
};

#include "duoram.tcc"

#endif
