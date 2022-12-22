#ifndef __DUORAM_HPP__
#define __DUORAM_HPP__

#include "types.hpp"

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
// new context, and call context() on it.  Be sure not to call context()
// on a Shape shared with other threads or coroutines.

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

    // The different Shapes are subclasses of this inner class
    class Shape;
    // These are the different Shapes that exist
    class Flat;

    // Pass the player number and desired size
    Duoram(int player, size_t size);

    // Get the size
    inline size_t size() { return oram_size; }

    // Get the basic Flat shape for this Duoram
    Flat flat(MPCTIO &tio, yield_t &yield, size_t start = 0,
            size_t len = 0) {
        return Flat(*this, tio, yield, start, len);
    }

    // For debugging; print the contents of the Duoram to stdout
    void dump() const;
};

// The parent class of all Shapes.  This is an abstract class that
// cannot itself be instantiated.

template <typename T>
class Duoram<T>::Shape {
    // Subclasses should be able to access _other_ Shapes' indexmap
    friend class Flat;

    // When you index into a shape (A[x]), you get one of these types,
    // depending on the type of x (the index), _not_ on the type T (the
    // underlying type of the Duoram).  That is, you can have an
    // additive-shared index (x) into an XOR-shared database (T), for
    // example.

    // When x is unshared explicit value
    class MemRefExpl;
    // When x is additively shared
    class MemRefAS;
    // When x is XOR shared
    class MemRefXS;

protected:
    // A reference to the parent shape.  As with ".." in the root
    // directory of a filesystem, the topmost shape is indicated by
    // having parent = *this.
    const Shape &parent;

    // A reference to the backing physical storage
    Duoram &duoram;

    // The size of this shape
    size_t shape_size;

    // The number of bits needed to address this shape (the number of
    // bits in shape_size-1)
    nbits_t addr_size;

    // And a mask with the low addr_size bits set
    address_t addr_mask;

    // The Shape's context (MPCTIO and yield_t)
    MPCTIO &tio;
    yield_t &yield;

    // A function to set the shape_size and compute addr_size and
    // addr_mask
    void set_shape_size(size_t sz);

    // We need a constructor because we hold non-static references; this
    // constructor is called by the subclass constructors
    Shape(const Shape &parent, Duoram &duoram, MPCTIO &tio,
        yield_t &yield) : parent(parent), duoram(duoram), shape_size(0),
        tio(tio), yield(yield) {}

    // The index-mapping function. Input the index relative to this
    // shape, and output the corresponding physical address.  The
    // strategy is to map the index relative to this shape to the index
    // relative to the parent shape, call the parent's indexmap function
    // on that (unless this is the topmost shape), and return what it
    // returns.  If this is the topmost shape, just return what you
    // would have passed to the parent's indexmap.
    //
    // This is a pure virtual function; all subclasses of Shape must
    // implement it, and of course Shape itself therefore cannot be
    // instantiated.
    virtual size_t indexmap(size_t idx) const = 0;

    // Get a pair (for the server) of references to the underlying
    // Duoram entries at share virtual index idx.  (That is, it gets
    // duoram.p0_blind[indexmap(idx)], etc.)
    inline std::tuple<T&,T&> get_server(size_t idx) const {
        size_t physaddr = indexmap(idx);
        return std::tie(
            duoram.p0_blind[physaddr],
            duoram.p1_blind[physaddr]);
    }

    // Get a triple (for the computational players) of references to the
    // underlying Duoram entries at share virtual index idx.  (That is,
    // it gets duoram.database[indexmap(idx)], etc.)
    inline std::tuple<T&,T&,T&> get_comp(size_t idx) const {
        size_t physaddr = indexmap(idx);
        return std::tie(
            duoram.database[physaddr],
            duoram.blind[physaddr],
            duoram.peer_blinded_db[physaddr]);
    }

public:
    // Get the size
    inline size_t size() { return shape_size; }

    // Update the context (MPCTIO and yield if you've started a new
    // thread, or just yield if you've started a new coroutine in the
    // same thread)
    void context(MPCTIO &new_tio, yield_t &new_yield) {
        tio = new_tio;
        yield = new_yield;
    }
    void context(yield_t &new_yield) { yield = new_yield; }

    // Index into this Shape in various ways
    MemRefAS operator[](const RegAS &idx) { return MemRefAS(*this, idx); }
    MemRefXS operator[](const RegXS &idx) { return MemRefXS(*this, idx); }
    MemRefExpl operator[](address_t idx) { return MemRefExpl(*this, idx); }

    // For debugging or checking your answers (using this in general is
    // of course insecure)

    // This one reconstructs the whole database
    std::vector<T> reconstruct() const;

    // This one reconstructs a single database value
    T reconstruct(const T& share) const;
};

// The most basic shape is Flat.  It is almost always the topmost shape,
// and serves to provide MPCTIO and yield_t context to a Duoram without
// changing the indices or size (but can specify a subrange if desired).

template <typename T>
class Duoram<T>::Flat : public Duoram<T>::Shape {
    // If this is a subrange, start may be non-0, but it's usually 0
    size_t start;

    inline size_t indexmap(size_t idx) const {
        size_t paridx = idx + start;
        if (&(this->parent) == this) {
            return paridx;
        } else {
            return this->parent.indexmap(paridx);
        }
    }

public:
    // Constructor.  len=0 means the maximum size (the parent's size
    // minus start).
    Flat(Duoram &duoram, MPCTIO &tio, yield_t &yield, size_t start = 0,
        size_t len = 0);
};

// An additively shared memory reference.  You get one of these from a
// Shape A and an additively shared RegAS index x with A[x].  Then you
// perform operations on this object, which do the Duoram operations.

template <typename T>
class Duoram<T>::Shape::MemRefAS {
    const Shape &shape;
    RegAS idx;

public:
    MemRefAS(const Shape &shape, const RegAS &idx) :
        shape(shape), idx(idx) {}

    // Oblivious read from an additively shared index of Duoram memory
    operator T();

    // Oblivious update to an additively shared index of Duoram memory
    MemRefAS &operator+=(const T& M);
};

// An XOR shared memory reference.  You get one of these from a Shape A
// and an XOR shared RegXS index x with A[x].  Then you perform
// operations on this object, which do the Duoram operations.

template <typename T>
class Duoram<T>::Shape::MemRefXS {
    const Shape &shape;
    RegXS idx;

public:
    MemRefXS(const Shape &shape, const RegXS &idx) :
        shape(shape), idx(idx) {}
};

// An explicit memory reference.  You get one of these from a Shape A
// and an address_t index x with A[x].  Then you perform operations on
// this object, which update the Duoram state without performing Duoram
// operations.

template <typename T>
class Duoram<T>::Shape::MemRefExpl {
    const Shape &shape;
    address_t idx;

public:
    MemRefExpl(const Shape &shape, address_t idx) :
        shape(shape), idx(idx) {}

    // Explicit read from a given index of Duoram memory
    operator T();

    // Explicit update to a given index of Duoram memory
    MemRefExpl &operator+=(const T& M);
};

#include "duoram.tcc"

#endif
