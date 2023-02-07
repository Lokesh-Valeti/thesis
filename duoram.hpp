#ifndef __DUORAM_HPP__
#define __DUORAM_HPP__

#include "types.hpp"
#include "mpcio.hpp"
#include "coroutine.hpp"

// Implementation of the 3-party protocols described in:
// Adithya Vadapalli, Ryan Henry, Ian Goldberg, "Duoram: A
// Bandwidth-Efficient Distributed ORAM for 2- and 3-Party Computation".

// A Duoram object is like physical memory: it's just a flat address
// space, and you can't access it directly.  Instead, you need to access
// it through a "Shape", such as Flat, Tree, Path, etc.  Shapes can be
// nested, so you can have a Path of a Subtree of a Tree sitting on the
// base Duoram.  Each Shape's parent must remain in scope (references to
// it must remain valid) for the lifetime of the child Shape.  Each
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

    // When x is additively or XOR shared
    // U is the sharing type of the indices, while T is the sharing type
    // of the data in the database.  If we are referencing an entire
    // entry of type T, then the field type FT will equal T, and the
    // field selector type FST will be nullopt_t.  If we are referencing
    // a particular field of T, then FT will be the type of the field
    // (RegAS or RegXS) and FST will be a pointer-to-member T::* type
    // pointing to that field.  Sh is the specific Shape subtype used to
    // create the MemRefS.
    template <typename U, typename FT, typename FST, typename Sh>
    class MemRefS;
    // When x is unshared explicit value.  FT and FST are as above.
    template <typename FT, typename FST>
    class MemRefExpl;
    // When x is a vector or array of values of type U, used to denote a
    // collection of independent memory operations that can be performed
    // simultaneously.  Sh is the specific Shape subtype used to create
    // the MemRefInd.
    template <typename U, typename Sh>
    class MemRefInd;

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

    // If you enable explicit-only mode, sending updates of your blind
    // to the server and of your blinded database to your peer will be
    // temporarily disabled.  When you disable it (which will happen
    // automatically at the next ORAM read or write, or you can do it
    // explicitly), new random blinds will be chosen for the whole
    // Shape, and the blinds sent to the server, and the blinded
    // database sent to the peer.
    bool explicitmode;

    // A function to set the shape_size and compute addr_size and
    // addr_mask
    void set_shape_size(size_t sz);

    // We need a constructor because we hold non-static references; this
    // constructor is called by the subclass constructors
    Shape(const Shape &parent, Duoram &duoram, MPCTIO &tio,
        yield_t &yield) : parent(parent), duoram(duoram), shape_size(0),
        tio(tio), yield(yield), explicitmode(false) {}

    // Copy the given Shape except for the tio and yield
    Shape(const Shape &copy_from, MPCTIO &tio, yield_t &yield) :
        parent(copy_from.parent), duoram(copy_from.duoram),
        shape_size(copy_from.shape_size),
        addr_size(copy_from.addr_size), addr_mask(copy_from.addr_mask),
        tio(tio), yield(yield),
        explicitmode(copy_from.explicitmode) {}

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
    inline std::tuple<T&,T&> get_server(size_t idx,
        std::nullopt_t null = std::nullopt) const {
        size_t physaddr = indexmap(idx);
        return std::tie(
            duoram.p0_blind[physaddr],
            duoram.p1_blind[physaddr]);
    }

    // Get a triple (for the computational players) of references to the
    // underlying Duoram entries at share virtual index idx.  (That is,
    // it gets duoram.database[indexmap(idx)], etc.)
    inline std::tuple<T&,T&,T&> get_comp(size_t idx,
        std::nullopt_t null = std::nullopt) const {
        size_t physaddr = indexmap(idx);
        return std::tie(
            duoram.database[physaddr],
            duoram.blind[physaddr],
            duoram.peer_blinded_db[physaddr]);
    }

    // Get a pair (for the server) of references to a particular field
    // of the underlying Duoram entries at share virtual index idx.
    // (That is, it gets duoram.p0_blind[indexmap(idx)].field, etc.)
    template <typename FT>
    inline std::tuple<FT&,FT&> get_server(size_t idx, FT T::*field) const {
        size_t physaddr = indexmap(idx);
        return std::tie(
            duoram.p0_blind[physaddr].*field,
            duoram.p1_blind[physaddr].*field);
    }

    // Get a triple (for the computational players) of references to a
    // particular field to the underlying Duoram entries at share
    // virtual index idx.  (That is, it gets
    // duoram.database[indexmap(idx)].field, etc.)
    template <typename FT>
    inline std::tuple<FT&,FT&,FT&> get_comp(size_t idx, FT T::*field) const {
        size_t physaddr = indexmap(idx);
        return std::tie(
            duoram.database[physaddr].*field,
            duoram.blind[physaddr].*field,
            duoram.peer_blinded_db[physaddr].*field);
    }

public:
    // Get the size
    inline size_t size() { return shape_size; }

    // Enable or disable explicit-only mode.  Only using [] with
    // explicit (address_t) indices are allowed in this mode.  Using []
    // with RegAS or RegXS indices will automatically turn off this
    // mode, or you can turn it off explicitly.  In explicit-only mode,
    // updates to the memory in the Shape will not induce communication
    // to the server or peer, but when it turns off, a message of the
    // size of the entire Shape will be sent to each of the server and
    // the peer.  This is useful if you're going to be doing multiple
    // explicit writes to every element of the Shape before you do your
    // next oblivious read or write.  Bitonic sort is a prime example.
    void explicitonly(bool enable);

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
    size_t len;

    inline size_t indexmap(size_t idx) const {
        size_t paridx = idx + start;
        if (&(this->parent) == this) {
            return paridx;
        } else {
            return this->parent.indexmap(paridx);
        }
    }

    // Internal function to aid bitonic_sort
    void butterfly(address_t start, nbits_t depth, bool dir);

public:
    // Constructor.  len=0 means the maximum size (the parent's size
    // minus start).
    Flat(Duoram &duoram, MPCTIO &tio, yield_t &yield, size_t start = 0,
        size_t len = 0);

    // Copy the given Flat except for the tio and yield
    Flat(const Flat &copy_from, MPCTIO &tio, yield_t &yield) :
        Shape(copy_from, tio, yield), start(copy_from.start),
        len(copy_from.len) {}

    // Update the context (MPCTIO and yield if you've started a new
    // thread, or just yield if you've started a new coroutine in the
    // same thread).  Returns a new Shape with an updated context.
    Flat context(MPCTIO &new_tio, yield_t &new_yield) const {
        return Flat(*this, new_tio, new_yield);
    }
    Flat context(yield_t &new_yield) const {
        return Flat(*this, this->tio, new_yield);
    }

    // Index into this Flat in various ways
    typename Duoram::Shape::template MemRefS<RegAS,T,std::nullopt_t,Flat>
            operator[](const RegAS &idx) {
        typename Duoram<T>::Shape::
            template MemRefS<RegAS,T,std::nullopt_t,Flat>
            res(*this, idx, std::nullopt);
        return res;
    }
    typename Duoram::Shape::template MemRefS<RegXS,T,std::nullopt_t,Flat>
            operator[](const RegXS &idx) {
        typename Duoram<T>::Shape::
            template MemRefS<RegXS,T,std::nullopt_t,Flat>
            res(*this, idx, std::nullopt);
        return res;
    }
    typename Duoram::Shape::template MemRefExpl<T,std::nullopt_t>
            operator[](address_t idx) {
        typename Duoram<T>::Shape::
            template MemRefExpl<T,std::nullopt_t>
            res(*this, idx, std::nullopt);
        return res;
    }
    template <typename U>
    Duoram::Shape::MemRefInd<U, Flat>
            operator[](const std::vector<U> &indcs) {
        typename Duoram<T>::Shape::
            template MemRefInd<U,Flat>
            res(*this, indcs);
        return res;
    }
    template <typename U, size_t N>
    Duoram::Shape::MemRefInd<U, Flat>
            operator[](const std::array<U,N> &indcs) {
        typename Duoram<T>::Shape::
            template MemRefInd<U,Flat>
            res(*this, indcs);
        return res;
    }

    // Oblivious sort the elements indexed by the two given indices.
    // Without reconstructing the values, if dir=0, this[idx1] will
    // become a share of the smaller of the reconstructed values, and
    // this[idx2] will become a share of the larger.  If dir=1, it's the
    // other way around.
    //
    // Note: this only works for additively shared databases
    template<typename U,typename V>
    void osort(const U &idx1, const V &idx2, bool dir=0);

    // Bitonic sort the elements from start to start+(1<<depth)-1, in
    // increasing order if dir=0 or decreasing order if dir=1. Note that
    // the elements must be at most 63 bits long each for the notion of
    // ">" to make consistent sense.
    void bitonic_sort(address_t start, nbits_t depth, bool dir=0);

    // Assuming the memory is already sorted, do an oblivious binary
    // search for the largest index containing the value at most the
    // given one.  (The answer will be 0 if all of the memory elements
    // are greate than the target.) This Flat must be a power of 2 size.
    // Only available for additive shared databases for now.
    RegAS obliv_binary_search(RegAS &target);
};

// An additive or XOR shared memory reference.  You get one of these
// from a Shape A and an additively shared RegAS index x, or an XOR
// shared RegXS index x, with A[x].  Then you perform operations on this
// object, which do the Duoram operations.  As above, T is the sharing
// type of the data in the database, while U is the sharing type of the
// index used to create this memory reference.  If we are referencing an
// entire entry of type T, then the field type FT will equal T, and the
// field selector type FST will be nullopt_t.  If we are referencing a
// particular field of T, then FT will be the type of the field (RegAS
// or RegXS) and FST will be a pointer-to-member T::* type pointing to
// that field.  Sh is the specific Shape subtype used to create the
// MemRefS.

template <typename T>
template <typename U, typename FT, typename FST, typename Sh>
class Duoram<T>::Shape::MemRefS {
    Sh &shape;
    U idx;
    FST fieldsel;

private:
    // Oblivious update to a shared index of Duoram memory, only for
    // FT = RegAS or RegXS
    MemRefS<U,FT,FST,Sh> &oram_update(const FT& M, const prac_template_true&);
    // Oblivious update to a shared index of Duoram memory, for
    // FT not RegAS or RegXS
    MemRefS<U,FT,FST,Sh> &oram_update(const FT& M, const prac_template_false&);

public:
    MemRefS<U,FT,FST,Sh>(Sh &shape, const U &idx, FST fieldsel) :
        shape(shape), idx(idx), fieldsel(fieldsel) {}

    // Create a MemRefExpl for accessing a partcular field of T
    template <typename SFT>
    MemRefS<U,SFT,SFT T::*,Sh> field(SFT T::*subfieldsel) {
        auto res = MemRefS<U,SFT,SFT T::*,Sh>(this->shape, idx, subfieldsel);
        return res;
    }

    // Oblivious read from a shared index of Duoram memory
    operator FT();

    // Oblivious update to a shared index of Duoram memory
    MemRefS<U,FT,FST,Sh> &operator+=(const FT& M);

    // Oblivious write to a shared index of Duoram memory
    MemRefS<U,FT,FST,Sh> &operator=(const FT& M);
};

// An explicit memory reference.  You get one of these from a Shape A
// and an address_t index x with A[x].  Then you perform operations on
// this object, which update the Duoram state without performing Duoram
// operations.  If we are referencing an entire entry of type T, then
// the field type FT will equal T, and the field selector type FST will
// be nullopt_t.  If we are referencing a particular field of T, then FT
// will be the type of the field (RegAS or RegXS) and FST will be a
// pointer-to-member T::* type pointing to that field.

template <typename T> template <typename FT, typename FST>
class Duoram<T>::Shape::MemRefExpl {
    Shape &shape;
    address_t idx;
    FST fieldsel;

public:
    MemRefExpl(Shape &shape, address_t idx, FST fieldsel) :
        shape(shape), idx(idx), fieldsel(fieldsel) {}

    // Create a MemRefExpl for accessing a partcular field of T
    template <typename SFT>
    MemRefExpl<SFT,SFT T::*> field(SFT T::*subfieldsel) {
        auto res = MemRefExpl<SFT,SFT T::*>(this->shape, idx, subfieldsel);
        return res;
    }

    // Explicit read from a given index of Duoram memory
    operator FT();

    // Explicit update to a given index of Duoram memory
    MemRefExpl &operator+=(const FT& M);

    // Explicit write to a given index of Duoram memory
    MemRefExpl &operator=(const FT& M);

    // Convenience function
    MemRefExpl &operator-=(const FT& M) { *this += (-M); return *this; }
};

// A collection of independent memory references that can be processed
// simultaneously.  You get one of these from a Shape A (of specific
// subclass Sh) and a vector or array of indices v with each element of
// type U.

template <typename T> template <typename U, typename Sh>
class Duoram<T>::Shape::MemRefInd {
    Sh &shape;
    std::vector<U> indcs;

public:
    MemRefInd(Sh &shape, std::vector<U> indcs) :
        shape(shape), indcs(indcs) {}
    template <size_t N>
    MemRefInd(Sh &shape, std::array<U,N> aindcs) :
        shape(shape) { for ( auto &i : aindcs ) { indcs.push_back(i); } }

    // Independent reads from shared or explicit indices of Duoram memory
    operator std::vector<T>();

    // Independent updates to shared or explicit indices of Duoram memory
    MemRefInd &operator+=(const std::vector<T>& M);
    template <size_t N>
    MemRefInd &operator+=(const std::array<T,N>& M);

    // Independent writes to shared or explicit indices of Duoram memory
    MemRefInd &operator=(const std::vector<T>& M);
    template <size_t N>
    MemRefInd &operator=(const std::array<T,N>& M);

    // Convenience function
    MemRefInd &operator-=(const std::vector<T>& M) { *this += (-M); return *this; }
    template <size_t N>
    MemRefInd &operator-=(const std::array<T,N>& M) { *this += (-M); return *this; }
};

#include "duoram.tcc"

#endif
