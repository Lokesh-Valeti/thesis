// Templated method implementations for duoram.hpp

#include <stdio.h>

#include "cdpf.hpp"

// Pass the player number and desired size
template <typename T>
Duoram<T>::Duoram(int player, size_t size) : player(player),
        oram_size(size), p0_blind(blind), p1_blind(peer_blinded_db) {
    if (player < 2) {
        database.resize(size);
        blind.resize(size);
        peer_blinded_db.resize(size);
    } else {
        p0_blind.resize(size);
        p1_blind.resize(size);
    }
}

// For debugging; print the contents of the Duoram to stdout
template <typename T>
void Duoram<T>::dump() const
{
    for (size_t i=0; i<oram_size; ++i) {
        if (player < 2) {
            printf("%04lx %016lx %016lx %016lx\n",
                i, database[i].share(), blind[i].share(),
                peer_blinded_db[i].share());
        } else {
            printf("%04lx %016lx %016lx\n",
                i, p0_blind[i].share(), p1_blind[i].share());
        }
    }
    printf("\n");
}

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
template <typename T>
void Duoram<T>::Shape::explicitonly(bool enable)
{
    if (enable == true)  {
        explicitmode = true;
    } else if (explicitmode == true) {
        explicitmode = false;
        // Reblind the whole Shape
        int player = tio.player();
        if (player < 2) {
            for (size_t i=0; i<shape_size; ++i) {
                auto [ DB, BL, PBD ] = get_comp(i);
                BL.randomize();
                tio.iostream_server() << BL;
                tio.iostream_peer() << (DB + BL);
            }
            yield();
            for (size_t i=0; i<shape_size; ++i) {
                auto [ DB, BL, PBD ] = get_comp(i);
                tio.iostream_peer() >> PBD;
            }
        } else {
            yield();
            for (size_t i=0; i<shape_size; ++i) {
                auto [BL0, BL1] = get_server(i);
                tio.iostream_p0() >> BL0;
                tio.iostream_p1() >> BL1;
            }
        }
    }

}

// For debugging or checking your answers (using this in general is
// of course insecure)
// This one reconstructs the whole database
template <typename T>
std::vector<T> Duoram<T>::Shape::reconstruct() const
{
    int player = tio.player();
    std::vector<T> res;
    res.resize(duoram.size());
    // Player 1 sends their share of the database to player 0
    if (player == 1) {
        tio.queue_peer(duoram.database.data(), duoram.size()*sizeof(T));
    } else if (player == 0) {
        tio.recv_peer(res.data(), duoram.size()*sizeof(T));
        for(size_t i=0;i<duoram.size();++i) {
            res[i] += duoram.database[i];
        }
    }
    // The server (player 2) does nothing

    // Players 1 and 2 will get an empty vector here
    return res;
}

// This one reconstructs a single database value
template <typename T>
T Duoram<T>::Shape::reconstruct(const T& share) const
{
    int player = tio.player();
    T res;

    // Player 1 sends their share of the value to player 0
    if (player == 1) {
        tio.queue_peer(&share, sizeof(T));
    } else if (player == 0) {
        tio.recv_peer(&res, sizeof(T));
        res += share;
    }
    // The server (player 2) does nothing

    // Players 1 and 2 will get 0 here
    return res;
}

// Function to set the shape_size of a shape and compute the number of
// bits you need to address a shape of that size (which is the number of
// bits in sz-1).  This is typically called by subclass constructors.
template <typename T>
void Duoram<T>::Shape::set_shape_size(size_t sz)
{
    shape_size = sz;
    // Compute the number of bits in (sz-1)
    // But use 0 if sz=0 for some reason (though that should never
    // happen)
    if (sz > 1) {
        addr_size = 64-__builtin_clzll(sz-1);
        addr_mask = address_t((size_t(1)<<addr_size)-1);
    } else {
        addr_size = 0;
        addr_mask = 0;
    }
}

// Constructor for the Flat shape.  len=0 means the maximum size (the
// parent's size minus start).
template <typename T>
Duoram<T>::Flat::Flat(Duoram &duoram, MPCTIO &tio, yield_t &yield,
    size_t start, size_t len) : Shape(*this, duoram, tio, yield)
{
    size_t parentsize = duoram.size();
    if (start > parentsize) {
        start = parentsize;
    }
    this->start = start;
    size_t maxshapesize = parentsize - start;
    if (len > maxshapesize || len == 0) {
        len = maxshapesize;
    }
    this->len = len;
    this->set_shape_size(len);
}

// Bitonic sort the elements from start to start+(1<<depth)-1, in
// increasing order if dir=0 or decreasing order if dir=1. Note that
// the elements must be at most 63 bits long each for the notion of
// ">" to make consistent sense.
template <typename T>
void Duoram<T>::Flat::bitonic_sort(address_t start, nbits_t depth, bool dir)
{
    if (depth == 0) return;
    if (depth == 1) {
        osort(start, start+1, dir);
        return;
    }
    // Recurse on the first half (increasing order) and the second half
    // (decreasing order) in parallel
    run_coroutines(this->yield,
        [this, start, depth](yield_t &yield) {
            Flat Acoro = context(yield);
            Acoro.bitonic_sort(start, depth-1, 0);
        },
        [this, start, depth](yield_t &yield) {
            Flat Acoro = context(yield);
            Acoro.bitonic_sort(start+(1<<(depth-1)), depth-1, 1);
        });
    // Merge the two into the desired order
    butterfly(start, depth, dir);
}

// Internal function to aid bitonic_sort
template <typename T>
void Duoram<T>::Flat::butterfly(address_t start, nbits_t depth, bool dir)
{
    if (depth == 0) return;
    if (depth == 1) {
        osort(start, start+1, dir);
        return;
    }
    // Sort pairs of elements half the width apart in parallel
    address_t halfwidth = address_t(1)<<(depth-1);
    std::vector<coro_t> coroutines;
    for (address_t i=0; i<halfwidth;++i) {
        coroutines.emplace_back(
            [this, start, halfwidth, dir, i](yield_t &yield) {
                Flat Acoro = context(yield);
                Acoro.osort(start+i, start+i+halfwidth, dir);
            });
    }
    run_coroutines(this->yield, coroutines);
    // Recurse on each half in parallel
    run_coroutines(this->yield,
        [this, start, depth, dir](yield_t &yield) {
            Flat Acoro = context(yield);
            Acoro.butterfly(start, depth-1, dir);
        },
        [this, start, halfwidth, depth, dir](yield_t &yield) {
            Flat Acoro = context(yield);
            Acoro.butterfly(start+halfwidth, depth-1, dir);
        });
}

// Assuming the memory is already sorted, do an oblivious binary
// search for the largest index containing the value at most the
// given one.  (The answer will be 0 if all of the memory elements
// are greate than the target.) This Flat must be a power of 2 size.
// Only available for additive shared databases for now.
template <>
RegAS Duoram<RegAS>::Flat::obliv_binary_search(RegAS &target)
{
    nbits_t depth = this->addr_size;
    // Start in the middle
    RegAS index;
    index.set(this->tio.player() ? 0 : 1<<(depth-1));
    // Invariant: index points to the first element of the right half of
    // the remaining possible range
    while (depth > 0) {
        // Obliviously read the value there
        RegAS val = operator[](index);
        // Compare it to the target
        CDPF cdpf = tio.cdpf(this->yield);
        auto [lt, eq, gt] = cdpf.compare(this->tio, this->yield,
            val-target, tio.aes_ops());
        if (depth > 1) {
            // If val > target, the answer is strictly to the left
            // and we should subtract 2^{depth-2} from index
            // If val <= target, the answer is here or to the right
            // and we should add 2^{depth-2} to index
            // So we unconditionally subtract 2^{depth-2} from index, and
            // add (lt+eq)*2^{depth-1}.
            RegAS uncond;
            uncond.set(tio.player() ? 0 : address_t(1)<<(depth-2));
            RegAS cond;
            cond.set(tio.player() ? 0 : address_t(1)<<(depth-1));
            RegAS condprod;
            RegBS le = lt ^ eq;
            mpc_flagmult(this->tio, this->yield, condprod, le, cond);
            index -= uncond;
            index += condprod;
        } else {
            // If val > target, the answer is strictly to the left
            // If val <= target, the answer is here or to the right
            // so subtract gt from index
            RegAS cond;
            cond.set(tio.player() ? 0 : 1);
            RegAS condprod;
            mpc_flagmult(this->tio, this->yield, condprod, gt, cond);
            index -= condprod;
        }
        --depth;
    }

    return index;
}

// Helper functions to specialize the read and update operations for
// RegAS and RegXS shared indices
template <typename U>
inline address_t IfRegAS(address_t val);
template <typename U>
inline address_t IfRegXS(address_t val);

template <>
inline address_t IfRegAS<RegAS>(address_t val) { return val; }
template <>
inline address_t IfRegAS<RegXS>(address_t val) { return 0; }
template <>
inline address_t IfRegXS<RegAS>(address_t val) { return 0; }
template <>
inline address_t IfRegXS<RegXS>(address_t val) { return val; }

// Oblivious read from an additively or XOR shared index of Duoram memory
// T is the sharing type of the _values_ in the database; U is the
// sharing type of the _indices_ in the database.
template <typename T> template <typename U>
Duoram<T>::Shape::MemRefS<U>::operator T()
{
    T res;
    Shape &shape = this->shape;
    shape.explicitonly(false);
    int player = shape.tio.player();
    if (player < 2) {
        // Computational players do this

        RDPFTriple dt = shape.tio.rdpftriple(shape.yield, shape.addr_size);

        // Compute the index offset
        U indoffset = dt.target<U>();
        indoffset -= idx;

        // We only need two of the DPFs for reading
        RDPFPair dp(std::move(dt), 0, player == 0 ? 2 : 1);
        // The RDPFTriple dt is now broken, since we've moved things out
        // of it.

        // Send it to the peer and the server
        shape.tio.queue_peer(&indoffset, BITBYTES(shape.addr_size));
        shape.tio.queue_server(&indoffset, BITBYTES(shape.addr_size));

        shape.yield();

        // Receive the above from the peer
        U peerindoffset;
        shape.tio.recv_peer(&peerindoffset, BITBYTES(shape.addr_size));

        // Reconstruct the total offset
        auto indshift = combine(indoffset, peerindoffset, shape.addr_size);

        // Evaluate the DPFs and compute the dotproducts
        ParallelEval pe(dp, IfRegAS<U>(indshift), IfRegXS<U>(indshift),
            shape.shape_size, shape.tio.cpu_nthreads(),
            shape.tio.aes_ops());
        T init;
        res = pe.reduce(init, [&dp, &shape] (int thread_num, address_t i,
                const RDPFPair::node &leaf) {
            // The values from the two DPFs
            auto [V0, V1] = dp.unit<T>(leaf);
            // References to the appropriate cells in our database, our
            // blind, and our copy of the peer's blinded database
            auto [DB, BL, PBD] = shape.get_comp(i);
            return (DB + PBD) * V0.share() - BL * (V1-V0).share();
        });

        shape.yield();

        // Receive the cancellation term from the server
        T gamma;
        shape.tio.iostream_server() >> gamma;
        res += gamma;
    } else {
        // The server does this

        RDPFPair dp = shape.tio.rdpfpair(shape.yield, shape.addr_size);
        U p0indoffset, p1indoffset;

        shape.yield();

        // Receive the index offset from the computational players and
        // combine them
        shape.tio.recv_p0(&p0indoffset, BITBYTES(shape.addr_size));
        shape.tio.recv_p1(&p1indoffset, BITBYTES(shape.addr_size));
        auto indshift = combine(p0indoffset, p1indoffset, shape.addr_size);

        // Evaluate the DPFs to compute the cancellation terms
        std::tuple<T,T> init, gamma;
        ParallelEval pe(dp, IfRegAS<U>(indshift), IfRegXS<U>(indshift),
            shape.shape_size, shape.tio.cpu_nthreads(),
            shape.tio.aes_ops());
        gamma = pe.reduce(init, [&dp, &shape] (int thread_num, address_t i,
                const RDPFPair::node &leaf) {
            // The values from the two DPFs
            auto [V0, V1] = dp.unit<T>(leaf);

            // shape.get_server(i) returns a pair of references to the
            // appropriate cells in the two blinded databases
            auto [BL0, BL1] = shape.get_server(i);
            return std::make_tuple(-BL0 * V1.share(), -BL1 * V0.share());
        });

        // Choose a random blinding factor
        T rho;
        rho.randomize();

        std::get<0>(gamma) += rho;
        std::get<1>(gamma) -= rho;

        // Send the cancellation terms to the computational players
        shape.tio.iostream_p0() << std::get<0>(gamma);
        shape.tio.iostream_p1() << std::get<1>(gamma);

        shape.yield();
    }
    return res;  // The server will always get 0
}

// Oblivious update to an additively or XOR shared index of Duoram memory
template <typename T> template <typename U>
typename Duoram<T>::Shape::template MemRefS<U>
    &Duoram<T>::Shape::MemRefS<U>::operator+=(const T& M)
{
    Shape &shape = this->shape;
    shape.explicitonly(false);
    int player = shape.tio.player();
    if (player < 2) {
        // Computational players do this

        RDPFTriple dt = shape.tio.rdpftriple(shape.yield, shape.addr_size);

        // Compute the index and message offsets
        U indoffset = dt.target<U>();
        indoffset -= idx;
        auto Moffset = std::make_tuple(M, M, M);
        Moffset -= dt.scaled_value<T>();

        // Send them to the peer, and everything except the first offset
        // to the server
        shape.tio.queue_peer(&indoffset, BITBYTES(shape.addr_size));
        shape.tio.iostream_peer() << Moffset;
        shape.tio.queue_server(&indoffset, BITBYTES(shape.addr_size));
        shape.tio.iostream_server() << std::get<1>(Moffset) <<
            std::get<2>(Moffset);

        shape.yield();

        // Receive the above from the peer
        U peerindoffset;
        std::tuple<T,T,T> peerMoffset;
        shape.tio.recv_peer(&peerindoffset, BITBYTES(shape.addr_size));
        shape.tio.iostream_peer() >> peerMoffset;

        // Reconstruct the total offsets
        auto indshift = combine(indoffset, peerindoffset, shape.addr_size);
        auto Mshift = combine(Moffset, peerMoffset);

        // Evaluate the DPFs and add them to the database
        ParallelEval pe(dt, IfRegAS<U>(indshift), IfRegXS<U>(indshift),
            shape.shape_size, shape.tio.cpu_nthreads(),
            shape.tio.aes_ops());
        int init = 0;
        pe.reduce(init, [&dt, &shape, &Mshift, player] (int thread_num,
                address_t i, const RDPFTriple::node &leaf) {
            // The values from the three DPFs
            auto [V0, V1, V2] = dt.scaled<T>(leaf) + dt.unit<T>(leaf) * Mshift;
            // References to the appropriate cells in our database, our
            // blind, and our copy of the peer's blinded database
            auto [DB, BL, PBD] = shape.get_comp(i);
            DB += V0;
            if (player == 0) {
                BL -= V1;
                PBD += V2-V0;
            } else {
                BL -= V2;
                PBD += V1-V0;
            }
            return 0;
        });
    } else {
        // The server does this

        RDPFPair dp = shape.tio.rdpfpair(shape.yield, shape.addr_size);
        U p0indoffset, p1indoffset;
        std::tuple<T,T> p0Moffset, p1Moffset;

        // Receive the index and message offsets from the computational
        // players and combine them
        shape.tio.recv_p0(&p0indoffset, BITBYTES(shape.addr_size));
        shape.tio.iostream_p0() >> p0Moffset;
        shape.tio.recv_p1(&p1indoffset, BITBYTES(shape.addr_size));
        shape.tio.iostream_p1() >> p1Moffset;
        auto indshift = combine(p0indoffset, p1indoffset, shape.addr_size);
        auto Mshift = combine(p0Moffset, p1Moffset);

        // Evaluate the DPFs and subtract them from the blinds
        ParallelEval pe(dp, IfRegAS<U>(indshift), IfRegXS<U>(indshift),
            shape.shape_size, shape.tio.cpu_nthreads(),
            shape.tio.aes_ops());
        int init = 0;
        pe.reduce(init, [&dp, &shape, &Mshift] (int thread_num,
                address_t i, const RDPFPair::node &leaf) {
            // The values from the two DPFs
            auto V = dp.scaled<T>(leaf) + dp.unit<T>(leaf) * Mshift;
            // shape.get_server(i) returns a pair of references to the
            // appropriate cells in the two blinded databases, so we can
            // subtract the pair directly.
            shape.get_server(i) -= V;
            return 0;
        });
    }
    return *this;
}

// Oblivious write to an additively or XOR shared index of Duoram memory
template <typename T> template <typename U>
typename Duoram<T>::Shape::template MemRefS<U>
    &Duoram<T>::Shape::MemRefS<U>::operator=(const T& M)
{
    T oldval = *this;
    T update = M - oldval;
    *this += update;
    return *this;
}

// Oblivious sort with the provided other element.  Without
// reconstructing the values, *this will become a share of the
// smaller of the reconstructed values, and other will become a
// share of the larger.
//
// Note: this only works for additively shared databases
template <> template <typename U,typename V>
void Duoram<RegAS>::Flat::osort(const U &idx1, const V &idx2, bool dir)
{
    // Load the values in parallel
    RegAS val1, val2;
    run_coroutines(yield,
        [this, &idx1, &val1](yield_t &yield) {
            Flat Acoro = context(yield);
            val1 = Acoro[idx1];
        },
        [this, &idx2, &val2](yield_t &yield) {
            Flat Acoro = context(yield);
            val2 = Acoro[idx2];
        });
    // Get a CDPF
    CDPF cdpf = tio.cdpf(yield);
    // Use it to compare the values
    RegAS diff = val1-val2;
    auto [lt, eq, gt] = cdpf.compare(tio, yield, diff, tio.aes_ops());
    RegBS cmp = dir ? lt : gt;
    // Get additive shares of cmp*diff
    RegAS cmp_diff;
    mpc_flagmult(tio, yield, cmp_diff, cmp, diff);
    // Update the two locations in parallel
    run_coroutines(yield,
        [this, &idx1, &cmp_diff](yield_t &yield) {
            Flat Acoro = context(yield);
            Acoro[idx1] -= cmp_diff;
        },
        [this, &idx2, &cmp_diff](yield_t &yield) {
            Flat Acoro = context(yield);
            Acoro[idx2] += cmp_diff;
        });
}

// Explicit read from a given index of Duoram memory
template <typename T>
Duoram<T>::Shape::MemRefExpl::operator T()
{
    Shape &shape = this->shape;
    T res;
    int player = shape.tio.player();
    if (player < 2) {
        res = std::get<0>(shape.get_comp(idx));
    }
    return res;  // The server will always get 0
}

// Explicit update to a given index of Duoram memory
template <typename T>
typename Duoram<T>::Shape::MemRefExpl
    &Duoram<T>::Shape::MemRefExpl::operator+=(const T& M)
{
    Shape &shape = this->shape;
    int player = shape.tio.player();
    // In explicit-only mode, just update the local DB; we'll sync the
    // blinds and the blinded DB when we leave explicit-only mode.
    if (shape.explicitmode) {
        if (player < 2) {
            auto [ DB, BL, PBD ] = shape.get_comp(idx);
            DB += M;
        }
        return *this;
    }
    if (player < 2) {
        // Computational players do this

        // Pick a blinding factor
        T blind;
        blind.randomize();

        // Send the blind to the server, and the blinded value to the
        // peer
        shape.tio.iostream_server() << blind;
        shape.tio.iostream_peer() << (M + blind);

        shape.yield();

        // Receive the peer's blinded value
        T peerblinded;
        shape.tio.iostream_peer() >> peerblinded;

        // Our database, our blind, the peer's blinded database
        auto [ DB, BL, PBD ] = shape.get_comp(idx);
        DB += M;
        BL += blind;
        PBD += peerblinded;
    } else if (player == 2) {
        // The server does this

        // Receive the updates to the blinds
        T p0blind, p1blind;
        shape.tio.iostream_p0() >> p0blind;
        shape.tio.iostream_p1() >> p1blind;

        // The two computational parties' blinds
        auto [ BL0, BL1 ] = shape.get_server(idx);
        BL0 += p0blind;
        BL1 += p1blind;
    }
    return *this;
}

// Explicit write to a given index of Duoram memory
template <typename T>
typename Duoram<T>::Shape::MemRefExpl
    &Duoram<T>::Shape::MemRefExpl::operator=(const T& M)
{
    T oldval = *this;
    T update = M - oldval;
    *this += update;
    return *this;
}

// Independent U-shared reads into a Shape of subtype Sh on a Duoram
// with values of sharing type T
template <typename T> template <typename U, typename Sh>
Duoram<T>::Shape::MemRefInd<U,Sh>::operator std::vector<T>()
{
    std::vector<T> res;
    size_t size = indcs.size();
    res.resize(size);
    std::vector<coro_t> coroutines;
    for (size_t i=0;i<size;++i) {
        coroutines.emplace_back([this, &res, i] (yield_t &yield) {
            Sh Sh_coro = shape.context(yield);
            res[i] = Sh_coro[indcs[i]];
        });
    }
    run_coroutines(shape.yield, coroutines);

    return res;
}

// Independent U-shared updates into a Shape of subtype Sh on a Duoram
// with values of sharing type T (vector version)
template <typename T> template <typename U, typename Sh>
typename Duoram<T>::Shape::template MemRefInd<U,Sh>
    &Duoram<T>::Shape::MemRefInd<U,Sh>::operator+=(const std::vector<T>& M)
{
    size_t size = indcs.size();
    assert(M.size() == size);

    std::vector<coro_t> coroutines;
    for (size_t i=0;i<size;++i) {
        coroutines.emplace_back([this, &M, i] (yield_t &yield) {
            Sh Sh_coro = shape.context(yield);
            Sh_coro[indcs[i]] += M[i];
        });
    }
    run_coroutines(shape.yield, coroutines);

    return *this;
}

// Independent U-shared updates into a Shape of subtype Sh on a Duoram
// with values of sharing type T (array version)
template <typename T> template <typename U, typename Sh> template <size_t N>
typename Duoram<T>::Shape::template MemRefInd<U,Sh>
    &Duoram<T>::Shape::MemRefInd<U,Sh>::operator+=(const std::array<T,N>& M)
{
    size_t size = indcs.size();
    assert(N == size);

    std::vector<coro_t> coroutines;
    for (size_t i=0;i<size;++i) {
        coroutines.emplace_back([this, &M, i] (yield_t &yield) {
            Sh Sh_coro = shape.context(yield);
            Sh_coro[indcs[i]] += M[i];
        });
    }
    run_coroutines(shape.yield, coroutines);

    return *this;
}

// Independent U-shared writes into a Shape of subtype Sh on a Duoram
// with values of sharing type T (vector version)
template <typename T> template <typename U, typename Sh>
typename Duoram<T>::Shape::template MemRefInd<U,Sh>
    &Duoram<T>::Shape::MemRefInd<U,Sh>::operator=(const std::vector<T>& M)
{
    size_t size = indcs.size();
    assert(M.size() == size);

    std::vector<coro_t> coroutines;
    for (size_t i=0;i<size;++i) {
        coroutines.emplace_back([this, &M, i] (yield_t &yield) {
            Sh Sh_coro = shape.context(yield);
            Sh_coro[indcs[i]] = M[i];
        });
    }
    run_coroutines(shape.yield, coroutines);

    return *this;
}

// Independent U-shared writes into a Shape of subtype Sh on a Duoram
// with values of sharing type T (array version)
template <typename T> template <typename U, typename Sh> template <size_t N>
typename Duoram<T>::Shape::template MemRefInd<U,Sh>
    &Duoram<T>::Shape::MemRefInd<U,Sh>::operator=(const std::array<T,N>& M)
{
    size_t size = indcs.size();
    assert(N == size);

    std::vector<coro_t> coroutines;
    for (size_t i=0;i<size;++i) {
        coroutines.emplace_back([this, &M, i] (yield_t &yield) {
            Sh Sh_coro = shape.context(yield);
            Sh_coro[indcs[i]] = M[i];
        });
    }
    run_coroutines(shape.yield, coroutines);

    return *this;
}
