// Templated method implementations for duoram.hpp

#include <stdio.h>

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
    this->set_shape_size(len);
}

// Oblivious read from an additively shared index of Duoram memory
template <typename T>
Duoram<T>::Shape::MemRefAS::operator T()
{
    T res;
    int player = shape.tio.player();
    if (player < 2) {
        // Computational players do this

        RDPFTriple dt = shape.tio.rdpftriple(shape.addr_size);

        // Compute the index offset
        RegAS indoffset = dt.as_target;
        indoffset -= idx;

        // We only need two of the DPFs for reading
        RDPFPair dp(std::move(dt), 0, player == 0 ? 2 : 1);

        // Send it to the peer and the server
        shape.tio.queue_peer(&indoffset, BITBYTES(shape.addr_size));
        shape.tio.queue_server(&indoffset, BITBYTES(shape.addr_size));

        shape.yield();

        // Receive the above from the peer
        RegAS peerindoffset;
        shape.tio.recv_peer(&peerindoffset, BITBYTES(shape.addr_size));

        // Reconstruct the total offset
        auto indshift = combine(indoffset, peerindoffset, shape.addr_size);

        // Evaluate the DPFs and compute the dotproducts
        StreamEval ev(dp, indshift, 0, shape.tio.aes_ops());
        for (size_t i=0; i<shape.shape_size; ++i) {
            auto L = ev.next();
            // The values from the two DPFs
            auto [V0, V1] = dp.unit<T>(L);
            // References to the appropriate cells in our database, our
            // blind, and our copy of the peer's blinded database
            auto [DB, BL, PBD] = shape.get_comp(i);
            res += (DB + PBD) * V0.share() - BL * (V1-V0).share();
        }

        // Receive the cancellation term from the server
        T gamma;
        shape.tio.iostream_server() >> gamma;
        res += gamma;
    } else {
        // The server does this

        RDPFPair dp = shape.tio.rdpfpair(shape.addr_size);
        RegAS p0indoffset, p1indoffset;

        // Receive the index offset from the computational players and
        // combine them
        shape.tio.recv_p0(&p0indoffset, BITBYTES(shape.addr_size));
        shape.tio.recv_p1(&p1indoffset, BITBYTES(shape.addr_size));
        auto indshift = combine(p0indoffset, p1indoffset, shape.addr_size);

        // Evaluate the DPFs to compute the cancellation terms
        T gamma0, gamma1;
        StreamEval ev(dp, indshift, 0, shape.tio.aes_ops());
        for (size_t i=0; i<shape.shape_size; ++i) {
            auto L = ev.next();

            // The values from the two DPFs
            auto [V0, V1] = dp.unit<T>(L);

            // shape.get_server(i) returns a pair of references to the
            // appropriate cells in the two blinded databases
            auto [BL0, BL1] = shape.get_server(i);
            gamma0 -= BL0 * V1.share();
            gamma1 -= BL1 * V0.share();
        }

        // Choose a random blinding factor
        T rho;
        rho.randomize();

        gamma0 += rho;
        gamma1 -= rho;

        // Send the cancellation terms to the computational players
        shape.tio.iostream_p0() << gamma0;
        shape.tio.iostream_p1() << gamma1;

        shape.yield();
    }
    return res;  // The server will always get 0
}

// Oblivious update to an additively shared index of Duoram memory
template <typename T>
typename Duoram<T>::Shape::MemRefAS
    &Duoram<T>::Shape::MemRefAS::operator+=(const T& M)
{
    int player = shape.tio.player();
    if (player < 2) {
        // Computational players do this

        RDPFTriple dt = shape.tio.rdpftriple(shape.addr_size);

        // Compute the index and message offsets
        RegAS indoffset = dt.as_target;
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
        RegAS peerindoffset;
        std::tuple<T,T,T> peerMoffset;
        shape.tio.recv_peer(&peerindoffset, BITBYTES(shape.addr_size));
        shape.tio.iostream_peer() >> peerMoffset;

        // Reconstruct the total offsets
        auto indshift = combine(indoffset, peerindoffset, shape.addr_size);
        auto Mshift = combine(Moffset, peerMoffset);

        // Evaluate the DPFs and add them to the database
        StreamEval ev(dt, indshift, 0, shape.tio.aes_ops());
        for (size_t i=0; i<shape.shape_size; ++i) {
            auto L = ev.next();
            // The values from the three DPFs
            auto [V0, V1, V2] = dt.scaled<T>(L) + dt.unit<T>(L) * Mshift;
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
        }
    } else {
        // The server does this

        RDPFPair dp = shape.tio.rdpfpair(shape.addr_size);
        RegAS p0indoffset, p1indoffset;
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
        StreamEval ev(dp, indshift, 0, shape.tio.aes_ops());
        for (size_t i=0; i<shape.shape_size; ++i) {
            auto L = ev.next();
            // The values from the two DPFs
            auto V = dp.scaled<T>(L) + dp.unit<T>(L) * Mshift;
            // shape.get_server(i) returns a pair of references to the
            // appropriate cells in the two blinded databases, so we can
            // subtract the pair directly.
            shape.get_server(i) -= V;
        }
    }
    return *this;
}

// Explicit read from a given index of Duoram memory
template <typename T>
Duoram<T>::Shape::MemRefExpl::operator T()
{
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
    int player = shape.tio.player();
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
