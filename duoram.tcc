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
    return res;
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
        RegAS indoffset = idx;
        indoffset -= dt.as_target;
        auto Moffset = std::make_tuple(M, M, M);
        Moffset -= dt.scaled_sum();

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
        std::tuple<RegAS,RegAS,RegAS> peerMoffset;
        shape.tio.recv_peer(&peerindoffset, BITBYTES(shape.addr_size));
        shape.tio.iostream_peer() >> peerMoffset;

        // Reconstruct the total offsets
        auto indshift = combine(indoffset, peerindoffset, shape.addr_size);
        auto Mshift = combine(Moffset, peerMoffset);

        // Evaluate the DPFs and add them to the database
        StreamEval ev(dt, indshift, shape.tio.aes_ops());
        for (size_t i=0; i<shape.shape_size; ++i) {
            auto L = ev.next();
            shape.get_comp(i) += dt.scaled_as(L) + dt.unit_as(L) * Mshift;
        }
    } else {
        // The server does this
        RDPFPair dp = shape.tio.rdpfpair(shape.addr_size);
        RegAS p0indoffset, p1indoffset;
        RegAS p0Moffset[2];
        RegAS p1Moffset[2];
        shape.tio.recv_p0(&p0indoffset, BITBYTES(shape.addr_size));
        shape.tio.iostream_p0() >> p0Moffset[0] >> p0Moffset[1];
        shape.tio.recv_p1(&p1indoffset, BITBYTES(shape.addr_size));
        shape.tio.iostream_p1() >> p1Moffset[0] >> p1Moffset[1];
        p0indoffset += p1indoffset;
        p0Moffset[0] += p1Moffset[0];
        p0Moffset[1] += p1Moffset[1];
    }
    return *this;
}
