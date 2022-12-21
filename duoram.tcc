// Templated method implementations for duoram.hpp

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
        RegAS indoffset = idx;
        indoffset -= dt.as_target;
        RegAS Moffset[3];
        Moffset[0] = Moffset[1] = Moffset[2] = M;
        Moffset[0] -= dt.dpf[0].scaled_sum;
        Moffset[1] -= dt.dpf[1].scaled_sum;
        Moffset[2] -= dt.dpf[2].scaled_sum;
        shape.tio.queue_peer(&indoffset, BITBYTES(shape.addr_size));
        shape.tio.iostream_peer() << Moffset[0] << Moffset[1] <<
            Moffset[2];
        shape.tio.queue_server(&indoffset, BITBYTES(shape.addr_size));
        shape.tio.iostream_server() << Moffset[1] << Moffset[2];
        shape.yield();
        RegAS peerindoffset;
        RegAS peerMoffset[3];
        shape.tio.recv_peer(&peerindoffset, BITBYTES(shape.addr_size));
        shape.tio.iostream_peer() >> peerMoffset[0] >> peerMoffset[1] >>
            peerMoffset[2];
        indoffset += peerindoffset;
        indoffset &= shape.addr_mask;
        Moffset[0] += peerMoffset[0];
        Moffset[1] += peerMoffset[1];
        Moffset[2] += peerMoffset[2];
        StreamEval<RDPFTriple> ev(dt, indoffset.ashare,
            shape.tio.aes_ops());
        for (size_t i=0; i<shape.shape_size; ++i) {
            auto [L0, L1, L2] = ev.next();

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
