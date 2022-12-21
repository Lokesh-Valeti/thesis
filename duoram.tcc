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
    this->shape_size = len;
}

template <typename T>
Duoram<T>::MemRefAS::operator T()
{
    T res;
    return res;
}
