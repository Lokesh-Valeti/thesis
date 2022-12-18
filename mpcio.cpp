#include "mpcio.hpp"

template<typename T>
PreCompStorage<T>::PreCompStorage(unsigned player, bool preprocessing,
        const char *filenameprefix, unsigned thread_num) {
    if (preprocessing) return;
    std::string filename(filenameprefix);
    char suffix[20];
    sprintf(suffix, ".p%d.t%u", player%10, thread_num);
    filename.append(suffix);
    storage.open(filename);
    if (storage.fail()) {
        std::cerr << "Failed to open " << filename << "\n";
        exit(1);
    }
    count = 0;
}

template<typename T>
void PreCompStorage<T>::get(T& nextval) {
    storage.read((char *)&nextval, sizeof(T));
    if (storage.gcount() != sizeof(T)) {
        std::cerr << "Failed to read precomputed value from storage\n";
        exit(1);
    }
    ++count;
}

void MPCSingleIO::async_send_from_msgqueue()
{
#ifdef SEND_LAMPORT_CLOCKS
    std::vector<boost::asio::const_buffer> tosend;
    tosend.push_back(boost::asio::buffer(messagequeue.front().header));
    tosend.push_back(boost::asio::buffer(messagequeue.front().message));
#endif
    boost::asio::async_write(sock,
#ifdef SEND_LAMPORT_CLOCKS
        tosend,
#else
        boost::asio::buffer(messagequeue.front()),
#endif
        [&](boost::system::error_code ec, std::size_t amt){
            messagequeuelock.lock();
            messagequeue.pop();
            if (messagequeue.size() > 0) {
                async_send_from_msgqueue();
            }
            messagequeuelock.unlock();
        });
}

size_t MPCSingleIO::queue(const void *data, size_t len, lamport_t lamport)
{
    // Is this a new message?
    size_t newmsg = 0;

    dataqueue.append((const char *)data, len);

    // If this is the first queue() since the last explicit send(),
    // which we'll know because message_lamport will be nullopt, set
    // message_lamport to the current Lamport clock.  Note that the
    // boolean test tests whether message_lamport is nullopt, not
    // whether its value is zero.
    if (!message_lamport) {
        message_lamport = lamport;
        newmsg = 1;
    }

    // If we already have some full packets worth of data, may as
    // well send it.
    if (dataqueue.size() > 28800) {
        send(true);
    }

    return newmsg;
}

void MPCSingleIO::send(bool implicit_send)
{
    size_t thissize = dataqueue.size();
    // Ignore spurious calls to send(), except for resetting
    // message_lamport if this was an explicit send().
    if (thissize == 0) {
#ifdef SEND_LAMPORT_CLOCKS
        // If this was an explicit send(), reset the message_lamport so
        // that it gets updated at the next queue().
        if (!implicit_send) {
            message_lamport.reset();
        }
#endif
        return;
    }

#ifdef RECORD_IOTRACE
    iotrace.push_back(thissize);
#endif

    messagequeuelock.lock();
    // Move the current message to send into the message queue (this
    // moves a pointer to the data, not copying the data itself)
#ifdef SEND_LAMPORT_CLOCKS
    messagequeue.emplace(std::move(dataqueue),
        message_lamport.value());
    // If this was an explicit send(), reset the message_lamport so
    // that it gets updated at the next queue().
    if (!implicit_send) {
        message_lamport.reset();
    }
#else
    messagequeue.emplace(std::move(dataqueue));
#endif
    // If this is now the first thing in the message queue, launch
    // an async_write to write it
    if (messagequeue.size() == 1) {
        async_send_from_msgqueue();
    }
    messagequeuelock.unlock();
}

size_t MPCSingleIO::recv(void *data, size_t len, lamport_t &lamport)
{
#ifdef SEND_LAMPORT_CLOCKS
    char *cdata = (char *)data;
    size_t res = 0;
    while (len > 0) {
        while (recvdataremain == 0) {
            // Read a new header
            char hdr[sizeof(uint32_t) + sizeof(lamport_t)];
            uint32_t datalen;
            lamport_t recv_lamport;
            boost::asio::read(sock, boost::asio::buffer(hdr, sizeof(hdr)));
            memmove(&datalen, hdr, sizeof(datalen));
            memmove(&recv_lamport, hdr+sizeof(datalen), sizeof(lamport_t));
            lamport_t new_lamport = recv_lamport + 1;
            if (lamport < new_lamport) {
                lamport = new_lamport;
            }
            if (datalen > 0) {
                recvdata.resize(datalen, '\0');
                boost::asio::read(sock, boost::asio::buffer(recvdata));
                recvdataremain = datalen;
            }
        }
        size_t amttoread = len;
        if (amttoread > recvdataremain) {
            amttoread = recvdataremain;
        }
        memmove(cdata, recvdata.data()+recvdata.size()-recvdataremain,
            amttoread);
        cdata += amttoread;
        len -= amttoread;
        recvdataremain -= amttoread;
        res += amttoread;
    }
#else
    size_t res = boost::asio::read(sock, boost::asio::buffer(data, len));
#endif
#ifdef RECORD_IOTRACE
    iotrace.push_back(-(ssize_t(res)));
#endif
    return res;
}

#ifdef RECORD_IOTRACE
void MPCSingleIO::dumptrace(std::ostream &os, const char *label)
{
    if (label) {
        os << label << " ";
    }
    os << "IO trace:";
    for (auto& s: iotrace) {
        os << " " << s;
    }
    os << "\n";
}
#endif

void MPCIO::reset_stats()
{
    msgs_sent.clear();
    msg_bytes_sent.clear();
    aes_ops.clear();
    for (size_t i=0; i<num_threads; ++i) {
        msgs_sent.push_back(0);
        msg_bytes_sent.push_back(0);
        aes_ops.push_back(0);
    }
    steady_start = boost::chrono::steady_clock::now();
    cpu_start = boost::chrono::process_cpu_clock::now();
}

void MPCIO::dump_stats(std::ostream &os)
{
    size_t tot_msgs_sent = 0;
    size_t tot_msg_bytes_sent = 0;
    size_t tot_aes_ops = 0;
    for (auto& n : msgs_sent) {
        tot_msgs_sent += n;
    }
    for (auto& n : msg_bytes_sent) {
        tot_msg_bytes_sent += n;
    }
    for (auto& n : aes_ops) {
        tot_aes_ops += n;
    }
    auto steady_elapsed =
        boost::chrono::steady_clock::now() - steady_start;
    auto cpu_elapsed =
        boost::chrono::process_cpu_clock::now() - cpu_start;

    os << tot_msgs_sent << " messages sent\n";
    os << tot_msg_bytes_sent << " message bytes sent\n";
    os << tot_aes_ops << " local AES operations\n";
    os << lamport << " Lamport clock (latencies)\n";
    os << boost::chrono::duration_cast
        <boost::chrono::milliseconds>(steady_elapsed) <<
        " wall clock time\n";
    os << cpu_elapsed << " {real;user;system}\n";
}

MPCPeerIO::MPCPeerIO(unsigned player, bool preprocessing,
        std::deque<tcp::socket> &peersocks,
        std::deque<tcp::socket> &serversocks) :
    MPCIO(player, preprocessing, peersocks.size())
{
    unsigned num_threads = unsigned(peersocks.size());
    for (unsigned i=0; i<num_threads; ++i) {
        triples.emplace_back(player, preprocessing, "triples", i);
    }
    for (unsigned i=0; i<num_threads; ++i) {
        halftriples.emplace_back(player, preprocessing, "halves", i);
    }
    for (auto &&sock : peersocks) {
        peerios.emplace_back(std::move(sock));
    }
    for (auto &&sock : serversocks) {
        serverios.emplace_back(std::move(sock));
    }
}

void MPCPeerIO::dump_precomp_stats(std::ostream &os)
{
    for (size_t i=0; i<triples.size(); ++i) {
        if (i > 0) {
            os << " ";
        }
        os << "T" << i << " t:" << triples[i].get_stats() <<
            " h:" << halftriples[i].get_stats();
    }
    os << "\n";
}

void MPCPeerIO::reset_precomp_stats()
{
    for (size_t i=0; i<triples.size(); ++i) {
        triples[i].reset_stats();
        halftriples[i].reset_stats();
    }
}

void MPCPeerIO::dump_stats(std::ostream &os)
{
    MPCIO::dump_stats(os);
    os << "Precomputed values used: ";
    dump_precomp_stats(os);
}

MPCServerIO::MPCServerIO(bool preprocessing,
        std::deque<tcp::socket> &p0socks,
        std::deque<tcp::socket> &p1socks) :
    MPCIO(2, preprocessing, p0socks.size())
{
    for (auto &&sock : p0socks) {
        p0ios.emplace_back(std::move(sock));
    }
    for (auto &&sock : p1socks) {
        p1ios.emplace_back(std::move(sock));
    }
}

// Sync our per-thread lamport clock with the master one in the
// mpcio.  You only need to call this explicitly if your MPCTIO
// outlives your thread (in which case call it after the join), or
// if your threads do interthread communication amongst themselves
// (in which case call it in the sending thread before the send, and
// call it in the receiving thread after the receive).
void MPCTIO::sync_lamport()
{
    // Update the mpcio Lamport time to be max of the thread Lamport
    // time and what we thought it was before.  We use this
    // compare_exchange construction in order to atomically
    // do the comparison, computation, and replacement
    lamport_t old_lamport = mpcio.lamport;
    lamport_t new_lamport = thread_lamport;
    do {
        if (new_lamport < old_lamport) {
            new_lamport = old_lamport;
        }
    // The next line atomically checks if lamport still has
    // the value old_lamport; if so, it changes its value to
    // new_lamport and returns true (ending the loop).  If
    // not, it sets old_lamport to the current value of
    // lamport, and returns false (continuing the loop so
    // that new_lamport can be recomputed based on this new
    // value).
    } while (!mpcio.lamport.compare_exchange_weak(
        old_lamport, new_lamport));
    thread_lamport = new_lamport;
}

// Queue up data to the peer or to the server

void MPCTIO::queue_peer(const void *data, size_t len)
{
    if (mpcio.player < 2) {
        MPCPeerIO &mpcpio = static_cast<MPCPeerIO&>(mpcio);
        size_t newmsg = mpcpio.peerios[thread_num].queue(data, len, thread_lamport);
        mpcpio.msgs_sent[thread_num] += newmsg;
        mpcpio.msg_bytes_sent[thread_num] += len;
    }
}

void MPCTIO::queue_server(const void *data, size_t len)
{
    if (mpcio.player < 2) {
        MPCPeerIO &mpcpio = static_cast<MPCPeerIO&>(mpcio);
        size_t newmsg = mpcpio.serverios[thread_num].queue(data, len, thread_lamport);
        mpcpio.msgs_sent[thread_num] += newmsg;
        mpcpio.msg_bytes_sent[thread_num] += len;
    }
}

// Receive data from the peer or to the server

size_t MPCTIO::recv_peer(void *data, size_t len)
{
    if (mpcio.player < 2) {
        MPCPeerIO &mpcpio = static_cast<MPCPeerIO&>(mpcio);
        return mpcpio.peerios[thread_num].recv(data, len, thread_lamport);
    }
    return 0;
}

size_t MPCTIO::recv_server(void *data, size_t len)
{
    if (mpcio.player < 2) {
        MPCPeerIO &mpcpio = static_cast<MPCPeerIO&>(mpcio);
        return mpcpio.serverios[thread_num].recv(data, len, thread_lamport);
    }
    return 0;
}

// Queue up data to p0 or p1

void MPCTIO::queue_p0(const void *data, size_t len)
{
    if (mpcio.player == 2) {
        MPCServerIO &mpcsrvio = static_cast<MPCServerIO&>(mpcio);
        size_t newmsg = mpcsrvio.p0ios[thread_num].queue(data, len, thread_lamport);
        mpcsrvio.msgs_sent[thread_num] += newmsg;
        mpcsrvio.msg_bytes_sent[thread_num] += len;
    }
}

void MPCTIO::queue_p1(const void *data, size_t len)
{
    if (mpcio.player == 2) {
        MPCServerIO &mpcsrvio = static_cast<MPCServerIO&>(mpcio);
        size_t newmsg = mpcsrvio.p1ios[thread_num].queue(data, len, thread_lamport);
        mpcsrvio.msgs_sent[thread_num] += newmsg;
        mpcsrvio.msg_bytes_sent[thread_num] += len;
    }
}

// Receive data from p0 or p1

size_t MPCTIO::recv_p0(void *data, size_t len)
{
    if (mpcio.player == 2) {
        MPCServerIO &mpcsrvio = static_cast<MPCServerIO&>(mpcio);
        return mpcsrvio.p0ios[thread_num].recv(data, len, thread_lamport);
    }
    return 0;
}

size_t MPCTIO::recv_p1(void *data, size_t len)
{
    if (mpcio.player == 2) {
        MPCServerIO &mpcsrvio = static_cast<MPCServerIO&>(mpcio);
        return mpcsrvio.p1ios[thread_num].recv(data, len, thread_lamport);
    }
    return 0;
}

// Send all queued data for this thread

void MPCTIO::send()
{
    if (mpcio.player < 2) {
        MPCPeerIO &mpcpio = static_cast<MPCPeerIO&>(mpcio);
        mpcpio.peerios[thread_num].send();
        mpcpio.serverios[thread_num].send();
    } else {
        MPCServerIO &mpcsrvio = static_cast<MPCServerIO&>(mpcio);
        mpcsrvio.p0ios[thread_num].send();
        mpcsrvio.p1ios[thread_num].send();
    }
}

// Functions to get precomputed values.  If we're in the online
// phase, get them from PreCompStorage.  If we're in the
// preprocessing phase, read them from the server.
MultTriple MPCTIO::triple()
{
    MultTriple val;
    if (mpcio.player < 2) {
        MPCPeerIO &mpcpio = static_cast<MPCPeerIO&>(mpcio);
        if (mpcpio.preprocessing) {
            recv_server(&val, sizeof(val));
        } else {
            mpcpio.triples[thread_num].get(val);
        }
    } else if (mpcio.preprocessing) {
        // Create triples (X0,Y0,Z0),(X1,Y1,Z1) such that
        // (X0*Y1 + Y0*X1) = (Z0+Z1)
        value_t X0, Y0, Z0, X1, Y1, Z1;
        arc4random_buf(&X0, sizeof(X0));
        arc4random_buf(&Y0, sizeof(Y0));
        arc4random_buf(&Z0, sizeof(Z0));
        arc4random_buf(&X1, sizeof(X1));
        arc4random_buf(&Y1, sizeof(Y1));
        Z1 = X0 * Y1 + X1 * Y0 - Z0;
        MultTriple T0, T1;
        T0 = std::make_tuple(X0, Y0, Z0);
        T1 = std::make_tuple(X1, Y1, Z1);
        queue_p0(&T0, sizeof(T0));
        queue_p1(&T1, sizeof(T1));
    }
    return val;
}

HalfTriple MPCTIO::halftriple()
{
    HalfTriple val;
    if (mpcio.player < 2) {
        MPCPeerIO &mpcpio = static_cast<MPCPeerIO&>(mpcio);
        if (mpcpio.preprocessing) {
            recv_server(&val, sizeof(val));
        } else {
            mpcpio.halftriples[thread_num].get(val);
        }
    } else if (mpcio.preprocessing) {
        // Create half-triples (X0,Z0),(Y1,Z1) such that
        // X0*Y1 = Z0 + Z1
        value_t X0, Z0, Y1, Z1;
        arc4random_buf(&X0, sizeof(X0));
        arc4random_buf(&Z0, sizeof(Z0));
        arc4random_buf(&Y1, sizeof(Y1));
        Z1 = X0 * Y1 - Z0;
        HalfTriple H0, H1;
        H0 = std::make_tuple(X0, Z0);
        H1 = std::make_tuple(Y1, Z1);
        queue_p0(&H0, sizeof(H0));
        queue_p1(&H1, sizeof(H1));
    }
    return val;
}

AndTriple MPCTIO::andtriple()
{
    AndTriple val;
    if (mpcio.player < 2) {
        MPCPeerIO &mpcpio = static_cast<MPCPeerIO&>(mpcio);
        if (mpcpio.preprocessing) {
            recv_server(&val, sizeof(val));
        } else {
            std::cerr << "Attempted to read AndTriple in online phase\n";
        }
    } else if (mpcio.preprocessing) {
        // Create triples (X0,Y0,Z0),(X1,Y1,Z1) such that
        // (X0&Y1 ^ Y0&X1) = (Z0^Z1)
        DPFnode X0, Y0, Z0, X1, Y1, Z1;
        arc4random_buf(&X0, sizeof(X0));
        arc4random_buf(&Y0, sizeof(Y0));
        arc4random_buf(&Z0, sizeof(Z0));
        arc4random_buf(&X1, sizeof(X1));
        arc4random_buf(&Y1, sizeof(Y1));
        Z1 = ((X0 & Y1) ^ (X1 & Y0)) ^ Z0;
        AndTriple T0, T1;
        T0 = std::make_tuple(X0, Y0, Z0);
        T1 = std::make_tuple(X1, Y1, Z1);
        queue_p0(&T0, sizeof(T0));
        queue_p1(&T1, sizeof(T1));
    }
    return val;
}

// The port number for the P1 -> P0 connection
static const unsigned short port_p1_p0 = 2115;

// The port number for the P2 -> P0 connection
static const unsigned short port_p2_p0 = 2116;

// The port number for the P2 -> P1 connection
static const unsigned short port_p2_p1 = 2117;

void mpcio_setup_computational(unsigned player,
    boost::asio::io_context &io_context,
    const char *p0addr,  // can be NULL when player=0
    int num_threads,
    std::deque<tcp::socket> &peersocks,
    std::deque<tcp::socket> &serversocks)
{
    if (player == 0) {
        // Listen for connections from P1 and from P2
        tcp::acceptor acceptor_p1(io_context,
            tcp::endpoint(tcp::v4(), port_p1_p0));
        tcp::acceptor acceptor_p2(io_context,
            tcp::endpoint(tcp::v4(), port_p2_p0));

        peersocks.clear();
        serversocks.clear();
        for (int i=0;i<num_threads;++i) {
            peersocks.emplace_back(io_context);
            serversocks.emplace_back(io_context);
        }
        for (int i=0;i<num_threads;++i) {
            tcp::socket peersock = acceptor_p1.accept();
            // Read 2 bytes from the socket, which will be the thread
            // number
            unsigned short thread_num;
            boost::asio::read(peersock,
                boost::asio::buffer(&thread_num, sizeof(thread_num)));
            if (thread_num >= num_threads) {
                std::cerr << "Received bad thread number from peer\n";
            } else {
                peersocks[thread_num] = std::move(peersock);
            }
        }
        for (int i=0;i<num_threads;++i) {
            tcp::socket serversock = acceptor_p2.accept();
            // Read 2 bytes from the socket, which will be the thread
            // number
            unsigned short thread_num;
            boost::asio::read(serversock,
                boost::asio::buffer(&thread_num, sizeof(thread_num)));
            if (thread_num >= num_threads) {
                std::cerr << "Received bad thread number from server\n";
            } else {
                serversocks[thread_num] = std::move(serversock);
            }
        }
    } else if (player == 1) {
        // Listen for connections from P2, make num_threads connections to P0
        tcp::acceptor acceptor_p2(io_context,
            tcp::endpoint(tcp::v4(), port_p2_p1));

        tcp::resolver resolver(io_context);
        boost::system::error_code err;
        peersocks.clear();
        serversocks.clear();
        for (int i=0;i<num_threads;++i) {
            serversocks.emplace_back(io_context);
        }
        for (unsigned short thread_num = 0; thread_num < num_threads; ++thread_num) {
            tcp::socket peersock(io_context);
            while(1) {
                boost::asio::connect(peersock,
                    resolver.resolve(p0addr, std::to_string(port_p1_p0)), err);
                if (!err) break;
                std::cerr << "Connection to p0 refused, will retry.\n";
                sleep(1);
            }
            // Write 2 bytes to the socket indicating which thread
            // number this socket is for
            boost::asio::write(peersock,
                boost::asio::buffer(&thread_num, sizeof(thread_num)));
            peersocks.push_back(std::move(peersock));
        }
        for (int i=0;i<num_threads;++i) {
            tcp::socket serversock = acceptor_p2.accept();
            // Read 2 bytes from the socket, which will be the thread
            // number
            unsigned short thread_num;
            boost::asio::read(serversock,
                boost::asio::buffer(&thread_num, sizeof(thread_num)));
            if (thread_num >= num_threads) {
                std::cerr << "Received bad thread number from server\n";
            } else {
                serversocks[thread_num] = std::move(serversock);
            }
        }
    } else {
        std::cerr << "Invalid player number passed to mpcio_setup_computational\n";
    }
}

void mpcio_setup_server(boost::asio::io_context &io_context,
    const char *p0addr, const char *p1addr, int num_threads,
    std::deque<tcp::socket> &p0socks,
    std::deque<tcp::socket> &p1socks)
{
    // Make connections to P0 and P1
    tcp::resolver resolver(io_context);
    boost::system::error_code err;
    p0socks.clear();
    p1socks.clear();
    for (unsigned short thread_num = 0; thread_num < num_threads; ++thread_num) {
        tcp::socket p0sock(io_context);
        while(1) {
            boost::asio::connect(p0sock,
                resolver.resolve(p0addr, std::to_string(port_p2_p0)), err);
            if (!err) break;
            std::cerr << "Connection to p0 refused, will retry.\n";
            sleep(1);
        }
        // Write 2 bytes to the socket indicating which thread
        // number this socket is for
        boost::asio::write(p0sock,
            boost::asio::buffer(&thread_num, sizeof(thread_num)));
        p0socks.push_back(std::move(p0sock));
    }
    for (unsigned short thread_num = 0; thread_num < num_threads; ++thread_num) {
        tcp::socket p1sock(io_context);
        while(1) {
            boost::asio::connect(p1sock,
                resolver.resolve(p1addr, std::to_string(port_p2_p1)), err);
            if (!err) break;
            std::cerr << "Connection to p1 refused, will retry.\n";
            sleep(1);
        }
        // Write 2 bytes to the socket indicating which thread
        // number this socket is for
        boost::asio::write(p1sock,
            boost::asio::buffer(&thread_num, sizeof(thread_num)));
        p1socks.push_back(std::move(p1sock));
    }
}
