#ifndef __MCPIO_HPP__
#define __MCPIO_HPP__

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <queue>
#include <string>
#include <atomic>
#include <optional>
#include <bsd/stdlib.h> // arc4random_buf

#include <boost/asio.hpp>
#include <boost/thread.hpp>
#include <boost/chrono.hpp>

#include "types.hpp"

using boost::asio::ip::tcp;

// Classes to represent stored precomputed data (e.g., multiplication triples)

template<typename T>
class PreCompStorage {
public:
    PreCompStorage(unsigned player, bool preprocessing,
        const char *filenameprefix, unsigned thread_num);
    void get(T& nextval);

    inline size_t get_stats() { return count; }
    inline void reset_stats() { count = 0; }
private:
    std::ifstream storage;
    size_t count;
};

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

// If we want to send Lamport clocks in messages, define this.  It adds
// an 8-byte header to each message (length and Lamport clock), so it
// has a small network cost.  We always define and pass the Lamport
// clock member of MPCIO to the IO functions for simplicity, but they're
// ignored if this isn't defined
#define SEND_LAMPORT_CLOCKS
using lamport_t = uint32_t;
using atomic_lamport_t = std::atomic<lamport_t>;
using opt_lamport_t = std::optional<lamport_t>;
#ifdef SEND_LAMPORT_CLOCKS
struct MessageWithHeader {
    std::string header;
    std::string message;

    MessageWithHeader(std::string &&msg, lamport_t lamport) :
        message(std::move(msg)) {
            char hdr[sizeof(uint32_t) + sizeof(lamport_t)];
            uint32_t msglen = uint32_t(message.size());
            memmove(hdr, &msglen, sizeof(msglen));
            memmove(hdr+sizeof(msglen), &lamport, sizeof(lamport));
            header.assign(hdr, sizeof(hdr));
    }
};
#endif

// A class to wrap a socket to another MPC party.  This wrapping allows
// us to do some useful logging, and perform async_writes transparently
// to the application.

class MPCSingleIO {
    tcp::socket sock;
    size_t totread, totwritten;
#ifdef RECORD_IOTRACE
    std::vector<ssize_t> iotrace;
#endif

    // To avoid blocking if both we and our peer are trying to send
    // something very large, and neither side is receiving, we will send
    // with async_write.  But this has a number of implications:
    // - The data to be sent has to be copied into this MPCSingleIO,
    //   since asio::buffer pointers are not guaranteed to remain valid
    //   after the end of the coroutine that created them
    // - We have to keep a queue of messages to be sent, in case
    //   coroutines call send() before the previous message has finished
    //   being sent
    // - This queue may be accessed from the async_write thread as well
    //   as the work thread that uses this MPCSingleIO directly (there
    //   should be only one of the latter), so we need some locking

    // This is where we accumulate data passed in queue()
    std::string dataqueue;

    // When send() is called, the above dataqueue is appended to this
    // messagequeue, and the dataqueue is reset.  If messagequeue was
    // empty before this append, launch async_write to write the first
    // thing in the messagequeue.  When async_write completes, it will
    // delete the first thing in the messagequeue, and see if there are
    // any more elements.  If so, it will start another async_write.
    // The invariant is that there is an async_write currently running
    // iff messagequeue is nonempty.
#ifdef SEND_LAMPORT_CLOCKS
    std::queue<MessageWithHeader> messagequeue;
#else
    std::queue<std::string> messagequeue;
#endif

    // If a single message is broken into chunks in order to get the
    // first part of it out on the wire while the rest of it is still
    // being computed, we want the Lamport clock of all the chunks to be
    // that of when the message is first created.  This value will be
    // nullopt when there has been no queue() since the last explicit
    // send() (as opposed to the implicit send() called by queue()
    // itself if it wants to get a chunk on its way), and will be set to
    // the current lamport clock when that first queue() after each
    // explicit send() happens.
    opt_lamport_t message_lamport;

#ifdef SEND_LAMPORT_CLOCKS
    // If Lamport clocks are being sent, then the data stream is divided
    // into chunks, each with a header containing the length of the
    // chunk and the Lamport clock.  So when we read, we'll read a whole
    // chunk, and store it here.  Then calls to recv() will read pieces
    // of this buffer until it has all been read, and then read the next
    // header and chunk.
    std::string recvdata;
    size_t recvdataremain;
#endif

    // Never touch the above messagequeue without holding this lock (you
    // _can_ touch the strings it contains, though, if you looked one up
    // while holding the lock).
    boost::mutex messagequeuelock;

    // Asynchronously send the first message from the message queue.
    // * The messagequeuelock must be held when this is called! *
    // This method may be called from either thread (the work thread or
    // the async_write handler thread).
    void async_send_from_msgqueue() {
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

public:
    MPCSingleIO(tcp::socket &&sock) :
        sock(std::move(sock)), totread(0), totwritten(0) {}

    // Returns 1 if a new message is started, 0 otherwise
    size_t queue(const void *data, size_t len, lamport_t lamport) {
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

    void send(bool implicit_send = false) {
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

    size_t recv(void *data, size_t len, lamport_t &lamport) {
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
        return res;
#else
        size_t res = boost::asio::read(sock, boost::asio::buffer(data, len));
#ifdef RECORD_IOTRACE
        iotrace.push_back(-(ssize_t(res)));
#endif
        return res;
#endif
    }

#ifdef RECORD_IOTRACE
    void dumptrace(std::ostream &os, const char *label = NULL) {
        if (label) {
            os << label << " ";
        }
        os << "IO trace:";
        for (auto& s: iotrace) {
            os << " " << s;
        }
        os << "\n";
    }

    void resettrace() {
        iotrace.clear();
    }
#endif
};

// A base class to represent all of a computation peer or server's IO,
// either to other parties or to local storage (the computation and
// server cases are separate subclasses below).

struct MPCIO {
    int player;
    bool preprocessing;
    size_t num_threads;
    atomic_lamport_t lamport;
    std::vector<size_t> msgs_sent;
    std::vector<size_t> msg_bytes_sent;
    boost::chrono::steady_clock::time_point steady_start;
    boost::chrono::process_cpu_clock::time_point cpu_start;

    MPCIO(int player, bool preprocessing, size_t num_threads) :
        player(player), preprocessing(preprocessing),
        num_threads(num_threads), lamport(0)
    {
        reset_stats();
    }

    void reset_stats()
    {
        msgs_sent.clear();
        msg_bytes_sent.clear();
        for (size_t i=0; i<num_threads; ++i) {
            msgs_sent.push_back(0);
            msg_bytes_sent.push_back(0);
        }
        steady_start = boost::chrono::steady_clock::now();
        cpu_start = boost::chrono::process_cpu_clock::now();
    }

    void dump_stats(std::ostream &os)
    {
        size_t tot_msgs_sent = 0;
        size_t tot_msg_bytes_sent = 0;
        for (auto& n : msgs_sent) {
            tot_msgs_sent += n;
        }
        for (auto& n : msg_bytes_sent) {
            tot_msg_bytes_sent += n;
        }
        auto steady_elapsed =
            boost::chrono::steady_clock::now() - steady_start;
        auto cpu_elapsed =
            boost::chrono::process_cpu_clock::now() - cpu_start;

        os << tot_msgs_sent << " messages sent\n";
        os << tot_msg_bytes_sent << " message bytes sent\n";
        os << lamport << " Lamport clock (latencies)\n";
        os << boost::chrono::duration_cast
            <boost::chrono::milliseconds>(steady_elapsed) <<
            " wall clock time\n";
        os << cpu_elapsed << " {real;user;system}\n";
    }
};

// A class to represent all of a computation peer's IO, either to other
// parties or to local storage

struct MPCPeerIO : public MPCIO {
    // We use a deque here instead of a vector because you can't have a
    // vector of a type without a copy constructor (tcp::socket is the
    // culprit), but you can have a deque of those for some reason.
    std::deque<MPCSingleIO> peerios;
    std::deque<MPCSingleIO> serverios;
    std::vector<PreCompStorage<MultTriple>> triples;
    std::vector<PreCompStorage<HalfTriple>> halftriples;

    MPCPeerIO(unsigned player, bool preprocessing,
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

    void dump_precomp_stats(std::ostream &os)
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

    void reset_precomp_stats()
    {
        for (size_t i=0; i<triples.size(); ++i) {
            triples[i].reset_stats();
            halftriples[i].reset_stats();
        }
    }
};

// A class to represent all of the server party's IO, either to
// computational parties or to local storage

struct MPCServerIO : public MPCIO {
    std::deque<MPCSingleIO> p0ios;
    std::deque<MPCSingleIO> p1ios;

    MPCServerIO(bool preprocessing,
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
};

// A handle to one thread's sockets and streams in a MPCIO

class MPCTIO {
    int thread_num;
    lamport_t thread_lamport;
    MPCIO &mpcio;

public:
    MPCTIO(MPCIO &mpcio, int thread_num):
        thread_num(thread_num), thread_lamport(mpcio.lamport),
        mpcio(mpcio) {}

    // Sync our per-thread lamport clock with the master one in the
    // mpcio.  You only need to call this explicitly if your MPCTIO
    // outlives your thread (in which case call it after the join), or
    // if your threads do interthread communication amongst themselves
    // (in which case call it in the sending thread before the send, and
    // call it in the receiving thread after the receive).
    void sync_lamport() {
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

    // The normal case, where the MPCIO is created inside the thread,
    // and so destructed when the thread ends, is handles automatically
    // here.
    ~MPCTIO() {
        sync_lamport();
    }

    // Queue up data to the peer or to the server

    void queue_peer(const void *data, size_t len) {
        if (mpcio.player < 2) {
            MPCPeerIO &mpcpio = static_cast<MPCPeerIO&>(mpcio);
            size_t newmsg = mpcpio.peerios[thread_num].queue(data, len, thread_lamport);
            mpcpio.msgs_sent[thread_num] += newmsg;
            mpcpio.msg_bytes_sent[thread_num] += len;
        }
    }

    void queue_server(const void *data, size_t len) {
        if (mpcio.player < 2) {
            MPCPeerIO &mpcpio = static_cast<MPCPeerIO&>(mpcio);
            size_t newmsg = mpcpio.serverios[thread_num].queue(data, len, thread_lamport);
            mpcpio.msgs_sent[thread_num] += newmsg;
            mpcpio.msg_bytes_sent[thread_num] += len;
        }
    }

    // Receive data from the peer or to the server

    size_t recv_peer(void *data, size_t len) {
        if (mpcio.player < 2) {
            MPCPeerIO &mpcpio = static_cast<MPCPeerIO&>(mpcio);
            return mpcpio.peerios[thread_num].recv(data, len, thread_lamport);
        }
        return 0;
    }

    size_t recv_server(void *data, size_t len) {
        if (mpcio.player < 2) {
            MPCPeerIO &mpcpio = static_cast<MPCPeerIO&>(mpcio);
            return mpcpio.serverios[thread_num].recv(data, len, thread_lamport);
        }
        return 0;
    }

    // Queue up data to p0 or p1

    void queue_p0(const void *data, size_t len) {
        if (mpcio.player == 2) {
            MPCServerIO &mpcsrvio = static_cast<MPCServerIO&>(mpcio);
            size_t newmsg = mpcsrvio.p0ios[thread_num].queue(data, len, thread_lamport);
            mpcsrvio.msgs_sent[thread_num] += newmsg;
            mpcsrvio.msg_bytes_sent[thread_num] += len;
        }
    }

    void queue_p1(const void *data, size_t len) {
        if (mpcio.player == 2) {
            MPCServerIO &mpcsrvio = static_cast<MPCServerIO&>(mpcio);
            size_t newmsg = mpcsrvio.p1ios[thread_num].queue(data, len, thread_lamport);
            mpcsrvio.msgs_sent[thread_num] += newmsg;
            mpcsrvio.msg_bytes_sent[thread_num] += len;
        }
    }

    // Receive data from p0 or p1

    size_t recv_p0(void *data, size_t len) {
        if (mpcio.player == 2) {
            MPCServerIO &mpcsrvio = static_cast<MPCServerIO&>(mpcio);
            return mpcsrvio.p0ios[thread_num].recv(data, len, thread_lamport);
        }
        return 0;
    }

    size_t recv_p1(void *data, size_t len) {
        if (mpcio.player == 2) {
            MPCServerIO &mpcsrvio = static_cast<MPCServerIO&>(mpcio);
            return mpcsrvio.p1ios[thread_num].recv(data, len, thread_lamport);
        }
        return 0;
    }

    // Send all queued data for this thread
    void send() {
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
    MultTriple triple() {
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

    HalfTriple halftriple() {
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

    // Accessors
    inline int player() { return mpcio.player; }
    inline bool preprocessing() { return mpcio.preprocessing; }
    inline bool is_server() { return mpcio.player == 2; }
};

// Set up the socket connections between the two computational parties
// (P0 and P1) and the server party (P2).  For each connection, the
// lower-numbered party does the accept() and the higher-numbered party
// does the connect().

// Computational parties call this version with player=0 or 1

void mpcio_setup_computational(unsigned player,
    boost::asio::io_context &io_context,
    const char *p0addr,  // can be NULL when player=0
    int num_threads,
    std::deque<tcp::socket> &peersocks,
    std::deque<tcp::socket> &serversocks);

// Server calls this version

void mpcio_setup_server(boost::asio::io_context &io_context,
    const char *p0addr, const char *p1addr, int num_threads,
    std::deque<tcp::socket> &p0socks,
    std::deque<tcp::socket> &p1socks);

#endif
