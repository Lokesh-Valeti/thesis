#ifndef __MCPIO_HPP__
#define __MCPIO_HPP__

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <queue>
#include <string>

#include <boost/asio.hpp>
#include <boost/coroutine2/all.hpp>
#include <boost/thread.hpp>

#include "types.hpp"

using boost::asio::ip::tcp;

// Classes to represent stored precomputed data (e.g., multiplication triples)

template<typename T>
class PreCompStorage {
public:
    PreCompStorage(unsigned player, bool preprocessing,
        const char *filenameprefix, unsigned thread_num);
    void get(T& nextval);
private:
    std::ifstream storage;
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
}

template<typename T>
void PreCompStorage<T>::get(T& nextval) {
    storage.read((char *)&nextval, sizeof(T));
    if (storage.gcount() != sizeof(T)) {
        std::cerr << "Failed to read precomputed value from storage\n";
        exit(1);
    }
}

// A class to wrap a socket to another MPC party.  This wrapping allows
// us to do some useful logging, and perform async_writes transparently
// to the application.

class MPCSingleIO {
    tcp::socket sock;
    size_t totread, totwritten;
    std::vector<ssize_t> iotrace;

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
    std::queue<std::string> messagequeue;

    // Never touch the above messagequeue without holding this lock (you
    // _can_ touch the strings it contains, though, if you looked one up
    // while holding the lock).
    boost::mutex messagequeuelock;

    // Asynchronously send the first message from the message queue.
    // * The messagequeuelock must be held when this is called! *
    // This method may be called from either thread (the work thread or
    // the async_write handler thread).
    void async_send_from_msgqueue() {
        boost::asio::async_write(sock,
            boost::asio::buffer(messagequeue.front()),
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

    void queue(const void *data, size_t len) {
        dataqueue.append((const char *)data, len);

        // If we already have some full packets worth of data, may as
        // well send it.
        if (dataqueue.size() > 28800) {
            send();
        }
    }

    void send() {
        size_t thissize = dataqueue.size();
        // Ignore spurious calls to send()
        if (thissize == 0) return;

        iotrace.push_back(thissize);

        messagequeuelock.lock();
        // Move the current message to send into the message queue (this
        // moves a pointer to the data, not copying the data itself)
        messagequeue.emplace(std::move(dataqueue));
        // If this is now the first thing in the message queue, launch
        // an async_write to write it
        if (messagequeue.size() == 1) {
            async_send_from_msgqueue();
        }
        messagequeuelock.unlock();
    }

    size_t recv(const std::vector<boost::asio::mutable_buffer>& buffers) {
        size_t res = boost::asio::read(sock, buffers);
        iotrace.push_back(-(ssize_t(res)));
        return res;
    }

    size_t recv(const boost::asio::mutable_buffer& buffer) {
        size_t res = boost::asio::read(sock, buffer);
        iotrace.push_back(-(ssize_t(res)));
        return res;
    }

    size_t recv(void *data, size_t len) {
        size_t res = boost::asio::read(sock, boost::asio::buffer(data, len));
        iotrace.push_back(-(ssize_t(res)));
        return res;
    }

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
};

// A class to represent all of a computation party's IO, either to other
// parties or to local storage

struct MPCIO {
    int player;
    // We use a deque here instead of a vector because you can't have a
    // vector of a type without a copy constructor (tcp::socket is the
    // culprit), but you can have a deque of those for some reason.
    std::deque<MPCSingleIO> peerios;
    MPCSingleIO serverio;
    std::vector<PreCompStorage<MultTriple>> triples;
    std::vector<PreCompStorage<HalfTriple>> halftriples;

    MPCIO(unsigned player, bool preprocessing,
            std::deque<tcp::socket> &peersocks, tcp::socket &&serversock) :
        player(player), serverio(std::move(serversock)) {
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
    }
};

// A class to represent all of the server party's IO, either to
// computational parties or to local storage

struct MPCServerIO {
    MPCSingleIO p0io;
    MPCSingleIO p1io;

    MPCServerIO(bool preprocessing, tcp::socket &&p0sock,
            tcp::socket &&p1sock) :
        p0io(std::move(p0sock)), p1io(std::move(p1sock)) {}
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
    std::deque<tcp::socket> &peersocks, tcp::socket &serversock);

// Server calls this version with player=2

void mpcio_setup_server(boost::asio::io_context &io_context,
    const char *p0addr, const char *p1addr,
    tcp::socket &p0sock, tcp::socket &p1sock);

#endif
