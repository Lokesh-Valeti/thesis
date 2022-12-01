#include "mpcio.hpp"

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
    std::deque<tcp::socket> &peersocks, tcp::socket &serversock)
{
    if (player == 0) {
        // Listen for connections from P1 and from P2
        tcp::acceptor acceptor_p1(io_context,
            tcp::endpoint(tcp::v4(), port_p1_p0));
        tcp::acceptor acceptor_p2(io_context,
            tcp::endpoint(tcp::v4(), port_p2_p0));

        for (int i=0;i<num_threads;++i) {
            peersocks.emplace_back(io_context);
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
        serversock = acceptor_p2.accept();
    } else if (player == 1) {
        // Listen for connections from P2, make num_threads connections to P0
        tcp::acceptor acceptor_p2(io_context,
            tcp::endpoint(tcp::v4(), port_p2_p1));

        tcp::resolver resolver(io_context);
        boost::system::error_code err;
        peersocks.clear();
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
        serversock = acceptor_p2.accept();
    } else {
        std::cerr << "Invalid player number passed to mpcio_setup_computational\n";
    }
}

void mpcio_setup_server(boost::asio::io_context &io_context,
    const char *p0addr, const char *p1addr,
    tcp::socket &p0sock, tcp::socket &p1sock)
{
    // Make connections to P0 and P1
    tcp::resolver resolver(io_context);
    boost::system::error_code err;
    while(1) {
        boost::asio::connect(p0sock,
            resolver.resolve(p0addr, std::to_string(port_p2_p0)), err);
        if (!err) break;
        std::cerr << "Connection to p0 refused, will retry.\n";
        sleep(1);
    }
    while(1) {
        boost::asio::connect(p1sock,
            resolver.resolve(p1addr, std::to_string(port_p2_p1)), err);
        if (!err) break;
        std::cerr << "Connection to p1 refused, will retry.\n";
        sleep(1);
    }
}
