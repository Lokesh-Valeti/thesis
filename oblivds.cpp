#include <iostream>

#include "mpcio.hpp"
#include "preproc.hpp"

static void usage(const char *progname)
{
    std::cerr << "Usage: " << progname << " [-p] player_num player_addrs args ...\n";
    std::cerr << "player_num = 0 or 1 for the computational players\n";
    std::cerr << "player_num = 2 for the server player\n";
    std::cerr << "player_addrs is omitted for player 0\n";
    std::cerr << "player_addrs is p0's hostname for player 1\n";
    std::cerr << "player_addrs is p0's hostname followed by p1's hostname for player 2\n";
    exit(1);
}

static void comp_player_main(boost::asio::io_context &io_context,
    unsigned player, bool preprocessing, const char *p0addr, char **args)
{
    tcp::socket peersock(io_context), serversock(io_context);
    mpcio_setup_computational(player, io_context,
        p0addr, peersock, serversock);
    MPCIO mpcio(player, preprocessing, std::move(peersock),
        std::move(serversock));

    // Queue up the work to be done
    boost::asio::post(io_context, [&]{
        if (preprocessing) {
            preprocessing_comp(mpcio, args);
        }
    });

    // Start another thread; one will perform the work and the other
    // will execute the async_write handlers
    boost::thread t([&]{io_context.run();});
    io_context.run();
    t.join();
}

static void server_player_main(boost::asio::io_context &io_context,
    bool preprocessing, const char *p0addr, const char *p1addr, char **args)
{
    tcp::socket p0sock(io_context), p1sock(io_context);
    mpcio_setup_server(io_context, p0addr, p1addr, p0sock, p1sock);
    MPCServerIO mpcserverio(preprocessing, std::move(p0sock),
        std::move(p1sock));

    // Queue up the work to be done
    boost::asio::post(io_context, [&]{
        if (preprocessing) {
            preprocessing_server(mpcserverio, args);
        }
    });

    // Start another thread; one will perform the work and the other
    // will execute the async_write handlers
    boost::thread t([&]{io_context.run();});
    io_context.run();
    t.join();
}

int main(int argc, char **argv)
{
    char **args = argv+1; // Skip argv[0] (the program name)
    bool preprocessing = false;
    unsigned player = 0;
    const char *p0addr = NULL;
    const char *p1addr = NULL;
    if (argc > 1 && !strcmp("-p", *args)) {
        preprocessing = true;
        ++args;
    }
    if (*args == NULL) {
        // No arguments?
        usage(argv[0]);
    } else {
        player = atoi(*args);
        ++args;
    }
    if (player > 2) {
        usage(argv[0]);
    }
    if (player > 0) {
        if (*args == NULL) {
            usage(argv[0]);
        } else {
            p0addr = *args;
            ++args;
        }
    }
    if (player > 1) {
        if (*args == NULL) {
            usage(argv[0]);
        } else {
            p1addr = *args;
            ++args;
        }
    }

    /*
    std::cout << "Preprocessing = " <<
            (preprocessing ? "true" : "false") << "\n";
    std::cout << "Player = " << player << "\n";
    if (p0addr) {
        std::cout << "Player 0 addr = " << p0addr << "\n";
    }
    if (p1addr) {
        std::cout << "Player 1 addr = " << p1addr << "\n";
    }
    std::cout << "Args =";
    for (char **a = args; *a; ++a) {
        std::cout << " " << *a;
    }
    std::cout << "\n";
    */

    // Make the network connections
    boost::asio::io_context io_context;

    if (player < 2) {
        comp_player_main(io_context, player, preprocessing, p0addr, args);
    } else {
        server_player_main(io_context, preprocessing, p0addr, p1addr, args);
    }

    return 0;
}
