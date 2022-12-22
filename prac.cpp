#include <iostream>
#include <deque>

#include "mpcio.hpp"
#include "preproc.hpp"
#include "online.hpp"
#include "options.hpp"

static void usage(const char *progname)
{
    std::cerr << "Usage: " << progname << " [-p] [-t num] player_num player_addrs args ...\n";
    std::cerr << "-p: preprocessing mode\n";
    std::cerr << "-t num: use num threads\n";
    std::cerr << "-c: store DPFs compressed (default is expanded)\n";
    std::cerr << "-x: use XOR-shared database (default is additive)\n";
    std::cerr << "player_num = 0 or 1 for the computational players\n";
    std::cerr << "player_num = 2 for the server player\n";
    std::cerr << "player_addrs is omitted for player 0\n";
    std::cerr << "player_addrs is p0's hostname for player 1\n";
    std::cerr << "player_addrs is p0's hostname followed by p1's hostname for player 2\n";
    exit(1);
}

static void comp_player_main(boost::asio::io_context &io_context,
    unsigned player, const PRACOptions &opts, const char *p0addr,
    char **args)
{
    std::deque<tcp::socket> peersocks, serversocks;
    mpcio_setup_computational(player, io_context, p0addr,
        opts.num_threads, peersocks, serversocks);
    MPCPeerIO mpcio(player, opts.preprocessing, peersocks, serversocks);

    // Queue up the work to be done
    boost::asio::post(io_context, [&]{
        if (opts.preprocessing) {
            preprocessing_comp(mpcio, opts, args);
        } else {
            online_main(mpcio, opts, args);
        }
    });

    // Start another thread; one will perform the work and the other
    // will execute the async_write handlers
    boost::thread t([&]{io_context.run();});
    io_context.run();
    t.join();
    mpcio.dump_stats(std::cout);
}

static void server_player_main(boost::asio::io_context &io_context,
    const PRACOptions &opts, const char *p0addr,
    const char *p1addr, char **args)
{
    std::deque<tcp::socket> p0socks, p1socks;
    mpcio_setup_server(io_context, p0addr, p1addr,
        opts.num_threads, p0socks, p1socks);
    MPCServerIO mpcserverio(opts.preprocessing, p0socks, p1socks);

    // Queue up the work to be done
    boost::asio::post(io_context, [&]{
        if (opts.preprocessing) {
            preprocessing_server(mpcserverio, opts, args);
        } else {
            online_main(mpcserverio, opts, args);
        }
    });

    // Start another thread; one will perform the work and the other
    // will execute the async_write handlers
    boost::thread t([&]{io_context.run();});
    io_context.run();
    t.join();
    mpcserverio.dump_stats(std::cout);
}

int main(int argc, char **argv)
{
    char **args = argv+1; // Skip argv[0] (the program name)
    PRACOptions opts;
    unsigned player = 0;
    const char *p0addr = NULL;
    const char *p1addr = NULL;
    // Get the options
    while (*args && *args[0] == '-') {
        if (!strcmp("-p", *args)) {
            opts.preprocessing = true;
            ++args;
        } else if (!strcmp("-t", *args)) {
            if (args[1]) {
                opts.num_threads = atoi(args[1]);
                if (opts.num_threads < 1) {
                    usage(argv[0]);
                }
                args += 2;
            } else {
                usage(argv[0]);
            }
        } else if (!strcmp("-c", *args)) {
            opts.expand_rdpfs = false;
            ++args;
        } else if (!strcmp("-x", *args)) {
            opts.use_xor_db = true;
            ++args;
        }
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
    std::cout << "Thread count = " << num_threads << "\n";
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
        comp_player_main(io_context, player, opts, p0addr, args);
    } else {
        server_player_main(io_context, opts, p0addr, p1addr, args);
    }

    return 0;
}
