#include <vector>

#include "types.hpp"
#include "coroutine.hpp"
#include "preproc.hpp"
#include "rdpf.hpp"

// Open a file for writing with name the given prefix, and ".pX.tY"
// suffix, where X is the (one-digit) player number and Y is the thread
// number
static std::ofstream openfile(const char *prefix, unsigned player,
    unsigned thread_num)
{
    std::string filename(prefix);
    char suffix[20];
    sprintf(suffix, ".p%d.t%u", player%10, thread_num);
    filename.append(suffix);
    std::ofstream f;
    f.open(filename);
    if (f.fail()) {
        std::cerr << "Failed to open " << filename << "\n";
        exit(1);
    }
    return f;
}

// The server-to-computational-peer protocol for sending precomputed
// data is:
//
// One byte: type
//   0x80: Multiplication triple
//   0x81: Multiplication half-triple
//   0x01 to 0x30: RAM DPF of that depth
//   0x40: Comparison DPF
//   0x00: End of preprocessing
//
// Four bytes: number of objects of that type (not sent for type == 0x00)
//
// Then that number of objects
//
// Repeat the whole thing until type == 0x00 is received

void preprocessing_comp(MPCIO &mpcio, int num_threads, char **args)
{
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, thread_num] {
            MPCTIO tio(mpcio, thread_num);
            std::vector<coro_t> coroutines;
            while(1) {
                unsigned char type = 0;
                unsigned int num = 0;
                size_t res = tio.recv_server(&type, 1);
                if (res < 1 || type == 0) break;
                tio.recv_server(&num, 4);
                if (type == 0x80) {
                    // Multiplication triples
                    std::ofstream tripfile = openfile("triples",
                        mpcio.player, thread_num);

                    MultTriple T;
                    for (unsigned int i=0; i<num; ++i) {
                        T = tio.triple();
                        tripfile.write((const char *)&T, sizeof(T));
                    }
                    tripfile.close();
                } else if (type == 0x81) {
                    // Multiplication half triples
                    std::ofstream halffile = openfile("halves",
                        mpcio.player, thread_num);

                    HalfTriple H;
                    for (unsigned int i=0; i<num; ++i) {
                        H = tio.halftriple();
                        halffile.write((const char *)&H, sizeof(H));
                    }
                    halffile.close();
                } else if (type >= 0x01 && type <= 0x30) {
                    // RAM DPFs
                    for (unsigned int i=0; i<num; ++i) {
                        coroutines.emplace_back(
                            [&](yield_t &yield) {
                                RegXS ri;
                                ri.randomize(type);
                                RDPF rdpf(tio, yield, ri, type);
                            });
                    }
                }
            }
            run_coroutines(tio, coroutines);
        });
    }
    pool.join();
}

void preprocessing_server(MPCServerIO &mpcsrvio, int num_threads, char **args)
{
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcsrvio, thread_num, args] {
            char **threadargs = args;
            MPCTIO stio(mpcsrvio, thread_num);
            std::vector<coro_t> coroutines;
            if (*threadargs && threadargs[0][0] == 'T') {
                // Per-thread initialization.  The args look like:
                // T0 t:50 h:10 T1 t:20 h:30 T2 h:20

                // Skip to the arg marking our thread
                char us[20];
                sprintf(us, "T%u", thread_num);
                while (*threadargs && strcmp(*threadargs, us)) {
                    ++threadargs;
                }
                // Now skip to the next arg if there is one
                if (*threadargs) {
                    ++threadargs;
                }
            }
            // Stop scanning for args when we get to the end or when we
            // get to another per-thread initialization marker
            while (*threadargs && threadargs[0][0] != 'T') {
                char *arg = strdup(*threadargs);
                char *colon = strchr(arg, ':');
                if (!colon) {
                    std::cerr << "Args must be type:num\n";
                    ++threadargs;
                    free(arg);
                    continue;
                }
                unsigned num = atoi(colon+1);
                *colon = '\0';
                char *type = arg;
                if (!strcmp(type, "t")) {
                    unsigned char typetag = 0x80;
                    stio.queue_p0(&typetag, 1);
                    stio.queue_p0(&num, 4);
                    stio.queue_p1(&typetag, 1);
                    stio.queue_p1(&num, 4);

                    for (unsigned int i=0; i<num; ++i) {
                        stio.triple();
                    }
                } else if (!strcmp(type, "h")) {
                    unsigned char typetag = 0x81;
                    stio.queue_p0(&typetag, 1);
                    stio.queue_p0(&num, 4);
                    stio.queue_p1(&typetag, 1);
                    stio.queue_p1(&num, 4);

                    for (unsigned int i=0; i<num; ++i) {
                        stio.halftriple();
                    }
                } else if (type[0] == 'r') {
                    int depth = atoi(type+1);
                    if (depth < 1 || depth > 48) {
                        std::cerr << "Invalid DPF depth\n";
                    } else {
                        unsigned char typetag = depth;
                        stio.queue_p0(&typetag, 1);
                        stio.queue_p0(&num, 4);
                        stio.queue_p1(&typetag, 1);
                        stio.queue_p1(&num, 4);

                        for (unsigned int i=0; i<num; ++i) {
                            coroutines.emplace_back(
                                [&](yield_t &yield) {
                                    RegXS ri;
                                    RDPF rdpf(stio, yield, ri, depth);
                                });
                        }
                    }
		}
                free(arg);
                ++threadargs;
            }
            // That's all
            unsigned char typetag = 0x00;
            stio.queue_p0(&typetag, 1);
            stio.queue_p1(&typetag, 1);
            run_coroutines(stio, coroutines);
        });
    }
    pool.join();
}
