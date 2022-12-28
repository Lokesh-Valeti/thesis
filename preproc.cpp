#include <vector>

#include "types.hpp"
#include "coroutine.hpp"
#include "preproc.hpp"
#include "rdpf.hpp"
#include "cdpf.hpp"

// Keep track of open files that coroutines might be writing into
class Openfiles {
    std::vector<std::ofstream> files;

public:
    class Handle {
        Openfiles &parent;
        size_t idx;
    public:
        Handle(Openfiles &parent, size_t idx) :
            parent(parent), idx(idx) {}

        // Retrieve the ofstream from this Handle
        std::ofstream &os() const { return parent.files[idx]; }
    };

    Handle open(const char *prefix, unsigned player,
        unsigned thread_num, nbits_t depth = 0);

    void closeall();
};

// Open a file for writing with name the given prefix, and ".pX.tY"
// suffix, where X is the (one-digit) player number and Y is the thread
// number.  If depth D is given, use "D.pX.tY" as the suffix.
Openfiles::Handle Openfiles::open(const char *prefix, unsigned player,
    unsigned thread_num, nbits_t depth)
{
    std::string filename(prefix);
    char suffix[20];
    if (depth > 0) {
        sprintf(suffix, "%02d.p%d.t%u", depth, player%10, thread_num);
    } else {
        sprintf(suffix, ".p%d.t%u", player%10, thread_num);
    }
    filename.append(suffix);
    std::ofstream &f = files.emplace_back(filename);
    if (f.fail()) {
        std::cerr << "Failed to open " << filename << "\n";
        exit(1);
    }
    return Handle(*this, files.size()-1);
}

// Close all the open files
void Openfiles::closeall()
{
    for (auto& f: files) {
        f.close();
    }
    files.clear();
}

// The server-to-computational-peer protocol for sending precomputed
// data is:
//
// One byte: type
//   0x80: Multiplication triple
//   0x81: Multiplication half-triple
//   0x01 to 0x30: RAM DPF of that depth
//   0x40: Comparison DPF
//   0x82: Counter (for testing)
//   0x00: End of preprocessing
//
// Four bytes: number of objects of that type (not sent for type == 0x00)
//
// Then that number of objects
//
// Repeat the whole thing until type == 0x00 is received

void preprocessing_comp(MPCIO &mpcio, const PRACOptions &opts, char **args)
{
    int num_threads = opts.num_threads;
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcio, &opts, thread_num] {
            MPCTIO tio(mpcio, thread_num);
            Openfiles ofiles;
            std::vector<coro_t> coroutines;
            while(1) {
                unsigned char type = 0;
                unsigned int num = 0;
                size_t res = tio.recv_server(&type, 1);
                if (res < 1 || type == 0) break;
                tio.recv_server(&num, 4);
                if (type == 0x80) {
                    // Multiplication triples
                    auto tripfile = ofiles.open("triples",
                        mpcio.player, thread_num);

                    for (unsigned int i=0; i<num; ++i) {
                        coroutines.emplace_back(
                            [&, tripfile](yield_t &yield) {
                                yield();
                                MultTriple T = tio.triple(yield);
                                tripfile.os() << T;
                            });
                    }
                } else if (type == 0x81) {
                    // Multiplication half triples
                    auto halffile = ofiles.open("halves",
                        mpcio.player, thread_num);

                    for (unsigned int i=0; i<num; ++i) {
                        coroutines.emplace_back(
                            [&, halffile](yield_t &yield) {
                                yield();
                                HalfTriple H = tio.halftriple(yield);
                                halffile.os() << H;
                            });
                    }
                } else if (type >= 0x01 && type <= 0x30) {
                    // RAM DPFs
                    auto tripfile = ofiles.open("rdpf",
                        mpcio.player, thread_num, type);
                    for (unsigned int i=0; i<num; ++i) {
                        coroutines.emplace_back(
                            [&, tripfile, type](yield_t &yield) {
                                yield();
                                RDPFTriple rdpftrip =
                                    tio.rdpftriple(yield, type, opts.expand_rdpfs);
                                printf("dep  = %d\n", type);
                                printf("usi0 = %016lx\n", rdpftrip.dpf[0].unit_sum_inverse);
                                printf("sxr0 = %016lx\n", rdpftrip.dpf[0].scaled_xor.xshare);
                                printf("usi1 = %016lx\n", rdpftrip.dpf[1].unit_sum_inverse);
                                printf("sxr1 = %016lx\n", rdpftrip.dpf[1].scaled_xor.xshare);
                                printf("usi2 = %016lx\n", rdpftrip.dpf[2].unit_sum_inverse);
                                printf("sxr2 = %016lx\n", rdpftrip.dpf[2].scaled_xor.xshare);
                                tripfile.os() << rdpftrip;
                            });
                    }
                } else if (type == 0x40) {
                    // Comparison DPFs
                    auto cdpffile = ofiles.open("cdpf",
                        mpcio.player, thread_num);

                    for (unsigned int i=0; i<num; ++i) {
                        coroutines.emplace_back(
                            [&, cdpffile](yield_t &yield) {
                                yield();
                                CDPF C = tio.cdpf(yield);
                                cdpffile.os() << C;
                            });
                    }
                } else if (type == 0x82) {
                    coroutines.emplace_back(
                        [&, num](yield_t &yield) {
                            yield();
                            unsigned int istart = 0x31415080;
                            for (unsigned int i=istart; i<istart+num; ++i) {
                                tio.queue_peer(&i, sizeof(i));
                                tio.queue_server(&i, sizeof(i));
                                yield();
                                unsigned int peeri, srvi;
                                tio.recv_peer(&peeri, sizeof(peeri));
                                tio.recv_server(&srvi, sizeof(srvi));
                                if (peeri != i || srvi != i) {
                                    printf("Incorrect counter received: "
                                        "peer=%08x srv=%08x\n", peeri,
                                        srvi);
                                }
                            }
                        });
                }
            }
            run_coroutines(tio, coroutines);
            ofiles.closeall();
        });
    }
    pool.join();
}

void preprocessing_server(MPCServerIO &mpcsrvio, const PRACOptions &opts, char **args)
{
    int num_threads = opts.num_threads;
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcsrvio, &opts, thread_num, args] {
            char **threadargs = args;
            MPCTIO stio(mpcsrvio, thread_num);
            Openfiles ofiles;
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
                        coroutines.emplace_back(
                            [&](yield_t &yield) {
                                yield();
                                stio.triple(yield);
                            });
                    }
                } else if (!strcmp(type, "h")) {
                    unsigned char typetag = 0x81;
                    stio.queue_p0(&typetag, 1);
                    stio.queue_p0(&num, 4);
                    stio.queue_p1(&typetag, 1);
                    stio.queue_p1(&num, 4);

                    for (unsigned int i=0; i<num; ++i) {
                        coroutines.emplace_back(
                            [&](yield_t &yield) {
                                yield();
                                stio.halftriple(yield);
                            });
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

                        auto pairfile = ofiles.open("rdpf",
                            mpcsrvio.player, thread_num, depth);
                        for (unsigned int i=0; i<num; ++i) {
                            coroutines.emplace_back(
                                [&, pairfile, depth](yield_t &yield) {
                                    yield();
                                    RDPFPair rdpfpair = stio.rdpfpair(yield, depth);
                                printf("usi0 = %016lx\n", rdpfpair.dpf[0].unit_sum_inverse);
                                printf("sxr0 = %016lx\n", rdpfpair.dpf[0].scaled_xor.xshare);
                                printf("dep0 = %d\n", rdpfpair.dpf[0].depth());
                                printf("usi1 = %016lx\n", rdpfpair.dpf[1].unit_sum_inverse);
                                printf("sxr1 = %016lx\n", rdpfpair.dpf[1].scaled_xor.xshare);
                                printf("dep1 = %d\n", rdpfpair.dpf[1].depth());
                                    if (opts.expand_rdpfs) {
                                        rdpfpair.dpf[0].expand(stio.aes_ops());
                                        rdpfpair.dpf[1].expand(stio.aes_ops());
                                    }
                                    pairfile.os() << rdpfpair;
                                });
                        }
                    }
                } else if (!strcmp(type, "c")) {
                    unsigned char typetag = 0x40;
                    stio.queue_p0(&typetag, 1);
                    stio.queue_p0(&num, 4);
                    stio.queue_p1(&typetag, 1);
                    stio.queue_p1(&num, 4);

                    for (unsigned int i=0; i<num; ++i) {
                        coroutines.emplace_back(
                            [&](yield_t &yield) {
                                yield();
                                stio.cdpf(yield);
                            });
                    }
                } else if (!strcmp(type, "i")) {
                    unsigned char typetag = 0x82;
                    stio.queue_p0(&typetag, 1);
                    stio.queue_p0(&num, 4);
                    stio.queue_p1(&typetag, 1);
                    stio.queue_p1(&num, 4);

                    coroutines.emplace_back(
                        [&, num] (yield_t &yield) {
                            unsigned int istart = 0x31415080;
                            yield();
                            for (unsigned int i=istart; i<istart+num; ++i) {
                                stio.queue_p0(&i, sizeof(i));
                                stio.queue_p1(&i, sizeof(i));
                                yield();
                                unsigned int p0i, p1i;
                                stio.recv_p0(&p0i, sizeof(p0i));
                                stio.recv_p1(&p1i, sizeof(p1i));
                                if (p0i != i || p1i != i) {
                                    printf("Incorrect counter received: "
                                        "p0=%08x p1=%08x\n", p0i,
                                        p1i);
                                }
                            }
                        });

		}
                free(arg);
                ++threadargs;
            }
            // That's all
            unsigned char typetag = 0x00;
            stio.queue_p0(&typetag, 1);
            stio.queue_p1(&typetag, 1);
            run_coroutines(stio, coroutines);
            ofiles.closeall();
        });
    }
    pool.join();
}
