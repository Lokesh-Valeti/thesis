#include <vector>
#include <bsd/stdlib.h> // arc4random_buf

#include "types.hpp"
#include "preproc.hpp"

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
//   0x01 to 0x40: DPF of that depth
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
                        res = tio.recv_server(&T, sizeof(T));
                        if (res < sizeof(T)) break;
                        tripfile.write((const char *)&T, sizeof(T));
                    }
                    tripfile.close();
                } else if (type == 0x81) {
                    // Multiplication half triples
                    std::ofstream halffile = openfile("halves",
                        mpcio.player, thread_num);

                    HalfTriple H;
                    for (unsigned int i=0; i<num; ++i) {
                        res = tio.recv_server(&H, sizeof(H));
                        if (res < sizeof(H)) break;
                        halffile.write((const char *)&H, sizeof(H));
                    }
                    halffile.close();
                }
            }
        });
    }
    pool.join();
}

// Create triples (X0,Y0,Z0),(X1,Y1,Z1) such that
// (X0*Y1 + Y0*X1) = (Z0+Z1)
static void create_triples(MPCTIO &stio, unsigned num)
{
    for (unsigned int i=0; i<num; ++i) {
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
        stio.queue_p0(&T0, sizeof(T0));
        stio.queue_p1(&T1, sizeof(T1));
    }
}

// Create half-triples (X0,Z0),(Y1,Z1) such that
// X0*Y1 = Z0 + Z1
static void create_halftriples(MPCTIO &stio, unsigned num)
{
    for (unsigned int i=0; i<num; ++i) {
        value_t X0, Z0, Y1, Z1;
        arc4random_buf(&X0, sizeof(X0));
        arc4random_buf(&Z0, sizeof(Z0));
        arc4random_buf(&Y1, sizeof(Y1));
        Z1 = X0 * Y1 - Z0;
        HalfTriple H0, H1;
        H0 = std::make_tuple(X0, Z0);
        H1 = std::make_tuple(Y1, Z1);
        stio.queue_p0(&H0, sizeof(H0));
        stio.queue_p1(&H1, sizeof(H1));
    }
}

void preprocessing_server(MPCServerIO &mpcsrvio, int num_threads, char **args)
{
    boost::asio::thread_pool pool(num_threads);
    for (int thread_num = 0; thread_num < num_threads; ++thread_num) {
        boost::asio::post(pool, [&mpcsrvio, thread_num, args] {
            char **threadargs = args;
            MPCTIO stio(mpcsrvio, thread_num);
            while (*threadargs) {
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

                    create_triples(stio, num);
                } else if (!strcmp(type, "h")) {
                    unsigned char typetag = 0x81;
                    stio.queue_p0(&typetag, 1);
                    stio.queue_p0(&num, 4);
                    stio.queue_p1(&typetag, 1);
                    stio.queue_p1(&num, 4);

                    create_halftriples(stio, num);
                }
                free(arg);
                ++threadargs;
            }
            // That's all
            unsigned char typetag = 0x00;
            stio.queue_p0(&typetag, 1);
            stio.queue_p1(&typetag, 1);
            stio.send();
        });
    }
    pool.join();
}
