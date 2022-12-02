#include <vector>
#include <bsd/stdlib.h> // arc4random_buf

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
//
// The incoming objects are written into num_threads files in a
// round-robin manner

void preprocessing_comp(MPCIO &mpcio, int num_threads, char **args)
{
    while(1) {
        unsigned char type = 0;
        unsigned int num = 0;
        size_t res = mpcio.serverio.recv(&type, 1);
        if (res < 1 || type == 0) break;
        mpcio.serverio.recv(&num, 4);
        if (type == 0x80) {
            // Multiplication triples
            std::vector<std::ofstream> tripfiles;
            for (int i=0; i<num_threads; ++i) {
                tripfiles.push_back(openfile("triples", mpcio.player, i));
            }
            unsigned thread_num = 0;

            MultTriple T;
            for (unsigned int i=0; i<num; ++i) {
                res = mpcio.serverio.recv(&T, sizeof(T));
                if (res < sizeof(T)) break;
                tripfiles[thread_num].write((const char *)&T, sizeof(T));
                thread_num = (thread_num + 1) % num_threads;
            }
            for (int i=0; i<num_threads; ++i) {
                tripfiles[i].close();
            }
        } else if (type == 0x81) {
            // Multiplication half triples
            std::vector<std::ofstream> halffiles;
            for (int i=0; i<num_threads; ++i) {
                halffiles.push_back(openfile("halves", mpcio.player, i));
            }
            unsigned thread_num = 0;

            HalfTriple H;
            for (unsigned int i=0; i<num; ++i) {
                res = mpcio.serverio.recv(&H, sizeof(H));
                if (res < sizeof(H)) break;
                halffiles[thread_num].write((const char *)&H, sizeof(H));
                thread_num = (thread_num + 1) % num_threads;
            }
            for (int i=0; i<num_threads; ++i) {
                halffiles[i].close();
            }
        }
    }
}

// Create triples (X0,Y0,Z0),(X1,Y1,Z1) such that
// (X0*Y1 + Y0*X1) = (Z0+Z1)
static void create_triples(MPCServerIO &mpcsrvio, unsigned num)
{
    for (unsigned int i=0; i<num; ++i) {
        uint64_t X0, Y0, Z0, X1, Y1, Z1;
        arc4random_buf(&X0, sizeof(X0));
        arc4random_buf(&Y0, sizeof(Y0));
        arc4random_buf(&Z0, sizeof(Z0));
        arc4random_buf(&X1, sizeof(X1));
        arc4random_buf(&Y1, sizeof(Y1));
        Z1 = X0 * Y1 + X1 * Y0 - Z0;
        MultTriple T0, T1;
        T0 = std::make_tuple(X0, Y0, Z0);
        T1 = std::make_tuple(X1, Y1, Z1);
        mpcsrvio.p0io.queue(&T0, sizeof(T0));
        mpcsrvio.p1io.queue(&T1, sizeof(T1));
    }
}

// Create half-triples (X0,Z0),(Y1,Z1) such that
// X0*Y1 = Z0 + Z1
static void create_halftriples(MPCServerIO &mpcsrvio, unsigned num)
{
    for (unsigned int i=0; i<num; ++i) {
        uint64_t X0, Z0, Y1, Z1;
        arc4random_buf(&X0, sizeof(X0));
        arc4random_buf(&Z0, sizeof(Z0));
        arc4random_buf(&Y1, sizeof(Y1));
        Z1 = X0 * Y1 - Z0;
        HalfTriple H0, H1;
        H0 = std::make_tuple(X0, Z0);
        H1 = std::make_tuple(Y1, Z1);
        mpcsrvio.p0io.queue(&H0, sizeof(H0));
        mpcsrvio.p1io.queue(&H1, sizeof(H1));
    }
}

void preprocessing_server(MPCServerIO &mpcsrvio, char **args)
{
    while (*args) {
        char *colon = strchr(*args, ':');
        if (!colon) {
            std::cerr << "Args must be type:num\n";
            ++args;
            continue;
        }
        unsigned num = atoi(colon+1);
        *colon = '\0';
        char *type = *args;
        if (!strcmp(type, "t")) {
            unsigned char typetag = 0x80;
            mpcsrvio.p0io.queue(&typetag, 1);
            mpcsrvio.p0io.queue(&num, 4);
            mpcsrvio.p1io.queue(&typetag, 1);
            mpcsrvio.p1io.queue(&num, 4);

            create_triples(mpcsrvio, num);
        } else if (!strcmp(type, "h")) {
            unsigned char typetag = 0x81;
            mpcsrvio.p0io.queue(&typetag, 1);
            mpcsrvio.p0io.queue(&num, 4);
            mpcsrvio.p1io.queue(&typetag, 1);
            mpcsrvio.p1io.queue(&num, 4);

            create_halftriples(mpcsrvio, num);
        }
        ++args;
    }
    // That's all
    unsigned char typetag = 0x00;
    mpcsrvio.p0io.queue(&typetag, 1);
    mpcsrvio.p1io.queue(&typetag, 1);
    mpcsrvio.p0io.send();
    mpcsrvio.p1io.send();
}
