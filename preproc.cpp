#include <bsd/stdlib.h> // arc4random_buf

#include "preproc.hpp"

// Open a file for writing with name the given prefix, and ".pX" suffix,
// where X is the (one-digit) player number
static std::ofstream openfile(const char *prefix, unsigned player)
{
    std::string filename(prefix);
    char suffix[4];
    sprintf(suffix, ".p%d", player%10);
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

void preprocessing_comp(MPCIO &mpcio, char **args)
{
    while(1) {
        unsigned char type = 0;
        unsigned int num = 0;
        size_t res = mpcio.serverio.recv(&type, 1);
        if (res < 1 || type == 0) break;
        mpcio.serverio.recv(&num, 4);
        if (type == 0x80) {
            // Multiplication triples
            std::ofstream tripfile = openfile("triples", mpcio.player);
            
            MultTriple T;
            for (unsigned int i=0; i<num; ++i) {
                res = mpcio.serverio.recv(&T, sizeof(T));
                if (res < sizeof(T)) break;
                tripfile.write((const char *)&T, sizeof(T));
            }
            tripfile.close();
        }
    }
}

void preprocessing_server(MPCServerIO &mpcsrvio, char **args)
{
    unsigned int numtriples = 100;
    if (*args) {
        numtriples = atoi(*args);
        ++args;
    }
    unsigned char type = 0x80;
    mpcsrvio.p0io.queue(&type, 1);
    mpcsrvio.p0io.queue(&numtriples, 4);
    mpcsrvio.p1io.queue(&type, 1);
    mpcsrvio.p1io.queue(&numtriples, 4);

    for (unsigned int i=0; i<numtriples; ++i) {
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

    // That's all
    type = 0x00;
    mpcsrvio.p0io.queue(&type, 1);
    mpcsrvio.p1io.queue(&type, 1);
    mpcsrvio.p0io.send();
    mpcsrvio.p1io.send();
}
