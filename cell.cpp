#include <functional>

#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"

struct Cell {
    RegAS key;
    RegXS pointers;
    RegXS value;

    void dump() const {
        printf("[%016lx %016lx %016lx]", key.share(), pointers.share(),
            value.share());
    }

    // You'll need to be able to create a random element, and do the
    // operations +=, +, -=, - (binary and unary).  Note that for
    // XOR-shared fields, + and - are both really XOR.

    inline void randomize() {
        key.randomize();
        pointers.randomize();
        value.randomize();
    }

    inline Cell &operator+=(const Cell &rhs) {
        this->key += rhs.key;
        this->pointers += rhs.pointers;
        this->value += rhs.value;
        return *this;
    }

    inline Cell operator+(const Cell &rhs) const {
        Cell res = *this;
        res += rhs;
        return res;
    }

    inline Cell &operator-=(const Cell &rhs) {
        this->key -= rhs.key;
        this->pointers -= rhs.pointers;
        this->value -= rhs.value;
        return *this;
    }

    inline Cell operator-(const Cell &rhs) const {
        Cell res = *this;
        res -= rhs;
        return res;
    }

    inline Cell operator-() const {
        Cell res;
        res.key = -this->key;
        res.pointers = -this->pointers;
        res.value = -this->value;
        return res;
    }

    // Multiply each field by the local share of the corresponding field
    // in the argument
    inline Cell mulshare(const Cell &rhs) const {
        Cell res = *this;
        res.key.mulshareeq(rhs.key);
        res.pointers.mulshareeq(rhs.pointers);
        res.value.mulshareeq(rhs.value);
        return res;
    }

    // You need a method to turn a leaf node of a DPF into a share of a
    // unit element of your type.  Typically set each RegAS to
    // dpf.unit_as(leaf) and each RegXS or RegBS to dpf.unit_bs(leaf).
    // Note that RegXS will extend a RegBS of 1 to the all-1s word, not
    // the word with value 1.  This is used for ORAM reads, where the
    // same DPF is used for all the fields.
    inline void unit(const RDPF &dpf, DPFnode leaf) {
        key = dpf.unit_as(leaf);
        pointers = dpf.unit_bs(leaf);
        value = dpf.unit_bs(leaf);
    }

    // Perform an update on each of the fields, using field-specific
    // MemRefs constructed from the Shape shape and the index idx
    template <typename Sh, typename U>
    inline static void update(Sh &shape, yield_t &shyield, U idx,
            const Cell &M) {
        run_coroutines(shyield,
            [&shape, &idx, &M] (yield_t &yield) {
                Sh Sh_coro = shape.context(yield);
                Sh_coro[idx].field(&Cell::key) += M.key;
            },
            [&shape, &idx, &M] (yield_t &yield) {
                Sh Sh_coro = shape.context(yield);
                Sh_coro[idx].field(&Cell::pointers) += M.pointers;
            },
            [&shape, &idx, &M] (yield_t &yield) {
                Sh Sh_coro = shape.context(yield);
                Sh_coro[idx].field(&Cell::value) += M.value;
            });
    }
};

template <typename T>
T& operator>>(T& is, Cell &x)
{
    is >> x.key >> x.pointers >> x.value;
    return is;
}

template <typename T>
T& operator<<(T& os, const Cell &x)
{
    os << x.key << x.pointers << x.value;
    return os;
}

DEFAULT_TUPLE_IO(Cell)

void cell(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
{
    nbits_t depth=4;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }

    MPCTIO tio(mpcio, 0, opts.num_threads);
    run_coroutines(tio, [&tio, depth] (yield_t &yield) {
        size_t size = size_t(1)<<depth;
        Duoram<Cell> oram(tio.player(), size);
        auto A = oram.flat(tio, yield);
        Cell c;
        c.key.set(0x0102030405060708);
        c.pointers.set(0x1112131415161718);
        c.value.set(0x2122232425262728);
        A[0] = c;
        RegAS idx;
        Cell expl_read_c = A[0];
        printf("expl_read_c = ");
        expl_read_c.dump();
        printf("\n");
        Cell oram_read_c = A[idx];
        printf("oram_read_c = ");
        oram_read_c.dump();
        printf("\n");

        RegXS valueupdate;
        valueupdate.set(0x4040404040404040 * tio.player());
        RegXS pointersset;
        pointersset.set(0x123456789abcdef0 * tio.player());
        A[1].field(&Cell::value) += valueupdate;
        A[3].field(&Cell::pointers) = pointersset;
        RegXS pointval = A[0].field(&Cell::pointers);
        printf("pointval = ");
        pointval.dump();
        printf("\n");

        idx.set(1 * tio.player());
        RegXS oram_value_read = A[idx].field(&Cell::value);
        printf("oram_value_read = ");
        oram_value_read.dump();
        printf("\n");
        valueupdate.set(0x8080808080808080 * tio.player());
        A[idx].field(&Cell::value) += valueupdate;
        idx.set(2 * tio.player());
        A[idx].field(&Cell::value) = valueupdate;

        c.key.set(0x0102030405060708 * tio.player());
        c.pointers.set(0x1112131415161718 * tio.player());
        c.value.set(0x2122232425262728 * tio.player());
        A[idx] += c;
        idx.set(3 * tio.player());
        A[idx] = c;

        printf("\n");

        if (depth < 10) {
            oram.dump();
            auto R = A.reconstruct();
            if (tio.player() == 0) {
                for(size_t i=0;i<R.size();++i) {
                    printf("\n%04lx ", i);
                    R[i].dump();
                }
                printf("\n");
            }
        }
    });
}
