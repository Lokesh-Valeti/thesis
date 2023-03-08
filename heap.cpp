#include <functional>

#include "types.hpp"

#include "duoram.hpp"

#include "cell.hpp"

#include  "heap.hpp"


void MinHeap::initialize(int num_players, size_t size) {
    this -> MAX_SIZE = size;
    this -> num_items = 0;
    oram = new Duoram < RegAS > (num_players, size);
}

RegAS reconstruct_AS(MPCTIO & tio, yield_t & yield, RegAS AS) {
    RegAS peer_AS;
    RegAS reconstructed_AS = AS;
    if (tio.player() == 1) {
        tio.queue_peer( & AS, sizeof(AS));
    } else {
        RegAS peer_AS;
        tio.recv_peer( & peer_AS, sizeof(peer_AS));
        reconstructed_AS += peer_AS;
    }

    yield();

    if (tio.player() == 0) {
        tio.queue_peer( & AS, sizeof(AS));
    } else {
        RegAS peer_flag;
        tio.recv_peer( & peer_AS, sizeof(peer_AS));
        reconstructed_AS += peer_AS;
    }

    return reconstructed_AS;
}

RegXS reconstruct_XS(MPCTIO & tio, yield_t & yield, RegXS XS) {
    RegXS peer_XS;
    RegXS reconstructed_XS = XS;
    if (tio.player() == 1) {
        tio.queue_peer( & XS, sizeof(XS));
    } else {
        RegXS peer_XS;
        tio.recv_peer( & peer_XS, sizeof(peer_XS));
        reconstructed_XS ^= peer_XS;
    }

    yield();

    if (tio.player() == 0) {
        tio.queue_peer( & XS, sizeof(XS));
    } else {
        RegXS peer_flag;
        tio.recv_peer( & peer_XS, sizeof(peer_XS));
        reconstructed_XS ^= peer_XS;
    }

    return reconstructed_XS;
}
bool reconstruct_flag(MPCTIO & tio, yield_t & yield, RegBS flag) {
    RegBS peer_flag;
    RegBS reconstructed_flag = flag;
    if (tio.player() == 1) {
        tio.queue_peer( & flag, sizeof(flag));
    } else {
        RegBS peer_flag;
        tio.recv_peer( & peer_flag, sizeof(peer_flag));
        reconstructed_flag ^= peer_flag;
    }

    yield();

    if (tio.player() == 0) {
        tio.queue_peer( & flag, sizeof(flag));
    } else {
        RegBS peer_flag;
        tio.recv_peer( & peer_flag, sizeof(peer_flag));
        reconstructed_flag ^= peer_flag;
    }

    return reconstructed_flag.bshare;
}

int MinHeap::insert(MPCTIO tio, yield_t & yield, RegAS val) {
    auto HeapArray = oram -> flat(tio, yield);
    num_items++;
    std::cout << "num_items = " << num_items << std::endl;

    std::cout << "we are adding in: " << std::endl;
    reconstruct_AS(tio, yield, val);
    val.dump();
    yield();
    size_t childindex = num_items;
    size_t parent = childindex / 2;
    std::cout << "childindex = " << childindex << std::endl;
    std::cout << "parent = " << parent << std::endl;
    HeapArray[num_items] = val;
    RegAS tmp = HeapArray[num_items];
    reconstruct_AS(tio, yield, tmp);
    yield();
    while (parent > 0) {
        //std::cout << "while loop\n";
        RegAS sharechild = HeapArray[childindex];
        RegAS shareparent = HeapArray[parent];
        // RegAS sharechildrec = reconstruct_AS(tio, yield, sharechild);
        // yield();
        // RegAS shareparentrec = reconstruct_AS(tio, yield, shareparent);
        // yield();
        // std::cout << "\nchild reconstruct_AS = \n";
        // sharechildrec.dump();
        // std::cout << "\nparent reconstruct_AS = \n";
        // shareparentrec.dump();
        // std::cout << "\n----\n";
        CDPF cdpf = tio.cdpf(yield);
        RegAS diff = sharechild - shareparent;
        // std::cout << "diff = " << std::endl;
        // RegAS diff_rec = reconstruct_AS(tio, yield, diff);
        // diff_rec.dump();
        // std::cout << std::endl << std::endl;
        auto[lt, eq, gt] = cdpf.compare(tio, yield, diff, tio.aes_ops());
        //auto lteq = lt ^ eq;
        // bool lteq_rec = reconstruct_flag(tio, yield, lteq);
        // yield();   
        // bool lt_rec = reconstruct_flag(tio, yield, lt);
        // yield();
        // std::cout <<"lt_rec = " << (int) lt_rec << std::endl;
        // std::cout << std::endl;
        // bool eq_rec = reconstruct_flag(tio, yield, eq);
        // yield();
        // std::cout <<"eq_rec = " << (int) eq_rec << std::endl;
        // std::cout << std::endl;
        // yield();   
        // bool gt_rec = reconstruct_flag(tio, yield, gt);
        // yield();
        // std::cout <<"gt_rec = " << (int) gt_rec << std::endl;
        // std::cout << std::endl;
        // if(lteq_rec) {
        //   if(sharechildrec.ashare > shareparentrec.ashare) {
        //     std::cout << "\nchild reconstruct_AS = \n";
        //     sharechildrec.dump();
        //     std::cout << "\nchild reconstruct_AS = \n";
        //     shareparentrec.dump();
        //     std::cout << "\n----\n";
        //   }
        //   assert(sharechildrec.ashare <= shareparentrec.ashare);
        // }
        // if(gt_rec) {
        //   if(sharechildrec.ashare < shareparentrec.ashare) {
        //     std::cout << "\nchild reconstruct_AS = \n";
        //     sharechildrec.dump();
        //     std::cout << "\nchild reconstruct_AS = \n";
        //     shareparentrec.dump();
        //     std::cout << "\n----\n";
        //   }
        //     assert(sharechildrec.ashare > shareparentrec.ashare);
        //   }

        //   std::cout << "child = " << child << std::endl;
        //   sharechildrec  =  reconstruct_AS(tio, yield, sharechild);
        //   sharechildrec.dump();
        //   yield();
        //   std::cout << "parent = " << parent << std::endl;
        //   shareparentrec =  reconstruct_AS(tio, yield, shareparent);
        //   shareparentrec.dump();
        //   yield();

        // std::cout << "\n^^^ before mpc_oswap\n";
        mpc_oswap(tio, yield, sharechild, shareparent, lt, 64);

        HeapArray[childindex] = sharechild;
        HeapArray[parent] = shareparent;

        // std::cout << "child = " << child << std::endl;
        // sharechildrec  =  reconstruct_AS(tio, yield, sharechild);
        // sharechildrec.dump();
        // yield();
        // std::cout << "parent = " << parent << std::endl;
        // shareparentrec =  reconstruct_AS(tio, yield, shareparent);
        // shareparentrec.dump();
        // yield();
        // std::cout << "\n^^^after mpc_oswap\n";
        // assert(sharechildrec.ashare >= shareparentrec.ashare);

        // std::cout << "we asserted that: \n";
        // sharechildrec.dump();
        // std::cout << std::endl << " < " << std::endl;
        // shareparentrec.dump();
        // std::cout << "\n ----- \n";

        childindex = parent;
        parent = childindex / 2;
    }

    return 1;
}

int MinHeap::verify_heap_property(MPCTIO tio, yield_t & yield) {
    std::cout << std::endl << std::endl << "verify_heap_property is being called " << std::endl;
    auto HeapArray = oram -> flat(tio, yield);

    RegAS heapreconstruction[num_items];
    for (size_t j = 0; j <= num_items; ++j) {
        RegAS tmp = HeapArray[j];
        heapreconstruction[j] = reconstruct_AS(tio, yield, tmp);
        yield();
        //  heapreconstruction[j].dump();
        //  std::cout << std::endl;
    }

    for (size_t j = 1; j < num_items / 2; ++j) {
        if (heapreconstruction[j].ashare > heapreconstruction[2 * j].ashare) {
            std::cout << "heap property failure\n\n";
            std::cout << "j = " << j << std::endl;
            heapreconstruction[j].dump();
            std::cout << std::endl;
            std::cout << "2*j = " << 2 * j << std::endl;
            heapreconstruction[2 * j].dump();
            std::cout << std::endl;

        }
        assert(heapreconstruction[j].ashare <= heapreconstruction[2 * j].ashare);
        assert(heapreconstruction[j].ashare <= heapreconstruction[2 * j + 1].ashare);

    }

    return 1;
}

void verify_parent_children_heaps(MPCTIO tio, yield_t & yield, RegAS parent, RegAS leftchild, RegAS rightchild) {
    RegAS parent_reconstruction = reconstruct_AS(tio, yield, parent);
    yield();
    RegAS leftchild_reconstruction = reconstruct_AS(tio, yield, leftchild);
    yield();
    RegAS rightchild_reconstruction = reconstruct_AS(tio, yield, rightchild);
    yield();

    assert(parent_reconstruction.ashare <= leftchild_reconstruction.ashare);
    assert(parent_reconstruction.ashare <= rightchild_reconstruction.ashare);
}

//  Let "x" be the root, and let "y" and "z" be the left and right children
//  For an array, we have A[i] = x, A[2i] = y, A[2i + 1] = z.
//  We want x \le y, and x \le z.
//  The steps are as follows:
//  Step 1: compare(y,z);  (1st call to to MPC Compare)
//  Step 2: smaller = min(y,z); This is done with an mpcselect (1st call to mpcselect)
//  Step 3: if(smaller == y) then smallerindex = 2i else smalleindex = 2i + 1;
//  Step 4: compare(x,smaller); (2nd call to to MPC Compare)
//  Step 5: smallest = min(x, smaller);  (2nd call to mpcselect)
//  Step 6: otherchild = max(x, smaller)   
//  Step 7: A[i] \gets smallest   (1st Duoam Write)
//  Step 8: A[smallerindex] \gets otherchild (2nd Duoam Write)
//  Overall restore_heap_property takes 2 MPC Comparisons, 2 MPC Selects, and 2 Duoram Writes
RegXS MinHeap::restore_heap_property(MPCTIO tio, yield_t & yield, RegXS index) {
    RegAS smallest;
    auto HeapArray = oram -> flat(tio, yield);
    RegAS parent = HeapArray[index];
    RegXS leftchildindex = index;
    leftchildindex = index << 1;

    RegXS rightchildindex;
    rightchildindex.xshare = leftchildindex.xshare ^ (tio.player());

    RegAS leftchild = HeapArray[leftchildindex];
    RegAS rightchild = HeapArray[rightchildindex];
    RegAS sum = parent + leftchild + rightchild;
    CDPF cdpf = tio.cdpf(yield);
    auto[lt, eq, gt] = cdpf.compare(tio, yield, leftchild - rightchild, tio.aes_ops());

    RegXS smallerindex;
    //mpc_select(tio, yield, smallerindex, lt, rightchildindex, leftchildindex, 64); 

    smallerindex = leftchildindex ^ lt;
    //smallerindex stores either the index of the left or child (whichever has the smaller value) 

    RegAS smallerchild;
    mpc_select(tio, yield, smallerchild, lt, rightchild, leftchild, 64);
    // the value smallerchild holds smaller of left and right child 
    RegAS largerchild = sum - parent - smallerchild;
    CDPF cdpf0 = tio.cdpf(yield);
    auto[lt0, eq0, gt0] = cdpf0.compare(tio, yield, smallerchild - parent, tio.aes_ops());
    //comparison between the smallerchild and the parent

    mpc_select(tio, yield, smallest, lt0, parent, smallerchild, 64);
    // smallest holds smaller of left/right child and parent

    RegAS otherchild;
    //mpc_select(tio, yield, otherchild, gt0, parent, smallerchild, 64);
    otherchild = sum - smallest - largerchild;
    // otherchild holds max(min(leftchild, rightchild), parent)

    HeapArray[index] = smallest;
    HeapArray[smallerindex] = otherchild;

    //verify_parent_children_heaps(tio, yield, HeapArray[index], HeapArray[leftchildindex] , HeapArray[rightchildindex]);

    return smallerindex;
}

/*

*/
RegXS MinHeap::restore_heap_property_at_root(MPCTIO tio, yield_t & yield) {
    size_t index = 1;
    auto HeapArray = oram -> flat(tio, yield);
    RegAS parent = HeapArray[index];
    RegAS leftchild = HeapArray[2 * index];
    RegAS rightchild = HeapArray[2 * index + 1];
    CDPF cdpf = tio.cdpf(yield);
    auto[lt, eq, gt] = cdpf.compare(tio, yield, leftchild - rightchild, tio.aes_ops()); // c_1 in the paper
    RegAS smallerchild;
    mpc_select(tio, yield, smallerchild, lt, rightchild, leftchild, 64); // smallerchild holds smaller of left and right child
    CDPF cdpf0 = tio.cdpf(yield);
    auto[lt0, eq0, gt0] = cdpf0.compare(tio, yield, smallerchild - parent, tio.aes_ops()); //c_0 in the paper
    RegAS smallest;
    mpc_select(tio, yield, smallest, lt0, parent, smallerchild, 64); // smallest holds smaller of left/right child and parent
    RegAS larger_p;
    mpc_select(tio, yield, larger_p, gt0, parent, smallerchild, 64); // smallest holds smaller of left/right child and parent
    parent = smallest;
    leftchild = larger_p;
    HeapArray[index] = smallest;
    //verify_parent_children_heaps(tio, yield, parent, leftchild, rightchild);
    //RegAS smallerindex;
    RegXS smallerindex(lt);
    uint64_t leftchildindex = (2 * index);
    uint64_t rightchildindex = (2 * index) + 1;
    // smallerindex2 &= (leftchildindex ^ rightchildindex);
    // smallerindex2.xshare ^= leftchildindex;
    smallerindex = (RegXS(lt) & leftchildindex) ^ (RegXS(gt) & rightchildindex);
    // RegXS smallerindex_reconstruction2 = reconstruct_XS(tio, yield, smallerindex2);
    // yield();
    HeapArray[smallerindex] = larger_p;
    // std::cout << "smallerindex XOR (root) == \n";
    // smallerindex_reconstruction2.dump();
    // std::cout << "\n\n ---? \n";
    return smallerindex;
}

RegAS MinHeap::extract_min(MPCTIO tio, yield_t & yield) {

    RegAS minval;
    auto HeapArray = oram -> flat(tio, yield);
    minval = HeapArray[1];
    HeapArray[1] = RegAS(HeapArray[num_items]);

    num_items--;
    RegXS smaller = restore_heap_property_at_root(tio, yield);
    std::cout << "num_items = " << num_items << std::endl;
    size_t height = std::log2(num_items);
    std::cout << "height = " << height << std::endl;
    for (size_t i = 0; i < height; ++i) {
        smaller = restore_heap_property(tio, yield, smaller);
    }
    return minval;
}

void Heap(MPCIO & mpcio,
    const PRACOptions & opts, char ** args) {
    nbits_t depth = atoi(args[0]);
    std::cout << "print arguements " << std::endl;
    std::cout << args[0] << std::endl;

    if ( * args) {
        depth = atoi( * args);
        ++args;
    }

    size_t items = (size_t(1) << depth) - 1;

    if ( * args) {
        items = atoi( * args);
        ++args;
    }

    std::cout << "items = " << items << std::endl;

    MPCTIO tio(mpcio, 0, opts.num_threads);

    run_coroutines(tio, [ & tio, depth, items](yield_t & yield) {
        size_t size = size_t(1) << depth;
        std::cout << "size = " << size << std::endl;

        MinHeap tree(tio.player(), size);

        for (size_t j = 0; j < 60; ++j) {
            RegAS inserted_val;
            inserted_val.randomize(62);
            inserted_val.ashare = inserted_val.ashare;
            tree.insert(tio, yield, inserted_val);
        }

        std::cout << std::endl << "=============[Insert Done]================" << std::endl << std::endl;
        tree.verify_heap_property(tio, yield);

        for (size_t j = 0; j < 10; ++j) {
            RegAS minval = tree.extract_min(tio, yield);
            tree.verify_heap_property(tio, yield);
            RegAS minval_reconstruction = reconstruct_AS(tio, yield, minval);
            yield();
            std::cout << "minval_reconstruction = ";
            minval_reconstruction.dump();
            std::cout << std::endl;
        }

        std::cout << std::endl << "=============[Extract Min Done]================" << std::endl << std::endl;
        tree.verify_heap_property(tio, yield);
    });
}