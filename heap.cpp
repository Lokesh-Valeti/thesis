#include <functional>
#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "heap.hpp"   
/*  
 _Protocol 4_ from PRAC: Round-Efficient 3-Party MPC for Dynamic Data Structures
  Consider the following insertion path with:  (x0 < x1 < x2 < x3 < x4)

        x0                      x0                                 x0           
         \                        \                                 \ 
          x1                       x1                                x1
           \                        \                                 \
            x2                       x2                                x2
             \                        \                                 \         
              x3                                                       NewElement
               \                        \                                 \
                x4                       x3                                x3
                 \                        \                                 \
            NewElement                    x4                                 x4
            
      (Path with new element)       (binary search to determine             (After insertion)
                                     the point where New Element 
                                     should be and shift the elements 
                                     from that point down the path 
                                     from the point)     

 The protocol begins by adding an empty node at the end of the heap array.
 The key observation is that after \heapinsert is complete, the only entries
 that might change are the ones on the path from the root to this new node. 
 Further, since the number of entries in the heap is public, \emph{which} entries in $\database$ form this path is also public. 
 The Observation is that this path (from root to leaf) starts off sorted, and will end up with the new element NewElement inserted into the correct position so as to keep the path sorted. 
 The path is of length $\lg \heapsize$, so we use our binary search to find the appropriate insertion position with a single IDPF of height $\lg \lg \heapsize$.
 The advice bits of that IDPF will then be bit shares of a vector $\flag = [0,0,1,0,\dots,0]$ with the $1$ indicating the position at which the new value must be inserted.
 The shares of $\flag$ are (locally) converted to shares of $\u = [0,0,1,1,\dots,1]$ by taking running XORs.
 The bits of $\flag$ and $\u$ are used in $2 \lg \heapsize$ parallel Flag-Word multiplications to shift the elements greater than $\insertval$ down one position, and to write $\insertval$ into the resulting hole, with a single message of communication. 
*/
int MinHeap::insert_optimized(MPCTIO tio, yield_t & yield, RegAS val) {
    auto HeapArray = oram.flat(tio, yield);
    num_items++;
    typename Duoram<RegAS>::Path old_P(HeapArray, tio, yield, num_items);
    const RegXS foundidx = old_P.binary_search(val);
    size_t childindex = num_items;
    uint64_t height = std::ceil(std::log2(num_items  + 1)) + 1;
    RegAS zero;
    zero.ashare = 0;
    HeapArray[childindex] = zero;
    typename Duoram<RegAS>::Path P(HeapArray, tio, yield, num_items);
    
    #ifdef VERBOSE
    uint64_t val_reconstruction = mpc_reconstruct(tio, yield, val, 64);
    std::cout << "val_reconstruction = " << val_reconstruction << std::endl;
    #endif

    uint64_t  logheight = std::ceil(double(std::log2(height))) + 1;
    std::vector<RegBS> flag(height+1);
    std::vector<RegBS> u(height+1);
    typename Duoram<RegAS>::template OblivIndex<RegXS,1> oidx(tio, yield, foundidx, logheight);
    u = oidx.unit_vector(tio, yield, 1 << logheight, foundidx);
    
    #ifdef VERBOSE
    uint64_t foundidx_reconstruction = mpc_reconstruct(tio, yield, foundidx);
    std::cout << "foundidx_reconstruction = " << foundidx_reconstruction << std::endl;
    std::cout << std::endl << " =============== " << std::endl;
    for (size_t j = 0; j < height; ++j) {
        uint64_t reconstruction = mpc_reconstruct(tio, yield, u[j]);
        std::cout << " --->> u[" << j << "] = " << reconstruction  <<  std::endl;
    }
    #endif

    for (size_t j = 0; j < height; ++j) {
        if(tio.player() !=2) {
            flag[j] = u[j];
            if(j > 0) u[j] = u[j] ^ u[j-1];
        }
    }
    
    #ifdef VERBOSE
    for (size_t j = 0; j < height; ++j) {
        uint64_t reconstruction = mpc_reconstruct(tio, yield, u[j]);
        std::cout << " --->> [0000111111]][" << j << "] = " << reconstruction << std::endl;
    }
    #endif

    RegAS * path = new RegAS[height+1];
    RegAS * w   = new RegAS[height+1];
    RegAS * v = new RegAS[height+1];
    for (size_t j = 0; j < height+1; ++j) path[j] = P[j];

    std::vector<coro_t> coroutines;
    for (size_t j = 1; j < height+1; ++j) {
        coroutines.emplace_back( 
                [&tio, w, u, path, j](yield_t &yield) {
            mpc_flagmult(tio, yield, w[j], u[j-1], (path[j-1]-path[j]));
        }
        );
        coroutines.emplace_back( 
                [&tio, v, flag, val, path, j](yield_t &yield) {
            mpc_flagmult(tio, yield, v[j-1], flag[j-1], (val - path[j-1]));
        }
        );
    }
    run_coroutines(tio, coroutines);

    for (size_t j = 0; j < height; ++j) P[j] += (w[j] + v[j]);

    #ifdef VERBOSE
    std::cout << "\n\n=================Before===========\n\n";
    for (size_t j = 0; j < height-1; ++j) {
        auto path_rec = mpc_reconstruct(tio, yield, P[j]);
        std::cout << j << " --->: " << path_rec << std::endl;
    }
    std::cout << "\n\n============================\n\n";
    
    std::cout << "\n\n=================Aftter===========\n\n";
    for (size_t j = 0; j < height-1; ++j) {
        auto path_rec = mpc_reconstruct(tio, yield, P[j]);
        std::cout << j << " --->: " << path_rec << std::endl;
    }
    std::cout << "\n\n============================\n\n";
    #endif

    delete[] path;
    delete[] w;
    delete[] v;

    return 1;
}

// The insert protocol works as follows:
// It adds a new element in the last entry of the array
// From the leaf (the element added), compare with its parent (1 oblivious compare)
// This Protocol 3 from PRAC: Round-Efficient 3-Party MPC for Dynamic Data Structures
int MinHeap::insert(MPCTIO tio, yield_t & yield, RegAS val) {
    
    auto HeapArray = oram.flat(tio, yield);
    num_items++;
    
    size_t childindex = num_items;
    size_t parentindex = childindex / 2;
    
    #ifdef VERBOSE
    std::cout << "childindex = " << childindex << std::endl;
    std::cout << "parentindex = " << parentindex << std::endl;
    #endif
    
    HeapArray[num_items] = val;
    typename Duoram<RegAS>::Path P(HeapArray, tio, yield, childindex);
    
    while (parentindex > 0) {
        RegAS sharechild = HeapArray[childindex];
        RegAS shareparent = HeapArray[parentindex];
        CDPF cdpf = tio.cdpf(yield);
        RegAS diff = sharechild - shareparent;
        auto[lt, eq, gt] = cdpf.compare(tio, yield, diff, tio.aes_ops());
        auto lteq = lt ^ eq;
        mpc_oswap(tio, yield, sharechild, shareparent, lteq, 64);
        HeapArray[childindex]  = sharechild;
        HeapArray[parentindex] = shareparent;
        childindex = parentindex;
        parentindex = parentindex / 2;
    }

    return 1;
}

// This is only for debugging purposes
// This function just verifies that the heap property is satisfied
int MinHeap::verify_heap_property(MPCTIO tio, yield_t & yield) {

    #ifdef VERBOSE
    std::cout << std::endl << std::endl << "verify_heap_property is being called " << std::endl;
    #endif
    
    auto HeapArray = oram.flat(tio, yield);
    uint64_t * heapreconstruction = new uint64_t[num_items + 1];
    for (size_t j = 0; j < num_items; ++j) {
        heapreconstruction[j] = mpc_reconstruct(tio, yield,  HeapArray[j]);
        #ifdef VERBOSE
        if(tio.player() < 2) std::cout << j << " -----> heapreconstruction[" << j << "] = " << heapreconstruction[j] << std::endl;
        #endif
    }
    
    for (size_t j = 1; j < num_items / 2; ++j) {
        if (heapreconstruction[j] > heapreconstruction[2 * j]) {
            std::cout << "heap property failure\n\n";
            std::cout << "j = " << j << std::endl;
            std::cout << heapreconstruction[j] << std::endl;
            std::cout << "2*j = " << 2 * j << std::endl;
            std::cout << heapreconstruction[2 * j] << std::endl;
        }
        if (heapreconstruction[j] > heapreconstruction[2 * j + 1]) {
            std::cout << "heap property failure\n\n";
            std::cout << "j = " << j << std::endl;
            std::cout << heapreconstruction[j] << std::endl;
            std::cout << "2*j + 1 = " << 2 * j + 1<< std::endl;
            std::cout << heapreconstruction[2 * j + 1] << std::endl;
        }
        assert(heapreconstruction[j] <= heapreconstruction[2 * j]);
        assert(heapreconstruction[j] <= heapreconstruction[2 * j + 1]);
    }

    delete [] heapreconstruction;

    return 1;
}

// This is only for debugging purposes
// The function asserts the fact that both the left and child's reconstruction is <= parent's reconstruction
void verify_parent_children_heaps(MPCTIO tio, yield_t & yield, RegAS parent, RegAS leftchild, RegAS rightchild) {
    uint64_t parent_reconstruction = mpc_reconstruct(tio, yield, parent);
    uint64_t leftchild_reconstruction = mpc_reconstruct(tio, yield, leftchild);
    uint64_t rightchild_reconstruction = mpc_reconstruct(tio, yield, rightchild);
    #ifdef VERBOSE
    std::cout << "parent_reconstruction = " << parent_reconstruction << std::endl;
    std::cout << "leftchild_reconstruction = " << leftchild_reconstruction << std::endl;
    std::cout << "rightchild_reconstruction = " << rightchild_reconstruction << std::endl << std::endl << std::endl;
    #endif
    assert(parent_reconstruction <= leftchild_reconstruction);
    assert(parent_reconstruction <= rightchild_reconstruction);
}

/*
    Protocol 6 from PRAC: Round-Efficient 3-Party MPC for Dynamic Data Structures
    Basic restore heap property has the following functionality:


    Before restoring heap property:                      z
                                                        /  \ 
                                                       y    x


    After restoring heap property:        if(y < x AND z < y)       if(y < x AND z > y)        if(y > x AND z < x)           if(y > x AND z > x)

                                                    z                         x                        y                              x    
                                                   /  \                      / \                      / \                            / \ 
                                                  y    x                    y   z                    x    z                         y   z

                                               
    The protocol works as follows:
    In the first step we compare left and right children
    Next, we compare the smaller child with the parent
    If the smaller child is smaller than the parent, we swap the smallerchild it with the root
    Requires three DORAM reads (in parallel) to read the parent, leftchild, and rightchild
    Requires two comparisons a) comparing left and right child b) comparing the smaller child and parent
    Next, we compute the offsets by which the parent and children need to be updated
    Offset computation requires 
      - one flag-flag multiplication
      - two flag-word multiplications
    Requires three DORAM update operations (in parallel) to update the parent, leftchild, and rightchild
    In total we need 3 DORAM reads, 2 Comparisons, 1 Flag-Flag Multiplication, 2-Flag-Word-Multiplication, and 3 DORAM Updates
    
    The function returns the XOR-share of the smallerchild's index
*/

RegXS MinHeap::restore_heap_property(MPCIO & mpcio, MPCTIO tio, yield_t & yield, RegXS index) {
    RegAS smallest;
    auto HeapArray = oram.flat(tio, yield);
    mpcio.reset_stats();
    tio.reset_lamport();
    RegXS leftchildindex = index;
    leftchildindex = index << 1;
    RegXS rightchildindex;
    rightchildindex.xshare = leftchildindex.xshare ^ (tio.player());
    
    RegAS parent, leftchild, rightchild;
    
    #ifdef VERBOSE
    auto index_reconstruction = mpc_reconstruct(tio, yield, index);
    auto leftchildindex_reconstruction = mpc_reconstruct(tio, yield, leftchildindex);
    auto rightchildindex_reconstruction = mpc_reconstruct(tio, yield, rightchildindex);
    std::cout << "index_reconstruction               =  " << index_reconstruction << std::endl;
    std::cout << "leftchildindex_reconstruction      =  " << leftchildindex_reconstruction << std::endl;
    std::cout << "rightchildindex_reconstruction     =  " << rightchildindex_reconstruction << std::endl;
    #endif
    
    std::vector<coro_t> coroutines_read;
    coroutines_read.emplace_back( 
        [&tio, &parent, &HeapArray, index](yield_t &yield) {
        auto Acoro = HeapArray.context(yield);
        parent = Acoro[index];
    }
    );
    coroutines_read.emplace_back( 
        [&tio, &HeapArray, &leftchild, leftchildindex](yield_t &yield) {
        auto Acoro = HeapArray.context(yield);
        leftchild  = Acoro[leftchildindex];
    }
    );
    coroutines_read.emplace_back( 
        [&tio, &rightchild, &HeapArray, rightchildindex](yield_t &yield) {
        auto Acoro = HeapArray.context(yield);
        rightchild = Acoro[rightchildindex];
    }
    );
    
    run_coroutines(tio, coroutines_read);
    CDPF cdpf = tio.cdpf(yield);
    auto[lt_c, eq_c, gt_c] = cdpf.compare(tio, yield, leftchild - rightchild, tio.aes_ops());
    auto lteq = lt_c ^ eq_c;
    RegXS smallerindex;
    RegAS smallerchild;
    
    run_coroutines(tio, [&tio, &smallerindex, lteq, rightchildindex, leftchildindex](yield_t &yield) {
    mpc_select(tio, yield, smallerindex, lteq, rightchildindex, leftchildindex, 64);
    },  [&tio, &smallerchild, lteq, rightchild, leftchild](yield_t &yield) {
        mpc_select(tio, yield, smallerchild, lteq, rightchild, leftchild, 64);
    }
    );
    
    CDPF cdpf0 = tio.cdpf(yield);
    auto[lt_p, eq_p, gt_p] = cdpf0.compare(tio, yield, smallerchild - parent, tio.aes_ops());
    auto lt_p_eq_p = lt_p ^ eq_p;
    
    RegBS ltlt1;
    
    mpc_and(tio, yield, ltlt1, lteq, lt_p_eq_p);
    
    RegAS update_index_by, update_leftindex_by;
    
    run_coroutines(tio, [&tio, &update_leftindex_by, ltlt1, parent, leftchild](yield_t &yield) {
    mpc_flagmult(tio, yield, update_leftindex_by, ltlt1, (parent - leftchild), 64);
    },  [&tio, &update_index_by, lt_p, parent, smallerchild](yield_t &yield) {
        mpc_flagmult(tio, yield, update_index_by, lt_p, smallerchild - parent, 64);
    }
    );
    
    std::vector<coro_t> coroutines;
    coroutines.emplace_back( 
        [&tio, &HeapArray, index, update_index_by](yield_t &yield) {
        auto Acoro = HeapArray.context(yield);
        Acoro[index] += update_index_by;
    }
    );
    coroutines.emplace_back( 
        [&tio, &HeapArray, leftchildindex, update_leftindex_by](yield_t &yield) {
        auto Acoro = HeapArray.context(yield);
        Acoro[leftchildindex] += update_leftindex_by;
    }
    );
    coroutines.emplace_back( 
        [&tio, &HeapArray, rightchildindex, update_index_by, update_leftindex_by](yield_t &yield) {
        auto Acoro = HeapArray.context(yield);
        Acoro[rightchildindex] += -(update_index_by + update_leftindex_by);
    }
    );
    run_coroutines(tio, coroutines);

    #ifdef DEBUG 
            verify_parent_children_heaps(tio, yield, HeapArray[index], HeapArray[leftchildindex] , HeapArray[rightchildindex]);
    #endif
    
    return smallerindex;
}

// This Protocol 7 from PRAC: Round-Efficient 3-Party MPC for Dynamic Data Structures
// This is an optimized version of restoring the heap property
// The onlydifference between the optimized and the basic versions is that, the reads and writes in the optimized are performed using a wide DPF
// Unlike the basic heap property, the function also returns (leftchild > rightchild)
// (leftchild > rightchild) is used in the extract_min to increment the oblivindx by 
auto MinHeap::restore_heap_property_optimized(MPCTIO tio, yield_t & yield, RegXS index, size_t layer, size_t depth, typename Duoram < RegAS > ::template OblivIndex < RegXS, 3 > (oidx)) {
    
    auto HeapArray = oram.flat(tio, yield);

    RegXS leftchildindex = index;
    leftchildindex = index << 1;

    RegXS rightchildindex;
    rightchildindex.xshare = leftchildindex.xshare ^ (tio.player());

    typename Duoram < RegAS > ::Flat P(HeapArray, tio, yield, 1 << layer, 1 << layer);
    typename Duoram < RegAS > ::Flat C(HeapArray, tio, yield, 2 << layer, 2 << layer);
    typename Duoram < RegAS > ::Stride L(C, tio, yield, 0, 2);
    typename Duoram < RegAS > ::Stride R(C, tio, yield, 1, 2);

    RegAS parent_tmp, leftchild_tmp, rightchild_tmp; 

    std::vector<coro_t> coroutines_read;
    coroutines_read.emplace_back( 
    [&tio, &parent_tmp, &P, &oidx](yield_t &yield) { 
            auto Acoro = P.context(yield); 
            parent_tmp = Acoro[oidx];
     });             
   
    coroutines_read.emplace_back( 
    [&tio, &L, &leftchild_tmp, &oidx](yield_t &yield) { 
            auto Acoro = L.context(yield); 
            leftchild_tmp  = Acoro[oidx];  
     }); 

    coroutines_read.emplace_back( 
    [&tio, &R, &rightchild_tmp, &oidx](yield_t &yield) { 
           auto Acoro = R.context(yield); 
           rightchild_tmp = Acoro[oidx];
     }); 

    run_coroutines(tio, coroutines_read);

    CDPF cdpf = tio.cdpf(yield);

    auto[lt, eq, gt] = cdpf.compare(tio, yield, leftchild_tmp - rightchild_tmp, tio.aes_ops());
    auto lteq = lt ^ eq;

    RegXS smallerindex;
    RegAS smallerchild;

    run_coroutines(tio, [&tio, &smallerindex, lteq, rightchildindex, leftchildindex](yield_t &yield)
            { mpc_select(tio, yield, smallerindex, lteq, rightchildindex, leftchildindex, 64);;},
            [&tio, &smallerchild, lt, rightchild_tmp, leftchild_tmp](yield_t &yield)
            { mpc_select(tio, yield, smallerchild, lt, rightchild_tmp, leftchild_tmp, 64);;});

    CDPF cdpf0 = tio.cdpf(yield);
    auto[lt1, eq1, gt1] = cdpf0.compare(tio, yield, smallerchild - parent_tmp, tio.aes_ops());
    
    auto lt1eq1 = lt1 ^ eq1;

    RegBS ltlt1;    

    mpc_and(tio, yield, ltlt1, lteq, lt1eq1);

    RegAS update_index_by, update_leftindex_by;

    run_coroutines(tio, [&tio, &update_leftindex_by, ltlt1, parent_tmp, leftchild_tmp](yield_t &yield)
            { mpc_flagmult(tio, yield, update_leftindex_by, ltlt1, (parent_tmp - leftchild_tmp), 64);},
            [&tio, &update_index_by, lt1eq1, parent_tmp, smallerchild](yield_t &yield)
            {mpc_flagmult(tio, yield, update_index_by, lt1eq1, smallerchild - parent_tmp, 64);} 
            );
    
    std::vector<coro_t> coroutines;

    coroutines.emplace_back( 
    [&tio, &P, &oidx, update_index_by](yield_t &yield) {
            auto Acoro = P.context(yield); 
            Acoro[oidx] += update_index_by;
     });             
   
    coroutines.emplace_back( 
    [&tio, &L,  &oidx, update_leftindex_by](yield_t &yield) {
            auto Acoro = L.context(yield);
            Acoro[oidx] += update_leftindex_by;
     });

    coroutines.emplace_back( 
    [&tio, &R,  &oidx, update_leftindex_by, update_index_by](yield_t &yield) {
            auto Acoro = R.context(yield);
            Acoro[oidx] += -(update_leftindex_by + update_index_by);
     }); 

    run_coroutines(tio, coroutines);

    return std::make_pair(smallerindex, gt);
}


// Intializes the heap array with 0x7fffffffffffff
void MinHeap::initialize(MPCTIO tio, yield_t & yield) {
    auto HeapArray = oram.flat(tio, yield);
    HeapArray.init(0x7fffffffffffff);
}


// This function simply initializes a heap with values 1,2,...,n
// We use this function only to setup our heap 
// to do timing experiments on insert and extractmins
void MinHeap::initialize_heap(MPCTIO tio, yield_t & yield) {
    auto HeapArray = oram.flat(tio, yield);
    std::vector<coro_t> coroutines;
    for (size_t j = 1; j <= num_items; ++j) {
        coroutines.emplace_back(
            [&tio, &HeapArray, j](yield_t &yield) {
            auto Acoro = HeapArray.context(yield);
            RegAS v;
            v.ashare = j * tio.player();
            Acoro[j] = v;
        }
       );
    }
  run_coroutines(tio, coroutines);
}


// Prints the heap
void MinHeap::print_heap(MPCTIO tio, yield_t & yield) {
    auto HeapArray = oram.flat(tio, yield);
    uint64_t * Pjreconstruction = new uint64_t[num_items + 1];
    for (size_t j = 0; j <= num_items; ++j)  Pjreconstruction[j] = mpc_reconstruct(tio, yield, HeapArray[j]);
    for (size_t j = 0; j <= num_items; ++j) {
        if(2 * j < num_items) {
            std::cout << j << "-->> HeapArray[" << j << "] = " << std::dec << Pjreconstruction[j] << ", children are: " << Pjreconstruction[2 * j] << " and " << Pjreconstruction[2 * j + 1] <<  std::endl;
        } else {
            std::cout << j << "-->> HeapArray[" << j << "] = " << std::dec << Pjreconstruction[j] << " is a LEAF " <<  std::endl;
        }
    }

    delete[] Pjreconstruction;
}


/*
    Restore the head property at the root.
       
                root
                /  \ 
       leftchild    rightchild

After restoring heap property:        
if(leftchild < rightchild AND root < leftchild)       if(leftchild < rightchild AND root > leftchild)     if(leftchild > rightchild AND root < rightchild)     if(leftchild > rightchild AND root > rightchild)


               root                                                     rightchild                                              leftchild                                             rightchild    
               /  \                                                        / \                                                     / \                                                     /  \ 
         leftchild  rightchild                                      leftchild  root                                       rightchild  root                                          leftchild   root


    In the first step we compare left and right children
    Next, we compare the smaller child with the root
    If the smaller child is smaller than the root, we swap the smallerchild it with the root
    
    Unlike restore_heap_property, restore_heap_property_at_root begins with three regular (non-DORAM) read operations
    Requires two comparisons a) comparing left and right child b) comparing the smaller child and parent
    Next, we compute the offsets by which the parent and children need to be updated
    Offset computation requires 
      - one flag-flag multiplication
      - two flag-word multiplications
    Requires three DORAM update operations (in parallel) to update the parent, leftchild, and rightchild
    In total we need  2 Comparisons, 1 Flag-Flag Multiplication, 2-Flag-Word-Multiplication, and 3 DORAM Updates
    
    The function returns the XOR-share of the smallerchild's index

*/

auto MinHeap::restore_heap_property_at_root(MPCTIO tio, yield_t & yield, size_t index = 1) {
    auto HeapArray = oram.flat(tio, yield);
    RegAS parent = HeapArray[index];
    RegAS leftchild = HeapArray[2 * index];
    RegAS rightchild = HeapArray[2 * index + 1];
    CDPF cdpf = tio.cdpf(yield);
    auto[lt, eq, gt] = cdpf.compare(tio, yield, leftchild - rightchild, tio.aes_ops());
   
    auto lteq = lt ^ eq;
    RegAS smallerchild;
    mpc_select(tio, yield, smallerchild, lteq, rightchild, leftchild);
    RegXS smallerindex(lt);
   
    uint64_t leftchildindex = (2 * index);
    uint64_t rightchildindex = (2 * index) + 1;
    smallerindex = (RegXS(lteq) & leftchildindex) ^ (RegXS(gt) & rightchildindex);
    CDPF cdpf0 = tio.cdpf(yield);
    auto[lt1, eq1, gt1] = cdpf0.compare(tio, yield, smallerchild - parent, tio.aes_ops());
    auto lt1eq1 = lt1 ^ eq1;
    RegBS ltlt1;
   
    mpc_and(tio, yield, ltlt1, lteq, lt1eq1);
    RegAS update_index_by, update_leftindex_by;
   
    run_coroutines(tio, [&tio, &update_leftindex_by, ltlt1, parent, leftchild](yield_t &yield) {
        mpc_flagmult(tio, yield, update_leftindex_by, ltlt1, (parent - leftchild), 64);
    }, [&tio, &update_index_by, lt1eq1, parent, smallerchild](yield_t &yield) {
        mpc_flagmult(tio, yield, update_index_by, lt1eq1, smallerchild - parent, 64);
    }
    );
   
    std::vector<coro_t> coroutines;
    coroutines.emplace_back( 
        [&tio, &HeapArray, index, update_index_by](yield_t &yield) {
        auto Acoro = HeapArray.context(yield);
        Acoro[index] += update_index_by;
    }
    );
    coroutines.emplace_back( 
        [&tio, &HeapArray, leftchildindex, update_leftindex_by](yield_t &yield) {
        auto Acoro = HeapArray.context(yield);
        Acoro[leftchildindex] += update_leftindex_by;
    }
    );
    coroutines.emplace_back( 
        [&tio, &HeapArray, rightchildindex, update_index_by, update_leftindex_by](yield_t &yield) {
        auto Acoro = HeapArray.context(yield);
        Acoro[rightchildindex] += -(update_index_by + update_leftindex_by);
    }
    );

    run_coroutines(tio, coroutines);
    
    #ifdef VERBOSE
    RegAS new_parent = HeapArray[index];
    RegAS new_left   = HeapArray[leftchildindex];
    RegAS new_right  = HeapArray[rightchildindex];
    uint64_t parent_R  = mpc_reconstruct(tio, yield, new_parent);
    uint64_t left_R    = mpc_reconstruct(tio, yield, new_left);
    uint64_t right_R   = mpc_reconstruct(tio, yield, new_right);
    std::cout << "parent_R = " << parent_R << std::endl;
    std::cout << "left_R = " << left_R << std::endl;
    std::cout << "right_R = " << right_R << std::endl;
    #endif
    
    #ifdef DEBUG
    verify_parent_children_heaps(tio, yield, HeapArray[index], HeapArray[leftchildindex] , HeapArray[rightchildindex]);
    #endif
    
    return std::make_pair(smallerindex, gt);
}


// This Protocol 5 from PRAC: Round-Efficient 3-Party MPC for Dynamic Data Structures
// Like in the paper, there is only version of extract_min
// the optimized version calls the optimized restore_heap_property
RegAS MinHeap::extract_min(MPCIO & mpcio, MPCTIO tio, yield_t & yield, int is_optimized) {

    size_t height = std::log2(num_items);
    RegAS minval;
    auto HeapArray = oram.flat(tio, yield);
    minval = HeapArray[1];
    HeapArray[1] = RegAS(HeapArray[num_items]);
    num_items--;
    auto outroot = restore_heap_property_at_root(tio, yield);
    RegXS smaller = outroot.first;
    
    if(is_optimized > 0) {
     typename Duoram < RegAS > ::template OblivIndex < RegXS, 3 > oidx(tio, yield, height);
     oidx.incr(outroot.second);

     for (size_t i = 0; i < height-1; ++i) {
         auto out = restore_heap_property_optimized(tio, yield, smaller, i + 1, height, typename Duoram < RegAS > ::template OblivIndex < RegXS, 3 > (oidx));
         smaller = out.first;
         oidx.incr(out.second);
     }
    }

    if(is_optimized == 0) {
      for (size_t i = 0; i < height - 1; ++i) {
        smaller = restore_heap_property(mpcio, tio, yield, smaller);
      }
    }
    
    return minval;
}

/*
 This function is *NOT USED* in the evaluation in PRAC: Round-Efficient 3-Party MPC for Dynamic Data Structures
*/  
void MinHeap::heapify_at_level(MPCIO & mpcio, MPCTIO tio, yield_t & yield, size_t index = 1) {
    auto outroot = restore_heap_property_at_root(tio, yield, index);
    RegXS smaller = outroot.first;
    
    #ifdef VERBOSE
     uint64_t smaller_rec = mpc_reconstruct(tio, yield, smaller, 64);
     std::cout << "smaller_rec = " << smaller_rec << std::endl;
     std::cout << "num_items = " << num_items << std::endl;
     std::cout << "index = " << index << std::endl;
    #endif
    
    size_t height =  std::log2(num_items) - std::floor(log2(index)) ;
    
    #ifdef VERBOSE
    std::cout << "height = " << height << std::endl << "===================" << std::endl;
    #endif
    
    for (size_t i = 0; i < height - 1; ++i) {     
    #ifdef VERBOSE
     std::cout << "index = " << index <<  ",  i = " << i << std::endl;
     uint64_t smaller_rec = mpc_reconstruct(tio, yield, smaller, 64);
     std::cout << "[inside loop] smaller_rec = " << smaller_rec << std::endl;
    #endif
    
    smaller = restore_heap_property(mpcio, tio, yield, smaller);
   }
}

// This function is *NOT USED* in the evaluation in PRAC: Round-Efficient 3-Party MPC for Dynamic Data Structures
// This function takes in a random array turns into a heap
void MinHeap::heapify(MPCIO & mpcio, MPCTIO tio, yield_t & yield) {
    size_t startIdx = ((num_items + 1) / 2) - 1;
    for (size_t i = startIdx; i >= 1; i--) {
        heapify_at_level(mpcio, tio, yield, i);
    }
}

void Heap(MPCIO & mpcio,
    
    const PRACOptions & opts, char ** args) {
    int argc = 12;

    int maxdepth = 0;
    int heapdepth = 0;
    size_t n_inserts = 0;
    size_t n_extracts = 0;
    int is_optimized = 0;
    int run_sanity = 0;

    // Process command line arguments
    for (int i = 0; i < argc; i += 2) {
        std::string option = args[i];
        if (option == "-m" && i + 1 < argc) {
            maxdepth = std::atoi(args[i + 1]);
        } else if (option == "-d" && i + 1 < argc) {
            heapdepth = std::atoi(args[i + 1]);
        } else if (option == "-i" && i + 1 < argc) {
            n_inserts = std::atoi(args[i + 1]);
        } else if (option == "-e" && i + 1 < argc) {
            n_extracts = std::atoi(args[i + 1]);
        } else if (option == "-opt" && i + 1 < argc) {
            is_optimized = std::atoi(args[i + 1]);
        } else if (option == "-s" && i + 1 < argc) {
            run_sanity = std::atoi(args[i + 1]);
        }
    }

    std::cout << "maxdepth: " << maxdepth << std::endl;
    std::cout << "heapdepth: " << heapdepth << std::endl;
    std::cout << "n_inserts: " << n_inserts << std::endl;
    std::cout << "n_extracts: " << n_extracts << std::endl;
    std::cout << "is_optimized: " << is_optimized << std::endl;
    std::cout << "run_sanity: " << run_sanity << std::endl;

    
    MPCTIO tio(mpcio, 0, opts.num_threads);
    
    run_coroutines(tio, [ & tio, maxdepth, heapdepth, n_inserts, n_extracts, is_optimized, run_sanity, &mpcio](yield_t & yield) {
        size_t size = size_t(1) << maxdepth;
        MinHeap tree(tio.player(), size);
        tree.initialize(tio, yield);
        tree.num_items = (size_t(1) << heapdepth) - 1;
        tree.initialize_heap(tio, yield);
        std::cout << "\n===== Init Stats =====\n";
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();
        for (size_t j = 0; j < n_inserts; ++j) {
            
            RegAS inserted_val;
            inserted_val.randomize(10);
           
            #ifdef VERBOSE
            inserted_val.ashare = inserted_val.ashare;
            uint64_t inserted_val_rec = mpc_reconstruct(tio, yield, inserted_val, 64);
            std::cout << "inserted_val_rec = " << inserted_val_rec << std::endl << std::endl;
            #endif
            
            if(is_optimized > 0)  tree.insert_optimized(tio, yield, inserted_val);
            if(is_optimized == 0) tree.insert(tio, yield, inserted_val);
        }

        std::cout << "\n===== Insert Stats =====\n";
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        
        if(run_sanity == 1 && n_inserts != 0) tree.verify_heap_property(tio, yield);
         
        mpcio.reset_stats();
        tio.reset_lamport();
        

        
        #ifdef VERBOSE
        tree.print_heap(tio, yield);
        #endif
        
        for (size_t j = 0; j < n_extracts; ++j) {
        tree.extract_min(mpcio, tio, yield, is_optimized);

        #ifdef VERBOSE
        RegAS minval = tree.extract_min(mpcio, tio, yield, is_optimized);
        uint64_t minval_reconstruction = mpc_reconstruct(tio, yield, minval, 64);
        std::cout << "minval_reconstruction = " << minval_reconstruction << std::endl;
        #endif
        
         #ifdef DEBUG
         tree.verify_heap_property(tio, yield);
         #endif 
         
         #ifdef VERBOSE
         tree.print_heap(tio, yield);
         #endif

        }

        std::cout << "\n===== Extract Min Stats =====\n";
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        
        #ifdef VERBOSE
        tree.print_heap(tio, yield);
        #endif

        if(run_sanity == 1 && n_extracts != 0) tree.verify_heap_property(tio, yield);
    }
    );
}