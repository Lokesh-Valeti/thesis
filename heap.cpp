#include <functional>

#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "heap.hpp"   


// The optimized insertion protocol works as follows
// Do a binary search from the root of the rightmost child
// Binary search returns the index in the path where the new element should go to
// We next do a trickle down. 
// Consdier a path P = [a1, a2, a3, a4, ..., a_n]
// Consider a standard-basis vector [0, ... 0, 1 , 0, ... , 0] (at i)
// we get a new path P = [a1, ... , ai,ai ..., a_n]
// P[i] = newelement   
int MinHeap::insert_optimized(MPCTIO tio, yield_t & yield, RegAS val) {
    
    auto HeapArray = oram.flat(tio, yield);
    num_items++;
    //std::cout << "num_items = " << num_items << std::endl;

    uint64_t height = std::ceil(std::log2(num_items)) + 1;
    
    //std::cout << "height = " << height << std::endl;

    size_t childindex = num_items;
    RegAS zero;    
    zero.ashare = 0;
    HeapArray[childindex] = zero;

    typename Duoram<RegAS>::Path P(HeapArray, tio, yield, childindex);
 
    const RegXS foundidx = P.binary_search(val);
 
    uint64_t  logheight = std::ceil(double(std::log2(height))) + 1; 
  //  std::cout << "logheight = " << logheight << std::endl;
    // RDPF<1> dpf2(tio, yield, foundidx, logheight, false, false);
//    RegBS * flags_array = new RegBS[height];

       
  
    std::vector<RegBS> standard_basis_vector(height+1);
    typename Duoram<RegAS>::template OblivIndex<RegXS,1> oidx(tio, yield, foundidx, logheight);

    auto flags_array = oidx.unit_vector(tio, yield, 1 << logheight, foundidx);
    
    

    //  for(size_t j = 0; j < height; ++j)
    // {
    //      uint64_t reconstruction = mpc_reconstruct(tio, yield, new_stand[j], 64); 
    //     if(reconstruction != 0) std::cout << j << " --->> reconstruction from OblivIndex [new_stand] = " << reconstruction << std::endl;
    // }    
   

 
    for(size_t j = 0; j < height; ++j)
    {
        if(tio.player() !=2) 
        {
         // RDPF<1>::LeafNode leafval = dpf2.leaf(j, tio.aes_ops());
         // flags_array[j] =  dpf2.unit_bs(leafval);
         standard_basis_vector[j] = flags_array[j];        
         if(j > 0) flags_array[j] = flags_array[j] ^ flags_array[j-1];
        }
    }

    // //  #ifdef VERBOSE
    //     for(size_t j = 0; j < height; ++j)
    //     {
    //         uint64_t reconstruction = mpc_reconstruct(tio, yield, standard_basis_vector[j], 64);
    //         if(reconstruction != 0) std::cout << j << " --->> reconstruction Explicitly Calling the DPF [standard_basis_vector] = " << reconstruction << std::endl;
    //     }    
    // //   #endif






    RegAS * z_array2 = new RegAS[height];
    RegAS * z2_tmp   = new RegAS[height]; 
    RegAS * standard_basis_vector_time_value   = new RegAS[height];
    for(size_t j = 0; j < height; ++j) z_array2[j] = P[j];            
    
    //print_heap(tio, yield);

    std::vector<coro_t> coroutines;
    for(size_t j = 1; j < height; ++j)
    {
        coroutines.emplace_back( 
        [&tio, z2_tmp, flags_array, z_array2, j](yield_t &yield) { 
            mpc_flagmult(tio, yield, z2_tmp[j], flags_array[j-1], (z_array2[j-1]-z_array2[j]), 64); 
         });

        coroutines.emplace_back( 
        [&tio, standard_basis_vector_time_value, standard_basis_vector, val, z_array2, j](yield_t &yield) { 
            mpc_flagmult(tio, yield, standard_basis_vector_time_value[j-1], standard_basis_vector[j-1], (val - z_array2[j-1]) , 64); 
         });             
    }    
    
    run_coroutines(tio, coroutines);

    // //#ifdef VERBOSE
    //     for(size_t j = 0; j < height; ++j)
    //     {
    //         int64_t reconstruction = mpc_reconstruct(tio, yield, z2_tmp[j], 64);
    //         std::cout << j << " --->> reconstruction [z2_tmp] = " << reconstruction << std::endl;
    //     }    
    //     std::cout << std::endl << " =============== " << std::endl;
    //     for(size_t j = 0; j < height; ++j)
    //     {
    //         int64_t reconstruction = mpc_reconstruct(tio, yield, flags_array[j], 64);
    //         std::cout << j << " --->> reconstruction [flags_array] = " << reconstruction << std::endl;
    //     }   
    //            std::cout << std::endl << " =============== " << std::endl;
    //     for(size_t j = 0; j < height; ++j)
    //     {
    //         int64_t reconstruction = mpc_reconstruct(tio, yield, standard_basis_vector[j], 64);
    //         std::cout << j << " --->> reconstruction [standard_basis_vector] = " << reconstruction << std::endl;
    //     }   
    //     std::cout << std::endl << " =============== " << std::endl;
    //     for(size_t j = 0; j < height; ++j)
    //     {
    //         int64_t reconstruction = mpc_reconstruct(tio, yield, z_array2[j], 64);
    //         std::cout << j << " --->> reconstruction [z_array2] = " << reconstruction << std::endl;
    //     }   
    // //#endif

    // for(size_t j = 0; j < height; ++j) P[j] += (z2_tmp[j] + standard_basis_vector_time_value[j]);   
    
    // //#ifdef VERBOSE
    //    std::cout << std::endl << " =============== " << std::endl;
    //     for(size_t j = 0; j < height; ++j)
    //     {
    //         int64_t reconstruction = mpc_reconstruct(tio, yield, P[j], 64);
    //         std::cout << j << " --->> reconstruction [P] = " << reconstruction << std::endl;
    //     }  
    //  print_heap(tio, yield);
    // //#endif

    

    // for(size_t j = 1; j < height; ++j) P[j] += z2_tmp[j];
    // typename Duoram<RegAS>::template OblivIndex<RegXS,1> oidx(tio, yield, foundidx, height);
        
     //P[oidx] = val;

    return 1;
}


// The insert protocol works as follows:
// It adds a new element in the last entry of the array
// From the leaf (the element added), compare with its parent (1 oblivious compare)
// If the child is larger, then we do an OSWAP.
int MinHeap::insert(MPCTIO tio, yield_t & yield, RegAS val) {
    auto HeapArray = oram.flat(tio, yield);

    num_items++;
    //std::cout << "num_items = " << num_items << std::endl;

    // uint64_t val_reconstruct = mpc_reconstruct(tio, yield, val);
    // std::cout << "val_reconstruct = " << val_reconstruct << std::endl;

    size_t childindex = num_items;
    size_t parentindex = childindex / 2;
    
    #ifdef VERBOSE
        std::cout << "childindex = " << childindex << std::endl;
        std::cout << "parentindex = " << parentindex << std::endl;
    #endif

    HeapArray[num_items] = val;
    typename Duoram<RegAS>::Path P(HeapArray, tio, yield, childindex);
    //RegXS foundidx = P.binary_search(val);   

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



int MinHeap::verify_heap_property(MPCTIO tio, yield_t & yield) {
    std::cout << std::endl << std::endl << "verify_heap_property is being called " << std::endl;
    auto HeapArray = oram.flat(tio, yield);

    uint64_t heapreconstruction[num_items];
    for (size_t j = 0; j <= num_items; ++j) {
        heapreconstruction[j] = mpc_reconstruct(tio, yield,  HeapArray[j]);
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
        //assert(heapreconstruction[j] <= heapreconstruction[2 * j]);
        //assert(heapreconstruction[j] <= heapreconstruction[2 * j + 1]);

    }

    return 1;
}

void verify_parent_children_heaps(MPCTIO tio, yield_t & yield, RegAS parent, RegAS leftchild, RegAS rightchild) {

    std::cout << "calling this ... \n";
    
    uint64_t parent_reconstruction = mpc_reconstruct(tio, yield, parent);
    
    uint64_t leftchild_reconstruction = mpc_reconstruct(tio, yield, leftchild);
    
    uint64_t rightchild_reconstruction = mpc_reconstruct(tio, yield, rightchild);
    

    assert(parent_reconstruction <= leftchild_reconstruction);
    assert(parent_reconstruction <= rightchild_reconstruction);
}

 

RegXS MinHeap::restore_heap_property(MPCTIO tio, yield_t & yield, RegXS index) {
    RegAS smallest;
    auto HeapArray = oram.flat(tio, yield);


    RegXS leftchildindex = index;
    leftchildindex = index << 1;

    RegXS rightchildindex;
    rightchildindex.xshare = leftchildindex.xshare ^ (tio.player());

    RegAS parent;     //   = HeapArray[index];
    RegAS leftchild;  //  = HeapArray[leftchildindex];
    RegAS rightchild; // = HeapArray[rightchildindex];
 
    std::vector<coro_t> coroutines_read;
    coroutines_read.emplace_back( 
    [&tio, &parent, &HeapArray, index](yield_t &yield) { 
            auto Acoro = HeapArray.context(yield); 
            parent = Acoro[index]; //inserted_val;
     });             
   
    coroutines_read.emplace_back( 
    [&tio, &HeapArray, &leftchild, leftchildindex](yield_t &yield) { 
            auto Acoro = HeapArray.context(yield); 
            leftchild  = Acoro[leftchildindex]; //inserted_val;
     }); 

    coroutines_read.emplace_back( 
    [&tio, &rightchild, &HeapArray, rightchildindex](yield_t &yield) { 
            auto Acoro = HeapArray.context(yield); 
            rightchild = Acoro[rightchildindex];
     }); 

    run_coroutines(tio, coroutines_read);

    //RegAS sum = parent + leftchild + rightchild;

    CDPF cdpf = tio.cdpf(yield);
    auto[lt_c, eq_c, gt_c] = cdpf.compare(tio, yield, leftchild - rightchild, tio.aes_ops());
    auto lteq = lt_c ^ eq_c;    
    RegXS smallerindex;
    RegAS smallerchild;


    #ifdef VERBOSE
        uint64_t LC_rec = mpc_reconstruct(tio, yield, leftchildindex);
        std::cout << "LC_rec = " << LC_rec << std::endl;
    #endif


    // mpc_select(tio, yield, smallerindex, lteq, rightchildindex, leftchildindex, 64);
    // mpc_select(tio, yield, smallerchild, lt_c, rightchild, leftchild, 64);

    run_coroutines(tio, [&tio, &smallerindex, lteq, rightchildindex, leftchildindex](yield_t &yield)
            { mpc_select(tio, yield, smallerindex, lteq, rightchildindex, leftchildindex, 64);},
            [&tio, &smallerchild, lt_c, rightchild, leftchild](yield_t &yield)
            { mpc_select(tio, yield, smallerchild, lt_c, rightchild, leftchild, 64);});


    #ifdef VERBOSE
        uint64_t smallerindex_rec = mpc_reconstruct(tio, yield, smallerindex);
        std::cout << "smallerindex_rec = " << smallerindex_rec << std::endl; 
    #endif

    


    CDPF cdpf0 = tio.cdpf(yield);

    auto[lt_p, eq_p, gt_p] = cdpf0.compare(tio, yield, smallerchild - parent, tio.aes_ops());
    
    auto lt_p_eq_p = lt_p ^ eq_p;

    RegBS ltlt1;
    
    mpc_and(tio, yield, ltlt1, lteq, lt_p_eq_p);

    RegAS update_index_by, update_leftindex_by;

    run_coroutines(tio, [&tio, &update_leftindex_by, ltlt1, parent, leftchild](yield_t &yield)
            { mpc_flagmult(tio, yield, update_leftindex_by, ltlt1, (parent - leftchild), 64);},
            [&tio, &update_index_by, lt_p, parent, smallerchild](yield_t &yield)
            {mpc_flagmult(tio, yield, update_index_by, lt_p, smallerchild - parent, 64);}
            );

    std::vector<coro_t> coroutines;

    coroutines.emplace_back( 
    [&tio, &HeapArray, index, update_index_by](yield_t &yield) { 
            auto Acoro = HeapArray.context(yield); 
            Acoro[index] += update_index_by; //inserted_val;
     });             
   
    coroutines.emplace_back( 
    [&tio, &HeapArray, leftchildindex, update_leftindex_by](yield_t &yield) { 
            auto Acoro = HeapArray.context(yield); 
            Acoro[leftchildindex] += update_leftindex_by; //inserted_val;
     }); 

    coroutines.emplace_back( 
    [&tio, &HeapArray, rightchildindex, update_index_by, update_leftindex_by](yield_t &yield) { 
            auto Acoro = HeapArray.context(yield); 
            Acoro[rightchildindex] += -(update_index_by + update_leftindex_by);
     }); 

    run_coroutines(tio, coroutines);

    
   // verify_parent_children_heaps(tio, yield, HeapArray[index], HeapArray[leftchildindex] , HeapArray[rightchildindex]);

    return smallerindex;
}



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

    // RegAS parent_tmp = P[oidx];
    // RegAS leftchild_tmp = L[oidx];
    // RegAS rightchild_tmp = R[oidx];

    RegAS parent_tmp;     //   = HeapArray[index];
    RegAS leftchild_tmp;  //  = HeapArray[leftchildindex];
    RegAS rightchild_tmp; // = HeapArray[rightchildindex];
 
    std::vector<coro_t> coroutines_read;
    coroutines_read.emplace_back( 
    [&tio, &parent_tmp, &P, &oidx](yield_t &yield) { 
           // auto Acoro = P.context(yield); 
            parent_tmp = P[oidx]; //inserted_val;
     });             
   
    coroutines_read.emplace_back( 
    [&tio, &L, &leftchild_tmp, &oidx](yield_t &yield) { 
            //auto Acoro = L.context(yield); 
            leftchild_tmp  = L[oidx]; //inserted_val;
     }); 

    coroutines_read.emplace_back( 
    [&tio, &R, &rightchild_tmp, &oidx](yield_t &yield) { 
         //  auto Acoro = R.context(yield); 
           rightchild_tmp = R[oidx];
     }); 

    run_coroutines(tio, coroutines_read);




    //RegAS sum = parent_tmp + leftchild_tmp + rightchild_tmp;

    CDPF cdpf = tio.cdpf(yield);

    auto[lt, eq, gt] = cdpf.compare(tio, yield, leftchild_tmp - rightchild_tmp, tio.aes_ops());
    auto lteq = lt ^ eq;

    RegXS smallerindex;
    RegAS smallerchild;

    // mpc_select(tio, yield, smallerindex, lteq, rightchildindex, leftchildindex, 64);
    // mpc_select(tio, yield, smallerchild, lt, rightchild_tmp, leftchild_tmp, 64);

    run_coroutines(tio, [&tio, &smallerindex, lteq, rightchildindex, leftchildindex](yield_t &yield)
            { mpc_select(tio, yield, smallerindex, lteq, rightchildindex, leftchildindex, 64);;},
            [&tio, &smallerchild, lt, rightchild_tmp, leftchild_tmp](yield_t &yield)
            { mpc_select(tio, yield, smallerchild, lt, rightchild_tmp, leftchild_tmp, 64);;});

    CDPF cdpf0 = tio.cdpf(yield);
    auto[lt1, eq1, gt1] = cdpf0.compare(tio, yield, smallerchild - parent_tmp, tio.aes_ops());

    // RegAS z;
    // mpc_flagmult(tio, yield, z, lt1, smallerchild - parent_tmp, 64);

    
    auto lt1eq1 = lt1 ^ eq1;

    RegBS ltlt1;
    // RegAS zz;
    
    mpc_and(tio, yield, ltlt1, lteq, lt1eq1);
    // // mpc_flagmult(tio, yield, zz, ltlt1, (parent_tmp - leftchild_tmp), 64);

    // run_coroutines(tio, [&tio, &ltlt1, lteq, lt1eq1](yield_t &yield)
    //         { mpc_and(tio, yield, ltlt1, lteq, lt1eq1);},
    //         [&tio, &zz, ltlt1, parent_tmp, leftchild_tmp](yield_t &yield)
    //         { mpc_flagmult(tio, yield, zz, ltlt1, (parent_tmp - leftchild_tmp), 64);});

    RegAS update_index_by, update_leftindex_by;

    run_coroutines(tio, [&tio, &update_leftindex_by, ltlt1, parent_tmp, leftchild_tmp](yield_t &yield)
            { mpc_flagmult(tio, yield, update_leftindex_by, ltlt1, (parent_tmp - leftchild_tmp), 64);},
            [&tio, &update_index_by, lt1eq1, parent_tmp, smallerchild](yield_t &yield)
            {mpc_flagmult(tio, yield, update_index_by, lt1eq1, smallerchild - parent_tmp, 64);} 
            );
    

    // RegAS leftchildplusparent = RegAS(HeapArray[index]) + RegAS(HeapArray[leftchildindex]);
    // RegAS tmp = (sum - leftchildplusparent);

    std::vector<coro_t> coroutines;

    coroutines.emplace_back( 
    [&tio, &P, &oidx, update_index_by](yield_t &yield) { 
            auto Acoro = P.context(yield); 
            Acoro[oidx] += update_index_by; //inserted_val;
     });             
   
    coroutines.emplace_back( 
    [&tio, &L,  &oidx, update_leftindex_by](yield_t &yield) { 
            auto Acoro = L.context(yield); 
            Acoro[oidx] += update_leftindex_by; //inserted_val;
     }); 

    coroutines.emplace_back( 
    [&tio, &R,  &oidx, update_leftindex_by, update_index_by](yield_t &yield) { 
            auto Acoro = R.context(yield); 
            Acoro[oidx] += -(update_leftindex_by + update_index_by);
     }); 

    run_coroutines(tio, coroutines);



    // P[oidx] += z;
    // L[oidx] += zz; 
    // R[oidx] += zzz; 

    return std::make_pair(smallerindex, gt);
}

void MinHeap::initialize(MPCTIO tio, yield_t & yield) {
    auto HeapArray = oram.flat(tio, yield);
    HeapArray.init(0x7fffffffffffff);
}

void MinHeap::initialize_random(MPCTIO tio, yield_t & yield) {
    auto HeapArray = oram.flat(tio, yield);
    std::cout << "initialize_random " << num_items << std::endl;

    std::vector<coro_t> coroutines;

   //  RegAS v[num_items+1];
 
   // for(size_t j = 1; j < num_items; ++j) v[j].ashare = j * tio.player();

   for(size_t j = 1; j <= num_items; ++j)
    {
        coroutines.emplace_back( 
        [&tio, &HeapArray, j](yield_t &yield) { 
                auto Acoro = HeapArray.context(yield); 
                RegAS v;
                v.ashare = j * tio.player();
                Acoro[j] = v; //inserted_val;
         });             
    }    
    
    run_coroutines(tio, coroutines);

//    for(size_t j = 1; j <= num_items; ++j)
    // {
    //     RegAS v;
    //     RegAS inserted_val;
    //     inserted_val.randomize(6);
    //     v.ashare = j * tio.player();
    //     HeapArray[j] = v; //inserted_val;
       
    //     // HeapArray.init([v] (size_t j) { return v; });
    // }
 //   HeapArray.init(0x7fffffffffffff);
}
void MinHeap::print_heap(MPCTIO tio, yield_t & yield) {
    auto HeapArray = oram.flat(tio, yield);
    uint64_t * Pjreconstruction = new uint64_t[num_items + 1];
    for(size_t j = 0; j <= num_items; ++j)  Pjreconstruction[j] = mpc_reconstruct(tio, yield, HeapArray[j]);
    for(size_t j = 0; j <= num_items; ++j) 
    {   
        if(2 * j < num_items) {

          std::cout << j << "-->> HeapArray[" << j << "] = " << std::dec << Pjreconstruction[j] << ", children are: " << Pjreconstruction[2 * j] << " and " << Pjreconstruction[2 * j + 1] <<  std::endl;  
        } 
        else
        {
          std::cout << j << "-->> HeapArray[" << j << "] = " << std::dec << Pjreconstruction[j] << " is a LEAF " <<  std::endl;     
        }
    }
}

auto MinHeap::restore_heap_property_at_root(MPCTIO tio, yield_t & yield, size_t index = 1) {

    //size_t index = 1;
    //std::cout << "index = " << index << std::endl;
    
    auto HeapArray = oram.flat(tio, yield);
    RegAS parent = HeapArray[index];
    RegAS leftchild = HeapArray[2 * index];
    RegAS rightchild = HeapArray[2 * index + 1];
    RegAS sum = parent + leftchild + rightchild;
    CDPF cdpf = tio.cdpf(yield);
    auto[lt, eq, gt] = cdpf.compare(tio, yield, leftchild - rightchild, tio.aes_ops()); // c_1 in the paper
    RegAS smallerchild;
    mpc_select(tio, yield, smallerchild, lt, rightchild, leftchild, 64); // smallerchild holds smaller of left and right child
    auto lteq = lt ^ eq;
    RegXS smallerindex(lt);
    uint64_t leftchildindex = (2 * index);
    uint64_t rightchildindex = (2 * index) + 1;
    smallerindex = (RegXS(lteq) & leftchildindex) ^ (RegXS(gt) & rightchildindex);    
  

    CDPF cdpf0 = tio.cdpf(yield);
    auto[lt1, eq1, gt1] = cdpf0.compare(tio, yield, smallerchild - parent, tio.aes_ops());
    
     auto lt_p_eq_p = lt1 ^ eq1;

    RegBS ltlt1;
    
    mpc_and(tio, yield, ltlt1, lteq, lt_p_eq_p);

    RegAS z, zz;

    run_coroutines(tio, [&tio, &zz, ltlt1, parent, leftchild](yield_t &yield)
            { mpc_flagmult(tio, yield, zz, ltlt1, (parent - leftchild), 64);},
            [&tio, &z, lt1, parent, smallerchild](yield_t &yield)
            {mpc_flagmult(tio, yield, z, lt1, smallerchild - parent, 64);}
            );

    

    HeapArray[index]           += z;
    HeapArray[leftchildindex]  += zz;

    RegAS leftchildplusparent = RegAS(HeapArray[index]) + RegAS(HeapArray[leftchildindex]);
    RegAS tmp = (sum - leftchildplusparent);

    HeapArray[rightchildindex] += tmp - rightchild;

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

    //verify_parent_children_heaps(tio, yield, HeapArray[index], HeapArray[leftchildindex] , HeapArray[rightchildindex]);

    return std::make_pair(smallerindex, gt);
}

    


RegAS MinHeap::extract_min(MPCTIO tio, yield_t & yield, int is_optimized) {

    RegAS minval;
    auto HeapArray = oram.flat(tio, yield);

    minval = HeapArray[1];
    HeapArray[1] = RegAS(HeapArray[num_items]);

    num_items--;

    auto outroot = restore_heap_property_at_root(tio, yield);

    RegXS smaller = outroot.first;

    size_t height = std::log2(num_items);
    
    if(is_optimized > 0)
    {
            typename Duoram < RegAS > ::template OblivIndex < RegXS, 3 > oidx(tio, yield, height);
            oidx.incr(outroot.second);

            for (size_t i = 0; i < height; ++i) {
            auto out = restore_heap_property_optimized(tio, yield, smaller, i + 1, height, typename Duoram < RegAS > ::template OblivIndex < RegXS, 3 > (oidx));;
            smaller = out.first;
            oidx.incr(out.second);
        }
    }


    if(is_optimized == 0)
    {
        for (size_t i = 0; i < height; ++i) {
            smaller = restore_heap_property(tio, yield, smaller);    
        }
    }

    return minval;
}


void MinHeap::heapify2(MPCTIO tio, yield_t & yield, size_t index = 1) {

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

        smaller = restore_heap_property(tio, yield, smaller);    
    }
}

void MinHeap::heapify(MPCTIO tio, yield_t & yield) {

    size_t startIdx = ((num_items + 1) / 2) - 1;

    //std::cout << "startIdx " << startIdx << std::endl;

    for (size_t i = startIdx; i >= 1; i--) {
         heapify2(tio, yield, i);

         //print_heap(tio, yield);
     }
}   

void Heap(MPCIO & mpcio,
    const PRACOptions & opts, char ** args) {
    nbits_t depth     = atoi(args[0]);
    nbits_t depth2    = atoi(args[1]);
    size_t n_inserts  = atoi(args[2]);
    size_t n_extracts = atoi(args[3]);
    int is_optimized  = atoi(args[4]);
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

    //
    std::cout << "items = " << items << std::endl;

    MPCTIO tio(mpcio, 0, opts.num_threads);

    run_coroutines(tio, [ & tio, depth, depth2, items, n_inserts, n_extracts, is_optimized, &mpcio](yield_t & yield) {



        size_t size = size_t(1) << depth;
       // std::cout << "size = " << size << std::endl;

        MinHeap tree(tio.player(), size);
        tree.initialize(tio, yield); 

        tree.num_items = (size_t(1) << depth2) - 1;
        //  std::cout << "num_items " << tree.num_items << std::endl;
        tree.initialize_random(tio, yield);
 


       std::cout << "\n===== Init Stats =====\n";
       


        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
 

        mpcio.reset_stats();
        tio.reset_lamport();


        // tree.heapify(tio, yield);
        // tree.print_heap(tio, yield);

        for (size_t j = 0; j < n_inserts; ++j) {
            RegAS inserted_val;
            inserted_val.randomize(6);
            

            #ifdef VERBOSE
            inserted_val.ashare = inserted_val.ashare;
            uint64_t inserted_val_rec = mpc_reconstruct(tio, yield, inserted_val, 64);
            std::cout << "inserted_val_rec = " << inserted_val_rec << std::endl << std::endl;
            #endif

            if(is_optimized > 0)  tree.insert_optimized(tio, yield, inserted_val);
            if(is_optimized == 0) tree.insert(tio, yield, inserted_val);
            //tree.print_heap(tio, yield);
        }

        std::cout << "\n===== Insert Stats =====\n";
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);

         mpcio.reset_stats();
         tio.reset_lamport();



        // tree.verify_heap_property(tio, yield);
        // tree.print_heap(tio, yield);

        for (size_t j = 0; j < n_extracts; ++j) {
           
            tree.extract_min(tio, yield, is_optimized);
            
            //RegAS minval = tree.extract_min(tio, yield, is_optimized);
            // uint64_t minval_reconstruction = mpc_reconstruct(tio, yield, minval, 64);
            // std::cout << "minval_reconstruction = " << minval_reconstruction << std::endl;
            // tree.verify_heap_property(tio, yield);
            // tree.print_heap(tio, yield);
        }
        std::cout << "\n===== Extract Min Stats =====\n";
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);

        //tree.print_heap(tio, yield);
        //tree.verify_heap_property(tio, yield);
        
      
    });
}