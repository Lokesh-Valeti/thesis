#include <functional>

#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include  "heap.hpp"


 
void HEAP::initialize(int num_players, size_t size) {
    this->MAX_SIZE = size;
    this->num_items = 0;
    oram = new Duoram<RegAS>(num_players, size);
}

RegAS reconstruct_AS(MPCTIO &tio, yield_t &yield, RegAS AS) {
    RegAS peer_AS;
    RegAS reconstructed_AS = AS;
    if (tio.player() == 1) {
        tio.queue_peer(&AS, sizeof(AS));
    } else {
        RegAS peer_AS;
        tio.recv_peer(&peer_AS, sizeof(peer_AS));
         reconstructed_AS += peer_AS;
    }

    yield();

    if (tio.player() == 0) {
        tio.queue_peer(&AS, sizeof(AS));
    } else {
        RegAS peer_flag;
        tio.recv_peer(&peer_AS, sizeof(peer_AS));
        reconstructed_AS += peer_AS;
    }
    
 
   return reconstructed_AS;
}

bool reconstruct_flag(MPCTIO &tio, yield_t &yield, RegBS flag) {
    RegBS peer_flag;
    RegBS reconstructed_flag = flag;
    if (tio.player() == 1) {
        tio.queue_peer(&flag, sizeof(flag));
    } else {
        RegBS peer_flag;
        tio.recv_peer(&peer_flag, sizeof(peer_flag));
        reconstructed_flag ^= peer_flag;
    }

    yield();

    if (tio.player() == 0) {
        tio.queue_peer(&flag, sizeof(flag));
    } else {
        RegBS peer_flag;
        tio.recv_peer(&peer_flag, sizeof(peer_flag));
        reconstructed_flag ^= peer_flag;
    }
    
    return reconstructed_flag.bshare;
}

int HEAP::insert(MPCTIO tio, yield_t &yield, RegAS val) {
  auto HeapArray = oram->flat(tio, yield);
  num_items++;
  std::cout << "num_items = " << num_items << std::endl;
  val.randomize();
  val.ashare = val.ashare/1000;
  std::cout << "we are adding in: " << std::endl;
  reconstruct_AS(tio, yield, val);
  yield();
  size_t child  = num_items;
  size_t parent = child/2;
  std::cout << "child = " << child << std::endl;
  std::cout << "parent = " << parent << std::endl;
  HeapArray[num_items] = val;
  RegAS tmp = HeapArray[num_items];
  reconstruct_AS(tio, yield, tmp);
  yield();
  while(parent != 0) {
    std::cout << "while loop\n";
    RegAS sharechild  = HeapArray[child];
    RegAS shareparent = HeapArray[parent];
    RegAS sharechildrec = reconstruct_AS(tio, yield, sharechild);
    yield();
    RegAS shareparentrec = reconstruct_AS(tio, yield, shareparent);
    yield();
    std::cout << "\nchild reconstruct_AS = \n";
    sharechildrec.dump();
    std::cout << "\nparent reconstruct_AS = \n";
    shareparentrec.dump();
    std::cout << "\n----\n";
    CDPF cdpf = tio.cdpf(yield);
    RegAS diff =  sharechild-shareparent;
    std::cout << "diff = " << std::endl;
    RegAS diff_rec = reconstruct_AS(tio, yield, diff);
    diff_rec.dump();
    std::cout << std::endl << std::endl;
    auto [lt, eq, gt] = cdpf.compare(tio, yield, diff, tio.aes_ops());
    auto lteq = lt ^ eq;
    bool lteq_rec = reconstruct_flag(tio, yield, lteq);
    yield();   
    bool lt_rec = reconstruct_flag(tio, yield, lt);
    yield();
    std::cout <<"lt_rec = " << (int) lt_rec << std::endl;
    std::cout << std::endl;
    bool eq_rec = reconstruct_flag(tio, yield, eq);
    yield();
    std::cout <<"eq_rec = " << (int) eq_rec << std::endl;
    std::cout << std::endl;
    yield();   
    bool gt_rec = reconstruct_flag(tio, yield, gt);
    yield();
    std::cout <<"gt_rec = " << (int) gt_rec << std::endl;
    std::cout << std::endl;
    if(lteq_rec) {
      if(sharechildrec.ashare > shareparentrec.ashare) {
        std::cout << "\nchild reconstruct_AS = \n";
        sharechildrec.dump();
        std::cout << "\nchild reconstruct_AS = \n";
        shareparentrec.dump();
        std::cout << "\n----\n";
      }
      assert(sharechildrec.ashare <= shareparentrec.ashare);
    }
    if(gt_rec) {
      if(sharechildrec.ashare < shareparentrec.ashare) {
        std::cout << "\nchild reconstruct_AS = \n";
        sharechildrec.dump();
        std::cout << "\nchild reconstruct_AS = \n";
        shareparentrec.dump();
        std::cout << "\n----\n";
      }
        assert(sharechildrec.ashare > shareparentrec.ashare);
      }

      std::cout << "child = " << child << std::endl;
      sharechildrec  =  reconstruct_AS(tio, yield, sharechild);
      sharechildrec.dump();
      yield();
      std::cout << "parent = " << parent << std::endl;
      shareparentrec =  reconstruct_AS(tio, yield, shareparent);
      shareparentrec.dump();
      yield();

      std::cout << "\n^^^ before mpc_oswap\n";
      mpc_oswap(tio, yield, sharechild, shareparent, lt, 64);

      HeapArray[child]       = sharechild;
      HeapArray[parent]      = shareparent;

      std::cout << "child = " << child << std::endl;
      sharechildrec  =  reconstruct_AS(tio, yield, sharechild);
      sharechildrec.dump();
      yield();
      std::cout << "parent = " << parent << std::endl;
      shareparentrec =  reconstruct_AS(tio, yield, shareparent);
      shareparentrec.dump();
      yield();
      std::cout << "\n^^^after mpc_oswap\n";
      assert(sharechildrec.ashare >= shareparentrec.ashare);

      std::cout << "we asserted that: \n";
      sharechildrec.dump();
      std::cout << std::endl << " < " << std::endl;
      shareparentrec.dump();
      std::cout << "\n ----- \n";

      child = parent;
      parent = child/2;
    }

    return 1;
}

int HEAP::verify_heap_property(MPCTIO tio, yield_t &yield) {

std::cout << std::endl << std::endl << "verify_heap_property \n\n\n";
auto HeapArray = oram->flat(tio, yield);
 
    RegAS heapreconstruction[num_items];
    for(size_t j = 0; j <= num_items; ++j)
    {
        RegAS tmp = HeapArray[j];
        heapreconstruction[j] = reconstruct_AS(tio, yield, tmp);
        yield();
    }

    for(size_t j = 1; j < num_items/2; ++j)
    {
        if(heapreconstruction[j].ashare > heapreconstruction[2*j].ashare)
        {
            std::cout << "j = " << j << std::endl;
            heapreconstruction[j].dump();
            std::cout << std::endl;
            std::cout << "2*j = " << 2*j << std::endl;
            heapreconstruction[2*j].dump();
            std::cout << std::endl;

        }
        assert(heapreconstruction[j].ashare <= heapreconstruction[2*j].ashare);
    }
 

return 1;
}


void restore_heap_property(MPCTIO tio, yield_t &yield, RegAS& parent, RegAS& leftchild, RegAS& rightchild)
{

  RegAS sum = parent + leftchild + rightchild;
  CDPF cdpf = tio.cdpf(yield);
  auto [lt, eq, gt] = cdpf.compare(tio, yield, leftchild-rightchild, tio.aes_ops());  

  RegAS smallerchild;
  mpc_select(tio, yield, smallerchild, lt, leftchild, rightchild, 64); // smallerchild holds smaller of left and right child 


  CDPF cdpf0 = tio.cdpf(yield);
  auto [lt0, eq0, gt0] = cdpf.compare(tio, yield, smallerchild-parent, tio.aes_ops());

  RegAS smallest;
  mpc_select(tio, yield, smallest, lt0, smallerchild, parent, 64); // smallerchild holds smaller of left and right child   

  parent     = smallest;
  leftchild  = smallerchild;
  rightchild = (sum - smallerchild - smallest);
}

// RegAS HEAP::extract_min(MPCTIO tio, yield_t &yield) {
//     // this->MAX_SIZE = size;
//     // oram = new Duoram<RegXS>(num_players, size);
//     std::cout << "extract_min" << std::endl;
//     std::cout << "num_items = " << num_items << std::endl;
//     // if(num_items==0)
//     //     return RegAS(0);
    

// //    restore_heap_property(tio, yield, A[0], A[1], A[2]);
//     //RegAS minval = A[1];    
//     //A[1] = A[num_items-1];
//     //if(num_items!=0) 
//     {
//         //Delete root

//         //Node zero;
//         // for(size_t j = 0; j < num_items; ++j)   
//         // {
//         //     RegAS val;
//         //     val.set((0+j)*tio.player());
//         //     A[j] = val;
//         // }

//        // num_items--;

//      //    CDPF cdpf = tio.cdpf(yield);
//      //    RegAS tmp = A[10];
//      //    RegAS tmp1 = A[12];

//      //    auto [lt, eq, gt] = cdpf.compare(tio, yield, tmp-tmp1, tio.aes_ops());    
    
//      // //   reconstruct_flag(tio, yield, lt);
//      //    RegAS selected_val;
//      //    mpc_select(tio, yield, selected_val, gt, tmp, tmp1, 64);

//      //    printf("selected_val is: \n");
//      //    reconstruct_AS(tio, yield, selected_val);
       
//      //    yield();
//      //    printf("first option is: \n");
//      //    reconstruct_AS(tio, yield, tmp);
//      //    yield();
//      //    printf("second option is: \n");
//      //    reconstruct_AS(tio, yield, tmp1);
 

//         //return 1; 
//     } 
//      return 1;
// }
 void Heap(MPCIO &mpcio, const PRACOptions &opts, char **args)
{

    nbits_t depth=3;

    if (*args) {
        depth = atoi(*args);
        ++args;
    }
    size_t items = (size_t(1)<<depth)-1;
    if (*args) {
        items = atoi(*args);
        ++args;
    }

       

        MPCTIO tio(mpcio, 0, opts.num_threads);

       
       
        run_coroutines(tio, [&tio, depth, items] (yield_t &yield) {
        size_t size = size_t(1)<<depth;
        std::cout << "size = " << size << std::endl;
        


 
        HEAP tree(tio.player(), size);
 
        
        for(size_t j = 0; j < 100; ++j)
        {
         RegAS inserted_val;
         tree.insert(tio, yield, inserted_val);
         tree.verify_heap_property(tio, yield);
        }

     
      
         //tree.extract_min(tio, yield);
        
    });
}
