#include <functional>

#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include  "heap.hpp"


 
void HEAP::initialize(int num_players, size_t size) {
    this->MAX_SIZE = size;
    this->num_items = size;
    oram = new Duoram<RegAS>(num_players, size);
}

void reconstruct_AS(MPCTIO &tio, yield_t &yield, RegAS AS) {
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
    
   printf("reconstructed_AS = ");
   reconstructed_AS.dump();
   printf("\n\n");
}

bool reconstruct_flag(MPCTIO &tio, yield_t &yield, RegBS flag) {
    RegBS peer_flag;
    RegBS reconstructed_flag;
    if (tio.player() == 1) {
        tio.queue_peer(&flag, sizeof(flag));
    } else {
        RegBS peer_flag;
        tio.recv_peer(&peer_flag, sizeof(peer_flag));
        reconstructed_flag ^= peer_flag;
    }

    if (tio.player() == 0) {
        tio.queue_peer(&flag, sizeof(flag));
    } else {
        RegBS peer_flag;
        tio.recv_peer(&peer_flag, sizeof(peer_flag));
        reconstructed_flag ^= peer_flag;
    }
    
    return reconstructed_flag.bshare;
}

int HEAP::insert(MPCTIO tio, yield_t &yield, std::vector<RegAS>& A, RegAS val) {

    auto HeapArray = oram->flat(tio, yield);
    num_items++;
    //auto A = oram->flat(tio, yield);
    std::cout << "num_items = " << num_items << std::endl;
 
    val.randomize();
    size_t child  = num_items-1;
    size_t parent = child/2;
    std::cout << "child = " << child << std::endl;
    HeapArray[num_items] = val;
    while(parent != 0)
    {
     CDPF cdpf = tio.cdpf(yield);
 
     auto [lt, eq, gt] = cdpf.compare(tio, yield, A[child]-A[parent], tio.aes_ops());   
     std::cout << "child = " << child << std::endl;
 
     

     mpc_oswap(tio, yield, A[child], A[parent], lt, 64);
 
 
     child = parent;
     parent = parent/2;
    }
    return 1;
}

int HEAP::verify_heap_property(MPCTIO tio, yield_t &yield, std::vector<RegAS> A) {

for(size_t i = 0; i < num_items; ++i)
{

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

RegAS HEAP::extract_min(MPCTIO tio, yield_t &yield, std::vector<RegAS> A) {
    // this->MAX_SIZE = size;
    // oram = new Duoram<RegXS>(num_players, size);
    std::cout << "extract_min" << std::endl;
    std::cout << "num_items = " << num_items << std::endl;
    // if(num_items==0)
    //     return RegAS(0);
    

    restore_heap_property(tio, yield, A[0], A[1], A[2]);
    RegAS minval = A[1];    
    //A[1] = A[num_items-1];
    //if(num_items!=0) 
    {
        //Delete root

        //Node zero;
        // for(size_t j = 0; j < num_items; ++j)   
        // {
        //     RegAS val;
        //     val.set((0+j)*tio.player());
        //     A[j] = val;
        // }

       // num_items--;

     //    CDPF cdpf = tio.cdpf(yield);
     //    RegAS tmp = A[10];
     //    RegAS tmp1 = A[12];

     //    auto [lt, eq, gt] = cdpf.compare(tio, yield, tmp-tmp1, tio.aes_ops());    
    
     // //   reconstruct_flag(tio, yield, lt);
     //    RegAS selected_val;
     //    mpc_select(tio, yield, selected_val, gt, tmp, tmp1, 64);

     //    printf("selected_val is: \n");
     //    reconstruct_AS(tio, yield, selected_val);
       
     //    yield();
     //    printf("first option is: \n");
     //    reconstruct_AS(tio, yield, tmp);
     //    yield();
     //    printf("second option is: \n");
     //    reconstruct_AS(tio, yield, tmp1);
 

        //return 1; 
    } 
     return minval;
}
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
        


        std::vector<RegAS> HeapArray;
        HeapArray.resize(size);
        for(size_t i = 0; i < size; ++i)
        {
          HeapArray[i].dump();
          std::cout << std::endl;  
        }

        HEAP tree(tio.player(), size);
        for(size_t i = 0; i < HeapArray.size(); ++i)
        {
          std::cout << i << " : \n";
          HeapArray[i].dump();
          std::cout << std::endl;  
        }

        std::cout << "\n--------\n\n";
        
        for(size_t j = 0; j < 10; ++j)
        {
         RegAS inserted_val;
         tree.insert(tio, yield, HeapArray, inserted_val);
        }

        std::cout << "\n--------\n\n";

        for(size_t i = 0; i < HeapArray.size(); ++i)
        {
          std::cout << i << " : \n";
          reconstruct_AS(tio, yield, HeapArray[i]);
          yield();
          std::cout << std::endl;  
        }
      
         tree.extract_min(tio, yield, HeapArray);
        // tree.verify_heap_property(tio, yield, HeapArray);
    });
}
