#include <functional>
#include "types.hpp"
#include "duoram.hpp"
#include "cell.hpp"
#include "rdpf.hpp"
#include "shapes.hpp"
#include "heap.hpp"
#include "trie.hpp"


// void TrieClass::insert(MPCTIO &tio, yield_t &yield, RegXS val, RegXS y,unsigned player) {
//     auto TrieArray = oram.flat(tio, yield);
//     num_items++;
//     RegXS b = TrieArray[val];
//     //mpc_reconstruct()
//     RegXS x;
//     value_t check = mpc_reconstruct(tio,yield,b,64);
//     std::cout<< check << " ";
//     // std::cout<< b.xshare <<"  ";     
//     // RegXS z;  
//     // RegXS b_not;
//     // mpc_not(b_not, b, 1);
//     // mpc_and(tio, yield,z, b_not, y, 1);
//     // std::cout<< mpc_reconstruct(tio,yield,z,1)<<"  ";

//     if(check==0)
//     b.xshare = 0;
//     if(check==1)
//     b.xshare=player;
//     mpc_xor_if(tio,yield,x,y,b,player);
//     TrieArray[val] = x; 
// }


void TrieClass::insert(MPCTIO &tio, yield_t &yield, RegXS index, RegXS &insert_value, unsigned player) {
    auto TrieArray = oram.flat(tio, yield);
    num_items++;
    std::cout<< mpc_reconstruct(tio,yield,index,64)<<" " ;

    RegXS b = TrieArray[index];
    //std::cout<<b<<" ";
    
    // Get reconstructed value
    value_t check = mpc_reconstruct(tio, yield, b, 64);
    //std::cout << check << " "<<b.xshare<<" ";
    RegXS input;
    if(check == 0){
        input.xshare = player;
    }
    else{
        input.xshare = 0;
    }
    //std::cout << input.xshare<<" ";
    RegXS x;  
    mpc_xor_if(tio, yield, x, insert_value,b,input, player);

    TrieArray[index] = x;  
}


void TrieClass::search(MPCTIO &tio, yield_t & yield, RegXS index,RegBS &Z,unsigned player){
    auto TrieArray = oram.flat(tio, yield);
    RegXS val =  TrieArray[index];
    RegBS bval;
    if(mpc_reconstruct(tio,yield,val,64)==1){
        bval.bshare = player;
    }
    else bval.bshare = 0;
    //std::cout<<"val reconstruct"<<mpc_reconstruct(tio,yield,val,64);
    RegBS temp;
    mpc_and(tio,yield,temp,bval,Z);
    Z = temp;
    //std::cout<<"Z share "<<Z.bshare<<" z recobstruct "<<mpc_reconstruct(tio,yield,Z)<<" ";
         
}


void TrieClass::init(MPCTIO &tio, yield_t & yield) {
    auto TrieArray = oram.flat(tio, yield);
    TrieArray.init(0x7fffffffffffffff);
    num_items = 0;
}

void TrieClass::init(MPCTIO &tio, yield_t &yield, size_t n) {
    auto TrieArray = oram.flat(tio, yield);  // Get a flat representation of the ORAM.

    num_items = n;  // Store the number of elements in the heap.

    // Initialize all elements with zero
    TrieArray.init([](size_t i) {
        return size_t(0);  // Always return 0 for all elements.
    });
}

void TrieClass::print_trie(MPCTIO &tio, yield_t &yield, size_t size){
    auto HeapArray = oram.flat(tio, yield);
    auto Pjreconstruction = HeapArray.reconstruct();
    for(size_t i = 0 ; i< size ; i++){
        std::cout <<i<<"->"<< Pjreconstruction[i].share()<<"   ";
    }

}

#define HEAP_VERBOSE

int preIndex[10];
int letterToIndex(char x, int pos, int alphasize){
    char ch = x-'a'+1;
    if(pos==0){
        preIndex[pos]=ch;
        return ch;
    }
    else{
        // here 4 is the length of the aplhabet size.
        preIndex[pos] = preIndex[pos-1]*alphasize+ch;
        return preIndex[pos];
    }
}

size_t sumOfPowers(size_t n, size_t m) {
    if (n == 1) return m;  // Special case when n = 1

    return (n * (std::pow(n, m) - 1)) / (n - 1);
}


void Trie(unsigned p,MPCIO & mpcio,  const PRACOptions & opts, char ** args) {


    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);

    int nargs = 0;

    while (args[nargs] != nullptr) {
        ++nargs;
    }
    unsigned player = p;
    int alphasize = 0;
    int triedepth = 0;
    size_t n_inserts = 0;
    size_t n_searches = 0;
    int is_optimized = 0;
    int run_sanity = 0;

    for (int i = 0; i < nargs; i += 2) {
        std::string option = args[i];
        if (option == "-m" && i + 1 < nargs) {
            alphasize = std::atoi(args[i + 1]);
        } else if (option == "-d" && i + 1 < nargs) {
            triedepth = std::atoi(args[i + 1]);
        } else if (option == "-i" && i + 1 < nargs) {
            n_inserts = std::atoi(args[i + 1]);
        } else if (option == "-e" && i + 1 < nargs) {
            n_searches = std::atoi(args[i + 1]);
        } else if (option == "-opt" && i + 1 < nargs) {
            is_optimized = std::atoi(args[i + 1]);
        } else if (option == "-s" && i + 1 < nargs) {
            run_sanity = std::atoi(args[i + 1]);
        }
    }

    run_coroutines(tio, [ & tio, alphasize, triedepth, n_inserts, n_searches, is_optimized, run_sanity,player, &mpcio](yield_t & yield) {
        size_t size = sumOfPowers(alphasize,triedepth)+1;
        std::cout << size<<"\n";
        TrieClass tree(tio.player(), size);
        tree.init(tio, yield, size);
        tree.print_trie(tio, yield,size);
        std::cout << "\n===== Init Stats=====\n";
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();

        std::string insertArray[] = {"dad","aab","aca","dca"};
        std::string searchArray[] = {"ddd","aab","aca","dca"};
        for(size_t i=0;i<n_inserts;i++){

            for(size_t j = 0 ; j < insertArray[i].length() ; j++){
                RegXS share;
                RegXS insert_value;
                insert_value.xshare = 1;
                share.xshare = 1000;

                size_t inserted_index =  letterToIndex(insertArray[i][j],j,alphasize);
                RegXS i_index;
                i_index.xshare = inserted_index;
                //std::cout<< i_index.xshare<<" ";
                //std::cout<< share.xshare<<" ";

                if(player==0){
                    share = i_index^share;
                    insert_value.xshare = 0;
                }
                //std::cout<< mpc_reconstruct(tio,yield,share,64)<<" " ;

                //if(is_optimized > 1)   tree.insert_optimized(tio, yield, share);
                //if(is_optimized == 1)  tree.insert_semi_optimized(tio, yield, share);
                if(is_optimized == 0) tree.insert(tio, yield,share,insert_value,player);

            }
            
            // RegXS input;
            // RegXS output;
            // RegXS cond;
            // if(player == 0) cond.xshare = 0;
            // if(player == 1) cond.xshare = 1;  
            
            // if(player == 0) input.xshare = 0;
            // if(player == 1) input.xshare = 1; 


            // if(player == 0) output.xshare = 0;
            // if(player == 1) output.xshare = 1; 
            
            // mpc_xor_if(tio, yield, input, output, cond, player);

            std::cout << "\ninserted value is " << insertArray[i] << std::endl;
            tree.print_trie(tio, yield,size);
            std::cout<<"\n";
        }        
        
        

        std::cout << "\n===== Insert Stats =====\n";
        tio.sync_lamport();
        mpcio.dump_stats(std::cout);

        mpcio.reset_stats();
        tio.reset_lamport();

        #ifdef HEAP_VERBOSE
        //tree.print_heap(tio, yield);
        #endif
        tree.print_trie(tio,yield,size);
        std::cout<<"\n";
        for(size_t i = 0 ; i< n_searches; i++){
            RegBS Z;
            if(player==0){
                Z.bshare=false;
            }
            else{
                Z.bshare=true;
            }
            //std::cout<<mpc_reconstruct(tio,yield,Z)<< "  "<<  Z.bshare<<"    ";
            for(size_t j = 0;j<searchArray[i].length();j++){
                RegXS share;
                
                share.xshare = 2000;

                size_t inserted_index =  letterToIndex(searchArray[i][j],j,alphasize);
                RegXS i_index;
                i_index.xshare = inserted_index;
                

                if(player==0){
                    share = i_index^share;
                    
                }
                std::cout<<mpc_reconstruct(tio,yield,share,64)<<" ";
                tree.search(tio,yield,share,Z,player);
                               
        }

        //mpc_reconstruct(tio,yield,Z,64);
        if(mpc_reconstruct(tio,yield,Z))
        std::cout << "\nThe value  " << searchArray[i] << " is present" << std::endl;
        else
        std::cout << "\nthe value " << searchArray[i] << " is not present" << std::endl;
       }

    

       std::cout << "\n=====  Search Stats =====\n";
       tio.sync_lamport();
       mpcio.dump_stats(std::cout);

    }
    );
}
