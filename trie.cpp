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
    std::cout<< mpc_reconstruct(tio,yield,index,64)<<"  --  " ;

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
    auto StringArray = second_oram.flat(tio,yield);
    num_items = n;  // Store the number of elements in the heap.

    // Initialize all elements with zero
    TrieArray.init([](size_t i) {
        return size_t(0);  // Always return 0 for all elements.
    });
    StringArray.init([](size_t i) {
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

void TrieClass::print_trie_stringcheck(MPCTIO &tio, yield_t &yield, size_t size){
    auto HeapArray = second_oram.flat(tio, yield);
    auto Pjreconstruction = HeapArray.reconstruct();
    for(size_t i = 0 ; i< size ; i++){
        std::cout <<i<<"->"<< Pjreconstruction[i].share()<<"   ";
    }

}

#define HEAP_VERBOSE

int preIndex[10];

int letterToIndex(char x, int pos, int alphasize,int is_optimized){
    int ch;
    if(is_optimized>0){
        ch = x - 'a';
    }
    else{
        ch = x - 'a' + 1;
    }
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

// size_t PowerOfLevel(size_t n, size_t level, is_optimized){
//     return 
// }


size_t Power(size_t x, size_t y) {
    return static_cast<size_t>(std::pow(x, (y+1)));
}


void basic(MPCIO &mpcio, yield_t &yield, int alphasize, int triedepth, size_t n_inserts, size_t n_searches , int is_optimized, unsigned player, MPCTIO &tio){
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
            
            RegXS share;
            for(size_t j = 0; j < insertArray[i].length() ; j++){
                
                RegXS insert_value;
                insert_value.xshare = 1;
                share.xshare = 1000;
                size_t inserted_index =  letterToIndex(insertArray[i][j],j,alphasize,is_optimized);
                RegXS i_index;
                i_index.xshare = inserted_index;

                if(player==0){
                    share = i_index^share;
                    insert_value.xshare = 0;
                }
               
                tree.insert(tio, yield,share,insert_value,player);

                if(is_optimized==1){
                    tree.print_trie(tio,yield,size);
                    std::cout<<"\n";
                }
            }
            RegXS check;
            check.xshare=player;
            //std::cout<<"  ----- "<< mpc_reconstruct(tio, yield,check)<<"  ------ ";

            auto End_String = tree.second_oram.flat(tio, yield);
            End_String[share] = check;

            //std::cout<<"  ----- "<< mpc_reconstruct(tio, yield,End_String[share])<<"  ------ ";

            std::cout << "\ninserted value is " << insertArray[i] << std::endl;
            tree.print_trie(tio, yield,size);
            std::cout<<"\n";
            std::cout<<"String presence array \n";
            tree.print_trie_stringcheck(tio,yield,size);
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
            RegXS share;
            for(size_t j = 0;j<searchArray[i].length();j++){
                
                
                share.xshare = 2000;

                size_t inserted_index =  letterToIndex(searchArray[i][j],j,alphasize,is_optimized);


    
                inserted_index = letterToIndex(searchArray[i][j],j,alphasize,is_optimized);
                RegXS i_index;
                i_index.xshare = inserted_index;
                

                if(player==0){
                    share = i_index^share;
                    
                }
                std::cout<<mpc_reconstruct(tio,yield,share,64)<<" ";
                tree.search(tio,yield,share,Z,player);
                               
        }
        auto End_String = tree.second_oram.flat(tio, yield);
        RegXS check = End_String[share];
        RegBS temp;

        //std::cout<<"  ----- "<< mpc_reconstruct(tio, yield,check)<<"  ------ ";
        // mpc and between check and Z but we are reconstructing the check value because as of now we dont have mpc_ and
        if(mpc_reconstruct(tio, yield,check)==1){
            temp.bshare = player;
        }
        else{
            temp.bshare=0;
        }
        RegBS value;
        mpc_and(tio,yield,value,temp,Z);
        Z = value;

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


void semi_optimized(MPCIO &mpcio, yield_t &yield, int alphasize, int triedepth, size_t n_inserts, size_t n_searches , int is_optimized, unsigned player, MPCTIO &tio){
    TrieClass* trieArray[triedepth];  // Declaring an array of pointers
        for(int i = 0 ; i< triedepth; i++){
            size_t size  = Power(alphasize,i);
            trieArray[i] = new TrieClass(tio.player(),size);

        }

        std::string insertArray[] = {"dad","aab","aca","dca"};
        std::string searchArray[] = {"ddd","aab","aca","dca"};

        for(size_t i=0;i<n_inserts;i++){
            RegXS share;
            size_t j = 0; 
            for(; j < insertArray[i].length() ; j++){

                size_t size  = Power(alphasize,j);
                
                RegXS insert_value;
                insert_value.xshare = 1;
                share.xshare = 1000;
                size_t inserted_index =  letterToIndex(insertArray[i][j],j,alphasize,is_optimized);
                RegXS i_index;
                i_index.xshare = inserted_index;

                if(player==0){
                    share = i_index^share;
                    insert_value.xshare = 0;
                }
               
                trieArray[j]->insert(tio, yield,share,insert_value,player);
                    trieArray[j]->print_trie(tio,yield,size);
                    std::cout<<"\n";

            }
            RegXS check;
            check.xshare=player;
            //std::cout<<"  ----- "<< mpc_reconstruct(tio, yield,check)<<"  ------ ";

            auto End_String = trieArray[j-1]->second_oram.flat(tio, yield);
            End_String[share] = check;

            std::cout << "inserted value is " << insertArray[i] << std::endl;
            
            
        }

        for(size_t i = 0 ; i< n_searches; i++){
            RegBS Z;
            if(player==0){
                Z.bshare=false;
            }
            else{
                Z.bshare=true;
            }
            size_t j = 0;
            RegXS share;
            //std::cout<<mpc_reconstruct(tio,yield,Z)<< "  "<<  Z.bshare<<"    ";
            for(;j<searchArray[i].length();j++){
                
                
                share.xshare = 2000;

                size_t inserted_index =  letterToIndex(searchArray[i][j],j,alphasize,is_optimized);


    
                inserted_index = letterToIndex(searchArray[i][j],j,alphasize,is_optimized);
                RegXS i_index;
                i_index.xshare = inserted_index;
                

                if(player==0){
                    share = i_index^share;
                    
                }
                std::cout<<mpc_reconstruct(tio,yield,share,64)<<" ";
                trieArray[j]->search(tio,yield,share,Z,player);
                               
        }
        auto End_String = trieArray[j-1]->second_oram.flat(tio, yield);
        RegXS check = End_String[share];
        RegBS temp;

        //std::cout<<"  ----- "<< mpc_reconstruct(tio, yield,check)<<"  ------ ";
        // mpc and between check and Z but we are reconstructing the check value because as of now we dont have mpc_ and
        if(mpc_reconstruct(tio, yield,check)==1){
            temp.bshare = player;
        }
        else{
            temp.bshare=0;
        }
        RegBS value;
        mpc_and(tio,yield,value,temp,Z);
        Z = value;

        //mpc_reconstruct(tio,yield,Z,64);
        if(mpc_reconstruct(tio,yield,Z))
        std::cout << "\nThe value  " << searchArray[i] << " is present" << std::endl;
        else
        std::cout << "\nthe value " << searchArray[i] << " is not present" << std::endl;
       }
}


void optimized(MPCIO &mpcio, yield_t &yield, int alphasize, int triedepth, size_t n_inserts, size_t n_searches , int is_optimized, unsigned player, MPCTIO &tio){
    std::string insertArray[] = {"aaa","aab","aca","dca"};
    std::string searchArray[] = {"ddd","aab","aca","dca"};

    size_t max  =  std::pow(alphasize,triedepth);
    std::cout<< "  .."<< max<<"    ...   ";

    size_t current = 1;
    size_t size=0;
    while(current <= max){
        size+=current;
        current=current<<1;
    }
    
    TrieClass tree(tio.player(), size);
    tree.init(tio,yield,size);
    Duoram<RegXS> :: Flat F(tree.oram,tio,yield);

    for(size_t i = 0 ; i < n_inserts; i++){
        RegXS share;
            size_t j = 0;
            size_t tree_layer=1,trie_layer=alphasize,sum=0; 
            for(; j < insertArray[i].length() ; j++){

               // size_t size  = Power(alphasize,j);
                
                RegXS insert_value;
                insert_value.xshare = 1;
                share.xshare = 1000;
                while(tree_layer<trie_layer){
                    sum+=tree_layer;
                    tree_layer=tree_layer<<1;
                }

                size_t inserted_index =  letterToIndex(insertArray[i][j],j,alphasize,is_optimized);
                //std::cout<<"   "<< inserted_index<<"   ";
                std::cout<< "..."<<inserted_index;
                RegXS i_index;
                i_index.xshare = (inserted_index);

                if(player==0){
                    share = i_index^share;
                    insert_value.xshare = 0;
                }

                std::cout<<"  ..  "<<mpc_reconstruct(tio,yield,share,64)<<" ..  ";
                
                
                trie_layer = trie_layer*trie_layer;
                Duoram<RegXS>::Stride S(F,tio,yield,sum,1);
                //auto Trie_Array = tree.oram.Stride(this,tio,yield,sum,1);
                //Trie_Array[share] = insert_value;
                //share.xshare = sum + share.xshare;
                //doubt...needs to be verified
                //S[share]=insert_value;

            }std::cout<<"\n";
            tree.print_trie(tio,yield,size);

    }
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
        
        
        if(is_optimized==0){
            basic(mpcio,yield,alphasize,triedepth,n_inserts,n_searches,is_optimized,player,tio);        
        }
        else if(is_optimized==1){
            semi_optimized(mpcio,yield,alphasize,triedepth,n_inserts,n_searches,is_optimized,player,tio);
        }
        else if(is_optimized==2){
            std::cout<<"hello";
            optimized(mpcio,yield,alphasize,triedepth,n_inserts,n_searches,is_optimized,player,tio);
        }
    }
    
    );
}
