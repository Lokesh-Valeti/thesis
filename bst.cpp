#include <functional>

#include "bst.hpp"

// This file demonstrates how to implement custom ORAM wide cell types.
// Such types can be structures of arbitrary numbers of RegAS and RegXS
// fields.  The example here imagines a node of a binary search tree,
// where you would want the key to be additively shared (so that you can
// easily do comparisons), the pointers field to be XOR shared (so that
// you can easily do bit operations to pack two pointers and maybe some
// tree balancing information into one field) and the value doesn't
// really matter, but XOR shared is usually slightly more efficient.

std::tuple<RegBS, RegBS> compare_keys(MPCTIO tio, yield_t &yield, Node n1, Node n2) {
    CDPF cdpf = tio.cdpf(yield);
    auto [lt, eq, gt] = cdpf.compare(tio, yield, n2.key - n1.key, tio.aes_ops());
    RegBS lteq = lt^eq;
    return {lteq, gt};
}

// Assuming pointer of 64 bits is split as:
// - 32 bits Left ptr
// - 32 bits Right ptr
// < Left, Right>

inline RegXS extractLeftPtr(RegXS pointer){ 
    return ((pointer&(0xFFFFFFFF00000000))>>32); 
}

inline RegXS extractRightPtr(RegXS pointer){
    return (pointer&(0x00000000FFFFFFFF)); 
}

inline void setLeftPtr(RegXS &pointer, RegXS new_ptr){ 
    pointer&=(0x00000000FFFFFFFF);
    pointer+=(new_ptr<<32);
}

inline void setRightPtr(RegXS &pointer, RegXS new_ptr){
    pointer&=(0xFFFFFFFF00000000);
    pointer+=(new_ptr);
}

// Pretty-print a reconstructed BST, rooted at node. is_left_child and
// is_right_child indicate whether node is a left or right child of its
// parent.  They cannot both be true, but the root of the tree has both
// of them false.
void BST::pretty_print(const std::vector<Node> &R, value_t node,
    const std::string &prefix = "", bool is_left_child = false,
    bool is_right_child = false)
{
    if (node == 0) {
        // NULL pointer
        if (is_left_child) {
            printf("%s\xE2\x95\xA7\n", prefix.c_str()); // ╧
        } else if (is_right_child) {
            printf("%s\xE2\x95\xA4\n", prefix.c_str()); // ╤
        } else {
            printf("%s\xE2\x95\xA2\n", prefix.c_str()); // ╢
        }
        return;
    }
    const Node &n = R[node];
    value_t left_ptr = extractLeftPtr(n.pointers).xshare;
    value_t right_ptr = extractRightPtr(n.pointers).xshare;
    std::string rightprefix(prefix), leftprefix(prefix),
        nodeprefix(prefix);
    if (is_left_child) {
        rightprefix.append("\xE2\x94\x82"); // │
        leftprefix.append(" ");
        nodeprefix.append("\xE2\x94\x94"); // └
    } else if (is_right_child) {
        rightprefix.append(" ");
        leftprefix.append("\xE2\x94\x82"); // │
        nodeprefix.append("\xE2\x94\x8C"); // ┌
    } else {
        rightprefix.append(" ");
        leftprefix.append(" ");
        nodeprefix.append("\xE2\x94\x80"); // ─
    }
    pretty_print(R, right_ptr, rightprefix, false, true);
    printf("%s\xE2\x94\xA4", nodeprefix.c_str()); // ┤
    n.dump();
    printf("\n");
    pretty_print(R, left_ptr, leftprefix, true, false);
}

bool reconstruct_flag(MPCTIO &tio, yield_t &yield, RegBS flag) {
    RegBS peer_flag;
    RegBS reconstructed_flag = flag;
    tio.queue_peer(&flag, 1);
    yield();
    tio.recv_peer(&peer_flag, 1);
    reconstructed_flag ^= peer_flag;
    //return 1;
    /* 
    //Opt 1:
    if(reconstructed_flag.bshare)
      return 1;
    else
      return 0;
    */
    //Opt 2:
    return reconstructed_flag.bshare;
    
}

void BST::pretty_print(MPCTIO &tio, yield_t &yield) {
    RegXS peer_root;
    RegXS reconstructed_root = root;
    if (tio.player() == 1) {
        tio.queue_peer(&root, sizeof(root));
    } else {
        RegXS peer_root;
        tio.recv_peer(&peer_root, sizeof(peer_root));
        reconstructed_root += peer_root;
    }

    auto A = oram->flat(tio, yield);
    auto R = A.reconstruct();
    if(tio.player()==0) {
        pretty_print(R, reconstructed_root.xshare);
    }
}

// Check the BST invariant of the tree (that all keys to the left are
// less than or equal to this key, all keys to the right are strictly
// greater, and this is true recursively).  Returns a
// tuple<bool,address_t>, where the bool says whether the BST invariant
// holds, and the address_t is the height of the tree (which will be
// useful later when we check AVL trees).
std::tuple<bool, address_t> BST::check_bst(const std::vector<Node> &R,
    value_t node, value_t min_key = 0, value_t max_key = ~0)
{
    if (node == 0) {
        return { true, 0 };
    }
    const Node &n = R[node];
    value_t key = n.key.ashare;
    value_t left_ptr = extractLeftPtr(n.pointers).xshare;
    value_t right_ptr = extractRightPtr(n.pointers).xshare;
    auto [leftok, leftheight ] = check_bst(R, left_ptr, min_key, key);
    auto [rightok, rightheight ] = check_bst(R, right_ptr, key+1, max_key);
    address_t height = leftheight;
    if (rightheight > height) {
        height = rightheight;
    }
    height += 1;
    return { leftok && rightok && key >= min_key && key <= max_key,
        height };
}

void BST::check_bst(MPCTIO &tio, yield_t &yield) {
    auto A = oram->flat(tio, yield);
    auto R = A.reconstruct();

    auto [ ok, height ] = check_bst(R, root.xshare);
    printf("BST structure %s\nBST height = %u\n",
        ok ? "ok" : "NOT OK", height);
}  

void newnode(Node &a) {
    a.key.randomize(8);
    a.pointers.set(0);
    a.value.randomize();
}

void BST::initialize(int num_players, size_t size) {
    this->MAX_SIZE = size;
    oram = new Duoram<Node>(num_players, size);
}


std::tuple<RegXS, RegBS> BST::insert(MPCTIO &tio, yield_t &yield, RegXS ptr,
    const Node &new_node, Duoram<Node>::Flat &A, int TTL, RegBS isDummy) {
    if(TTL==0) {
        RegBS zero;
        return {ptr, zero};
    }

    RegBS isNotDummy = isDummy ^ (tio.player());
    Node cnode = A[ptr];
    // Compare key
    auto [lteq, gt] = compare_keys(tio, yield, cnode, new_node);

    // Depending on [lteq, gt] select the next ptr/index as
    // upper 32 bits of cnode.pointers if lteq
    // lower 32 bits of cnode.pointers if gt 
    RegXS left = extractLeftPtr(cnode.pointers);
    RegXS right = extractRightPtr(cnode.pointers);
    
    RegXS next_ptr;
    mpc_select(tio, yield, next_ptr, gt, left, right, 32);

    CDPF dpf = tio.cdpf(yield);
    size_t &aes_ops = tio.aes_ops();
    // F_z: Check if this is last node on path
    RegBS F_z = dpf.is_zero(tio, yield, next_ptr, aes_ops);
    RegBS F_i;

    // F_i: If this was last node on path (F_z), and isNotDummy insert.
    mpc_and(tio, yield, F_i, (isNotDummy), F_z);
     
    isDummy^=F_i;
    auto [wptr, direction] = insert(tio, yield, next_ptr, new_node, A, TTL-1, isDummy);
    
    RegXS ret_ptr;
    RegBS ret_direction;
    // If we insert here (F_i), return the ptr to this node as wptr
    // and update direction to the direction taken by compare_keys
    mpc_select(tio, yield, ret_ptr, F_i, wptr, ptr);
    //ret_direction = direction + F_p(direction - gt)
    mpc_and(tio, yield, ret_direction, F_i, direction^gt);
    ret_direction^=direction;  

    return {ret_ptr, ret_direction};
}


// Insert(root, ptr, key, TTL, isDummy) -> (new_ptr, wptr, wnode, f_p)
void BST::insert(MPCTIO &tio, yield_t &yield, const Node &node, Duoram<Node>::Flat &A) {
    bool player0 = tio.player()==0;
    // If there are no items in tree. Make this new item the root.
    if(num_items==0) {
        Node zero;
        A[0] = zero;
        A[1] = node;
        (root).set(1*tio.player());
        num_items++;
        //printf("num_items == %ld!\n", num_items);
        return;
    } else {
        // Insert node into next free slot in the ORAM
        int new_id = 1 + num_items;
        int TTL = num_items++;
        A[new_id] = node;
        RegXS new_addr;
        new_addr.set(new_id * tio.player());
        RegBS isDummy;

        //Do a recursive insert
        auto [wptr, direction] = insert(tio, yield, root, node, A, TTL, isDummy);

        //Complete the insertion by reading wptr and updating its pointers
        RegXS pointers = A[wptr].NODE_POINTERS;
        RegXS left_ptr = extractLeftPtr(pointers);
        RegXS right_ptr = extractRightPtr(pointers);
        RegXS new_right_ptr, new_left_ptr;
        mpc_select(tio, yield, new_right_ptr, direction, right_ptr, new_addr);
        if(player0) {
            direction^=1;
        }
        mpc_select(tio, yield, new_left_ptr, direction, left_ptr, new_addr);
        setLeftPtr(pointers, new_left_ptr);
        setRightPtr(pointers, new_right_ptr);
        A[wptr].NODE_POINTERS = pointers;
        //printf("num_items == %ld!\n", num_items);
    } 
}


void BST::insert(MPCTIO &tio, yield_t &yield, Node &node) {
    auto A = oram->flat(tio, yield);
    auto R = A.reconstruct();

    insert(tio, yield, node, A);
    /*
    // To visualize database and tree after each insert:
    if (tio.player() == 0) {
        for(size_t i=0;i<R.size();++i) {
            printf("\n%04lx ", i);
            R[i].dump();
        }
        printf("\n");
    }
    pretty_print(R, 1);
    */
}

/*
// Compute in MPC a | b 
void mpc_or(MPCTIO &tio, yield_t &yield, RegBS &result, RegBS a, RegBS b) {
    int player0 = tio.player();
    if(player0) {
        a^=1;
        b^=1;    
    } 

    mpc_and(tio, yield, result, a, b);
    if(player0)
        result^=1;
}
*/

bool BST::del(MPCTIO &tio, yield_t &yield, RegXS ptr, RegAS del_key,
     Duoram<Node>::Flat &A, RegBS af, RegBS fs, int TTL, 
    del_return &ret_struct) {
    printf("TTL = %d\n", TTL);
    if(TTL==0) {
        //Reconstruct and return af
        bool success = reconstruct_flag(tio, yield, af);
        printf("Reconstructed flag = %d\n", success);
        return success;
    } else {
        bool player0 = tio.player()==0;
        Node node = A[ptr];
        // Compare key

        CDPF cdpf = tio.cdpf(yield);
        auto [lt, eq, gt] = cdpf.compare(tio, yield, del_key - node.key, tio.aes_ops());
        // c is the direction bit for next_ptr 
        // (c=0: go left or c=1: go right)
        RegBS c = gt;
        // lf = local found. We found the key to delete in this level.
        RegBS lf = eq;

        // Depending on [lteq, gt] select the next ptr/index as
        // upper 32 bits of cnode.pointers if lteq
        // lower 32 bits of cnode.pointers if gt 
        RegXS left = extractLeftPtr(node.pointers);
        RegXS right = extractRightPtr(node.pointers);
        
        CDPF dpf = tio.cdpf(yield);
        size_t &aes_ops = tio.aes_ops();
        // Check if left and right children are 0, and compute F_0, F_1, F_2
        RegBS l0 = dpf.is_zero(tio, yield, left, aes_ops);
        RegBS r0 = dpf.is_zero(tio, yield, right, aes_ops);
        RegBS F_0, F_1, F_2;
        // F_0 = l0 & r0
        mpc_and(tio, yield, F_0, l0, r0);
        // F_1 = l0 \xor r0
        F_1 = l0 ^ r0;
        // F_2 = !(F_0 + F_1) (Only 1 of F_0, F_1, and F_2 can be true)
        F_2 = F_0 ^ F_1;
        if(player0)
            F_2^=1;

        // We set next ptr based on c, but we need to handle three 
        // edge cases where we do not go by just the comparison result
        RegXS next_ptr;
        RegBS c_prime;
        // Case 1: found the node here (lf) or we are finding successor (fs)
        // and there is only one child. We traverse down the lone child path.
        RegBS F_c1a, F_c1b, F_c2, F_c3;
        // Case 1a: lf & F_1
        mpc_and(tio, yield, F_c1a, lf, F_1);
        // Case 1b: fs & F_1
        mpc_and(tio, yield, F_c1b, fs, F_1);
        // Set c_prime for Case 1a and 1b
        mpc_select(tio, yield, c_prime, F_c1a, c, l0);
        mpc_select(tio, yield, c_prime, F_c1b, c, l0);

        // s1: shares of 1 bit, s0: shares of 0 bit
        RegBS s1, s0;
        s1.set(tio.player()==1);
        // Case 2: found the node here (lf) and node has both children (F_2)
        // In find successor case, so find inorder successor
        // (Go right and then find leftmost child.)
        mpc_and(tio, yield, F_c2, lf, F_2);
        mpc_select(tio, yield, c_prime, F_c2, c, s1);

        // Case 3: finding successor (fs) and node has both children (F_2)
        // Go left. 
        mpc_and(tio, yield, F_c3, fs, F_2);
        mpc_select(tio, yield, c_prime, F_c3, c, s0);

        // Set next_ptr
        mpc_select(tio, yield, next_ptr, c_prime, left, right, 32);
        
        RegBS af_prime, fs_prime;
        mpc_or(tio, yield, af_prime, af, lf);

        //A[ptr].dump();
        printf("af = %d, lf = %d\n", af.bshare, lf.bshare);
        // If in Case 2, set fs. We are now finding successor
        mpc_or(tio, yield, fs_prime, fs, F_c2); 
        bool key_found = del(tio, yield, next_ptr, del_key, A, af_prime, fs_prime, TTL-1, ret_struct);

        // If we didn't find the key, we can end here.
        if(!key_found)
          return 0;

        printf("TTL = %d\n", TTL);
        /*
        RegBS F_rs;
        // Flag here should be direction (c_prime) and F_r i.e. we need to swap return ptr in,
        // F_r needs to be returned in ret_struct
        mpc_and(tio, yield, F_rs, c_prime, ret_struct.F_r);
        mpc_select(tio, yield, right, F_rs, right, ret_struct.ret_ptr);
        if(player0)
            c_prime^=1; 
        mpc_and(tio, yield, F_rs, c_prime, ret_struct.F_r);
        mpc_select(tio, yield, left, F_rs, left, ret_struct.ret_ptr); 

        // Update the return structure        
        RegBS F_nd, F_ns, F_r, F_rp, F_rp0;
        mpc_or(tio, yield, ret_struct.F_ss, ret_struct.F_ss, F_c2);
        if(player0)
            af^=1; 
        mpc_and(tio, yield, F_nd, lf, af);
        mpc_and(tio, yield, F_ns, fs, F_0);
        // F_r = F_d.(!F_2)
        if(player0)
            F_2^=1;
        // If we have to delete here, and it doesn't have two children we have to
        // update child pointer in parent with the returned pointer
        mpc_and(tio, yield, F_r, F_nd, F_2);
        mpc_and(tio, yield, F_rp, F_nd, F_1);
        mpc_and(tio, yield, F_rp0, F_nd, F_0);

        mpc_select(tio, yield, ret_struct.N_d, F_nd, ret_struct.N_s, ptr);
        mpc_select(tio, yield, ret_struct.N_s, F_ns, ret_struct.N_s, ptr);
        
        mpc_select(tio, yield, ret_struct.ret_ptr, F_rp, ptr, ret_struct.ret_ptr);
        mpc_select(tio, yield, ret_struct.ret_ptr, F_rp0, ptr, s0);
        ret_struct.F_r = F_r;
        */
        return 1;
    }
}


bool BST::del(MPCTIO &tio, yield_t &yield, RegAS del_key) {
    if(num_items==0)
        return 0;
    if(num_items==1) {
        //Delete root
        auto A = oram->flat(tio, yield);
        Node zero;
        A[0] = zero;
        num_items--;
        return 1; 
    } else {
        int TTL = num_items;
        // Flags for already found (af) item to delete and find successor (fs)
        // if this deletion requires a successor swap
        RegBS af;
        RegBS fs;
        del_return ret_struct; 
        auto A = oram->flat(tio, yield);
        int success = del(tio, yield, root, del_key, A, af, fs, TTL, ret_struct); 
        printf ("Success =  %d\n", success); 
        if(!success){
            return 0;
        }
        else{
            //Fix up the actual deletion and succesor swap (if needed) here
            Node del_node = A.reconstruct(A[ret_struct.N_d]);
            Node suc_ptr = A.reconstruct(A[ret_struct.N_s]);
            printf("del_node key = %ld, suc_node key = %ld\n", 
                del_node.key.ashare, suc_ptr.key.ashare); 
            //print("flag_s = %d\n", rec_struct.F_ss);
        }

      return 1;
    }
}

// Now we use the node in various ways.  This function is called by
// online.cpp.
void bst(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
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
        BST tree(tio.player(), size);

        
        Node node; 
        for(size_t i = 1; i<=items; i++) {
          newnode(node);
          node.key.set(i * tio.player());
          tree.insert(tio, yield, node);
        }
       
        tree.pretty_print(tio, yield);
        
        
        RegAS del_key;
        del_key.set(1);
        tree.del(tio, yield, del_key);

        tree.pretty_print(tio, yield);
        /*
        if (depth < 10) {
            //oram.dump();
            auto R = A.reconstruct();
            // Reconstruct the root
            if (tio.player() == 1) {
                tio.queue_peer(&root, sizeof(root));
            } else {
                RegXS peer_root;
                tio.recv_peer(&peer_root, sizeof(peer_root));
                root += peer_root;
            }
            if (tio.player() == 0) {
                for(size_t i=0;i<R.size();++i) {
                    printf("\n%04lx ", i);
                    R[i].dump();
                }
                printf("\n");
                pretty_print(R, root.xshare);
                auto [ ok, height ] = check_bst(R, root.xshare);
                printf("BST structure %s\nBST height = %u\n",
                    ok ? "ok" : "NOT OK", height);
            }
        }
        */ 
    });
}
