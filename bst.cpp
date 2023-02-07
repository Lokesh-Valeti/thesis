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

std::tuple<RegBS, RegBS> compare_keys(Node n1, Node n2, MPCTIO tio, yield_t &yield) {
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

std::tuple<RegXS, RegBS> BST::insert(MPCTIO &tio, yield_t &yield, RegXS ptr,
    const Node &new_node, Duoram<Node>::Flat &A, int TTL, RegBS isDummy) {
    if(TTL==0) {
        RegBS zero;
        return {ptr, zero};
    }

    RegBS isNotDummy = isDummy ^ (tio.player());
    Node cnode = A[ptr];
    // Compare key
    auto [lteq, gt] = compare_keys(cnode, new_node, tio, yield);

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
