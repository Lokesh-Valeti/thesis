#include <functional>

#include "avl.hpp"

void print_green(std::string line) {
    printf("%s%s%s", KGRN, line.c_str(), KNRM);
}

void print_red(std::string line) {
    printf("%s%s%s", KRED, line.c_str(), KNRM);
}

// Pretty-print a reconstructed BST, rooted at node. is_left_child and
// is_right_child indicate whether node is a left or right child of its
// parent.  They cannot both be true, but the root of the tree has both
// of them false.
void AVL::pretty_print(const std::vector<Node> &R, value_t node,
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
    value_t left_ptr = getAVLLeftPtr(n.pointers).xshare;
    value_t right_ptr = getAVLRightPtr(n.pointers).xshare;
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
    dumpAVL(n);
    printf("\n");
    pretty_print(R, left_ptr, leftprefix, true, false);
}

void AVL::print_oram(MPCTIO &tio, yield_t &yield) {
    auto A = oram.flat(tio, yield);
    auto R = A.reconstruct();

    for(size_t i=0;i<R.size();++i) {
        printf("\n%04lx ", i);
        R[i].dump();
    }
    printf("\n");
}

void AVL::pretty_print(MPCTIO &tio, yield_t &yield) {
    RegXS peer_root;
    RegXS reconstructed_root = root;
    if (tio.player() == 1) {
        tio.queue_peer(&root, sizeof(root));
    } else {
        RegXS peer_root;
        tio.recv_peer(&peer_root, sizeof(peer_root));
        reconstructed_root += peer_root;
    }

    auto A = oram.flat(tio, yield);
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
std::tuple<bool, bool, address_t> AVL::check_avl(const std::vector<Node> &R,
    value_t node, value_t min_key = 0, value_t max_key = ~0)
{
    if (node == 0) {
        return { true, true, 0 };
    }
    const Node &n = R[node];
    value_t key = n.key.ashare;
    value_t left_ptr = getAVLLeftPtr(n.pointers).xshare;
    value_t right_ptr = getAVLRightPtr(n.pointers).xshare;
    auto [leftok, leftavlok, leftheight ] = check_avl(R, left_ptr, min_key, key);
    auto [rightok, rightavlok, rightheight ] = check_avl(R, right_ptr, key+1, max_key);
    address_t height = leftheight;
    if (rightheight > height) {
        height = rightheight;
    }
    height += 1;
    int heightgap = leftheight - rightheight;
    bool avlok = (abs(heightgap)<2);
    //printf("node = %ld, leftok = %d, rightok = %d\n", node, leftok, rightok);
    return { leftok && rightok && key >= min_key && key <= max_key,
        avlok && leftavlok && rightavlok, height};
}

void AVL::check_avl(MPCTIO &tio, yield_t &yield) {
    auto A = oram.flat(tio, yield);
    auto R = A.reconstruct();

    RegXS rec_root = this->root;
    if (tio.player() == 1) {
        tio.queue_peer(&(this->root), sizeof(this->root));
    } else {
        RegXS peer_root;
        tio.recv_peer(&peer_root, sizeof(peer_root));
        rec_root+= peer_root;
    }
    if (tio.player() == 0) {
      auto [ bst_ok, avl_ok, height ] = check_avl(R, rec_root.xshare);
      printf("BST structure %s\nAVL structure %s\nTree height = %u\n",
          bst_ok ? "ok" : "NOT OK", avl_ok ? "ok" : "NOT OK", height);
    }
}


/*
  Rotate: (gp = grandparent (if exists), p = parent, c = child)

  This rotates the p -> c link.

   gp                           gp
     \                            \
      p   --- Left rotate  --->    c
       \                          /
        c                        p


   gp                           gp
     \                            \
      p  --- Right rotate  --->    c
     /                              \
    c                                p

*/
void AVL::rotate(MPCTIO &tio, yield_t &yield, RegXS &gp_pointers, RegXS p_ptr,
      RegXS &p_pointers, RegXS c_ptr, RegXS &c_pointers, RegBS dir_gpp,
      RegBS dir_pc, RegBS isReal, RegBS F_gp) {

    bool player0 = tio.player()==0;
    RegXS gp_left = getAVLLeftPtr(gp_pointers);
    RegXS gp_right = getAVLRightPtr(gp_pointers);
    RegXS p_left = getAVLLeftPtr(p_pointers);
    RegXS p_right = getAVLRightPtr(p_pointers);
    RegXS c_left = getAVLLeftPtr(c_pointers);
    RegXS c_right = getAVLRightPtr(c_pointers);
    RegXS ptr_upd;

    // F_gpp: Flag to update gp -> p link, F_pc: Flag to update p -> c link
    RegBS F_gpp, F_pc_l, F_pc_r, F_gppr, F_gppl;
    // We care about !F_gp. If !F_gp, then we do the gp->p link updates.
    // Otherwise, we do NOT do any updates to gp-> p link;
    // since F_gp==1, implies gp does not exist and parent is root.
    if(player0)
        F_gp^=1;
    mpc_and(tio, yield, F_gpp, F_gp, isReal);

    // i) gp[dir_gpp] <-- c_ptr
    RegBS nt_dir_gpp = dir_gpp;
    if(player0)
        nt_dir_gpp^=1;
    mpc_select(tio, yield, ptr_upd, F_gpp, p_ptr, c_ptr);

    RegBS not_dir_pc_l = dir_pc, not_dir_pc_r = dir_pc;
    if(player0)
        not_dir_pc_r^=1;
    RegXS c_not_dir_pc; //c[!dir_pc]
    // ndpc_right: if not_dir_pc is right
    // ndpc_left: if not_dir_pc is left
    RegBS F_ndpc_right, F_ndpc_left;
    RegBS nt_dir_pc = dir_pc;
    if(player0)
        nt_dir_pc^=1;

    std::vector<coro_t> coroutines;
    coroutines.emplace_back( 
        [&tio, &F_gppr, F_gpp, dir_gpp](yield_t &yield) { 
            mpc_and(tio, yield, F_gppr, F_gpp, dir_gpp); 
    });
    coroutines.emplace_back( 
        [&tio, &F_gppl, F_gpp, nt_dir_gpp](yield_t &yield) { 
            mpc_and(tio, yield, F_gppl, F_gpp, nt_dir_gpp);
    });
    // ii) p[dir_pc] <-- c[!dir_pc] and iii) c[!dir_pc] <-- p_ptr
    coroutines.emplace_back(
        [&tio, &F_ndpc_right, isReal, not_dir_pc_r](yield_t &yield) { 
            mpc_and(tio, yield, F_ndpc_right, isReal, not_dir_pc_r);
    });
    coroutines.emplace_back(
        [&tio, &F_ndpc_left, isReal, not_dir_pc_l](yield_t &yield) {
            mpc_and(tio, yield, F_ndpc_left, isReal, not_dir_pc_l);
    });
    coroutines.emplace_back(
        [&tio, &F_pc_l, dir_pc, isReal](yield_t &yield) { 
            mpc_and(tio, yield, F_pc_l, dir_pc, isReal);
    });
    coroutines.emplace_back(
        [&tio, &F_pc_r, nt_dir_pc, isReal](yield_t &yield) {
            mpc_and(tio, yield, F_pc_r, nt_dir_pc, isReal);
    });
    run_coroutines(tio, coroutines);

    run_coroutines(tio, [&tio, &gp_right, F_gppr, ptr_upd](yield_t &yield)
        { mpc_select(tio, yield, gp_right, F_gppr, gp_right, ptr_upd);},
        [&tio, &gp_left, F_gppl, ptr_upd](yield_t &yield)
        { mpc_select(tio, yield, gp_left, F_gppl, gp_left, ptr_upd);},
        [&tio, &c_not_dir_pc, F_ndpc_right, c_right](yield_t &yield)
        { mpc_select(tio, yield, c_not_dir_pc, F_ndpc_right, c_not_dir_pc, c_right, AVL_PTR_SIZE);});


     //[&tio, &c_not_dir_pc, F_ndpc_left, c_left](yield_t &yield) 
     mpc_select(tio, yield, c_not_dir_pc, F_ndpc_left, c_not_dir_pc, c_left, AVL_PTR_SIZE);


    // ii) p[dir_pc] <-- c[!dir_pc]
    // iii): c[!dir_pc] <-- p_ptr
    run_coroutines(tio, [&tio, &p_left, F_ndpc_right, c_not_dir_pc](yield_t &yield) 
        { mpc_select(tio, yield, p_left, F_ndpc_right, p_left, c_not_dir_pc, AVL_PTR_SIZE);},
        [&tio, &p_right, F_ndpc_left, c_not_dir_pc](yield_t &yield)
        { mpc_select(tio, yield, p_right, F_ndpc_left, p_right, c_not_dir_pc, AVL_PTR_SIZE);},
        [&tio, &ptr_upd, isReal, c_not_dir_pc, p_ptr](yield_t &yield)
        { mpc_select(tio, yield, ptr_upd, isReal, c_not_dir_pc, p_ptr, AVL_PTR_SIZE);});

    run_coroutines(tio, [&tio, &c_left, F_pc_l, ptr_upd](yield_t &yield)
        { mpc_select(tio, yield, c_left, F_pc_l, c_left, ptr_upd, AVL_PTR_SIZE);},
        [&tio, &c_right, F_pc_r, ptr_upd](yield_t &yield)
        { mpc_select(tio, yield, c_right, F_pc_r, c_right, ptr_upd, AVL_PTR_SIZE);});

    setAVLLeftPtr(gp_pointers, gp_left);
    setAVLRightPtr(gp_pointers, gp_right);
    setAVLLeftPtr(p_pointers, p_left);
    setAVLRightPtr(p_pointers, p_right);
    setAVLLeftPtr(c_pointers, c_left);
    setAVLRightPtr(c_pointers, c_right);
}

/*
  In updateBalanceDel, the position of imbalance, and shift direction for both
  cases are inverted, since a bal_upd on a child implies it reduced height.
  If F_rs: (bal_upd & right_child)
    imbalance, bal_l, balanced, bal_r
    And then left shift to get imbalance bit, and new bal_l, bal_r bits
  else if F_ls: (bal_upd & left_child)
    bal_l, balanced, bal_r, imbalance
    And then right shift to get imbalance bit, and new bal_l, bal_r bits

*/
std::tuple<RegBS, RegBS, RegBS, RegBS> AVL::updateBalanceDel(MPCTIO &tio, yield_t &yield,
        RegBS bal_l, RegBS bal_r, RegBS bal_upd, RegBS child_dir) {
    bool player0 = tio.player()==0;
    RegBS s0;
    RegBS F_rs, F_ls, balanced, imbalance;
    RegBS nt_child_dir = child_dir;
    if(player0) {
        nt_child_dir^=1;
    }

    // balanced = is the node currently balanced
    balanced = bal_l ^ bal_r;
    if(player0) {
        balanced^=1;
    }

    //F_ls (Flag left shift) <- child_dir & bal_upd
    //F_rs (Flag right shift) <- !child_dir & bal_upd
    run_coroutines(tio, [&tio, &F_ls, child_dir, bal_upd](yield_t &yield)
        { mpc_and(tio, yield, F_ls, child_dir, bal_upd);},
        [&tio, &F_rs, nt_child_dir, bal_upd](yield_t &yield)
        { mpc_and(tio, yield, F_rs, nt_child_dir, bal_upd);});

    // Left shift if F_ls
    run_coroutines(tio, [&tio, &imbalance, F_ls, bal_l](yield_t &yield)
        { mpc_select(tio, yield, imbalance, F_ls, imbalance, bal_l);},
        [&tio, &bal_l, F_ls, balanced](yield_t &yield)
        { mpc_select(tio, yield, bal_l, F_ls, bal_l, balanced);},
        [&tio, &balanced, F_ls, bal_r](yield_t &yield)
        { mpc_select(tio, yield, balanced, F_ls, balanced, bal_r);},
        [&tio, &bal_r, F_ls, s0](yield_t &yield)
        { mpc_select(tio, yield, bal_r, F_ls, bal_r, s0);});

    // Right shift if F_rs
    run_coroutines(tio, [&tio, &imbalance, F_rs, bal_r](yield_t &yield)
        { mpc_select(tio, yield, imbalance, F_rs, imbalance, bal_r);},
        [&tio, &bal_r, F_rs, balanced](yield_t &yield)
        { mpc_select(tio, yield, bal_r, F_rs, bal_r, balanced);},
        [&tio, &balanced, F_rs, bal_l](yield_t &yield)
        { mpc_select(tio, yield, balanced, F_rs, balanced, bal_l);},
        [&tio, &bal_l, F_rs, s0](yield_t &yield) 
        { mpc_select(tio, yield, bal_l, F_rs, bal_l, s0);});

    // if(bal_upd) and not imbalance bal_upd<-0
    RegBS bu0;
    mpc_and(tio, yield, bu0, bal_upd, balanced);
    mpc_select(tio, yield, bal_upd, bu0, bal_upd, s0);

    // Any bal_upd, propogates all the way up to root
    return {bal_l, bal_r, bal_upd, imbalance};
}


/*
  If F_rs: (bal_upd & right_child)
    bal_l, balanced, bal_r, imbalance
    And then right shift to get imbalance bit, and new bal_l, bal_r bits
  else if F_ls: (bal_upd & left_child)
    imbalance, bal_l, balanced, bal_r
    And then left shift to get imbalance bit, and new bal_l, bal_r bits

*/
std::tuple<RegBS, RegBS, RegBS, RegBS> AVL::updateBalanceIns(MPCTIO &tio, yield_t &yield,
        RegBS bal_l, RegBS bal_r, RegBS bal_upd, RegBS child_dir) {
    bool player0 = tio.player()==0;
    RegBS s1, s0;
    s1.set(tio.player()==1);
    RegBS F_rs, F_ls, balanced, imbalance, nt_child_dir;

    // balanced = is the node currently balanced
    balanced = bal_l ^ bal_r;
    nt_child_dir = child_dir;
    if(player0){
        nt_child_dir^=1;
    }
    if(player0) {
        balanced^=1;
    }
    run_coroutines(tio, [&tio, &F_rs, child_dir, bal_upd](yield_t &yield) 
        { //F_rs (Flag right shift) <- child_dir & bal_upd
          mpc_and(tio, yield, F_rs, child_dir, bal_upd);},
        [&tio, &F_ls, nt_child_dir, bal_upd](yield_t &yield)
        { //F_ls (Flag left shift) <- !child_dir & bal_upd
          mpc_and(tio, yield, F_ls, nt_child_dir, bal_upd);});


    std::vector<coro_t> coroutines;
    // Right shift if child_dir = 1 & bal_upd = 1
    coroutines.emplace_back(
        [&tio, &imbalance, F_rs, bal_r, balanced](yield_t &yield) {
            mpc_select(tio, yield, imbalance, F_rs, imbalance, bal_r);
        });
    coroutines.emplace_back(
        [&tio, &bal_r, F_rs, balanced](yield_t &yield) {
            mpc_select(tio, yield, bal_r, F_rs, bal_r, balanced);
        });
    coroutines.emplace_back(
        [&tio, &balanced, F_rs, bal_l](yield_t &yield) {
            mpc_select(tio, yield, balanced, F_rs, balanced, bal_l);
        });
    coroutines.emplace_back(
        [&tio, &bal_l, F_rs, s0](yield_t &yield) {
            mpc_select(tio, yield, bal_l, F_rs, bal_l, s0);
        });
    run_coroutines(tio, coroutines);
    coroutines.clear();

    // Left shift if child_dir = 0 & bal_upd = 1
    coroutines.emplace_back(
        [&tio, &imbalance, F_ls, bal_l] (yield_t &yield) {
            mpc_select(tio, yield, imbalance, F_ls, imbalance, bal_l);
        });
    coroutines.emplace_back(
        [&tio, &bal_l, F_ls, balanced] (yield_t &yield) {
            mpc_select(tio, yield, bal_l, F_ls, bal_l, balanced);
        });
    coroutines.emplace_back(
        [&tio, &balanced, F_ls, bal_r] (yield_t &yield) {
            mpc_select(tio, yield, balanced, F_ls, balanced, bal_r);
        });
    coroutines.emplace_back(
        [&tio, &bal_r, F_ls, s0](yield_t &yield) {
            mpc_select(tio, yield, bal_r, F_ls, bal_r, s0);
        });
    run_coroutines(tio, coroutines);

    // bal_upd' <- bal_upd ^ imbalance
    RegBS F_bu0;
    mpc_and(tio, yield, F_bu0, bal_upd, balanced);
    mpc_select(tio, yield, bal_upd, F_bu0, bal_upd, s0);
    mpc_select(tio, yield, bal_upd, imbalance, bal_upd, s0);
    return {bal_l, bal_r, bal_upd, imbalance};
}

std::tuple<RegBS, RegBS, RegXS, RegBS> AVL::insert(MPCTIO &tio, yield_t &yield, RegXS ptr, RegXS ins_addr,
    RegAS insert_key, Duoram<Node>::Flat &A, int TTL, RegBS isDummy, avl_insert_return *ret) {
    if(TTL==0) {
        RegBS z;
        return {z, z, z, z};
    }

    RegBS isReal = isDummy ^ (tio.player());
    Node cnode;
    #ifdef OPT_ON
        typename Duoram<Node>::template OblivIndex<RegXS,1> oidx(tio, yield, ptr, MAX_DEPTH);
        cnode = A[oidx];
    #else
        cnode = A[ptr];
    #endif
    RegXS old_pointers = cnode.pointers;

    // Compare key
    auto [lteq, gt] = compare_keys(tio, yield, cnode.key, insert_key);

    // Depending on [lteq, gt] select the next_ptr
    RegXS next_ptr;
    RegXS left = getAVLLeftPtr(cnode.pointers);
    RegXS right = getAVLRightPtr(cnode.pointers);
    RegBS bal_l = getLeftBal(cnode.pointers);
    RegBS bal_r = getRightBal(cnode.pointers);
    /*
    size_t rec_left = reconstruct_RegXS(tio, yield, left);
    size_t rec_right = reconstruct_RegXS(tio, yield, right);
    size_t rec_key = reconstruct_RegAS(tio, yield, cnode.key);
    printf("\n\nKey = %ld\n", rec_key);
    printf("rec_left = %ld, rec_right = %ld\n", rec_left, rec_right);
    */

    mpc_select(tio, yield, next_ptr, gt, left, right, AVL_PTR_SIZE);
    /*
    size_t rec_next_ptr = reconstruct_RegXS(tio, yield, next_ptr);
    printf("rec_next_ptr = %ld\n", rec_next_ptr);
    */

    CDPF dpf = tio.cdpf(yield);
    size_t &aes_ops = tio.aes_ops();
    // F_z: Check if this is last node on path
    RegBS F_z = dpf.is_zero(tio, yield, next_ptr, aes_ops);
    RegBS F_i;

    // F_i: If this was last node on path (F_z), and isReal insert.
    mpc_and(tio, yield, F_i, (isReal), F_z);

    isDummy^=F_i;
    auto [bal_upd, F_gp, prev_node, prev_dir] = insert(tio, yield,
          next_ptr, ins_addr, insert_key, A, TTL-1, isDummy, ret);
    /*
    rec_bal_upd = reconstruct_RegBS(tio, yield, bal_upd);
    rec_F_gp = reconstruct_RegBS(tio, yield, F_gp);
    printf("Insert returns: rec_bal_upd = %d, rec_F_gp = %d\n",
          rec_bal_upd, rec_F_gp);
    size_t rec_ptr = reconstruct_RegXS(tio, yield, ptr);
    printf("\nrec_ptr = %ld\n", rec_ptr);
    */

    // Save insertion pointer and direction
    /*
    mpc_select(tio, yield, ret->i_node, F_i, ret->i_node, ptr, AVL_PTR_SIZE);
    mpc_select(tio, yield, ret->dir_i, F_i, ret->dir_i, gt);
    */
    

    // Update balance
    // If we inserted at this level (F_i), bal_upd is set to 1
    mpc_or(tio, yield, bal_upd, bal_upd, F_i);
    auto [new_bal_l, new_bal_r, new_bal_upd, imbalance] = updateBalanceIns(tio, yield, bal_l, bal_r, bal_upd, gt);

    // Store if this insert triggers an imbalance
    ret->imbalance ^= imbalance;

    std::vector<coro_t> coroutines;
    // Save grandparent pointer
    coroutines.emplace_back(
        [&tio, &ret, F_gp, ptr](yield_t &yield) {
            mpc_select(tio, yield, ret->gp_node, F_gp, ret->gp_node, ptr, AVL_PTR_SIZE);
        });
    coroutines.emplace_back(
        [&tio, &ret, F_gp, gt](yield_t &yield) {
            mpc_select(tio, yield, ret->dir_gpp, F_gp, ret->dir_gpp, gt);
        });
    // Save parent pointer
    coroutines.emplace_back(
        [&tio, &ret, imbalance, ptr](yield_t &yield) {
            mpc_select(tio, yield, ret->p_node, imbalance, ret->p_node, ptr, AVL_PTR_SIZE);
        });
    coroutines.emplace_back(
        [&tio, &ret, imbalance, gt](yield_t &yield) {
            mpc_select(tio, yield, ret->dir_pc, imbalance, ret->dir_pc, gt);
        });
    // Save child pointer
    coroutines.emplace_back(
        [&tio, &ret, imbalance, prev_node](yield_t &yield) {
            mpc_select(tio, yield, ret->c_node, imbalance, ret->c_node, prev_node, AVL_PTR_SIZE);
        });
    coroutines.emplace_back(
        [&tio, &ret, imbalance, prev_dir](yield_t &yield) {
            mpc_select(tio, yield, ret->dir_cn, imbalance, ret->dir_cn, prev_dir);
        });
    run_coroutines(tio, coroutines);


    // Store new_bal_l and new_bal_r for this node
    setLeftBal(cnode.pointers, new_bal_l);
    setRightBal(cnode.pointers, new_bal_r);
    // We have to write the node pointers anyway to resolve balance updates
    RegBS F_ir, F_il;
    run_coroutines(tio, [&tio, &F_ir, F_i, gt](yield_t &yield) 
        { mpc_and(tio, yield, F_ir, F_i, gt); },
        [&tio, &F_il, F_i, lteq](yield_t &yield)
        { mpc_and(tio, yield, F_il, F_i, lteq); });

    run_coroutines(tio, [&tio, &left, F_il, ins_addr](yield_t &yield) 
        { mpc_select(tio, yield, left, F_il, left, ins_addr);},
        [&tio, &right, F_ir, ins_addr](yield_t &yield)
        { mpc_select(tio, yield, right, F_ir, right, ins_addr);});

    setAVLLeftPtr(cnode.pointers, left);
    setAVLRightPtr(cnode.pointers, right);
    #ifdef OPT_ON
        A[oidx].NODE_POINTERS+=(cnode.pointers - old_pointers);
    #else
        A[ptr].NODE_POINTERS = cnode.pointers;
    #endif
    // s0 = shares of 0
    RegBS s0;

    // Update F_gp flag: If there was an imbalance then we set this to store
    // the grandparent node (node in the level above) into the ret_struct
    mpc_select(tio, yield, F_gp, imbalance, s0, imbalance);

    return {new_bal_upd, F_gp, ptr, gt};
}


// Insert(root, ptr, key, TTL, isDummy) -> (new_ptr, wptr, wnode, f_p)
void AVL::insert(MPCTIO &tio, yield_t &yield, const Node &node) {
    bool player0 = tio.player()==0;
    auto A = oram.flat(tio, yield);
    // If there are no items in tree. Make this new item the root.
    if(num_items==0) {
        Node zero;
        A[0] = zero;
        A[1] = node;
        (root).set(1*tio.player());
        num_items++;
        return;
    } else {
        // Insert node into next free slot in the ORAM
        int new_id;
        RegXS insert_address;
        num_items++;
        int TTL = AVL_TTL(num_items);
        bool insertAtEmptyLocation = (numEmptyLocations() > 0);
        if(insertAtEmptyLocation) {
            insert_address = empty_locations.back();
            empty_locations.pop_back();
            A[insert_address] = node;
        } else {
            new_id = num_items;
            A[new_id] = node;
            insert_address.set(new_id * tio.player());
        }

        RegBS isDummy;
        avl_insert_return ret;
        RegAS insert_key = node.key;
        // Recursive insert function
        auto [bal_upd, F_gp, prev_node, prev_dir] = insert(tio, yield, root, 
            insert_address, insert_key, A, TTL, isDummy, &ret);
        /*
        // Debug code
        bool rec_bal_upd, rec_F_gp, ret_dir_pc, ret_dir_cn;
        rec_bal_upd = reconstruct_RegBS(tio, yield, bal_upd);
        rec_F_gp = reconstruct_RegBS(tio, yield, F_gp);
        ret_dir_pc = reconstruct_RegBS(tio, yield, ret.dir_pc);
        ret_dir_cn = reconstruct_RegBS(tio, yield, ret.dir_cn);
        printf("(Top level) Insert returns: rec_bal_upd = %d, rec_F_gp = %d\n",
              rec_bal_upd, rec_F_gp);
        printf("(Top level) Insert returns: ret.dir_pc = %d, rt.dir_cn = %d\n",
              ret_dir_pc, ret_dir_cn);
        */

        // Perform balance procedure
        RegXS gp_pointers, parent_pointers, child_pointers;
        #ifdef OPT_ON
            typename Duoram<Node>::template OblivIndex<RegXS,1> oidx_gp(tio, yield, ret.gp_node, TTL+1);
            typename Duoram<Node>::template OblivIndex<RegXS,1> oidx_p(tio, yield, ret.p_node, TTL+1);
            typename Duoram<Node>::template OblivIndex<RegXS,1> oidx_c(tio, yield, ret.c_node, TTL+1); 
            gp_pointers = A[oidx_gp].NODE_POINTERS;
            parent_pointers = A[oidx_p].NODE_POINTERS;
            child_pointers = A[oidx_c].NODE_POINTERS;
        #else
            RegXS gp_pointers = A[ret.gp_node].NODE_POINTERS;
            RegXS parent_pointers = A[ret.p_node].NODE_POINTERS;
            RegXS child_pointers = A[ret.c_node].NODE_POINTERS;
        #endif
        // n_node (child's next node)
        RegXS child_left = getAVLLeftPtr(child_pointers);
        RegXS child_right = getAVLRightPtr(child_pointers);
        RegXS n_node, n_pointers;
        mpc_select(tio, yield, n_node, ret.dir_cn, child_left, child_right, AVL_PTR_SIZE);

        #ifdef OPT_ON
            typename Duoram<Node>::template OblivIndex<RegXS,1> oidx_n(tio, yield, n_node, TTL+1);
            n_pointers = A[oidx_n].NODE_POINTERS;
        #else
            n_pointers = A[n_node].NODE_POINTERS;  
        #endif

        RegXS old_gp_pointers, old_parent_pointers, old_child_pointers, old_n_pointers;
        #ifdef OPT_ON
            old_gp_pointers = gp_pointers;
            old_parent_pointers = parent_pointers;
            old_child_pointers = child_pointers;
            old_n_pointers = n_pointers;
        #endif

        // F_dr = (dir_pc != dir_cn) : i.e., double rotation case if
        // (parent->child) and (child->new_node) are not in the same direction
        RegBS F_dr = (ret.dir_pc) ^ (ret.dir_cn);

        /* Flags:  F_cn_rot = child->node rotate
                   F_ur = update root.

           In case of an imbalance we have to always rotate p->c link. (L or R case)
           In case of an imbalance where p->c and c->n are in different directions, we have
           to perform a double rotation (LR or RL case). In such cases, first rotate
           c->n link, and then p->c link
           (Note: in the second rotation c is actually n, since the the first rotation
            swaps their positions)
        */
        RegBS F_cn_rot, F_ur, s0;
        run_coroutines(tio, [&tio, &F_ur, F_gp, ret](yield_t &yield) 
            {mpc_and(tio, yield, F_ur, F_gp, ret.imbalance);},
            [&tio, &F_cn_rot, ret, F_dr](yield_t &yield) 
            {mpc_and(tio, yield, F_cn_rot, ret.imbalance, F_dr);});

        // Get the n children information for 2nd rotate fix before rotations happen.
        RegBS n_bal_l, n_bal_r;
        RegXS n_l = getAVLLeftPtr(n_pointers);
        RegXS n_r = getAVLRightPtr(n_pointers);
        n_bal_l = getLeftBal(n_pointers);
        n_bal_r = getRightBal(n_pointers);

        // First rotation: c->n link
        rotate(tio, yield, parent_pointers, ret.c_node, child_pointers, n_node,
            n_pointers, ret.dir_pc, ret.dir_cn, F_cn_rot, s0);

        // If F_cn_rot, i.e. we did first rotation. Then c and n need to swap before the second rotate.
        RegXS new_child_pointers, new_child;

        run_coroutines(tio, [&tio, &new_child_pointers, F_cn_rot, child_pointers, n_pointers] (yield_t &yield)
            {mpc_select(tio, yield, new_child_pointers, F_cn_rot, child_pointers, n_pointers);}, 
            [&tio, &new_child, F_cn_rot, ret, n_node](yield_t &yield)
            {mpc_select(tio, yield, new_child, F_cn_rot, ret.c_node, n_node, AVL_PTR_SIZE);});

        // Second rotation: p->c link
        rotate(tio, yield, gp_pointers, ret.p_node, parent_pointers, new_child,
            new_child_pointers, ret.dir_gpp, ret.dir_pc, ret.imbalance, F_gp);

        // Set parent and child balances to 0 if there was an imbalance.
        // parent balances are already set to 0 from updateBalanceIns
        RegBS temp_bal, p_bal_l, p_bal_r, p_bal_ndpc;
        RegBS c_bal_l, c_bal_r, c_bal_dpc, n_bal_dpc, n_bal_ndpc;
        p_bal_l = getLeftBal(parent_pointers);
        p_bal_r = getRightBal(parent_pointers);

        run_coroutines(tio, [&tio, &child_pointers, F_cn_rot, new_child_pointers] (yield_t &yield)
            {mpc_select(tio, yield, child_pointers, F_cn_rot, new_child_pointers, child_pointers);},
            [&tio, &n_pointers, F_cn_rot, new_child_pointers] (yield_t &yield) 
            {mpc_select(tio, yield, n_pointers, F_cn_rot, n_pointers, new_child_pointers);});

        c_bal_l = getLeftBal(child_pointers);
        c_bal_r = getRightBal(child_pointers);

        run_coroutines(tio, [&tio, &c_bal_l, ret, s0] (yield_t &yield)
            {mpc_select(tio, yield, c_bal_l, ret.imbalance, c_bal_l, s0);},
            [&tio, &c_bal_r, ret, s0] (yield_t &yield) 
            {mpc_select(tio, yield, c_bal_r, ret.imbalance, c_bal_r, s0);});

        /*   In the double rotation case: balance of c and p have a tweak
           p_bal_ndpc <- !(n_bal_ndpc)
           c_bal_dpc <- !(n_bal_dpc) */
        CDPF cdpf = tio.cdpf(yield);
        size_t &aes_ops = tio.aes_ops();
      
        RegBS n_l0, n_r0;
        run_coroutines(tio, [&tio, &n_l0, n_l, &cdpf, &aes_ops] (yield_t &yield) 
            {n_l0 = cdpf.is_zero(tio, yield, n_l, aes_ops);},
            [&tio, &n_r0, n_r, &cdpf, &aes_ops] (yield_t &yield)
            {n_r0 = cdpf.is_zero(tio, yield, n_r, aes_ops);});

        RegBS p_c_update, n_has_children;
        // n_has_children = !(n_l0 & n_r0)
        mpc_and(tio, yield, n_has_children, n_l0, n_r0);
        if(player0) {
            n_has_children^=1;
        }
      
        run_coroutines(tio, [&tio, &p_c_update, F_cn_rot, n_has_children] (yield_t &yield) 
            {mpc_and(tio, yield, p_c_update, F_cn_rot, n_has_children);}, 
            [&tio, &n_bal_ndpc, ret, n_bal_l, n_bal_r] (yield_t &yield)
            {mpc_select(tio, yield, n_bal_ndpc, ret.dir_pc, n_bal_r, n_bal_l);},
            [&tio, &n_bal_dpc, ret, n_bal_l, n_bal_r] (yield_t &yield)
            {mpc_select(tio, yield, n_bal_dpc, ret.dir_pc, n_bal_l, n_bal_r);},
            [&tio, &p_bal_ndpc, ret, p_bal_r, p_bal_l] (yield_t &yield)
            {mpc_select(tio, yield, p_bal_ndpc, ret.dir_pc, p_bal_r, p_bal_l);});

        // !n_bal_ndpc, !n_bal_dpc
        if(player0) {
            n_bal_ndpc^=1;
            n_bal_dpc^=1;
        }

        run_coroutines(tio, [&tio, &p_bal_ndpc, p_c_update, n_bal_ndpc] (yield_t &yield)
            {mpc_select(tio, yield, p_bal_ndpc, p_c_update, p_bal_ndpc, n_bal_ndpc);},
            [&tio, &c_bal_dpc, p_c_update, n_bal_dpc] (yield_t &yield)
            {mpc_select(tio, yield, c_bal_dpc, p_c_update, c_bal_dpc, n_bal_dpc);});

        std::vector<coro_t> coroutines;
        coroutines.emplace_back([&tio, &p_bal_r, ret, p_bal_ndpc] (yield_t &yield)
            {mpc_select(tio, yield, p_bal_r, ret.dir_pc, p_bal_ndpc, p_bal_r);});
        coroutines.emplace_back([&tio, &p_bal_l, ret, p_bal_ndpc] (yield_t &yield)
            {mpc_select(tio, yield, p_bal_l, ret.dir_pc, p_bal_l, p_bal_ndpc);});
        coroutines.emplace_back([&tio, &c_bal_r, ret, c_bal_dpc] (yield_t &yield)
            {mpc_select(tio, yield, c_bal_r, ret.dir_pc, c_bal_r, c_bal_dpc);});
        coroutines.emplace_back([&tio, &c_bal_l, ret, c_bal_dpc] (yield_t &yield)
            {mpc_select(tio, yield, c_bal_l, ret.dir_pc, c_bal_dpc, c_bal_l);});
        // If double rotation (LR/RL) case, n ends up with 0 balance.
        // In all other cases, n's balance remains unaffected by rotation during insertion.
        coroutines.emplace_back([&tio, &n_bal_l, F_cn_rot, s0] (yield_t &yield)
            {mpc_select(tio, yield, n_bal_l, F_cn_rot, n_bal_l, s0);});
        coroutines.emplace_back([&tio, &n_bal_r, F_cn_rot, s0] (yield_t &yield)
            {mpc_select(tio, yield, n_bal_r, F_cn_rot, n_bal_r, s0);});
        run_coroutines(tio, coroutines);

        setLeftBal(parent_pointers, p_bal_l);
        setRightBal(parent_pointers, p_bal_r);
        setLeftBal(child_pointers, c_bal_l);
        setRightBal(child_pointers, c_bal_r);
        setLeftBal(n_pointers, n_bal_l);
        setRightBal(n_pointers, n_bal_r);

        // Write back update pointers and balances into gp, p, c, and n
        #ifdef OPT_ON
            A[oidx_n].NODE_POINTERS+=(n_pointers - old_n_pointers);
            A[oidx_c].NODE_POINTERS+=(child_pointers - old_child_pointers); 
            A[oidx_p].NODE_POINTERS+=(parent_pointers - old_parent_pointers); 
            A[oidx_gp].NODE_POINTERS+=(gp_pointers - old_gp_pointers); 
        #else
            A[ret.c_node].NODE_POINTERS = child_pointers;
            A[ret.p_node].NODE_POINTERS = parent_pointers;
            A[ret.gp_node].NODE_POINTERS = gp_pointers;
            A[n_node].NODE_POINTERS = n_pointers;
        #endif

        // Handle root pointer update (if F_ur is true)
        // If F_ur and we did a double rotation: root <-- new node
        // If F_ur and we did a single rotation: root <-- child node
        RegXS temp_root = root;
        run_coroutines(tio, [&tio, &temp_root, F_ur, ret] (yield_t &yield)
            {mpc_select(tio, yield, temp_root, F_ur, temp_root, ret.c_node, AVL_PTR_SIZE);},
            [&tio, &F_ur, F_gp, F_dr] (yield_t &yield)
            {mpc_and(tio, yield, F_ur, F_gp, F_dr);});

        mpc_select(tio, yield, temp_root, F_ur, temp_root, n_node, AVL_PTR_SIZE);
        root = temp_root;
    }
}


bool AVL::lookup(MPCTIO &tio, yield_t &yield, RegXS ptr, RegAS key, Duoram<Node>::Flat &A,
    int TTL, RegBS isDummy, Node *ret_node) {
    if(TTL==0) {
        // Reconstruct and return isDummy
        // If we found the key, then isDummy will be true
        bool found = reconstruct_RegBS(tio, yield, isDummy);
        return found;
    }

    RegBS isNotDummy = isDummy ^ (tio.player());
    Node cnode = A[ptr];
    // Compare key
    CDPF cdpf = tio.cdpf(yield);
    auto [lt, eq, gt] = cdpf.compare(tio, yield, key - cnode.key, tio.aes_ops());

    // Depending on [lteq, gt] select the next ptr/index as
    // upper 32 bits of cnode.pointers if lteq
    // lower 32 bits of cnode.pointers if gt
    RegXS left = getAVLLeftPtr(cnode.pointers);
    RegXS right = getAVLRightPtr(cnode.pointers);

    RegXS next_ptr;
    mpc_select(tio, yield, next_ptr, gt, left, right, 32);

    RegBS F_found;
    // If we haven't found the key yet, and the lookup matches the current node key,
    // then we found the node to return
    mpc_and(tio, yield, F_found, isNotDummy, eq);
    mpc_select(tio, yield, ret_node->key, eq, ret_node->key, cnode.key);
    mpc_select(tio, yield, ret_node->value, eq, ret_node->value, cnode.value);

    isDummy^=F_found;
    bool found = lookup(tio, yield, next_ptr, key, A, TTL-1, isDummy, ret_node);

    return found;
}

bool AVL::lookup(MPCTIO &tio, yield_t &yield, RegAS key, Node *ret_node) {
    auto A = oram.flat(tio, yield);
    RegBS isDummy;
    bool found = lookup(tio, yield, root, key, A, num_items, isDummy, ret_node);
    return found;
}

void AVL::updateChildPointers(MPCTIO &tio, yield_t &yield, RegXS &left, RegXS &right,
          RegBS c_prime, avl_del_return ret_struct) {
    bool player0 = tio.player()==0;
    RegBS F_rr; // Flag to resolve F_r by updating right child ptr
    RegBS F_rl; // Flag to resolve F_r by updating left child ptr
    RegBS nt_c_prime = c_prime;
    if(player0)
        nt_c_prime^=1;

    run_coroutines(tio, [&tio, &F_rr, c_prime, ret_struct](yield_t &yield)
        { mpc_and(tio, yield, F_rr, c_prime, ret_struct.F_r);},
        [&tio, &F_rl, nt_c_prime, ret_struct](yield_t &yield)
        { mpc_and(tio, yield, F_rl, nt_c_prime, ret_struct.F_r);});

    run_coroutines(tio, [&tio, &right, F_rr, ret_struct](yield_t &yield)
        { mpc_select(tio, yield, right, F_rr, right, ret_struct.ret_ptr);},
        [&tio, &left, F_rl, ret_struct](yield_t &yield)
        { mpc_select(tio, yield, left, F_rl, left, ret_struct.ret_ptr);});
}


// Perform rotations if imbalance (else dummy rotations)
/*
   For capturing both the symmetric L and R cases of rotations, we'll capture directions with
   dpc  = dir_pc = direction from parent to child, and
   ndpc = not(dir_pc)
   When we travelled down the stack, we went from p->c. But in deletions to handle any imbalance
   we look at c's sibling cs (child's sibling). And the rotation is between p and cs if there
   was an imbalance at p, and perhaps even cs and it's child (the child in dir_pc, as that's the
   only case that results in a double rotation when deleting).

   In case of an imbalance we have to always rotate p->cs link. (L or R case)
   If cs.bal_(dir_pc), then we have a double rotation (LR or RL) case.
   In such cases, first rotate cs->gcs link, and then p->cs link. gcs = grandchild on cs path

   Layout: In the R (or LR) case:

         p
       /   \
      cs    c
     /  \
    a   gcs
        /  \
       x    y

   - One of x or y must exist for it to be an LR case,
     since then cs.bal_(dir_pc) = cs.bal_r = 1

   Layout: In the L (or RL) case:

         p
       /   \
      c     cs
           /  \
         gcs   a
        /   \
       x     y

   - One of x or y must exist for it to be an RL case,
     since then cs.bal_(dir_pc) = cs.bal_l = 1

   (Note: if double rotation case, in the second rotation cs is actually gcs,
    since the the first rotation swaps their positions)
*/

void AVL::fixImbalance(MPCTIO &tio, yield_t &yield, Duoram<Node>::Flat &A, 
        Duoram<Node>::OblivIndex<RegXS, 1> oidx, RegXS oidx_oldptrs, RegXS ptr, 
        RegXS nodeptrs, RegBS new_p_bal_l, RegBS new_p_bal_r, RegBS &bal_upd, 
        RegBS c_prime, RegXS cs_ptr, RegBS imb, RegBS &F_ri, 
        avl_del_return &ret_struct) {
    bool player0 = tio.player()==0;
    RegBS s0, s1;
    s1.set(tio.player()==1);

    RegXS old_cs_ptr;
    Node cs_node;
    #ifdef OPT_ON
        typename Duoram<Node>::template OblivIndex<RegXS,1> oidx_cs(tio, yield, cs_ptr, MAX_DEPTH);
        cs_node = A[oidx_cs];
        old_cs_ptr = cs_node.pointers;
    #else
        cs_node = A[cs_ptr];
    #endif
    //dirpc = dir_pc = dpc = c_prime
    RegBS cs_bal_l, cs_bal_r, cs_bal_dpc, cs_bal_ndpc, F_dr, not_c_prime;
    RegXS gcs_ptr, cs_left, cs_right, cs_dpc, cs_ndpc, null;
    // child's sibling node's balances in dir_pc (dpc), and not_dir_pc (ndpc)
    cs_bal_l = getLeftBal(cs_node.pointers);
    cs_bal_r = getRightBal(cs_node.pointers);
    cs_left = getAVLLeftPtr(cs_node.pointers);
    cs_right = getAVLRightPtr(cs_node.pointers);

    run_coroutines(tio, [&tio, &cs_bal_dpc, c_prime, cs_bal_l, cs_bal_r](yield_t &yield)
        { mpc_select(tio, yield, cs_bal_dpc, c_prime, cs_bal_l, cs_bal_r);},
        [&tio, &cs_bal_ndpc, c_prime, cs_bal_r, cs_bal_l](yield_t &yield)
        { mpc_select(tio, yield, cs_bal_ndpc, c_prime, cs_bal_r, cs_bal_l);},
        [&tio, &cs_dpc, c_prime, cs_left, cs_right](yield_t &yield)
        { mpc_select(tio, yield, cs_dpc, c_prime, cs_left, cs_right);},
        [&tio, &cs_ndpc, c_prime, cs_right, cs_left](yield_t &yield)
        { mpc_select(tio, yield, cs_ndpc, c_prime, cs_right, cs_left);});

    // We need to double rotate (LR or RL case) if cs_bal_dpc is 1
    run_coroutines(tio, [&tio, &F_dr, imb, cs_bal_dpc] (yield_t &yield)
        { mpc_and(tio, yield, F_dr, imb, cs_bal_dpc);},
        [&tio, &gcs_ptr, cs_bal_dpc, cs_ndpc, cs_dpc](yield_t &yield)
        { mpc_select(tio, yield, gcs_ptr, cs_bal_dpc, cs_ndpc, cs_dpc, AVL_PTR_SIZE);});

    Node gcs_node;
    RegXS old_gcs_ptr;
    #ifdef OPT_ON
        typename Duoram<Node>::template OblivIndex<RegXS,1> oidx_gcs(tio, yield, gcs_ptr, MAX_DEPTH);
        gcs_node = A[oidx_gcs];
        old_gcs_ptr = gcs_node.pointers;
    #else
        gcs_node = A[gcs_ptr];
    #endif

    not_c_prime = c_prime;
    if(player0) {
      not_c_prime^=1;
    }
    // First rotation: cs->gcs link
    rotate(tio, yield, nodeptrs, cs_ptr, cs_node.pointers, gcs_ptr,
        gcs_node.pointers, not_c_prime, c_prime, F_dr, s0);

    // If F_dr, we did first rotation. Then cs and gcs need to swap before the second rotate.
    RegXS new_cs_pointers, new_cs, new_ptr;
    run_coroutines(tio, [&tio, &new_cs_pointers, F_dr, cs_node, gcs_node](yield_t &yield)
        { mpc_select(tio, yield, new_cs_pointers, F_dr, cs_node.pointers, gcs_node.pointers);},
        [&tio, &new_cs, F_dr, cs_ptr, gcs_ptr](yield_t &yield)
        { mpc_select(tio, yield, new_cs, F_dr, cs_ptr, gcs_ptr, AVL_PTR_SIZE);},
        [&tio, &new_ptr, F_dr, cs_ptr, gcs_ptr](yield_t &yield) 
        { mpc_select(tio, yield, new_ptr, F_dr, cs_ptr, gcs_ptr);});

    // Second rotation: p->cs link
    // Since we don't have access to gp node here we just send a null and s0
    // for gp_pointers and dir_gpp. Instead this pointer fix is handled by F_r
    // and ret_struct.ret_ptr.
    rotate(tio, yield, null, ptr, nodeptrs, new_cs,
        new_cs_pointers, s0, not_c_prime, imb, s1);

    // If imb (we do some rotation), then update F_r, and ret_ptr, to
    // fix the gp->p link (The F_r clauses later, and this are mutually
    // exclusive events. They will never trigger together.)

    std::vector<coro_t> coroutines;
    coroutines.emplace_back([&tio, &F_ri, imb, s0, s1](yield_t &yield) {
        mpc_select(tio, yield, F_ri, imb, s0, s1);
    });
    coroutines.emplace_back([&tio, &ret_struct, imb, new_ptr](yield_t &yield) {
        mpc_select(tio, yield, ret_struct.ret_ptr, imb, ret_struct.ret_ptr, new_ptr);
    });
    // Write back new_cs_pointers correctly to (cs_node/gcs_node).pointers
    // and then balance the nodes
    coroutines.emplace_back([&tio, &cs_node, F_dr, new_cs_pointers](yield_t &yield) {
        mpc_select(tio, yield, cs_node.pointers, F_dr, new_cs_pointers, cs_node.pointers);
    });
    coroutines.emplace_back([&tio, &gcs_node, F_dr, new_cs_pointers](yield_t &yield) {
        mpc_select(tio, yield, gcs_node.pointers, F_dr, gcs_node.pointers, new_cs_pointers);
    });
    run_coroutines(tio, coroutines);
    coroutines.clear();
    
    /*
       Update balances based on imbalance and type of rotations that happen.
       In the case of an imbalance, updateBalance() sets bal_l and bal_r of p to 0.

    */
    RegBS IC1, IC2, IC3; // Imbalance Case 1, 2 or 3
    run_coroutines(tio, [&tio, &IC1, imb, cs_bal_ndpc] (yield_t &yield) {
        // IC1 = Single rotation (L/R). L/R = dpc
        mpc_and(tio, yield, IC1, imb, cs_bal_ndpc);
        },
        [&tio, &IC3, imb, cs_bal_dpc](yield_t &yield) {
        // IC3 = Double rotation (LR/RL). 1st rotate direction = ndpc, 2nd direction = dpc
        mpc_and(tio, yield, IC3, imb, cs_bal_dpc);
        });

    // IC2 = Single rotation (L/R).
    IC2 = IC1 ^ IC3;
    if(player0) {
      IC2^=1;
    }

    RegBS p_bal_dpc, p_bal_ndpc;
    RegBS IC2_ndpc_l, IC2_ndpc_r, IC2_dpc_l, IC2_dpc_r;

    run_coroutines(tio, [&tio, &IC2, imb] (yield_t &yield)
        { mpc_and(tio, yield, IC2, imb, IC2);},
        [&tio, &cs_bal_dpc, imb, s0](yield_t &yield) 
        { // IC1, IC2, IC3: CS.bal = 0 0
          mpc_select(tio, yield, cs_bal_dpc, imb, cs_bal_dpc, s0);},
        [&tio, &cs_bal_ndpc, c_prime, imb, s0](yield_t &yield) {
        mpc_select(tio, yield, cs_bal_ndpc, imb, cs_bal_ndpc, s0);});

    run_coroutines(tio, [&tio, &cs_bal_r, c_prime, cs_bal_ndpc, cs_bal_dpc] (yield_t &yield) 
        { mpc_select(tio, yield, cs_bal_r, c_prime, cs_bal_ndpc, cs_bal_dpc);},
        [&tio, &cs_bal_l, c_prime, cs_bal_dpc, cs_bal_ndpc](yield_t &yield)
        { mpc_select(tio, yield, cs_bal_l, c_prime, cs_bal_dpc, cs_bal_ndpc);});

    // IC2: p.bal_ndpc = 1, cs.bal_dpc = 1
    // (IC2 & not_c_prime)
    coroutines.emplace_back([&tio, &p_bal_ndpc, c_prime, new_p_bal_r, new_p_bal_l](yield_t &yield)
        { mpc_select(tio, yield, p_bal_ndpc, c_prime, new_p_bal_r, new_p_bal_l);});
    coroutines.emplace_back([&tio, &IC2_ndpc_l, c_prime, IC2] (yield_t &yield)
        { mpc_and(tio, yield, IC2_ndpc_l, IC2, c_prime);});
    coroutines.emplace_back([&tio, &IC2_ndpc_r, IC2, not_c_prime](yield_t &yield)
        { mpc_and(tio, yield, IC2_ndpc_r, IC2, not_c_prime);});
    coroutines.emplace_back([&tio, &IC2_dpc_l, IC2, not_c_prime](yield_t &yield)
        { mpc_and(tio, yield, IC2_dpc_l, IC2, not_c_prime);});
    coroutines.emplace_back([&tio, &IC2_dpc_r, IC2, c_prime](yield_t &yield)
        { mpc_and(tio, yield, IC2_dpc_r, IC2, c_prime);});
    run_coroutines(tio, coroutines); 
    coroutines.clear();

    cs_bal_dpc^=IC2;
    p_bal_ndpc^=IC2;
  
    coroutines.emplace_back([&tio, &new_p_bal_l, IC2_ndpc_l, p_bal_ndpc](yield_t &yield)
        { mpc_select(tio, yield, new_p_bal_l, IC2_ndpc_l, new_p_bal_l, p_bal_ndpc);});
    coroutines.emplace_back([&tio, &new_p_bal_r, IC2_ndpc_r, p_bal_ndpc](yield_t &yield)
        { mpc_select(tio, yield, new_p_bal_r, IC2_ndpc_r, new_p_bal_r, p_bal_ndpc);});
    coroutines.emplace_back([&tio, &cs_bal_l, IC2_dpc_l, cs_bal_dpc](yield_t &yield)
        { mpc_select(tio, yield, cs_bal_l, IC2_dpc_l, cs_bal_l, cs_bal_dpc);});
    coroutines.emplace_back([&tio, &cs_bal_r, IC2_dpc_r, cs_bal_dpc](yield_t &yield)
        { mpc_select(tio, yield, cs_bal_r, IC2_dpc_r, cs_bal_r, cs_bal_dpc);});
    coroutines.emplace_back([&tio, &bal_upd, IC2, s0](yield_t &yield)
        {
        // In the IC2 case bal_upd = 0 (The rotation doesn't end up
        // decreasing height of this subtree.
        mpc_select(tio, yield, bal_upd, IC2, bal_upd, s0);});
    run_coroutines(tio, coroutines);
    coroutines.clear();

    // IC3:
    // To set balance in this case we need to know if gcs.dpc child exists
    // and similarly if gcs.ndpc child exitst.
    // if(gcs.ndpc child exists): cs.bal_ndpc = 1
    // if(gcs.dpc child exists): p.bal_dpc = 1
    RegBS gcs_dpc_exists, gcs_ndpc_exists;
    RegXS gcs_l = getAVLLeftPtr(gcs_node.pointers);
    RegXS gcs_r = getAVLRightPtr(gcs_node.pointers);
    RegBS gcs_bal_l = getLeftBal(gcs_node.pointers);
    RegBS gcs_bal_r = getRightBal(gcs_node.pointers);
    RegXS gcs_dpc, gcs_ndpc;

    run_coroutines(tio, [&tio, &gcs_dpc, c_prime, gcs_l, gcs_r] (yield_t &yield)
        { mpc_select(tio, yield, gcs_dpc, c_prime, gcs_l, gcs_r);},
        [&tio, &gcs_ndpc, not_c_prime, gcs_l, gcs_r] (yield_t &yield)
        { mpc_select(tio, yield, gcs_ndpc, not_c_prime, gcs_l, gcs_r);});

    CDPF cdpf = tio.cdpf(yield);
    run_coroutines(tio, [&tio, &gcs_dpc_exists, gcs_dpc, &cdpf](yield_t &yield) 
        { gcs_dpc_exists = cdpf.is_zero(tio, yield, gcs_dpc, tio.aes_ops());},
        [&tio, &gcs_ndpc_exists, gcs_ndpc, &cdpf](yield_t &yield)
        { gcs_ndpc_exists = cdpf.is_zero(tio, yield, gcs_ndpc, tio.aes_ops());});

    cs_bal_ndpc^=IC3;
    RegBS IC3_ndpc_l, IC3_ndpc_r, IC3_dpc_l, IC3_dpc_r;
    
    run_coroutines(tio, [&tio, &IC3_ndpc_l, IC3, c_prime](yield_t &yield)
        { mpc_and(tio, yield, IC3_ndpc_l, IC3, c_prime);},
        [&tio, &IC3_ndpc_r, IC3, not_c_prime](yield_t &yield)
        { mpc_and(tio, yield, IC3_ndpc_r, IC3, not_c_prime);},
        [&tio, &IC3_dpc_l, IC3, not_c_prime](yield_t &yield)
        { mpc_and(tio, yield, IC3_dpc_l, IC3, not_c_prime);},
        [&tio, &IC3_dpc_r, IC3, c_prime](yield_t &yield)
        { mpc_and(tio, yield, IC3_dpc_r, IC3, c_prime);});

    RegBS f0, f1, f2, f3;
    run_coroutines(tio, [&tio, &f0, IC3_dpc_l, gcs_dpc_exists] (yield_t &yield)
        { mpc_and(tio, yield, f0, IC3_dpc_l, gcs_dpc_exists);},
        [&tio, &f1, IC3_dpc_r, gcs_dpc_exists] (yield_t &yield)
        { mpc_and(tio, yield, f1, IC3_dpc_r, gcs_dpc_exists);},
        [&tio, &f2, IC3_ndpc_l, gcs_ndpc_exists] (yield_t &yield)
        { mpc_and(tio, yield, f2, IC3_ndpc_l, gcs_ndpc_exists);},
        [&tio, &f3, IC3_ndpc_r, gcs_ndpc_exists] (yield_t &yield) 
        { mpc_and(tio, yield, f3, IC3_ndpc_r, gcs_ndpc_exists);});

    
    coroutines.emplace_back([&tio, &new_p_bal_l, f0, IC3](yield_t &yield) {
        mpc_select(tio, yield, new_p_bal_l, f0, new_p_bal_l, IC3);});
    coroutines.emplace_back([&tio, &new_p_bal_r, f1, IC3](yield_t &yield) {
        mpc_select(tio, yield, new_p_bal_r, f1, new_p_bal_r, IC3);});
    coroutines.emplace_back([&tio, &cs_bal_l, f2, IC3](yield_t &yield) {
        mpc_select(tio, yield, cs_bal_l, f2, cs_bal_l, IC3);});
    coroutines.emplace_back([&tio, &cs_bal_r, f3, IC3](yield_t &yield) {
        mpc_select(tio, yield, cs_bal_r, f3, cs_bal_r, IC3);});
    // In IC3 gcs.bal = 0 0
    coroutines.emplace_back([&tio, &gcs_bal_l, IC3, s0](yield_t &yield) {
        mpc_select(tio, yield, gcs_bal_l, IC3, gcs_bal_l, s0);});
    coroutines.emplace_back([&tio, &gcs_bal_r, IC3, s0](yield_t &yield) {
        mpc_select(tio, yield, gcs_bal_r, IC3, gcs_bal_r, s0);});
    run_coroutines(tio, coroutines); 

    // Write back <cs_bal_dpc, cs_bal_ndpc> and <gcs_bal_l, gcs_bal_r>
    setLeftBal(gcs_node.pointers, gcs_bal_l);
    setRightBal(gcs_node.pointers, gcs_bal_r);
    setLeftBal(cs_node.pointers, cs_bal_l);
    setRightBal(cs_node.pointers, cs_bal_r);

    A[oidx_cs].NODE_POINTERS+= (cs_node.pointers - old_cs_ptr);
    A[oidx_gcs].NODE_POINTERS+= (gcs_node.pointers - old_gcs_ptr);

    // Write back updated pointers correctly accounting for rotations
    setLeftBal(nodeptrs, new_p_bal_l);
    setRightBal(nodeptrs, new_p_bal_r);
    #ifdef OPT_ON
        A[oidx].NODE_POINTERS +=(nodeptrs - oidx_oldptrs);
    #else
        A[ptr].NODE_POINTERS = nodeptrs;
    #endif
}

/* Update the return structure
   F_dh = Delete Here flag,
   F_sf = successor found (no more left children while trying to find successor)
   F_rs is a subflag for F_r (update children pointers with ret ptr)
   F_rs: Flag for updating the correct child pointer of this node
   This happens if F_r is set in ret_struct. F_r indicates if we need
   to update a child pointer at this level by skipping the current
   child in the direction of traversal. We do this in two cases:
     i) F_d & (!F_2) : If we delete here, and this node does not have
        2 children (;i.e., we are not in the finding successor case)
    ii) F_sf: Found the successor (no more left children while
        traversing to find successor)
   In cases i and ii we skip the next node, and make the current node
   point to the node after the next node.

   The third case for F_r:
   iii) We did rotation(s) at the lower level, changing the child in
        that position. So we update it to the correct node in that
        position now.
   Whether skip happens or just update happens is handled by how
   ret_struct.ret_ptr is set.
*/

void AVL::updateRetStruct(MPCTIO &tio, yield_t &yield, RegXS ptr, RegBS F_2, RegBS F_c2, 
        RegBS F_c4, RegBS lf, RegBS F_ri, RegBS &found, RegBS &bal_upd, 
        avl_del_return &ret_struct) {
    bool player0 = tio.player()==0;
    RegBS s0, s1;
    s1.set(tio.player()==1);
    RegBS F_dh, F_sf, F_rs;
    RegBS not_found = found;
    if(player0) {
        not_found^=1;
    }

    run_coroutines(tio, [&tio, &ret_struct, F_c2](yield_t &yield)
        { mpc_or(tio, yield, ret_struct.F_ss, ret_struct.F_ss, F_c2);},
        [&tio, &F_dh, lf, not_found] (yield_t &yield)
        { mpc_and(tio, yield, F_dh, lf, not_found);},
        [&tio, &ret_struct, F_dh, ptr] (yield_t &yield)
        { mpc_select(tio, yield, ret_struct.N_d, F_dh, ret_struct.N_d, ptr);});

    // F_sf = Successor found = F_c4 = Finding successor & no more left child
    F_sf = F_c4;
    if(player0)
        F_2^=1;
    // If we have to i) delete here, and it doesn't have two children
    // we have to update child pointer in parent with the returned pointer
    mpc_and(tio, yield, F_rs, F_dh, F_2);
    // ii) if we found successor here
    run_coroutines(tio, [&tio, &F_rs, F_sf](yield_t &yield)
        { mpc_or(tio, yield, F_rs, F_rs, F_sf);},
        [&tio, &ret_struct, F_sf, ptr] (yield_t &yield)
        { mpc_select(tio, yield, ret_struct.N_s, F_sf, ret_struct.N_s, ptr);});

    // F_rs and F_ri will never trigger together. So the line below
    // set ret_ptr to the correct pointer to handle either case
    // If neither F_rs nor F_ri, we set the ret_ptr to current ptr.
    RegBS F_nr;
    mpc_or(tio, yield, F_nr, F_rs, F_ri);
    // F_nr = F_rs || F_ri
    ret_struct.F_r = F_nr;

    if(player0) {
      F_nr^=1;
    }
    // F_nr = !(F_rs || F_ri)
    run_coroutines(tio, [&tio, &ret_struct, F_nr, ptr](yield_t &yield)
        { mpc_select(tio, yield, ret_struct.ret_ptr, F_nr, ret_struct.ret_ptr, ptr);},
        [&tio, &bal_upd, F_rs, s1](yield_t &yield)
        { // If F_rs, we skipped a node, so update bal_upd to 1
          mpc_select(tio, yield, bal_upd, F_rs, bal_upd, s1);});
}

std::tuple<bool, RegBS> AVL::del(MPCTIO &tio, yield_t &yield, RegXS ptr, RegAS del_key,
      Duoram<Node>::Flat &A, RegBS found, RegBS find_successor, int TTL,
      avl_del_return &ret_struct) {
    bool player0 = tio.player()==0;
    if(TTL==0) {
        //Reconstruct and return found
        bool success = reconstruct_RegBS(tio, yield, found);
        RegBS zero;
        return {success, zero};
    } else {
        Node node;
        RegXS oldptrs;
        #ifdef OPT_ON
            typename Duoram<Node>::template OblivIndex<RegXS,1> oidx(tio, yield, ptr, MAX_DEPTH);
            node = A[oidx];
            oldptrs = node.pointers;
        #else
            node = A[ptr];
        #endif

        // Compare key
        CDPF cdpf = tio.cdpf(yield);
        auto [lt, eq, gt] = cdpf.compare(tio, yield, del_key - node.key, tio.aes_ops());

        // c is the direction bit for next_ptr
        // (c=0: go left or c=1: go right)
        RegBS c = gt;
        // lf = local found. We found the key to delete in this level.
        RegBS lf = eq;

        // Select the next ptr
        RegXS left = getAVLLeftPtr(node.pointers);
        RegXS right = getAVLRightPtr(node.pointers);

        size_t &aes_ops = tio.aes_ops();
        RegBS l0, r0;
        // Check if left and right children are 0, and compute F_0, F_1, F_2
        run_coroutines(tio, [&tio, &l0, left, &aes_ops, &cdpf](yield_t &yield)
            { l0 = cdpf.is_zero(tio, yield, left, aes_ops);},
            [&tio, &r0, right, &aes_ops, &cdpf](yield_t &yield)
            { r0 = cdpf.is_zero(tio, yield, right, aes_ops);});

        RegBS F_0, F_1, F_2;
        RegBS F_c1, F_c2, F_c3, F_c4;
        RegXS next_ptr, cs_ptr;
        RegBS c_prime;

        // F_1 = l0 \xor r0
        F_1 = l0 ^ r0;

        // F_0 = l0 & r0
        // Case 1: lf & F_1
        run_coroutines(tio, [&tio, &F_0, l0, r0](yield_t &yield)
            { mpc_and(tio, yield, F_0, l0, r0);},
            [&tio, &F_c1, lf, F_1](yield_t &yield)
            { mpc_and(tio, yield, F_c1, lf, F_1);});

        // F_2 = !(F_0 + F_1) (Only 1 of F_0, F_1, and F_2 can be true)
        F_2 = F_0 ^ F_1;
        if(player0)
            F_2^=1;
        // s1: shares of 1 bit, s0: shares of 0 bit
        RegBS s1, s0;
        s1.set(tio.player()==1);

        // We set next ptr based on c, but we need to handle three
        // edge cases where we do not pick next_ptr by just the comparison result
        // Case 1: found the node here (lf), and node has only one child.
        // Then we iterate down the only child.
        // Set c_prime for Case 1
        run_coroutines(tio, [&tio, &c_prime, F_c1, c, l0](yield_t &yield)
            { mpc_select(tio, yield, c_prime, F_c1, c, l0);},
            [&tio, &F_c2, lf, F_2](yield_t &yield)
            { mpc_and(tio, yield, F_c2, lf, F_2);});

        // Case 2: found the node here (lf) and node has both children (F_2)
        // In find successor case, so we find inorder successor for node to be deleted
        // (inorder successor = go right and then find leftmost child.)

        // Case 3: finding successor (find_successor) and node has both children (F_2)
        // Go left.

        run_coroutines(tio, [&tio, &c_prime, F_c2, s1](yield_t &yield)
            { mpc_select(tio, yield, c_prime, F_c2, c_prime, s1);},
            [&tio, &F_c3, find_successor, F_2](yield_t &yield)
            { mpc_and(tio, yield, F_c3, find_successor, F_2);});

        // Case 4: finding successor (find_successor) and node has no more left children (l0)
        // This is the successor node then.
        // Go right (since no more left)
        run_coroutines(tio, [&tio, &c_prime, F_c3, s0](yield_t &yield) 
            { mpc_select(tio, yield, c_prime, F_c3, c_prime, s0);},
            [&tio, &F_c4, find_successor, l0](yield_t &yield)
            { mpc_and(tio, yield, F_c4, find_successor, l0);});

        RegBS found_prime, find_successor_prime;
        mpc_select(tio, yield, c_prime, F_c4, c_prime, l0);
        // Set next_ptr
        mpc_select(tio, yield, next_ptr, c_prime, left, right, AVL_PTR_SIZE);
        // cs_ptr: child's sibling pointer

        run_coroutines(tio, [&tio, &cs_ptr, c_prime, right, left](yield_t &yield)
            { mpc_select(tio, yield, cs_ptr, c_prime, right, left, AVL_PTR_SIZE);},
            [&tio, &found_prime, found, lf](yield_t &yield)
            { mpc_or(tio, yield, found_prime, found, lf);},
            // If in Case 2, set find_successor. We are now finding successor
            [&tio, &find_successor_prime, find_successor, F_c2](yield_t &yield)
            { mpc_or(tio, yield, find_successor_prime, find_successor, F_c2);});

        // If in Case 4. Successor found here already. Toggle find_successor off
        find_successor_prime=find_successor_prime^F_c4;

        TTL-=1;
        auto [key_found, bal_upd] = del(tio, yield, next_ptr, del_key, A, found_prime, find_successor_prime, TTL, ret_struct);

        // If we didn't find the key, we can end here.
        if(!key_found) {
          return {0, s0};
        }

        updateChildPointers(tio, yield, left, right, c_prime, ret_struct);
        setAVLLeftPtr(node.pointers, left);
        setAVLRightPtr(node.pointers, right);
        // Delay storing pointers back until balance updates are done as well.
        // Since we resolved the F_r flag returned with updateChildPointers(),
        // we set it back to 0.
        ret_struct.F_r = s0;

        RegBS p_bal_l, p_bal_r;
        p_bal_l = getLeftBal(node.pointers);
        p_bal_r = getRightBal(node.pointers);
        auto [new_p_bal_l, new_p_bal_r, new_bal_upd, imb] =
            updateBalanceDel(tio, yield, p_bal_l, p_bal_r, bal_upd, c_prime);
        
        // F_ri: subflag for F_r. F_ri = returned flag set to 1 from imbalance fix.
        RegBS F_ri;
        fixImbalance(tio, yield, A, oidx, oldptrs, ptr, node.pointers, new_p_bal_l, new_p_bal_r, bal_upd, 
              c_prime, cs_ptr, imb, F_ri, ret_struct);

        updateRetStruct(tio, yield, ptr, F_2, F_c2, F_c4, lf, F_ri, found, bal_upd, ret_struct); 

        return {key_found, bal_upd};
    }
}


bool AVL::del(MPCTIO &tio, yield_t &yield, RegAS del_key) {
    if(num_items==0)
        return 0;

    auto A = oram.flat(tio, yield);
    if(num_items==1) {
        //Delete root if root's key = del_key
        Node zero;
        typename Duoram<Node>::template OblivIndex<RegXS,1> oidx(tio, yield, root, MAX_DEPTH);
        Node node = A[oidx];
        // Compare key
        CDPF cdpf = tio.cdpf(yield);
        auto [lt, eq, gt] = cdpf.compare(tio, yield, del_key - node.key, tio.aes_ops());
        bool success = reconstruct_RegBS(tio, yield, eq);
        if(success) {
            empty_locations.emplace_back(root);
            A[oidx] = zero;
            num_items--;
            return 1;
        } else {
            return 0;
        } 
    } else {
        int TTL = AVL_TTL(num_items);
        // Flags for already found (found) item to delete and find successor (find_successor)
        // if this deletion requires a successor swap
        RegBS found, find_successor;
        avl_del_return ret_struct;
        auto [success, bal_upd] = del(tio, yield, root, del_key, A, found, find_successor, TTL, ret_struct);
        printf ("Success =  %d\n", success);
        if(!success){
            return 0;
        }
        else{
            num_items--;
            /*
            printf("In delete's swap portion\n");
            Node rec_del_node = A.reconstruct(A[ret_struct.N_d]);
            Node rec_suc_node = A.reconstruct(A[ret_struct.N_s]);
            printf("del_node key = %ld, suc_node key = %ld\n",
                rec_del_node.key.ashare, rec_suc_node.key.ashare);
            printf("flag_s = %d\n", ret_struct.F_ss.bshare);
            */
            Node del_node, suc_node;
            typename Duoram<Node>::template OblivIndex<RegXS,2> oidx_nd(tio, yield, ret_struct.N_d, MAX_DEPTH);
            typename Duoram<Node>::template OblivIndex<RegXS,2> oidx_ns(tio, yield, ret_struct.N_s, MAX_DEPTH);
            #ifdef OPT_ON
                del_node = A[oidx_nd];
                suc_node = A[oidx_ns];
            #else
                del_node = A[ret_struct.N_d];
                suc_node = A[ret_struct.N_s];
            #endif
            RegAS zero_as; RegXS zero_xs;
            // Update root if needed
            mpc_select(tio, yield, root, ret_struct.F_r, root, ret_struct.ret_ptr);

            /*
            bool rec_F_ss = reconstruct_RegBS(tio, yield, ret_struct.F_ss);
            size_t rec_del_key = reconstruct_RegAS(tio, yield, del_node.key);
            size_t rec_suc_key = reconstruct_RegAS(tio, yield, suc_node.key);
            printf("rec_F_ss = %d, del_node.key = %lu, suc_nod.key = %lu\n",
                rec_F_ss, rec_del_key, rec_suc_key);
            */            
            RegXS old_del_value;
            RegAS old_del_key;
            #ifdef OPT_ON
                old_del_value = del_node.value;
                old_del_key = del_node.key;
            #endif
            RegXS empty_loc;

            run_coroutines(tio, [&tio, &del_node, ret_struct, suc_node](yield_t &yield)
                { mpc_select(tio, yield, del_node.key, ret_struct.F_ss, del_node.key, suc_node.key);},
                [&tio, &del_node, ret_struct, suc_node] (yield_t &yield)
                { mpc_select(tio, yield, del_node.value, ret_struct.F_ss, del_node.value, suc_node.value);},
                [&tio, &empty_loc, ret_struct](yield_t &yield)
                { mpc_select(tio, yield, empty_loc, ret_struct.F_ss, ret_struct.N_d, ret_struct.N_s);});

            #ifdef OPT_ON
                A[oidx_nd].NODE_KEY+=(del_node.key - old_del_key);
                A[oidx_nd].NODE_VALUE+=(del_node.value - old_del_value);
                A[oidx_ns].NODE_KEY+=(-suc_node.key);
                A[oidx_ns].NODE_VALUE+=(suc_node.value);
            #else
                A[ret_struct.N_d].NODE_KEY = del_node.key;
                A[ret_struct.N_d].NODE_VALUE = del_node.value;
                A[ret_struct.N_s].NODE_KEY = zero_as;
                A[ret_struct.N_s].NODE_VALUE = zero_xs;
            #endif

            //Add deleted (empty) location into the empty_locations vector for reuse in next insert()
            empty_locations.emplace_back(empty_loc);
        }

      return 1;
    }
}

void AVL::initialize(MPCTIO &tio, yield_t &yield, size_t depth) {
    size_t init_size = (size_t(1)<<depth) - 1;
    auto A = oram.flat(tio, yield);

    for(size_t i=1; i<=depth; i++) {
        size_t start = size_t(1)<<(i-1);
        size_t gap = size_t(1)<<i;
        size_t current = start;
        for(size_t j=1; j<=(size_t(1)<<(depth-i)); j++) {
            //printf("current = %ld ", current);
            Node node;
            node.key.set(current * tio.player());
            if(i!=1) {
                //Set left and right child pointers and balance bits
                size_t ptr_gap = start/2;
                RegXS lptr, rptr;
                lptr.set(tio.player() * (current-(ptr_gap)));
                rptr.set(tio.player() * (current+(ptr_gap)));
                setAVLLeftPtr(node.pointers, lptr);
                setAVLRightPtr(node.pointers, rptr);
            }
            A[current] = node;
            current+=gap;
        }
    }

    // Set num_items to init_size after they have been initialized;
    num_items = init_size;
    // Set root correctly
    root.set(tio.player() * size_t(1)<<(depth-1));
} 

// Now we use the AVL class in various ways.  This function is called by
// online.cpp.
void avl(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
{
    size_t depth=4, n_inserts=0, n_deletes=0;
    if (*args) {
        depth = atoi(args[0]);
        n_inserts = atoi(args[1]);
        n_deletes = atoi(args[2]);
    }

    /* The ORAM will be initialized with 2^depth-1 items, but the 0 slot is reserved.
       So we initialize (initial inserts) with 2^depth-2 items.
       The ORAM size is set to 2^depth-1 + n_insert.
    */
    size_t init_size = (size_t(1)<<(depth));
    size_t oram_size = init_size + 1 + n_inserts; // +1 because init_size does not account for slot at 0.

    MPCTIO tio(mpcio, 0, opts.num_threads);
    run_coroutines(tio, [&tio, &mpcio, depth, oram_size, init_size, n_inserts, n_deletes] (yield_t &yield) {

        std::cout << "\n===== SETUP =====\n";
        AVL tree(tio.player(), oram_size);
        tree.initialize(tio, yield, depth);
        //tree.pretty_print(tio, yield);
        tio.sync_lamport();

        Node node;
        mpcio.dump_stats(std::cout);
        std::cout << "\n===== INSERTS =====\n";
        mpcio.reset_stats();
        tio.reset_lamport();

        for(size_t i = 1; i<=n_inserts; i++) {
            newnode(node);
            node.key.set((i+init_size) * tio.player());
            tree.insert(tio, yield, node);
        }

        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        std::cout << "\n===== DELETES =====\n";
        mpcio.reset_stats();
        tio.reset_lamport();
        for(size_t i = 1; i<=n_deletes; i++) {
            RegAS del_key;
            del_key.set((i+init_size) * tio.player());
            tree.del(tio, yield, del_key);
        }
    });
}


void avl_tests(MPCIO &mpcio,
    const PRACOptions &opts, char **args)
{
    // Not taking arguments for tests
    nbits_t depth=4;
    size_t items = (size_t(1)<<depth)-1;

    MPCTIO tio(mpcio, 0, opts.num_threads);
    run_coroutines(tio, [&tio, depth, items] (yield_t &yield) {
        size_t size = size_t(1)<<depth;
        bool player0 = tio.player()==0;
        AVL tree(tio.player(), size);

        // (T1) : Test 1 : L rotation (root modified)
        /*
            Operation:

              5                   7
               \                /   \
                7       --->   5     9
                 \
                  9

            T1 checks:
            - root is 7
            - 5,7,9 in correct positions
            - 5 and 9 have no children and 0 balances
        */
        {
            bool success = 1;
            int insert_array[] = {5, 7, 9};
            size_t insert_array_size = 2;
            Node node;

            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }
            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();

            Node root_node, left_node, right_node;
            size_t left_index, right_index;
            root_node = R[root];
            if((root_node.key).share()!=7) {
                success = false;
            }
            left_index = (getAVLLeftPtr(root_node.pointers)).share();
            right_index = (getAVLRightPtr(root_node.pointers)).share();
            left_node = R[left_index];
            right_node = R[right_index];

            if(left_node.key.share()!=5 || right_node.key.share()!=9) {
                success = false;
            }

            //To check that left and right have no children and 0 balances
            size_t sum = left_node.pointers.share() + right_node.pointers.share();
            if(sum!=0) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T1 : SUCCESS\n");
                } else {
                    print_red("T1 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }

        // (T2) : Test 2 : L rotation (root unmodified)
        /*
            Operation:
                 5                        5
               /   \                    /   \
              3     7                  3     9
                     \     --->             /  \
                      9          7         7    12
                       \
                        12


            T2 checks:
            - root is 5
            - 3, 7, 9, 12 in expected positions
            - Nodes 3, 7, 12 have 0 balance and no children
            - 5's bal = 0 1

        */
        {
            bool success = 1;
            int insert_array[] = {5, 3, 7, 9, 12};
            size_t insert_array_size = 4;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n3, n7, n9, n12;
            size_t n3_index, n7_index, n9_index, n12_index;
            root_node = R[root];
            if((root_node.key).share()!=5) {
                success = false;
            }
            n3_index = (getAVLLeftPtr(root_node.pointers)).share();
            n9_index = (getAVLRightPtr(root_node.pointers)).share();
            n3 = R[n3_index];
            n9 = R[n9_index];
            n7_index = getAVLLeftPtr(n9.pointers).share();
            n12_index = getAVLRightPtr(n9.pointers).share();
            n7 = R[n7_index];
            n12 = R[n12_index];

            // Node value checks
            if(n3.key.share()!=3 || n9.key.share()!=9) {
                success = false;
            }
            if(n7.key.share()!=7 || n12.key.share()!=12) {
                success = false;
            }

            // Node children and balance checks
            size_t zero = 0;
            zero+=(n3.pointers.share());
            zero+=(n7.pointers.share());
            zero+=(n12.pointers.share());
            zero+=(getLeftBal(root_node.pointers).share());
            zero+=(getLeftBal(n9.pointers).share());
            zero+=(getRightBal(n9.pointers).share());
            if(zero!=0) {
                success = false;
            }
            int one = (getRightBal(root_node.pointers).share());
            if(one!=1) {
                success = false;
            }

            if(player0) {
                if(success) {
                    print_green("T2 : SUCCESS\n");
                } else {
                    print_red("T2 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }


        // (T3) : Test 3 : R rotation (root modified)
        /*
            Operation:

                   9                 7
                  /                /   \
                 7         --->   5     9
                /
               5

            T3 checks:
            - root is 7
            - 5,7,9 in correct positions
            - 5 and 9 have no children

        */
        {
            bool success = 1;
            int insert_array[] = {9, 7, 5};
            size_t insert_array_size = 2;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();

            Node root_node, left_node, right_node;
            size_t left_index, right_index;
            root_node = R[root];
            if((root_node.key).share()!=7) {
                success = false;
            }
            left_index = (getAVLLeftPtr(root_node.pointers)).share();
            right_index = (getAVLRightPtr(root_node.pointers)).share();
            left_node = R[left_index];
            right_node = R[right_index];

            if(left_node.key.share()!=5 || right_node.key.share()!=9) {
                success = false;
            }

            //To check that left and right have no children and 0 balances
            size_t sum = left_node.pointers.share() + right_node.pointers.share();
            if(sum!=0) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T3 : SUCCESS\n");
                } else{
                    print_red("T3 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }


        // (T4) : Test 4 : R rotation (root unmodified)
        /*
            Operation:
                     9                        9
                   /   \                    /   \
                  7     12                 5     12
                 /              --->      / \
                5                    7   3   7
               /
              3


            T4 checks:
            - root is 9
            - 3,5,7,12 are in correct positions
            - Nodes 3,7,12 have 0 balance
            - Nodes 3,7,12 have no children
            - 9's bal = 1 0

        */
        {
            bool success = 1;
            int insert_array[] = {9, 12, 7, 5, 3};
            size_t insert_array_size = 4;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n3, n7, n5, n12;
            size_t n3_index, n7_index, n5_index, n12_index;
            root_node = R[root];
            if((root_node.key).share()!=9) {
                success = false;
            }
            n5_index = (getAVLLeftPtr(root_node.pointers)).share();
            n12_index = (getAVLRightPtr(root_node.pointers)).share();
            n5 = R[n5_index];
            n12 = R[n12_index];
            n3_index = getAVLLeftPtr(n5.pointers).share();
            n7_index = getAVLRightPtr(n5.pointers).share();
            n7 = R[n7_index];
            n3 = R[n3_index];

            // Node value checks
            if(n12.key.share()!=12 || n5.key.share()!=5) {
                success = false;
            }
            if(n3.key.share()!=3 || n7.key.share()!=7) {
                success = false;
            }

            // Node balance checks
            size_t zero = 0;
            zero+=(n3.pointers.share());
            zero+=(n7.pointers.share());
            zero+=(n12.pointers.share());
            zero+=(getRightBal(root_node.pointers).share());
            zero+=(getLeftBal(n5.pointers).share());
            zero+=(getRightBal(n5.pointers).share());
            if(zero!=0) {
                success = false;
            }
            int one = (getLeftBal(root_node.pointers).share());
            if(one!=1) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T4 : SUCCESS\n");
                } else {
                    print_red("T4 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }


        // (T5) : Test 5 : LR rotation (root modified)
        /*
            Operation:
                    9              9           7
                   /              /          /   \
                  5      -->     7     -->  5     9
                   \            /
                    7          5


            T5 checks:
            - root is 7
            - 9,5,7 are in correct positions
            - Nodes 5,7,9 have 0 balance
            - Nodes 5,9 have no children
        */

        {
            bool success = 1;
            int insert_array[] = {9, 5, 7};
            size_t insert_array_size = 2;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n9, n5;
            size_t n9_index, n5_index;
            root_node = R[root];
            if((root_node.key).share()!=7) {
                success = false;
            }
            n5_index = (getAVLLeftPtr(root_node.pointers)).share();
            n9_index = (getAVLRightPtr(root_node.pointers)).share();
            n5 = R[n5_index];
            n9 = R[n9_index];

            // Node value checks
            if(n9.key.share()!=9 || n5.key.share()!=5) {
                success = false;
            }

            // Node balance checks
            size_t zero = 0;
            zero+=(n5.pointers.share());
            zero+=(n9.pointers.share());
            zero+=(getRightBal(root_node.pointers).share());
            zero+=(getLeftBal(n5.pointers).share());
            zero+=(getRightBal(n5.pointers).share());
            zero+=(getLeftBal(n5.pointers).share());
            zero+=(getRightBal(n9.pointers).share());
            zero+=(getLeftBal(n9.pointers).share());
            if(zero!=0) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T5 : SUCCESS\n");
                } else {
                    print_red("T5 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }


        // (T6) : Test 6 : LR rotation (root unmodified)
        /*
            Operation:

                 9                      9                   9
               /   \                  /   \               /   \
              7     12               7     12            5     12
             /             --->     /           --->    / \
            3                      5                   3   7
             \                    /
              5                  3


            T6 checks:
            - root is 9
            - 3,5,7,12 are in correct positions
            - Nodes 3,7,12 have 0 balance
            - Nodes 3,7,12 have no children
            - 9's bal = 1 0

        */
        {
            bool success = 1;
            int insert_array[] = {9, 12, 7, 3, 5};
            size_t insert_array_size = 4;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n3, n7, n5, n12;
            size_t n3_index, n7_index, n5_index, n12_index;
            root_node = R[root];
            if((root_node.key).share()!=9) {
                success = false;
            }
            n5_index = (getAVLLeftPtr(root_node.pointers)).share();
            n12_index = (getAVLRightPtr(root_node.pointers)).share();
            n5 = R[n5_index];
            n12 = R[n12_index];
            n3_index = getAVLLeftPtr(n5.pointers).share();
            n7_index = getAVLRightPtr(n5.pointers).share();
            n7 = R[n7_index];
            n3 = R[n3_index];

            // Node value checks
            if(n5.key.share()!=5 || n12.key.share()!=12) {
                success = false;
            }
            if(n3.key.share()!=3 || n7.key.share()!=7) {
                success = false;
            }

            // Node balance checks
            size_t zero = 0;
            zero+=(n3.pointers.share());
            zero+=(n7.pointers.share());
            zero+=(n12.pointers.share());
            zero+=(getRightBal(root_node.pointers).share());
            zero+=(getLeftBal(n5.pointers).share());
            zero+=(getRightBal(n5.pointers).share());
            if(zero!=0) {
                success = false;
            }
            int one = (getLeftBal(root_node.pointers).share());
            if(one!=1) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T6 : SUCCESS\n");
                } else {
                    print_red("T6 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }


        // (T7) : Test 7 : RL rotation (root modified)
        /*
            Operation:
                5              5                7
                 \              \             /   \
                  9      -->     7      -->  5     9
                 /                \
                7                  9


            T7 checks:
            - root is 7
            - 9,5,7 are in correct positions
            - Nodes 5,7,9 have 0 balance
            - Nodes 5,9 have no children
        */

        {
            bool success = 1;
            int insert_array[] = {5, 9, 7};
            size_t insert_array_size = 2;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n9, n5;
            size_t n9_index, n5_index;
            root_node = R[root];
            if((root_node.key).share()!=7) {
                success = false;
            }
            n5_index = (getAVLLeftPtr(root_node.pointers)).share();
            n9_index = (getAVLRightPtr(root_node.pointers)).share();
            n5 = R[n5_index];
            n9 = R[n9_index];

            // Node value checks
            if(n9.key.share()!=9 || n5.key.share()!=5) {
                success = false;
            }

            // Node balance checks
            size_t zero = 0;
            zero+=(n5.pointers.share());
            zero+=(n9.pointers.share());
            zero+=(getRightBal(root_node.pointers).share());
            zero+=(getLeftBal(n5.pointers).share());
            zero+=(getRightBal(n5.pointers).share());
            zero+=(getLeftBal(n5.pointers).share());
            zero+=(getRightBal(n9.pointers).share());
            zero+=(getLeftBal(n9.pointers).share());
            if(zero!=0) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T7 : SUCCESS\n");
                } else {
                    print_red("T7 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }

        // (T8) : Test 8 : RL rotation (root unmodified)
        /*
            Operation:

                 5                     5                   5
               /   \                 /   \               /   \
              3     12              3     12            3     9
                   /      --->           /     --->          /  \
                  7                     9                   7    12
                   \                   /
                    9                 7


            T8 checks:
            - root is 5
            - 3,9,7,12 are in correct positions
            - Nodes 3,7,12 have 0 balance
            - Nodes 3,7,12 have no children
            - 5's bal = 0 1

        */
        {
            bool success = 1;
            int insert_array[] = {5, 3, 12, 7, 9};
            size_t insert_array_size = 4;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n3, n7, n9, n12;
            size_t n3_index, n7_index, n9_index, n12_index;
            root_node = R[root];
            if((root_node.key).share()!=5) {
                success = false;
            }
            n3_index = (getAVLLeftPtr(root_node.pointers)).share();
            n9_index = (getAVLRightPtr(root_node.pointers)).share();
            n3 = R[n3_index];
            n9 = R[n9_index];
            n7_index = getAVLLeftPtr(n9.pointers).share();
            n12_index = getAVLRightPtr(n9.pointers).share();
            n7 = R[n7_index];
            n12 = R[n12_index];

            // Node value checks
            if(n3.key.share()!=3 || n9.key.share()!=9) {
                success = false;
            }
            if(n7.key.share()!=7 || n12.key.share()!=12) {
                success = false;
            }

            // Node balance checks
            size_t zero = 0;
            zero+=(n3.pointers.share());
            zero+=(n7.pointers.share());
            zero+=(n12.pointers.share());
            zero+=(getLeftBal(root_node.pointers).share());
            zero+=(getLeftBal(n9.pointers).share());
            zero+=(getRightBal(n9.pointers).share());
            if(zero!=0) {
                success = false;
            }
            int one = (getRightBal(root_node.pointers).share());
            if(one!=1) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T8 : SUCCESS\n");
                } else {
                    print_red("T8 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }

        // Deletion Tests:

        // (T9) : Test 9 : L rotation (root modified)
        /*
            Operation:

              5                    7
             / \       Del 3     /   \
            3   7     ------>   5     9
                 \
                  9

            T9 checks:
            - root is 7
            - 5,7,9 in correct positions
            - 5 and 9 have no children and 0 balances
            - 7 has 0 balances
        */
        {
            bool success = 1;
            int insert_array[] = {5, 3, 7, 9};
            size_t insert_array_size = 3;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            RegAS del_key;
            del_key.set(3 * tio.player());
            tree.del(tio, yield, del_key);
            tree.check_avl(tio, yield);

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();

            Node root_node, left_node, right_node;
            size_t left_index, right_index;
            root_node = R[root];
            if((root_node.key).share()!=7) {
                success = false;
            }
            left_index = (getAVLLeftPtr(root_node.pointers)).share();
            right_index = (getAVLRightPtr(root_node.pointers)).share();
            left_node = R[left_index];
            right_node = R[right_index];

            if(left_node.key.share()!=5 || right_node.key.share()!=9) {
                success = false;
            }

            //To check that left and right have no children and 0 balances
            size_t sum = left_node.pointers.share() + right_node.pointers.share();
            if(sum!=0) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T9 : SUCCESS\n");
                } else {
                    print_red("T9 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }

        // (T10) : Test 10 : L rotation (root unmodified)
        /*
            Operation:
                 5                        5
               /   \                    /   \
              3     7      Del 6       3     9
             /     / \    ------>     /     /  \
            1     6   9              1     7    12
                       \
                        12


            T10 checks:
            - root is 5
            - 3, 7, 9, 12 in expected positions
            - Nodes 3, 7, 12 have 0 balance and no children
            - 5's bal = 0 1

        */
        {
            bool success = 1;
            int insert_array[] = {5, 3, 7, 9, 6, 1, 12};
            size_t insert_array_size = 6;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            RegAS del_key;
            del_key.set(6 * tio.player());
            tree.del(tio, yield, del_key);
            tree.check_avl(tio, yield);

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n1, n3, n7, n9, n12;
            size_t n1_index, n3_index, n7_index, n9_index, n12_index;
            root_node = R[root];
            if((root_node.key).share()!=5) {
                success = false;
            }
            n3_index = (getAVLLeftPtr(root_node.pointers)).share();
            n9_index = (getAVLRightPtr(root_node.pointers)).share();
            n3 = R[n3_index];
            n9 = R[n9_index];
            n7_index = getAVLLeftPtr(n9.pointers).share();
            n12_index = getAVLRightPtr(n9.pointers).share();
            n7 = R[n7_index];
            n12 = R[n12_index];
            n1_index = getAVLLeftPtr(n3.pointers).share();
            n1 = R[n1_index];

            // Node value checks
            if(n3.key.share()!=3 || n9.key.share()!=9) {
                success = false;
            }
            if(n7.key.share()!=7 || n12.key.share()!=12 || n1.key.share()!=1) {
                success = false;
            }

            // Node children and balance checks
            size_t zero = 0;
            zero+=(n1.pointers.share());
            zero+=(n7.pointers.share());
            zero+=(n12.pointers.share());
            zero+=(getLeftBal(root_node.pointers).share());
            zero+=(getRightBal(root_node.pointers).share());
            zero+=(getLeftBal(n9.pointers).share());
            zero+=(getRightBal(n9.pointers).share());
            zero+=(getRightBal(n3.pointers).share());
            if(zero!=0) {
                success = false;
            }
            int one = (getLeftBal(n3.pointers).share());
            if(one!=1) {
                success = false;
            }

            if(player0) {
                if(success) {
                    print_green("T10 : SUCCESS\n");
                } else {
                    print_red("T10 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }

        // (T11) : Test 11 : R rotation (root modified)
        /*
            Operation:

                   9                   7
                  / \     Del 12     /   \
                 7   12  ------->   5     9
                /
               5

            T11 checks:
            - root is 7
            - 5,7,9 in correct positions and balances to 0
            - 5 and 9 have no children

        */
        {
            bool success = 1;
            int insert_array[] = {9, 7, 12, 5};
            size_t insert_array_size = 3;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            RegAS del_key;
            del_key.set(12 * tio.player());
            tree.del(tio, yield, del_key);
            tree.check_avl(tio, yield);

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();

            Node root_node, left_node, right_node;
            size_t left_index, right_index;
            root_node = R[root];
            if((root_node.key).share()!=7) {
                success = false;
            }
            left_index = (getAVLLeftPtr(root_node.pointers)).share();
            right_index = (getAVLRightPtr(root_node.pointers)).share();
            left_node = R[left_index];
            right_node = R[right_index];

            if(left_node.key.share()!=5 || right_node.key.share()!=9) {
                success = false;
            }

            //To check that left and right have no children and 0 balances
            size_t zero = left_node.pointers.share() + right_node.pointers.share();
            zero+=(getLeftBal(left_node.pointers).share());
            zero+=(getRightBal(left_node.pointers).share());
            zero+=(getLeftBal(right_node.pointers).share());
            zero+=(getRightBal(right_node.pointers).share());
            if(zero!=0) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T11 : SUCCESS\n");
                } else{
                    print_red("T11 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }


        // (T12) : Test 12 : R rotation (root unmodified)
        /*
            Operation:
                     9                        9
                   /   \                    /   \
                  7     12      Del 8      5     12
                 / \     \     ------>    / \     \
                5   8     15             3   7     15
               /
              3


            T4 checks:
            - root is 9
            - 3,5,7,12,15 are in correct positions
            - Nodes 3,7,15 have 0 balance
            - Nodes 3,7,15 have no children
            - 9,5 bal = 0 0
            - 12 bal = 0 1
        */
        {
            bool success = 1;
            int insert_array[] = {9, 12, 7, 5, 8, 15, 3};
            size_t insert_array_size = 6;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            RegAS del_key;
            del_key.set(8 * tio.player());
            tree.del(tio, yield, del_key);
            tree.check_avl(tio, yield);

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n3, n7, n5, n12, n15;
            size_t n3_index, n7_index, n5_index, n12_index, n15_index;
            root_node = R[root];
            if((root_node.key).share()!=9) {
                success = false;
            }
            n5_index = (getAVLLeftPtr(root_node.pointers)).share();
            n12_index = (getAVLRightPtr(root_node.pointers)).share();
            n5 = R[n5_index];
            n12 = R[n12_index];
            n3_index = getAVLLeftPtr(n5.pointers).share();
            n7_index = getAVLRightPtr(n5.pointers).share();
            n7 = R[n7_index];
            n3 = R[n3_index];
            n15_index = getAVLRightPtr(n12.pointers).share();
            n15 = R[n15_index];

            // Node value checks
            if(n12.key.share()!=12 || n5.key.share()!=5) {
                success = false;
            }
            if(n3.key.share()!=3 || n7.key.share()!=7 || n15.key.share()!=15) {
                success = false;
            }

            // Node balance checks
            size_t zero = 0;
            zero+=(n3.pointers.share());
            zero+=(n7.pointers.share());
            zero+=(n15.pointers.share());
            zero+=(getRightBal(root_node.pointers).share());
            zero+=(getLeftBal(root_node.pointers).share());
            zero+=(getLeftBal(n5.pointers).share());
            zero+=(getRightBal(n5.pointers).share());
            if(zero!=0) {
                success = false;
            }
            int one = (getRightBal(n12.pointers).share());
            if(one!=1) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T12 : SUCCESS\n");
                } else {
                    print_red("T12 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }


        // (T13) : Test 13 : LR rotation (root modified)
        /*
            Operation:
                    9                   9           7
                   / \     Del 12      /          /   \
                  5   12  ------->    7     -->  5     9
                   \                 /
                    7               5


            T5 checks:
            - root is 7
            - 9,5,7 are in correct positions
            - Nodes 5,7,9 have 0 balance
            - Nodes 5,9 have no children
        */

        {
            bool success = 1;
            int insert_array[] = {9, 5, 12, 7};
            size_t insert_array_size = 3;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            RegAS del_key;
            del_key.set(12 * tio.player());
            tree.del(tio, yield, del_key);
            tree.check_avl(tio, yield);

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n9, n5;
            size_t n9_index, n5_index;
            root_node = R[root];
            if((root_node.key).share()!=7) {
                success = false;
            }
            n5_index = (getAVLLeftPtr(root_node.pointers)).share();
            n9_index = (getAVLRightPtr(root_node.pointers)).share();
            n5 = R[n5_index];
            n9 = R[n9_index];

            // Node value checks
            if(n9.key.share()!=9 || n5.key.share()!=5) {
                success = false;
            }

            // Node balance checks
            size_t zero = 0;
            zero+=(n5.pointers.share());
            zero+=(n9.pointers.share());
            zero+=(getRightBal(root_node.pointers).share());
            zero+=(getLeftBal(n5.pointers).share());
            zero+=(getRightBal(n5.pointers).share());
            zero+=(getLeftBal(n5.pointers).share());
            zero+=(getRightBal(n9.pointers).share());
            zero+=(getLeftBal(n9.pointers).share());
            if(zero!=0) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T13 : SUCCESS\n");
                } else {
                    print_red("T13 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }


        // (T14) : Test 14 : LR rotation (root unmodified)
        /*
            Operation:

                 9                      9                     9
               /   \                  /   \                 /   \
              7     12     Del 8     7     12              5     12
             / \          ------>   /             --->    / \
            3   8                  5                     3   7
             \                    /
              5                  3


            T6 checks:
            - root is 9
            - 3,5,7,12 are in correct positions
            - Nodes 3,7,12 have 0 balance
            - Nodes 3,7,12 have no children
            - 9's bal = 1 0

        */
        {
            bool success = 1;
            int insert_array[] = {9, 12, 7, 3, 5};
            size_t insert_array_size = 4;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            RegAS del_key;
            del_key.set(8 * tio.player());
            tree.del(tio, yield, del_key);
            tree.check_avl(tio, yield);

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n3, n7, n5, n12;
            size_t n3_index, n7_index, n5_index, n12_index;
            root_node = R[root];
            if((root_node.key).share()!=9) {
                success = false;
            }
            n5_index = (getAVLLeftPtr(root_node.pointers)).share();
            n12_index = (getAVLRightPtr(root_node.pointers)).share();
            n5 = R[n5_index];
            n12 = R[n12_index];
            n3_index = getAVLLeftPtr(n5.pointers).share();
            n7_index = getAVLRightPtr(n5.pointers).share();
            n7 = R[n7_index];
            n3 = R[n3_index];

            // Node value checks
            if(n5.key.share()!=5 || n12.key.share()!=12) {
                success = false;
            }
            if(n3.key.share()!=3 || n7.key.share()!=7) {
                success = false;
            }

            // Node balance checks
            size_t zero = 0;
            zero+=(n3.pointers.share());
            zero+=(n7.pointers.share());
            zero+=(n12.pointers.share());
            zero+=(getRightBal(root_node.pointers).share());
            zero+=(getLeftBal(n5.pointers).share());
            zero+=(getRightBal(n5.pointers).share());
            if(zero!=0) {
                success = false;
            }
            int one = (getLeftBal(root_node.pointers).share());
            if(one!=1) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T14 : SUCCESS\n");
                } else {
                    print_red("T14 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }

        // (T15) : Test 15 : RL rotation (root modified)
        /*
            Operation:
                5              5                7
               / \     Del 3    \             /   \
              3   9   ------->   7      -->  5     9
                 /                \
                7                  9


            T15 checks:
            - root is 7
            - 9,5,7 are in correct positions
            - Nodes 5,7,9 have 0 balance
            - Nodes 5,9 have no children
        */

        {
            bool success = 1;
            int insert_array[] = {5, 9, 3, 7};
            size_t insert_array_size = 3;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            RegAS del_key;
            del_key.set(3 * tio.player());
            tree.del(tio, yield, del_key);
            tree.check_avl(tio, yield);

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n9, n5;
            size_t n9_index, n5_index;
            root_node = R[root];
            if((root_node.key).share()!=7) {
                success = false;
            }
            n5_index = (getAVLLeftPtr(root_node.pointers)).share();
            n9_index = (getAVLRightPtr(root_node.pointers)).share();
            n5 = R[n5_index];
            n9 = R[n9_index];

            // Node value checks
            if(n9.key.share()!=9 || n5.key.share()!=5) {
                success = false;
            }

            // Node balance checks
            size_t zero = 0;
            zero+=(n5.pointers.share());
            zero+=(n9.pointers.share());
            zero+=(getRightBal(root_node.pointers).share());
            zero+=(getLeftBal(n5.pointers).share());
            zero+=(getRightBal(n5.pointers).share());
            zero+=(getLeftBal(n5.pointers).share());
            zero+=(getRightBal(n9.pointers).share());
            zero+=(getLeftBal(n9.pointers).share());
            if(zero!=0) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T15 : SUCCESS\n");
                } else {
                    print_red("T15 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }

        // (T16) : Test 16 : RL rotation (root unmodified)
        /*
            Operation:

                 5                     5                   5
               /   \                 /   \               /   \
              3     12    Del 1     3     12            3     9
             /      /    ------>          /     --->         / \
            1      7                     9                  7   12
                    \                   /
                     9                 7


            T8 checks:
            - root is 5
            - 3,9,7,12 are in correct positions
            - Nodes 3,7,12 have 0 balance
            - Nodes 3,7,12 have no children
            - 5's bal = 0 1

        */
        {
            bool success = 1;
            int insert_array[] = {5, 3, 12, 7, 1, 9};
            size_t insert_array_size = 5;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            RegAS del_key;
            del_key.set(1 * tio.player());
            tree.del(tio, yield, del_key);
            tree.check_avl(tio, yield);

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n3, n7, n9, n12;
            size_t n3_index, n7_index, n9_index, n12_index;
            root_node = R[root];
            if((root_node.key).share()!=5) {
                success = false;
            }
            n3_index = (getAVLLeftPtr(root_node.pointers)).share();
            n9_index = (getAVLRightPtr(root_node.pointers)).share();
            n3 = R[n3_index];
            n9 = R[n9_index];
            n7_index = getAVLLeftPtr(n9.pointers).share();
            n12_index = getAVLRightPtr(n9.pointers).share();
            n7 = R[n7_index];
            n12 = R[n12_index];

            // Node value checks
            if(n3.key.share()!=3 || n9.key.share()!=9) {
                success = false;
            }
            if(n7.key.share()!=7 || n12.key.share()!=12) {
                success = false;
            }

            // Node balance checks
            size_t zero = 0;
            zero+=(n3.pointers.share());
            zero+=(n7.pointers.share());
            zero+=(n12.pointers.share());
            zero+=(getLeftBal(root_node.pointers).share());
            zero+=(getLeftBal(n9.pointers).share());
            zero+=(getRightBal(n9.pointers).share());
            if(zero!=0) {
                success = false;
            }
            int one = (getRightBal(root_node.pointers).share());
            if(one!=1) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T16 : SUCCESS\n");
                } else {
                    print_red("T16 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }


        // (T17) : Test 17 : Double imbalance (root modified)
        /*
            Operation:

                      9                                     9
                   /     \                               /     \
                  5        12        Del 10             5        15
                /   \     /  \      -------->        /    \      /  \
               3     7   10   15                    3      7    12   20
              / \   / \        \                   / \    / \
             2   4 6   8        20                2   4  6   8
            /                                    /
           1                                    1


                                5
                             /     \
                           3          9
            ----->        / \       /   \
                         2   4     7      15
                        /         / \    /   \
                       1         6   8  10    20




            T17 checks:
            - root is 5
            - all other nodes are in correct positions
            - balances and children are correct
        */
        {
            bool success = 1;
            int insert_array[] = {9, 5, 12, 7, 3, 10, 15, 2, 4, 6, 8, 20, 1};
            size_t insert_array_size = 12;
            Node node;
            for(size_t i = 0; i<=insert_array_size; i++) {
              newnode(node);
              node.key.set(insert_array[i] * tio.player());
              tree.insert(tio, yield, node);
              tree.check_avl(tio, yield);
            }

            RegAS del_key;
            del_key.set(10 * tio.player());
            tree.del(tio, yield, del_key);
            tree.check_avl(tio, yield);

            Duoram<Node>* oram = tree.get_oram();
            RegXS root_xs = tree.get_root();
            size_t root = reconstruct_RegXS(tio, yield, root_xs);
            auto A = oram->flat(tio, yield);
            auto R = A.reconstruct();
            Node root_node, n3, n7, n9;
            Node n1, n2, n4, n6, n8, n12, n15, n20;
            size_t n3_index, n7_index, n9_index;
            size_t n1_index, n2_index, n4_index, n6_index;
            size_t n8_index, n12_index, n15_index, n20_index;
            root_node = R[root];
            if((root_node.key).share()!=5) {
                success = false;
            }
            n3_index = (getAVLLeftPtr(root_node.pointers)).share();
            n9_index = (getAVLRightPtr(root_node.pointers)).share();
            n3 = R[n3_index];
            n9 = R[n9_index];

            n2_index = getAVLLeftPtr(n3.pointers).share();
            n4_index = getAVLRightPtr(n3.pointers).share();
            n7_index = getAVLLeftPtr(n9.pointers).share();
            n15_index = getAVLRightPtr(n9.pointers).share();
            n2 = R[n2_index];
            n4 = R[n4_index];
            n7 = R[n7_index];
            n15 = R[n15_index];

            n1_index = getAVLLeftPtr(n2.pointers).share();
            n6_index = getAVLLeftPtr(n7.pointers).share();
            n8_index = getAVLRightPtr(n7.pointers).share();
            n12_index = getAVLLeftPtr(n15.pointers).share();
            n20_index = getAVLRightPtr(n15.pointers).share();
            n1 = R[n1_index];
            n6 = R[n6_index];
            n8 = R[n8_index];
            n12 = R[n12_index];
            n20 = R[n20_index];

            // Node value checks
            if(n3.key.share()!=3 || n9.key.share()!=9) {
                success = false;
            }
            if(n2.key.share()!=2 || n4.key.share()!=4) {
                success = false;
            }
            if(n7.key.share()!=7 || n15.key.share()!=15) {
                success = false;
            }
            if(n1.key.share()!=1 || n6.key.share()!=6 || n8.key.share()!=8) {
                success = false;
            }
            if(n12.key.share()!=12 || n20.key.share()!=20) {
                success = false;
            }

            // Node balance checks
            size_t zero = 0;
            zero+=(n1.pointers.share());
            zero+=(n4.pointers.share());
            zero+=(n6.pointers.share());
            zero+=(n8.pointers.share());
            zero+=(n12.pointers.share());
            zero+=(n20.pointers.share());
            zero+=(getLeftBal(n7.pointers).share());
            zero+=(getRightBal(n7.pointers).share());
            zero+=(getLeftBal(n9.pointers).share());
            zero+=(getRightBal(n9.pointers).share());
            zero+=(getLeftBal(n15.pointers).share());
            zero+=(getRightBal(n15.pointers).share());
            zero+=(getRightBal(n3.pointers).share());
            zero+=(getLeftBal(root_node.pointers).share());
            zero+=(getRightBal(root_node.pointers).share());
            if(zero!=0) {
                success = false;
            }
            int one = (getLeftBal(n3.pointers).share());
            if(one!=1) {
                success = false;
            }
            if(player0) {
                if(success) {
                    print_green("T17 : SUCCESS\n");
                } else {
                    print_red("T17 : FAIL\n");
                }
            }
            A.init();
            tree.init();
        }

    });
}
