"""
Script for route optimization for multiple searchers
"""
import pandas as pd
import numpy as np
from itertools import product
import gurobipy as gp
from gurobipy import GRB
import os
import platform
import time
import json

##################### Helper function ###########################
def is_nearby_cell(c, c_prime):
    if c == c_prime:
        return True
    elif c[0] == c_prime[0]:
        if c[1] + 1 == c_prime[1] or c[1] - 1 == c_prime[1]:
            return True
    elif c[1] == c_prime[1]:
        if c[0] + 1 == c_prime[0] or c[0] - 1 == c_prime[0]:
            return True
    else:
        return False
    
def is_corner_cell(c, grid_size):
    if c[0] == 1 and c[1] == 1:
        return True
    elif c[0] == 1 and c[1] == grid_size:
        return True
    elif c[0] == grid_size and c[1] == 1:
        return True
    elif c[0] == grid_size and c[1] == grid_size:
        return True
    else:
        return False
    
def is_side_cell(c, grid_size):
    if c[0] == 1 or c[0] == grid_size:
        if c[1] not in [1, grid_size]:
            return True
    elif c[1] == 1 or c[1] == grid_size:
        if c[0] not in [1, grid_size]:
            return True
    else:
        return False
    
def return_nearby_cell(c, grid_size):
    """ return nearby cells given a cell and grid_size """
    stay = c
    up = (c[0] - 1, c[1]) if (c[0] - 1) >= 1 else None
    down = (c[0] + 1, c[1]) if (c[0] + 1) <= grid_size else None
    left = (c[0], c[1] - 1) if (c[1] - 1) >= 1 else None
    right = (c[0], c[1] + 1) if (c[1] + 1) <= grid_size else None
    nearby =  [stay, up, down, left, right]
    return [direction for direction in nearby if direction != None]

def is_searcher_occ(C, T, grid_size):
    searcher_occ = {}
    for t in T:                
        for c in C:
            # when t = 1, searcher can only occupy (1, 1), (1, 2) and (2, 1) because searchers are at (1, 1) at time 0
            if t == 1: 
                nearby_cell = return_nearby_cell((1, 1), grid_size)
                for each in nearby_cell:
                    searcher_occ[each, 1] = 1
            # return searcher_occ
            else:
                already_occ = list(searcher_occ.keys())
                last_iter_occ = [c_t for c_t in already_occ if c_t[1] == t - 1] # check where the searcher could be at (t - 1)
                for c_t in last_iter_occ:
                    # if c_t[1] == t - 1:
                    nearby_cell = return_nearby_cell(c_t[0], grid_size)
                    for each in nearby_cell:
                        searcher_occ[each, t] = 1
    return searcher_occ

##################### End of helper function ####################

# ending_time_grid = list(range(7, 16))
# ending_time_grid = [7, 8 , 9 , 10]
ending_time_grid = [7, 8, 9, 10, 11, 12, 13, 14, 15]


""" Import data
"""
if platform.system() == 'Windows':
    data_folder = 'E:\\Research\\Route_Optimization_for_Multiple_Searchers\\Python\\'
else:
    data_folder = os.path.dirname(os.path.realpath(__file__))
# zeta_raw = pd.read_csv(data_folder + '/Zeta.csv', header = None, index_col = 0)
# q_raw = pd.read_csv(data_folder + '/q.csv', header = 0, index_col = 0)

time_log = {}
for ending_time in ending_time_grid:
    # ending_time = ending_
    # ending_time = 9
    print('===========================')
    print('ending time is', ending_time)
    print('===========================')

    grid_size = 9
    print('===== grid size is =====', )
    ending_time = ending_time
    # num_scenario = 1000
    
    """
    Creating set
    """    
    C = [(i, j) for i in range(1, grid_size + 1) for j in range(1, grid_size + 1)]
    T = list(range(1, ending_time + 1))
    T0 = [0] + T
    # Omega = list(range(1, num_scenario + 1))
    J = 15
    I = list(range(0, J * ending_time + 1))
    # print('i is', I)
    
    xx = {}
    for c in C:
        for t in T0:
            if c == (1, 1) and t == 0:
                xx[c, t] = J
            else:
                xx[c, t] = 0
        
    """
    Creating parameters
    """
    p = {}
    for c in C:
        if c == ((grid_size + 1) / 2, (grid_size + 1) / 2):
            p[c] = 1
        else:
            p[c] = 0
    
    alpha_sub = list(product(C, T))
    alpha_value = [-3 * np.log(0.4) / J] * grid_size * grid_size * ending_time
    alpha = dict(zip(alpha_sub, alpha_value))
    alpha = -3 * np.log(0.4) / J
    
    
    sub_gamma = list(product(C, C, T))
    gamma = {}
    stay_prob = 0.6
    for sub in sub_gamma:
        c, c_prime, t = sub
        if c == c_prime:
            gamma[sub] = stay_prob
        elif is_corner_cell(c, grid_size = grid_size):
            if is_nearby_cell(c, c_prime):
                gamma[sub] = (1 - stay_prob) / 2
            else:
                gamma[sub] = 0
        elif is_side_cell(c, grid_size = grid_size):
            if is_nearby_cell(c, c_prime):
                gamma[sub] = (1 - stay_prob) / 3
            else:
                gamma[sub] = 0
        else:
            if is_nearby_cell(c, c_prime):
                gamma[sub] = (1 - stay_prob) / 4
            else:
                gamma[sub] = 0
    
    sub_q = list(product(C, T)) 
    sub_q = sorted(sub_q, key = lambda x: x[1])
    q = {} # create param values for q
    for sub in sub_q:
        c_two_dim, t = sub
        # print(sub)
        if t == 1:
            q[c_two_dim, t] = p[c_two_dim]
        else:
            # c_one_dim = (c_two_dim[0] - 1) * grid_size + c_two_dim[1]
            q[c_two_dim, t] = sum([q[c_prime, t - 1] * gamma[c_prime, c_two_dim, t - 1] for c_prime in C if is_nearby_cell(c_two_dim, c_prime)])
    
    searcher_occ = is_searcher_occ(C, T, grid_size)
    occ_cell = list(searcher_occ.keys())
    
    W_param = {}
    for c in C:
        for t in T:
            if q[c, t] == 0 or (c, t) not in occ_cell: # no target or searchers cannot occupy at t
                W_param[c, t] = 0
            else:
                W_param[c, t] = q[c, t]
    
    W_param = dict(sorted(W_param.items(), key = lambda x: x[0][1]))
    
    searcher_reachable = {}
    for t in T:
        for c in C:
            if (c,t) not in occ_cell:
                searcher_reachable[c, t] = 0
            else:
                searcher_reachable[c, t] = 1
        
    T_no_detect = []
    for t in T:
        if sum(q[c,t] * searcher_reachable[c,t] for c in C) == 0:
            T_no_detect.append(t)
    
    part_T_no_detect = {}
    for t in set(T) - set(T_no_detect):
        for c in C:
            if q[c, t] * searcher_reachable[c, t] == 0:
                if t not in part_T_no_detect.keys():
                    part_T_no_detect[t] = [(c, t)]
                else:
                    part_T_no_detect[t].append((c, t))
# =============================================================================
#     num_groups = 11
#     # max_W_param = max(W_param.values())
#     
#     # W_param_excl_0 = [val for val in W_param.values() if val != 0]
#     # W_param_excl_0 = [val for val in W_param.values()]
#     # min_W_param = min(W_param_excl_0)
#     cat_group = {}
#     W_param_cut = [0, 0.02, 0.04, 0.06, 0.09, 0.12, 0.1240, 0.1455, 0.1757, 0.22, 0.2881]
#     for c in C:
#         for t in T:
#             if W_param[c, t] == 0:
#                 cat_group[c, t] = 1 
#             else:
#                 for ind in range(0, len(W_param_cut) - 1):
#                     if W_param_cut[ind] < W_param[c, t] <= W_param_cut[ind + 1]:
#                         cat_group[c, t] = ind + 2 # assign a group for each c,t
#     assert len(cat_group.keys()) == len(W_param.keys()), "Not all c,t pairs are grouped"
# =============================================================================
    # Z_0 = m.addVars(T_no_detect, lb = J, ub = J, name = 'Z_0', vtype = GRB.CONTINUOUS)
    # constr_Z_0 = m.a
    part_T_no_detect_pairs = sum(part_T_no_detect.values(), [])

    # num_groups = 5
    cat_group = {}
    # W_param_t = {}
    
    
    # =============================================================================
    #         W_param_t = {c_t: prob for c_t, prob in W_param.items() if c_t[1] == t}
    #         max_W_param = max(W_param_t.values())
    #         # W_param_t_prob = [val for val in W_param_t.values()]
    #         min_W_param = min(W_param_t.values())
    #         W_param_cut = list(np.linspace(min_W_param, max_W_param, num = num_groups))
    #         if min_W_param != 0:
    #             W_param_cut.append(0)
    #             W_param_cut = sorted(W_param_cut)
    #         for c in C:
    #             # for t in T:                
    #             if W_param_t[c, t] == 0:
    #                 # print('c,t is', c, t)
    #                 cat_group[c, t] = (1, t) 
    #             else:
    #                 # print('c,t is', c, t)
    #                 for ind in range(0, len(W_param_cut) - 1):
    #                     if W_param_cut[ind] < W_param_t[c, t] <= W_param_cut[ind + 1]:
    #                         # print('assigned')
    #                         cat_group[c, t] = (ind + 2, t) # assign a group for each c,t
    # =============================================================================
    
    
    for t in range(1, ending_time + 1):
        if t not in T_no_detect:
            for c in C:
                if (c, t) not in part_T_no_detect_pairs:
        # for ind, c in enumerate(C):
                    cat_group[c, t] = (c, t)
    print('c, t pairs in the individual c,t group', len((list(cat_group.keys()))))
    print('c, t pairs in the full no value t group', len(T_no_detect) * grid_size * grid_size)
    print('c, t pairs in the part no value g group', len(part_T_no_detect_pairs))
    ct_pairs_sum = len(list(cat_group.keys())) + len(T_no_detect) * grid_size * grid_size + len(part_T_no_detect_pairs)
    print('sum is', ct_pairs_sum)
    print('the total number of c t paris are', grid_size * grid_size * ending_time)
    assert ct_pairs_sum == grid_size * grid_size * ending_time, "Not all c,t pairs are properly bucketized"
    # assert len(cat_group.keys()) == len(W_param.keys()), "Not all c,t pairs are grouped"    
    
# =============================================================================
#     test_obj = {}
#     
#     for c_t, k_t in cat_group.items():
#         if c_t[1] not in test_obj.keys():
#             test_obj[c_t[1]] = 1
#         else:
#             test_obj[c_t[1]] += 1
# =============================================================================
            
    
    
    
    cat_group = dict(sorted(cat_group.items(), key = lambda x: x[0][1]))
    
# =============================================================================
#     group_prob = {}
#     for c in C:
#         for t in T:
#             group = cat_group[c, t]
#             if group in group_prob.keys():
#                 group_prob[group].append(W_param[c, t])
#             else:
#                 group_prob[group] = [W_param[c,t]]
#     
#     group_prob_diff = {}
#     
#     for group in group_prob.keys():
#         group_prob_diff[group] = [min(group_prob[group]), max(group_prob[group])]
#         group_prob_diff[group].append(group_prob_diff[group][1] - group_prob_diff[group][0])
#         
#     print("===== group probability difference =====", group_prob_diff)
# =============================================================================
# =============================================================================
#     np.random.seed(2022)
#     q = np.random.uniform(low = 0, high = 1, size = num_scenario)
#     q = q / sum(q) # normalize to a probablity distribution summing up to 1
#     q = dict(zip(Omega, q))
#     alpha = -3 * np.log(0.4) / J
# =============================================================================
    
    # =============================================================================
    # q = pd.read_csv(data_folder + 'q.csv')
    # q = dict(zip(q['index'], q.q))
    # =============================================================================
    
    """
    Defining decision variables
    """
    
    
    # create sets
    sub_X = list(product(C, C, T0))
    sub_Z = list(product(C, T))
    
    
    
    ############ step 0 ###################    
    delta = 1e-4

    
    Xi_lb = 0
    Xi_ub = 1
    counter = 1
    
    cuts = []
    
    model_name = 'sp1_algo2'
    m = gp.Model(model_name)
    m.setParam(GRB.Param.TimeLimit, 15 * 60)
    m.setParam(GRB.Param.Threads, 1)
    m.setParam(GRB.Param.LogFile, model_name)
    # m.NumStart = 0
    
    Z_0 = m.addVars(T_no_detect, lb = J, ub = J, name = 'Z_0', vtype = GRB.CONTINUOUS)
    part_Z_0 = m.addVars(list(part_T_no_detect.keys()), name = 'part_Z_0', vtype = GRB.CONTINUOUS)
    
    # add variables
    # X = m.addVars(sub_X, lb = 0, name = 'X', vtype = GRB.INTEGER)
    X = m.addVars(sub_X, lb = 0, name = 'X', vtype = GRB.CONTINUOUS)

    # Z = m.addVars(sub_Z, lb = 0, ub = J, vtype = GRB.INTEGER, name = 'Z')
    
    group_cnt = {}
    for group in cat_group.values():
        if group not in group_cnt.keys():
            group_cnt[group] = 1
        else:
            group_cnt[group] += 1
    
    # record 
    group_by_cell = {}
    for c_t, k_t in cat_group.items():
        if k_t not in group_by_cell.keys():
            group_by_cell[k_t] = [c_t]
        else:
            group_by_cell[k_t].append(c_t)
            
    group_by_cell = dict(sorted(group_by_cell.items(), key = lambda x: (x[0][1], x[0][0])))
# =============================================================================
#     print('group(k,t) and its component cell(c,t)')
#     for k_t, c_t in group_by_cell.items():
#         print(k_t, ':', c_t, '\n')
# =============================================================================
    
    
    print(' ===== number of cells per group =====', group_cnt)
    
    ZZZ = {} # dict to save ZZZ variable
    ZZZ_param = {}
    for group in group_cnt.keys():
        ZZZ[group] = m.addVar(lb = 0, ub = J, vtype = GRB.INTEGER)
        # ZZZ[group] = m.addVar(lb = 0, ub = 6, vtype = GRB.INTEGER)
        ZZZ[group].Start = 0
        ZZZ_param[group] = 0
    
    
    Z_param = {}
    for sub in sub_Z:
        Z_param[sub] = 0
    
    Xi = m.addVar(vtype = GRB.CONTINUOUS, name = 'Xi')
    # m.update()
    

    """
    Defining objective function
    """
    m.setObjective(Xi, GRB.MINIMIZE)    
    
    """
    Defining constraints
    """
    constr_14 = m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == sum(X[c, c_prime, t] for c_prime in C if is_nearby_cell(c, c_prime))  for c in C for t in T), name = '14') #2d
    constr_15 = m.addConstrs((sum(X[c, c_prime, 0] for c_prime in C if is_nearby_cell(c, c_prime)) == xx[c, 0] for c in C), name = '15') #2d
    
    group_by_t = {}
    for k_t in group_cnt.keys():
        if k_t[1] not in group_by_t.keys():
            group_by_t[k_t[1]] = [k_t[0]] # add group k for time t
        else:
            group_by_t[k_t[1]].append(k_t[0]) # add group k for time t
    
    # 2.22
    constr_Z_t_J = m.addConstrs((part_Z_0[t] + sum(ZZZ[k, t] for k in group_by_t[t] ) == J for t in group_by_t.keys()), name = 'Z_kt bound')
    # cell_by_group 
    
    # 2.23
    # constr_Z_ct_0 = m.addConstrs((ZZZ[group] == 0 for group in group_cnt.keys() if searcher_reachable[group] == 0), name = 'z_ct_0')
    
    
    constr_16 = {}
    for group in group_cnt.keys():
        # group = 1
        # constr_16[group] = (m.addConstr((sum(X[c_prime, c, t - 1] for c in C for t in T for c_prime in C if cat_group[c, t] == group and is_nearby_cell(c, c_prime)) == ZZZ[group]), name = '16_' + str(group))) #2d
        # if searcher_reachable[group] == 0:
        #    print(group)
        constr_16[group] = (m.addConstr((sum(X[c_prime, group[0], group[1] - 1] for c_prime in C if is_nearby_cell(group[0], c_prime)) == ZZZ[group]), name = '16_' + str(group))) #2d
    
    for t in T_no_detect:
        constr_link_x_Z_0 = m.addConstr((sum(X[c_prime, c, t - 1] for c in C for c_prime in C if is_nearby_cell(c, c_prime)) == Z_0[t]), name = '16_link_Z_0_' + str(t)) 
        
    for t in part_T_no_detect.keys():
        constr_link_x_part_Z_0 = m.addConstr((sum(X[c_prime, c_t[0], t - 1] for c_t in part_T_no_detect[t] for c_prime in C if searcher_reachable[c_t] > 0 and is_nearby_cell(c, c_prime)) == part_Z_0[t]), name = '16_link_part_Z_0_' + str(t)) 
# =============================================================================
#     constr_Z_ct_fd_0 = m.addConstrs((ZZZ[group] == 0 for group in group_cnt.keys() if searcher_reachable[group] == 1 and Remove_finite_diff[group] == 0), name = 'z_ct_fd_0')
#             
#             
#     constr_x_z = {}
#     for group in group_cnt.keys():
#         # group = 1
#         # constr_16[group] = (m.addConstr((sum(X[c_prime, c, t - 1] for c in C for t in T for c_prime in C if cat_group[c, t] == group and is_nearby_cell(c, c_prime)) == ZZZ[group]), name = '16_' + str(group))) #2d
#         if searcher_reachable[group] == 1 and Remove_finite_diff[group] == 0:
#             constr_x_z[group] = (m.addConstr((sum(X[c_prime, group[0], group[1] - 1] for c_prime in C if is_nearby_cell(group[0], c_prime)) == ZZZ[group]), name = 'x_z' + str(group))) #2d
# =============================================================================

    # m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == Z[c, t] for c in C for t in T), name = '16') #2d
    lhs = {}
    lhs_val = {}
    zzz_val = {}
    x_val = {}
    z_recov_val = {}
    finite_diff_val = {}
    
    start_time = time.time()
    # while abs(Xi_ub - Xi_lb) > delta * Xi_lb: #and counter <= 100:
    # while counter <= 100 and abs(Xi_ub - Xi_lb) > delta * Xi_lb and time.time() - start_time <= 900:
    while abs(Xi_ub - Xi_lb) > delta * Xi_lb and time.time() - start_time <= 900:
        print('============= iteration', counter, '===============')
# =============================================================================
#         Z_param = Z_param_new.copy()
#         ################ step 1 ################
#         print('============= iteration', counter, '===============')
#         
#         # Z_param = Z_param_new.copy()
#         for group in group_cnt.keys():
#             ZZZ[group].lb = 0
#             ZZZ[group].ub = 0
#         
#         z_ct_optimal = {
#         (1,1):3,
#         (1,2):3,
#         (1,3):3,
#         (1,4):3,
#         (5,5):3,
#         (5,6):3,
#         (5,7):3,
#         (4,8):2,
#         (5,8):1,
#         (4,9):3,
#         # (4,9):1
#         }
#         
# 
#         for group in z_ct_optimal.keys():
#             ZZZ[group].lb = z_ct_optimal[group]
#             ZZZ[group].ub = z_ct_optimal[group]
#         
#         for each in sub_X:
#             X[each].ub = 0
#             X[each].lb = 0
#         
#         x_assign = {((1, 1), (2, 1), 0):3.0,
#                     ((2, 1), (3, 1), 1):3.0,
#                     ((3, 1), (3, 2), 2):3.0,
#                     ((3, 2), (4, 2), 3):3.0,
#                     ((4, 2), (4, 3), 4):3.0,
#                     ((4, 3), (4, 4), 5):3.0,
#                     ((4, 4), (4, 5), 6):3.0,
#                     ((4, 5), (3, 5), 9):2.0,
#                     ((4, 5), (4, 5), 7):2.0,
#                     ((4, 5), (4, 5), 8):2.0,
#                     ((4, 5), (5, 5), 7):1.0,
#                     ((5, 5), (5, 6), 8):1.0,
#                     ((5, 6), (4, 6), 9):1.0,
#                    }    
#         
#         for each in x_assign.keys():
#             X[each].ub = x_assign[each]
#             X[each].lb = x_assign[each]
# =============================================================================
        
        # ZZZ[(3,9)].ub = 2 
        
        # Xi.ub = 0.45
        
        r = {}
        for t in T:
            for c in C:
                if t == 1:
                    r[c, t] = p[c]
                else:
                    r[c, t] = sum([r[c_prime, t - 1] * np.exp(-alpha * Z_param[c_prime, t - 1]) * gamma[c_prime, c, t - 1] for c_prime in C])
        
# =============================================================================
#         r_stand = {}
#         for t in T:
#             for c in C:
#                 if t == 1:
#                     r_stand[c, t] = p[c]
#                 else:
#                     r_stand[c, t] = sum([r_stand[c_prime, t - 1] * np.exp(-alpha * Z_param_new[c_prime, t - 1]) * gamma[c_prime, c, t - 1] for c_prime in C])
#         
# =============================================================================
# =============================================================================
#         for t in T:
#             for c in C:
#                 if r_stand[c, t] != r[c, t]:
#                     print('r[c,t]', c, t, 'not equal', r[c,t], r_stand[c,t])
#         
# =============================================================================
        s = {}
        for t in T[::-1]:
            for c in C:
                if t == ending_time:
                    s[c, t] = 1
                else:
                    s[c, t] = sum([s[c_prime, t + 1] * np.exp(-alpha * Z_param[c_prime, t + 1]) * gamma[c, c_prime, t] for c_prime in C])
        
# =============================================================================
#         s_stand = {}
#         for t in T[::-1]:
#             for c in C:
#                 if t == ending_time:
#                     s_stand[c, t] = 1
#                 else:
#                     s_stand[c, t] = sum([s_stand[c_prime, t + 1] * np.exp(-alpha * Z_param_new[c_prime, t + 1]) * gamma[c, c_prime, t] for c_prime in C])
#                     
# =============================================================================
# =============================================================================
#         for t in T:
#             for c in C:
#                 if s_stand[c, t] != s[c, t]:
#                     print('s[c,t]', c, t, 'not equal', s[c,t], s_stand[c,t])
# =============================================================================
# =============================================================================
#         for key, val in r.items():
#             #if val != 0:
#             print(key, val, '\n')
# =============================================================================
        
        # f_Z = sum([r[c, 1] * np.exp(-alpha * Z_param[c, 1]) * s[c, 1] for c in C])
        # f_Z_2 = sum([r[c, ending_time] * np.exp(-alpha * Z_param[c, ending_time]) * s[c, ending_time] for c in C])
        f_Z = sum([r[c, 5] * np.exp(-alpha * Z_param[c, 5]) * s[c, 5] for c in C])

        print('f(Z) equals to', f_Z)
        # print('f(Z) equals to', f_Z_2)
        # print('f(Z) equals to', f_Z_5)
        
        if f_Z < Xi_ub:
            Xi_ub = f_Z   
# =============================================================================
#         if abs(Xi_ub - Xi_lb) <= delta * Xi_lb:
#             break
# =============================================================================
        
        print('Before optimization')
        print('Xi upper', Xi_ub)
        print('Xi lower', Xi_lb)
        

        ################ step 2 ##############
        # def solve_p():

        if Xi_lb == 0:
            g = np.inf
        else:
            g = abs(Xi_ub - Xi_lb) / Xi_lb
        print('counter is', counter)
        print('g is', g)
        
        if counter == 1:
            mip_gap = 0
        else:
            mip_gap = min([0.03, g / 3])
        # mip_gap = 0.1 / 2 ** (counter - 1)
        mip_gap = 1e-4
        print('==== MIP Gap ====', mip_gap)
        m.setParam("MIPGap", mip_gap)
        
        # turn off warm start
# =============================================================================
#         for ind in sub_X:
#             X[ind].Start = 0
#         for ind in sub_Z:
#             Z[ind].Start = 0
#         Xi = 1
# =============================================================================
    
    

        # m.addConstrs((sum(X[c_prime, sub[0], sub[1] - 1] for c_prime in C if is_nearby_cell(sub[0], c_prime)) == Z_New[sub] for sub in sub_WW), name = '16') #2d


        # f(Z^j) + sum (f(Z^j + delta_ct) - f(Z^j))(Z_ct - Z^j_ct) <= Xi
        # cuts.append(f_Z + sum([r[c, t] * (np.exp(-alpha * (Z_param[c, t] + 1)) - np.exp(-alpha * Z_param[c, t])) * s[c, t] * (Z[c, t] - Z_param[c, t]) for c in C for t in T]) <= Xi)
        
        # m.addConstr(f_Z + sum([r[c, t] * (np.exp(-alpha * (Z_param[c, t] + 1)) - np.exp(-alpha * Z_param[c, t])) * s[c, t] * (Z[c, t] - Z_param[c, t]) for c in C for t in T]) <= Xi, name = 'cut_' + str(counter))
        # finite_diff_coef = {}
        Remove_finite_diff = {}
# =============================================================================
#         Remove_r = {}
#         Remove_s = {}
#         Remove_middle_term_1 = {}
#         Remove_middle_term_2 = {}
# =============================================================================
        adj_finite_diff_coef = {}

        for t in T:
            for c in C:
                # Remove_r[c, t] = r[c, t]
                # Remove_s[c, t] = s[c, t]
                # Remove_middle_term_1[c, t] = np.exp(-alpha * (Z_param[c, t] + 1))
                # Remove_middle_term_2[c, t] = np.exp(-alpha * (Z_param[c, t]))
                if t in T_no_detect or (c,t) in part_T_no_detect_pairs:
                    Remove_finite_diff[c, t] = 0
                    continue
                else:
                    Remove_finite_diff[c, t] = r[c, t] * (np.exp(-alpha * (Z_param[c, t] + 1)) - np.exp(-alpha * Z_param[c, t])) * s[c, t]
                group = cat_group[c, t]
                if group in adj_finite_diff_coef.keys():
                    # finite_diff_coef[group].append(r[c, t] * (np.exp(-alpha * (Z_param[c, t] + 1)) - np.exp(-alpha * Z_param[c, t])) * s[c, t])
                    adj_finite_diff_coef[group].append(searcher_reachable[c, t] * r[c, t] * (np.exp(-alpha * (Z_param[c, t] + 1)) - np.exp(-alpha * Z_param[c, t])) * s[c, t])
                else:
                    # finite_diff_coef[group] = [r[c, t] * (np.exp(-alpha * (Z_param[c, t] + 1)) - np.exp(-alpha * Z_param[c, t])) * s[c, t]]
                    adj_finite_diff_coef[group] = [searcher_reachable[c, t] * r[c, t] * (np.exp(-alpha * (Z_param[c, t] + 1)) - np.exp(-alpha * Z_param[c, t])) * s[c, t]]
        
        

        
        
        
# =============================================================================
#         set_q_zero = []
#         for t in T:
#             for c in C:
#                 if q[c,t] == 0:
#                     set_q_zero.append((c,t))
# =============================================================================
        
# =============================================================================
#         for c_t in set_q_zero:
#             print(c_t, Remove_finite_diff[c_t])
# =============================================================================
        
        
# =============================================================================
#         q_sorted = dict(sorted(q.items(), key = lambda x: x[1]))
#         Remove_finite_diff = dict(sorted(Remove_finite_diff.items(), key = lambda x: x[1], reverse = True))
#         
# =============================================================================
# =============================================================================
#         for t in T:
#             for c in C:
#                 if Remove_finite_diff[c,t] == 0 and q[c,t] * W_param[c,t] > 0:
#                     print(c,t)
# =============================================================================
        
        
        # 0 + sum([min(adj_finite_diff_coef[group]) * (ZZZ[group] - ZZZ_param[group]) for group in group_cnt.keys()]) 
        # lhs[counter] = 0.732658674027209 + sum([min(adj_finite_diff_coef[group]) * (ZZZ[group] - ZZZ_param[group]) for group in group_cnt.keys()]) 
        lhs[counter] = f_Z + sum([min(adj_finite_diff_coef[group]) * (ZZZ[group] - ZZZ_param[group]) for group in group_cnt.keys()]) 
        # finite_diff_by_k = {group: min(adj_finite_diff_coef[group]) for group in group_cnt.keys()}
# =============================================================================
#         max_diff_by_k = {group: max(adj_finite_diff_coef[group]) for group in group_cnt.keys()}
# =============================================================================

# =============================================================================
#         print('========= finite difference by group ===========')
#         for key,val in finite_diff_by_k.items():
#             print(key, val)
#         print(finite_diff_by_k)
# =============================================================================
        
# =============================================================================
#         print('===== finite difference among different gorups =====')
#         print('median', {group: np.median(finite_diff_coef[group]) for group in group_cnt.keys()})
#         print('average', {group: np.mean(finite_diff_coef[group]) for group in group_cnt.keys()})
#         print('max', {group: max(finite_diff_coef[group]) for group in group_cnt.keys()})        
#         print('min', {group: min(finite_diff_coef[group]) for group in group_cnt.keys()})
#         
# =============================================================================
# =============================================================================
#         for t in T:
#             for c in C:
#                 coef = r[c, t] * (np.exp(-alpha * (Z_param[c, t] + 1)) - np.exp(-alpha * Z_param[c, t])) * s[c, t]
#                 if coef != 0:
#                     finite_diff_val[str((counter, c, t))] = r[c, t] * (np.exp(-alpha * (Z_param[c, t] + 1)) - np.exp(-alpha * Z_param[c, t])) * s[c, t]
# =============================================================================
        
        # finite_diff_val = dict(sorted(finite_diff_val.items(), key = lambda item: item[1]))
        
        m.addConstr(lhs[counter] <= Xi, name = 'cut_' + str(counter))

        
# =============================================================================
#         print('number of cuts', len(cuts))
#         for ind, cut in enumerate(cuts):       
#             m.addConstr(cut, name = 'cut_' + str(ind + 1))
# =============================================================================
            
        """ Solving
        """
        # m.NumStart = 0
        m.optimize()

# =============================================================================
#         m.computeIIS()
#         if m.IISMinimal:
#             print('IIS is minimal \n')
#         else:
#             print('IIS is not minimal \n')
#         print('\n The following constraint(s) cannot be satisfied:')
#         
#         infeasible = []
#         for c in m.getConstrs():
#             if c.IISConstr:
#                 print('%s' % c.constrName)
#                 infeasible.append(c.constrName)
#         
#         contr_14_list = [((4, 5),8), ((5, 4),8), ((5, 5),8)]
#         for each in contr_14_list:
#             print('constr_14', each)
#             print(m.getRow(contr_14[each]))
#             print('\n')
#             
#         contr_16_list = [(4,8), (5,8), (3,9), (4,9), (5,9)]
#         for each in contr_16_list:
#             print('constr_16', each)
#             print(m.getRow(contr_16[each]))
#             print('\n')
# =============================================================================

                
                
            
# =============================================================================
#         modified_model = m.copy()
#         orig_num_vars = modified_model.NumVars
#         modified_model.feasRelaxS(0, False, True, False)
#         modified_model.optimize()
# =============================================================================
        
        
# =============================================================================
#         slacks = modified_model.getVars()[orig_num_vars:]
#         for sv in slacks:
#             if sv.X > 1e-9:
#                 print('%s = %g' % (sv.VarName, sv.X))
#         
#         lhs_val[counter]  = lhs[counter].getValue()
#         print('====== lhs after optimization is', lhs_val[counter])
#     
# =============================================================================
# =============================================================================
#         for group in group_cnt.keys():
#             zzz_val[str((counter, group))] = ZZZ[group].X
# =============================================================================
        
# =============================================================================
#         print('========== Z_wt ==========')
#         for group in group_cnt.keys():
#             print(group, ZZZ[group].X)        
#         print('========== Z_wt hat ======')
#         for group in group_cnt.keys():
#             print(group, ZZZ_param[group])
# =============================================================================
        
# =============================================================================
#         for key, val in X.items():
#             if val.X != 0:
#                 x_val[str((counter, key))] = val.X
# =============================================================================
        
        Xi_iter_lb = m.poolObjBound
        Xi_iter_ub = m.objVal

        print('Pi result after optimization')
        print('Xi iter upper', Xi_iter_ub)
        print('Xi iter lower', Xi_iter_lb)

        
        if Xi_iter_lb > Xi_lb:
            Xi_lb = Xi_iter_lb
        
        
        subs = []
        print('====== recalculated Z_ct ======')
        Z_param_prior_opt = Z_param.copy()
        for sub in sub_Z: # for the grouping algorithm, we need to recalculate Z
            # print(sub)
            # if Z_param[sub] != Z[sub].X:
            #    print(sub, Z_param[sub], Z[sub])
            #    subs.append(sub)
            c, t = sub
            Z_param[c, t] = sum(X[c_prime, c, t - 1].X for c_prime in C if is_nearby_cell(c, c_prime))
# =============================================================================
#             if Z_param[c, t] != 0:
#                 z_recov_val[str((counter, c, t))] = Z_param[c, t]
# =============================================================================
                # print(c,t, Z_param[c, t])
        
# =============================================================================
#         second_lhs = f_Z + sum(r[c, t] * (np.exp(-alpha * (Z_param_prior_opt[c, t] + 1)) - np.exp(-alpha * Z_param_prior_opt[c, t])) * s[c, t] * (Z_param[c, t] - Z_param_prior_opt[c, t]) for c in C for t in T)        
# =============================================================================
# =============================================================================
#         print('===== lhs from recalculated Z_ct =====',second_lhs)
# =============================================================================
        
        for group in group_cnt.keys():
            ZZZ_param[group] = ZZZ[group].X
            # print('===== iteration: =====', counter)
# =============================================================================
#             if ZZZ[group].X != 0:
#                 print('===== ZZZ value is', group, ZZZ[group].X)
# =============================================================================
            
# =============================================================================
#         for sub in sub_Z:
#             c,t = sub
#             if Z_param[c, t] != 0:
#                 print('===== recalculated Z_ct is', c,t, Z_param[c, t])
# =============================================================================
# =============================================================================
#         print('checking if Z_param is updated')
#         for sub in subs:
#             print(sub, Z_param[sub])
# =============================================================================
        
            # Z_param[sub] = Z[sub].X
            
            # print(Z_param[sub], Z[sub].X)
        
        print('After optimization')
        print('Xi upper', Xi_ub)
        print('Xi lower', Xi_lb)
        
        counter += 1
        
# =============================================================================
#         Z_param_new = {}
#         for sub in sub_Z:
#             c,t = sub
#             Z_param_new[c, t] = 0
#         stand_algo_2 = [((1,2),1), ((1,3),2), ((1,4),3), ((2,4),4),
#                         ((3,4),5), ((4,4),6), ((5,4),7), ((5,5),8), ((5,5),9)]
# =============================================================================
        
# =============================================================================
#         for c_t in stand_algo_2:
#             Z_param_new[c_t] = 3
#         
#         for sub in sub_Z:
#             c,t = sub
#             if Z_param_new[c, t] != 0:
#                 print('===== standard algo 2 solution is', c,t, Z_param_new[c, t])
# =============================================================================

        
        
        # =============================================================================
        #     for U in m.getVars():
        #         print('%s %g' % (U.varName, U.x))
        # =============================================================================
        # print(U)
        # print(X)
        # =============================================================================
        #     for key, value in Z.items():    
        #         if value.X != 0:
        #             print(key, value.X)
        # =============================================================================
    gap = abs(Xi_ub - Xi_lb) / Xi_lb
    print('Final MIPGap is', gap)
    end_time = time.time()
    running_time = end_time - start_time
    print("Running time is", running_time)
    time_log[ending_time] = [gap, running_time, Xi_ub]
    
print(time_log)
with open('time_log_T8.txt', 'w') as log_result:
    log_result.write(json.dumps(time_log))
# =============================================================================
#     log_result.write(json.dumps(time_log))
#     with open('Grouping_ZZZ.txt', 'w') as log_result:
#         log_result.write(json.dumps(zzz_val))
#     with open('Grouping_Recov_Z.txt', 'w') as log_result:
#         log_result.write(json.dumps(z_recov_val))
#     with open('Grouping_LHS.txt', 'w') as log_result:
#         log_result.write(json.dumps(lhs_val))
#     with open('Grouping_X.txt', 'w') as log_result:
#         log_result.write(json.dumps(x_val))
# =============================================================================
# =============================================================================
#     with open('Grouping_finite_diff.txt', 'w') as log_result:
#         log_result.write(json.dumps(finite_diff_val))
# =============================================================================
