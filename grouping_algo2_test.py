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

ending_time_grid = list(range(7, 16))
ending_time_grid = [5]

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
    ending_time = 5
    print('===========================')
    print('ending time is', ending_time)
    print('===========================')

    grid_size = 5
    ending_time = ending_time
    num_scenario = 1000
    
    """
    Creating set
    """    
    C = [(i, j) for i in range(1, grid_size + 1) for j in range(1, grid_size + 1)]
    T = list(range(1, ending_time + 1))
    T0 = [0] + T
    Omega = list(range(1, num_scenario + 1))
    J = 3
    I = list(range(0, J * ending_time + 1))
    # print('i is', I)
    
    xx = {}
    for c in C:
        for t in T0:
            if c == (1, 1) and t == 0:
                xx[c, t] = J
            else:
                xx[c, t] = 0
        
# =============================================================================
#     q_raw.columns = [int(col) for col in q_raw.columns]
#     sub_q = list(product(C, T))
#     q = {}
#     for sub in sub_q:
#         c_two_dim, t = sub
#         c_one_dim = (c_two_dim[0] - 1) * grid_size + c_two_dim[1]
#         q[c_two_dim, t] = q_raw.loc[c_one_dim, t]
# =============================================================================
    
    
    
    
    
# =============================================================================
#     # zeta_raw = pd.read_csv(data_folder + '\Zeta.csv', header = None, index_col = 0)
#     
#     Zeta = {}
#     for path in range(1, zeta_raw.shape[0] + 1):
#         for t in range(1, ending_time + 1):
#             all_cells = C
#             for cell in all_cells:
#                 Zeta[(cell, t, path)] = 0 # set Zeta equal to 0
#             cell_one_dim = zeta_raw.loc[path, 3 * (t - 1) + 1] # extract the occupied cell loc from Zeyu's path
#             cell_two_dim = (cell_one_dim // grid_size + 1, np.mod(cell_one_dim, grid_size))
#             Zeta[(cell_two_dim, t, path)] = 1 # set Zeta equal to 1 for occupied celll
# =============================================================================
    
# =============================================================================
#     W = {}
#     for c in C:
#         for t in T:
#             if sum(Zeta[c, t, omega] for omega in Omega) >= 1:
#                 W[c, t] = 1
#             else:
#                 W[c, t] = 0
# =============================================================================
        
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
    
    num_groups = 11
    # max_W_param = max(W_param.values())
    
    # W_param_excl_0 = [val for val in W_param.values() if val != 0]
    # W_param_excl_0 = [val for val in W_param.values()]
    # min_W_param = min(W_param_excl_0)
    cat_group = {}
    W_param_cut = [0, 0.02, 0.04, 0.06, 0.09, 0.12, 0.1240, 0.1455, 0.1757, 0.22, 0.2881]
    for c in C:
        for t in T:
            if W_param[c, t] == 0:
                cat_group[c, t] = 1 
            else:
                for ind in range(0, len(W_param_cut) - 1):
                    if W_param_cut[ind] < W_param[c, t] <= W_param_cut[ind + 1]:
                        cat_group[c, t] = ind + 2 # assign a group for each c,t
    assert len(cat_group.keys()) == len(W_param.keys()), "Not all c,t pairs are grouped"
    
    group_prob = {}
    for c in C:
        for t in T:
            group = cat_group[c, t]
            if group in group_prob.keys():
                group_prob[group].append(W_param[c, t])
            else:
                group_prob[group] = [W_param[c,t]]
    
    group_prob_diff = {}
    
    for group in group_prob.keys():
        group_prob_diff[group] = [min(group_prob[group]), max(group_prob[group])]
        group_prob_diff[group].append(group_prob_diff[group][1] - group_prob_diff[group][0])
        
    print("===== group probability difference =====", group_prob_diff)
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
    
    # add variables
    X = m.addVars(sub_X, lb = 0, name = 'X')
    # Z = m.addVars(sub_Z, lb = 0, ub = J, vtype = GRB.INTEGER, name = 'Z')
    
    group_cnt = {}
    for group in cat_group.values():
        if group not in group_cnt.keys():
            group_cnt[group]  = 1
        else:
            group_cnt[group] += 1
            
    print(' ===== number of cells per group =====', group_cnt)
    
    ZZZ = {} # dict to save ZZZ variable
    ZZZ_param = {}
    for group in group_cnt.keys():
        ZZZ[group] = m.addVar(lb = 0, ub = group_cnt[group], vtype = GRB.INTEGER)
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
    m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == sum(X[c, c_prime, t] for c_prime in C if is_nearby_cell(c, c_prime))  for c in C for t in T), name = '14') #2d
    m.addConstrs((sum(X[c, c_prime, 0] for c_prime in C if is_nearby_cell(c, c_prime)) == xx[c, 0] for c in C), name = '15') #2d
    
    # cell_by_group 
    
    for group in group_cnt.keys():
        # group = 1
        m.addConstr((sum(X[c_prime, c, t - 1] for c in C for t in T for c_prime in C if cat_group[c, t] == group and is_nearby_cell(c, c_prime)) == ZZZ[group]), name = '16_' + str(group)) #2d
    
    # m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == Z[c, t] for c in C for t in T), name = '16') #2d
    lhs = {}
    lhs_val = {}
    zzz_val = {}
    x_val = {}
    z_recov_val = {}
    
    start_time = time.time()
    # while Xi_ub - Xi_lb > delta * Xi_lb and counter <= 100:
    while counter <= 3 and Xi_ub - Xi_lb > delta * Xi_lb and time.time() - start_time <= 900:
        
        ################ step 1 ################
        print('=============', counter, '===============')
        
        r = {}
        for t in T:
            for c in C:
                if t == 1:
                    r[c, t] = p[c]
                else:
                    r[c, t] = sum([r[c_prime, t - 1] * np.exp(-alpha * Z_param[c_prime, t - 1]) * gamma[c_prime, c, t - 1] for c_prime in C])
        
        s = {}
        for t in T[::-1]:
            for c in C:
                if t == ending_time:
                    s[c, t] = 1
                else:
                    s[c, t] = sum([s[c_prime, t + 1] * np.exp(-alpha * Z_param[c_prime, t + 1]) * gamma[c, c_prime, t] for c_prime in C])
    
        # f_Z = sum([r[c, 1] * np.exp(-alpha * Z_param[c, 1]) * s[c, 1] for c in C])
        # f_Z_2 = sum([r[c, ending_time] * np.exp(-alpha * Z_param[c, ending_time]) * s[c, ending_time] for c in C])
        f_Z = sum([r[c, 5] * np.exp(-alpha * Z_param[c, 5]) * s[c, 5] for c in C])

        print('f(Z) equals to', f_Z)
        # print('f(Z) equals to', f_Z_2)
        # print('f(Z) equals to', f_Z_5)
        
        if f_Z < Xi_ub:
            Xi_ub = f_Z   
        if Xi_ub - Xi_lb <= delta * Xi_lb:
            break
        
        print('Before optimization')
        print('Xi upper', Xi_ub)
        print('Xi lower', Xi_lb)
        

        ################ step 2 ##############
        # def solve_p():

        if Xi_lb == 0:
            g = np.inf
        else:
            g = (Xi_ub - Xi_lb) / Xi_lb
        print('counter is', counter)
        print('g is', g)
        
        if counter == 1:
            mip_gap = 0
        else:
            mip_gap = min([0.03, g / 3])
        # mip_gap = 0.1 / 2 ** (counter - 1)
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
        lhs[counter] = f_Z + sum([r[c, t] * (np.exp(-alpha * (Z_param[c, t] + 1)) - np.exp(-alpha * Z_param[c, t])) * s[c, t] * (ZZZ[cat_group[c, t]] - Z_param[c, t]) for c in C for t in T]) 
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
        
        lhs_val[counter]  = lhs[counter].getValue()
        
        for group in range(1, num_groups + 1):
            zzz_val[counter, group] = ZZZ[group].X
        
        for key, val in X.items():
            x_val[counter, key] = val.X
        
        Xi_iter_lb = m.poolObjBound
        Xi_iter_ub = m.objVal

        print('Pi result after optimization')
        print('Xi iter upper', Xi_iter_ub)
        print('Xi iter lower', Xi_iter_lb)

        
        if Xi_iter_lb > Xi_lb:
            Xi_lb = Xi_iter_lb
        
        
        subs = []
        for sub in sub_Z: # for the grouping algorithm, we need to recalculate Z
            # print(sub)
            # if Z_param[sub] != Z[sub].X:
            #    print(sub, Z_param[sub], Z[sub])
            #    subs.append(sub)
            c, t = sub
            Z_param[c, t] = sum(X[c_prime, c, t - 1].X for c_prime in C if is_nearby_cell(c, c_prime))
        z_recov_val[counter, c, t] = Z_param[c, t]
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
    gap = (Xi_ub - Xi_lb) / Xi_lb
    print('Final MIPGap is', gap)
    end_time = time.time()
    running_time = end_time - start_time
    print("Running time is", running_time)
    time_log[ending_time] = [gap, running_time, Xi_ub]
    
print(time_log)
with open('time_log_T8.txt', 'w') as log_result:
    log_result.write(json.dumps(time_log))
with open('ZZZ.txt', 'w') as log_result:
    log_result.write(json.dumps(zzz_val))
with open('Recov_Z.txt', 'w') as log_result:
    log_result.write(json.dumps(z_recov_val))
with open('LHS.txt', 'w') as log_result:
    log_result.write(json.dumps(lhs_val))
with open('X.txt', 'w') as log_result:
    log_result.write(json.dumps(x_val))
    