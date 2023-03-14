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
#
"""
Define helper functions
"""

def is_forward_cell(c_0, c_1, 
                    starting_c = (-100, -100), ending_c = (-100, 100),
                    on_map_start = (1, 1), on_map_end = (1, 5)):
    """ function to check if c_1 is a forward cell for c_0 """ 
    
    """ add special condition for starting state """
    if c_0 == on_map_end: # add additional ending_c
        if c_1 == ending_c:
            return True
        elif c_0 == c_1:
            return True
        elif c_0[0] == c_1[0]:
            if c_0[1] + 1 == c_1[1] or c_0[1] - 1 == c_1[1]:
                return True
        elif c_0[1] == c_1[1]:
            if c_0[0] + 1 == c_1[0] or c_0[0] - 1 == c_1[0]:
                return True
        else:
            return False
    elif c_0 == starting_c: # only itself and on_map_start
        if c_1 == on_map_start or c_1 == starting_c: # the first cell has to be (1, 1)
            return True
        else:
            return False
    elif c_0 == ending_c: # only itself
        if c_1 == ending_c:
            return True
        else:
            return False
    else: # normall cells
        if c_0 == c_1:
            return True
        elif c_0[0] == c_1[0]:
            if c_0[1] + 1 == c_1[1] or c_0[1] - 1 == c_1[1]:
                return True
        elif c_0[1] == c_1[1]:
            if c_0[0] + 1 == c_1[0] or c_0[0] - 1 == c_1[0]:
                return True
        else:
            return False
        
def is_backward_cell(c_0, c_1, 
                     starting_c = (-100, -100), ending_c = (-100, 100),
                     on_map_start = (1, 1), on_map_end = (1, 5)):
    """ function to check if c_0 is a backward cell for c_1 """ 
    
    """ add special condition for starting state """
    if c_1 == on_map_start:
        if c_0 == starting_c: # add starting_c as the backward
            return True
        if c_0 == c_1:
            return True
        elif c_0[0] == c_1[0]:
            if c_0[1] + 1 == c_1[1] or c_0[1] - 1 == c_1[1]:
                return True
        elif c_0[1] == c_1[1]:
            if c_0[0] + 1 == c_1[0] or c_0[0] - 1 == c_1[0]:
                return True
        return False           
    elif c_1 == ending_c: # the backward state of the ending cell is the (1, 5)
        if c_0 == ending_c:
            return True
        elif c_0 == on_map_end:
            return True
        else:
            return False
    elif c_1 == starting_c: # off map starting, only itself is the backward 
        if c_0 == starting_c:
            return True
        else:
            return False        
    else:
        if c_0 == c_1:
            return True
        elif c_0[0] == c_1[0]:
            if c_0[1] + 1 == c_1[1] or c_0[1] - 1 == c_1[1]:
                return True
        elif c_0[1] == c_1[1]:
            if c_0[0] + 1 == c_1[0] or c_0[0] - 1 == c_1[0]:
                return True
        else:
            return False        

def is_backward_state_with_cam(s_c_0, s_c_1):
    """ function to check if s_c_0 is a backward state for s_c_1 """ 
    """ this function is for the target which state s_c contains s and cam binary indicator """
    s_0, c_0 = s_c_0
    s_1, c_1 = s_c_1
    
    if s_0 == s_1: # no matter if c_0 is camou or not, s_0 and s_1 are connected
        return True
    else:
        if c_0 == 0 and c_1 == 0: # current and next cam are 0;
            if s_0[0] == s_1[0]: # row is the same
                if s_0[1] + 1 == s_1[1] or s_0[1] - 1 == s_1[1]:
                    return True
            elif s_0[1] == s_1[1]: # column is the same
                if s_0[0] + 1 == s_1[0] or s_0[0] - 1 == s_1[0]:
                    return True
        else:
            return False       



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

grid_size = 9
# ending_time_grid = [10, 12, 14, 15, 16, 17, 18, 20]
ending_time_grid = [10, 12, 14, 15, 16, 17, 18, 20]
# ending_time = 10
# ending_time = 15
# num_scenario = 1000
J = 3
J_2 = int(J * 0.7)
J_1 = J - J_2

# ending_time = 10
for ending_time in ending_time_grid:
    ending_time = ending_time
    # ending_time = 9
    print('=============================')
    print('Ending time is', ending_time)
    print('=============================')

    """
    Creating set
    """
    
    S = [(i, j) for i in range(1, grid_size + 1) for j in range(1, grid_size + 1)] # a state for target
    T = list(range(1, ending_time + 1)) # time period during which a detection could happen
    T0 = [0] + T # time period add t0
    # Omega_num = list(range(1, num_scenario + 1)) # numerical list with each element represent a path number
    C = [0, 1]
    S_C = list(product(S, C)) # state contains (c, camouflage_state)
    L = [1, 2] # set of searchers' type
    # L = [1]
    n_L = {1: J_1, 2: J_2} # number of searchers for each searcher type
    # n_L = {1: 1}
    # alpha_l = {1: 0.1, 2: 0.2} # detection rate for each searcher type
    # I = list(range(0, J * ending_time + 1))
    # print('i is', I)
    J_ind = {1: list(range(1, J_1 + 1)), 2: list(range(1, J_2 + 1))}
    
    
    total_J = sum(n_L.values())
    s_init = (-100, -100)
    s_end = (-100, 100)
    on_map_init = (grid_size // 2, 1)
    on_map_end = (grid_size, grid_size // 2)
    S_expand = S + [s_init] + [s_end]
    tau = {1: ending_time * 0.8, 2: ending_time * 0.6} # operation duration limit for searcher of type l
    
    n = {} # number of seachers per state per time
    for s in S_expand:
        for t in T:
            n[s, t] = total_J
            
    """ Import data
    """
# =============================================================================
#     if platform.system() == 'Windows':
#         data_folder = 'E:\\Research\\Route_Optimization_for_Multiple_Searchers\\Python\\'
#     else:
#         data_folder = os.path.dirname(os.path.realpath(__file__))
# =============================================================================
    

    
    """
    Creating parameters
    """
     
    """ p is where target is at t = 1"""
    p = {}
    for s_c in S_C:
        if s_c == (((grid_size + 1) / 2, (grid_size + 1) / 2), 0):
            p[s_c] = 1
        else:
            p[s_c] = 0
    
    alpha = {}
    for l in L:
        for c in C:
            if c == 1:
                alpha[l, c] = 0
            else:
                alpha[l, c] = -3 * np.log(0.4) / n_L[l]
# =============================================================================
#     alpha_sub = list(product(S, T))
#     alpha_value = [-3 * np.log(0.4) / J] * grid_size * grid_size * ending_time
#     alpha = dict(zip(alpha_sub, alpha_value))
# =============================================================================
    
    sub_gamma = list(product(S_C, S_C, T))
    gamma = {}
    stay_prob = 0.5
    to_cam = 0.1
    from_cam = 5 / 6
    stay_cam = 1 / 6
    
    for sub in sub_gamma:
        s, s_prime, t = sub
        
        s_c, s_cam = s[0], s[1]
        s_prime_c, s_prime_cam = s_prime[0], s_prime[1]
        
        if s_cam == 1: # target is in camouflage state
            if s_c == s_prime_c and s_prime_cam == 1: # stay in the same cell and remain in camouflage
                gamma[sub] = stay_cam
            elif s_c == s_prime_c and s_prime_cam == 0: # stay in the same cell and comes out of camouflage
                gamma[sub] = from_cam
            # do not define gamma for other state
            else:
                gamma[sub] = 0 # define other prob as 0; we might be able to remove this.
        else: # target is not in camoflage state
            if s_prime_cam == 1: # target move from non-camouflage to camouflage
                if s_c == s_prime_c:
                    gamma[sub] = to_cam
            else: # target remains in a non-camouflage state        
                non_cam_prob = 1 - to_cam
                if s_c == s_prime_c:
                    gamma[sub] = stay_prob
                elif is_corner_cell(s_c, grid_size = grid_size):
                    if is_nearby_cell(s_c, s_prime_c):
                        gamma[sub] = (non_cam_prob - stay_prob) / 2
                    else:
                        gamma[sub] = 0
                elif is_side_cell(s_c, grid_size = grid_size):
                    if is_nearby_cell(s_c, s_prime_c):
                        gamma[sub] = (non_cam_prob - stay_prob) / 3
                    else:
                        gamma[sub] = 0
                else:
                    if is_nearby_cell(s_c, s_prime_c):
                        gamma[sub] = (non_cam_prob - stay_prob) / 4
                    else:
                        gamma[sub] = 0
                
    """
    the following code is to help with debugging only
    To check if gamma is correctly calculated
    """
# =============================================================================
#     test_t = 5
#     test_state = ((3, 3), 0)
#     end_state_set = list(product(S, C))           
#     
#     for key in gamma.keys():
#         state_0, state_1, t = key
#         if state_0 == test_state and t == test_t:
#             if gamma[state_0, state_1, t] != 0:
#                 print(state_0, ',', state_1, ',', t, 'prob', gamma[state_0, state_1, t])
# =============================================================================
    
# =============================================================================
#     start_state = ((1, 1), 0)
#     for s_prime in S:
#         if start_state[1] == 1:
#             if gamma[(start_state, s_prime, 1)] != 0:
#                 print(start_state, s_prime, 1, 'gamma value', gamma[(start_state, s_prime, 1)])
#         else:
#             s_prime_c, s_prime_cam = s_prime[0], s_prime[1]
#             s_c, s_cam = start_state[0], start_state[1]
#             if is_nearby_cell(s_c, s_prime_c) and s_prime[1] == 0:
#                 print(start_state, s_prime, 1, 'gamma value', gamma[(start_state, s_prime, 1)])
#             elif s_c == s_prime_c and s_prime[1] == 1:
#                 print(start_state, s_prime, 1, 'gamma value', gamma[(start_state, s_prime, 1)])
# =============================================================================
# =============================================================================
#         for s_c in S_C:
#             if p[s_c] != 0:
#                 print(s_c, p[s_c])
# =============================================================================
    
    sub_q = list(product(S_C, T)) 
    sub_q = sorted(sub_q, key = lambda x: x[1])
    q = {} # create param values for q
    for sub in sub_q:
        s_c, t = sub
        # print(sub)
        if t == 1:
            q[s_c, t] = p[s_c]
        else:
            # c_one_dim = (c_two_dim[0] - 1) * grid_size + c_two_dim[1]
            s, cam = s_c[0], s_c[1]
            if cam == 1: # if target is in camouflage at t, it either in the same state at t - 1 or change into camouflage from the same cell
                q[s_c, t] = q[(s, 1), t - 1] * stay_cam + q[(s, 0), t - 1] * to_cam
            else: # if target is not in camouflage at t, it either moves from a non-camouflage cell or changes back a camouflage state at the same cell
                q[s_c, t] = q[(s, 1), t - 1] * from_cam + sum([q[(s_prime, 0), t - 1] * gamma[(s_prime, 0), (s, 0), t - 1] for s_prime in S if is_nearby_cell(s, s_prime)])
    
    """ code to check if a state is reachable for a searcher """
    searcher_occ = is_searcher_occ(C = S, T = T, grid_size = grid_size)
    occ_cell = list(searcher_occ.keys())
    
    """ init and occ cell are occupiable """
    init_occ = list(product([s_init], T))
    end_occ = list(product([s_end], T))
    # occ_cell.extend(init_occ)
    # occ_cell.extend(end_occ)
    
    searcher_reachable = {}
    for t in T:
        for s in S_expand:
            if s == s_init or s == s_end:
                searcher_reachable[s, t] = 1
            elif (s, t) not in occ_cell:
                searcher_reachable[s, t] = 0
            else:
                searcher_reachable[s, t] = 1
    
    """ creating set T that has no detection across all states
    """
    
    T_no_detect = []
    for t in T:
        if sum(q[(s, c), t] * searcher_reachable[s, t] for s, c in S_C if c == 0) == 0:
            T_no_detect.append(t)
    
    
    """ creating set (s, t) that has no detection
    """
    part_T_no_detect = {}
    for t in set(T) - set(T_no_detect):
        for s in S:
            if q[(s, 0), t] * searcher_reachable[s, t] == 0:
                if t not in part_T_no_detect.keys():
                    part_T_no_detect[t] = [(s, t)]
                else:
                    part_T_no_detect[t].append((s, t))
    
    """ s_init and s_end for t should be non-detectable """
    for t in set(T) - set(T_no_detect):
        part_T_no_detect[t].append((s_init, t))
        part_T_no_detect[t].append((s_end, t))
    
    part_T_no_detect_pairs = sum(part_T_no_detect.values(), [])

# =============================================================================
#     xx = {} 
#     for c in C:
#         if c == (1, 1):
#             xx[c] = J
#         else:
#             xx[c] = 0
# =============================================================================
    x = {}
    for l in L:
        for s in S_expand:
            for t in T0:
                if s in [s_init] and t == 0:
                    x[l, s, t] = n_L[l]
                else:
                    x[l, s, t] = 0
    
    # searcher index based on searcher type
    J_l = {}
    for l in L:
        J_l[l] = list(range(1, n_L[l] + 1))
    
    

    
    Xi_lb = 0
    Xi_ub = 1
    counter = 1
    
    cuts = []
    
    """
    Defining decision variables
    """
    """
    """
    model_name = 'sp1_lm'
    m = gp.Model(model_name)
    m.setParam(GRB.Param.TimeLimit, 15 * 60)
    m.setParam(GRB.Param.Threads, 1)
    m.setParam(GRB.Param.LogFile, model_name)
    # m.setParam(GRB.Param.MIPGap, 1e-5)
    
    sub_X = list(product(L, S_expand, S_expand, T0))
    sub_Z = list(product(L, S_expand, T))
    sub_O = list(product(L, T))
    
    
    delta = 1e-4

    
    X = m.addVars(sub_X, lb = 0, name = 'X')
    # Z = m.addVars(sub_Z, vtype = GRB.INTEGER, lb = 0, name = 'Z')
# =============================================================================
#     for l in L:
#         for s_expand in S_expand:
#             for t in T:
#                 Z[l, s_expand, t].ub = min(n_L[l], n[s, t])
# =============================================================================
    
    O = m.addVars(sub_O, vtype = GRB.INTEGER, lb = 0, name = 'O')
    # O = m.addVars(sub_O, lb = 0, name = 'O')    
    for l in L:
        for t in T:
            O[l, t].ub = n_L[l]
    
    
    Z_param = {} # Z_param is for recover the Z value post optimization
    for sub in sub_Z:
        Z_param[sub] = 0
    
    """ ZZZ is only defined for 'detectable' state """
    sub_ZZZ = [(l, s, t) for l in L for s in S_expand for t in T if t not in T_no_detect and (s, t) not in part_T_no_detect_pairs]
    sub_ZZZ = sorted(sub_ZZZ, key = lambda x: (x[0], x[2]))
    
    
    ZZZ = {} # dict to save ZZZ variable
    ZZZ_param = {}
    for sub in sub_ZZZ:
        l, s, t = sub
        ZZZ[sub] = m.addVar(lb = 0, ub = min(n_L[l], n[s, t]), vtype = GRB.INTEGER)
        ZZZ[sub].Start = 0
        ZZZ_param[sub] = 0
    
    Xi = m.addVar(vtype = GRB.CONTINUOUS, name = 'Xi')

    """
    Defining objective function
    """
    m.setObjective(Xi, GRB.MINIMIZE)    
    # m.setObjective(3, GRB.MINIMIZE)    
  
    """
    Defining constraints
    """
# =============================================================================
#     m.addConstrs((sum(X[l, s_prime, s, t - 1] for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
#                   == Z[l, s, t] for l in L for s in S_expand for t in T), name = 'link_x_z') #2d
# =============================================================================
    
    m.addConstrs((sum(X[l, s_prime, s, t - 1] for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
                  == ZZZ[l, s, t] for l in L for s in S_expand for t in T if t not in T_no_detect and (s, t) not in part_T_no_detect_pairs), name = 'link_x_z') #2d
    
    """ variables and constraints with t that has no detection on any cell """
    Z_0 = m.addVars(list(product(L, T_no_detect)), name = 'Z_0', vtype = GRB.CONTINUOUS)
    for l in L:
        for t in T_no_detect:
            Z_0[l, t].ub = n_L[l]
            Z_0[l, t].lb = n_L[l]
            
    for l in L:
        for t in T_no_detect:
            constr_link_x_Z_0 = m.addConstr((sum(X[l, s_prime, s, t - 1] for s in S_expand for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end)) == Z_0[l, t]), name = '16_link_Z_0_' + str(l) + str(t)) 
    
    """ variables and constraints with no detection pairs """
    part_Z_0 = m.addVars(list(part_T_no_detect.keys()), name = 'part_Z_0', vtype = GRB.CONTINUOUS)
    
    part_detect_T = set(T) - set(T_no_detect)
    for t_loop in part_detect_T:
        # print(t_loop)
        constr_link_x_part_Z_0 = m.addConstr((sum(X[l, s_prime, s, t_loop - 1] for s, t in part_T_no_detect[t_loop] for s_prime in S_expand if searcher_reachable[s, t_loop] > 0 and is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end)) == part_Z_0[t_loop]), name = '16_link_part_Z_0_' + str(t_loop)) 
    
# =============================================================================
#     constr_Z_t_J = {}
#     for l in L:
#         for t in part_detect_T:
#             constr_Z_t_J[l, t] = m.addConstr((part_Z_0[t] + sum(ZZZ[l, s, t] for l, s, t_inner in sub_ZZZ if t_inner == t) >= n_L[l]), name = 'constr_' + str(l) + '_' + str(t))
# 
# =============================================================================
# =============================================================================
#         for l in L:
#             for t in part_T_no_detect.keys():
#                 m.addConstr((part_Z_0[l, t] + sum([ZZZ[l, group] for group in group_cnt.keys() if group[1] == t]) == n_L[l]), name = '100_' + str(l) + str(group)) #2d
# =============================================================================
        

    """ adding 2.4b to 2.4f     
    """
    # 2.4b
    m.addConstrs((sum(X[l, s_prime, s, t - 1] for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
                  == sum(X[l, s, s_prime, t] for s_prime in S_expand if is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))  for l in L for s in S_expand for t in T), name = 'continum') #2d
    
    # 2.4c
    """
    For this specific constraint, s is the default starting state at t = 0, which is defined by the number of searchers
    through xx[l, s, 0]. And s_prime (which is where the searchers could be at t = 1) need to be defined in a separate
    way from the general is_nearby_cell function
    """
    m.addConstrs((sum(X[l, s, s_prime, 0] for s_prime in S_expand if is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
                  == x[l, s, 0] for l in L for s in S_expand), name = 'continum_start') #2d
    
# =============================================================================
#     for s in S_expand:
#         if x[1, s, 0] == 1:
#                 # print(1, s, 0, x[1, s, 0])
#                 # print('s is', s, 's_prime is', s_prime)
#             print(s)
#             for s_prime in S_expand:
#                 if X[1, s, s_prime, 0].X == 1 and is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end):
#                     print(1, s, s_prime, 0, X[1, s, s_prime, 0].X)  
#                 
#     for s in S_expand:
#         if x[1, s, 0] != 0:
#             print(1, s, 0, x[1, s, 0])
#     
#     for s in S_expand:
#         for s_prime in S_expand:
#             if X[1, s, s_prime, 0].X != 0:
#                 print(1, s, s_prime, 0, ':', X[1, s, s_prime, 0].X)
# =============================================================================
    # m.addConstrs((sum(X[l, s, s_prime, 0] for s_prime in S_expand if is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end)) == x[l, s, 0] for l in L for s in S), name = 'continum_start') #2d

    # m.addConstrs((sum(X[l, s, s_prime, 0] for s_prime in S if s == s_prime) == x[l, s, 0] for l in L for s in S), name = 'continum_start') #2d

    """ 2.4d """
    m.addConstrs((sum(X[l, s_expand_prime, s_expand, t - 1] for s_expand_prime in S_expand for l in L if is_backward_cell(s_expand_prime, s_expand, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
                  <= n[s, t] for s_expand in S_expand for t in T), name = 'state_capacity') #2d
    
    # 2.4e
    """ 
    this s can only be at the searcher's init position 
    
    """
    m.addConstrs((sum(X[l, s, s_prime, t - 1] for s in [s_init] for s_prime in S_expand if s_prime != s and is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
                  == O[l, t] for l in L for t in T), name = 'duration_start') #2d

    # 2.4f
    """ should we consider those move off the map next time period """
    m.addConstrs((sum(X[l, s, s_prime, t] for s in S for s_prime in S_expand if is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end)) 
                  <= sum(O[l, t_prime] for t_prime in T if t_prime <= t and t_prime >= t - tau[l] + 1)  for l in L for t in T), name = 'duration') #2d

    # m.optimize()
# =============================================================================
#     m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == sum(X[c, c_prime, t] for c_prime in C if is_nearby_cell(c, c_prime))  for c in C for t in T), name = '14') #2d
#     
#     m.addConstrs((sum(X[c, c_prime, 0] for c_prime in C if is_nearby_cell(c, c_prime)) == xx[c] for c in C), name = '15') #2d
# =============================================================================
    
    start_time = time.time()
    # while Xi_ub - Xi_lb > delta * Xi_lb and counter <= 100:
    while Xi_ub - Xi_lb > delta * Xi_lb and time.time() - start_time <= 900:
    # while counter <= 5 and Xi_ub - Xi_lb > delta * Xi_lb and time.time() - start_time <= 900:

        ################ step 1 ################
        print('=============', counter, '===============')
        
        r = {}
        for t in T:
            for s_c in S_C:
                if t == 1:
                    r[s_c, t] = p[s_c]
                else:
                    r[s_c, t] = sum([r[s_c_prime, t - 1] * np.exp(-sum(alpha[l, s_c[1]] * Z_param[l, s_c_prime[0], t - 1] for l in L)) * gamma[s_c_prime, s_c, t - 1] for s_c_prime in S_C if is_backward_state_with_cam(s_c_prime, s_c)])
                    # r[s_c, t] = sum([r[s_c_prime, t - 1] * np.exp(-sum(alpha[l, s_c[1]] * Z_param[l, s_c_prime[0], t - 1] for l in L)) * gamma[s_c_prime, s_c, t - 1] for s_c_prime in S_C])

        s = {}
        for t in T[::-1]:
            for s_c in S_C:
                if t == ending_time:
                    s[s_c, t] = 1
                else:
                    s[s_c, t] = sum([s[s_c_prime, t + 1] * np.exp(-sum(alpha[l, s_c[1]] * Z_param[l, s_c_prime[0], t + 1] for l in L)) * gamma[s_c, s_c_prime, t] for s_c_prime in S_C if is_backward_state_with_cam(s_c, s_c_prime)]) # if s_c_prime is s_c's forward
                    # s[s_c, t] = sum([s[s_c_prime, t + 1] * np.exp(-sum(alpha[l, s_c[1]] * Z_param[l, s_c_prime[0], t + 1] for l in L)) * gamma[s_c, s_c_prime, t] for s_c_prime in S_C]) # if s_c_prime is s_c's forward

        # f_Z = sum([r[c, 1] * np.exp(-alpha * Z_param[c, 1]) * s[c, 1] for c in C])
        # f_Z_2 = sum([r[c, ending_time] * np.exp(-alpha * Z_param[c, ending_time]) * s[c, ending_time] for c in C])
        f_Z = sum([r[s_c, 5] * np.exp(-sum(alpha[l, s_c[1]] * Z_param[l, s_c[0], 5] for l in L)) * s[s_c, 5] for s_c in S_C])

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

        
        g = (Xi_ub - Xi_lb) / Xi_lb if Xi_lb != 0 else np.inf
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
        # m.addConstr(f_Z + sum([r[s_c, t] * (np.exp(-sum(alpha[l, s_c[1]] * (Z_param[l, s_c[0], t] + 1) for l in L)) - np.exp(-sum(alpha[l, s_c[1]] * (Z_param[l, s_c[0], t]) for l in L))) * s[s_c, t] * (Z[s_c, t] - Z_param[s_c, t]) for s_c in S_C for t in T]) <= Xi, name = 'cut_' + str(counter))
        
        # adj_finite_diff = {}
        """
        use 1 - s_c[1] to adjust the finite difference: when the searcher is in the camouflage mode, the finite difference will be 0
        """
        adj_finite_diff = {}
        for l in L:
            for s_c in S_C:
                for t in T:
                  adj_finite_diff[l, s_c, t] = (1 - s_c[1]) * r[s_c, t] * (np.exp(-alpha[l, s_c[1]] * (Z_param[l, s_c[0], t] + 1)) - np.exp(-alpha[l, s_c[1]] * (Z_param[l, s_c[0], t]))) * s[s_c, t]
                  
        m.addConstr(f_Z + sum([adj_finite_diff[l, s_c, t] * (ZZZ[l, s_c[0], t] - ZZZ_param[l, s_c[0], t]) for l in L for s_c in S_C for t in T if t not in T_no_detect and (s_c[0], t) not in part_T_no_detect_pairs]) <= Xi, name = 'cut_' + str(counter))

        # m.addConstr(f_Z + sum([(1 - s_c[1]) * r[s_c, t] * (np.exp(-alpha[l, s_c[1]] * (Z_param[l, s_c[0], t] + 1)) - np.exp(-alpha[l, s_c[1]] * (Z_param[l, s_c[0], t]))) * s[s_c, t] * (Z[l, s_c[0], t] - Z_param[l, s_c[0], t]) for l in L for s_c in S_C for t in T]) <= Xi, name = 'cut_' + str(counter))

        
# =============================================================================
#         print('number of cuts', len(cuts))
#         for ind, cut in enumerate(cuts):       
#             m.addConstr(cut, name = 'cut_' + str(ind + 1))
# =============================================================================
            
        """ Solving
        """
        # m.NumStart = 0
        m.optimize()
        
        Xi_iter_lb = m.poolObjBound
        Xi_iter_ub = m.objVal

        print('Pi result after optimization')
        print('Xi iter upper', Xi_iter_ub)
        print('Xi iter lower', Xi_iter_lb)

        
        if Xi_iter_lb > Xi_lb:
            Xi_lb = Xi_iter_lb
        
        
        """ updating defined ZZZ """
        """ ZZZ_param is used to calculate the new cut """
        for sub in sub_ZZZ:
            ZZZ_param[sub] = ZZZ[sub].X
        
        """ updating Z param """
        """ Z param is used to calculate f(Z), r and s """
        """ Z param is defined at all s, t """
        
        # subs = []
        print('====== recalculated Z_ct ======')
        Z_param_prior_opt = Z_param.copy()
        for sub in sub_Z: # for the grouping algorithm, we need to recalculate Z
            # print(sub)
            # if Z_param[sub] != Z[sub].X:
            #    print(sub, Z_param[sub], Z[sub])
            #    subs.append(sub)
            l, s, t = sub
            # Z_param[l, s, t] = sum(X[l, s_prime, s, t - 1].X for s_prime in S_expand if is_nearby_cell(s, s_prime))
            Z_param[l, s, t] = sum(X[l, s_prime, s, t - 1].X for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))

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
    # time_log[grid_size] = [gap, running_time, Xi_ub]
# print(time_log)
# with open('time_log_T10.txt', 'w') as log_result:
#    log_result.write(json.dumps(time_log)) 


    

