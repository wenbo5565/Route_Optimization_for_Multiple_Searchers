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
    Defining several helper functions
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
                nearby_cell = return_nearby_cell((4, 1), grid_size)
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

""" End of helper function definition """

""" Start of the optimization problem  formulation """

grid_size = 9
# ending_time_grid = [10, 12, 14, 15, 16, 17, 18, 20]
# ending_time_grid = [10, 12, 14, 15] #, 16, 17, 18, 20] # , 15] # , 12, 14]
ending_time_grid = [10, 12, 14, 15]
# ending_time_grid = [16, 17, 18, 20]
# ending_time_grid = [16]

# ending_time = 10
# ending_time = 16
# num_scenario = 1000
J = 15
#
J_2 = int(J * 0.7)
J_1 = J - J_2

# ending_time = 10
# ending_time = 14
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
    C = [0, 1]
    S_C = list(product(S, C)) # state contains (c, camouflage_state)
    L = [1, 2] # set of searchers' type
    J_L = {1: J_1, 2: J_2} # number of searchers for each searcher type
    J_ind = {1: list(range(1, J_1 + 1)), 2: list(range(1, J_2 + 1))} # index for each type l in L
    
    
    total_J = sum(J_L.values())
    s_init = (-100, -100)
    s_end = (-100, 100)
    on_map_init = (grid_size // 2, 1)
    on_map_end = (grid_size, grid_size // 2)
    S_expand = S + [s_init] + [s_end]
    tau = {1: int(ending_time * 0.8), 2: int(ending_time * 0.6)} # operation duration limit for searcher of type l
    
    n = {} # number of seachers per state per time
    for s in S_expand:
        for t in T:
            n[s, t] = total_J
        
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
                alpha[l, c] = -3 * np.log(0.4) / total_J
    
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
                        
    """ calculating value for parameter q """
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
    
    """ calculating value for parameter x """
    x = {}
    for l in L:
        for s in S_expand:
            for t in T0:
                if s in [s_init] and t == 0:
                    x[l, s, t] = J_L[l]
                else:
                    x[l, s, t] = 0
                    
    """ calculating if a (s, t) pair is detectable """
    """ code to check if a state is reachable for a searcher """
    searcher_occ = is_searcher_occ(C = S, T = T, grid_size = grid_size)
    occ_cell = list(searcher_occ.keys())
    
    """ init and occ cell are occupiable """
    init_occ = list(product([s_init], T))
    end_occ = list(product([s_end], T))                

    searcher_reachable = {}
    for t in T:
        for s in S_expand:
            if s == s_init or s == s_end: # init and occ cell are occupiable
                searcher_reachable[s, t] = 1
            elif (s, t) not in occ_cell:
                searcher_reachable[s, t] = 0
            else:
                searcher_reachable[s, t] = 1    
                
    """ creating set T that has no detection across all states
    """
    
    T_nd = []
    for t in T:
        if sum(q[(s, c), t] * searcher_reachable[s, t] for s, c in S_C if c == 0) == 0:
            T_nd.append(t)
    
    """ creating set (s, t) that has no detection
    """
    T_d = set(T) - set(T_nd)
    
    V_nd = {}
    for t in T_d:
        for s in S:
            if q[(s, 0), t] * searcher_reachable[s, t] == 0:
                if t not in V_nd.keys():
                    V_nd[t] = [(s, t)]
                else:
                    V_nd[t].append((s, t))
    
    for t in T_d:
        if t not in V_nd.keys():
            V_nd[t] = []
    
    """ setting parameter for the outer optimization problem """
    delta = 1e-4
    Xi_lb = 0
    Xi_ub = 1
    counter = 1
    
    
    """
    Defining model parameters
    """

    model_name = 'sp1_lm'
    m = gp.Model(model_name)
    m.setParam(GRB.Param.TimeLimit, 15 * 60)
    m.setParam(GRB.Param.Threads, 1)
    m.setParam(GRB.Param.LogFile, model_name)
    
    """ create index set """
    sub_X = list(product(L, S_expand, S_expand, T0))
    sub_Z = []
    for t in T_d:
        for s in S_expand:
            for l in L:
                if (s, t) not in V_nd[t]:
                    sub_Z = sub_Z + [(l, s, t)]
    sub_O = list(product(L, T))
    sub_recov_Z = list(product(L, S_expand, T))
    
    """ create variable """
    X = m.addVars(sub_X, vtype = GRB.CONTINUOUS, lb = 0, name = 'X')
    Z = m.addVars(sub_Z, vtype = GRB.INTEGER, lb = 0, name = 'Z')
    O = m.addVars(sub_O, vtype = GRB.CONTINUOUS, lb = 0, name = 'O')
    # O = m.addVars(sub_O, lb = 0, name = 'O')    
    for l in L:
        for t in T:
            O[l, t].ub = J_L[l]
    Xi = m.addVar(vtype = GRB.CONTINUOUS, name = 'Xi')
    
    
    """ set a Z_param to track value of decision variable Z """
    Z_param = {}
    for sub in sub_Z:
        Z_param[sub] = 0
    
    """ set a Z_recov_param to track Z value by recalculating X for every l, s, t """
    Z_recov_param = {} # Z_param is for recover the Z value post optimization
    for sub in sub_recov_Z:
        Z_recov_param[sub] = 0
    
    """
    Defining objective function
    """
    m.setObjective(Xi, GRB.MINIMIZE)    
    # m.setObjective(3, GRB.MINIMIZE)    
  
    """
    Defining constraints
    """    

    # 2.4b

    m.addConstrs((sum(X[l, s_prime, s, t - 1] for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
                  == sum(X[l, s, s_prime, t] for s_prime in S_expand if is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))  for l in L for s in S_expand for t in T), name = '24b') #2d
    
    # 2.4c
    """
    """
    m.addConstrs((sum(X[l, s, s_prime, 0] for s_prime in S_expand if is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
                  == x[l, s, 0] for l in L for s in S_expand), name = '24c') #2d
    
    # 2.4e
    """ 
    this s can only be at the searcher's init position 
    
    """
    m.addConstrs((sum(X[l, s, s_prime, t - 1] for s in [s_init] for s_prime in S_expand if s_prime != s and is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
                  == O[l, t] for l in L for t in T), name = '24e') #2d

    # 2.4f
    """ should we consider those move off the map next time period """
    m.addConstrs((sum(X[l, s, s_prime, t] for s in S for s_prime in S_expand if is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end)) 
                  <= sum(O[l, t_prime] for t_prime in T if t_prime <= t and t_prime >= t - tau[l] + 1)  for l in L for t in T), name = '24f') #2d

    # 2.5c
    m.addConstrs((sum(Z[l, s, t] for l in L )
                  <= n[s, t] for t in T_d for s in S_expand if (s, t) not in V_nd[t]), name = '25c') #2d
    
    # 2.5d
    for l in L:
        for s_expand in S_expand:
            for s_expand_prime in S_expand:
                for t in T:
                    X[l, s_expand, s_expand_prime, t].ub = min(J_L[l], n[s, t])
                    X[l, s_expand, s_expand_prime, t].lb = 0
                    
                    
# =============================================================================
#     # 2.5b
#     m.addConstrs((sum(X[l, s_prime, s, t - 1] for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
#                   == Z[l, s, t] for l in L for s in S_expand for t in T), name = '25b') #2d
# =============================================================================
           
# =============================================================================
#     # 2.5e
#     for l in L:
#         for s_expand in S_expand:
#             for t in T:
#                 Z[l, s_expand, t].ub = min(J_L[l], n[s, t])
#                 Z[l, s_expand, t].lb = 0
# =============================================================================
    
    # 4.14c
    m.addConstrs((sum(X[l, s_prime, s, t - 1] for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
                  == Z[l, s, t] for l, s, t in sub_Z), name = '414c') #2d
    
    # 4.14d
    m.addConstrs((sum(O[l, t] for t in T) == J_L[l] for l in L), name = '414d') #2d
        
    # 4.14e
    for l, s, t in sub_Z:
        Z[l, s, t].ub = min(J_L[l], n[s, t])
    
    start_time = time.time()
    
    # while Xi_ub - Xi_lb > delta * Xi_lb and counter <= 100:
    while Xi_ub - Xi_lb > delta * Xi_lb and time.time() - start_time <= 900:
    # while counter <= 5 and Xi_ub - Xi_lb > delta * Xi_lb and time.time() - start_time <= 900:
        
        """ for t in T_nd, print all states with non-zero searchers """
        for t in T_nd:
            for s in S_expand:
                for l in L:
                    if Z_recov_param[l, s, t] != 0:
                        print(l, s, t, Z_recov_param[l, s, t])
                        # print('q[(s, 0), t] is', q[(s, 0), t], 'q[(s, 1), t] is', q[(s, 1), t])
                        # print('searcher reachable', searcher_reachable[s, t])
                        
        """ for t in T_d but in V_nd[t], print all states with non-zero searchers """
        for t in T_d:
            for s in S_expand:
                if (s, t) in V_nd[t]:
                    for l in L:
                        if Z_recov_param[l, s, t] != 0:
                            print(l, s, t, Z_recov_param[l, s, t])
        
        ################ step 1 ################
        print('=============', counter, '===============')
        
        adj_Z_recov_param = Z_recov_param.copy()
        for l, s_temp, t in sub_recov_Z:
            if t in T_nd:
                adj_Z_recov_param[l, s_temp, t] = 0
            elif (s, t) in V_nd[t]:
                adj_Z_recov_param[l, s_temp, t] = 0
                

                
# =============================================================================
#         r_adj = {}
#         for t in T:
#             for s_c in S_C:
#                 if t == 1:
#                     r_adj[s_c, t] = p[s_c]
#                 else:
#                     r_adj[s_c, t] = sum([r_adj[s_c_prime, t - 1] * np.exp(-sum(alpha[l, s_c[1]] * adj_Z_recov_param[l, s_c_prime[0], t - 1] for l in L)) * gamma[s_c_prime, s_c, t - 1] for s_c_prime in S_C if is_backward_state_with_cam(s_c_prime, s_c)])
#                     # r[s_c, t] = sum([r[s_c_prime, t - 1] * np.exp(-sum(alpha[l, s_c[1]] * Z_param[l, s_c_prime[0], t - 1] for l in L)) * gamma[s_c_prime, s_c, t - 1] for s_c_prime in S_C])
#         
#         s_adj = {}
#         for t in T[::-1]:
#             for s_c in S_C:
#                 if t == ending_time:
#                     s_adj[s_c, t] = 1
#                 else:
#                     s_adj[s_c, t] = sum([s_adj[s_c_prime, t + 1] * np.exp(-sum(alpha[l, s_c[1]] * adj_Z_recov_param[l, s_c_prime[0], t + 1] for l in L)) * gamma[s_c, s_c_prime, t] for s_c_prime in S_C if is_backward_state_with_cam(s_c, s_c_prime)]) # if s_c_prime is s_c's forward
#                     # s[s_c, t] = sum([s[s_c_prime, t + 1] * np.exp(-sum(alpha[l, s_c[1]] * Z_param[l, s_c_prime[0], t + 1] for l in L)) * gamma[s_c, s_c_prime, t] for s_c_prime in S_C]) # if s_c_prime is s_c's forward
# 
# =============================================================================
        
        r = {}
        for t in T:
            for s_c in S_C:
                if t == 1:
                    r[s_c, t] = p[s_c]
                else:
                    r[s_c, t] = sum([r[s_c_prime, t - 1] * np.exp(-sum(alpha[l, s_c[1]] * Z_recov_param[l, s_c_prime[0], t - 1] for l in L)) * gamma[s_c_prime, s_c, t - 1] for s_c_prime in S_C if is_backward_state_with_cam(s_c_prime, s_c)])
                    # r[s_c, t] = sum([r[s_c_prime, t - 1] * np.exp(-sum(alpha[l, s_c[1]] * Z_param[l, s_c_prime[0], t - 1] for l in L)) * gamma[s_c_prime, s_c, t - 1] for s_c_prime in S_C])

        s = {}
        for t in T[::-1]:
            for s_c in S_C:
                if t == ending_time:
                    s[s_c, t] = 1
                else:
                    s[s_c, t] = sum([s[s_c_prime, t + 1] * np.exp(-sum(alpha[l, s_c[1]] * Z_recov_param[l, s_c_prime[0], t + 1] for l in L)) * gamma[s_c, s_c_prime, t] for s_c_prime in S_C if is_backward_state_with_cam(s_c, s_c_prime)]) # if s_c_prime is s_c's forward
                    # s[s_c, t] = sum([s[s_c_prime, t + 1] * np.exp(-sum(alpha[l, s_c[1]] * Z_param[l, s_c_prime[0], t + 1] for l in L)) * gamma[s_c, s_c_prime, t] for s_c_prime in S_C]) # if s_c_prime is s_c's forward

        # f_Z = sum([r[c, 1] * np.exp(-alpha * Z_param[c, 1]) * s[c, 1] for c in C])
        # f_Z_2 = sum([r[c, ending_time] * np.exp(-alpha * Z_param[c, ending_time]) * s[c, ending_time] for c in C])
        
        
        
        
        test_time = 5
        f_Z = sum([r[s_c, test_time] * np.exp(-sum(alpha[l, s_c[1]] * Z_recov_param[l, s_c[0], test_time] for l in L)) * s[s_c, test_time] for s_c in S_C])
        # f_Z_reduced = sum([r[s_c, test_time] * np.exp(-sum(alpha[l, s_c[1]] * adj_Z_recov_param[l, s_c[0], test_time] for l in L)) * s[s_c, test_time] for s_c in S_C])
        
        for l, s_temp, t in sub_recov_Z:
            if adj_Z_recov_param[l, s_temp, t] != Z_recov_param[l, s_temp, t]:
                print('stop', l, s_temp, t)
                print('adj', adj_Z_recov_param[l, s_temp, t])
                print('org', Z_recov_param[l, s_temp, t])
                if s_temp != s_init and s_temp != s_end:
                    print('r, s is', r[(s_temp, 0), t], s[(s_temp, 0), t])
                # print('r, s is', r_adj[(s_temp, 0), t], s_adj[(s_temp, 0), t])
        
# =============================================================================
#         if f_Z_reduced != f_Z:
#             print('f_Z_reduced', f_Z_reduced)
#             print('f_Z', f_Z)
#             break
# =============================================================================


        print('f(Z) equals to', f_Z)    
        if f_Z < Xi_ub:
            Xi_ub = f_Z   
        if Xi_ub - Xi_lb <= delta * Xi_lb:
            break
        
        print('Before optimization')
        print('Xi upper', Xi_ub)
        print('Xi lower', Xi_lb)
        

        ################ step 2 ##############
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
        
        """
        use 1 - s_c[1] to adjust the finite difference: when the searcher is in the camouflage mode, the finite difference will be 0
        """
        # m.addConstr(f_Z + sum([(1 - s_c[1]) * r[s_c, t] * (np.exp(-alpha[l, s_c[1]] * (Z_param[l, s_c[0], t] + 1)) - np.exp(-alpha[l, s_c[1]] * (Z_param[l, s_c[0], t]))) * s[s_c, t] * (Z[l, s_c[0], t] - Z_param[l, s_c[0], t]) for l in L for s_c in S_C for t in T]) <= Xi, name = 'cut_' + str(counter))
        # m.addConstr(f_Z + sum([searcher_reachable[s_c[0], t] * r[s_c, t] * (np.exp(-alpha[l, s_c[1]] * (Z_param[l, s_c[0], t] + 1)) - np.exp(-alpha[l, s_c[1]] * (Z_param[l, s_c[0], t]))) * s[s_c, t] * (Z[l, s_c[0], t] - Z_param[l, s_c[0], t]) for l in L for s_c in S_C for t in T_d if (s_c[0], t) not in V_nd[t]]) <= Xi, name = 'cut_' + str(counter))
        # m.addConstr(f_Z + sum([r[s_c, t] * (np.exp(-alpha[l, s_c[1]] * (Z_param[l, s_c[0], t] + 1)) - np.exp(-alpha[l, s_c[1]] * (Z_param[l, s_c[0], t]))) * s[s_c, t] * (Z[l, s_c[0], t] - Z_param[l, s_c[0], t]) for l in L for s_c in S_C for t in T_d if (s_c[0], t) not in V_nd[t]]) <= Xi, name = 'cut_' + str(counter))
        m.addConstr(f_Z + sum([r[s_c, t] * (np.exp(-alpha[l, s_c[1]] * (Z_param[l, s_c[0], t] + 1)) - np.exp(-alpha[l, s_c[1]] * (Z_param[l, s_c[0], t]))) * s[s_c, t] * (Z[l, s_c[0], t] - Z_param[l, s_c[0], t]) for l in L for s_c in S_C for t in T_d if (s_c[0], t) not in V_nd[t] and s_c[1] == 0]) <= Xi, name = 'cut_' + str(counter))

        
        """ Solving
        """

        m.optimize()
        
        Xi_iter_lb = m.poolObjBound
        Xi_iter_ub = m.objVal

        print('Pi result after optimization')
        print('Xi iter upper', Xi_iter_ub)
        print('Xi iter lower', Xi_iter_lb)

        
        if Xi_iter_lb > Xi_lb:
            Xi_lb = Xi_iter_lb
        
        """ updating Z_param that tracking the value of Z """
        for sub in sub_Z:
            Z_param[sub] = Z[sub].X
        
        print('====== recalculated Z_ct ======')
        # Z_recov_prior = Z_recov_param.copy()
        for sub in sub_recov_Z: # for the grouping algorithm, we need to recalculate Z)
            l, s, t = sub
            Z_recov_param[l, s, t] = sum(X[l, s_prime, s, t - 1].X for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
            
        print('After optimization')
        print('Xi upper', Xi_ub)
        print('Xi lower', Xi_lb)
        
        counter += 1
        print('sub problem MIP gap is', mip_gap)

    print('********** alpha is **********', alpha[1, 0], alpha[2, 0])
    print('********** T is', ending_time)
    print('********** Duration is', tau)
    print('********** J_l is', J_L)
    print('********** Grid size is', grid_size)    
    gap = (Xi_ub - Xi_lb) / Xi_lb
    print('Final MIPGap is', gap)
    end_time = time.time()
    running_time = end_time - start_time
    print("Running time is", running_time)
    # time_log[grid_size] = [gap, running_time, Xi_ub]
# print(time_log)
# with open('time_log_T10.txt', 'w') as log_result:
#    log_result.write(json.dumps(time_log)) 
    print('********* optimal solution for X ***********')
    sub_X = sorted(sub_X, key = lambda x: (x[0], x[3]))    
    for sub in sub_X:
        if X[sub].X != 0:
            print(sub, X[sub].X)
            
    print('********* optimal solution for O ***********')
    sub_O = sorted(sub_O, key = lambda x: (x[0], x[1]))    
    for sub in sub_O:
        if O[sub].X != 0:
            print(sub, O[sub].X)

    

