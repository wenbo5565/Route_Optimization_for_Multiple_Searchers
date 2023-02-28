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
ending_time_grid = [10, 12, 14, 15, 16, 17, 18, 20]
# ending_time = 15
# num_scenario = 1000
J = 15
J_2 = int(J * 0.7)
J_1 = J - J_2

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
    """
    model_name = 'sp1_lm'
    m = gp.Model(model_name)
    m.setParam(GRB.Param.TimeLimit, 15 * 60)
    m.setParam(GRB.Param.Threads, 1)
    m.setParam(GRB.Param.LogFile, model_name)
    # m.setParam(GRB.Param.MIPGap, 1e-5)
    
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
    
    
    """
    Defining decision variables
    """

    sub_X = list(product(L, S_expand, S_expand, T0))
    sub_P = list(product(S_C, T)) # corresponds to the target
    sub_V = []
    for l in L:
        sub_V = sub_V + list(product([l], S_expand, T, J_l[l]))
    sub_Q = []
    for l in L:
        sub_Q = sub_Q + list(product([l], S_C, T, J_l[l]))
    sub_W = list(product(S_C, T))
# =============================================================================
#     for l in L:
#         sub_W = sub_W + list(product([l], S_C, T, J_l[l]))
# =============================================================================
    sub_O = list(product(L, T))
    
    X = m.addVars(sub_X, lb = 0, name = 'X')
    # Z = m.addVars(sub_ct, lb = 0, ub = J, vtype = GRB.INTEGER, name = 'Z')
    # U = m.addVars(Omega, lb = 0, name = 'U')
    P = m.addVars(sub_P, lb = 0, name = 'P')
    Q = m.addVars(sub_Q, lb = 0, name = 'Q')
    V = m.addVars(sub_V, vtype = GRB.BINARY, name = 'V')
    W = m.addVars(sub_W, lb = 0, name = 'W')
    O = m.addVars(sub_O, vtype = GRB.INTEGER, lb = 0, name = 'O')
    for l in L:
        for t in T:
            O[l, t].ub = n_L[l]
    
    """
    Defining objective function
    """
    m.setObjective(1 - sum(Q[l, s_c, t, j] for l, s_c, t, j in sub_Q), GRB.MINIMIZE)    
    
    """
    Defining constraints
    """
    
    
    # =============================================================================
    # def is_forward_cell(c, c_next):
    #     if c == c_next:
    #         return True
    #     elif c[0] == c_next[0]:
    #         if c[1] + 1 == c_next[1] or c[1] - 1 == c_next[1]:
    #             return True
    #     elif c[1] == c_next[1]:
    #         if c[0] + 1 == c_next[0] or c[0] - 1 == c_next[0]:
    #             return True
    #     else:
    #         return False
    #     
    # def is_reverse_cell(c, c_last):
    #     if c == c_last:
    #         return True
    #     elif c[0] == c_last[0]:
    #         if c[1] + 1 == c_last[1] or c[1] - 1 == c_last[1]:
    #             return True
    #     elif c[1] == c_last[1]:
    #         if c[0] + 1 == c_last[0] or c[0] - 1 == c_last[0]:
    #             return True
    #     else:
    #         return False
    # =============================================================================
    """ check this constraints """
    # m.addConstrs((Q[l, s_c, t, j] <= q[s_c, t] * (1 - np.exp(sum(-j_l * alpha[l, s_c[1]] for j_l in J_ind[l]))) * V[l, s_c[0], t, j] for l, s_c, t, j in sub_Q), name = '23')
    
    m.addConstrs((Q[l, s_c, t, j] <= q[s_c, t] * (1 - np.exp(-j * alpha[l, s_c[1]])) * V[l, s_c[0], t, j] for l in L for s_c in S_C for t in T for j in J_l[l]), name = '23')

    
    m.addConstrs((Q[l, s_c, t, j] <= (1 - np.exp(-j * alpha[l, s_c[1]])) * P[s_c, t] for l, s_c, t, j in sub_Q), name = '24')
    
    m.addConstrs((P[s_c, t + 1] == sum(gamma[s_c_prime, s_c, t] * W[s_c_prime, t] for s_c_prime in S_C if is_backward_state_with_cam(s_c_prime, s_c)) for s_c in S_C for t in T[:-1]), name = '25')
    
    m.addConstrs((W[s_c, t] <= P[s_c, t] for s_c in S_C for t in T), name = '26')
    
    # """ check this constraint """
    m.addConstrs((W[s_c, t] <= np.exp(-j * alpha[l, s_c[1]]) * P[s_c, t] + q[s_c, t] * (1 - np.exp(-j * alpha[l, s_c[1]])) * (1 - V[l, s_c[0], t, j]) for l in L for s_c in S_C for t in T for j in J_ind[l]), name = '27')
    
    m.addConstrs((P[s_c, 1] == p[s_c] for s_c in S_C), name = '28')
    
    m.addConstrs((P[s_c, t] <= q[s_c, t] for s_c in S_C for t in T), name = '29')
    
    m.addConstrs((sum(X[l, s_expand_prime, s_expand, t - 1] for s_expand_prime in S_expand if is_backward_cell(s_expand_prime, s_expand, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end)) 
                  == sum(j * V[l, s_expand, t, j] for j in J_ind[l]) for l in L for s_expand in S_expand for t in T), name = '30')
    
    m.addConstrs((sum(V[l, s_expand, t, j] for j in J_ind[l]) <= 1 for l in L for s_expand in S_expand for t in T), name = '31')
    
    """ adding 2.4b to 2.4f     
    """
    # 2.4b

    m.addConstrs((sum(X[l, s_prime, s, t - 1] for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
                  == sum(X[l, s, s_prime, t] for s_prime in S_expand if is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))  for l in L for s in S_expand for t in T), name = 'continum') #2d
    
    # m.addConstrs((sum(X[l, s_prime, s, t - 1] for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end)) == sum(X[l, s, s_prime, t] for s_prime in S if is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end))  for l in L for s in (S + s_init) for t in T), name = 'continum') #2d

    # m.addConstrs((sum(X[l, s_prime, s, t - 1] for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end)) == sum(X[l, s, s_prime, t] for s_prime in S if is_forward_cell(s, s_prime, starting_c = s_init, ending_c = s_end))  for l in L for s in S for t in T), name = 'continum') #2d


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

    
    
# =============================================================================
#     m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == sum(X[c, c_prime, t] for c_prime in C if is_nearby_cell(c, c_prime))  for c in C for t in T), name = '14') #2d
#     
#     m.addConstrs((sum(X[c, c_prime, 0] for c_prime in C if is_nearby_cell(c, c_prime)) == xx[c] for c in C), name = '15') #2d
# =============================================================================
    
    m.optimize()

    """ Solving
    """
    #
    print("********** optimal solution for V **********")
    sub_V = sorted(sub_V, key = lambda x: (x[0], x[2]))    
    for sub in sub_V:
        if V[sub].X != 0:
            print(sub, V[sub].X)
            
    print('********** optimal solution for O **********')
    sub_O = sorted(sub_O, key = lambda x: x[1])
    for sub in sub_O:
        if O[sub].X != 0:
            print(sub, O[sub].X)
            
    print("********** optimal solution for Q **********")
    sub_Q = sorted(sub_Q, key = lambda x: (x[0], x[2]))    
    for sub in sub_Q:
        if Q[sub].X != 0:
            print(sub, Q[sub].X)
            
# J = 15