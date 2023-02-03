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
Helper function
"""

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

grid_size = 9
ending_time_grid = range(7, 8)
# num_scenario = 1000



for ending_time in ending_time_grid:
    ending_time = ending_time
    print('=============================')
    print('Ending time is', ending_time)
    print('=============================')

    """
    Creating set
    """
    
    C = [(i, j) for i in range(1, grid_size + 1) for j in range(1, grid_size + 1)]
    camouflage_state = [0, 1]
    S = list(product(C, camouflage_state)) # state contains (c, camouflage_state)
    S_searcher = [s for s in S if s[1] == 0]
    
    
    T = list(range(1, ending_time + 1))
    T0 = [0] + T
    # Omega = list(range(1, num_scenario + 1))
    # J = 3
    # J_set = list(range(1, J + 1))
    L = [1, 2] # set of searchers' type
    n_L = {1: 2, 2: 1} # number of searchers for each searcher type
    alpha_l = {1: 0.1, 2: 0.2} # detection rate for each searcher type
    
    J_l = {}
    for l in L:
        J_l[l] = list(range(1, n_L[l] + 1))
    
    
    # I = list(range(0, J * ending_time + 1))
    


    
    
    """ Import data
    """
    if platform.system() == 'Windows':
        data_folder = 'E:\\Research\\Route_Optimization_for_Multiple_Searchers\\Python\\'
    else:
        data_folder = os.path.dirname(os.path.realpath(__file__))
    
    
    
# =============================================================================
#     zeta_raw = pd.read_csv(data_folder + '/Zeta.csv', header = None, index_col = 0)
#     q_raw = pd.read_csv(data_folder + '/q.csv', header = 0, index_col = 0)
#     q_raw.columns = [int(col) for col in q_raw.columns]
#     sub_q = list(product(C, T))
#     q = {}
#     for sub in sub_q:
#         c_two_dim, t = sub
#         c_one_dim = (c_two_dim[0] - 1) * grid_size + c_two_dim[1]
#         q[c_two_dim, t] = q_raw.loc[c_one_dim, t]
# =============================================================================
    # zeta_raw = pd.read_csv(data_folder + '\Zeta.csv', header = None, index_col = 0)
        

                            
# =============================================================================
#     Zeta = {}
#     for path in range(1, zeta_raw.shape[0] + 1):
#         for t in range(1, ending_time + 1):
#             all_cells = C
#             for cell in all_cells:
#                 Zeta[(cell, t, path)] = 0 # set Zeta equal to 0
#             cell_one_dim = zeta_raw.loc[path, 3 * (t - 1) + 1] # extract the occupied cell loc from Zeyu's path
#             cell_two_dim = (cell_one_dim // grid_size + 1, np.mod(cell_one_dim, grid_size))
#             Zeta[(cell_two_dim, t, path)] = 1 # set Zeta equal to 1 for occupied cell
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
    
    # =============================================================================
    # q = pd.read_csv('')
    # q = np.random.uniform(low = 0, high = 1, size = num_scenario)
    # q = q / sum(q) # normalize to a probablity distribution summing up to 1
    # q = dict(zip(Omega, q))
    # =============================================================================
    
    
    p = {} # prob target is in s at t = 1
    # The target could exist at any cell in row 1 with equal probability
    # The target cannot be in a camouflage state at t = 1?
    for s in S:
        c, cam_s = s[0], s[1]
        prob = 1 / grid_size
        if c[0] == 1 and cam_s == 0: # row 1 and not in camouflage mode
            p[s] = prob
        else:
            p[s] = 0
    
    # define detection rate by l, s, t
    alpha_sub = list(product(L, S, T))
    alpha = {} # detection rate
    for ind in alpha_sub:
        if ind[1][1] == 1: # target is in camouflage mode
            alpha[ind] = 0
        else:
            alpha[ind] = alpha_l[ind[0]] # use the detection rate from type l
    
# =============================================================================
#     alpha_value = [-3 * np.log(0.4) / J] * grid_size * grid_size * ending_time
#     alpha = dict(zip(alpha_sub, alpha_value))
# =============================================================================
    
    sub_gamma = list(product(S, S, T))
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
    start_state = ((1, 1), 0)
    for s_prime in S:
        if start_state[1] == 1:
            if gamma[(start_state, s_prime, 1)] != 0:
                print(start_state, s_prime, 1, 'gamma value', gamma[(start_state, s_prime, 1)])
        else:
            s_prime_c, s_prime_cam = s_prime[0], s_prime[1]
            s_c, s_cam = start_state[0], start_state[1]
            if is_nearby_cell(s_c, s_prime_c) and s_prime[1] == 0:
                print(start_state, s_prime, 1, 'gamma value', gamma[(start_state, s_prime, 1)])
            elif s_c == s_prime_c and s_prime[1] == 1:
                print(start_state, s_prime, 1, 'gamma value', gamma[(start_state, s_prime, 1)])
    
    
    sub_q = list(product(S, T)) 
    sub_q = sorted(sub_q, key = lambda x: x[1])
    q = {} # create param values for q
    for sub in sub_q:
        s, t = sub
        # print(sub)
        if t == 1:
            q[s, t] = p[s]
        else:
            # c_one_dim = (c_two_dim[0] - 1) * grid_size + c_two_dim[1]
            c, cam = s[0], s[1]
            if cam == 1: # if target is in camouflage at t, it either in the same state at t - 1 or change into camouflage from the same cell
                q[s, t] = q[(c, 1), t - 1] * stay_cam + q[(c, 0), t - 1] * to_cam
            else: # if target is not in camouflage at t, it either moves from a non-camouflage cell or changes back a camouflage state at the same cell
                q[s, t] = q[(c, 1), t - 1] * from_cam + sum([q[(c_prime, 0), t - 1] * gamma[(c_prime, 0), (c, 0), t - 1] for c_prime in C if is_nearby_cell(c, c_prime)])
    
    """
    The following code is to check if q is calculated correctly
    """
    test_state = (((1, 1), 0), 1)
    test_state = (((1, 1), 1), 1)
    test_state = (((2, 1), 0), 1)
    test_state = (((2, 1), 0), 2)
    test_state = (((2, 1), 1), 2)
    q[test_state]
    t_test = 2
    for s in S:
        print(s, t_test, q[(s, t_test)])
    print(sum([q[(s, t_test)] for s in S]))
    
    # =============================================================================
    # xx = {} 
    # for c in C:
    #     for t in T0:
    #         if c == (1, 1) and t == 0:
    #             xx[c, t] = J
    #         else:
    #             xx[c, t] = 0
    #             
    # =============================================================================
    searcher_init_loc = (0, 0)
    searcher_init_state = (searcher_init_loc, 0)
    
    xx = {} # number of searchers at t = 0
    for l in L:
        for s in S_searcher:
            if s == searcher_init_state:
                xx[l, s] = n_L[l] # set the number equal to the number of searchers of type l
            else:
                xx[l, s] = 0
    
# =============================================================================
#     xx = {} 
#     for c in C:
#         if c == (1, 1):
#             xx[c] = J
#         else:
#             xx[c] = 0
# =============================================================================
    
# =============================================================================
#     w_param = {}
#     for c in C:
#         for t in T:
#             if q[c, t] > 0:
#                 w_param[c, t] = 1
#             else:
#                 w_param[c, t] = 0
# =============================================================================
    
    w_param = {}
    for s in S: # this s is for the target
        for t in T:
            if q[s, t] > 0:
                w_param[s, t] = 1
            else:
                w_param[s, t] = 1 # do not reduce any pair for now


    """
    Defining decision variables
    """
    sub_X = list(product(L, S_searcher, S_searcher, T0))
    sub_ct = list(product(S, T)) # P and W is indexed by sub_ct
    
    sub_ctj = []
    for l in L:
        sub_ctj = sub_ctj + list(product([l], S, T, range(1, n_L[l] + 1))) # Q is indexed by sub_ctj
    
    
    
    sub_zset = list(product(S, T))
    sub_zset = [each for each in sub_zset if w_param[each] == 1]
    
    sub_zset_j = []
    for l in L:
        sub_zset_j = sub_zset_j + list(product([l], S_searcher, T, J_l[l])) # used by V
    
    
    X = m.addVars(sub_X, lb = 0, name = 'X') # related to searcher; should be indexed by S_searcher
    # Z = m.addVars(sub_ct, lb = 0, ub = J, vtype = GRB.INTEGER, name = 'Z')
    # U = m.addVars(Omega, lb = 0, name = 'U')
    P = m.addVars(sub_ct, lb = 0, name = 'P') # prob that the target is in states in period t and was not detected prior to t
    Q = m.addVars(sub_ctj, lb = 0, name = 'Q') # auxiliary variable related to P and alpha
    # V = m.addVars(sub_ctj, vtype = GRB.BINARY, name = 'V')
    V = m.addVars(sub_zset_j, vtype = GRB.BINARY, name = 'V') # = 1 if J_l searchers of class l in state s in period t
    W = m.addVars(sub_ct, lb = 0, name = 'W') # auxiliary variable related to P and alpha
    
    """
    Defining objective function
    """
    m.setObjective(1 - sum(Q[c, t, j] for t in T for c in C for j in J_set), GRB.MINIMIZE)    
    
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
        
    # m.addConstrs((Q[c, t, j] <= q[c, t] * (1 - np.exp(-j * alpha[c, t])) * V[c, t, j] for c in C for t in T for j in J_set), name = '23')
    m.addConstrs((Q[reduce_ct[0], reduce_ct[1], j] <= q[reduce_ct[0], reduce_ct[1]] * (1 - np.exp(-j * alpha[reduce_ct[0], reduce_ct[1]])) * V[(reduce_ct[0], reduce_ct[1]), j] for reduce_ct in sub_zset for j in J_set), name = '23')

    
    m.addConstrs((Q[c, t, j] <= (1 - np.exp(-j * alpha[c, t])) * P[c, t] for c in C for t in T for j in J_set), name = '24')
    m.addConstrs((P[c, t + 1] == sum(gamma[c_prime, c, t] * W[c_prime, t] for c_prime in C) for c in C for t in T[:-1]), name = '25')
    m.addConstrs((W[c, t] <= P[c, t] for c in C for t in T), name = '26')
    # m.addConstrs((W[c, t] <= np.exp(-j * alpha[c, t]) * P[c, t] + q[c, t] * (1 - np.exp(-j * alpha[c, t])) * (1 - V[c, t, j]) for c in C for t in T for j in J_set), name = '27')
    m.addConstrs((W[reduce_ct[0], reduce_ct[1]] <= np.exp(-j * alpha[reduce_ct[0], reduce_ct[1]]) * P[reduce_ct[0], reduce_ct[1]] + q[reduce_ct[0], reduce_ct[1]] * (1 - np.exp(-j * alpha[reduce_ct[0], reduce_ct[1]])) * (1 - V[(reduce_ct[0], reduce_ct[1]), j]) for reduce_ct in sub_zset for j in J_set), name = '27')
    
    m.addConstrs((P[c, 1] == p[c] for c in C), name = '28')
    m.addConstrs((P[c, t] <= q[c, t] for c in C for t in T), name = '29')
    m.addConstrs((sum(X[c_prime, reduce_ct[0], reduce_ct[1] - 1] for c_prime in C if is_nearby_cell(reduce_ct[0], c_prime)) == sum(j * V[(reduce_ct[0], reduce_ct[1]), j] for j in J_set) for reduce_ct in sub_zset), name = '30')
    m.addConstrs((sum(V[(reduce_ct[0], reduce_ct[1]), j] for j in J_set) <= 1 for reduce_ct in sub_zset), name = '31')
    
    
    m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == sum(X[c, c_prime, t] for c_prime in C if is_nearby_cell(c, c_prime))  for c in C for t in T), name = '14') #2d
    m.addConstrs((sum(X[c, c_prime, 0] for c_prime in C if is_nearby_cell(c, c_prime)) == xx[c] for c in C), name = '15') #2d
    
    """ Solving
    """
    m.optimize()
