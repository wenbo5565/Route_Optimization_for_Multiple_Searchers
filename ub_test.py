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

grid_size = 9
# ending_time_grid = range(7, 16)
ending_time_grid = [15]
# num_scenario = 1000


for ending_time in ending_time_grid:
    ending_time = ending_time
    # ending_time = ending_time
    print('=============================')
    print('Ending time is', ending_time)
    print('=============================')

    """
    Creating set
    """
    
    C = [(i, j) for i in range(1, grid_size + 1) for j in range(1, grid_size + 1)]
    T = list(range(1, ending_time + 1))
    T0 = [0] + T
    # Omega = list(range(1, num_scenario + 1))
    J = 15
    J_set = list(range(1, J + 1))
    I = list(range(0, J * ending_time + 1))
    
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
#     
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
    
    
    p = {}
    for c in C:
        if c == ((grid_size + 1) / 2, (grid_size + 1) / 2):
            p[c] = 1
        else:
            p[c] = 0
    
    alpha_sub = list(product(C, T))
    alpha_value = [-3 * np.log(0.4) / J] * grid_size * grid_size * ending_time
    alpha = dict(zip(alpha_sub, alpha_value))
    
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
    
    xx = {} 
    for c in C:
        if c == (1, 1):
            xx[c] = J
        else:
            xx[c] = 0
    
    
    
    """
    Defining decision variables
    """
    sub_X = list(product(C, C, T0))
    sub_ct = list(product(C, T))
    sub_ctj = list(product(C, T, J_set))
    
    X = m.addVars(sub_X, lb = 0, name = 'X')
    # Z = m.addVars(sub_ct, lb = 0, ub = J, vtype = GRB.INTEGER, name = 'Z')
    # U = m.addVars(Omega, lb = 0, name = 'U')
    P = m.addVars(sub_ct, lb = 0, name = 'P')
    Q = m.addVars(sub_ctj, lb = 0, name = 'Q')
    V = m.addVars(sub_ctj, vtype = GRB.BINARY, name = 'V')
    W = m.addVars(sub_ct, lb = 0, name = 'W')
    
    # sol = [[[1, 2], 1, 1.0], [[2, 1], 1, 2.0], [[2, 2], 2, 1.0], [[3, 1], 2, 2.0], [[2, 3], 3, 1.0], [[3, 2], 3, 2.0], [[2, 4], 4, 1.0], [[4, 2], 4, 2.0], [[3, 4], 5, 1.0], [[4, 3], 5, 2.0], [[3, 5], 6, 1.0], [[4, 4], 6, 1.0], [[5, 3], 6, 1.0], [[4, 4], 7, 1.0], [[4, 5], 7, 1.0], [[5, 4], 7, 1.0], [[4, 5], 8, 1.0], [[4, 6], 8, 1.0], [[5, 5], 8, 1.0], [[5, 5], 9, 1.0], [[5, 6], 9, 1.0], [[6, 5], 9, 1.0], [[5, 7], 10, 0.55], [[6, 4], 10, 1.0], [[6, 5], 10, 1.0], [[6, 6], 10, 0.45], [[5, 4], 11, 1.0], [[5, 7], 11, 0.5], [[6, 4], 11, 0.15], [[6, 6], 11, 0.5], [[6, 7], 11, 0.05], [[7, 5], 11, 0.35], [[7, 6], 11, 0.45], [[4, 4], 12, 1.0], [[5, 6], 12, 1.0], [[6, 4], 12, 0.15], [[6, 6], 12, 0.5], [[7, 5], 12, 0.35], [[4, 5], 13, 1.0], [[5, 5], 13, 1.0], [[6, 5], 13, 1.0]]
    sol = [[[1, 2], 1, 13.0], [[2, 1], 1, 2.0], [[2, 2], 2, 13.0], [[3, 1], 2, 2.0], [[3, 2], 3, 13.0], [[4, 1], 3, 2.0], [[3, 3], 4, 13.0], [[5, 1], 4, 2.0], [[3, 4], 5, 7.0], [[4, 3], 5, 6.0], [[5, 2], 5, 2.0], [[3, 5], 6, 5.0], [[4, 4], 6, 6.0], [[5, 3], 6, 4.0], [[3, 5], 7, 1.0], [[4, 4], 7, 3.0], [[4, 5], 7, 5.0], [[5, 3], 7, 1.0], [[5, 4], 7, 5.0], [[4, 5], 8, 3.0], [[4, 6], 8, 3.0], [[5, 4], 8, 2.0], [[5, 5], 8, 4.0], [[6, 4], 8, 3.0], [[4, 5], 9, 1.0], [[4, 6], 9, 2.17], [[5, 4], 9, 1.0], [[5, 5], 9, 3.0], [[5, 6], 9, 4.0], [[6, 4], 9, 0.83], [[6, 5], 9, 3.0], [[4, 5], 10, 1.0], [[4, 7], 10, 0.34], [[5, 4], 10, 1.0], [[5, 5], 10, 1.0], [[5, 6], 10, 1.0], [[5, 7], 10, 1.76], [[6, 4], 10, 0.83], [[6, 5], 10, 2.0], [[6, 6], 10, 3.11], [[7, 4], 10, 0.83], [[7, 5], 10, 2.13], 
           [[4, 5], 11, 1.0], [[4, 6], 11, 1.0], [[4, 7], 11, 0.44], [[5, 4], 11, 1.0], [[5, 5], 11, 1.0], [[5, 6], 11, 2.0], [[5, 7], 11, 0.42], [[6, 4], 11, 1.55], [[6, 5], 11, 1.0], [[6, 6], 11, 1.56], [[6, 7], 11, 0.55], [[7, 4], 11, 0.79], [[7, 5], 11, 0.97], [[7, 6], 11, 1.71], [[4, 4], 12, 1.0], [[4, 5], 12, 2.0], [[4, 6], 12, 0.9], [[4, 7], 12, 1.06], [[5, 5], 12, 1.0], [[5, 6], 12, 1.0], [[5, 7], 12, 1.54], [[6, 3], 12, 0.52], [[6, 4], 12, 0.75], [[6, 5], 12, 2.0], [[6, 6], 12, 0.08], [[6, 7], 12, 1.47], [[7, 4], 12, 0.05], [[7, 5], 12, 1.55], [[7, 6], 12, 0.08], [[4, 4], 13, 0.3], [[4, 5], 13, 1.0], [[4, 6], 13, 0.48], [[5, 3], 13, 0.45], [[5, 4], 13, 2.0], [[5, 5], 13, 2.0], [[5, 6], 13, 2.0], [[5, 7], 13, 1.17], [[6, 4], 13, 1.05], [[6, 5], 13, 1.0], [[6, 6], 13, 1.0], [[6, 7], 13, 0.55], [[7, 5], 13, 0.67], [[7, 6], 13, 1.33]]
    
    sol_ctj = [[(ctj[0][0], ctj[0][1]), ctj[1], ctj[2]] for ctj in sol]
    sol_ct = [[(ctj[0][0], ctj[0][1]), ctj[1]] for ctj in sol]
    
    for t in T:
        for c in C:
            for j in J_set:
                if [c, t, j] in sol_ctj: 
                    V[c, t, j].ub = 1
                    V[c, t, j].lb = 1
                elif [c, t] in sol_ct:
                    continue
                else:
                    V[c, t, j].ub = 0
                    V[c, t, j].lb = 0
                
    
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
        
    m.addConstrs((Q[c, t, j] <= q[c, t] * (1 - np.exp(-j * alpha[c, t])) * V[c, t, j] for c in C for t in T for j in J_set), name = '23')
    m.addConstrs((Q[c, t, j] <= (1 - np.exp(-j * alpha[c, t])) * P[c, t] for c in C for t in T for j in J_set), name = '24')
    m.addConstrs((P[c, t + 1] == sum(gamma[c_prime, c, t] * W[c_prime, t] for c_prime in C) for c in C for t in T[:-1]), name = '25')
    m.addConstrs((W[c, t] <= P[c, t] for c in C for t in T), name = '26')
    m.addConstrs((W[c, t] <= np.exp(-j * alpha[c, t]) * P[c, t] + q[c, t] * (1 - np.exp(-j * alpha[c, t])) * (1 - V[c, t, j]) for c in C for t in T for j in J_set), name = '27')
    m.addConstrs((P[c, 1] == p[c] for c in C), name = '28')
    m.addConstrs((P[c, t] <= q[c, t] for c in C for t in T), name = '29')
    m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == sum(j * V[c, t, j] for j in J_set) for c in C for t in T), name = '30')
    m.addConstrs((sum(V[c, t, j] for j in J_set) <= 1 for c in C for t in T), name = '31')
    
    
    m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == sum(X[c, c_prime, t] for c_prime in C if is_nearby_cell(c, c_prime))  for c in C for t in T), name = '14') #2d
    m.addConstrs((sum(X[c, c_prime, 0] for c_prime in C if is_nearby_cell(c, c_prime)) == xx[c] for c in C), name = '15') #2d
    
    """ Solving
    """
    m.optimize()
    
    for t in T:
        for c in C:
            for j in J_set:
                if V[c, t, j].X != 0:
                    print(c, t, j, V[c, t, j].X)