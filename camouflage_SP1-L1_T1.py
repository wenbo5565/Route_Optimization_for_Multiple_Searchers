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
import ast

"""
Define helper functions
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

def is_forward_cell(c, c_next):
    if c == c_next:
        return True
    elif c[0] == c_next[0]:
        if c[1] + 1 == c_next[1] or c[1] - 1 == c_next[1]:
            return True
    elif c[1] == c_next[1]:
        if c[0] + 1 == c_next[0] or c[0] - 1 == c_next[0]:
            return True
    else:
        return False
    
def is_reverse_cell(c, c_last):
    if c == c_last:
        return True
    elif c[0] == c_last[0]:
        if c[1] + 1 == c_last[1] or c[1] - 1 == c_last[1]:
            return True
    elif c[1] == c_last[1]:
        if c[0] + 1 == c_last[0] or c[0] - 1 == c_last[0]:
            return True
    else:
        return False

def is_reverse_state(s_prime, s):
    """
        Function to check if s_prime is a reverse (R) state of s
    """
    c, cam = s
    c_prime, cam_prime = s_prime
    if cam == 1: # current state is a camouflage state
        if c == c_prime:
            return True
    elif cam == 0: # current state is a non-camouflage state
        if cam_prime == 1 and c == c_prime: # case when it switches from camouflage to non camouflage
            return True
        elif cam_prime == 0 and is_nearby_cell(c, c_prime): # case when it only changes cell but not camouflage state
            return True
    else:
        return False
    
def is_forward_state(s, s_prime):
    """
        Function to check if s_prime is a forward (F) state of s
    """
    # This is equivalent to say if s is a reverse state of s_prime
    return is_reverse_state(s, s_prime)


ending_time_grid = list(range(7, 16))

ending_time = 7
for ending_time in ending_time_grid:
    print('===========================')
    print('ending time is', ending_time)
    print('===========================')

    grid_size = 9
    ending_time = ending_time
    num_scenario = 500
    
    """
    Creating set
    """
    
    S = [(i, j) for i in range(1, grid_size + 1) for j in range(1, grid_size + 1)] # a state for searcher
    T = list(range(1, ending_time + 1)) # time period during which a detection could happen
    T0 = [0] + T # time period add t0
    Omega_num = list(range(1, num_scenario + 1)) # numerical list with each element represent a path number
    L = [1, 2] # set of searchers' type
    n_L = {1: 2, 2: 1} # number of searchers for each searcher type
    alpha_l = {1: 0.1, 2: 0.2} # detection rate for each searcher type
    # I = list(range(0, J * ending_time + 1))
    # print('i is', I)
    total_J = sum(n_L.values())
    searcher_init = {1: (1, 1), 2: (grid_size, 1)} # searcher type 1 starts from 1,1 # searcher type 2 starts from last row, left-most column
    
    tau = {1: 8, 2: 7} # operation duration limit for searcher of type l
    """ taking-off states """
# =============================================================================
#     S_plus = {}
#     for l in L:
#         S_plus[l] = [(0, 0)]
# =============================================================================
    
    """ landing states """
# =============================================================================
#     S_minus = {}
#     for l in L:
#         S_minus[l] = [(10_000, 10_000)]
# =============================================================================
    
    """ terminal state """
# =============================================================================
#     s_term = {}
#     for l in L:
#         s_term[l] = [(50_000, 50_000)]
# =============================================================================
    
# =============================================================================
#     S = S + S_plus + S_minus + s_term
# =============================================================================
    


    x = {}
    for l in L:
        for s in S:
            for t in T0:
                if s in searcher_init[l] and t == 0:
                    x[s, t] = n_L[l]
                else:
                    x[s, t] = 0
        

    """ Import data
    """
    if platform.system() == 'Windows':
        data_folder = 'E:\\Research\\Route_Optimization_for_Multiple_Searchers\\Python\\'
    else:
        data_folder = os.path.dirname(os.path.realpath(__file__))
    
    # zeta_raw = pd.read_csv(r'C:\Users\Wenbo Ma\Desktop\Route Optimization\Python\SP1-L\Zeta.csv', header = None, index_col = 0)
    zeta_raw = pd.read_csv(data_folder + '/zeta_1000_csv.csv', header = None, index_col = 0)
    for col in zeta_raw.columns:
        zeta_raw[col] = zeta_raw[col].apply(ast.literal_eval)
    
    Zeta = {}
    for path in range(1, zeta_raw.shape[0] + 1):
        for t in range(1, ending_time + 1):
            all_states = S
            for state in all_states:
                Zeta[(state, t, path)] = 0 # set Zeta equal to 0
            # cell_one_dim = zeta_raw.loc[path, 3 * (t - 1) + 1] # extract the occupied cell loc from Zeyu's path
            # cell_two_dim = (cell_one_dim // grid_size + 1, np.mod(cell_one_dim, grid_size))
            s, cam = zeta_raw.loc[path, t]
            if cam == 0: # not in camouflage mode
                Zeta[(s, t, path)] = 1 # set Zeta equal to 1 for occupied cell
    
    """ debugging """
# =============================================================================
#     test_path = 500
#     for t in T:
#         for s in S:
#             if Zeta[(s, t, test_path)] != 0:
#                 print(s, t)
# =============================================================================
    """ end of debugging """
    """ This is where we can set a subset to reduce c,t pair """
# =============================================================================
#     W = {}
#     for c in C:
#         for t in T:
#             if sum(Zeta[c, t, omega] for omega in Omega) >= 1:
#                 W[c, t] = 1
#             else:
#                 W[c, t] = 0
# =============================================================================
    # =============================================================================
    # Zeta = {}
    # for path in range(1, zeta_raw.shape[0] + 1):
    #     for t in range(1, ending_time + 1):
    #         all_cells = C
    #         for cell in all_cells:
    #             Zeta[(cell, t, path)] = 0 # set Zeta equal to 0
    #         cell_one_dim = zeta_raw.loc[path, 3 * (t - 1) + 1] # extract the occupied cell loc from Zeyu's path
    #         # cell_two_dim = (cell_one_dim // grid_size + 1, np.mod(cell_one_dim, grid_size))
    #         Zeta[(cell_one_dim, t, path)] = 1 # set Zeta equal to 1 for occupied cell
    # 
    # =============================================================================
    """
    """
    model_name = 'sp1_l_new_formul'
    m = gp.Model(model_name)
    m.setParam(GRB.Param.TimeLimit, 15 * 60)
    m.setParam(GRB.Param.Threads, 1)
    m.setParam(GRB.Param.LogFile, model_name)
    # m.setParam(GRB.Param.NumericFocus,1)
    
    n = {} # number of seachers per state per time
    for s in S:
        for t in T:
            n[s, t] = 2
            
    N = sum(n.values())
    
    """
    Creating parameters
    """
    np.random.seed(2022)
    q = np.random.uniform(low = 0, high = 1, size = num_scenario)
    q = q / sum(q) # normalize to a probablity distribution summing up to 1
    q = dict(zip(Omega_num, q))
    alpha = -3 * np.log(0.4) / 3
    
    # =============================================================================
    # q = pd.read_csv(data_folder + 'q.csv')
    # q = dict(zip(q['index'], q.q))
    # =============================================================================
    
    """
    Defining decision variables
    """
    sub_X = list(product(L, S, S, T0))
    sub_Z = list(product(L, S, T))
    sub_O = list(product(L, T))
    # sub_WW = list(product(S, T))
    # sub_WW = [each for each in sub_WW if W[each] == 1]
    
    X = m.addVars(sub_X, lb = 0, name = 'X')
    Z = m.addVars(sub_Z, lb = 0, ub = max(n_L.values()), vtype = GRB.INTEGER, name = 'Z')
    # Z_New = m.addVars(sub_WW, lb = 0, ub = J, vtype = GRB.INTEGER, name = 'Z')
    U = m.addVars(Omega_num, lb = 0, name = 'U')
    O = m.addVars(sub_O, lb = 0, name = 'O')
    
    """
    Defining objective function
    """
    m.setObjective(sum(q[omega_num] * U[omega_num] for omega_num in Omega_num), GRB.MINIMIZE)    
    

        
    coef_scale = 1
    
    """ check this constraint """
    """ constraint 3.1 """
    m.addConstrs((np.exp(-i * alpha) * (1 + i - i * np.exp(-alpha)) 
                  - np.exp(-i * alpha) * (1 - np.exp(-alpha)) * sum(Zeta[s, t, omega_num] * Z[l, s, t] for l in L for s in S for t in T) 
                  <= U[omega_num] for omega_num in Omega_num for i in range(1, N + 1)), name = 'cut') #2d
    # m.addConstrs((coef_scale * np.exp(-i * alpha) * (1 + i - i * np.exp(-alpha)) + coef_scale * np.exp(-i * alpha) * (np.exp(-alpha) - 1) * sum(Z_New[sub] for sub in sub_WW if Zeta[sub[0], sub[1], omega] == 1) <= coef_scale * U[omega] for omega in Omega for i in I), name = '19') #2d
    
    # 2.4b
    
    m.addConstrs((sum(X[l, s_prime, s, t - 1] for s_prime in S if is_nearby_cell(s, s_prime)) == sum(X[l, s, s_prime, t] for s_prime in S if is_nearby_cell(s, s_prime))  for l in L for s in S for t in T), name = 'continum') #2d
    
    # 2.4c
    """
    For this specific constraint, s is the default starting state at t = 0, which is defined by the number of searchers
    through xx[l, s, 0]. And s_prime (which is where the searchers could be at t = 1) need to be defined in a separate
    way from the general is_nearby_cell function
    """
    m.addConstrs((sum(X[l, s, s_prime, 0] for s_prime in S if is_nearby_cell(s, s_prime)) == x[s, 0] for s in S), name = 'continum_start') #2d
    
    # 2.4e
    """ 
    this s can only be at the searcher's init position 
    
    """
    m.addConstrs((sum(X[l, s, s_prime, t] for s in [searcher_init[l]] for s_prime in S  if is_nearby_cell(s, s_prime)) == O[l, t] for l in L for t in T), name = 'duration_start') #2d

    # 2.4f
    m.addConstrs((sum(X[l, s, s_prime, t] for s in [searcher_init[l]] for s_prime in S if is_nearby_cell(s, s_prime))  <= sum(O[l, t_prime] for t_prime in T if t_prime <= t and t_prime >= t - tau[l] + 1)  for l in L for t in T), name = 'duration') #2d

    
    # 2.5b
    m.addConstrs((sum(X[l, s_prime, s, t - 1] for s_prime in S if is_nearby_cell(s, s_prime)) == Z[l, s, t] for l in L for s in S for t in T), name = 'link_x_z') #2d
    # m.addConstrs((sum(X[c_prime, sub[0], sub[1] - 1] for c_prime in C if is_nearby_cell(sub[0], c_prime)) == Z_New[sub] for sub in sub_WW), name = '16') #2d
    
    # 2.5c
    m.addConstrs((sum(Z[l, s, t] for s in S) <= n[s, t] for l in L for s in S for t in T), name = 'capacity') #2d

    
    # 2.5d
    m.addConstrs((Z[l, s, t] <= min(n[s, t], n_L[l]) for l in L for s in S for t in T), name = 'capacity') #2d

    
    """ Solving
    """
    m.optimize()
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

