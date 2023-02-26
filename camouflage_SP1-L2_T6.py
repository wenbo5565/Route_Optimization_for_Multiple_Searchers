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
# check

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

# ending_time_grid = list(range(7, 16))
# ending_time_grid = list(range(7, 16))
# ending_time= 7
J_total = [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50]


for J in J_total:
    J_2 = int(J * 0.7)
    J_1 = J - J_2
    print('===========================')
    print('J total is', J)
    print('J1 is', J_1, 'J2 is', J_2)
    print('===========================')

    grid_size = 9
    ending_time = 15
    num_scenario = 1000
    
    """
    Creating set
    """
    
    S = [(i, j) for i in range(1, grid_size + 1) for j in range(1, grid_size + 1)] # a state for searcher
    T = list(range(1, ending_time + 1)) # time period during which a detection could happen
    T0 = [0] + T # time period add t0
    Omega_num = list(range(1, num_scenario + 1)) # numerical list with each element represent a path number
    L = [1, 2] # set of searchers' type
    # L = [1]
    n_L = {1: J_1, 2: J_2} # number of searchers for each searcher type
    # n_L = {1: 1}
    # alpha_l = {1: 0.1, 2: 0.2} # detection rate for each searcher type
    # I = list(range(0, J * ending_time + 1))
    # print('i is', I)
    total_J = sum(n_L.values())
    
    s_init = (-100, -100)
    s_end = (-100, 100)
    on_map_init = (grid_size // 2, 1)
    on_map_end = (grid_size, grid_size // 2)
    
    # searcher_init = {1: (1, 1), 2: (3, 1)} # searcher type 1 starts from 1,1 # searcher type 2 starts from last row, left-most column
    
    S_expand = S + [s_init] + [s_end]
    
    tau = {1: ending_time * 0.8, 2: ending_time * 0.6} # operation duration limit for searcher of type l
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
        for s in S_expand:
            for t in T0:
                if s in [s_init] and t == 0:
                    x[l, s, t] = n_L[l]
                else:
                    x[l, s, t] = 0
    
    """ debug """
# =============================================================================
#     for l in L:
#         for t in T0:
#             for s in S_expand:
#                 if x[l, s, t] != 0:
#                     print(l, s, t, x[l, s, t])
# =============================================================================
    """ end of debug """
    
    
        

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
            all_states = S # not use S expand because this is for searcher
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
    D = {}
    for s in S:
        for t in T:
            if sum(Zeta[s, t, omega_num] for omega_num in Omega_num) >= 1: # if s is in cam under omega_num, Zeta = 0 due to the setup
                D[s, t] = 1
            else:
                D[s, t] = 0
                
    for t in T:
        D[s_init, t] = 1
        D[s_end, t] = 1
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
    # m.setParam(GRB.Param.MIPGap, 1e-4)
    # m.setParam(GRB.Param.NumericFocus,1)
    
    n = {} # number of seachers per state per time
    for s in S_expand:
        for t in T:
            n[s, t] = total_J
            
    # N = sum(n.values())
    N = sum(n_L.values()) * ending_time
    
    """
    Creating parameters
    """
    np.random.seed(2022)
    q = np.random.uniform(low = 0, high = 1, size = num_scenario)
    q = q / sum(q) # normalize to a probablity distribution summing up to 1
    q = dict(zip(Omega_num, q))
    alpha = -3 * np.log(0.4) / total_J
    
    # =============================================================================
    # q = pd.read_csv(data_folder + 'q.csv')
    # q = dict(zip(q['index'], q.q))
    # =============================================================================
    
    """
    Defining decision variables
    """
    sub_X = list(product(L, S_expand, S_expand, T0))
    sub_Z = list(product(L, S_expand, T))
    sub_Z = [each for each in sub_Z if D[each[1], each[2]] == 1]
    sub_O = list(product(L, T))
    # sub_WW = list(product(S, T))
    # sub_WW = [each for each in sub_WW if W[each] == 1]
    
    X = m.addVars(sub_X, lb = 0, name = 'X')
    Z = m.addVars(sub_Z, lb = 0, ub = max(n_L.values()), vtype = GRB.INTEGER, name = 'Z')
    # Z_New = m.addVars(sub_WW, lb = 0, ub = J, vtype = GRB.INTEGER, name = 'Z')
    U = m.addVars(Omega_num, lb = 0, name = 'U')
    O = m.addVars(sub_O, lb = 0, name = 'O')
    
# =============================================================================
#     test_state = (1, 1)
#     test_time = [t for t in range(1, 8)]
#     # for t in test_time:
#     Z[l, test_state, 1].ub = 1
#     Z[l, test_state, 1].lb = 1
# =============================================================================

    
    """
    Defining objective function
    """
    m.setObjective(sum(q[omega_num] * U[omega_num] for omega_num in Omega_num), GRB.MINIMIZE)    
    

        
    coef_scale = 1
    
    """ check this constraint """
    """ constraint 3.1 """
    print("**************** Adding the cut **************************")
    m.addConstrs((np.exp(-i * alpha) * (1 + i - i * np.exp(-alpha)) 
                  - np.exp(-i * alpha) * (1 - np.exp(-alpha)) * sum(Zeta[s, t, omega_num] * Z[l, s, t] for l in L for s in S for t in T if D[s, t] == 1) # D[s, t] is only defined in S not S_expand
                  <= U[omega_num] for omega_num in Omega_num for i in range(1, N + 1)), name = 'cut') #2d
    
    print("**************** End *************************************")
    # m.addConstrs((coef_scale * np.exp(-i * alpha) * (1 + i - i * np.exp(-alpha)) + coef_scale * np.exp(-i * alpha) * (np.exp(-alpha) - 1) * sum(Z_New[sub] for sub in sub_WW if Zeta[sub[0], sub[1], omega] == 1) <= coef_scale * U[omega] for omega in Omega for i in I), name = '19') #2d
    
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

    
    # 3.4c
    """ need to modify """
    m.addConstrs((sum(X[l, s_prime, s, t - 1] for s_prime in S_expand if is_backward_cell(s_prime, s, starting_c = s_init, ending_c = s_end, on_map_start = on_map_init, on_map_end = on_map_end))
                  == Z[l, s, t] for l in L for s in S_expand for t in T if D[s, t] == 1), name = 'link_x_z') #2d
    
# =============================================================================
#     for s in S_expand:
#         for s_prime in S_expand:
#             if X[1, s_prime, s, 0].X != 0:
#                 print(1, s_prime, s, 0, ':', X[1, s_prime, s, 0].X)
#     
#     Z[1, (1, 1), 1].X
#     
#     for s_prime in S_expand:
#         if is_backward_cell(s_prime, (1, 1), starting_c = s_init, ending_c = s_end,
#                             on_map_start = on_map_init, on_map_end = on_map_end):
#             print('backward', s_prime)
#         if is_backward_cell(s_prime, (1, 1), starting_c = s_init, ending_c = s_end,
#                             on_map_start = on_map_init, on_map_end = on_map_end):
#             print(s_prime, (1, 1), X[1, s_prime, (1, 1), 0].X) 
#     
#     for s in S_expand:
#         if Z[1, s, 1].X != 0:
#             print(l, s, 1, ':', Z[l, s, 1].X)
# =============================================================================
    
    # m.addConstrs((sum(X[c_prime, sub[0], sub[1] - 1] for c_prime in C if is_nearby_cell(sub[0], c_prime)) == Z_New[sub] for sub in sub_WW), name = '16') #2d
    
    # 3.4d
    m.addConstrs((sum(Z[l, s, t] for l in L) <= n[s, t] for s in S_expand for t in T if D[s, t] == 1), name = 'loc_capacity') #2d

    
    # 3.4e
    """ remove Z """
    m.addConstrs((Z[l, s, t] <= min(n[s, t], n_L[l]) for l in L for s in S_expand for t in T if D[s, t] == 1), name = 'capacity') #2d

    
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
    print("********** optimal solution for Z **********")
    sub_Z = sorted(sub_Z, key = lambda x: (x[0], x[2]))    
    for sub in sub_Z:
        if Z[sub].X != 0:
            print(sub, Z[sub].X)
            
# =============================================================================
#     print("********** optimal solution for X **********")
#     sub_X = sorted(sub_X, key = lambda x: (x[0], x[3]))    
#     for sub in sub_X:
#         if X[sub].X != 0:
#             print(sub, X[sub].X)
# =============================================================================
    print("********** number of possible looks *********")
    M = {}
    for omega_num in Omega_num:
        M[omega_num] = sum(Zeta[s, t, omega_num] * Z[l, s, t].X for l in L for s in S for t in T if D[s, t] == 1) 
    print('M by path is', M)
    print('M total is', sum(M.values()))
    
    
    print('********** optimal solution for O **********')
    sub_O = sorted(sub_O, key = lambda x: x[1])
    for sub in sub_O:
        if O[sub].X != 0:
            print(sub, O[sub].X)
