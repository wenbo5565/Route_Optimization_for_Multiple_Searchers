"""
Script for route optimization for multiple searchers
"""

"""
######### Import libraries
"""

import pandas as pd
import numpy as np
from itertools import product
import gurobipy as gp
from gurobipy import GRB
import os
import platform

"""
########## Define some helper functions
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

"""
    Setting parameters for the searching environment
"""

grid_size = 9
ending_time = 7
num_scenario = 200

if platform.system() == 'Windows':
    data_folder = 'E:\\Research\\Route_Optimization_for_Multiple_Searchers\\Python\\'
else:
    data_folder = os.path.dirname(os.path.realpath(__file__))


"""
Creating set for the optimization model
"""
C = [(i, j) for i in range(1, grid_size + 1) for j in range(1, grid_size + 1)]
T = list(range(1, ending_time + 1))
T0 = [0] + T
Omega = list(range(1, num_scenario + 1))
J = 3
I = list(range(0, J * ending_time + 1))
    
xx = {} # this is the variable x in the formulation
for c in C:
    for t in T0:
        if c == (1, 1) and t == 0:
            xx[c, t] = J
        else:
            xx[c, t] = 0
        

""" Import data
"""

# zeta_raw = pd.read_csv(r'C:\Users\Wenbo Ma\Desktop\Route Optimization\Python\SP1-L\Zeta.csv', header = None, index_col = 0)
zeta_raw = pd.read_csv(data_folder + '/Zeta.csv', header = None, index_col = 0)
Zeta = {}
for path in range(1, zeta_raw.shape[0] + 1):
    for t in range(1, ending_time + 1):
        all_cells = C
        for cell in all_cells:
            Zeta[(cell, t, path)] = 0 # set Zeta equal to 0
        cell_one_dim = zeta_raw.loc[path, 3 * (t - 1) + 1] # extract the occupied cell loc from Zeyu's path
        cell_two_dim = (cell_one_dim // grid_size + 1, np.mod(cell_one_dim, grid_size))
        Zeta[(cell_two_dim, t, path)] = 1 # set Zeta equal to 1 for occupied cell

np.random.seed(2022)
q = np.random.uniform(low = 0, high = 1, size = num_scenario)
q = q / sum(q) # normalize to a probablity distribution summing up to 1
q = dict(zip(Omega, q))
alpha = -3 * np.log(0.4) / J    
    
    
    
# =============================================================================
# W = {}
# for c in C:
#     for t in T:
#         if sum(Zeta[c, t, omega] for omega in Omega) >= 1:
#             W[c, t] = 1
#         else:
#             W[c, t] = 0
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
    
    
"""
Creating parameters
"""

    
    # =============================================================================
    # q = pd.read_csv(data_folder + 'q.csv')
    # q = dict(zip(q['index'], q.q))
    # =============================================================================
    
"""
Defining decision variables
"""
sub_X = list(product(C, C, T0))
sub_Z = list(product(C, T))
# =============================================================================
#     sub_WW = list(product(C, T))
#     sub_WW = [each for each in sub_WW if W[each] == 1]
# =============================================================================
    
X = m.addVars(sub_X, lb = 0, name = 'X')
Z = m.addVars(sub_Z, lb = 0, ub = J, vtype = GRB.INTEGER, name = 'Z')
# Z_New = m.addVars(sub_WW, lb = 0, ub = J, vtype = GRB.INTEGER, name = 'Z')
U = m.addVars(Omega, lb = 0, name = 'U')
    
"""
Defining objective function
"""
m.setObjective(sum(q[omega] * U[omega] for omega in Omega), GRB.MINIMIZE)    

"""
Defining constraints
"""

m.addConstrs((np.exp(-i * alpha) * (1 + i - i * np.exp(-alpha)) + np.exp(-i * alpha) * (np.exp(-alpha) - 1) * sum(Zeta[c, t, omega] * Z[c, t] for c in C for t in T) <= U[omega] for omega in Omega for i in I), name = '19') #2d
m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == sum(X[c, c_prime, t] for c_prime in C if is_nearby_cell(c, c_prime))  for c in C for t in T), name = '14') #2d
m.addConstrs((sum(X[c, c_prime, 0] for c_prime in C if is_nearby_cell(c, c_prime)) == xx[c, 0] for c in C), name = '15') #2d
m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == Z[c, t] for c in C for t in T), name = '16') #2d
# m.addConstrs((sum(X[c_prime, sub[0], sub[1] - 1] for c_prime in C if is_nearby_cell(sub[0], c_prime)) == Z_New[sub] for sub in sub_WW), name = '16') #2d

""" Solving
"""
m.optimize()

""" Extracting optimal solution, bound, objective value
"""
for key, val in Z.items():
    if val.X != 0:
        print('key is', key, 'value is', val.X)
        
print('best objective value is', m.objVal)
print('best bound is', m.objBound)
