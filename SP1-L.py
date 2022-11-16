"""
Script for route optimization for multiple searchers
"""
import pandas as pd
import numpy as np
from itertools import product
import gurobipy as gp
from gurobipy import GRB
import os

grid_size = 9
ending_time = 15
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

""" Import data
"""
data_folder = os.path.dirname(os.path.realpath(__file__))

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

"""
"""
model_name = 'sp1_l'
m = gp.Model(model_name)
m.setParam(GRB.Param.TimeLimit, 15 * 60)
m.setParam(GRB.Param.Threads, 4)
m.setParam(GRB.Param.LogFile, model_name)



"""
Creating parameters
"""
q = np.random.uniform(low = 0, high = 1, size = num_scenario)
q = q / sum(q) # normalize to a probablity distribution summing up to 1
q = dict(zip(Omega, q))
alpha = -3 * np.log(0.4) / J


"""
Defining decision variables
"""
sub_X = list(product(C, C, T0))
sub_Z = list(product(C, T))

X = m.addVars(sub_X, lb = 0, name = 'X')
Z = m.addVars(sub_Z, lb = 0, ub = J, vtype = GRB.INTEGER, name = 'Z')
U = m.addVars(Omega, lb = 0, name = 'U')

"""
Defining objective function
"""
m.setObjective(sum(q[omega] * U[omega] for omega in Omega), GRB.MINIMIZE)    

"""
Defining constraints
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
    

m.addConstrs((np.exp(-i * alpha) * (1 + i - i * np.exp(-alpha)) + np.exp(-i * alpha) * (np.exp(-alpha) - 1) * sum(Zeta[c, t, omega] * Z[c, t] for c in C for t in T) <= U[omega] for omega in Omega for i in I), name = '19') #2d
m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == sum(X[c, c_prime, t - 1] for c_prime in C if is_nearby_cell(c, c_prime))  for c in C for t in T), name = '14') #2d
m.addConstrs((sum(X[c, c_prime, 0] for c_prime in C if is_nearby_cell(c, c_prime)) == J for c in C), name = '15') #2d
m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == Z[c, t] for c in C for t in T), name = '16') #2d

""" Solving
"""
m.optimize()