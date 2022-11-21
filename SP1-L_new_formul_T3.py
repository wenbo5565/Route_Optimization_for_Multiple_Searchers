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

scenario_grid = [1000, 2000, 4000, 8000, 16000, 32000]

for scenario in scenario_grid:
    print('===========================')
    print('num of scenarios', scenario)
    print('===========================')

    grid_size = 9
    ending_time = 8
    num_scenario = scenario
    
    """
    Creating set
    """
    
    
    C = [(i, j) for i in range(1, grid_size + 1) for j in range(1, grid_size + 1)]
    # C = list(range(1, grid_size * grid_size + 1))
    # =============================================================================
    # E = [(1,1)	,(1,2)	,(1,10)	,(2,1)	,(2,2)	,(2,3)	,(2,11)	,(3,2)	,(3,3)	,
    #      (3,4)	,(3,12)	,(4,3)	,(4,4)	,(4,5)	,(4,13)	,(5,4)	,(5,5)	,(5,6)	,
    #      (5,14)	,(6,5)	,(6,6)	,(6,7)	,(6,15)	,(7,6)	,(7,7)	,(7,8)	,(7,16)	,
    #      (8,7)	,(8,8)	,(8,9)	,(8,17)	,(9,8)	,(9,9)	,(9,18)	,(10,1)	,(10,10)	,
    #      (10,11)	,(10,19)	,(11,2)	,(11,10)	,(11,11)	,(11,12)	,(11,20)	,
    #      (12,3)	,(12,11)	,(12,12)	,(12,13)	,(12,21)	,(13,4)	,(13,12)	,
    #      (13,13)	,(13,14)	,(13,22)	,(14,5)	,(14,13)	,(14,14)	,(14,15)	,
    #      (14,23)	,(15,6)	,(15,14)	,(15,15)	,(15,16)	,(15,24)	,(16,7)	,
    #      (16,15)	,(16,16)	,(16,17)	,(16,25)	,(17,8)	,(17,16)	,(17,17)	,
    #      (17,18)	,(17,26)	,(18,9)	,(18,17)	,(18,18)	,(18,27)	,(19,10)	,
    #      (19,19)	,(19,20)	,(19,28)	,(20,11)	,(20,19)	,(20,20)	,(20,21)	,
    #      (20,29)	,(21,12)	,(21,20)	,(21,21)	,(21,22)	,(21,30)	,(22,13)	,
    #      (22,21)	,(22,22)	,(22,23)	,(22,31)	,(23,14)	,(23,22)	,(23,23)	,(23,24)	,
    #      (23,32)	,(24,15)	,(24,23)	,(24,24)	,(24,25)	,(24,33)	,(25,16)	,(25,24)	,
    #      (25,25)	,(25,26)	,(25,34)	,(26,17)	,(26,25)	,(26,26)	,(26,27)	,(26,35)	,
    #      (27,18)	,(27,26)	,(27,27)	,(27,36)	,(28,19)	,(28,28)	,(28,29)	,(28,37)	,
    #      (29,20)	,(29,28)	,(29,29)	,(29,30)	,(29,38)	,(30,21)	,(30,29)	,(30,30)	,
    #      (30,31)	,(30,39)	,(31,22)	,(31,30)	,(31,31)	,(31,32)	,(31,40)	,(32,23)	,(32,31)	,(32,32)	,(32,33)	,(32,41)	,
    #      (33,24)	,(33,32)	,(33,33)	,(33,34)	,(33,42)	,(34,25)	,(34,33)	,(34,34)	,
    #      (34,35)	,(34,43)	,(35,26)	,(35,34)	,(35,35)	,(35,36)	,(35,44)	,(36,27)	,(36,35)	,(36,36)	,(36,45)	,(37,28)	,
    #      (37,37)	,(37,38)	,(37,46)	,(38,29)	,(38,37)	,(38,38)	,(38,39)	,(38,47)	,(39,30)	,(39,38)	,(39,39)	,(39,40)	,
    #      (39,48)	,(40,31)	,(40,39)	,(40,40)	,(40,41)	,(40,49)	,(41,32)	,(41,40)	,(41,41)	,(41,42)	,(41,50)	,(42,33)	,
    #      (42,41)	,(42,42)	,(42,43)	,(42,51)	,(43,34)	,(43,42)	,(43,43)	,(43,44)	,(43,52)	,(44,35)	,(44,43)	,(44,44)	,
    #      (44,45)	,(44,53)	,(45,36)	,(45,44)	,(45,45)	,(45,54)	,(46,37)	,(46,46)	,(46,47)	,(46,55)	,(47,38)	,(47,46)	,
    #      (47,47)	,(47,48)	,(47,56)	,(48,39)	,(48,47)	,(48,48)	,(48,49)	,(48,57)	,(49,40)	,(49,48)	,(49,49)	,(49,50)	,
    #      (49,58)	,(50,41)	,(50,49)	,(50,50)	,(50,51)	,(50,59)	,(51,42)	,(51,50)	,(51,51)	,(51,52)	,(51,60)	,(52,43)	,
    #      (52,51)	,(52,52)	,(52,53)	,(52,61)	,(53,44)	,(53,52)	,(53,53)	,(53,54)	,(53,62)	,(54,45)	,(54,53)	,(54,54)	,
    #      (54,63)	,(55,46)	,(55,55)	,(55,56)	,(55,64)	,(56,47)	,(56,55)	,(56,56)	,(56,57)	,(56,65)	,(57,48)	,(57,56)	,
    #      (57,57)	,(57,58)	,(57,66)	,(58,49)	,(58,57)	,(58,58)	,(58,59)	,(58,67)	,(59,50)	,(59,58)	,(59,59)	,(59,60)	,
    #      (59,68)	,(60,51)	,(60,59)	,(60,60)	,(60,61)	,(60,69)	,(61,52)	,(61,60)	,(61,61)	,(61,62)	,(61,70)	,(62,53)	,
    #      (62,61)	,(62,62)	,(62,63)	,(62,71)	,(63,54)	,(63,62)	,(63,63)	,(63,72)	,(64,55)	,(64,64)	,(64,65)	,(64,73)	,
    #      (65,56)	,(65,64)	,(65,65)	,(65,66)	,(65,74)	,(66,57)	,(66,65)	,(66,66)	,(66,67)	,(66,75)	,(67,58)	,(67,66)	,
    #      (67,67)	,(67,68)	,(67,76)	,(68,59)	,(68,67)	,(68,68)	,(68,69)	,(68,77)	,(69,60)	,(69,68)	,(69,69)	,(69,70)	,
    #      (69,78)	,(70,61)	,(70,69)	,(70,70)	,(70,71)	,(70,79)	,(71,62)	,(71,70)	,(71,71)	,(71,72)	,(71,80)	,(72,63)	,
    #      (72,71)	,(72,72)	,(72,81)	,(73,64)	,(73,73)	,(73,74)	,(74,65)	,(74,73)	,(74,74)	,(74,75)	,(75,66)	,(75,74)	,
    #      (75,75)	,(75,76)	,(76,67)	,(76,75)	,(76,76)	,(76,77)	,(77,68)	,(77,76)	,(77,77)	,(77,78)	,(78,69)	,(78,77)	,
    #      (78,78)	,(78,79)	,(79,70)	,(79,78)	,(79,79)	,(79,80)	,(80,71)	,(80,79)	,(80,80)	,(80,81)	,(81,72)	,(81,80)	,
    #      (81,81)]
    # =============================================================================
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
        

    """ Import data
    """
    if platform.system() == 'Windows':
        data_folder = 'E:\\Research\\Route_Optimization_for_Multiple_Searchers\\Python\\'
    else:
        data_folder = os.path.dirname(os.path.realpath(__file__))
    
    # zeta_raw = pd.read_csv(r'C:\Users\Wenbo Ma\Desktop\Route Optimization\Python\SP1-L\Zeta.csv', header = None, index_col = 0)
    zeta_raw = pd.read_csv(data_folder + '/Zeta_' + str(num_scenario) + '.csv', header = None, index_col = 0)
    
    Zeta = {}
    for path in range(1, zeta_raw.shape[0] + 1):
        for t in range(1, ending_time + 1):
            all_cells = C
            for cell in all_cells:
                Zeta[(cell, t, path)] = 0 # set Zeta equal to 0
            cell_one_dim = zeta_raw.loc[path, 3 * (t - 1) + 1] # extract the occupied cell loc from Zeyu's path
            cell_two_dim = (cell_one_dim // grid_size + 1, np.mod(cell_one_dim, grid_size))
            Zeta[(cell_two_dim, t, path)] = 1 # set Zeta equal to 1 for occupied cell
    
    W = {}
    for c in C:
        for t in T:
            if sum(Zeta[c, t, omega] for omega in Omega) >= 1:
                W[c, t] = 1
            else:
                W[c, t] = 0
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
    np.random.seed(2022)
    q = np.random.uniform(low = 0, high = 1, size = num_scenario)
    q = q / sum(q) # normalize to a probablity distribution summing up to 1
    q = dict(zip(Omega, q))
    alpha = -3 * np.log(0.4) / J
    
    # =============================================================================
    # q = pd.read_csv(data_folder + 'q.csv')
    # q = dict(zip(q['index'], q.q))
    # =============================================================================
    
    """
    Defining decision variables
    """
    sub_X = list(product(C, C, T0))
    sub_Z = list(product(C, T))
    sub_WW = list(product(C, T))
    sub_WW = [each for each in sub_WW if W[each] == 1]
    
    X = m.addVars(sub_X, lb = 0, name = 'X')
    # Z = m.addVars(sub_Z, lb = 0, ub = J, vtype = GRB.INTEGER, name = 'Z')
    Z_New = m.addVars(sub_WW, lb = 0, ub = J, vtype = GRB.INTEGER, name = 'Z')
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
        
    coef_scale = 1
    # m.addConstrs((coef_scale * np.exp(-i * alpha) * (1 + i - i * np.exp(-alpha)) + coef_scale * np.exp(-i * alpha) * (np.exp(-alpha) - 1) * sum(Zeta[c, t, omega] * Z[c, t] for c in C for t in T) <= coef_scale * U[omega] for omega in Omega for i in I), name = '19') #2d
    m.addConstrs((coef_scale * np.exp(-i * alpha) * (1 + i - i * np.exp(-alpha)) + coef_scale * np.exp(-i * alpha) * (np.exp(-alpha) - 1) * sum(Z_New[sub] for sub in sub_WW if Zeta[sub[0], sub[1], omega] == 1) <= coef_scale * U[omega] for omega in Omega for i in I), name = '19') #2d
    
    m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == sum(X[c, c_prime, t] for c_prime in C if is_nearby_cell(c, c_prime))  for c in C for t in T), name = '14') #2d
    m.addConstrs((sum(X[c, c_prime, 0] for c_prime in C if is_nearby_cell(c, c_prime)) == xx[c, 0] for c in C), name = '15') #2d
    
    # m.addConstrs((sum(X[c_prime, c, t - 1] for c_prime in C if is_nearby_cell(c, c_prime)) == Z[c, t] for c in C for t in T), name = '16') #2d
    m.addConstrs((sum(X[c_prime, sub[0], sub[1] - 1] for c_prime in C if is_nearby_cell(sub[0], c_prime)) == Z_New[sub] for sub in sub_WW), name = '16') #2d
    
    
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

