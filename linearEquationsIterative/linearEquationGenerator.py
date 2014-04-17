'''
Created on Apr 16, 2014

@author: NASSAR
'''

import pickle 
import numpy as np 
path_local="local/"

A = np.array([[3,1], [1,2]])
b = np.array([9,8])
x = np.linalg.solve(A, b)

np.save(path_local+"matrixA" , A)
np.save(path_local+"vectorb", b)
print A
print x
print b
