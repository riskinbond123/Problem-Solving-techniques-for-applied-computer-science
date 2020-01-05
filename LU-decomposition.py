# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 19:11:04 2019

@author: risha
"""

"""
The first part solves for a LU decomposition using Crout's algorithm
"""

import numpy as np

import pprint

A=[[11,2,-5,6,48],[1,0,17,29,-21],[-3,4,55,-61,0],[41,97,-32,47,23],[-6,9,-4,-8,50]]

b1 = [4,0,-7,-2,-11]
b2 = [2, 77,-1003,-7,10]

def LUdcmp(A):
    global L1
    global U1 
    global P1
    n = len(A)
    P = np.eye(5).tolist()      #Initialize identity matrix and coonvert to list
    for j in range(len(A)):
        row = max(range(j, n), key=lambda i: abs(A[i][j]))
        if j != row:
            P[j], P[row] = P[row], P[j]
        n = len(A)
    L = np.zeros([5,5]).tolist()
    U = np.zeros([5,5]).tolist()
    A2 = np.matmul(P, A)
    for j in range(n):
        L[j][j] = 1.0
        for i in range(j+1):
            s1 = sum(U[k][j] * L[i][k] for k in range(i))
            U[i][j] = A2[i][j] - s1
        for i in range(j, n):
            s2 = sum(U[k][j] * L[i][k] for k in range(j))
            L[i][j] = (A2[i][j] - s2) / U[j][j]
    
    print("Programming assignment 4 Question 2a:")
    print(" ")
    print("             A                  ")
    pprint.pprint(A)
    print(" ")
    print("             L                ")
    pprint.pprint(L)
    print(" ")
    print("             P                  ")
    pprint.pprint(P)
    print(" ")
    print("             U                 ")
    pprint.pprint(U)
    print(" ")
    
    L1 , U1, P1 = L,U,P
    return 0

LUdcmp(A)


def fwd_sub(L, PB):
    y = [0] * len(PB)
    for i in range(len(PB)):
        y[i] = PB[i] - sum([(L[i][j]*y[j]) for j in range(i)])
    return np.array(y)


def back_sub(U, y):
    x = [0] * len(y)
    for i in range(len(y)):
        i = len(y) - i - 1
        x[i] = (y[i] - sum([(U[i][j]*x[j]) for j in range(i, len(y))])) / U[i][ i]
    return np.array(x)


def LUbksub(A, b):
    y = fwd_sub(L1, np.matmul(P1,b))
    x = back_sub(U1, y)
    print(x)
    return x
print("Programming assignment 4 Question 2b part a , x value::")
LUbksub(A,b1)


print(" ")
print("Programming assignment 4 Question 2b part b , x value:")
LUbksub(A,b2)
# np.matmul(A,LUbksub(A,b2)) ... cross-checked


































