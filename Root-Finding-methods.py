# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 12:00:02 2019

@author: risha
"""
import math
delta = 0.00005
arb_iterations = 100
def func(x): 
    return math.e**(-x) - math.sin(2*x)

def diff_func(x):
    return -1*math.e**(-x) - math.cos(2*x) * 2
print(" ")
print("Bisection method to calculate the roots between 0 and 1")
print(" ")
def bisection(a,b): 
    print("i  [ai,bi]  f(ai)  f(bi) (ai+bi)/2   f((ai+bi)/2))")
    for i in range(arb_iterations):
        m = (a+b) / 2
        if func(a)*func(m) <= 0:
            b = m
        else:
            a = m
        print(i, (a, b), func(a), func(b), (a+b)/2 , func((a+b)/2))
        if abs(a - b) < delta:
            return 0

bisection(0,1)
    
print(" ")
print("Modified regular Falsi method to calculate the roots between 0 and 1")
print(" ")
def MRF(a,b):
    print("i  [ai,bi]  F  G   wi+1   f(wi+1)")
    F, G, w, a_1, b_1 = func(a), func(b), a, a, b
    for i in range(arb_iterations):
        w_1 = (G*a_1 - F*b_1)/(G-F)
        if func(a_1)*func(w_1) <=0:
            b_1 = w_1
            G = func(w_1)
            if func(w)*func(w_1) > 0:
                F = F/2
        else:
            a_1 = w_1
            F = func(w_1)
            if func(w)*func(w_1) > 0:
                G = G/2
        print(i, (a_1, b_1), F , G , w_1 , func(w_1))
        if (abs(func(b_1))<delta):
            return 0

MRF(0,1)
print(" ")
print("Secant method to calculate the roots between 0 and 1")
print(" ")
def Secant(a,b):
    a_1 = a
    b_1 = b
    print("i    xi     f(xi)")
    for i in range(arb_iterations):
        c = b_1 - func(b_1)*(b_1 - a_1)/(func(b_1) - func(a_1))
        a_1, b_1 = b_1, c
        print(i,c,func(c))
        if (abs(func(c)) < delta):
            return 0
Secant(0,1)

print(" ")
print("Newton's method to calculate the roots between 0 and 1")
print(" ")
def Newtons(a):
    x = a
    print("i    xi     f(xi)")
    for i in range(arb_iterations):
        x_n = x - (func(x) / diff_func(x) )
        print(i,x_n, func(x_n))
        if (abs(x_n - x) <= delta):
            return x_n
        x = x_n
        
Newtons(0.1)















