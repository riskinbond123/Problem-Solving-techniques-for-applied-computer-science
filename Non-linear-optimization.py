# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 12:50:12 2019

@author: risha
"""

x1=0.2
x2=0.2
def gT_dash(x1, x2):
    a = 3*(-grad_f(x1,x2)[1]**2*grad_f(x1,x2)[0] - (1/3.0)*grad_f(x1,x2)[0]**3)
    b = 2*(grad_f(x1,x2)[0]**2*x1 + grad_f(x1,x2)[1]**2*x1 + 4*x1*x2**2*grad_f(x1,x2)[0])
    c = -grad_f(x1,x2)[0] * x1**2 - grad_f(x1,x2)[1]**2 - x2**2 *grad_f(x1,x2)[0] + 3*grad_f(x1,x2)[0]
    t1 = (-b + (b**2 - 4*a*c) ** (1/2.0)) / (2*a)
    t2 = (-b - (b**2 - 4*a*c) ** (1/2.0)) / (2*a)
    return (t1, t2)

def f(x1,x2):
    return x1**3/3.0 + x2**2 * x1 - 3*x1

def grad_f(x1,x2):
    return(x1**2 + x2**2 - 3, 2*x1*x2);
print("Performing a steepest descent:")
print("x1   x2")
m = 1
fv = f(x1, x2)
while abs(m) > 10**(-6):
    t = gT_dash(x1, x2)
    grad = grad_f(x1, x2)
    x1 -= t[0] * grad[0]
    x2 -= t[0] * grad[1]
#   x1 -= t[1] * grad[0]       // Descent
#   x2 -= t[1] * grad[1]       // Descent
    m = fv - f(x1, x2)
    fv = f(x1, x2)
    print(x1,x2)
print(" ")
print("Final objective value is:" )
print(fv)
print("The values of x1 and x2 are-:")
print(x1,x2)
print(" ")
x1=0.2
x2=0.2

print("Performing a steepest ascent:")
print("x1   x2")
m = 1
fv = f(x1, x2)
while abs(m) > 10**(-6):
    t = gT_dash(x1, x2)
    grad = grad_f(x1, x2)
    x1 -= t[1] * grad[0]
    x2 -= t[1] * grad[1]
#   x1 -= t[1] * grad[0]
#   x2 -= t[1] * grad[1]
    m = fv - f(x1, x2)
    fv = f(x1, x2)
    print(x1,x2)
print(" ")
print("Final objective value is:" )
print(fv)
print("The values of x1 and x2 are-:")
print(x1,x2)
