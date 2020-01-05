# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 16:53:37 2019

@author: risha
"""
import math
import numpy as np
def evalPoly(n, a, t, v):
    eval1 = 0 
    eval2 =0
    for i in range(0,n):
        eval2 = eval2*t + eval1
        eval1 = eval1*t + a[i]
    v[0]= eval1 
    v[1]= eval2
#    print("Polynomial evaluation is: "+ str(v[0]) + ". The evaluated derivative of the polynomial at that point is " + str(v[1]))
    return v[0],v[1]
q_1b = [1,0,-170,0,7392,0,-39712,0,51200]
t_1b_1 = 1.414214
t_1b_2 = 1+2j


print("Question 1b :")
print(evalPoly(9, q_1b, t_1b_1,[0,0]))
print("              ")
print(evalPoly(9, q_1b, t_1b_2,[0,0]))
print("              ")


def DFT(a,n):
    ahat = [complex(0, 0)] *n
    if(n==1):
        return a
    omega_n = math.e**(1j*2*math.pi/n)
    omega = 1
    a_even = a[0::2]
    a_odd = a[1::2]
    ahat_even = DFT(a_even,int(n/2))
    ahat_odd = DFT(a_odd,int(n/2))
    for k in range(0,int(n/2)):
        ahat[k]= ahat_even[k] + omega*ahat_odd[k]
        ahat[k+int(n/2)] = ahat_even[k] - omega*ahat_odd[k]
        omega = omega*omega_n
    return ahat
# Changed sign
def inverse_DFT(a,n):
    ahat = [0]*n
    if(n==1):
        return a
    omega_n = math.e**(-1j*2*math.pi/n)
    omega = 1
    a_even = a[0::2]
    a_odd = a[1::2]
    ahat_even = DFT(a_even,int(n/2))
    ahat_odd = DFT(a_odd,int(n/2))
    for k in range(0,int(n/2)):
        ahat[k]= ahat_even[k] + omega*ahat_odd[k]
        ahat[k+int(n/2)] = ahat_even[k] - omega*ahat_odd[k]
        omega = omega*omega_n
    return [d/n for d in ahat]

# We need to pad a 0 coefficant to make the degree equal

def multiplyPolys(a, b):
    align = 1
    while align < 2* (max(len(a), len(b))):
        align *= 2
    a.extend([0]*(align - len(a)))
    b.extend([0]*(align - len(b)))

    a_hat = DFT(a,len(a))
    b_hat = DFT(b,len(b))
    c_hat = [a_hat[i] * b_hat[i] for i in range(align)]
    c = inverse_DFT(c_hat,len(c_hat))

    return [round(t.real, 2) for t in c[:(len(a) + len(b) - 1)]]


q_2d_p = [-6.8, 10.8, -10.8, 7.4, -3.7, 2.4, -70.1, 1]
q_2d_q = [51200, 0, -39712, 104.2, 7392, 0.614, -170, 0, 1]
print("              ")
print("Question 2d) The coefficants of the polynomial multiplication in ascending order are:")
print(multiplyPolys(q_2d_p, q_2d_q))  
print("  ")

def cbrt(x):
    if x == 0:
        return 0
    return x/abs(x) * (abs(x) ** (1./ 3.))

def cubic(*args):
    if len(args) == 4:
        a3,a2,a1,a0 = args[0],args[1],args[2],args[3]
        p,q,r = float(a2) / a3 ,float(a1) / a3,  float(a0) / a3
        a,b = 1/3 * (3*q - p**2), 1/27 * (2*p**3 - 9*p*q + 27*r)
        A, B = cbrt((-b/2 + (b**2/4 + a**3/27)**(1./2))),cbrt((-b/2 - (b**2/4 + a**3/27)**(1./2)))
    elif len(args) == 3:
        p,q,r = args[0],args[1],args[2]
        a,b = 1/3 * (3*q - p**2), 1/27 * (2*p**3 - 9*p*q + 27*r)
        A,B = cbrt((-b/2 + (b**2/4 + a**3/27)**(1./2))), cbrt((-b/2 - (b**2/4 + a**3/27)**(1./2)))
    elif len(args) == 2:
        p = q = 0
        a,b = args[0], args[1]
        A,B = cbrt((-b/2 + (b**2/4 + a**3/27)**(1./2))),cbrt((-b/2 - (b**2/4 + a**3/27)**(1./2)))
        
    y1 = A + B
    y2 = complex(-1/2*(A + B), 3**(1/2)/2 * (A - B))
    y3 = complex(-1/2*(A + B), -3**(1/2)/2 * (A - B))

    return y1 - p/3, y2 - p/3, y3 - p/3


def quartic(*args):
    if len(args) == 5:
        a4,a3,a2,a1,a0 = args[0],args[1],args[2],args[3],args[4]
        p,q,r,s = a3/ a4, a2/ a4, a1/ a4, a0 / a4
    elif len(args) == 4:
        p,q,r,s = args[0],args[1],args[2],args[3]
    
    rr = cubic(-q, p*r - 4*s, 4*q*s - r**2 - p**2*s)[0]

    R = (1/4*p**2 - q + rr)**(1./2)
    if R == 0:
        D = (3/4*p**2 - 2*q + 2*(rr**2 - 4*s)**(1./2))**(1./2)
        E = (3/4*p**2 - 2*q - 2*(rr**2 - 4*s)**(1./2))**(1./2)
    else:
        D = (3/4*p**2 - R**2 - 2*q + 1/4*(4*p*q - 8*r - p**3)/R)**(1./2)
        E = (3/4*p**2 - R**2 - 2*q - 1/4*(4*p*q - 8*r - p**3)/R)**(1./2)

    x1 = -p/4 + 1/2 * (R + D)
    x2 = -p/4 + 1/2 * (R - D)
    x3 = -p/4 - 1/2 * (R - E)
    x4 = -p/4 - 1/2 * (R + E)

    return x1, x2, x3, x4

print("Question 3b. part A has the following roots for the cubic polynomial:")
print(cubic(110, -23, 87, 4))
print("  ")
print("Question 3b. part B has the following roots for the quartic polynomial:")
print(quartic(43, 1.34, -7, 0, -3400))


accuracy = 0.0000001

def mullerRoot(p0):
    p = list(p0)
    n = len(p) - 1
    a0 = p[0]
    a1 = p[1] + 0.00000001
    an = p[-1]
    radius = min(n*abs(a0/a1), (abs(a0/an))**(1./n))
    r = [0.0] * 3
    r[0] = radius * 0.25
    r[1] = radius * 0.5
    r[2] = radius
    while abs(r[-1] - r[-2]) > accuracy:
        pr = list(p)
        pr.reverse()
        c = evalPoly(len(pr),pr, r[-1],[0,0])[0]
        f23 = (evalPoly(len(pr),pr, r[-1],[0,0])[0] - evalPoly(len(pr),pr, r[-2],[0,0])[0]) / (r[-1] - r[-2])
        f12 = (evalPoly(len(pr),pr, r[-2],[0,0])[0] - evalPoly(len(pr),pr, r[-3],[0,0])[0]) / (r[-2] - r[-3])
        f13 = (f12 - f23) / (r[-3] - r[-1])
        b = f23 + f13*(r[-1] - r[-2])
        a = f13
        x1 = -2*c / (b + (b**2 - 4*a*c)**(1./2))
        x2 = -2*c / (b - (b**2 - 4*a*c)**(1./2))
        if (abs(x1) > abs(x2)):
            res = x2 + r[-1]
        else:
            res = x1 + r[-1]
        r.append(res)
    return r[-1]


def Newton(a, r):
    ar = list(a)
    ar.reverse()
    v = evalPoly(len(a),ar, r,[0,0])
    rn = r - v[0]/v[1]
    while abs(rn - r) > accuracy:
        r = rn
        v = evalPoly(len(a),a, r,[0,0])
        rn = r - v[0]/v[1]
    return rn


def mullerMethod(a):
    r = []
    if len(a) - 1 == 4:
        r.extend(quartic(a[0], a[1], a[2], a[3], a[4]))
    elif len(a) - 1 == 3:
        r.extend(cubic(a[0], a[1], a[2], a[3]))
    elif len(a) - 1 > 4:
        p0 = list(a)
        p1 = list(a)
        l = 0
        while len(p1) - 1 > 4:
            r.append(complex(mullerRoot(p1)))
            if abs(r[l].imag) < accuracy:
                p1 = np.polydiv(np.array(p0), np.array([1, -r[l]]))[0].tolist()
                p0 = p1
                l = l + 1
                p1r = list(p1)
                p1r.reverse()
                while abs(evalPoly(len(p1r), p1r, r[l-1], [0,0])[0]) < accuracy:
                    r.append(r[l-1])
                    p1 = np.polydiv(np.array(p0), np.array([1, -r[l]]))[0].tolist()
                    p0 = p1
                    l = l + 1
            else:
                conjugate_r = r[l].conjugate()
                r.append(conjugate_r)
                p1 = np.polydiv(np.array(p0), np.array([1, -2*r[l].real, r[l].real**2 + r[l].imag**2]))[0].tolist()
                p0 = p1
                l = l + 2
                p1r = list(p1)
                p1r.reverse()
                while abs(evalPoly(len(p1r),p1r, r[l-1])[0]) < accuracy:
                    r.append(r[l-2])
                    r.append(r[l-2])
                    p1 = np.polydiv(np.array(p0), np.array([1, -2*r[l].real, r[l].real**2 + r[l].imag**2]))[0].tolist()
                    p0 = p1
                    l = l + 2
        if len(p1) - 1 == 4:
            r.extend(quartic(p1[0], p1[1], p1[2], p1[3], p1[4]))
        elif len(p1) -1 == 3:
            r.extend(cubic(p1[0], p1[1], p1[2], p1[3]))

    r_final = []
    for i in r:
        r_final.append(complex(round(Newton(a, i).real, 8), round(Newton(a, i).imag, 8)))

    return r_final

p = [1, -3.7, 7.4, -10.8, 10.8, -6.8]
q = [1, -0.843121, -8.35979, 10.1887, 14.6196, -25.7634, 9.15636, -0.360995, -0.180591, 0.00787276]
print("  ")
print("  ")
print("  ")
print("Question 4b, polynomial p has the following roots:")
print(mullerMethod(p))
print("  ")
print("  ")
print("  ")
print("Question 4b, polynomial q has the following roots:")
print(mullerMethod(q))
































