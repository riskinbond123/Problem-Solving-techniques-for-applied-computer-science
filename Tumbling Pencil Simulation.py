# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 14:03:00 2019

@author: risha
"""
import math
import numpy as np

g = np.array([[0,0,0,-9.8]]).transpose()

"Input values"
pi = math.pi
m = 1
r = 0.5
h1 = 3
h2 = 0.5
t1 = 0.4
t2 = 0.8

h = (6*(h1*h1) + 12*(h1)*(h2) + 3*(h2*h2))/((12*h1)+(4*h2))
p1 = h*np.array([[0,math.cos(pi/3),0,math.sin(pi/3)]]).transpose()

v1_neg = -5*np.array([[0,math.cos(pi/6),0,math.sin(pi/6)]]).transpose()
w1_neg = np.array([[0,1,5,0.5]]).transpose()       # write it as a quaternion
v1_pos = np.array([[0,-1.80954,-0.546988,1.2076]]).transpose()
w1_pos = np.array([[0,0.09957,-0.04174,0.5]]).transpose()

l = (h1/2) + h2 - h

"""
m1 and m2 need to be calculated wwhere m1 is mass of cylinder and m2 is mass
of cone.  m1 = (vol(cylinder)/(total volume))* m
          m2 = (vol(cone)/(total volume))*m
"""
m1 = (pi*r*r*h1/((1/3)*pi*r*r*h2 + pi*r*r*h1))*m     # mass of cylinder
m2 = ((((1/3)*pi*r*r*h2))/((1/3)*pi*r*r*h2 + pi*r*r*h1))*m  #mass of cone



" I will construct matrix Q below"

Q = np.zeros((3,3))

Q[0,0]= (m/(h1 + h2/3))*(h1*(((3*(r*r) + (h1*h1))/12) + l*l) + (h2/3)*((3/5)*((r*r)/4 + h2*h2 ) + h*h))   
Q[1,1]= (m/(h1 + h2/3))*(h1*(((3*(r*r) + (h1*h1))/12) + l*l) + (h2/3)*((3/5)*((r*r)/4 + h2*h2 ) + h*h))
Q[2,2]= ((1/2)*m1 + (3/10)*m2)*(r*r)


def angular_acceleration(w):
    return -1*np.matmul((np.linalg.inv(Q)),(np.cross(w,np.matmul(Q,w),axis=0)))

"""
I am going to represent my quaternions in 1-d arrays.
The quaternion operations such as multiply, 
"""
q1 = np.array([math.cos(pi/12), 0 , math.sin(pi/12), 0 ])

def quaternion_conjugate(quaternion):
    return [1,-1,-1,-1]*quaternion

def quaternion_multiply(quaternion1, quaternion0):
    w0, x0, y0, z0 = quaternion0
    w1, x1, y1, z1 = quaternion1
    return np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
                     x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                     -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                     x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0], dtype=np.float64)

def Lq(q,v):
    return quaternion_multiply(quaternion_multiply(q,v), quaternion_conjugate(q))

"""
Quaternion code tested on example given in handout.
Lq(np.array([0.5,0.5,0.5,0.5]),np.array([0,1,0,0]))
Out[97]: array([0., 0., 1., 0.])

"""

def norm(omega):
    return math.sqrt(quaternion_multiply(quaternion_conjugate(q1).reshape(4,),q1.reshape(4,))[0])

def tumble(p,q,v,w,T):
    t = 0.0
    v = v
    w = w
    q = q
    p = p
    count = 0
    print(p[1:])
    print("\n")
    print(v[1:])
    print("\n")
    print(w[1:])
    print("\n")
    print("\n")
    delta_t = 0.00001              # step size
    while(abs(t) < abs(T)):
        omega = Lq(q1, w.reshape(4,))
        omega_normalized = omega / norm(omega)
        if(T>0):
            delta_phi = norm(omega)*delta_t
        else:
            delta_phi = (-1 )* norm(omega)*delta_t
        r = np.array([math.cos(delta_phi/2),0,0,0]) + omega_normalized*math.sin(delta_phi/2)
        q = quaternion_multiply(r, q)
        if(T>0):
            v[1:] = v[1:] + g[1:]*(delta_t)                # step size
            w[1:] = w[1:] - (delta_t)*angular_acceleration(w[1:])
            print(angular_acceleration(w[1:]))
            t = t + delta_t
        else:
            v[1:] = v[1:] - g[1:]*(delta_t)                # step size
            w[1:] = w[1:] + (delta_t)*angular_acceleration(w[1:])
            print(angular_acceleration(w[1:]))
            t = t - delta_t

        p[1:] = p[1:] + v[1:]*delta_t + ((delta_t*delta_t)/2)*(g[1:])
        count = count + 1

        if(count % 10000 == 0):
            print(p[1:])
            print("\n")
            print(v[1:])
            print("\n")
            print(w[1:])
            print("\n")
            print("\n")

tumble(p1,q1,v1_neg,w1_neg,-0.4)
tumble(p1,q1,v1_pos,w1_pos,0.4)







  
        
    