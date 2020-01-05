# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 12:44:42 2019

@author: risha
"""


import math


def rho(phi):
    return(math.cos(phi) + math.sqrt(1-0.25*math.sin(phi)*math.sin(phi)) + math.sqrt(1-(361/400)*math.sin(phi)*math.sin(phi)))

def rho_dash(phi):
    return((-math.cos(phi)*math.sin(phi)/(4*math.sqrt(1-0.25*math.sin(phi)*math.sin(phi)))) - (361*math.cos(phi)*math.sin(phi)/(400*math.sqrt(1-(361/400)*math.sin(phi)*math.sin(phi)))) - math.sin(phi))                     
def rho_ddash(phi):
    return(math.sin(phi)**2/(4*math.sqrt(1-math.sin(phi)**2/4))-math.cos(phi)**2/(4*math.sqrt(1-math.sin(phi)**2/4))-(math.cos(phi)**2*math.sin(phi)**2)/(16*(1-math.sin(phi)**2/4)**(3/2))+(361*math.sin(phi)**2)/(400*math.sqrt(1-(361*math.sin(phi)**2)/400))-(361*math.cos(phi)**2)/(400*math.sqrt(1-(361*math.sin(phi)**2)/400))-(130321*math.cos(phi)**2*math.sin(phi)**2)/(160000*(1-(361*math.sin(phi)**2)/400)**(3/2))-math.cos(phi))
def curvature(phi):
    return abs(rho(phi)**2 + 2*rho_dash(phi)**2 - rho(phi)*rho_ddash(phi))/ (rho(phi)**2+rho_dash(phi)**2)**(1.5)

def curvature_dash(phi):
    h = 0.000001
    return((curvature(phi+h)-curvature(phi))/h)

def curvature_ddash(phi):
    h = 0.000001
    return((curvature_dash(phi+h)-curvature_dash(phi))/h)

t = 0.00001
print("Simple inflections for the curve are:")
print("Points                         Value")
while(t < math.pi * 2):
    if(abs(curvature(t)) < 0.00005):
        if(abs(curvature_dash(t))!=0.0):
            print(str(t) + "             " + str(rho(t)))
            t = t + 0.2
    t = t + 0.00001

t = 0.1
print("")
print("Vertices for the curve are:")
print("Points                          Value")
while(t < math.pi * 2):
    if(abs(curvature_dash(t)) < 0.0005):
        if(abs(curvature_ddash(t))!=0.0):
            print(str(t) + "             " + str(rho(t)))
            t = t+0.2
    t = t + 0.00001