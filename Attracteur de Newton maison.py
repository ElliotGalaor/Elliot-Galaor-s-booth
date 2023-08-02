# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 14:14:31 2022

@author: ellio
"""

import numpy as np
import matplotlib.pyplot as plt

pi=np.pi
cos=lambda x:np.cos(x)
sin=lambda x:np.sin(x)

def Polynome1(a):
    return lambda x:x-a

def Polynome2(a,b):
    return lambda x:(x-a)*(x-b)

def Polynome3(a,b,c):
    return lambda x:(x-a)*(x-b)*(x-c)

def Polynome4(a,b,c,d):
    return lambda x:(x-a)*(x-b)*(x-c)*(x-d)

def Polynome5(a,b,c,d,e):
    return lambda x:(x-a)*(x-b)*(x-c)*(x-d)*(x-e),

def dPolynome1(a):
    return lambda x:x-a

def dPolynome2(a,b):
    f,g=Polynome1(a),Polynome1(b)
    return lambda x:f(x)+g(x)

def dPolynome3(a,b,c):
    f,g,h=Polynome2(a,b),Polynome2(b,c),Polynome2(a,c)
    return lambda x:f(x)+g(x)+h(x)

def dPolynome4(a,b,c,d):
    f,g,h,i=Polynome3(a,b,c),Polynome3(b,c,d),Polynome3(a,c,d),Polynome3(a,b,d)
    return lambda x:f(x)+g(x)+h(x)+i(x)

def dPolynome5(a,b,c,d,e):
    f,g,h,i,j=Polynome4(a,b,c,d),Polynome4(b,c,d,e),Polynome4(a,c,d,e),Polynome4(a,b,d,e),Polynome4(a,b,c,e)
    return lambda x:f(x)+g(x)+h(x)+i(x)+j(x)

#           Choix des racines
a=cos(pi/3)+sin(pi/3)*1j
b=cos(pi/3)-sin(pi/3)*1j
c=-1
P=Polynome3(a,b,c)
dP=dPolynome3(a,b,c)
i=10
def newton1d(F,DF,x0,eps=1e-15,i):
    x=x0
    for i in range(N):
        Fx=F(x)
        DFx=DF(x)
        if abs(Fx)<eps:
            return x
        if abs(DFx)<eps:
            raise Exception(f"La dérivée DF = {DFx} est trop petite")
        x-=Fx/DFx
    raise Exception(f"L'erreur après {N} itération est {abs(Fx)} > {eps}")



