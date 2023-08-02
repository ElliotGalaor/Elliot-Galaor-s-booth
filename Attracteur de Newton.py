# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 15:49:34 2021

@author: ellio
"""
import numpy as np
import matplotlib.pyplot as plt
def newton1d(F,DF,x0,eps=1e-15,N=100):       #Fonction appliquant la méthode de Newton en chaque point du tableau
    x=x0
    for i in range(N):
        Fx=F(x)
        DFx=DF(x)
        if abs(Fx)<eps:
            return x
        if abs(DFx)<eps:
            raise Exception(f"La dérivée DF = {DFx} est trop petite")
        x-=Fx/DFx
    if abs(Fx)>eps:
        return 0        #raise Exception(f"L'erreur après {N} itération est {abs(Fx)} > {eps}")
 
def parallel(z0,eps=1e-3,N=100):
    z=z0.copy()
    cond=np.abs(z**3-1)>eps
    for i in range(N):
        cond[cond]=np.abs(z[cond]**3-1)>eps
        if cond.any()==False:
            return z
        z[cond]=1/(3*z[cond]**2)+2*z[cond]/3
    raise Exception(f"Certains z n'ont pas convergés")
    
   
F=lambda z:z**5-1       #Polynome
DF=lambda z:5*z**4        #Polynome dérivé

[newton1d(F,DF,z0) for z0 in [-1,-1j,1j,1]]
nb=8000        #Nb de points pris dans le tableau
lst=np.linspace(-1,1,nb)
x0,y0=np.meshgrid(lst,lst)
z0=x0+1j*y0
out=parallel(z0)


out=z0.copy()
for i in range(nb):
    for j in range(nb):
        out[i,j]=newton1d(F,DF,z0[i,j])
plt.figure(figsize=(10,10))
plt.imshow(np.angle(out),extent=[-3,3,-3,3])