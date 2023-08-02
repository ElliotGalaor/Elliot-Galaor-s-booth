# -*- coding: utf-8 -*-
"""
Created on Mon May 23 16:35:20 2022

@author: ellio
"""

import matplotlib.pyplot as plt
import numpy as np

"""
P0=1
N=10000
n=133
t1=1666
t2=2083
t3=5833
a=n*P0/(t2-t1)
b=-a*t1/n
tau=-(t3-t2)/(n*np.log(1/2))
f1=lambda x:a*x+b
f2=lambda x:P0*np.exp(-(x-t2/n)/tau)
x=[]
t=[i*n/N for i in range(N)]
for i in range(0,t1):
    x.append(0)
for i in range(t1,t2):
    x.append(f1(i/n))
for i in range(t2,N):
    x.append(f2(i/n))

barx=np.array([t2*n/N,t3*n/N])
bary=np.array([1,1/2])


plt.plot(t,x)
plt.bar(barx,bary,width=0.1,color="grey")
plt.barh(bary,barx,height=3*10**(-3),color="grey")
plt.xlim(0,n)
plt.ylim(0,1.05)
plt.xlabel("Temps en µs")
plt.ylabel("P/P0")
plt.annotate('décroissance exponentielle',xy=(1, 1), xycoords='axes fraction',xytext=(-40, -160), textcoords='offset pixels',horizontalalignment='right',verticalalignment='bottom')
plt.show()
"""
""
plt.close()
N=10000
Q=4
w0=0.05
a=((w0/Q)**2+4*(w0)**4)**(1/2)
tau=w0/Q
A=0.03/a
n=100
t=[i*n/N for i in range(N)]
f=lambda x:A*(1-np.exp(-x/tau)*(np.cos(a*x)))
y=[f(i) for i in t]
plt.plot(t,y)
plt.xlim(0,100)
plt.ylim(0,2.5)
plt.grid()
plt.show()
print(y[0],max(y),min(y))
"""
plt.close()
N=10000
omega=1
omega0=1.25
Q=3
a=omega*omega0*(abs(1/Q-4))**(1/2)
tau1=16
tau2=1
A=1.05
B=1
n=100
t=[i*n/N for i in range(N)]
f=lambda x:A*np.exp(-x/tau1)*(np.sin(x*a))
g=lambda x:B*(1-np.exp(-x/tau2))
y=[f(i)+g(i) for i in t]
plt.plot(t,y)
plt.xlim(0,100)
plt.grid()
plt.show()
print(y[0],max(y),min(y),2*A/a)
"""