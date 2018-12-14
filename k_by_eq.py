# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 13:30:08 2018

@author: Emily
"""
import numpy as np
T = []
k = []
t = 250.0
print('please revise "kk" to the appropriate equation and "name" to the appropriate material')
while t < 2000.0:
    T.append(t)
    tau = t/1000.0
    kk = 1.158*(100/(7.5408+17.692*tau+3.6142*tau**2) + 6400/tau**(5/2)*np.exp(-16.35/tau)) #Uranium Equation
    k.append(round(kk,2))
    t = t + 1
    
    
name = 'UO2'
name = name + '.txt'
print('file name is:', name)

f = open(name, 'w+')
f.write('Temp(K)	    k(W/m-K)\n')
for i in range(0, len(T)):
    if len(str(T[i])) == 5:
        f.write(str(T[i]) + '	      ' + str(round(k[i],2)) + '\n')
    elif len(str(T[i])) == 6:
        f.write(str(T[i]) + '	     ' + str(round(k[i],2)) + '\n')
    else:
        f.write(str(T[i]) + '	   ' + str(round(k[i],2)) + '\n')