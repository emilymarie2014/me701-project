# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import re
import sys


def k_array(a):
    f = open('%s' % a , 'r')
    s = f.readlines()
    f.close()
    T = []
    k = []
    for j in range(0,len(s)):
        p = '(\d+)'
        l = re.search(p, s[j])
        if l != None:
            break 
            
    for i in range(j,len(s)):
        T.append(float(s[i].split()[0]))
        k.append(float(s[i].split()[1]))
        
    return T, k

def k_value(T, T_array, k_array):
    for i in range(0,len(k_array)):
        if T_array[i] == T:
            k = k_array[i]
            return k
        elif T_array[i] > T:
            k = (k_array[i]-k_array[i-1])/(T_array[i]-T_array[i-1])*(T-T_array[i-1])
            return k
    if T_array[-1] < T:
        k = k_array[-1]
        return k

def q_gen(n):
    
    if n == 'Y':
        q = '200*t**2 - 2*t'
    else:
        q = raw_input("What is qgen: ")
    qgen = []
    for i in range(0,len(time)):
        t = time[i]
        qgen.append(eval(q))
    return qgen
    


def iterate(t, k_T, hr_T):
    it = 0
    g = 0
    if g > 10000:
        print('*********BREAK***********')
        return
    while it == 0:

        g = g + 1
        print(g)    
        A = np.zeros((node, node), dtype='float64')
        h_r = eps*sig*(hr_T[1] + hr_T[0])*(hr_T[1]**2 + hr_T[0]**2)
        A[0, 0] = 1 + delta_t/(rho_gas*delta_x*c_gas)*h_r
        A[0, 1] = -delta_t/(rho_gas*delta_x*c_gas)*h_r
        

        #node 0
        k = k_value(k_T[1], Ta_array, ka_array)
        A[1, 0] = -delta_t/(rho_gas*delta_x*c_gas)*h_r
        A[1, 1] = 1 + delta_t/(rho_a*delta_x*c_a)*k/delta_x + delta_t/(rho_a*delta_x*c_a)*h_r
        A[1, 2] = -delta_t/(rho_a*delta_x*c_a)*k/delta_x
    
    
        for m in range(2,node):
            if m < 2 + nodes_per: #material a
                k = k_value(k_T[m], Ta_array, ka_array)
                A[m, m-1] = A[m, m+1] = -delta_t/(rho_a*delta_x*c_a)*k/delta_x
                A[m, m] = 1 + 2*delta_t/(rho_a*delta_x*c_a)*k/delta_x
            elif m == 2 + nodes_per: #material a & b
                k_a = k_value(k_T[m], Ta_array, ka_array)
                k_b = k_value(k_T[m], Tb_array, kb_array)
                A[m, m-1] = -delta_t/(rho_a*delta_x*c_a)*k_a/delta_x
                A[m, m+1] = -delta_t/(rho_b*delta_x*c_b)*k_b/delta_x
                A[m, m] = 1 + delta_t/(rho_a*delta_x*c_a)*k_a/delta_x + delta_t/(rho_b*delta_x*c_b)*k_b/delta_x
            
            elif 2 + nodes_per <= m <= 2 + 2*nodes_per: #material b
                k = k_value(k_T[m], Tb_array, kb_array)
                A[m, m-1] = A[m, m+1] = -delta_t/(rho_b*delta_x*c_b)*k/delta_x
                A[m, m] = 1 + 2*delta_t/(rho_b*delta_x*c_b)*k/delta_x
            elif m == 2 + 2*nodes_per:
                k_b = k_value(k_T[m], Tb_array, kb_array)
                k_c = k_value(k_T[m], Tc_array, kc_array)
                A[m, m-1] = -delta_t/(rho_b*delta_x*c_b)*k_b/delta_x
                A[m, m+1] = -delta_t/(rho_c*delta_x*c_c)*k_c/delta_x
                A[m, m] = 1 + delta_t/(rho_b*delta_x*c_b)*k_b/delta_x + delta_t/(rho_c*delta_x*c_b)*k_c/delta_x
                
            elif 2 + 2*nodes_per <= m <= 2 + 3*nodes_per: #material c
                k = k_value(k_T[m], Tc_array, kc_array)
                A[m, m-1] = A[m, m+1] = -delta_t/(rho_c*delta_x*c_c)*k/delta_x
                A[m, m] = 1 + 2*delta_t/(rho_c*delta_x*c_c)*k/delta_x
            elif m == 2 + 3*nodes_per:
                k_c = k_value(k_T[m], Tc_array, kc_array)
                k_d = k_value(k_T[m], Td_array, kd_array)
                A[m, m-1] = -delta_t/(rho_c*delta_x*c_c)*k_c/delta_x
                A[m, m+1] = -delta_t/(rho_d*delta_x*c_d)*k_d/delta_x
                A[m, m] = 1 + delta_t/(rho_c*delta_x*c_c)*k_c/delta_x + delta_t/(rho_d*delta_x*c_d)*k_d/delta_x
            elif 2 + 3*nodes_per <= m <= node-2:
                k = k_value(k_T[m], Td_array, kd_array)
                A[m, m-1] = A[m, m+1] = -delta_t/(rho_d*delta_x_d*c_d)*k/delta_x_d
                A[m, m] = 1 + 2*delta_t/(rho_d*delta_x_d*c_d)*k/delta_x_d
            elif m == node-1:
                k = k_value(k_T[m], Td_array, kd_array)
                A[m, m-1] = -delta_t/(rho_d*delta_x_d*c_d)*k/delta_x_d
                A[m, m] = 1 + delta_t/(rho_d*delta_x_d*c_d)*k/delta_x_d
                
        b = T[i-1].copy()
    
        Area = 1
        b[0] = b[0] + delta_t/(rho_gas*delta_x*c_gas)*qgen[i]/Area
        b[1] = b[1] + delta_t/(rho_a*delta_x_d*c_a)*qgen[i]/Area
    
        x = np.linalg.solve(A,b)
    
        
        for j in range(0,node):
            T[i,j] = x[j].copy()
            
        test = k_T - T[i] <= 1e-3
        
        
        
        for w in range(0,len(test)):
            if test[w] == False:
                k_T[w] = (T[i, w] - k_T[w])/2 + k_T[w]
                hr_T[w] = (T[i, w] - hr_T[w])/2 + hr_T[w]
                break
            elif test[-1] == True:
                it = 1
       
        





def transient(t, T):
    #Ax = b
    
    A = np.zeros((node, node), dtype='float64')
    k_T = np.zeros((node))
    hr_T = np.zeros((node))
    

    for i in range(1, len(time)): #go through all time steps
        t = time[i]

        for j in range(0,node): #initial
            k_T[j] = T[i-1,j]
            hr_T[j] = T[i-1,j]
            
        #node inf
        k = k_value(k_T[0], Ta_array, ka_array)
        h_r = eps*sig*(hr_T[1] + hr_T[0])*(hr_T[1]**2 + hr_T[0]**2)
        A[0, 0] = 1 + delta_t/(rho_gas*delta_x*c_gas)*h_r
        A[0, 1] = -delta_t/(rho_gas*delta_x*c_gas)*h_r


        #node 0
        k = k_value(k_T[1], Ta_array, ka_array)
        A[1, 0] = -delta_t/(rho_gas*delta_x*c_gas)*h_r
        A[1, 1] = 1 + delta_t/(rho_a*delta_x*c_a)*k/delta_x + delta_t/(rho_a*delta_x*c_a)*h_r
        A[1, 2] = -delta_t/(rho_a*delta_x*c_a)*k/delta_x
    
    
        for m in range(2,node):
            if m < 2 + nodes_per: #material a
                k = k_value(k_T[m], Ta_array, ka_array)
                A[m, m-1] = A[m, m+1] = -delta_t/(rho_a*delta_x*c_a)*k/delta_x
                A[m, m] = 1 + 2*delta_t/(rho_a*delta_x*c_a)*k/delta_x
            elif m == 2 + nodes_per: #material a & b
                k_a = k_value(k_T[m], Ta_array, ka_array)
                k_b = k_value(k_T[m], Tb_array, kb_array)
                A[m, m-1] = -delta_t/(rho_a*delta_x*c_a)*k_a/delta_x
                A[m, m+1] = -delta_t/(rho_b*delta_x*c_b)*k_b/delta_x
                A[m, m] = 1 + delta_t/(rho_a*delta_x*c_a)*k_a/delta_x + delta_t/(rho_b*delta_x*c_b)*k_b/delta_x
            
            elif 2 + nodes_per <= m <= 2 + 2*nodes_per: #material b
                k = k_value(k_T[m], Tb_array, kb_array)
                A[m, m-1] = A[m, m+1] = -delta_t/(rho_b*delta_x*c_b)*k/delta_x
                A[m, m] = 1 + 2*delta_t/(rho_b*delta_x*c_b)*k/delta_x
            elif m == 2 + 2*nodes_per:
                k_b = k_value(k_T[m], Tb_array, kb_array)
                k_c = k_value(k_T[m], Tc_array, kc_array)
                A[m, m-1] = -delta_t/(rho_b*delta_x*c_b)*k_b/delta_x
                A[m, m+1] = -delta_t/(rho_c*delta_x*c_c)*k_c/delta_x
                A[m, m] = 1 + delta_t/(rho_b*delta_x*c_b)*k_b/delta_x + delta_t/(rho_c*delta_x*c_b)*k_c/delta_x
                
            elif 2 + 2*nodes_per <= m <= 2 + 3*nodes_per: #material c
                k = k_value(k_T[m], Tc_array, kc_array)
                A[m, m-1] = A[m, m+1] = -delta_t/(rho_c*delta_x*c_c)*k/delta_x
                A[m, m] = 1 + 2*delta_t/(rho_c*delta_x*c_c)*k/delta_x
            elif m == 2 + 3*nodes_per:
                k_c = k_value(k_T[m], Tc_array, kc_array)
                k_d = k_value(k_T[m], Td_array, kd_array)
                A[m, m-1] = -delta_t/(rho_c*delta_x*c_c)*k_c/delta_x
                A[m, m+1] = -delta_t/(rho_d*delta_x*c_d)*k_d/delta_x
                A[m, m] = 1 + delta_t/(rho_c*delta_x*c_c)*k_c/delta_x + delta_t/(rho_d*delta_x*c_d)*k_d/delta_x
                
            elif 2 + 3*nodes_per <= m <= node-2:
                k = k_value(k_T[m], Td_array, kd_array)
                A[m, m-1] = A[m, m+1] = -delta_t/(rho_d*delta_x_d*c_d)*k/delta_x_d
                A[m, m] = 1 + 2*delta_t/(rho_d*delta_x_d*c_d)*k/delta_x_d
            elif m == node-1:
                k = k_value(k_T[m], Td_array, kd_array)
                A[m, m-1] = -delta_t/(rho_d*delta_x_d*c_d)*k/delta_x_d
                A[m, m] = 1 + delta_t/(rho_d*delta_x_d*c_d)*k/delta_x_d
                
        b = T[i-1].copy()
        
        Area = 1
        b[0] = b[0] + delta_t/(rho_gas*delta_x*c_gas)*qgen[i]/Area
        b[1] = b[1] + delta_t/(rho_a*delta_x_d*c_a)*qgen[i]/Area
        if what_end == 'A':
            b[-1] = T_int
        elif what_end == 'B':
            pass
        else:
            print("End condition wrong: input 'A' or 'B'")
            sys.exit()
            
            
            
        x = np.linalg.solve(A,b)
        
        
        for j in range(0,node):
            T[i,j] = x[j].copy()
            k_T[j] = x[j].copy()
            hr_T[j] = x[j].copy()
        
        iterate(t, k_T, hr_T)

        
    return A


if __name__ == '__main__':
    sig = 5.67e-8
    t = 0
    print(t)
    #default values?
    what_values = raw_input("Use default values? (Y/N) : ")
    what_values = what_values.capitalize()
    if what_values == 'Y':
        L = 1e-6 #m
        L_d = float(2) #m
        t_end = 20 #sec
        delta_t = 1 #sec        
        T_int = float(35 + 273.15) #K
        nodes_per = float(10)
        nodes_per_d = float(20)
        t_max = 360#sec
        rho_d = 1
        rho_c = 1
        rho_b = 1
        rho_a = 1
        rho_gas = 1
        c_d = 1
        c_c = 1
        c_b = 1
        c_a = 1
        c_gas = 1
        eps = 1
        time = []
        while t < t_max:
            time.append(t)
            t = t + delta_t
    else:
        L = float(input("The thickness of the layers: "))
        L_d = float(input("The thickness of the base layer: "))
        t_max = input("Run until t (sec) = ")
        T_int = float(input("Initial temperature (K): "))

    what_q = raw_input("Use default qgen? (Y/N) : ")
    what_q = what_q.capitalize()
    if what_q == 'Y':
        qgen = q_gen('Y')
    else:
        qgen = q_gen('n')
        
    delta_x = L/(nodes_per-1) 
    delta_x_d = L_d/(nodes_per_d-1)
    node = 1 #node in gas
    x = 0
    while x <= L:
        node = node + 1
        x = x + delta_x
    node = node*4 
    
    x = 0
    
    while x < L_d:
        node = node + 1
        x = x + delta_x_d
    node = node + 1
    x_d = L_d - x
        
        
    
    #retrieve k values
    what_k = raw_input("Use default materials? (Y/N) : ")
    what_k = what_k.capitalize()
    if what_k == 'Y':
        mat_a = 'Uranium.txt'
        mat_b = 'Platinum.txt'
        mat_c = 'Titanium.txt'
        mat_d = 'Alumina.txt'
        Ta_array, ka_array = k_array(mat_b)
        Tb_array, kb_array = k_array(mat_b)
        Tc_array, kc_array = k_array(mat_c)
        Td_array, kd_array = k_array(mat_b)
    else:
        print("code later")
    
    what_end = raw_input("Use (A) Constant end temperature or (B) Insulated end condition: ")
    what_end = what_end.capitalize()
    
    
    #Initial Conditions
    T = np.zeros((t_max/delta_t, node), dtype='float64')
    for i in range(0, node):
        T[0,i] = T_int
    
    
   
    A = transient(t, T)
        