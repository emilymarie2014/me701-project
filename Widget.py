#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 17:37:48 2018

@author: emilystallbaumer
"""

import sys
import platform
import numpy as np
import re
#Order of imports is important

from PyQt5.QtWidgets import (QMainWindow, QApplication, QWidget, QLineEdit, 
                             QVBoxLayout, QHBoxLayout, QAction, QMessageBox,QFileDialog,
                             QSizePolicy, QComboBox, QSplitter, QPushButton, QLabel)
from PyQt5.QtCore import QT_VERSION_STR, PYQT_VERSION_STR, Qt
from PyQt5.QtGui import QIcon, QPixmap

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure


class MainWindow(QMainWindow) :
    
    def __init__(self, parent=None) :
        super(MainWindow, self).__init__(parent)

        ########################################################################
        # ADD MENU ITEMS
        ########################################################################
        print('w')
        # Create the File menu
        self.menuFile = self.menuBar().addMenu("&File")
        self.actionSaveAs = QAction("&Save As", self)
        #self.actionSaveAs.triggered.connect(self.saveas)
        self.actionQuit = QAction("&Quit", self)
        self.actionQuit.triggered.connect(self.close)
        self.menuFile.addActions([self.actionSaveAs, self.actionQuit])
        
        # Create the Help menu
        #self.menuHelp = self.menuBar().addMenu("&Help")
        #self.actionAbout = QAction("&About",self)
        #self.actionAbout.triggered.connect(self.about)
        #self.menuHelp.addActions([self.actionAbout])
        




        #create widgets
        self.widget = QWidget()
        #self.plot = MatplotlibCanvas()
        self.box1 = QComboBox(self)
        self.box1.addItem('Uranium')
        self.box1.addItem('Platinum')
        self.box1.addItem('Titanium')
        self.box1.addItem('Alumina')
        self.box1.addItem('')
        self.box1.setEditable(True)
        self.box1.setInsertPolicy(1)
        
        self.box2 = QComboBox(self)
        self.box2.addItem('Platinum')
        self.box2.addItem('Titanium')
        self.box2.addItem('Alumina')
        self.box2.addItem('Uranium')
        self.box2.addItem('')
        self.box2.setEditable(True)
        self.box2.setInsertPolicy(1)

        self.box3 = QComboBox(self)
        self.box3.addItem('Titanium')
        self.box3.addItem('Alumina')
        self.box3.addItem('Uranium')
        self.box3.addItem('Platinum')
        self.box3.addItem('')
        self.box3.setEditable(True)
        self.box3.setInsertPolicy(1)

        self.box4 = QComboBox(self)
        self.box4.addItem('Alumina')
        self.box4.addItem('Uranium')
        self.box4.addItem('Platinum')
        self.box4.addItem('Titanium')
        self.box4.addItem('')
        self.box4.setEditable(True)
        self.box4.setInsertPolicy(1)
        
        self.box5 = QComboBox(self)
        self.box5.addItem('Argon')
        self.box5.addItem('')
        self.box5.setEditable(True)
        self.box5.setInsertPolicy(1)
        
        
        vlayout = QVBoxLayout()
        self.l1 = QLabel()
        self.l2 = QLabel()
        self.l3 = QLabel()
        self.l4 = QLabel()
        self.l5 = QLabel()
        
        self.l1.setText("Material A")
        self.l2.setText("Material B")
        self.l3.setText("Material C")
        self.l4.setText("Material D")
        self.l5.setText("Gas")
        
        
        #vlayout.addWidget(self.plot)
        vlayout.addWidget(self.l1)
        vlayout.addWidget(self.box1)
        vlayout.addWidget(self.l2)
        vlayout.addWidget(self.box2) 
        vlayout.addWidget(self.l3)
        vlayout.addWidget(self.box3)
        vlayout.addWidget(self.l4)
        vlayout.addWidget(self.box4)
        vlayout.addWidget(self.l5)
        vlayout.addWidget(self.box5)
        
        vlayout2 = QVBoxLayout()
        
        self.L_a = QLineEdit('30e-9') #m
        self.L_b = QLineEdit('50e-9') #m
        self.L_c = QLineEdit('5e-9')  #m
        self.L_d = QLineEdit('1e-2')  #m
        self.eps = QLineEdit('1.0')   #m
        
        self.lLa = QLabel()
        self.lLb = QLabel()
        self.lLc = QLabel()
        self.lLd = QLabel()
        self.leps = QLabel()
        
        self.lLa.setText("Length of A (m)")
        self.lLb.setText("Length of B (m)")
        self.lLc.setText("Length of C (m)")
        self.lLd.setText("Length of D (m)")  
        self.leps.setText("Assumed Emmisivity")
        
        vlayout2.addWidget(self.lLa)
        vlayout2.addWidget(self.L_a)
        vlayout2.addWidget(self.lLb)
        vlayout2.addWidget(self.L_b) 
        vlayout2.addWidget(self.lLc)
        vlayout2.addWidget(self.L_c)
        vlayout2.addWidget(self.lLd)
        vlayout2.addWidget(self.L_d)
        vlayout2.addWidget(self.leps)
        vlayout2.addWidget(self.eps)
        
        vlayout3 = QVBoxLayout()
        
        self.delta_t = QLineEdit("1.0") #sec        
        self.T_int = QLineEdit("308.15") #K
        self.nodes_per = QLineEdit("10")
        self.nodes_per_d = QLineEdit("20")
        self.t_max = QLineEdit("360.0")#sec
        
        self.ldt = QLabel()
        self.lTint = QLabel()  
        self.lnp = QLabel() 
        self.lnpd = QLabel() 
        self.ltmax = QLabel()  
        
        self.ldt.setText("time step")
        self.lTint.setText("initial temp")
        self.lnp.setText("# of nodes per materials A, B, C")
        self.lnpd.setText("# of nodes per material D")
        self.ltmax.setText("max time")
        
        
        vlayout3.addWidget(self.ldt)
        vlayout3.addWidget(self.delta_t)
        vlayout3.addWidget(self.ltmax)
        vlayout3.addWidget(self.t_max)
        vlayout3.addWidget(self.lTint)
        vlayout3.addWidget(self.T_int)
        vlayout3.addWidget(self.lnp)
        vlayout3.addWidget(self.nodes_per)
        vlayout3.addWidget(self.lnpd)
        vlayout3.addWidget(self.nodes_per_d)

        
        self.button1 = QPushButton('Run', self)
        #self.button2 = QPushButton('Clear', self)
#        
#        hlayout = QHBoxLayout()
#        hlayout.addWidget(self.button1)
#        hlayout.addWidget(self.button2)

        
        
        hlayout = QHBoxLayout()
#        
        hlayout.addLayout(vlayout)
        hlayout.addLayout(vlayout2)
        hlayout.addLayout(vlayout3)
#        
        vlayout4 = QVBoxLayout()
        
        hlayout2 = QHBoxLayout()
        
        hlayout2.addWidget(self.button1)
        #hlayout2.addWidget(self.button2)

        vlayout4.addLayout(hlayout)
        vlayout4.addLayout(hlayout2)
        
        
        self.widget.setLayout(vlayout4)
         
        self.setCentralWidget(self.widget)
        self.show()
        
        self.button1.clicked.connect(self.run)
        #self.button2.clicked.connect(self.clear)
#        self.box2.returnPressed.connect(self.run)
        

        
        
        
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
            b[0] = b[0] + delta_t/(rho_gas*delta_x*c_gas)*qgen/Area
            b[1] = b[1] + delta_t/(rho_a*delta_x_d*c_a)*qgen/Area
        
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
        
    
        for i in range(1, t_max/delta_t): #go through all time steps
            t = t + delta_t
    
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
            b[0] = b[0] + delta_t/(rho_gas*delta_x*c_gas)*qgen/Area
            b[1] = b[1] + delta_t/(rho_a*delta_x_d*c_a)*qgen/Area
    
            x = np.linalg.solve(A,b)
            
            
            for j in range(0,node):
                T[i,j] = x[j].copy()
            
            iterate(k_T, hr_T)
    
            
        return A
       
        
        
    def run(self):
        sig = 5.67e-8
        t = 0
        print(t)
        #default values?
        nodes_per = float(eval(self.lnp.text()))
        L_a = .000004
        L_b = .0002
        L_c = .000003124
        L_d = .1
        nodes_per_d = float(eval(self.lnod.text()))
        print(type(L_a))
        delta_x_a = L_a/(nodes_per-1) 
        delta_x_b = L_b/(nodes_per-1)
        delta_x_c = L_c/(nodes_per-1)
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


        mat_a = box1 + '.txt'
        mat_b = box2 + '.txt'
        mat_c = box3 + '.txt'
        mat_d = box4 + '.txt'
        Ta_array, ka_array = k_array(mat_b)
        Tb_array, kb_array = k_array(mat_b)
        Tc_array, kc_array = k_array(mat_c)
        Td_array, kd_array = k_array(mat_b)
        
        
        #Initial Conditions
        T = np.zeros((t_max/delta_t, node), dtype='float64')
        for i in range(0, node):
            T[0,i] = T_int
        
        
       
        A = transient(t, T)
        
        #self.update()
app = QApplication(sys.argv)
form = MainWindow()
form.show()
app.exec_()
        