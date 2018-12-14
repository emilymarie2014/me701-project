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
import os
import matplotlib.pyplot as plt
#Order of imports is important

from PyQt5.QtWidgets import (QMainWindow, QApplication, QWidget, QLineEdit, 
                             QVBoxLayout, QHBoxLayout, QAction, QMessageBox,QFileDialog,
                             QSizePolicy, QComboBox, QSplitter, QPushButton, QLabel)
from PyQt5.QtCore import QT_VERSION_STR, PYQT_VERSION_STR, Qt
from PyQt5.QtGui import QIcon, QPixmap

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

class ThirdWindow(QMainWindow): #results for steady state
    def __init__(self, x_array, T, h, melt):
        super(ThirdWindow, self).__init__()
        self.titleS = 'Steady State Results'
        self.setMinimumSize(600,500)
        self.menuFile = self.menuBar().addMenu("&File")
        self.actionSaveAs_s = QAction("&Steady State Save As", self)
        self.actionSaveAs_s.triggered.connect(self.saveas)
        self.actionQuit = QAction("&Quit", self)
        self.actionQuit.triggered.connect(self.close)
        self.menuFile.addActions([self.actionSaveAs_s, self.actionQuit])
        self.setMinimumSize(500,500)
        #self.setGeometry(300,300,700,700)
        self.widg = QWidget()
        self.p = MatplotlibCanvas()
        x = x_array
        y = T
        
        
        vlayoutp = QVBoxLayout()
        vlayoutp.addWidget(self.p)

        vlayout = QVBoxLayout()

        self.t = QLabel()
        self.t.setText('Hotest Temperature in each material')
        self.hA = QLabel()
        self.hA.setText('Material A: %0.2f K' % h[0])
        self.hB = QLabel()
        self.hB.setText('Material B: %0.2f K' % h[1])
        self.hC = QLabel()
        self.hC.setText('Material C: %0.2f K' % h[2])
        self.hD = QLabel()
        self.hD.setText('Material D: %0.2f K' % h[3])
        vlayout.addWidget(self.t)
        vlayout.addWidget(self.hA)
        vlayout.addWidget(self.hB)
        vlayout.addWidget(self.hC)
        vlayout.addWidget(self.hD)


        self.tt = QLabel()
        self.tt.setText("Does the material melt?")
        
        difA = float(melt[0]) - float(h[0])        
        self.meltA = QLabel()
        if difA > 0:
            self.meltA.setText('No, Material A reaches %0.2f K below the melting point' % difA)
        else: 
            difA = difA*-1
            self.meltA.setText('WARNING!!!! Material A reaches %0.2f ABOVE the melting point' % difA)
        difB = float(melt[1]) - float(h[1])  
        self.meltB = QLabel()
        if difB > 0:
            self.meltB.setText('No, Material B reaches %0.2f K below the melting point' % difB)
        else: 
            difB = difB*-1
            self.meltB.setText('WARNING!!!! Material B reaches %0.2f ABOVE the melting point' % difB)
        difC = float(melt[2]) - float(h[2])
        self.meltC = QLabel()
        if difC > 0:
            self.meltC.setText('No, Material C reaches %0.2f K below the melting point' % difC)
        else: 
            difC =difC*-1
            self.meltC.setText('WARNING!!!! Material C reaches %0.2f ABOVE the melting point' % difC)
        difD = float(melt[3]) - float(h[3])
        self.meltD = QLabel()
        if difD > 0:
            self.meltD.setText('No, Material D reaches %0.2f K below the melting point' % difD)
        else: 
            difD = difD*-1
            self.meltD.setText('WARNING!!!! Material D reaches %0.2f ABOVE the melting point' % difD)
        
        vlayout2 = QVBoxLayout()
        vlayout2.addWidget(self.tt)
        vlayout2.addWidget(self.meltA)
        vlayout2.addWidget(self.meltB)
        vlayout2.addWidget(self.meltC)
        vlayout2.addWidget(self.meltD)
        
        hlayout = QHBoxLayout()
        hlayout.addLayout(vlayout)
        hlayout.addLayout(vlayout2)
        
        vlayoutp.addLayout(hlayout)
        self.widg.setLayout(vlayoutp)
        self.titleS = 'Steady State Results'
        self.setWindowTitle(self.titleS)
        self.setCentralWidget(self.widg)
        
        

        self.widg.setLayout(vlayout)
        self.setCentralWidget(self.widg)
        self.show()
        
        self.p = self.p.ddraw(x, y, 0, 0, 1)
        
    def saveas(self):
        """Save the steady state Temperature values"""
        fname = QFileDialog.getSaveFileName(self,'save steady state file')[0]
        os.rename("temporary_steady_state.txt", fname)

class SecondWindow(QMainWindow): #results for transient
    def __init__(self, x_array, T, h, melt):
        super(SecondWindow, self).__init__()
        self.titleT = 'Transient Results'
        self.setMinimumSize(600,500)
        self.menuFile = self.menuBar().addMenu("&File")
        self.actionSaveAs_t = QAction("&Transient Save As", self)
        self.actionSaveAs_t.triggered.connect(self.saveas)
        self.actionQuit = QAction("&Quit", self)
        self.actionQuit.triggered.connect(self.close)
        self.menuFile.addActions([self.actionSaveAs_t, self.actionQuit])
        #self.setGeometry(300,300,700,700)
        self.widg = QWidget()
        self.p = MatplotlibCanvas()
        x = x_array[1:]
        y = T[0, 1:]
        z = T[-1, 1:]
        
        vlayoutp = QVBoxLayout()
        vlayoutp.addWidget(self.p)

        vlayout = QVBoxLayout()
        #TODO: find the hotest of each material
        self.t = QLabel()
        self.t.setText('Hotest Temperature in each material')
        self.hA = QLabel()
        self.hA.setText('Material A: %0.2f K' % h[0])
        self.hB = QLabel()
        self.hB.setText('Material B: %0.2f K' % h[1])
        self.hC = QLabel()
        self.hC.setText('Material C: %0.2f K' % h[2])
        self.hD = QLabel()
        self.hD.setText('Material D: %0.2f K' % h[3])
        vlayout.addWidget(self.t)
        vlayout.addWidget(self.hA)
        vlayout.addWidget(self.hB)
        vlayout.addWidget(self.hC)
        vlayout.addWidget(self.hD)
        

       
        self.tt = QLabel()
        self.tt.setText("Does the material melt?")
        
        difA = float(melt[0]) - float(h[0])        
        self.meltA = QLabel()
        if difA > 0:
            self.meltA.setText('No, Material A reaches %0.2f K below the melting point' % difA)
        else: 
            difA = difA*-1
            self.meltA.setText('WARNING!!!! Material A reaches %0.2f ABOVE the melting point' % difA)
        difB = float(melt[1]) - float(h[1])  
        self.meltB = QLabel()
        if difB > 0:
            self.meltB.setText('No, Material B reaches %0.2f K below the melting point' % difB)
        else: 
            difB = difB*-1
            self.meltB.setText('WARNING!!!! Material B reaches %0.2f ABOVE the melting point' % difB)
        difC = float(melt[2]) - float(h[2])
        self.meltC = QLabel()
        if difC > 0:
            self.meltC.setText('No, Material C reaches %0.2f K below the melting point' % difC)
        else: 
            difC =difC*-1
            self.meltC.setText('WARNING!!!! Material C reaches %0.2f ABOVE the melting point' % difC)
        difD = float(melt[3]) - float(h[3])
        self.meltD = QLabel()
        if difD > 0:
            self.meltD.setText('No, Material D reaches %0.2f K below the melting point' % difD)
        else: 
            difD = difD*-1
            self.meltD.setText('WARNING!!!! Material D reaches %0.2f ABOVE the melting point' % difD)
        
        vlayout2 = QVBoxLayout()
        vlayout2.addWidget(self.tt)
        vlayout2.addWidget(self.meltA)
        vlayout2.addWidget(self.meltB)
        vlayout2.addWidget(self.meltC)
        vlayout2.addWidget(self.meltD)
        
        hlayout = QHBoxLayout()
        hlayout.addLayout(vlayout)
        hlayout.addLayout(vlayout2)
        
        vlayoutp.addLayout(hlayout)
        self.widg.setLayout(vlayoutp)
        self.setWindowTitle(self.titleT)
        self.setCentralWidget(self.widg)
        
        
        self.show()
        self.p = self.p.ddraw(x, y, x, z, 2)
        #self.close()
        
    def saveas(self):
        """Save the transient Temperature values"""
        fname = QFileDialog.getSaveFileName(self,'save transient file')[0]
        os.rename("temporary_transient.txt", fname)
        
class MainWindow(QMainWindow) : #the GUI
    
    def __init__(self, parent=None) :
        """create the main window"""
        super(MainWindow, self).__init__(parent)
        self.TPlot = None
        self.SPlot = None
        ########################################################################
        # ADD MENU ITEMS
        ########################################################################
        # Create the File menu
        self.menuFile = self.menuBar().addMenu("&File")
        self.actionSaveAs_t = QAction("&Transient Save As", self)
        self.actionSaveAs_t.triggered.connect(self.saveas_t)
        self.actionSaveAs_ss = QAction("&Steady State Save As", self)
        self.actionSaveAs_ss.triggered.connect(self.saveas_ss)
        self.actionQuit = QAction("&Quit", self)
        self.actionQuit.triggered.connect(self.close)
        self.menuFile.addActions([self.actionSaveAs_t, self.actionSaveAs_ss, self.actionQuit])
        
        # Create the Help menu
        #self.menuHelp = self.menuBar().addMenu("&Help")
        #self.actionAbout = QAction("&About",self)
        #self.actionAbout.triggered.connect(self.about)
        #self.menuHelp.addActions([self.actionAbout])
        

        #create widgets
        self.widget = QWidget()
        self.plot = MatplotlibCanvas()
        self.box1 = QComboBox(self)
        self.box1.addItem('UO2')
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
        self.box2.addItem('UO2')
        self.box2.addItem('')
        self.box2.setEditable(True)
        self.box2.setInsertPolicy(1)

        self.box3 = QComboBox(self)
        self.box3.addItem('Titanium')
        self.box3.addItem('Alumina')
        self.box3.addItem('UO2')
        self.box3.addItem('Platinum')
        self.box3.addItem('')
        self.box3.setEditable(True)
        self.box3.setInsertPolicy(1)

        self.box4 = QComboBox(self)
        self.box4.addItem('Alumina')
        self.box4.addItem('UO2')
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
        
        self.La = QLineEdit('30e-9') #m
        self.Lb = QLineEdit('50e-9') #m
        self.Lc = QLineEdit('5e-9')  #m
        self.Ld = QLineEdit('1e-2')  #m
        self.ep = QLineEdit('1.0')   #m
        
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
        vlayout2.addWidget(self.La)
        vlayout2.addWidget(self.lLb)
        vlayout2.addWidget(self.Lb) 
        vlayout2.addWidget(self.lLc)
        vlayout2.addWidget(self.Lc)
        vlayout2.addWidget(self.lLd)
        vlayout2.addWidget(self.Ld)
        vlayout2.addWidget(self.leps)
        vlayout2.addWidget(self.ep)
        
        
        
        vlayout3 = QVBoxLayout()
        
        self.qto = QLineEdit('1')
        self.qt1 = QLineEdit('1.015')
        self.qgen_eq = QLineEdit('2362.5') 
        self.cycle = QLineEdit('1')
        self.tolerance = QLineEdit('1e-2')
        
        self.lqto = QLabel()
        self.lqt1 = QLabel()
        self.lqeq = QLabel()
        self.lcycle = QLabel()
        self.ltol = QLabel()
        
        self.lqto.setText("qgen starts (sec)")
        self.lqt1.setText("qgen ends (sec)")
        self.lqeq.setText("qgen (W/m^3)")
        self.lcycle.setText("Number of qgen cylces")
        self.ltol.setText("Tolerance")
        
        vlayout3.addWidget(self.lqto)
        vlayout3.addWidget(self.qto)
        vlayout3.addWidget(self.lqt1)
        vlayout3.addWidget(self.qt1)
        vlayout3.addWidget(self.lqeq)
        vlayout3.addWidget(self.qgen_eq)
        vlayout3.addWidget(self.lcycle)
        vlayout3.addWidget(self.cycle)
        vlayout3.addWidget(self.ltol)
        vlayout3.addWidget(self.tolerance) 
        
        
        vlayout4 = QVBoxLayout()  
        
        self.np = QLineEdit("10")
        self.npd = QLineEdit("20")
        self.endCon = QComboBox(self)
        self.endCon.addItem("Constant Temp")
        self.endCon.addItem("Insulated End")
        self.runCon = QComboBox(self)
        self.runCon.addItem("Transient")
        self.runCon.addItem("Steady State")
        self.runCon.addItem("Both")
        self.are = QLineEdit("4.90625e-14")
        
        self.lnp = QLabel() 
        self.lnpd = QLabel()
        self.lec = QLabel()
        self.lrc = QLabel()
        self.lare = QLabel()
        
        self.lnp.setText("# of nodes per materials A, B, C")
        self.lnpd.setText("# of nodes per material D")
        self.lec.setText("Condition at x_d")
        self.lrc.setText("Run as:")
        self.lare.setText("Cross section area, material A (m^2)")
        
        vlayout4.addWidget(self.lnp)
        vlayout4.addWidget(self.np)
        vlayout4.addWidget(self.lnpd)
        vlayout4.addWidget(self.npd)
        vlayout4.addWidget(self.lec)
        vlayout4.addWidget(self.endCon)
        vlayout4.addWidget(self.lrc)
        vlayout4.addWidget(self.runCon)
        vlayout4.addWidget(self.lare)
        vlayout4.addWidget(self.are)
        
        
        
        vlayout5 = QVBoxLayout()
        
        self.tra = QLabel()
        self.tra.setText("TRANSIENT VARIABLES")
        
        self.Tint = QLineEdit("308.15") #K
        self.dt = QLineEdit("0.005") #sec   
        
        self.lTint = QLabel() 
        self.ldt = QLabel()
        
        self.lTint.setText("initial temp (k)")
        self.ldt.setText("time step (sec)")
        
        vlayout5.addWidget(self.tra)
        vlayout5.addWidget(self.lTint)
        vlayout5.addWidget(self.Tint)
        vlayout5.addWidget(self.ldt)
        vlayout5.addWidget(self.dt)
        
        
        self.ste = QLabel()
        self.ste.setText("STEADY STATE VARIABLES")
        
        self.Tgas = QLineEdit('2000') #K
        self.Tend = QLineEdit('308.15') #K
        
        self.lTg = QLabel()
        self.lTe = QLabel()
        
        self.lTg.setText('temperature of gas (k)')
        self.lTe.setText('temperature of Constant Temp Boundry (k)')
        vlayout5.addWidget(self.ste)
        vlayout5.addWidget(self.lTg)
        vlayout5.addWidget(self.Tgas)
        vlayout5.addWidget(self.lTe)
        vlayout5.addWidget(self.Tend)
        
        #self.tmax = QLineEdit("360.0")#sec
        #self.ltmax = QLabel()
        #self.ltmax.setText("max time")
#        vlayout5.addWidget(self.ltmax)
#        vlayout5.addWidget(self.tmax)
        
        hlayout = QHBoxLayout() 
        hlayout.addLayout(vlayout)
        hlayout.addLayout(vlayout2)
        hlayout.addLayout(vlayout3)
        hlayout.addLayout(vlayout4) 
        hlayout.addLayout(vlayout5)
        
        
        self.button1 = QPushButton('Run', self)
        self.button2 = QPushButton('Reset', self)
        self.button3 = QPushButton('Clear', self)
                                   
        vlayout6 = QVBoxLayout()
        
        hlayout2 = QHBoxLayout()        
        hlayout2.addWidget(self.button1)
        hlayout2.addWidget(self.button2)
        hlayout2.addWidget(self.button3)
        vlayout6.addLayout(hlayout)
        vlayout6.addLayout(hlayout2)
        
        vlayoutplot = QVBoxLayout()
        vlayoutplot.addWidget(self.plot)
        vlayoutplot.addLayout(vlayout6)
        
        
        self.widget.setLayout(vlayoutplot)
         
        self.setCentralWidget(self.widget)
        self.show()
        
        self.button1.clicked.connect(self.run)
        self.button2.clicked.connect(self.reset)
        self.button3.clicked.connect(self.clear)
        self.name = QLineEdit()
        

        self.title = "Set Variables"
        self.setWindowTitle(self.title)
        

    def Hot_transient(self, text):
        """determines what the hottest temperature is in each material and how close 
        to the melting point the temperature is, for the transient problem"""
        a = 1 + self.nodes_per 
        b = a + self.nodes_per
        c = b + self.nodes_per
        d = c + self.nodes_per_d
        self.hot = [0, 0, 0, 0]
        for i in range(len(self.T)):
            for x in range(1,a):
                if self.T[i, x] > self.hot[0]:
                    self.hot[0] = self.T[i,x]
            for x in range(a-1, b):
                if self.T[i, x] > self.hot[1]:
                    self.hot[1] = self.T[i,x]
            for x in range(b-1, c):
                if self.T[i,x] > self.hot[2]:
                    self.hot[2] = self.T[i,x]
            for x in range(c-1, d):
                if self.T[i,x] > self.hot[3]:
                    self.hot[3] = self.T[i,x]
        
        return self.hot
    
    def Hot_steady(self, text):
        """determines what the hottest temperature is in each material and how close to 
        the melting point the temperature is, for the steady state problem"""
        a = self.nodes_per 
        b = a + self.nodes_per
        c = b + self.nodes_per
        d = c + self.nodes_per_d
        self.hot = [0, 0, 0, 0]
        for i in range(0,a):
                if self.T_ss[i] > self.hot[0]:
                    self.hot[0] = self.T_ss[i]
        for i in range(a-1, b):
                if self.T_ss[i] > self.hot[1]:
                    self.hot[1] = self.T_ss[i]
        for i in range(b-1, c):
                if self.T_ss[i] > self.hot[2]:
                    self.hot[2] = self.T_ss[i]
        for i in range(c-1, d):
                if self.T_ss[i] > self.hot[3]:
                    self.hot[3] = self.T_ss[i]
        
        return self.hot    
    
    def clear(self):
        """Clear the values and results as well as the plot"""
        self.La.setText('')
        self.Lb.setText('')
        self.Lc.setText('')
        self.Ld.setText('')
        self.ep.setText('')
        self.Tint.setText('')
        self.dt.setText('')
#        self.tmax.setText('')
        self.np.setText('')
        self.npd.setText('')
        self.qgen_eq.setText('')
        self.tolerance.setText('')
        self.qto.setText('')
        self.qt1.setText('')
        self.cycle.setText('')
        self.Tgas.setText('')
        self.Tend.setText('')
        self.are.setText('')
        
        x = [0,0]
        y = [0,0]
        self.plot.redraw(x,y,0)
        if self.TPlot != None:
            self.TPlot.close()
        if self.SPlot != None:
            self.SPlot.close()
            
    def reset(self):
        """resets the values and results as well as the plot"""
        self.La.setText('30e-9')
        self.Lb.setText('50e-9')
        self.Lc.setText('5e-9')
        self.Ld.setText('1e-2')
        self.ep.setText('1.0')
        self.Tint.setText('308.15')
        self.dt.setText('0.005')
#        self.tmax.setText('')
        self.np.setText('10')
        self.npd.setText('20')
        self.qgen_eq.setText('2362.5')
        self.tolerance.setText('1e-2')
        self.qto.setText('1')
        self.qt1.setText('1.015')
        self.cycle.setText('1')
        self.Tgas.setText('2000')
        self.Tend.setText('308.15')
        self.are.setText('4.90625e-14')
        
        x = [0,0]
        y = [0,0]
        self.plot.redraw(x,y, 0)
        if self.TPlot != None:
            self.TPlot.close()
        if self.SPlot != None:
            self.SPlot.close()


    def saveas_t(self):
        """Save the transient Temperature values"""
        fname = QFileDialog.getSaveFileName(self,'save transient file')[0]
        os.rename("temporary_transient.txt", fname)
    def saveas_ss(self):
        """Save the steady state Temperature values"""
        fname = QFileDialog.getSaveFileName(self,'save steady state file')[0]
        os.rename("temporary_steady_state.txt", fname)
    
    
    def closeEvent(self, event):
        """removes the temporary transient and steady state files, if they exist, as 
        well as closes the transient and steady state result windows, if they are open, 
        when the main window is closed"""
        if os.path.exists("temporary_transient.txt"):
            os.remove("temporary_transient.txt")
        if os.path.exists("temporary_steady_state.txt"):
            os.remove("temporary_steady_state.txt")
        if self.TPlot != None:
            self.TPlot.close()
        if self.SPlot != None:
            self.SPlot.close()
    
    
    def prop(self, p):
        """requires the material and reads the “Properties.txt” file and returns 
        the appropriate density, specific heat, and melting temperature for the material
        by a regular expression"""
        f = open('Properties.txt')
        s = f.readlines()
        f.close()
        
        for i in range(0,len(s)):
            match = re.search(p, s[i])
            if match != None:
                spl = s[i].split()
                rho = float(spl[1])
                c = float(spl[2])
                if len(spl) == 4:
                    melt = float(spl[3])
                    return rho, c, melt
                else:
                    return rho, c
     
        
    def k_array(self, a):
        """requires the material and reads the appropriate material file and returns 
        an array of the temperature and an array of the corresponding thermal conductivities"""
#        print('k_array')
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
    
    def k_value(self, T, T_array, k_array):
        """requires the desired temperature, the temperature array, created in “k_array,” and the thermal 
        conductivity array, created in “k_array.” The thermal conductivity for the temperature passed in is 
        solved for and returned"""
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
    
    def iterate(self, t, i, k_T, hr_T, qgen, A):
        """takes in the time, the time step, the array of temperatures to be used to determine the thermal conductivity,
        the array of temperatures to be used to determine the heat transfer  coefficient for radiation, the heat 
        generation, and matrix ‘A.’ It then solves for the new temperatures and compares them to the temperatures 
        used to solve for the thermal conductivity and the radiation heat transfer coefficient, until they are within 
        the tolerance and writes these values to a temporary file"""
        it = 0
        g = 0
        if g > 10000:
            print('*********BREAK: number of iterarations is excessive*********** \n')
            print('error is: ', k_T - self.T[i])
            return
        while it == 0:
            g = g + 1  
            h_r = self.eps*self.sig*(hr_T[1] + hr_T[0])*(hr_T[1]**2 + hr_T[0]**2) #get hr value
            #create array A
            A[0, 0] = 1 + self.delta_t/(self.rho_gas*self.delta_x_a*self.c_gas)*h_r
            A[0, 1] = -self.delta_t/(self.rho_gas*self.delta_x_a*self.c_gas)*h_r
            
    
            #node 0
            k = self.k_value(k_T[1], self.Ta_array, self.ka_array)
            A[1, 0] = -self.delta_t/(self.rho_gas*self.delta_x_a*self.c_gas)*h_r
            A[1, 1] = 1 + self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*k/self.delta_x_a + self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*h_r
            A[1, 2] = -self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*k/self.delta_x_a
        
        
            for m in range(2,self.node):
                if m < 2 + self.nodes_per: #material a
                    k = self.k_value(k_T[m], self.Ta_array, self.ka_array)
                    A[m, m-1] = A[m, m+1] = -self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*k/self.delta_x_a
                    A[m, m] = 1 + 2*self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*k/self.delta_x_a
                elif m == 2 + self.nodes_per: #material a & b
                    k_a = self.k_value(k_T[m], self.Ta_array, self.ka_array)
                    k_b = self.k_value(k_T[m], self.Tb_array, self.kb_array)
                    A[m, m-1] = -self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*k_a/self.delta_x_a
                    A[m, m+1] = -self.delta_t/(self.rho_b*self.delta_x_b*self.c_b)*k_b/self.delta_x_b
                    A[m, m] = 1 + self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*k_a/self.delta_x_a + self.delta_t/(self.rho_b*self.delta_x_b*self.c_b)*k_b/self.delta_x_b
                
                elif 2 + self.nodes_per <= m <= 2 + 2*self.nodes_per: #material b
                    k = self.k_value(k_T[m], self.Tb_array, self.kb_array)
                    A[m, m-1] = A[m, m+1] = -self.delta_t/(self.rho_b*self.delta_x_b*self.c_b)*k/self.delta_x_b
                    A[m, m] = 1 + 2*self.delta_t/(self.rho_b*self.delta_x_b*self.c_b)*k/self.delta_x_b
                elif m == 2 + 2*self.nodes_per:
                    k_b = self.k_value(k_T[m], self.Tb_array, self.kb_array)
                    k_c = self.k_value(k_T[m], self.Tc_array, self.kc_array)
                    A[m, m-1] = -self.delta_t/(self.rho_b*self.delta_x_b*self.c_b)*k_b/self.delta_x_b
                    A[m, m+1] = -self.delta_t/(self.rho_c*self.delta_x_c*self.c_c)*k_c/self.delta_x_c
                    A[m, m] = 1 + self.delta_t/(self.rho_b*self.delta_x_b*self.c_b)*k_b/self.delta_x_b + self.delta_t/(self.rho_c*self.delta_x_c*self.c_c)*k_c/self.delta_x_c
                    
                elif 2 + 2*self.nodes_per <= m <= 2 + 3*self.nodes_per: #material c
                    k = self.k_value(k_T[m], self.Tc_array, self.kc_array)
                    A[m, m-1] = A[m, m+1] = -self.delta_t/(self.rho_c*self.delta_x_c*self.c_c)*k/self.delta_x_c
                    A[m, m] = 1 + 2*self.delta_t/(self.rho_c*self.delta_x_c*self.c_c)*k/self.delta_x_c
                elif m == 2 + 3*self.nodes_per: #material c&d
                    k_c = self.k_value(k_T[m], self.Tc_array, self.kc_array)
                    k_d = self.k_value(k_T[m], self.Td_array, self.kd_array)
                    A[m, m-1] = -self.delta_t/(self.rho_c*self.delta_x_c*self.c_c)*k_c/self.delta_x_c
                    A[m, m+1] = -self.delta_t/(self.rho_d*self.delta_x_d*self.c_d)*k_d/self.delta_x_d
                    A[m, m] = 1 + self.delta_t/(self.rho_c*self.delta_x*self.c_c)*k_c/self.delta_x_c + self.delta_t/(self.rho_d*self.delta_x*self.c_d)*k_d/self.delta_x
                elif 2 + 3*self.nodes_per <= m <= self.node-2: #material d
                    k = self.k_value(k_T[m], self.Td_array, self.kd_array)
                    A[m, m-1] = A[m, m+1] = -self.delta_t/(self.rho_d*self.delta_x_d*self.c_d)*k/self.delta_x_d
                    A[m, m] = 1 + 2*self.delta_t/(self.rho_d*self.delta_x_d*self.c_d)*k/self.delta_x_d
                elif m == self.node-1:
                    if self.endCon.currentText() == "Constant Temp":
                        A[m, m] = 1
                    else:
                        k = self.k_value(k_T[m], self.Td_array, self.kd_array)
                        A[m, m-1] = -self.delta_t/(self.rho_d*self.delta_x_d*self.c_d)*k/self.delta_x_d
                        A[m, m] = 1 + self.delta_t/(self.rho_d*self.delta_x_d*self.c_d)*k/self.delta_x_d
            #create array b      
            b = self.T[i-1].copy()
        
            Area = float(eval(self.are.text()))
            b[0] = b[0] + self.delta_t/(self.rho_gas*self.delta_x_a*self.c_gas)*qgen/Area
            b[1] = b[1] + self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*qgen/Area
      
            x = np.linalg.solve(A,b)
        
            
            for j in range(0,self.node):
                self.T[i,j] = x[j].copy()
                
            test = k_T - self.T[i] <= float(eval(self.tolerance.text()))


            for w in range(0,len(test)): #check that temps are in tolerance
                if test[w] == False:
                    k_T[w] = (self.T[i, w] - k_T[w])/2 + k_T[w]
                    hr_T[w] = (self.T[i, w] - hr_T[w])/2 + hr_T[w]
                    break
                elif test[-1] == True:
                    f = open("temporary_transient.txt", 'a+')
                    f.write('  ')
                    if len(str(t)) == 1:
                        f.write(str(t) + '       ')
                    if len(str(t)) == 2:
                        f.write(str(t) + '      ')
                    if len(str(t)) == 3:
                        f.write(str(t) + '     ')
                    if len(str(t)) == 4:
                        f.write(str(t) + '    ')
                    if len(str(t)) == 5:
                        f.write(str(t) + '    ')
                    for j in range(0, len(test)):
                        text = str(round(self.T[i,j], 2))
                        if len(text) == 5:                            
                            f.write(text +  '      ')
                        elif len(text) == 6:                            
                            f.write(text +  '     ')
                        elif len(text) == 7:
                            f.write(text +  '    ')
                        elif len(text) == 8:
                            f.write(text +  '   ')
                        elif len(text) == 9:
                            f.write(text +  '  ')
                        else:
                            f.write(text +  ' ')
                    f.write("\n")
                    f.close()
                    it = 1

                    return 
    
    def transient(self, t):
        """creates matrix ‘A’ and ‘b’ and solves the transient problem for initial temperatures"""
        #Ax = b
        f = open("temporary_transient.txt", 'w+')
        f.write("Transient Temperatures \n  ")
        f.write('time    xgas       ')
        for i in range(0, self.node-1):
            text = 'x' + str(i) 
            if len(text) == 2:
                text = text + '         '
            elif len(text) == 3:
                text = text + '        '
            elif len(text) == 4:
                text = text + '       '
            elif len(text) == 5:
                text = text + '      '
            f.write(text)
        f.write("\n")
        f.close()

        A = np.zeros((self.node, self.node), dtype='float64')
        k_T = np.zeros((self.node))
        hr_T = np.zeros((self.node))
        self.cycle_start = []
        self.cycle_end = []
        j = float(eval(self.qto.text()))
        k = float(eval(self.qt1.text()))
        while j <= self.t_max:
            self.cycle_start.append(j)
            self.cycle_end.append(k)
            j = j + float(eval(self.qto.text()))
            k = k + float(eval(self.qt1.text()))
                
        for i in range(1, int(round(self.t_max/self.delta_t, 0))): #go through all time steps
            t = t + self.delta_t
            
            for m in range(0,len(self.cycle_start)):
                if self.cycle_start[m] <= t <= self.cycle_end[m]:
                    qgen = float(eval(self.qgen_eq.text()))
                else:
                    qgen = 0
                    
            for j in range(0,self.node): #initial
                k_T[j] = self.T[i-1,j]
                hr_T[j] = self.T[i-1,j]
                
            #node inf
            k = self.k_value(k_T[0], self.Ta_array, self.ka_array)
            h_r = self.eps*self.sig*(hr_T[1] + hr_T[0])*(hr_T[1]**2 + hr_T[0]**2)
            A[0, 0] = 1 + self.delta_t/(self.rho_gas*self.delta_x_a*self.c_gas)*h_r
            A[0, 1] = -self.delta_t/(self.rho_gas*self.delta_x_a*self.c_gas)*h_r
    
    
            #node 0
            k = self.k_value(k_T[1], self.Ta_array, self.ka_array)
            A[1, 0] = -self.delta_t/(self.rho_gas*self.delta_x_a*self.c_gas)*h_r
            A[1, 1] = 1 + self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*k/self.delta_x_a + self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*h_r
            A[1, 2] = -self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*k/self.delta_x_a
        
            for m in range(2,self.node):
                if m < 2 + self.nodes_per: #material a
                    k = self.k_value(k_T[m], self.Ta_array, self.ka_array)
                    A[m, m-1] = A[m, m+1] = -self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*k/self.delta_x_a
                    A[m, m] = 1 + 2*self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*k/self.delta_x_a
                elif m == 2 + self.nodes_per: #material a & b
                    k_a = self.k_value(k_T[m], self.Ta_array, self.ka_array)
                    k_b = self.k_value(k_T[m], self.Tb_array, self.kb_array)
                    A[m, m-1] = -self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*k_a/self.delta_x_a
                    A[m, m+1] = -self.delta_t/(self.rho_b*self.delta_x_b*self.c_b)*k_b/self.delta_x_b
                    A[m, m] = 1 + self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*k_a/self.delta_x_a + self.delta_t/(self.rho_b*self.delta_x_b*self.c_b)*k_b/self.delta_x_b
                    
                elif 2 + self.nodes_per <= m <= 2 + 2*self.nodes_per: #material b
                    k = self.k_value(k_T[m], self.Tb_array, self.kb_array)
                    A[m, m-1] = A[m, m+1] = -self.delta_t/(self.rho_b*self.delta_x_b*self.c_b)*k/self.delta_x_b
                    A[m, m] = 1 + 2*self.delta_t/(self.rho_b*self.delta_x_b*self.c_b)*k/self.delta_x_b
                elif m == 2 + 2*self.nodes_per: #material b & c
                    k_b = self.k_value(k_T[m], self.Tb_array, self.kb_array)
                    k_c = self.k_value(k_T[m], self.Tc_array, self.kc_array)
                    A[m, m-1] = -self.delta_t/(self.rho_b*self.delta_x_b*self.c_b)*k_b/self.delta_x_b
                    A[m, m+1] = -self.delta_t/(self.rho_c*self.delta_x_c*self.c_c)*k_c/self.delta_x_c
                    A[m, m] = 1 + self.delta_t/(self.rho_b*self.delta_x_b*self.c_b)*k_b/self.delta_x_b + self.delta_t/(self.rho_c*self.delta_x_c*self.c_c)*k_c/self.delta_x_c
                    
                elif 2 + 2*self.nodes_per <= m <= 2 + 3*self.nodes_per: #material c
                    k = self.k_value(k_T[m], self.Tc_array, self.kc_array)
                    A[m, m-1] = A[m, m+1] = -self.delta_t/(self.rho_c*self.delta_x_c*self.c_c)*k/self.delta_x_c
                    A[m, m] = 1 + 2*self.delta_t/(self.rho_c*self.delta_x_c*self.c_c)*k/self.delta_x_c
                elif m == 2 + 3*self.nodes_per:
                    k_c = self.k_value(k_T[m], self.Tc_array, self.kc_array)
                    k_d = self.k_value(k_T[m], self.Td_array, self.kd_array)
                    A[m, m-1] = -self.delta_t/(self.rho_c*self.delta_x_c*self.c_c)*k_c/self.delta_x_c
                    A[m, m+1] = -self.delta_t/(self.rho_d*self.delta_x_d*self.c_d)*k_d/self.delta_x_d
                    A[m, m] = 1 + self.delta_t/(self.rho_c*self.delta_x_c*self.c_c)*k_c/self.delta_x_c + self.delta_t/(self.rho_d*self.delta_x_d*self.c_d)*k_d/self.delta_x_d
                    
                elif 2 + 3*self.nodes_per <= m <= self.node-2:
                    k = self.k_value(k_T[m], self.Td_array, self.kd_array)
                    A[m, m-1] = A[m, m+1] = -self.delta_t/(self.rho_d*self.delta_x_d*self.c_d)*k/self.delta_x_d
                    A[m, m] = 1 + 2*self.delta_t/(self.rho_d*self.delta_x_d*self.c_d)*k/self.delta_x_d
                elif m == self.node-1:
                    if self.endCon.currentText() == "Constant Temp":
                        A[m, m] = 1
                    else:
                        k = self.k_value(k_T[m], self.Td_array, self.kd_array)
                        A[m, m-1] = -self.delta_t/(self.rho_d*self.delta_x_d*self.c_d)*k/self.delta_x_d
                        A[m, m] = 1 + self.delta_t/(self.rho_d*self.delta_x_d*self.c_d)*k/self.delta_x_d
                    
            b = self.T[i-1].copy()
            
            
            Area = float(eval(self.are.text()))
            b[0] = b[0] + self.delta_t/(self.rho_gas*self.delta_x_a*self.c_gas)*qgen/Area
            b[1] = b[1] + self.delta_t/(self.rho_a*self.delta_x_a*self.c_a)*qgen/Area
            
            x = np.linalg.solve(A,b)

            
            
            for j in range(0,self.node):
                self.T[i,j] = x[j].copy()
                k_T[j] = x[j].copy()
                hr_T[j] = x[j].copy()
            
            self.update(self.T[i], t)
            
            
            
        return A
    def ssIterate(self, k, T_temp):
        """takes in the array of temperatures to be used to determine the thermal conductivities
        and the radiation heat transfer coefficient, the old temperatures, and the heat generation. 
        It then solves for the constants and the temperatures until the old temperatures and the new 
        temperatures are within the tolerance"""
        it = 0
        g = 0 
        while it == 0:
            g = g + 1
            #print('ss', g)
            self.hr = self.eps*self.sig*(k[0] + self.T_gas)*(k[0]**2 + self.T_gas**2)
            

            
            k_a_0 = self.k_value(k[0], self.Ta_array, self.ka_array)
            k_a_b = self.k_value(k[1], self.Ta_array, self.ka_array)
            k_b_a = self.k_value(k[1], self.Tb_array, self.kb_array)
            k_b_c = self.k_value(k[2], self.Tb_array, self.kb_array)
            k_c_b = self.k_value(k[2], self.Tc_array, self.kc_array)
            k_c_d = self.k_value(k[3], self.Tc_array, self.kc_array)
            k_d_c = self.k_value(k[3], self.Td_array, self.kd_array)
#            k_d = self.k_value(k[4], self.Td_array, self.kd_array)
            self.A_ss[3, 6] = -k_a_0
            self.A_ss[3, 7] = self.hr
            self.A_ss[4, 4] = k_b_a * -1
            self.A_ss[4, 6] = k_a_b 
            self.A_ss[5, 2] = k_c_b * -1
            self.A_ss[5, 4] = k_b_c
            self.A_ss[6, 0] = k_d_c
            self.A_ss[6, 2] = k_c_d * -1
    
            self.b_ss[0] = self.L_a**2
            self.b_ss[3] =self.hr * self.T_gas
            
            if self.endCon.currentText() == "Constant Temp":
                self.A_ss[7, 0] = self.L_a + self.L_b + self.L_c + self.L_d
                self.A_ss[7, 1] = 1
                self.b_ss[7] =self.T_const
            
            else:
                self.A_ss[7, 0] = 1
                self.b_ss[7] = 0 
    
                
            c = np.linalg.solve(self.A_ss, self.b_ss)
                        
            x = 0.0
            for i in range(0, len(self.T_ss)):
                if i <= self.nodes_per:
                    #k_a = self.k_value(T_temp[i], self.Ta_array, self.ka_array)
                    T_temp[i] = c[6]*x + c[7]
                    x = x + self.delta_x_a
                elif self.nodes_per <= i <= 2*self.nodes_per:
                    T_temp[i] = c[4]*x + c[5]
                    x = x + self.delta_x_b
                elif 2*self.nodes_per <= i <= 3*self.nodes_per:
                    T_temp[i] = c[2]*x + c[3]
                    x = x + self.delta_x_c
                elif 3*self.nodes_per <= i:
                    T_temp[i] = c[0]*x + c[1]  
                    x = x + self.delta_x_d   
            
            test = []
            for i in range(0, len(T_temp)):
                test.append(self.T_ss[i] - T_temp[i] <= float(eval(self.tolerance.text())))
            for w in range(0,len(test)):
                if test[w] == False:
                    self.T_ss = T_temp
                    k[0] = self.T_ss[0]
                    k[1] = self.T_ss[self.nodes_per]
                    k[2] = self.T_ss[2*self.nodes_per]
                    k[3] = self.T_ss[3*self.nodes_per]
                    k[4] = self.T_ss[-1]
                    break
                elif test[-1] == True:
                    self.T_ss = T_temp
                    f = open("temporary_steady_state.txt", 'a+')
                    f.write('  ')
                    for j in range(0, len(test)):
                        text = str(round(self.T_ss[j], 2))
                        if len(text) == 5:                            
                            f.write(text +  '      ')
                        elif len(text) == 6:                            
                            f.write(text +  '     ')
                        elif len(text) == 7:
                            f.write(text +  '    ')
                        elif len(text) == 8:
                            f.write(text +  '   ')
                        elif len(text) == 9:
                            f.write(text +  '  ')
                        else:
                            f.write(text +  ' ')
                    f.write("\n")
                    f.close()
                    it = 1
                    return

        
    def SteadyState(self):
        """creates matrix ‘A’ and ‘b’ for the steady state problem and solves for 
        the constants and the temperatures initially"""
        f = open("temporary_steady_state.txt", 'w+')
        f.write("Steady State Temperatures \n")

        for i in range(0, self.node-1):
            text = '  x' + str(i) 
            if len(text) == 2:
                text = text + '         '
            elif len(text) == 3:
                text = text + '        '
            elif len(text) == 4:
                text = text + '       '
            elif len(text) == 5:
                text = text + '      '
            f.write(text)
        f.write("\n")
        f.close()
        
        
        self.A_ss = np.zeros((8,8), dtype='float64')
        self.b_ss = np.zeros((8), dtype='float64')

        qgen = 0 #b/c steady state

        
        self.T_gas = float(eval(self.Tgas.text()))
        if self.endCon.currentText() == "Constant Temp":
            self.T_const = float(eval(self.Tend.text()))
        
        k = [self.T_gas/4, self.T_gas/4, self.T_gas/4, self.T_gas/4, self.T_gas/4]
        self.hr = self.eps*self.sig*(k[0] + self.T_gas)*(k[0]**2 + self.T_gas**2)
        
        
        self.A_ss[0, 4] = self.L_a
        self.A_ss[0, 5] = 1
        self.A_ss[0, 6] = self.L_a * -1
        self.A_ss[0, 7] = -1
        
        self.A_ss[1, 2] = self.L_a + self.L_b
        self.A_ss[1, 3] = 1
        self.A_ss[1, 4] = self.L_a + self.L_b * -1
        self.A_ss[1, 5] = -1
        
        self.A_ss[2, 0] = self.L_a + self.L_b + self.L_c
        self.A_ss[2, 1] = 1
        self.A_ss[2, 2] = self.L_a + self.L_b + self.L_c * -1
        self.A_ss[2, 3] = -1
        
        self.A_ss[3, 6] = -1
        self.A_ss[3, 7] = self.hr
        
        

       
        if self.endCon.currentText() == "Constant Temp":
            self.A_ss[7, 0] = self.L_a + self.L_b + self.L_c + self.L_d
            self.A_ss[7, 1] = 1
            self.b_ss[7] =self.T_const
            
        else:
            self.A_ss[7, 0] = 1
            self.b_ss[7] = 0      
            
            
#        k_a_0 = self.k_value(k[0], self.Ta_array, self.ka_array)
        k_a_b = self.k_value(k[1], self.Ta_array, self.ka_array)
        k_b_a = self.k_value(k[1], self.Tb_array, self.kb_array)
        k_b_c = self.k_value(k[2], self.Tb_array, self.kb_array)
        k_c_b = self.k_value(k[2], self.Tc_array, self.kc_array)
        k_c_d = self.k_value(k[3], self.Tc_array, self.kc_array)
        k_d_c = self.k_value(k[3], self.Td_array, self.kd_array)
#        k_d = self.k_value(k[4], self.Td_array, self.kd_array)
        self.A_ss[4, 4] = k_b_a * -1
        self.A_ss[4, 6] = k_a_b 
        self.A_ss[5, 2] = k_c_b * -1
        self.A_ss[5, 4] = k_b_c
        self.A_ss[6, 0] = k_d_c
        self.A_ss[6, 2] = k_c_d * -1

        self.b_ss[0] = self.L_a**2
        self.b_ss[3] =self.hr * self.T_gas

            
        c = np.linalg.solve(self.A_ss, self.b_ss)
        self.T_ss = np.zeros((self.node),dtype='float64')
                
        T_temp = []
        for i in range(0, self.node):
            
            if i <= self.nodes_per:
                x = self.delta_x_a*i
                self.T_ss[i] = c[6]*x + c[7]
                T_temp.append(self.T_ss[i])
            elif self.nodes_per <= i <= 2*self.nodes_per:
                x = self.L_a + self.delta_x_b*(i - self.nodes_per)
                self.T_ss[i] = c[4]*x + c[5]
                T_temp.append(self.T_ss[i])
            elif 2*self.nodes_per <= i <= 3*self.nodes_per:
                x = self.L_a + self.L_b + self.delta_x_c*(i - 2*self.nodes_per)
                self.T_ss[i] = c[2]*x + c[3]
                T_temp.append(self.T_ss[i])
            elif 3*self.nodes_per <= i:
                x = self.L_a + self.L_b + self.L_c + self.delta_x_d*(i - 3*self.nodes_per)
                self.T_ss[i] = c[0]*x + c[1]
                T_temp.append(self.T_ss[i])
        k[0] = self.T_ss[0]
        k[1] = self.T_ss[self.nodes_per]
        k[2] = self.T_ss[2*self.nodes_per]
        k[3] = self.T_ss[3*self.nodes_per]
        k[4] = self.T_ss[-1]
        self.ssIterate(k, T_temp)
        self.update(self.T_ss, 0)
        
    def update(self, y, t):
        """Update the figure.
        """
        self.plot.redraw(self.x_array, y, t)

        
    def run(self):
        """run the program. Pulls the variables and call steady state and transient"""
        #        print('run')
        self.sig = 5.67e-8
        t = 0
        #default values?
        self.nodes_per = int(eval(self.np.text()))
        self.L_a = float(eval(self.La.text()))
        self.L_b = float(eval(self.Lb.text()))
        self.L_c = float(eval(self.Lc.text()))
        self.L_d = float(eval(self.Ld.text()))
        self.nodes_per_d = int(eval(self.npd.text()))
        self.delta_x_a = self.L_a/(self.nodes_per-1) 
        self.delta_x_b = self.L_b/(self.nodes_per-1)
        self.delta_x_c = self.L_c/(self.nodes_per-1)
        self.delta_x_d = self.L_d/(self.nodes_per_d-1)
        self.t_max = float(eval(self.cycle.text()))*float(eval(self.qt1.text()))
        self.delta_t = float(eval(self.dt.text()))
        self.eps = float(eval(self.ep.text()))
        self.node = 2 #node in gas
        self.x_array = [0, 0]
        x = 0
        while x <= self.L_a:
            self.node = self.node + 1
            x = x + self.delta_x_a
            self.x_array.append(self.x_array[-1] + self.delta_x_a)

        x = 0
        while x <= self.L_b:
            self.node = self.node + 1
            x = x + self.delta_x_b
            self.x_array.append(self.x_array[-1] + self.delta_x_b)
            
        x = 0
        while x <= self.L_c:
            self.node = self.node + 1
            x = x + self.delta_x_c
            self.x_array.append(self.x_array[-1] + self.delta_x_c)
        
        x = 0
        while x <= self.L_d:
            self.node = self.node + 1
            x = x + self.delta_x_d
            self.x_array.append(self.x_array[-1] + self.delta_x_d)

        

        mat_a = str(self.box1.currentText()) + '.txt'
        mat_b = str(self.box2.currentText()) + '.txt'
        mat_c = str(self.box3.currentText()) + '.txt'
        mat_d = str(self.box4.currentText()) + '.txt'
        self.Ta_array, self.ka_array = self.k_array(mat_a)
        self.Tb_array, self.kb_array = self.k_array(mat_b)
        self.Tc_array, self.kc_array = self.k_array(mat_c)
        self.Td_array, self.kd_array = self.k_array(mat_d)
        self.melt = [0, 0, 0, 0]
        self.rho_a, self.c_a, self.melt[0] = self.prop(self.box1.currentText())
        self.rho_b, self.c_b, self.melt[1] = self.prop(self.box2.currentText())
        self.rho_c, self.c_c, self.melt[2] = self.prop(self.box3.currentText())
        self.rho_d, self.c_d, self.melt[3] = self.prop(self.box4.currentText())
        self.rho_gas, self.c_gas = self.prop(self.box5.currentText())
        
        
        #Initial Conditions
        self.T = np.zeros((int(round(self.t_max/self.delta_t,0)), self.node), dtype='float64')
        for i in range(0, self.node):
            self.T[0,i] = float(eval(self.Tint.text()))
            
            
        
        if self.runCon.currentText() == "Transient":
            #Initial Conditions
            self.T = np.zeros((int(round(self.t_max/self.delta_t,0)), self.node), dtype='float64')
            for i in range(0, self.node):
                self.T[0,i] = float(eval(self.Tint.text()))    
            self.update(self.T[0], 0 )
            self.transient(t)
            self.hot = self.Hot_transient('temporary_transient.txt')
            self.TPlot = SecondWindow(self.x_array, self.T, self.hot, self.melt)
            self.TPlot.show()
        elif self.runCon.currentText() == "Steady State":
            self.SteadyState()
            self.hotS = self.Hot_steady('temporary_steady_state.txt')
            self.SPlot = ThirdWindow(self.x_array, self.T_ss, self.hotS, self.melt)
            self.SPlot.show()
        elif self.runCon.currentText() == 'Both':
            self.T = np.zeros((int(round(self.t_max/self.delta_t,0)), self.node), dtype='float64')
            for i in range(0, self.node):
                self.T[0,i] = float(eval(self.Tint.text()))    
            self.update(self.T[0], 0)
            self.transient(t)
            self.hot = self.Hot_transient('temporary_transient.txt')
            self.TPlot = SecondWindow(self.x_array, self.T, self.hot, self.melt)
            self.TPlot.show()
            self.SteadyState()
            self.hotS = self.Hot_steady('temporary_steady_state.txt')
            self.SPlot = ThirdWindow(self.x_array, self.T_ss, self.hotS, self.melt)
            self.SPlot.show()
            


        
        
class MatplotlibCanvas(FigureCanvas) :
    """ This is borrowed heavily from the matplotlib documentation;
        specifically, see:
        http://matplotlib.org/examples/user_interfaces/embedding_in_qt5.html
    """
    def __init__(self):
        
        # Initialize the figure and axes
        self.fig = Figure()
        self.axes = self.fig.add_subplot(111)
        
        # Give it some default plot (if desired).  
        self.axes.set_xlabel('x (m)')
        self.axes.set_ylabel('Temperature (K)')   
        self.axes.set_title('Plot of Temperature through the solid')        
        # Now do the initialization of the super class
        FigureCanvas.__init__(self, self.fig)
        #self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        
    def ddraw(self, x, y, a, b, T_SS):
        """takes in the array of x values, temperature values for the results windows, 
        and if it is the steady state or transient problem.
        It then draws the plot in the appropriate result window"""
        self.axes.clear()
        self.axes.set_xlabel('x (m)')
        self.axes.set_ylabel('Temperature (K)') 
        self.axes.set_title('Plot of Temperature through the solid')
        if T_SS == 1:
            self.axes.plot(x, y)
            self.axes.legend(["Steady State"])
        elif T_SS == 2:
            self.axes.plot(x, y, a, b)
            self.axes.legend(["At Initial Time", "At End Time"])
        self.draw()
        app.processEvents()
        
    def redraw(self, x, y, t) :
        """ Redraw the figure with new x and y values.
        """
        # clear the old image (axes.hold is deprecated)
        self.axes.clear()
        self.axes.set_xlabel('x (m)')
        self.axes.set_ylabel('Temperature (K)') 
        self.axes.set_title('Plot of Temperature through the solid, t = %0.2f sec' % t)
        self.axes.plot(x[1:], y[1:])
        self.draw()
        app.processEvents()
app = QApplication(sys.argv)
form = MainWindow()
form.show()
app.exec_()
        