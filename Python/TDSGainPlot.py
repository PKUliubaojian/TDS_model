# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 15:27:55 2018

@author: Baojian Liu
Plot the TDS Gain
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm


filename='''data\\GainTDS.npy'''
Gain=np.load(filename)
Azim=np.linspace(0,359,360)
Atitu=np.linspace(0,90,91)
X,Y=np.meshgrid(Azim,Atitu)
area=np.power(Atitu,2)/100
fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')
s=ax.scatter(X,Y, c=Gain.T, cmap='jet', alpha=0.75)
ax.set_rmax(90)
cb=plt.colorbar(s,ax=ax)
cb.set_label('TDS1 Antenna Gain(DB)')
