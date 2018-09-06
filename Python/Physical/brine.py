#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Collection of subfunctions related to brine parameters: Brine salinity, 
brine conductivity, brine volumen etc.

@author: mwi
"""
import numpy as np

def salinity(T,method='Assur'):
    """
    Brine salinity as function of temperature.
    Empirical expression from Assur (1960) and Poe et al (1972).
    Version from Ulaby et al, 1986 (E63). Valid range: -43.2oC--1.8oC.
    
    [sal] = salinity(T,method)
        sal:    Salinity [psu (or o/oo)]
        T:      Temperature [K]
        method: 'Assur' (default), 'unknown'
    """
    # Temperature in oC:
    Tc = T-273.15 
    
    # Salinity:
    sal = np.ones(len(Tc))
    
    if method == 'Assur':
        mask = (Tc>=-1.8) # Using Tc=-1.8oC (not 2oC) to ensure continuity of salinity to value of 34.2o/oo.    
        sal[mask]=34.2 # Constant level above -1.8oC
        
        mask = (Tc<-1.8) & (Tc>=-8.2)
        sal[mask]=1.725-18.756*Tc[mask]-0.3964*Tc[mask]**2
    
        mask = (Tc<-8.2) & (Tc>=-22.9)
        sal[mask]=57.041-9.929*Tc[mask]-0.16204*Tc[mask]**2-0.002396*Tc[mask]**3
    
        mask = (Tc<-22.9) & (Tc>=-36.8)
        sal[mask]=242.94+1.5299*Tc[mask]+0.0429*Tc[mask]**2
        
        mask = (Tc<-36.8) & (Tc>=-43.2)
        sal[mask]=508.18+14.535*Tc[mask]+0.2018*Tc[mask]**2
    
        # Constant level for temperatures below -43.2oC:
        mask = (Tc<-43.2)
        sal[mask]=508.18+14.535*-43.2+0.2018*((-43.2)**2) 
        
    elif method == 'unknown':
        sal[T>=0] = 0
        
        mask = (Tc<0) & (Tc>=-8)
        sal[mask]= 1./(0.001-0.05411/Tc[mask]) 
        
        mask = (Tc<-8) & (Tc>=-22.9)
        sal[mask] = -1.20 - 21.8*Tc[mask] - 0.919*Tc[mask]**2 - 0.0178*Tc[mask]**3
    
        mask = (Tc<-22.9) & (Tc>=-36.8)
        sal[mask]=242.94+1.5299*Tc[mask]+0.0429*Tc[mask]**2
        
        mask = (Tc<-36.8) & (Tc>=-43.2)
        sal[mask]=508.18+14.535*Tc[mask]+0.2018*Tc[mask]**2
    
        # Constant level for temperatures below -43.2oC:
        mask = (Tc<-43.2)
        sal[mask]=508.18+14.535*-43.2+0.2018*((-43.2)**2) 
    return sal

def normality(sal):
    """
    Normality of brine solution. 
    
    [N] = normality(sal)
       N:   Normality [unit?]
       sal: Salinity [psu]
    """
    
    # Normality of brine solution:
    N = 0.9141*sal* (1.707e-2 + 1.205e-5*sal + 4.058e-9*sal**2)
    return N    
   
def conductivity(T,sal,method='StogrynDesargant1985'):
    """
    Ionic conductivity of brine as function of brine temperature and salinity. 
    Version from Ulaby et al, 1986 (E20). Valid range: -43.2oC--1.8oC.
    Note that this version is somewhat different from what is given by Stogryn 
    and Desargant, 1985, which is calculated based on temperature only. 
    
    [cond, N] = conductivity(T,N)
        cond:   Brine conductivity [in S/m?]
        N:      Normality
        T:      Temperature [K]
        sal:    Salinity [psu]
        method: 'Stogryn1971', 'StogrynDesargant1985' (default)
    """    
    # Temperature in oC:
    Tc = T-273.15
    
    if method == 'Stogryn1971':
        # Normality of brine solution:
        N = normality(sal)
        
        # Conductivity:
        D = 25-Tc
        sig = N* (10.39-2.378*N+0.683*(N**2)-0.135*N**3+1.01e-2*N**4)
        c = 1.0 - 1.96e-2*D + 8.08e-5*(D**2) - N*D*(3.02e-5+3.92e-5*D + N*(1.72e-5 - 6.58e-6*D))
        cond=c*sig
        
    elif method == 'StogrynDesargant1985':
        # Alternative solution from Stogryn and Desargant, 1985, Eq.7:
        cond=np.zeros(len(Tc))
        mask = (Tc>=-22.9)
        cond[mask] = -Tc[mask]*np.exp(0.5193 + 0.8755*0.1*Tc[mask])
        cond[~mask] = -Tc[~mask]*np.exp(1.0334 + 0.1100*Tc[~mask])
        
        # Calculated normality is not required:
        N = None
        
    # Always positive:
    cond = np.clip(cond,0,None) 
    return cond, N

def volume(T,sal,method='Frankenstein'):
    """
    Volumen fraction of brine in sea ice.
    volbrine = volume(T,sal)
        volbrine: Volumen fraction of brine
        T:        Temperature [K]
        sal:      Sea ice salinity [o/oo, or psu, or ppt] 
        method:   'Frankenstein' (default), 'Frankenstein_simple', 'original'
    """   
    # Temperature in oC:
    Tc = T-273.15

    # From Ulaby et al, 1986, E71, and Ulaby p. 136, E4.51
    # Empirical expression from Frankenstein and Garner (1967): 
    # Applicable for the temperature range: -22.9oC- -0.5oC 
    if method == 'Frankenstein_simple':
        volbrine = np.zeros(len(Tc))
        mask = (Tc<-0.1)
        volbrine[mask] = 0.001*sal[mask]*((-49.185/Tc[mask])+0.532) 

        # Frankenstein equation cannot be calculated for Tc=0oC.
        # Using a different relationship for temperatures above -0.1oC:
        volbrine[~mask] = 0.001*sal[~mask]*9.717 # Not sure where this one comes from
        
    if method == 'Frankenstein':
        # Equations taken from Ulaby's old book, p. 2048:
        volbrine = np.zeros(len(Tc))
        mask = (Tc>-0.5)
        volbrine[mask] = 0.001*sal[mask]*9.717 # Not sure where this one comes from
        mask = (Tc<=-0.5) & (Tc>-2.06)
        volbrine[mask] = 0.001*sal[mask]*(-52.56/Tc[mask]-2.28)
        mask = (Tc<=-2.06) & (Tc>-8.2)
        volbrine[mask] = 0.001*sal[mask]*(-45.917/Tc[mask]+0.93)
        mask = (Tc<=-8.2) & (Tc>=-22.9)
        volbrine[mask] = 0.001*sal[mask]*(-43.795/Tc[mask]+1.189)
        mask = (Tc<-22.9)
        volbrine[mask] = 0.001*sal[mask]*(-43.795/(-22.9)+1.189)
        
    elif method == 'original':
        sal_brine = salinity(T)
        volbrine = sal/sal_brine #*rho_ice/rho_brine  
        # rho_brine is also temperature dependent, but much weaker than sal_brine. 
        # Calculated values may be some factor off the true values??
    
    volbrine = np.clip(volbrine,0,1)
    return volbrine