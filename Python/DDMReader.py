# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 13:15:27 2018

@author: RS101
"""

import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os
#%%Read specified DDMs
def readDDM(DDMfolder,SNRrange=[0,50],SP_class='ocean'):
    "e.g. D:\\data\\GNSS\\TDS\\L1B\\2015-09\\09\\H00"
    files=os.listdir(DDMfolder)
    fileDDMs=os.path.join(DDMfolder,"DDMs.nc")
    filemetadata=os.path.join(DDMfolder,"metadata.nc")
    if os.path.exists(filemetadata):
        metadata=Dataset(filemetadata)
        DDMs=Dataset(filemetadata)
        N_groups=len(metadata)
        group_num=np.arange(0,N_groups,dtype=np.int32)
        group_str= np.array(['{:0>6}'.format(x) for x in group_num])
        for group_name in group_str:
            gp_metadata=metadata.groups[group_name].NumberOfDelayPixels
            gp_DDM=DDMs.groups[group_name]
            SNR=gp_metadata=metadata.groups[group_name].SNR
    else:
        return []
    


#%%
def xyz2latlon(x,y,z):
    #x,y,z to latitude , longitude and altitude
    f = 1/298.257223563;        #   WGS-84 Flattening.
    e = np.sqrt(f*(2 - f));        #   Eccentricity.
    omega_ie = 7.292115e5;      #   WGS-84 Earth rate (rad/s).
    R_0 = 6378137;              #   WGS-84 equatorial radius (m).                            
    R_P = R_0*(1 - f);          #   Polar radius (m).
    mu_E = 3.986004418e14;      #   WGS-84 Earth's gravitational
    lon = np.arctan2(y,x)*(180/np.pi);
    p=np.linalg.norm(np.c_[x,y],axis=1)
    E=np.sqrt(R_0**2-R_P**2)
    F = 54*(R_P*z)**2;
    G = p**2 + (1 - e**2)*z**2 - (e*E)**2;
    c = e**4*F*p**2/G**3;
    s = (1 + c + np.sqrt(c**2 + 2*c))**(1/3);
    P = (F/(3*G**2))/((s + (1/s) + 1)**2);
    Q = sqrt(1 + 2*e**4*P);
    k_1 = -P*e**2*p/(1 + Q);
    k_2 = 0.5*R_0**2*(1 + 1/Q);
    k_3 = -P*(1 - e**2)*z**2/(Q*(1 + Q));
    k_4 = -0.5*P*p**2;
    r_0 = k_1 + np.sqrt(k_2 + k_3 + k_4);
    k_5 = (p - e**2*r_0);
    U = np.sqrt(k_5**2 + z**2);
    V = np.sqrt(k_5**2 + (1 - e**2)*z**2);
    alt = U*(1 - (R_P**2/(R_0*V)));
    #   Compute additional values required for computing
    #   latitude
    z_0 = (R_P**2*z)/(R_0*V);
    e_p = (R_0/R_P)*e;
    lat = np.arctan((z + z_0*(e_p)**2)/p)*(180/np.pi);
    return lon,lat,alt


#%%