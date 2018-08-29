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
from pandas import DataFrame,Series

#%%Read specified DDMs
class DDM:
    def __init__(self,DDMfolder,SNR_Min=0,SNR_Max=50,Incidence_Min=0,Incidence_Max=90,
                 Gain_Min=-20,Gain_Max=30,SP_class='ocean'):
        "D:\\data\\GNSS\\TDS\\L1B\\2015-09\\09\\H00"
        fileDDMs=os.path.join(DDMfolder,"DDMs.nc")
        filemetadata=os.path.join(DDMfolder,"metadata.nc")
        metadata=Dataset(filemetadata,'r')
        DDM_nc=Dataset(fileDDMs,'r')
        N_groups=len(metadata.groups)
        group_num=np.arange(0,N_groups,dtype=np.int32)
        group_str= np.array(['{:0>6}'.format(x) for x in group_num])
        SNR=np.zeros(0)
        Gain=np.zeros(0)
        IncidenceAngle=np.zeros(0)
        SPlat=np.zeros(0)
        SPlon=np.zeros(0)
        ddm_count=0
        for (idx,group_name) in zip(group_num,group_str):
            g=metadata.groups[group_name].variables
            g_SNR=np.array(g["DDMSNRAtPeakSingleDDM"])            
            g_Gain=np.array(g["AntennaGainTowardsSpecularPoint"])
            g_IncidenceAngle=np.array(g["SPIncidenceAngle"])
            g_SPlat=np.array(g["SpecularPointLat"])
            g_SPlon=np.array(g["SpecularPointLon"])
            SNR=np.append(SNR,g_SNR)
            Gain=np.append(Gain,g_Gain)
            IncidenceAngle=np.append(IncidenceAngle,g_IncidenceAngle)
            SPlat=np.append(SPlat,g_SPlat)
            SPlon=np.append(SPlon,g_SPlon)
            
            mask=np.logical_and((g_SNR>=SNR_Min) ,(g_SNR<SNR_Max) ,
                  (g_Gain>=Gain_Min) , (g_Gain<=Gain_Max) , 
                  (g_IncidenceAngle>=Incidence_Min) ,
                  (g_IncidenceAngle<=Incidence_Max))
            g_N=mask.sum()
            ddm_count=ddm_count+g_N
            
            g_ddm=np.array(DDM_nc.groups[group_name].variables["DDM"])
            g_ddm=g_ddm[mask,:,:]
            if(ddm_count==0 and g_N!=0):
                DDMs=g_ddm
            else:
                DDMs=np.append(DDMs,g_ddm,axis=0)
        self.DDMs=DDMs
        self.SNR=SNR
        self.IncidenceAngle=IncidenceAngle
        self.Gain=Gain
        metadata.close()
        DDM_nc.close()
#%% x,y,z to latitude , longitude and altitude
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


