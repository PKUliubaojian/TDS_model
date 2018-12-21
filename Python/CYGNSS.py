# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 15:31:48 2018

@author: RS101
"""

import numpy as np
from netCDF4 import Dataset
import os
import geopandas as gpd
from shapely.geometry import Point
from pandas import DataFrame

class DDM_L1:
    '''
        CYGNSS Level 1 DDM bistatic radar cross section data reader.
    '''
    def __init__(self,DDMfilename):
        DDM=Dataset(DDMfilename,'r')
        ts=np.array(DDM['ddm_timestamp_utc']).repeat(4).reshape(-1,4).reshape(-1,order='F')
        Lat=np.array(DDM['sp_lat']).reshape(-1)
        Lon=np.array(DDM['sp_lon']).reshape(-1)
        Alt=np.array(DDM['sp_alt']).reshape(-1)
        LES=np.array(DDM['ddm_les']).reshape(-1)
        nbrcs=np.array(DDM['ddm_nbrcs']).reshape(-1)
        P_t=np.array(DDM['gps_tx_power_db_w']).reshape(-1)
        G_t=np.array(DDM['gps_ant_gain_db_i']).reshape(-1)
        G_r=np.array(DDM['sp_rx_gain']).reshape(-1)
        R_t=np.array(DDM['tx_to_sp_range']).reshape(-1)
        R_r=np.array(DDM['rx_to_sp_range']).reshape(-1)
        I_angle=np.array(DDM['sp_inc_angle']).reshape(-1)
        SNR=np.array(DDM['ddm_snr'][:]).reshape(-1)
        
        df=DataFrame({'ts':ts,'Lat':Lat,'Lon':Lon,'Alt':Alt,'SNR':SNR,'nbrcs':nbrcs,
                      'LES':LES,'P_t':P_t,'G_t':G_t,'G_r':G_r,'R_t':R_t,'R_r':R_r,
                      'I_angle':I_angle})
        mask=(df["LES"]!=-9999)&(df["nbrcs"]!=-9999)&(df["SNR"]!=-9999)
        ddms=np.array(DDM['brcs'][:]).reshape(-1,17,11)
        self.ddms=ddms[mask,:,:]
        self.df=df[mask]