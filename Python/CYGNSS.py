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
        Track_ID=np.array(DDM['track_id']).reshape(-1)
        I_angle=np.array(DDM['sp_inc_angle']).reshape(-1)
        SNR=np.array(DDM['ddm_snr'][:]).reshape(-1)
        nbrcs_area=np.array(DDM['nbrcs_scatter_area'][:]).reshape(-1)
        peak_row=np.array(DDM['brcs_ddm_peak_bin_delay_row'][:]).reshape(-1)
        peak_col=np.array(DDM['brcs_ddm_peak_bin_dopp_col'][:]).reshape(-1)
        df=DataFrame({'ts':ts,'Lat':Lat,'Lon':Lon,'Alt':Alt,'SNR':SNR,'nbrcs':nbrcs,
              'LES':LES,'P_t':P_t,'G_t':G_t,'G_r':G_r,'R_t':R_t,'R_r':R_r,
              'I_angle':I_angle,'Track_ID':Track_ID,'nbrcs_area':nbrcs_area,
              'peak_row':peak_row,'peak_col':peak_col})
        mask=(df["LES"]!=-9999)&(df["nbrcs"]!=-9999)&(df["SNR"]!=-9999)
        ddms=np.array(DDM['brcs'][:]).reshape(-1,17,11)
        areas=np.array(DDM['eff_scatter'][:]).reshape(-1,17,11)
        self.ddms=ddms[mask,:,:]
        self.df=df[mask]
        self.areas=areas[mask]
        
def DDMA_CY(ddms,areas,p_delay,p_doppler,delaybin=3,dopplerbin=5):
    """
        输入的ddms和area都是一个三维矩阵，第一个维度为ddm数量，第二个维度delay固定为17，第三个维度doppler固定为10
        ddms: CYGNSS L1b DDM
        area: CYGNSS L1b effective scattering area
        p_delay: CYGNSS L1b peak delay bin in pixel
        p_doppler: CYGNSS L1b peak doppler bin in pixel
        输出:DDMA
        建议把ddm串起来计算以提高计算效率
    """
    N=np.arange(ddms.shape[0])
    x_d=delaybin//2
    y_d=dopplerbin//2
    x_l,y_l=np.meshgrid(np.arange(-x_d,x_d+1,1),np.arange(-y_d,y_d+1,1))
    x_l=x_l.reshape(-1)
    y_l=y_l.reshape(-1)
    sum_brcs=[ddms[N,(p_delay+dx)%16,(p_doppler+dy)%10] for (dx,dy) in zip(x_l,y_l)]
    sum_area=[areas[N,(p_delay+dx)%16,(p_doppler+dy)%10] for (dx,dy) in zip(x_l,y_l)]
    c_brcs=np.stack(sum_brcs)
    c_area=np.stack(sum_area)
    DDMA=np.sum(c_brcs,axis=0)/np.sum(c_area,axis=0)
    return DDMA

def LES_CY(ddms,areas,p_delay,p_doppler,delaybin=3,dopplerbin=5):
    """
        输入的ddms和area都是一个三维矩阵，第一个维度为ddm数量，第二个维度delay固定为17，第三个维度doppler固定为10
        ddms: CYGNSS L1b DDM
        area: CYGNSS L1b effective scattering area
        p_delay: CYGNSS L1b peak delay bin in pixel
        p_doppler: CYGNSS L1b peak doppler bin in pixel
        输出：LES
        建议把ddm串起来计算以提高计算效率
    """
    N=np.arange(ddms.shape[0])
    x_d=delaybin//2
    y_d=dopplerbin//2
    x_series=np.arange(-x_d,x_d+1,1)
    y_series=np.arange(-y_d,y_d+1,1)
    x_l,y_l=np.meshgrid(x_series,y_series)
    x_l=x_l.reshape(-1)
    y_l=y_l.reshape(-1)
    sum_area=[areas[N,(p_delay+dx)%16,(p_doppler+dy)%10] for (dx,dy) in zip(x_l,y_l)]
    c_area=np.sum(np.stack(sum_area),axis=0)    
    peak_IDW=[ddms[N,p_delay,(p_doppler+dy)%10] for dy in y_series]
    leading_IDW=[ddms[N,(p_delay-delaybin),(p_doppler+dy)%10] for dy in y_series]
    peak_IDW=np.sum(np.stack(peak_IDW),axis=0)
    leading_IDW=np.sum(np.stack(leading_IDW),axis=0)
    leading=(peak_IDW-leading_IDW)/delaybin
    LES=leading/c_area
    return LES

def TES_CY(ddms,areas,p_delay,p_doppler,delaybin=5,dopplerbin=5):
    """
        输入的ddms和area都是一个三维矩阵，第一个维度为ddm数量，第二个维度delay固定为17，第三个维度doppler固定为10
        ddms: CYGNSS L1b DDM
        area: CYGNSS L1b effective scattering area
        p_delay: CYGNSS L1b peak delay bin in pixel
        p_doppler: CYGNSS L1b peak doppler bin in pixel
        输出：TES
        建议把ddm串起来计算以提高计算效率
    """
    N=np.arange(ddms.shape[0])
    x_d=delaybin//2
    y_d=dopplerbin//2
    x_series=np.arange(-x_d,x_d+1,1)
    y_series=np.arange(-y_d,y_d+1,1)
    x_l,y_l=np.meshgrid(x_series,y_series)
    x_l=x_l.reshape(-1)
    y_l=y_l.reshape(-1)
    sum_area=[areas[N,(p_delay+dx)%16,(p_doppler+dy)%10] for (dx,dy) in zip(x_l,y_l)]
    peak_IDW=[ddms[N,p_delay,(p_doppler+dy)%10] for dy in y_series]
    trailing_IDW=[ddms[N,(p_delay+delaybin)%17,(p_doppler+dy)%10] for dy in y_series]
    c_area=np.sum(np.stack(sum_area),axis=0)
    peak_IDW=np.sum(np.stack(peak_IDW),axis=0)
    trailing_IDW=np.sum(np.stack(trailing_IDW),axis=0)
    trailing=(peak_IDW-trailing_IDW)/delaybin
    TES=trailing/c_area
    return TES

