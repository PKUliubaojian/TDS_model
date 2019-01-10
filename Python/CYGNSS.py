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
import xarray as xr

class DDM_L1:
    '''
        CYGNSS Level 1 DDM bistatic radar cross section data reader.
    '''
    def __init__(self,DDMfilename):
        lookup={'Lat':'sp_lat','Lon':'sp_lon','Alt':'sp_alt','LES':'ddm_les',
            'nbrcs':'ddm_nbrcs','P_t':'gps_tx_power_db_w','G_t':'gps_ant_gain_db_i',
            'G_r':'sp_rx_gain','R_t':'tx_to_sp_range','R_r':'rx_to_sp_range','Track_ID':'track_id',
            'I_angle':'sp_inc_angle','SNR':'ddm_snr','nbrcs_area':'nbrcs_scatter_area',
            'peak_row':'brcs_ddm_peak_bin_delay_row','peak_col':'brcs_ddm_peak_bin_dopp_col'}
        d=Dataset(DDMfilename,'r')
        dicts=[[key,np.array(d[value]).reshape(-1)] for key,value in lookup.items()]
        df=DataFrame(dict(dicts))
        df['ts']=np.array(d['ddm_timestamp_utc']).repeat(4).reshape(-1,4).reshape(-1,order='F')
        mask=(df["LES"]!=-9999)&(df["nbrcs"]!=-9999)&(df["SNR"]!=-9999)
        ddms=np.array(d['brcs'][:]).reshape(-1,17,11)
        areas=np.array(d['eff_scatter'][:]).reshape(-1,17,11)
        self.ddms=ddms[mask,:,:]
        self.df=df[mask]
        self.areas=areas[mask]
        
class DDM_L2:
    '''
        CYGNSS Level 2 MSS & ocean wind speed data.
    '''
    def __init__(self,fname):
        lookup={'Lat':'lat','Lon':'lon','mss':'mean_square_slope','mss_uncertainty':'mean_square_slope_uncertainty',
            'time':'sample_time','windspeed':'wind_speed','ws_uncertainty':'wind_speed_uncertainty',
            'fresnel_coeff':'fresnel_coeff','les':'les_mean','ddma':'nbrcs_mean','RCG':'range_corr_gain','azimuth_angle':'azimuth_angle',
            'incidence_angle':'incidence_angle'
            }
        d=Dataset(fname,'r')
        dicts=[[key,np.array(d[value]).reshape(-1)] for key,value in lookup.items()]
        df=DataFrame(dict(dicts))
        self.df=df

        
    
def DDMA(ddms,areas,p_delay,p_doppler,delaybin=3,dopplerbin=5):
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

def LES(ddms,areas,p_delay,p_doppler,delaybin=3,dopplerbin=5):
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

def TES(ddms,areas,p_delay,p_doppler,delaybin=5,dopplerbin=5):
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

