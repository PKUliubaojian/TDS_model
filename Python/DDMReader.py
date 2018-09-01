# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 13:15:27 2018

@author: RS101
"""

import numpy as np
from netCDF4 import Dataset
import os
import geopandas as gpd
from shapely.geometry import Point
from pandas import DataFrame
#%%Read specified DDMs
class DDM:
    def __init__(self,DDMfolder,SNR_Min=0,SNR_Max=50,Incidence_Min=0,Incidence_Max=90,
                 Gain_Min=-20,Gain_Max=30,SP_class='all'):
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
        SpecularPathRangeOffset=np.zeros(0)
        ddm_count=0
        firstgroup=metadata.groups[group_str[0]]
        
        TrackingOffsetDelayInPixels=np.zeros(N_groups)
        TrackingOffsetDopplerHz=np.zeros(N_groups)
        DopplerResolution=np.zeros(N_groups)
        #Attributes
        self.DelayResolution=firstgroup.DelayResolution
        self.DopplerResolution=firstgroup.DopplerResolution
        self.NumberOfDelayPixels=firstgroup.NumberOfDelayPixels
        self.NumberOfDopplerPixels=firstgroup.NumberOfDopplerPixels
        self.TrackingOffsetDelayInPixels=firstgroup.TrackingOffsetDelayInPixels
        self.TrackingOffsetDopplerHz=firstgroup.TrackingOffsetDopplerHz
        self.TrackingOffsetDopplerInPixels=firstgroup.TrackingOffsetDopplerInPixels
        DDMs_stack=[]
        for (group_num,group_name) in zip(group_num,group_str):
            grp=metadata.groups[group_name]
            TrackingOffsetDelayInPixels[group_num]=grp.TrackingOffsetDelayInPixels
            TrackingOffsetDopplerHz[group_num]=grp.TrackingOffsetDopplerHz
            DopplerResolution[group_num]=grp.DopplerResolution
            g=grp.variables
            
            g_SNR=np.array(g["DDMSNRAtPeakSingleDDM"])            
            g_Gain=np.array(g["AntennaGainTowardsSpecularPoint"])
            g_IncidenceAngle=np.array(g["SPIncidenceAngle"])
            g_SPlat=np.array(g["SpecularPointLat"])
            g_SPlon=np.array(g["SpecularPointLon"])
            g_SpecularPathRangeOffset=np.array(g["SpecularPathRangeOffset"])
            
            mask=g_SNR>=SNR_Min
            mask=np.logical_and(mask,g_SNR<SNR_Max)
            mask=np.logical_and(mask,g_Gain>=Gain_Min)
            mask=np.logical_and(mask,g_Gain<=Gain_Max)
            mask=np.logical_and(mask,g_SNR<SNR_Max)
            mask=np.logical_and(mask,g_IncidenceAngle>=Incidence_Min)
            mask=np.logical_and(mask,g_IncidenceAngle<=Incidence_Max)
            if SP_class=='all':
                pass
            else:
                mask[mask]=SPwithin(g_SPlon[mask],g_SPlat[mask],SP_class)
            g_N=mask.sum()
            if(g_N>0):
                SNR=np.append(SNR,g_SNR[mask])
                Gain=np.append(Gain,g_Gain[mask])
                IncidenceAngle=np.append(IncidenceAngle,g_IncidenceAngle[mask])
                SPlat=np.append(SPlat,g_SPlat[mask])
                SPlon=np.append(SPlon,g_SPlon[mask])
                SpecularPathRangeOffset=np.append(SpecularPathRangeOffset,g_SpecularPathRangeOffset[mask])
                g_ddm=np.array(DDM_nc.groups[group_name].variables["DDM"])
                g_ddm=g_ddm[mask,:,:]
#                if(ddm_count==0):
#                    DDMs=g_ddm
#                else:
#                    DDMs=np.append(DDMs,g_ddm,axis=0)
                DDMs_stack.append(g_ddm)
            else:
                pass
            ddm_count+=g_N
        DDMs=np.concatenate(DDMs_stack,axis=0)
        self.N=ddm_count
        self.DDMs=DDMs
        self.SNR=SNR
        self.IncidenceAngle=IncidenceAngle
        self.Gain=Gain
        self.SpecularPathRangeOffset=SpecularPathRangeOffset
        metadata.close()
        DDM_nc.close()
    def DDMA(self,DopplerWinWidth,DelayWinWitdth,DopplerUnit='pixel',DelayUnit='pixel'):
        '''
        Calculte the Delay Doppler Map Average near the DDM centre
        DopplerWinWidth: Doppler window width
        DelayWinWitdth： width of Delay window
        '''
        NDoP=self.NumberOfDopplerPixels
        if DopplerUnit=='pixel':
            w_Do=DopplerWinWidth//2
            C_Do=NDoP//2

        elif DopplerUnit=='Hz':
            w_Do=DopplerWinWidth/self.DopplerResolution//2
            C_Do=NDoP//2
        else:
            raise RuntimeError('Illegal parameter value in DopplerUnit : {}'.format(DopplerUnit))
        min_Do=C_Do-w_Do
        max_Do=C_Do+w_Do+1
        if DelayUnit=='ns':
            w_De=DelayWinWitdth/self.DelayResolution//2
        elif DelayUnit=='pixel':
            w_De=DelayWinWitdth//2
           # C_De=TrackingOffsetDelayInPixels-
        else:
            raise RuntimeError('Illegal parameter value in Delay Unit : {}'.format(DelayUnit))
        DDMslice=self.DDMs[:,min_Do:max_Do,:]
        De_min=(self.TrackingOffsetDelayInPixels-self.SpecularPathRangeOffset/self.DelayResolution-w_De).astype(np.int32)
        De_max=(self.TrackingOffsetDelayInPixels-self.SpecularPathRangeOffset/self.DelayResolution+w_De+1).astype(np.int32)
        ddma=[np.mean(DDMslice[i,:,De_min[i]:De_max[i]]) for i in np.arange(self.N)]
        ddma=np.array(ddma)
        return ddma
    
#%% x,y,z to latitude , longitude and altitude
def xyz2latlon(x,y,z):
    #x,y,z to latitude , longitude and altitude
    f = 1/298.257223563;        #   WGS-84 Flattening.
    e = np.sqrt(f*(2 - f));        #   Eccentricity.
    R_0 = 6378137;              #   WGS-84 equatorial radius (m).                            
    R_P = R_0*(1 - f);          #   Polar radius (m).
    lon = np.arctan2(y,x)*(180/np.pi);
    p=np.linalg.norm(np.c_[x,y],axis=1)
    E=np.sqrt(R_0**2-R_P**2)
    F = 54*(R_P*z)**2;
    G = p**2 + (1 - e**2)*z**2 - (e*E)**2;
    c = e**4*F*p**2/G**3;
    s = (1 + c + np.sqrt(c**2 + 2*c))**(1/3);
    P = (F/(3*G**2))/((s + (1/s) + 1)**2);
    Q = np.sqrt(1 + 2*e**4*P);
    k_1 = -P*e**2*p/(1 + Q);
    k_2 = 0.5*R_0**2*(1 + 1/Q);
    k_3 = -P*(1 - e**2)*z**2/(Q*(1 + Q));
    k_4 = -0.5*P*p**2;
    r_0 = k_1 + np.sqrt(k_2 + k_3 + k_4);
    k_5 = (p - e**2*r_0);
    U = np.sqrt(k_5**2 + z**2);
    V = np.sqrt(k_5**2 + (1 - e**2)*z**2);
    alt = U*(1 - (R_P**2/(R_0*V)));
    z_0 = (R_P**2*z)/(R_0*V);
    e_p = (R_0/R_P)*e;
    lat = np.arctan((z + z_0*(e_p)**2)/p)*(180/np.pi);
    return lon,lat,alt
#%%If a secular point is in land,lake or ocean
def SPwithin(Lon,Lat,source="land"):
    '''
    Find specular Points in land,lake,or ocean
    Return a list of boolean
    source : 'land','lake',or'ocean',or shapefile path provided by user's
    '''
    if source in ['land','ocean']:
        path="shp\\ne_110m_land\\ne_110m_land.shp"

    elif source=='lake':
        path="shp\\ne_10m_lakes\\ne_10m_lakes.shp"
    else:
        path=source
    df=DataFrame({"Coordinates":list(zip(Lon,Lat))})
    df['Coordinates'] = df['Coordinates'].apply(Point)
    pts = gpd.GeoDataFrame(df, geometry='Coordinates')
    shapes=gpd.read_file(path)
    union=shapes["geometry"].unary_union
    pts_in=pts.geometry.within(union)
    if source=='ocean':
        pts_in=~pts_in
    return pts_in
        