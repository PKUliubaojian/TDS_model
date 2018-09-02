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
import PIL
#%%Read specified DDMs
class DDM:
    def __init__(self,DDMfolder,SNR_Min=0,SNR_Max=50,Incidence_Min=0,Incidence_Max=90,
                 Gain_Min=-20,Gain_Max=30,SP_class='all'):
        
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
        firstgroup=metadata.groups[group_str[0]]
        
        TrackingOffsetDelayInPixels=np.zeros(N_groups)
        TrackingOffsetDopplerHz=np.zeros(N_groups)
        DopplerResolution=np.zeros(N_groups)
        #Attributes
        self.DelayResolution=firstgroup.DelayResolution
        self.DopplerResolution=firstgroup.DopplerResolution
        self.NumberOfDelayPixels=firstgroup.NumberOfDelayPixels
        self.NumberOfDopplerPixels=firstgroup.NumberOfDopplerPixels
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
            if(any(mask)):
                SNR=np.append(SNR,g_SNR[mask])
                Gain=np.append(Gain,g_Gain[mask])
                IncidenceAngle=np.append(IncidenceAngle,g_IncidenceAngle[mask])
                SPlat=np.append(SPlat,g_SPlat[mask])
                SPlon=np.append(SPlon,g_SPlon[mask])
                SpecularPathRangeOffset=np.append(SpecularPathRangeOffset,g_SpecularPathRangeOffset[mask])
                g_ddm=np.array(DDM_nc.groups[group_name].variables["DDM"])
                g_ddm=g_ddm[mask,:,:]
                DDMs_stack.append(g_ddm)
            else:
                pass
        Wave=np.concatenate(DDMs_stack,axis=0)
        self.Wave=Wave
        self.SNR=SNR
        self.IncidenceAngle=IncidenceAngle
        self.Gain=Gain
        self.SpecularPathRangeOffset=SpecularPathRangeOffset
        self.Lat=SPlat
        self.Lon=SPlon
        self.N=self.Wave.shape[0]
        if SP_class=='all':
            pass
        else:
            mask=SPwithin(SPlon,SPlat,SP_class)
            self.UpdateBymask(mask)

        metadata.close()
        DDM_nc.close()
        
    def ExecuteFilter(self,DEMmax=600):
        mask=DEMfilter(self.Lon,self.Lat,DEMmax)
        self.UpdateBymask(mask)
    
    def UpdateBymask(self,mask):
        self.Wave=self.Wave[mask,:,:]
        self.SNR=self.SNR[mask]
        self.IncidenceAngle=self.IncidenceAngle[mask]
        self.Gain=self.Gain[mask]
        self.SpecularPathRangeOffset=self.SpecularPathRangeOffset[mask]
        self.Lat=self.Lat[mask]
        self.Lon=self.Lon[mask]
        self.N=self.Wave.shape[0]
        
    def __repr__(self):
        return "Dataset of %d Delay Doppler Maps"%self.N
    
#%% Delay Doppler Maximum Average
    def DDMA(self,DopplerWinWidth=3,DelayWinWitdth=3,DopplerUnit='pixel',DelayUnit='pixel'):
        '''
        Calculte the Delay Doppler Map Average near the DDM centre
        DopplerWinWidth: Doppler window width
        DelayWinWitdth： width of Delay window
        '''
        if DopplerUnit=='pixel':
            w_Do=DopplerWinWidth//2
        elif DopplerUnit=='Hz':
            w_Do=DopplerWinWidth/self.DopplerResolution//2
        else:
            raise RuntimeError('Illegal parameter value in DopplerUnit : {}'.format(DopplerUnit))
        if DelayUnit=='ns':
            w_De=DelayWinWitdth/self.DelayResolution//2
        elif DelayUnit=='mus':
            w_De=DelayWinWitdth*1000/self.DelayResolution//2
        
        elif DelayUnit=='pixel':
            w_De=DelayWinWitdth//2
        else:
            raise RuntimeError('Illegal parameter value in Delay Unit : {}'.format(DelayUnit))
        DDMmax=np.where(self.Wave==65535)
        unique,ix=np.unique(DDMmax[0],return_index=True)
        Doc=DDMmax[1][ix]
        Dec=DDMmax[2][ix]
        Do_min=Doc-w_Do
        Do_max=Doc+w_Do+1
        De_min=Dec-w_De
        De_max=Dec+w_De+1
        ddma=[np.mean(self.Wave[i,Do_min[i]:Do_max[i],De_min[i]:De_max[i]]) for i in np.arange(self.N)]
        ddma=np.array(ddma)
        return ddma/65535
    
    #%%   DDM slice in Delay axis,Leading edge slope, and Trailing edge slope
    def DelaySlice(self,option='Interpolate'):
        if option=='Max':
            return np.max(self.Wave,axis=1)
        elif option=='Interpolate':
            return np.sum(self.Wave,axis=1)
        elif option=='Center':
            return self.Wave[:,self.NumberOfDopplerPixels//2,:]
        else:
            raise RuntimeError('Illegal option : {}'.format(option))
    

    def LES(self,DelayWinWitdth,DelayUnit='pixel',option='Interpolate'):
        if DelayUnit=='pixel':
            w=DelayWinWitdth
        ddmslice=self.DelaySlice(option=option)
        xlocate=np.arange(self.N)
        maxlocation=np.argmax(ddmslice,axis=1)
        maxlocation=np.where(maxlocation<w,
                             w,maxlocation)
        dif=ddmslice[xlocate,maxlocation]-ddmslice[xlocate,maxlocation-w]
        return dif/65535

    def TES(self,DelayWinWitdth,DelayUnit='pixel',option='Interpolate'):
        if DelayUnit=='pixel':
            w=DelayWinWitdth
        ddmslice=self.DelaySlice(option=option)
        xlocate=np.arange(self.N)
        maxlocation=np.argmax(ddmslice,axis=1)
        maxlocation=np.where(maxlocation>ddmslice.shape[1]-w-1,
                             ddmslice.shape[1]-w-1,maxlocation)
        dif=ddmslice[xlocate,maxlocation]-ddmslice[xlocate,maxlocation+w]
        return dif/65535
    
    def cor_WAF(self):
        '''
        The correlation coefficient with woodward ambiguty function
        '''
        waf=np.array([0,0,0,0.25,0.5,0.75,1,0.75,0.5,0.25,0,0,0])
        ddmslice=self.DelaySlice()
        w=6
        xlocate=np.arange(self.N)
        maxlocation=np.argmax(ddmslice,axis=1)
        maxlocation=np.where(maxlocation<w,w,maxlocation)
        maxlocation=np.where(maxlocation>ddmslice.shape[1]-w-1,
                             ddmslice.shape[1]-w-1,maxlocation)
        corr_waf=np.zeros(self.N)
        for i,j in list(zip(xlocate,maxlocation)):
            signal=ddmslice[i,j-w:j+w+1]
            corr_waf[i]=np.corrcoef(signal,waf)[0,1]
        return corr_waf
    
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
    source : 'land','lake','water',or shapefile path provided by user
    'water' is recommended to use due its less calculation time
    '''
    Lon=np.asarray(Lon)
    Lat=np.asarray(Lat)
    if source in ['water','land']:
        #downloaded from https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/georeferenced_tiff/
         watermask=np.load("data/WaterMask.npy")
         Lonindex=((Lon+180)*30).astype(np.int32)
         Latindex=((Lat+90)*30).astype(np.int32)
         if source=='land':
             return ~watermask[Latindex,Lonindex]
         else:
             return watermask[Latindex,Lonindex]
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
    

    return pts_in

def DEMfilter(Lon,Lat,DEMmax=600):
    dempath='data/dem_10km.npy'
    DEMglobe=np.load(dempath)
    Lonindex=((Lon+180)*10).astype(np.int32)
    Latindex=((Lat+90)*10).astype(np.int32)
    return DEMglobe[Latindex,Lonindex]<=DEMmax