# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from pandas import DataFrame
import geopandas as gpd
from shapely.geometry import Point
import fiona.drvsupport
def grid(location,resolution):
    #按照resolution的精度纬度进行分网格，如果对分辨率不满意可以进行修改
    return location//resolution*resolution
def GridCYGNSS(finput,foutput,resolutionX,resolutionY):
    '''
    把CYGNSS的csv文件作为输入，输出到一个网格化的csv上
    '''
    f=open(finput)
    df=pd.read_csv(f)
    gp=df.groupby([grid(df["Lat"],resolutionX),grid(df["Lon"],resolutionY)])
    SR=gp["SR"].mean()
    SRdf=DataFrame(SR)
    SRdf.to_csv(foutput)
    f.close()
#示例使用，示例网格是0.05°*0.05°，放在某个文件夹下迭代就可以了

def Grid_Shp(finput,foutput,resolutionX,resolutionY):
    '''
    Gridding the data using a specified resolution
    Usage: Grid_Shp(finput,foutput,resolutionX,resolutionY)
        finput:     input file name(string)
        fonput:     input file name(string)
        resolutionX:     resolution of Longitude(m)
        resolutionX:     resolution of Latitude(m)
    '''
    df=gpd.read_file(finput)
    gp=df.groupby([grid(df["Lat"],resolutionX),grid(df["Lon"],resolutionY)])
    SR=gp["SR"].mean()
    SRdf=DataFrame(SR)
    SRdf["Geo"]=SRdf.index
    SRdf["Geo"]=SRdf["Geo"].apply(Point)
    gdf = gpd.GeoDataFrame(SRdf, geometry='Geo')
    gdf.to_file(foutput)
 
    
    
def Grid_boundary_Shp(finput,foutput,boundary_source,resolutionX,resolutionY):
    '''
    Gridding the data using a specified resolution
    Usage: Grid_Shp(finput,foutput,resolutionX,resolutionY)
        finput:     input file name(string)
        fonput:     input file name(string)
        resolutionX:     resolution of Longitude(m)
        resolutionX:     resolution of Latitude(m)
    '''
    df=gpd.read_file(finput)
    gp=df.groupby([grid(df["Lat"],resolutionX),grid(df["Lon"],resolutionY)])
    SR=gp["SR"].mean()
    SRdf=DataFrame(SR)
    SRdf["Geo"]=SRdf.index
    SRdf["Geo"]=SRdf["Geo"].apply(Point)
    gdf = gpd.GeoDataFrame(SRdf, geometry='Geo')
    boundary=gpd.read_file(boundary_source)
    union=boundary["geometry"].unary_union
    basinin=gdf.geometry.within(union)
    gdf=gdf.drop(~basinin)
    gdf.to_file(foutput)