# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 00:27:59 2018

@author: lenovo
"""

import geopandas as gpd
from shapely.geometry import Point,Polygon
import numpy as np
from pandas import DataFrame,Series


def pixelswithin(Lon,Lat,source="land"):
    '''Points in shapeflie'''
    if source in ['land']:
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
    return pts_in
        