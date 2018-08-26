# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 00:27:59 2018

@author: lenovo
"""

import geopandas as gpd
from shapely.geometry import Point,Polygon
import pandas
from pandas import DataFrame,Series


def pixelswithin(points,source="land"):
    '''Points in shapeflie'''
    if source == 'land':
        path="shp\\ne_110m_land\\ne_110m_land.shp"

    elif source=='lake':
        path="shp\\ne_10m_lakes\\ne_10m_lakes.shp"
    else:
        path=source
    shapes=gpd.read_file(path)
    pts=gpd.GeoDataFrame(geometry=points)
    union=shapes["geometry"].unary_union
    pts_in=pts.geometry.within(union)
    return pts_in
        