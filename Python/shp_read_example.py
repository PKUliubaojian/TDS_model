# -*- coding: utf-8 -*-
"""
Created on Sun Aug 26 00:27:59 2018

@author: lenovo
"""
import shapefile
import geopandas as gpd
from shapely.geometry import Point,Polygon
path_land="shp\\ne_110m_land\\ne_110m_land.shp"
path_lake="shp\\ne_10m_lakes\\ne_10m_lakes.shp"
sr=shapefile.Reader(path_lake)
shapes=sr.shapes()

s_land=gpd.read_file(path_land)
s_lake=gpd.read_file(path_lake)
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
world.plot()
pt=[Point(2,2),Point(2,3),Point(3,3)]
pts=gpd.GeoDataFrame(geometry=pt)
polygon1=Polygon([(2.5,1),(1,3),(3,3),(3,1)])
polygon2=Polygon([(-1,1),(1,3),(3,3),(3,1)])
pys=gpd.GeoDataFrame(geometry=[polygon1,polygon2])
union=pts["geometry"].unary_union
pts.geometry.within(pys.geometry)
