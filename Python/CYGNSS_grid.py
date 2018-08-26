# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from pandas import DataFrame,Series
def grid(location,resolution):
    #按照resolution的精度纬度进行分网格，如果对分辨率不满意可以进行修改
    return location//resolution*resolution
def GridCYGNSS(finput,foutput,resolution):
    #把CYGNSS的csv文件作为输入，输出到一个网格化的csv上
    f=open(finput)
    df=pd.read_csv(f)
    gp=df.groupby([grid(df["Lat"],resolution),grid(df["Lon"],resolution)])
    SR=gp["SR"].mean()
    SRdf=DataFrame(SR)
    SRdf.to_csv(foutput)
    f.close()
#示例使用，示例网格是0.05°，放在某个文件夹下迭代就可以了
fileintest="""E:\CYGNSS\csv\\2017230_china.csv"""
fileouttest="""E:\CYGNSS\csv\\2017230_china_mytest.csv"""
GridCYGNSS(fileintest,fileouttest,0.05)