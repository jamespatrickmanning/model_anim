# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 05:38:37 2022

@author: James.manning
"""
from matplotlib import pyplot as plt
import pandas as pd
#url='http://erddap.dataexplorer.oceanobservatories.org/erddap/tabledap/ooi-cp01cnsm-rid26-04-velpta000.csvp?time%2Clatitude%2Clongitude%2Ceastward_sea_water_velocity%2Cnorthward_sea_water_velocity&time%3E=2022-06-20T02%3A00%3A00Z'
url='http://erddap.dataexplorer.oceanobservatories.org/erddap/tabledap/ooi-cp03issm-rid26-04-velpta000.csvp?time%2Ceastward_sea_water_velocity%2Cnorthward_sea_water_velocity&time%3E=2022-06-20T02%3A00%3A00Z'
ooi=pd.read_csv(url)
ooi=ooi.set_index('time (UTC)')
ooi.index=pd.to_datetime(ooi.index)
#ooi=ooi.rename(columns={'eastward_sea_water_velocity (m.s-1)':'u','northward_sea_water_velocity (m.s-1)':'v'})

fig=plt.figure(figsize=(8,5))
ooi['northward_sea_water_velocity (m.s-1)'].plot()
plt.ylabel('northward velocity (m/s)')
fig.autofmt_xdate()
    