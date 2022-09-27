# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 12:46:40 2021
Routine to plot CFRF haul locations
@author: James.manning
"""
import pandas as pd
import os
os.environ['PROJ_LIB'] = 'C:\\Users\\james.manning\\Anaconda3\\pkgs\\proj-7.1.0-h7d85306_1\\Library\share'
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot as plt
import numpy as np

def getgbox(area):
  # gets geographic box based on area
  if area=='SNE':
    gbox=[-70.,-64.,39.,42.] # for SNE
  elif area=='OOI':
    gbox=[-72.,-69.5,39.5,41.5] # for OOI
  elif area=='GBANK':
    gbox=[-71.,-66.,40.,42.] # for GBANK
  elif area=='GBANK_RING':
    gbox=[-70.,-64.,38.,42.] # for typical GBANK Ring 
  elif area=='GS':           
    gbox=[-71.,-63.,38.,42.5] # for Gulf Stream
  elif area=='NorthShore':
    gbox=[-71.,-69.5,41.5,43.] # for north shore
  elif area=='IpswichBay':
    gbox=[-71.,-70.,42.5,43.] # for IpswitchBay
  elif area=='CCBAY':
    gbox=[-70.75,-69.8,41.5,42.23] # CCBAY
  elif area=='inside_CCBAY':
    gbox=[-70.75,-70.,41.7,42.15] # inside CCBAY
  elif area=='NEC':
    gbox=[-69.,-64.,39.,43.5] # NE Channel
  elif area=='NE':
    gbox=[-76.,-66.,35.,44.5] # NE Shelf 
  return gbox

files=['2014 to 2016','2016 to 2018']
area='NE'
'''
for k in range(len(files)):
    infile='CFRF_LobsterCrabResearchFleet_TempData_'+files[k]+'.xlsx'
    df1=pd.read_excel(infile)
    if k==0:
        df=df1
    else:
        df=pd.concat([df, pd.read_excel(infile)], axis=0)
'''
df=pd.read_excel('sample.xlsx')
gb=getgbox(area)
fig = plt.figure(figsize=(10,6))
m = Basemap(projection='stere',lon_0=(gb[0]+gb[1])/2.,lat_0=(gb[2]+gb[3])/2.,lat_ts=0,llcrnrlat=gb[2],urcrnrlat=gb[3],\
                    llcrnrlon=gb[0],urcrnrlon=gb[1],rsphere=6371200.,resolution='f',area_thresh=100)# JiM changed resolution to "c" for crude
        # draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.fillcontinents(color='gray',zorder=3)

ids=np.unique(df['id'])
# for each "id" session the header dataframe, plot data and save to the output .csv file
#for k in range(len(dfhead)):
for k in [ids[0]]:
    dfds=df[df['id']==k]
    x,y=m(dfds['longitude'][0],dfds['latitude'][0])
    m.scatter(x,y)
plt.show()