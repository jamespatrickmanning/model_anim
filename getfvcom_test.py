# routine to test "getfvcom" function
from datetime import datetime as dt
from datetime import timedelta
from dateutil import parser
import pandas as pd
import numpy as np
import netCDF4
import sys
from conversions import c2f,dd2dm
lati=sys.argv[1]
loni=sys.argv[2]
depth=sys.argv[3]
yr=int(sys.argv[4])
mth=int(sys.argv[5])
day=int(sys.argv[6])
hr=int(sys.argv[7])#99999 for bottom
mn=int(sys.argv[8])
se=0#seconds

def get_FVCOM_temp(lati,loni,dtime,depth): # gets modeled temp using GOM3 forecast
        '''
        Taken primarily from Rich's blog at: http://rsignell-usgs.github.io/blog/blog/2014/01/08/fvcom/ on July 30, 2018
        where lati and loni are the position of interest, dtime is the datetime, and layer is "-1" for bottom
        '''
        #urlfvcom = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
        #urlfvcom = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc'
        urlfvcom = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'
        nc = netCDF4.Dataset(urlfvcom).variables
        #first find the index of the grid 
        lat = nc['lat'][:]
        lon = nc['lon'][:]
        inode = nearlonlat(lon,lat,loni,lati)
        #second find the index of time
        time_var = nc['time']
        itime = netCDF4.date2index(dtime,time_var,select='nearest')# where startime in datetime
        # figure out layer from depth
        w_depth=nc['h'][inode]
        layer=int(round(depth/w_depth*45.))
        return nc['temp'][itime,layer,inode]

dtime=dt(yr,mth,day,hr,mn,se)
temp=get_FVCOM_temp(lati,loni,dtime,depth)
print temp
