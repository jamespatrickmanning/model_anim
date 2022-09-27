# -*- coding: utf-8 -*-
"""
@author: Vitalii Sheremet, FATE Project, 2012-2021, vsheremet@whoi.edu
"""

import numpy as np
from netCDF4 import Dataset        # netCDF4 version
import matplotlib.pyplot as plt
#import matplotlib.mlab as ml
from datetime import datetime

#Warning monthly mean files for years < 2007 have zero lon,lat values

URL0='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/'
# monthly mean files are here
#/thredds/dodsC/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/monthly_mean/gom3_monthly_mean_197801.nc
# catalog at 
#http://www.smast.umassd.edu:8080/thredds/catalog/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/catalog.html
#http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/catalog.html?dataset=models/NECOFS/Archive/Seaplan_33_Hindcast_v1/gom3_197801.nc

# sometimes fails 199407005959
# need to round to 1 hour
#tt=np.round(tRD*24.)/24.
#ti=datetime.fromordinal(int(tt))
#YEAR=str(ti.year)
#MO=str(ti.month).zfill(2)
#DA=str(ti.day).zfill(2)
#hr=(tt-int(tt))*24
#HR=str(int(np.round(hr))).zfill(2)            
#TS=YEAR+MO+DA+HR+'0000'
#FN2=TS
#T1=datetime.strftime(t1,'%Y%m%dT%H:%M') # '%Y%m%dT%H:%M'
LON1=-70.;LON2=-65.+1.e-6;LAT1=39.0;LAT2=42.+1.e-6# JiM added this up front to limit results
ss=5 # JiM's subsampling rate to quiver 
scale=5 # Jim's scale of quiver where it was originally 10

YEARMOS=np.array([])
for yr in range(1978,2017):
#for yr in range(2010,2011):
#for yr in [2016]:
#    for mo in range(1,1+12):
    for mo in [4]:
        YEARMO=datetime.strftime(datetime(yr,mo,1),'%Y%m') # # '%Y%m%dT%H:%M'
        YEARMOS=np.append(YEARMOS,YEARMO)

#YEARMOS=np.array(['197904']) # test single file 
YEARMOS=np.array(['201409']) # test single file  phenomenal Haddock year class

for k in range(0,len(YEARMOS)):
    YEARMO=YEARMOS[k]
    print (k,YEARMO) 
    # 1h fields
    #http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/gom3_197801.nc
    URL='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/gom3_'+YEARMO+'.nc'

# monthly mean files are here
#/thredds/dodsC/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/monthly_mean/gom3_monthly_mean_197801.nc
    URL='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/monthly_mean/gom3_monthly_mean_'+YEARMO+'.nc'

    #xxx=ds['xxx']; np.save('gom3.xxx.npy',np.array(xxx))
    ds = Dataset(URL,'r').variables   # netCDF4 version
    time=np.array(ds['time'][:])
#    lon=np.array(ds['lon'][:])
#    lat=np.array(ds['lat'][:])
#    lonc=np.array(ds['lonc'][:])
#    latc=np.array(ds['latc'][:])
    us=np.array(ds['u'][:,0,:]) # surface current
    vs=np.array(ds['v'][:,0,:])
#    ub=np.array(ds['u'][:,-1,:]) # bottom current
#    vb=np.array(ds['v'][:,-1,:])
    ua=np.array(ds['ua'][:,:]) # depth avg current
    va=np.array(ds['va'][:,:])
#    uwind_stress=np.array(ds['uwind_stress'][:,:]) # wind stress
#    vwind_stress=np.array(ds['vwind_stress'][:,:])
#    zeta=np.array(ds['zeta'][:,:])
#    temps=np.array(ds['temp'][:,0,:]) # surface temperature
#    tempb=np.array(ds['temp'][:,-1,:]) # bottom temperature


# All relevant data for dynamical analysis    
#    D={'time':time,'lon':lon,'lat':lat,'lonc':lonc,'latc':latc,
#'zeta':zeta,'ua':ua,'va':va,'us':us,'vs':vs,'ub':ub,'vb':vb,
#'uwind_stress':uwind_stress,'vwind_stress':vwind_stress,
#'temps':temps,'tempb':tempb
#}
    FOLDER='DATA/' 
    FOLDER='' #JiM
    FN=FOLDER+'gom3_time_monthly_mean_'+YEARMO+'.npy'
    np.save(FN,time)
    FN=FOLDER+'gom3_ua_monthly_mean_'+YEARMO+'.npy'
    np.save(FN,ua)
    FN=FOLDER+'gom3_va_monthly_mean_'+YEARMO+'.npy'
    np.save(FN,va)
    FN=FOLDER+'gom3_us_monthly_mean_'+YEARMO+'.npy'
    np.save(FN,us)
    FN=FOLDER+'gom3_vs_monthly_mean_'+YEARMO+'.npy'
    np.save(FN,vs)

# D=np.load(FN).item() # to load the dictionary, D['time'] to refer to time array 

# uv only for dtr
#    D={'time':time,'lon':lon,'lat':lat,'lonc':lonc,'latc':latc,
#'zeta':zeta,'ua':ua,'va':va,'us':us,'vs':vs,'ub':ub,'vb':vb,
#'uwind_stress':uwind_stress,'vwind_stress':vwind_stress,
#}
#    FNOUT='gom3_uvz_'+YEARMO+'.npy'
#    np.save(FNOUT,D)

# temp only for front analysis
#    D={'time':time,'temps':temps,'tempb':tempb}
#    FNOUT='gom3_tempsb_'+YEARMO+'.npy'
#    np.save(FNOUT,D)


    # test plot velocity field
    lonc=np.array(ds['lonc'][:])
    latc=np.array(ds['latc'][:])
    lon=np.array(ds['lon'][:])
    lat=np.array(ds['lat'][:])
    h=np.array(ds['h'][:])

    #Warning monthly mean files for years < 2007 have zero lon,lat values
    # load lon,lat from preloded file
    Grid=np.load('Grid.npy',allow_pickle=True).item()
    lon=Grid['lon']
    lat=Grid['lat']
    lonc=Grid['lonc']
    latc=Grid['latc']


    u=ua;v=va
    fig=plt.figure(figsize=(15,10))
    A=1./np.cos(41*np.pi/180.)
    plt.tricontour(lon,lat,-h,[-200,-100,-70])  
    # shows problematic area SouthEast of the Northeast Channel with grid scale spurious oscillations in ua,ua
    #plt.quiver(lonc,latc,u,v,scale=4); LON1=-66.;LON2=-64.5+1.e-6;LAT1=41.5;LAT2=42.5+1.e-6 
    plt.quiver(lonc[::ss],latc[::ss],u[::ss,::ss],v[::ss,::ss],scale=scale); #LON1=-71.;LON2=-64.+1.e-6;LAT1=39.0;LAT2=44.+1.e-6
    
    plt.axis([LON1,LON2,LAT1,LAT2])
    ax=plt.gca()
    ax.set_aspect(A)
    plt.title(YEARMO+' mean flow with vectors subsampled by '+str(ss),fontsize=16)
    plt.show()
    fig.savefig('surface_mean_flow_'+YEARMO+'.png')

