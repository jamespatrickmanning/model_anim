# -*- coding: utf-8 -*-
"""
Created on Fri May 24 13:58:18 2013

@author: jmanning
routine to generate a vector plot from user-specified time and place
reads a control file for input parameters

modified in Dec 2015 to test on both eMOLT and COMET machine
revisited in Dec 2018 to overlay on drifter tracks in CC Bay
revisited in Oct 2019 to make detided vectors using sh_rmtide function
"""

from pylab import *
from matplotlib.collections import PolyCollection
import matplotlib.tri as Tri
from mpl_toolkits.basemap import Basemap
import datetime as dt
import netCDF4
import sys
import numpy as np
from datetime import timedelta
#from pydap.client import open_url

# HARDCODES ##########
ctrl_filename='ctrl_uvmodel_ccbay_2014_11_14.csv' # reads some HARDCODES listed in special control file based on project 

#############################################

# read input variables from control file
urlname=open(ctrl_filename, "r").readlines()[0][58:-1]
depth=int(open(ctrl_filename, "r").readlines()[1][22:-1])
TIME=open(ctrl_filename, "r").readlines()[2][31:-1]
lon_bin_size=float(open(ctrl_filename, "r").readlines()[3][32:-1])
lat_bin_size=float(open(ctrl_filename, "r").readlines()[4][32:-1])
gbox=array(open(ctrl_filename, "r").readlines()[5][5:-1].split(','),float)
detide=open(ctrl_filename, "r").readlines()[6][7:-1]

def inconvexpolygon(xp,yp,xv,yv):
    """
	check if point is inside a convex polygon

    i=inconvexpolygon(xp,yp,xv,yv)
    
    xp,yp - arrays of points to be tested
    xv,yv - vertices of the convex polygon

    i - boolean, True if xp,yp inside the polygon, False otherwise
    
    """    
    N=len(xv)   
    j=np.arange(N)
    ja=(j+1)%N # next vertex in the sequence    
    NP=len(xp)
    i=np.zeros(NP,dtype=bool)
    for k in range(NP):
        # area of triangle p,j,j+1
        Aj=np.cross(np.array([xv[j]-xp[k],yv[j]-yp[k]]).T,np.array([xv[ja]-xp[k],yv[ja]-yp[k]]).T) 
    	# if a point is inside the convect polygon all these Areas should be positive 
    	# (assuming the area of polygon is positive, counterclockwise contour)
        Aj /= Aj.sum()
    	# Now there should be no negative Aj
    	# unless the point is outside the triangular mesh
        i[k]=(Aj>0.).all()
        
    return i

def sh_rmtide(f,dt=1.,ends=0.):
    """
    removes solar and lunar tidal signal by sequentially averaging 
    over their periods: of 24h and 24h50.47m. This is equivalent to
    applying convolution with a trapezoidal shaped window (almost triagular).
    
    f - input tidal sequence
    dt - uniform time step of input series, in hours, dt=1. hr default
    ends - fill value for the ends: 0. or nan
    
    fm = sh_rmtide(f,dt=1.,ends=np.nan)
    """
    TS=24. # Solar hour angle period 24h
    TL=24.+50.47/60. # Lunar hour angle period
    N=len(f)
    fm=np.zeros(N)*ends

    # remove solar period    
    T=TS/dt
    m=int(np.floor((T-1.)/2.))
    w=(T-(2*m+1))*0.5
    # all weights =1 except at ends      
    for i in range(m+1,N-m-1-1):
        #print len(f[i-m:i+m+1])
        fm[i]=(np.sum(f[i-m:i+m+1])+w*f[i-m-1]+w*f[i+m+1])/T
        
    # remove lunar period
    f=fm*1.0  # deep copy!  
    T=TL/dt
    m=int(np.floor((T-1.)/2.))
    w=(T-(2*m+1))*0.5
    # all weights =1 except at ends      
    #    print w
    for i in range(m+1,N-m-1-1):
        #print len(f[i-m:i+m+1])
        fm[i]=(np.sum(f[i-m:i+m+1])+w*f[i-m-1]+w*f[i+m+1])/T
    
    return fm

def sh_bindata(x, y, z, xbins, ybins):
    ix=np.digitize(x,xbins)
    iy=np.digitize(y,ybins)
    xb=0.5*(xbins[:-1]+xbins[1:]) # bin x centers
    yb=0.5*(ybins[:-1]+ybins[1:]) # bin y centers
    zb_mean=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    for iix in range(1,len(xbins)):
        for iiy in range(1,len(ybins)):
            k,=np.where((ix==iix) & (iy==iiy))
            zb_mean[iix-1,iiy-1]=np.mean(z[k])
    return xb,yb,zb_mean

#MAIN CODE 
num_days=2#number of days to make picture
if urlname[0:7]=="massbay":
    TIME=dt.datetime.strptime(TIME, "%Y-%m-%d %H:%M:%S") 
    now=dt.datetime.now()

    timeperiod=(TIME)-(now-timedelta(days=3))
    startrecord=(timeperiod.seconds)/60/60
    if urlname=='massbay_forecast':
      if TIME>now:
         diff=(TIME-now).days
      else:
         diff=(now-TIME).days
      if diff>3:
        print("please check your input start time,within 3 days both side form now on")
        sys.exit(0)
      url='http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc?lon,lat,lonc,latc,time,nv,h,siglay,v,u'   
    elif urlname=='massbay_archives': # case of archived forecasts
      #year=raw_input('What year? ')
      #mth=raw_input('What month? ')
      year=TIME.year
      mth=TIME.month
      url='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/NECOFS_MASS_BAY/'+str(year)+'/mbn_'+str(year)+str(mth)+'.nc'

if urlname=="30yr":
      stime=dt.datetime.strptime(TIME, "%Y-%m-%d %H:%M:%S")
      timesnum=stime.year-1981
      standardtime=dt.datetime.strptime(str(stime.year)+'-01-01 00:00:00', "%Y-%m-%d %H:%M:%S")
      timedeltaprocess=(stime-standardtime).days
      startrecord=26340+35112*(timesnum/4)+8772*(timesnum%4)+1+timedeltaprocess*24     
      url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3?lon,lat,lonc,latc,time,nv,h,siglay,v,u'

print(urlname)
nc = netCDF4.Dataset(url)
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
latc = nc.variables['latc'][:]
lonc = nc.variables['lonc'][:]
h = nc.variables['h'][:]
siglay=nc.variables['siglay']
u= nc.variables['u']
v= nc.variables['v']
nv = nc.variables['nv'][:].T - 1
print('we have the model data and now we are finding indices within the gbox')
gbox_poly=array([[gbox[0],gbox[2]],[gbox[1],gbox[2]],[gbox[1],gbox[3]],[gbox[0],gbox[3]]])
i=inconvexpolygon(lonc,latc,gbox_poly[:,0],gbox_poly[:,1])
i=np.argwhere(i).flatten()
llond0=lonc[i] # reduces lat/lon to inside the polygon
llatd0=latc[i]
num_model_pts_in_gbox=len(i)
print('number of model points inside gbox = ',num_model_pts_in_gbox)
#THIS ONLY CODED FOR THE DETIDED CASE!!!!!!!!!!!!!!!!!!!!!!!
for n in range(num_days): # though through one day at a time (increments by 24 hours at end of loop)
    utotal=[0]*num_model_pts_in_gbox
    vtotal=[0]*num_model_pts_in_gbox
    utotal=np.zeros(num_model_pts_in_gbox,dtype=np.real)
    vtotal=np.zeros(num_model_pts_in_gbox,dtype=np.real)

    if depth==-1: # case of surface flow
        if detide=='no':
            utotal=u[startrecord,0,:]
            vtotal=v[startrecord,0,:]
        else:
            print('detiding')
        for k in range(num_model_pts_in_gbox): # for each velocity inside the polygon, detide
	      		utotal[k]=sh_rmtide(u[startrecord-24:startrecord+24,0,i[k]])[24]# uses 48 hours to apply a filter around both solar and lunar periods
	      		vtotal[k]=sh_rmtide(v[startrecord-24:startrecord+24,0,i[k]])[24]
    else:
        for i in range(len(lon)):
            depthtotal=siglay[:,i]*h[i]
            layer=np.argmin(abs(depthtotal+depth))
            utotal.append(u[startrecord,layer,i])
            vtotal.append(v[startrecord,layer,i])
            utotal=np.array(utotal)
            vtotal=np.array(vtotal)

    print('now lets bin the data')
    xi = np.arange(gbox[0]-0.1,gbox[1]+0.1,lon_bin_size)
    yi = np.arange(gbox[2]-0.1,gbox[3]+0.1,lat_bin_size)
    xb,yb,ub_mean = sh_bindata(llond0[::-1], llatd0[::-1], utotal, xi, yi)
    xb,yb,vb_mean = sh_bindata(llond0[::-1], llatd0[::-1], vtotal, xi, yi)
    xxb,yyb = np.meshgrid(xb, yb)
    latsize=[gbox[2]-0.1,gbox[3]+0.1]
    lonsize=[gbox[0]-0.1,gbox[1]+0.1]

    print('and plot')
    plt.figure()
    m = Basemap(projection='cyl',llcrnrlat=min(latsize),urcrnrlat=max(latsize),\
            llcrnrlon=min(lonsize),urcrnrlon=max(lonsize),resolution='h')#,fix_aspect=False)
    m.drawparallels(np.arange(int(min(latsize)),int(ceil(max(latsize))),round(diff(latsize)[0]/4,1)),labels=[1,0,0,0])
    m.drawmeridians(np.arange(int(min(lonsize)),int(ceil(max(lonsize))),round(diff(lonsize)[0]/4,1)),labels=[0,0,0,1])
    #m.drawparallels(np.arange(min(latsize),max(latsize),5),labels=[1,0,0,0])
    #m.drawmeridians(np.arange(min(lonsize),max(lonsize),10),labels=[0,0,0,1])
    m.drawcoastlines()
    m.fillcontinents(color='grey')
    m.drawmapboundary()
    ub = np.ma.array(ub_mean, mask=np.isnan(ub_mean))
    vb = np.ma.array(vb_mean, mask=np.isnan(vb_mean))
    Q=m.quiver(xxb,yyb,ub,vb,scale=1)

    plt.quiverkey(Q,0.8,0.05,.2, '0.2 m/s', labelpos='W')
    plt.title(urlname+' Depth:'+str(depth)+' Time:'+str(TIME)) 
    plt.show()
    print('Saved figure as <grid>_vectors<TIME>')
    plt.savefig(urlname+'_vectors_'+str(TIME)[0:10]+'_'+str(TIME)[11:16]+'_binned_to_'+str(lat_bin_size)+'_deg_'+str(n)+'.png')
    #plt.close('all')
    startrecord=startrecord+24# jumps one day
