#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Apr 17 11:29:14 2019
@author: leizhao
Modified by JiM in Fall 2019
Modified by JiM in Summer 2020 with Jack Polentes help
Modified by JiM in Fall 2020 by adding drifter tracks and, in late Oct 2020, NCEP wind vector
Modified by JiM in Spring 2021 to add all model options including FVCOM
Modified by JiM in Fall 2021 to add OOI variables 
Stored in Github as "model_anim" repository
See hardcodes at top of code where there is, for example, "start_date" and "ndays" 
Note: You may need to adjust the "clevs" to get good range of colors.
Note: You may want to adjust the Basemap resolution to "c" for crude in experiments.
Note: You may need to "conda install -c conda-forge basemap-data-hires" in order to use the higher resolution coastlines
Note: 
"""
#hardcodes########################
area='SNE'#'SNW'#'NorthShore'#'SNW'#'GBANK_RING'#'Gloucester'
start_date='2022-08-21'#'2013-04-01'
#clevs=[65.,80.,.5]#gb ring surf June
clevs=[50.,72.,.5]#ns bottom June
clevs=[52.,78.,.5]#gb ring June
clevs=[58.,74.,.5]#SNE-W in July
clevs=[58.,80.,.5]#SNE in August
dtime=[]
units='degF'
ndays=2 # number of days to run
detide='n'# for FVCOM hindast only now
include_temp_obs='yes' # 'yes' overlays positions of observed bottom temps
include_temp_obs_LFA='no' # 'yes' overlays LFA positions of observed bottom temp
include_wind='no'
include_ooi='yes'
include_weather_balloon='no'
#cluster='ep_2022_1' # code of drifter cluster to include
cluster='shp_2022_1'
surf_or_bot=-1#0 for most but -1 for surface (opposite for FVCOM)
lat_w,lon_w=40.,-68.5 # base of wind vector legend (actual vector appears mid plot)
model='DOPPIO'#'DOPPIO'# FVCOM, DOPPIO, or GOMOFS ...
#########
import os,imageio
import conda
import pandas as pd
from scipy.interpolate import griddata as gd
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
#os.environ['PROJ_LIB'] = 'c:\\Users\\Joann\\anaconda3\\pkgs\\proj4-5.2.0-ha925a31_1\\Library\share'
os.environ['PROJ_LIB'] = 'C:\\Users\\james.manning\\Anaconda3\\pkgs\\proj-7.1.0-h7d85306_1\\Library\share'
from mpl_toolkits.basemap import Basemap
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
import time
import zlconversions as zl
import gomofs_modules
import sys
import warnings
warnings.filterwarnings("ignore") # gets rid of warnings at runtime but you may want to comment this out to see things
import netCDF4 # for wind access
from math import sqrt
try:
    import cPickle as pickle
except ImportError:
    import pickle
import glob
import ssl # added this in Oct 2021 when I got a "certificate_verify_failed" message
ssl._create_default_https_context = ssl._create_unverified_context

def getgbox(area):
  # gets geographic box based on area
  if area=='SNE':
    gbox=[-72.,-66.,39.,42.] # for SNE shelf east
  elif area=='SNW':
    gbox=[-71.5,-69.5,40.,41.75] # for SNw shelf west
  elif area=='MABN':
    gbox=[-73.,-68.,39.,42.] # for SNw shelf west  
  elif area=='OOI':
    gbox=[-71.5,-69.,39.5,41.6] # for OOI
  elif area=='GBANK':
    gbox=[-71.,-66.,40.,42.] # for GBANK
  elif area=='GBANK_RING':
    gbox=[-71.,-65.,39.,42.5] # for typical GBANK Ring 
  elif area=='GS':           
    gbox=[-71.,-63.,38.,42.5] # for Gulf Stream
  elif area=='NorthShore':
    gbox=[-71.,-69.5,41.75,43.25] # for north shore
  elif area=='Gloucester':
    gbox=[-71.,-70.,42.25,43.] # for north shore
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


#def get_wind_ncep(starttime,endtime,lat,lon):
def get_wind_ncep(year,lat,lon):    
        #function get a time series of u & v wind m/s from pre-downloaded ncep file
        #Note: there is a new version of this function inside "jamespatrickmanning/stick_plot" github repository in get_stick_plot.py
        url_input=""
        url_uwind=url_input+'uwnd.sig995.'+str(year)+'.nc'## where I downloaded these from ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/surface/
        url_vwind=url_input+'vwnd.sig995.'+str(year)+'.nc'
        #print(url_vwind)
        ncv=netCDF4.Dataset(url_vwind)
        ncu=netCDF4.Dataset(url_uwind)
        
        t=ncv.variables['time'][:]
        u_wind=ncu.variables['uwnd']
        v_wind=ncv.variables['vwnd']
        
        LAT=ncu.variables['lat'][:]
        LON=ncu.variables['lon'][:]
        
        for i in range(len(LON)):#transfer lon from (0,360) to (-180.180)
            if(LON[i]>180):
                LON[i]=-360+LON[i]
        ###########find index of nearest point###########
        index=[]
        d=[]
        for a in np.arange(len(LAT)):
            d1=[]
            for b in np.arange(len(LON)):
                d2=sqrt((LAT[a]-lat)*(LAT[a]-lat)+(LON[b]-lon)*(LON[b]-lon))
                d1.append(d2)
            d.append(d1) 
        #print(np.argmin(d)/len(LON), np.argmin(d)%len(LON),np.hstack(d)[np.argmin(d)],d[np.argmin(d)/len(LON)][np.argmin(d)%len(LON)]
        index.append(np.argmin(d)/len(LON))#index of LAT
        index.append(np.argmin(d)%len(LON))#index of LON
        #print(index
        cptime="%i,01,01,00,00"  %year
        cptimes=datetime.strptime(cptime, '%Y,%m,%d,%H,%M')
        moddate=[]
        for k in range(len(t)):
            moddate.append(cptimes+timedelta(hours=t[k]-t[0]))
        du= u_wind[:,index[0],index[1]].tolist()
        dv= v_wind[:,index[0],index[1]].tolist()
        du=list(map(float,du))
        dv=list(map(float,dv))
        df_w=pd.DataFrame([moddate,du,dv])
        #df_w=pd.DataFrame([moddate,u_wind[:,index[0],index[1]].tolist(),v_wind[:,index[0],index[1]].tolist()])  
        df_w=df_w.T.set_index(0) # takes teh transpose and then makes datetime column the index
        df_w.columns=['u','v']
        dfu=pd.to_numeric(df_w.u)
        dfuh=dfu.resample('H').mean().interpolate('linear') 
        dfv=pd.to_numeric(df_w.v)
        dfvh=dfv.resample('H').mean().interpolate('linear')
        df_wf=pd.concat([dfuh,dfvh],axis=1)
        # make it hourly
        return df_wf

def plotit(lons,lats,slons,slats,stemp,temp,depth,time_str,path_save,dpi=80,area='OOI',clevs=[39.,44.,0.5],lat_w=42.5,lon_w=-70.5,wind=pd.DataFrame([[datetime.now()-timedelta(1),datetime.now()],[.1,.1],[.1,.1]]),ooi=pd.DataFrame([[datetime.now()-timedelta(1),datetime.now()],[.1,.1],[.1,.1]]),dtime=dtime,model='DOPPIO',detide='n',ooiw=pd.DataFrame([[datetime.now()-timedelta(1),datetime.now()],[.1,.1],[.1,.1]])):
    '''dtimes,u,v are wind in m/s'''
    if clevs[1]<32.:
        units='degC'
    else:
        units='degF'
    fig = plt.figure(figsize=(10,6))
    #ax = fig.add_axes([0.01,0.05,0.98,0.87])
    ax = fig.add_axes([0.05,0.05,0.9,0.82])
    # create polar stereographic Basemap instance.
    gb=getgbox(area)
    m = Basemap(projection='stere',lon_0=(gb[0]+gb[1])/2.,lat_0=(gb[2]+gb[3])/2.,lat_ts=0,llcrnrlat=gb[2],urcrnrlat=gb[3],\
                llcrnrlon=gb[0],urcrnrlon=gb[1],rsphere=6371200.,resolution='f',area_thresh=100)# JiM changed resolution to "c" for crude
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.fillcontinents(color='gray',zorder=3)
    if include_wind=='yes':
       xw,yw=m(lon_w,lat_w)  
       ax.quiver(xw,yw,10.0,0.0,color='red',scale=50.,zorder=2)#wind legend in middle
       ax.text(xw,yw+(m.ymin+m.ymax)/20,'10 m/s (~20 knots) NCEP wind',fontsize=16,color='red',zorder=2) 
       idex=wind.index.get_loc(datetime.strptime(time_str[0:15],"%Y-%m-%d %H%M"),method='nearest')
       ax.quiver((m.xmin+m.xmax)/2,(m.ymin+m.ymax)/2,wind.u[idex],wind.v[idex],scale=50.0,color='red',zorder=2) # plot an arrow
    if include_ooi=='yes':
       #xo,yo=m(-70.77,39.943)  
       xo,yo=m(-70.78,40.13)
       ax.quiver(xo,yo+(m.ymax-m.ymin)/20,.1,0.0,color='b',scale=2.,zorder=2)#current legend in middle
       ax.text(xo,yo+(m.ymax-m.ymin)/15,'.1 m/s (~1/4 knot current) ',fontsize=14,color='b',zorder=2) 
       ooi.index=pd.to_datetime(ooi.index)
       idex=ooi.index.get_loc(datetime.strptime(time_str[0:15],"%Y-%m-%d %H%M"),method='nearest')
       ax.quiver(xo,yo+(m.ymax-m.ymin)/40,ooi.u[idex],ooi.v[idex],scale=2.,color='b',zorder=2) # plot an arrow
       #wind
       xo,yo=m(-70.78,39.7)
       ax.quiver(xo,yo+(m.ymax-m.ymin)/20,10.,0.0,color='grey',scale=100.,zorder=2)#current legend in middle
       ax.text(xo,yo+(m.ymax-m.ymin)/15,'10.0 m/s (~20 knot wind) ',fontsize=14,color='grey',zorder=2) 
       ooiw.index=pd.to_datetime(ooiw.index)
       idex=ooiw.index.get_loc(datetime.strptime(time_str[0:15],"%Y-%m-%d %H%M"),method='nearest')
       ax.quiver(xo,yo,ooiw.uw[idex],ooiw.vw[idex],scale=100.,color='grey',zorder=2) 
    if include_temp_obs_LFA=='yes':# special case where you have a set of points to post as in CC Bay 2020 study here
        dfLFoM=pd.read_csv('C:\\Users\Joann\Downloads\LFoM\LFoM_sites.csv')
        xx,yy=m(list(dfLFoM['lon']),list(dfLFoM['lat']))
        ax.plot(xx,yy,'ro')
    if len(slons)!=0:
        x1,y1=m(slons,slats)
        for jj in range(len(stemp)):
            ax.text(x1[jj],y1[jj],str.format('{0:.1f}',stemp[jj]*1.8+32),color='w',fontsize=10,fontweight='bold',horizontalalignment='center',verticalalignment='center')
    # draw parallels.
    if (area=='IpswichBay') or (area=='Gloucester'):
        labint=0.2
        dept_clevs=[30,50,100, 150]
    elif area[0:3]=='CCB':
        labint=0.5
        dept_clevs=[30,50,100]
    elif (area=='NE') or (area=='SNW'):
        labint=1.0
        dept_clevs=[50,100,1000]
        x,y=m(-69.5,40.5)    
        #plt.text(x,y,'Great South Channel',fontsize=12, rotation=0) 

    elif (area=='NorthShore'):
        labint=.50
        dept_clevs=[50,100,150]  
    elif (area=='GBANK') or (area=='GBANK_RING'):
        labint=1.0
        dept_clevs=[50,100,150]
        x,y=m(-68.3,40.65)    
        plt.text(x,y,' Georges Bank',fontsize=16, rotation=30) 
    else:
        labint=1.0
        dept_clevs=[30,50,100, 150,300,1000]
    parallels = np.arange(0.,90,labint)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=20)
    # draw meridians
    meridians = np.arange(180.,360.,labint)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=20)
    x, y = m(lons, lats) # compute map proj coordinates.
    dtthis=datetime(int(time_str[0:4]),int(time_str[5:7]),int(time_str[8:10]),int(time_str[11:13]),int(time_str[13:15]))
    clevs=np.arange(clevs[0],clevs[1],clevs[2])  #for all year:np.arange(34,84,1) or np.arange(34,68,1)
    if (model=='DOPPIO') or (model=='GOMOFS'):
            # draw filled contours.
            cs = m.contourf(x,y,temp,clevs,cmap=plt.get_cmap('rainbow'))
            # draw depth contours.
            dept_cs=m.contour(x,y,depth,dept_clevs,colors='black')
            plt.clabel(dept_cs, inline = True, fontsize =15,fmt="%1.0f")

    elif model=='FVCOM':
            if dtthis<datetime(2020,6,1,0,0,0):
                if detide=='y': # here we use the seaplan daily averages
                    url='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/daily_mean/gom3_daily_mean_'+str(dtthis.year)+str(dtthis.month).zfill(2)+'.nc' 
                else:
                    url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
            elif dtthis>datetime.now()-timedelta(days=2):    
                url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'
            else:
                print('FVCOM not available for this date/time.')
            #else:
            	
            nc = netCDF4.Dataset(url).variables
            time_var = nc['time']
            #itime = netCDF4.date2index(dtime,time_var,select='nearest')
            itime = netCDF4.date2index(dtthis,time_var,select='nearest')
            # Get lon,lat coordinates for nodes (depth)
            lats = nc['lat'][:]
            lons = nc['lon'][:]
            # Get lon,lat coordinates for cell centers (depth)
            latc = nc['latc'][:]
            lonc = nc['lonc'][:]
            # Get Connectivity array
            nv = nc['nv'][:].T - 1 
            # Get depth
            depth = nc['h'][:]  # depth
            dtime = netCDF4.num2date(time_var[itime],time_var.units)
            daystr = dtime.strftime('%Y-%b-%d %H:%M')
            if surf_or_bot==3:# vertically averaged
                m_temp = nc['temp'][itime,0,:]#if len(levels)==0:
                u = nc['ua'][itime,:]# surface fields
                v = nc['va'][itime,:]
            else:    
                m_temp = nc['temp'][itime, surf_or_bot, :]#if len(levels)==0:
                u = nc['u'][itime,surf_or_bot,:]# surface fields
                v = nc['v'][itime,surf_or_bot,:]
            if units=='degF':
                temp=m_temp*1.8+32
            else:
                temp=m_temp
            # transform lon / lat coordinates to map projection
            x,y = m(lons, lats)
            xc,yc=m(lonc,latc)
            
            #kfGB=np.argwhere(lonc>gb[0])&(lonc<gb[1])&(latc>gb[2])&(latc<gb[3]).flatten() # Georges Bank
            #loni=np.arange(gb[0],gb[1],0.02)
            #lati=np.arange(gb[2],gb[3],0.02)
            #lonii,latii=np.meshgrid(loni,lati)
            # grid datam
            #numcols, numrows = 1000, 1000
            numcols, numrows = 25,25
            #xi = np.linspace(x.min(), x.max(), numcols)
            #yi = np.linspace(y.min(), y.max(), numrows)
            xm,ym=m([gb[0],gb[1]],[gb[2],gb[3]])
            xi = np.linspace(xm[0], xm[1], numcols)
            yi = np.linspace(ym[0], ym[1], numrows)
            xi, yi = np.meshgrid(xi, yi)

            # interpolate
            #zi = griddata(x, y, z, xi, yi)
            zi=gd((x,y),temp,(xi,yi),method='linear')
            ui=gd((xc,yc),u,(xi,yi),method='linear')
            vi=gd((xc,yc),v,(xi,yi),method='linear')
            cs=m.contourf(xi,yi,zi,clevs,cmap=plt.get_cmap('rainbow'),zorder=0)
            plt.quiver(xi,yi,ui,vi,scale=20)
            #plt.quiver(xi[::100],yi[::100],ui[::100],vi[::100],scale=10)
            #strm=plt.streamplot(xi,yi,ui,vi,density=1,color=np.sqrt(ui*ui+vi*vi),cmap='GnBu',arrowsize=2,zorder=2)
            zid=gd((x,y),depth,(xi,yi),method='linear')
            #dept_cs=m.contour(xi,yi,zid,dept_clevs,colors='black',zorder=1)
            plt.tricontour(x,y,depth,[200.],colors='purple')
            #plt.clabel(dept_cs, inline = True, fontsize =15,fmt="%1.0f")
    

    # add colorbar.
    cbar = m.colorbar(cs,location='right',pad="2%",size="5%")
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    cbar.set_label(units,fontsize=25)

    #df=pd.read_csv('http://apps-nefsc.fisheries.noaa.gov/drifter/drift_'+cluster+'.csv')
    #df1=pd.read_csv('http://apps-nefsc.fisheries.noaa.gov/drifter/drift_ep_2022_1_ap3.csv')
    #df2=pd.read_csv('http://apps-nefsc.fisheries.noaa.gov/drifter/drift_wms_2022_1.csv')
    df1=pd.read_csv('http://apps-nefsc.fisheries.noaa.gov/drifter/drift_rlsa_2022_1.csv')
    #df2=pd.read_csv('http://apps-nefsc.fisheries.noaa.gov/drifter/drift_shp_2022_1.csv')
    df2=pd.read_csv('drift_shp_2022_1.csv') # had to make a local one with same header
    #df2=pd.read_csv('http://apps-nefsc.fisheries.noaa.gov/drifter/drift_ep_2022_1_ap3.csv')
    df=pd.concat([df1,df2],axis=0)
    ids=np.unique(df['ID'])
    
    for k in ids:
        #print(k)
        datett=[]
        df1=df[df['ID']==k]
        for kk in range(len(df1)):
            datett.append(datetime(int(start_date[0:4]),1,1)+timedelta(df1['YEARDAY'].values[kk]-1))
        df1['datet']=datett
        df1=df1[df1['datet']<dtthis]# only the track less than this time
        df1=df1[df1['datet']>dtthis-timedelta(2.0)]
        #make a plot for this hour and every hour previous with decreasing thickness
        if k==226400691:
            print(len(df1))
        for kk in range(len(df1)-1):
            x1,y1=m(df1['LON'].values[kk],df1['LAT'].values[kk])
            x2,y2=m(df1['LON'].values[kk+1],df1['LAT'].values[kk+1])
            m.plot([x1,x2],[y1,y2],'m',linewidth=kk/(len(df1)/4))#markersize=kk)#,linewidth=kk)
    if include_weather_balloon=='yes':
        # add weather balloon track    
        dfwb=pd.read_csv('ActivityReport.csv',skiprows=4)    
        dfwb[['lat','lon']]=dfwb['Lat/Lng'].str.split(',', 1, expand=True)
        dfwb['datet']=pd.to_datetime(dfwb['Date'])
        dfwb=dfwb[dfwb['datet']<dtthis]# only the track less than this time
        dfwb=dfwb[dfwb['datet']>dtthis-timedelta(2.0)]
        x,y=m(dfwb['lon'].values,dfwb['lat'].values)
        m.plot(x,y,'w',linewidth=4)
    
    # add title
    clayer=''# default current layer
    if (surf_or_bot==-1) & (model[1]=='O'): #case of both DOPPIO and GOMOFS
        layer='surface'
    elif (surf_or_bot==0) & (model[1]=='O'):
        layer='bottom'
    elif (surf_or_bot==-1) & (model=='FVCOM'):
        layer='bottom'
    elif (surf_or_bot==0) & (model=='FVCOM'):
        layer='surface'
    elif surf_or_bot==3:
        clayer='VA_'
        layer='surface'
    #plt.title(model+' '+layer+' temps (color),'+clayer+'current (black) & depth (m)',fontsize=12,fontweight='bold')
    #plt.title('eMOLT bottom temps (white#s) '+model+' '+layer+' temps (color) '+clayer+' & depth (meters)',fontsize=12,fontweight='bold')
    plt.title('eMOLT bottom temps (white#s) '+model+' '+layer+' temps (color) '+clayer+' depth(meters) & drifter (purple)',fontsize=12,fontweight='bold')
    if detide=='y':
        time_str=time_str[0:10]
    plt.suptitle('OOI obs & '+model+' at '+time_str, fontsize=24) 
    if not os.path.exists(path_save):
        os.makedirs(path_save)
    plt.savefig(os.path.join(path_save,time_str.replace(' ','t')+'.png'),dpi=dpi)
    plt.close()

def mean_temp(temps): #makes daily averages if needed for multi-week animations
    mean_temp=temps[0,0]
    for i in range(1,24):
        mean_temp+=temps[i,0]# note we are using the bottom level 0
    return mean_temp/24.0
        
def make_images(dpath,path,dt=datetime(2019,5,1,0,0,0),interval=31,area='OOI',clevs=[39.,44.,0.5],lat_w=42.5,lon_w=-70.3,wind=pd.DataFrame([[datetime.now()-timedelta(1),datetime.now()],[.1,.1],[.1,.1]]),ooi=pd.DataFrame([[datetime.now()-timedelta(1),datetime.now()],[.1,.1],[.1,.1]]),model=model,detide='n',ooiw=pd.DataFrame([[datetime.now()-timedelta(1),datetime.now()],[.1,.1],[.1,.1]])):
    '''dpath: the path of dictionary, use to store telemetered data
        path: use to store images
        dt: start time
        interval: how many days we need make 
        dtimes,u,v are wind in m/s
    '''
    with open(dpath,'rb') as fp:
         telemetered_dict=pickle.load(fp)
    if detide=='n':     
        interval=interval*24
        tdint=timedelta(hours=1)
    else:
        interval=interval
        tdint=timedelta(days=1)
        #tdint=timedelta(hours=1)
    for j in range(interval):
        #dtime=dt+timedelta(days=j)
        dtime=dt+tdint*j
        print(dtime)
        if model=='DOPPIO':
            #url=get_doppio_url(dtime)
            url='http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/History_Best'
            while True:
                if zl.isConnected(address=url):
                    break
                print('check if website is well or internet is connected?')
            time.sleep(5)
            skip=0
            #while True: 
            #        try:
            nc = NetCDFFile(url)
            lons=nc.variables['lon_rho'][:]
            lats=nc.variables['lat_rho'][:]
            temps=nc.variables['temp']
            depth=nc.variables['h'][:]
            time_var = nc['time']
            itime = netCDF4.date2index(dtime,time_var,select='nearest')
            #break
            #        except RuntimeError:
            #            print(str(url)+': need reread')
            #        except OSError:
            #                if zl.isConnected(address=url):
            #                    print(str(url)+': file not exist.')
            #                    skip=1
            #                    break
            #        except KeyboardInterrupt:
            #                    sys.exit()
            m_temp=temps[itime,surf_or_bot]#-1 for surface?
            temp=m_temp*1.8+32

            if skip==1:
                continue
            
        elif model=='GOMOFS':
            url=gomofs_modules.get_gomofs_url(dtime)
            nc = NetCDFFile(url)
            lons=nc.variables['lon_rho'][:]
            lats=nc.variables['lat_rho'][:]
            temps=nc.variables['temp']
            depth=nc.variables['h'][:]
            #time_var = nc['time']
            #itime = netCDF4.date2index(dtime,time_var,select='nearest')
            m_temp=temps[0,surf_or_bot,:,:]#-1 for surface
            temp=m_temp*1.8+32
        elif model=="FVCOM":
            #print('fvcom')
            lons=[]
            lats=[]
            depth=[]
            temp=[]


                
        ntime=dtime
        time_str=ntime.strftime('%Y-%m-%d %H%MUTC')
        
        Year=str(ntime.year)
        Month=str(ntime.month)
        Day=str(ntime.day)
        slons,slats=[],[]
        try:
            slons,slats,stemp=[],[],[]
            for i in telemetered_dict[Year][Month][Day].index:
                slons.append(telemetered_dict[Year][Month][Day]['lon'].iloc[i])
                slats.append(telemetered_dict[Year][Month][Day]['lat'].iloc[i])
                stemp.append(telemetered_dict[Year][Month][Day]['temp'].iloc[i])
        except:
            slons,slats,stemp=[],[],[]
        #print(slons,temp)            
        plotit(lons,lats,slons,slats,stemp,temp,depth,time_str,path,dpi=80,area=area,clevs=clevs,lat_w=lat_w,lon_w=lon_w,wind=wind,ooi=ooi,dtime=dtime,model=model,detide=detide,ooiw=ooiw)
            
        
def read_telemetry(path):
    """read the telemetered data and fix a standard format, the return the standard data"""
    tele_df=pd.read_csv(path,sep='\s+',names=['vessel_n','esn','month','day','Hours','minutes','fracyrday',\
                                          'lon','lat','dum1','dum2','depth','rangedepth','timerange','temp','stdtemp','year'])
    if len(tele_df)<6000:
        print('Warning: the emolt.dat file is not complete at this time.')
        #sys.exit()
        
    return tele_df

#def seperate(filepathsave,ptelemetered='https://www.nefsc.noaa.gov/drifter/emolt.dat'):
#def seperate(filepathsave,ptelemetered='https://apps-nefsc.fisheries.noaa.gov/drifter/emolt_ap3.dat'):
#def make_dict(filepathsave,ptelemetered='https://nefsc.noaa.gov/drifter/emolt.dat'):   
def make_dict(filepathsave,ptelemetered='http://emolt.org/emoltdata/emolt.dat'):   
    '''create a dictionary use to store the data from telemetered, index series is year, month, day and hour
    ptelemetered: the path of telemetered
    '''
    dfdict={}
    df=read_telemetry(ptelemetered)
    
    # JiM added the following 6/16/2021 to get mean value per vessel per day
    df=df[df['depth']>2.0]
    df = df.groupby(['vessel_n','year','month','day'], as_index=False).agg({'temp': 'mean', 'lat': 'mean', 'lon': 'mean','Hours':'first','minutes':'first'})
    for i in df.index:
        #if df['depth'][i]<2.0:
        #    continue
        #if df['minutes'].iloc[i]<=30:
        Ctime=datetime.strptime(str(df['year'].iloc[i])+'-'+str(df['month'].iloc[i])+'-'+str(df['day'].iloc[i])+' '+\
                                         str(df['Hours'].iloc[i])+':'+str(df['minutes'].iloc[i])+':'+'00','%Y-%m-%d %H:%M:%S')
        #else:
        #    Ctime=datetime.strptime(str(df['year'].iloc[i])+'-'+str(df['month'].iloc[i])+'-'+str(df['day'].iloc[i])+' '+\
        #                                 str(df['Hours'].iloc[i])+':'+str(df['minutes'].iloc[i])+':'+'00','%Y-%m-%d %H:%M:%S')+timedelta(seconds=1800)
        Year=str(Ctime.year)
        Month=str(Ctime.month)
        Day=str(Ctime.day)
        #Vessel=df['vessel_n'][i]
        if not Year in dfdict:
            dfdict[Year]={}
        if not Month in dfdict[Year]:
            dfdict[Year][Month]={}
        if not Day in dfdict[Year][Month]:
            dfdict[Year][Month][Day]={}
        #if not Vessel in dfdict[Year][Month][Day]:
        #    dfdict[Year][Month][Day][Vessel]={}

        if len(dfdict[Year][Month][Day])!=0:
            dfdict[Year][Month][Day]=dfdict[Year][Month][Day].append(pd.DataFrame(data=[[df['lat'].iloc[i],df['lon'].iloc[i],df['temp'].iloc[i]]],columns=['lat','lon','temp']).iloc[0])
            dfdict[Year][Month][Day].index=range(len(dfdict[Year][Month][Day]))
        else:
            dfdict[Year][Month][Day]=pd.DataFrame(data=[[df['lat'].iloc[i],df['lon'].iloc[i],df['temp'].iloc[i]]],columns=['lat','lon','temp'])
    with open(filepathsave,'wb') as fp:
        pickle.dump(dfdict,fp,protocol=pickle.HIGHEST_PROTOCOL)


def make_gif(gif_name,png_dir,start_time=False,end_time=False,frame_length = 0.2,end_pause = 4 ):
    '''use images to make the gif
    frame_length: seconds between frames
    end_pause: seconds to stay on last frame
    the format of start_time and end time is string, for example: %Y-%m-%d(YYYY-MM-DD)'''
    
    if not os.path.exists(os.path.dirname(gif_name)):
        os.makedirs(os.path.dirname(gif_name))
    allfile_list = glob.glob(os.path.join(png_dir,'*.png')) # Get all the pngs in the current directory
    print(allfile_list)
    file_list=[]
    if start_time:    
        for file in allfile_list:
            if start_time<=os.path.basename(file).split('.')[0]<=end_time:
                file_list.append(file)
    else:
        file_list=allfile_list
    list.sort(file_list, key=lambda x: x.split('/')[-1].split('t')[0]) # Sort the images by time, this may need to be tweaked for your use case
    images=[]
    # loop through files, join them to image array, and write to GIF called 'wind_turbine_dist.gif'
    for ii in range(0,len(file_list)):       
        file_path = os.path.join(png_dir, file_list[ii])
        if ii==len(file_list)-1:
            for jj in range(0,int(end_pause/frame_length)):
                images.append(imageio.imread(file_path))
        else:
            images.append(imageio.imread(file_path))
    # the duration is the time spent on each image (1/duration is frame rate)
    imageio.mimsave(gif_name, images,'GIF',duration=frame_length)
    
#MAINCODE###########################
start_date_datetime=datetime(int(start_date[0:4]),int(start_date[5:7]),int(start_date[8:10]),0,0,0)
end_date_datetime=datetime(int(start_date[0:4]),int(start_date[5:7]),int(start_date[8:10]),0,0,0)+timedelta(days=ndays)
end_date=str(end_date_datetime.year)+'-'+str(end_date_datetime.month).zfill(2)+'-'+str(end_date_datetime.day).zfill(2)
realpath=os.path.dirname(os.path.abspath(__file__))
dpath=realpath[::-1].replace('py'[::-1],'result/Doppio'[::-1],1)[::-1]  # the directory of the result
if not os.path.exists(dpath):
    os.makedirs(dpath)
dictionary=os.path.join(dpath,'dictionary_emolt.p')
gif_path=os.path.join(dpath,'gif')
map_save=os.path.join(dpath,'map')
gif_name =os.path.join(gif_path,start_date+area+'_'+model+'_'+str(surf_or_bot)+'.gif')
if include_wind=='yes': # get wind time series at one location
    df_w=get_wind_ncep(int(start_date[0:4]),lat_w,lon_w)# returns a dataframe of wind
else:
    df_w=np.nan
if include_ooi=='yes':
    #ooi=pd.read_csv('http://erddap.dataexplorer.oceanobservatories.org/erddap/tabledap/ooi-cp02pmuo-wfp01-01-vel3dk000.csvp?time%2Ceastward_sea_water_velocity_profiler_depth_enabled%2Cnorthward_sea_water_velocity_profiler_depth_enabled&time%3E=2021-10-04T00%3A00%3A00Z&time%3C=2021-10-06T15%3A04%3A00Z&z%3E=-430&z%3C=-17')
    #url='http://erddap.dataexplorer.oceanobservatories.org/erddap/tabledap/ooi-cp03issm-rid26-04-velpta000.csvp?time%2Clatitude%2Clongitude%2Ceastward_sea_water_velocity%2Cnorthward_sea_water_velocity&time%3E=2022-06-18T02%3A00%3A00Z'
    url='http://erddap.dataexplorer.oceanobservatories.org/erddap/tabledap/ooi-cp01cnsm-rid26-04-velpta000.csvp?time%2Clatitude%2Clongitude%2Ceastward_sea_water_velocity%2Cnorthward_sea_water_velocity&time%3E=2022-06-20T02%3A00%3A00Z'
    ooi=pd.read_csv(url)
    ooi=ooi.set_index('time (UTC)')
    ooi.index=pd.to_datetime(ooi.index)
    ooi=ooi[~ooi.index.duplicated()]
    #ooi=ooi.rename(columns={'eastward_sea_water_velocity_profiler_depth_enabled (m.s-1)':'u','northward_sea_water_velocity_profiler_depth_enabled (m.s-1)':'v'})
    ooi=ooi.rename(columns={'eastward_sea_water_velocity (m.s-1)':'u','northward_sea_water_velocity (m.s-1)':'v'})
    # now the wind
    urlw='http://erddap.dataexplorer.oceanobservatories.org/erddap/tabledap/ooi-cp01cnsm-sbd11-06-metbka000.csvp?time%2Ceastward_sea_water_velocity%2Cnorthward_sea_water_velocity%2Ceastward_wind%2Cnorthward_wind&time%3E=2021-10-25T00%3A00%3A00Z'
    ooiw=pd.read_csv(urlw)
    ooiw=ooiw.set_index('time (UTC)')
    ooiw.index=pd.to_datetime(ooiw.index)
    ooiw=ooiw[~ooiw.index.duplicated()]
    ooiw=ooiw.rename(columns={'eastward_wind (m.s-1)':'uw','northward_wind (m.s-1)':'vw'})

else:
    ooi=np.nan
    ooiw=np.nan
    

#############################
#run functions
if include_temp_obs=='yes':
    make_dict(filepathsave=dictionary)
make_images(dpath=dictionary,path=map_save,dt=start_date_datetime,interval=ndays,area=area,clevs=clevs,lat_w=lat_w,lon_w=lon_w,wind=df_w,ooi=ooi,model=model,detide=detide,ooiw=ooiw)
make_gif(gif_name,map_save,start_time=start_date,end_time=end_date)# sending datetimes
    


