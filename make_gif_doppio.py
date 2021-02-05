#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Apr 17 11:29:14 2019
@author: leizhao
Modified by JiM in Fall 2019
Modified by JiM in Summer 2020 with Jack Polentes help
Modified by JiM in Fall 2020 by adding drifter tracks and, in late Oct 2020, NCEP wind vector
See hardcodes at top of code where there is, for example, "start_date" and "ndays" 
Note: You may need to adjust the "clevs" to get good range of colors.
Note: You may want to adjust the Basemap resolution to "c" for crude in experiments.
Note: You may need to "conda install -c conda-forge basemap-data-hires" in order to use the higher resolution coastlines
"""
#hardcodes########################
area='NE'
#area='IpswichBay'
start_date='2020-12-25'
#clevs=[55.,80.,2.0]# min,max, and interval of temperature contour levels wanted
clevs=[44.,64.,2.0]
ndays=2 # number of days to run
include_temp_obs='no' # 'yes' overlays LFA positions of observed bottom temps
include_wind='no'
cluster='ga_2020_1'
#cluster='X'
surf_or_bot=0 # -1 for surface
lat_w,lon_w=40.,-68.5 # base of wind vector legend (actual vector appears mid plot)
#########
import os,imageio
import pytz
import conda
import pandas as pd
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
#os.environ["PROJ_LIB"] = proj_lib NOTE: JiM had to point to specific directory
os.environ['PROJ_LIB'] = 'c:\\Users\\Joann\\anaconda3\\pkgs\\proj4-5.2.0-ha925a31_1\\Library\share'
from mpl_toolkits.basemap import Basemap
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
from matplotlib.dates import date2num,num2date
import time
import zlconversions as zl
import sys
import warnings
warnings.filterwarnings("ignore") # gets rid of warnings at runtime but you may want to comment this out to see things
import pandas as pd
import netCDF4 # for wind access
from operator import itemgetter
from math import sqrt
try:
    import cPickle as pickle
except ImportError:
    import pickle
import glob

def get_doppio_url(dtime):
    '''dtime ids gmt time'''
    date=dtime.strftime('%Y-%m-%d')
    url='http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/runs/History_RUN_2018-11-12T00:00:00Z'
    #url='http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/History_Best'
    return url.replace('2018-11-12',date)
    #return url

def getgbox(area):
  # gets geographic box based on area
  if area=='SNE':
    gbox=[-70.,-64.,39.,42.] # for SNE
  elif area=='OOI':
    gbox=[-72.,-69.5,39.5,41.5] # for OOI
  elif area=='GBANK':
    gbox=[-71.,-66.,40.,42.] # for GBANK
  elif area=='GBANK_RING':
    #gbox=[-70.,-64.,39.,42.] # for typical GBANK Ring 
    gbox=[-70.,-64.,38.,42.] # for typical GBANK Ring 
  elif area=='GS':           
    gbox=[-71.,-63.,38.,42.5] # for Gulf Stream
  elif area=='NorthShore':
    gbox=[-71.,-69.5,41.5,43.] # for north shore
  elif area=='IpswichBay':
    gbox=[-71.,-70.,42.5,43.] # for north shore
  elif area=='CCBAY':
    gbox=[-70.75,-69.8,41.5,42.23] # CCBAY
  elif area=='inside_CCBAY':
    gbox=[-70.75,-70.,41.7,42.15] # CCBAY
  elif area=='NEC':
    gbox=[-69.,-64.,39.,43.5] # NE Channel
  elif area=='NE':
    gbox=[-76.,-66.,35.,44.5] # NE Shelf 
  return gbox


#def get_wind_ncep(starttime,endtime,lat,lon):
def get_wind_ncep(year,lat,lon):    
        #function get a time series of u & v wind m/s from pre-downloaded ncep file
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

def plotit(lons,lats,slons,slats,stemp,temp,depth,time_str,path_save,dpi=80,area='OOI',clevs=[39.,44.,0.5],lat_w=42.5,lon_w=-70.5,wind=pd.DataFrame([[datetime.now()-timedelta(1),datetime.now()],[.1,.1],[.1,.1]])):
    '''dtimes,u,v are wind in m/s'''
    fig = plt.figure(figsize=(10,6))
    #ax = fig.add_axes([0.01,0.05,0.98,0.87])
    ax = fig.add_axes([0.05,0.05,0.9,0.82])
    # create polar stereographic Basemap instance.
    gb=getgbox(area)
    m = Basemap(projection='stere',lon_0=(gb[0]+gb[1])/2.,lat_0=(gb[2]+gb[3])/2.,lat_ts=0,llcrnrlat=gb[2],urcrnrlat=gb[3],\
                llcrnrlon=gb[0],urcrnrlon=gb[1],rsphere=6371200.,resolution='f',area_thresh=100)# JiM changed resolution to "c" for crude
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    #m.drawstates()
    #m.drawcountries()
    if include_wind=='yes':
       xw,yw=m(lon_w,lat_w)  
       ax.quiver(xw,yw,10.0,0.0,color='red',scale=50.,zorder=2)#wind legend in middle
       ax.text(xw,yw+(m.ymin+m.ymax)/20,'10 m/s (~20 knots) NCEP wind',fontsize=16,color='red',zorder=2) 
       idex=wind.index.get_loc(datetime.strptime(time_str[0:15],"%Y-%m-%d %H%M"),method='nearest')
       #idex=min(range(len(dtimes)), key=lambda x: abs(x-datetime.strptime(time_str[0:15],"%Y-%m-%d %H%M").toordinal()))
       #print(str(wind.u[idex]),str(wind.v[idex]))
       ax.quiver((m.xmin+m.xmax)/2,(m.ymin+m.ymax)/2,wind.u[idex],wind.v[idex],scale=50.0,color='red',zorder=2) # plot an arrow
    if include_temp_obs=='yes':
        dfLFoM=pd.read_csv('C:\\Users\Joann\Downloads\LFoM\LFoM_sites.csv')
        xx,yy=m(list(dfLFoM['lon']),list(dfLFoM['lat']))
        ax.plot(xx,yy,'ro')
    if len(slons)!=0:
        x1,y1=m(slons,slats)
        for jj in range(len(stemp)):
            ax.text(x1[jj],y1[jj],str.format('{0:.1f}',stemp[jj]*1.8+32),color='w',fontsize=10,fontweight='bold',horizontalalignment='center',verticalalignment='center')
    # draw parallels.
    if area=='IpswichBay':
        labint=0.2
        dept_clevs=[30,50,100, 150]
    elif area=='NE':
        labint=2.0
        dept_clevs=[50,100,1000]
    else:
        labint=1.0
        dept_clevs=[30,50,100, 150,300,1000]
    parallels = np.arange(0.,90,labint)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=20)
    # draw meridians
    meridians = np.arange(180.,360.,labint)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=20)
    x, y = m(lons, lats) # compute map proj coordinates.
    # draw filled contours.
    
    
    
    dept_cs=m.contour(x,y,depth,dept_clevs,colors='black')
    plt.clabel(dept_cs, inline = True, fontsize =15,fmt="%1.0f")
    
    
    clevs=np.arange(clevs[0],clevs[1],clevs[2])  #for all year:np.arange(34,84,1) or np.arange(34,68,1)
    cs = m.contourf(x,y,temp,clevs,cmap=plt.get_cmap('rainbow'))
    # add colorbar.
    cbar = m.colorbar(cs,location='right',pad="2%",size="5%")
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    cbar.set_label('Fahrenheit',fontsize=25)

    #df=pd.read_csv('http://nefsc.noaa.gov/drifter/drift_'+cluster+'.csv')
    df=pd.read_csv('http://apps-nefsc.fisheries.noaa.gov/drifter/drift_'+cluster+'.csv')
    ids=np.unique(df['ID'])
    #datetd=[]#drifters datetime
    dtthis=datetime(2020,int(time_str[5:7]),int(time_str[8:10]),int(time_str[11:13]),int(time_str[13:15]))
    for k in ids:
        #print(k)
        datett=[]
        df1=df[df['ID']==k]
        for kk in range(len(df1)):
            #datett.append(datetime(2020,df1['MTH'][kk],df1['DAY'][kk],df1['HR_GMT'][kk],df1['MIN'][kk],0,pytz.UTC))
            #datetd.append(num2date(df1['YEARDAY'][kk]).replace(year=2020).replace(tzinfo=pytz.UTC))
            datett.append(datetime(2020,1,1)+timedelta(df1['YEARDAY'].values[kk]))
             #.replace(month=int(df1['MTH'][kk])).replace(day=int(df1['DAY'][kk]))
        df1['datet']=datett
        #df1=df1[df1['DAY']==datetime_wanted.day]
        #df1=df1[df1['MTH']==datetime_wanted.month]# JiM added sep 2020
        #df1=df1[df1['DAY']==int(time_str[8:10])]
        #df1=df1[df1['MTH']==int(time_str[5:7])]# JiM added sep 2020 
        #df1=df1[df1['HR_GMT']==int(time_str[11:13])]
        df1=df1[df1['datet']<dtthis]# only the track less than this time
        df1=df1[df1['datet']>dtthis-timedelta(1.0)]
        x,y=m(df1['LON'].values,df1['LAT'].values)
        m.plot(x,y,'m',linewidth=4)
    # add title
    if surf_or_bot==-1:
        layer='surface'
    else:
        layer='bottom'
    plt.title('Doppio '+layer+' temps (color), observations (white#s) & depth (m)',fontsize=12,fontweight='bold')
    plt.suptitle(time_str, fontsize=24)
    if not os.path.exists(path_save):
        os.makedirs(path_save)
    plt.savefig(os.path.join(path_save,time_str.replace(' ','t')+'.png'),dpi=dpi)
    plt.close()

def mean_temp(temps): #makes daily averages if needed for multi-week animations
    mean_temp=temps[0,0]
    for i in range(1,24):
        mean_temp+=temps[i,0]# note we are using the bottom level 0
    return mean_temp/24.0
        
def make_images(dpath,path,dt=datetime(2019,5,1,0,0,0),interval=31,area='OOI',clevs=[39.,44.,0.5],lat_w=42.5,lon_w=-70.3,wind=pd.DataFrame([[datetime.now()-timedelta(1),datetime.now()],[.1,.1],[.1,.1]])):
    '''dpath: the path of dictionary, use to store telemetered data
        path: use to store images
        dt: start time
        interval: how many days we need make 
        dtimes,u,v are wind in m/s
    '''
    with open(dpath,'rb') as fp:
         telemetered_dict=pickle.load(fp)
    interval=interval*24
    for j in range(interval):
        #dtime=dt+timedelta(days=j)
        dtime=dt+timedelta(hours=j)
        print(dtime)
        #url=get_doppio_url(dtime)
        url='http://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/History_Best'
        while True:
            if zl.isConnected(address=url):
                break
            print('check the website is well or internet is connected?')
            time.sleep(5)
        skip=0
        while True: 
            try:
                nc = NetCDFFile(url)
                lons=nc.variables['lon_rho'][:]
                lats=nc.variables['lat_rho'][:]
                temps=nc.variables['temp']
                depth=nc.variables['h'][:]
                time_var = nc['time']
                itime = netCDF4.date2index(dtime,time_var,select='nearest')
                break
            except RuntimeError:
                print(str(url)+': need reread')
            except OSError:
                if zl.isConnected(address=url):
                    print(str(url)+': file not exit.')
                    skip=1
                    break
            except KeyboardInterrupt:
                sys.exit()
        if skip==1:
            continue
                
        #m_temp=mean_temp(temps)# here we are taking a daily average
        #m_temp=temps[j,0] # JiM made this change 7/16/2020 since we are looking at hourly not daily images
        m_temp=temps[itime,surf_or_bot]#-1 for surface?
        #m_temp=temps[np.mod(j,24),0]
        ntime=dtime
        #time_str=ntime.strftime('%Y-%m-%d')
        time_str=ntime.strftime('%Y-%m-%d %H00UTC')
        temp=m_temp*1.8+32
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
        plotit(lons,lats,slons,slats,stemp,temp,depth,time_str,path,dpi=80,area=area,clevs=clevs,lat_w=lat_w,lon_w=lon_w,wind=wind)
        
def read_telemetry(path):
    """read the telemetered data and fix a standard format, the return the standard data"""
    tele_df=pd.read_csv(path,sep='\s+',names=['vessel_n','esn','month','day','Hours','minates','fracyrday',\
                                          'lon','lat','dum1','dum2','depth','rangedepth','timerange','temp','stdtemp','year'])
    if len(tele_df)<6000:
        print('Warning: the emolt.dat file is not complete at this time.')
        #sys.exit()
        
    return tele_df

#def seperate(filepathsave,ptelemetered='https://www.nefsc.noaa.gov/drifter/emolt.dat'):
#def seperate(filepathsave,ptelemetered='https://apps-nefsc.fisheries.noaa.gov/drifter/emolt_ap3.dat'):
def seperate(filepathsave,ptelemetered='https://nefsc.noaa.gov/drifter/emolt.dat'):   
    '''create a dictionary use to store the data from telemetered, index series is year, month, day and hour
    ptelemetered: the path of telemetered
    '''
    dfdict={}
    df=read_telemetry(ptelemetered)
    for i in df.index:
        if df['depth'][i]<2.0:
            continue
        if df['minates'].iloc[i]<=30:
            Ctime=datetime.strptime(str(df['year'].iloc[i])+'-'+str(df['month'].iloc[i])+'-'+str(df['day'].iloc[i])+' '+\
                                         str(df['Hours'].iloc[i])+':'+str(df['minates'].iloc[i])+':'+'00','%Y-%m-%d %H:%M:%S')
        else:
            Ctime=datetime.strptime(str(df['year'].iloc[i])+'-'+str(df['month'].iloc[i])+'-'+str(df['day'].iloc[i])+' '+\
                                         str(df['Hours'].iloc[i])+':'+str(df['minates'].iloc[i])+':'+'00','%Y-%m-%d %H:%M:%S')+timedelta(seconds=1800)
        Year=str(Ctime.year)
        Month=str(Ctime.month)
        Day=str(Ctime.day)
        if not Year in dfdict:
            dfdict[Year]={}
        if not Month in dfdict[Year]:
            dfdict[Year][Month]={}
        if not Day in dfdict[Year][Month]:
            dfdict[Year][Month][Day]={}

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
gif_name =os.path.join(gif_path,start_date+area+'_Doppio.gif')
if include_wind=='yes': # get wind time series at one location
    df_w=get_wind_ncep(int(start_date[0:4]),lat_w,lon_w)# returns a dataframe of wind
else:
    df_w=np.nan
    

#############################
#run functions
if include_temp_obs=='yes':
    seperate(filepathsave=dictionary)
make_images(dpath=dictionary,path=map_save,dt=start_date_datetime,interval=ndays,area=area,clevs=clevs,lat_w=lat_w,lon_w=lon_w,wind=df_w)
make_gif(gif_name,map_save,start_time=start_date,end_time=end_date)# sending datetimes
    


