#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 12:10:25 2019
GoMOFs
@author: leizhao
Modified by JiM in Sep 2019 with more comments
Modified by JiM in Oct 2019 with gbox
"""

import os,imageio
import conda
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib
from mpl_toolkits.basemap import Basemap
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
import time
import zlconversions as zl
import sys
import pandas as pd
try:
    import cPickle as pickle
except ImportError:
    import pickle
import glob
import math

def getgbox(area):
  # gets geographic box based on area
  if area=='SNE':
    gbox=[-70.,-64.,39.,42.] # for SNE
  elif area=='OOI':
    gbox=[-71.,-69.,39.,41.] # for OOI
  elif area=='GBANK':
    gbox=[-70.,-64.,39.,42.] # for GBANK
  elif area=='GS':           
    gbox=[-71.,-63.,38.,42.5] # for Gulf Stream
  elif area=='NorthShore':
    gbox=[-71.,-70.,42.,43.] # for north shore
  elif area=='CCBAY':
    gbox=[-70.75,-69.8,41.5,42.23] # CCBAY
  elif area=='inside_CCBAY':
    gbox=[-70.75,-70.,41.7,42.23] # CCBAY
  elif area=='NEC':
    gbox=[-69.,-64.,39.,43.5] # NE Channel
  return gbox

def get_gomofs_url(date):
    """
    the format of date is:datetime.datetime(2019, 2, 27, 11, 56, 51, 666857)
    the date is the american time 
    input date and return the url of data
    """
#    print('start calculate the url!') 
    date=date+timedelta(hours=4.5)
    date_str=date.strftime('%Y%m%d%H%M%S')
    hours=int(date_str[8:10])+int(date_str[10:12])/60.+int(date_str[12:14])/3600.
    tn=int(math.floor((hours)/6.0)*6)  ## for examole: t12z the number is 12
    if len(str(tn))==1:
        tstr='t0'+str(tn)+'z'   # tstr in url represent hour string :t00z
    else:
        tstr='t'+str(tn)+'z'
    if round((hours)/3.0-1.5,0)==tn/3:
        nstr='n006'       # nstr in url represent nowcast string: n003 or n006
    else:
        nstr='n003'
    url='http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/GOMOFS/MODELS/'\
    +date_str[:6]+'/nos.gomofs.fields.'+nstr+'.'+date_str[:8]+'.'+tstr+'.nc'
    return url

def plotit(lons,lats,slons,slats,temp,depth,time_str,path_save,dpi=80,area='OOI'):
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_axes([0.01,0.05,0.98,0.87])
    # create polar stereographic Basemap instance.
    gb=getgbox(area)
    m = Basemap(projection='stere',lon_0=(gb[0]+gb[1])/2.,lat_0=(gb[2]+gb[3])/2.,lat_ts=0,llcrnrlat=gb[2],urcrnrlat=gb[3],\
                llcrnrlon=gb[0],urcrnrlon=gb[1],rsphere=6371200.,resolution='f',area_thresh=100)
#    m = Basemap(projection='stere',lon_0=-67,lat_0=42,lat_ts=0,llcrnrlat=38,urcrnrlat=46,\
#                llcrnrlon=-73,urcrnrlon=-61,rsphere=6371200.,resolution='i',area_thresh=100)
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    if len(slons)!=0:
        x1,y1=m(slons,slats)
        ax.plot(x1,y1,'ro',markersize=10)
    # draw parallels.
    parallels = np.arange(0.,90,1.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=20)
    # draw meridians
    meridians = np.arange(180.,360.,1.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=20)
    x, y = m(lons, lats) # compute map proj coordinates.
    # draw filled contours.
    
    
    dept_clevs=[50,100, 150,300,1000]
    dept_cs=m.contour(x,y,depth,dept_clevs,colors='black')
    plt.clabel(dept_cs, inline = True, fontsize =15,fmt="%1.0f")
    
    
    clevs=np.arange(39,47,0.5)  #for all year:np.arange(34,84,1) or np.arange(34,68,1)
    cs = m.contourf(x,y,temp,clevs,cmap=plt.get_cmap('rainbow'))
    #plt.xlim([gbox[0],gbox[1]])
    #plt.ylim([gbox[2],gbox[3]])
    # add colorbar.
    cbar = m.colorbar(cs,location='right',pad="2%",size="5%")
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    cbar.set_label('Fahrenheit',fontsize=25)
    # add title
    plt.title('GoMOFS MODEL BOTTOM TEMP '+time_str,fontsize=30)
    if not os.path.exists(path_save):
        os.makedirs(path_save)
    #print(time_str)
    #print(os.path.join(path_save,time_str.replace(' ','t')+'.png'))
    #plt.show()
    plt.savefig(os.path.join(path_save,time_str.replace(' ','t')+'.png'),dpi=dpi)

def make_images(dpath,path,dt=datetime(2019,5,1,0,0,0),interval=31,area='OOI'):
#def make_images(dpath,path,dt=datetime(2019,5,1,0,0,0),area='OOI'):
    '''dpath: the path of dictionary, use to store telemetered data
        path: use to store images
        dt: start time
        interval: how many days we need make 
        gbox: area to focus on
    '''
    with open(dpath,'rb') as fp: #load the data of observation
         telemetered_dict=pickle.load(fp) 
    for j in range(interval): # loop every days files 
        dtime=dt+timedelta(days=j)
        print(dtime)
        count,skip=0,0  #count use to count how many files load successfully
        for i in range(0,24,3): #loop every file of day, every day have 8 files
            ntime=dtime+timedelta(hours=i)
            url=get_gomofs_url(ntime)
            print(url)
            while True:#check the internet
                if zl.isConnected(address=url):
                    break
                print('check the website is well or internet is connected?')
                time.sleep(5)
            while True:  #load data
                try:
                    nc = NetCDFFile(url)
                    lons=nc.variables['lon_rho'][:]
                    lats=nc.variables['lat_rho'][:]
                    temps=nc.variables['temp']
                    depth=nc.variables['h'][:]
                    break
                except KeyboardInterrupt:
                    sys.exit()
                except OSError:
                    if zl.isConnected(address=url):
                        print(str(url)+': file not exit.')
                        skip=1
                        break
                except:
                    print('reread data:'+str(url))
            if skip==1:  #if file is not exist   
                continue
            m_temp=temps[0,0] # JiM added this 2/19/2020
            '''
            if i==0: 
                count+=1
                m_temp=temps[0,0]
            else:
                m_temp+=temps[0,0]
                count+=1
            '''
            #m_temp=m_temp/float(count)
            #ntime=dtime
            time_str=ntime.strftime('%Y-%m-%d-%H')
            temp=m_temp*1.8+32
            Year=str(ntime.year)
            Month=str(ntime.month)
            Day=str(ntime.day)
            slons,slats=[],[]
            try:
                slons,slats=[],[]
                for i in telemetered_dict[Year][Month][Day].index:
                    slons.append(telemetered_dict[Year][Month][Day]['lon'].iloc[i])
                    slats.append(telemetered_dict[Year][Month][Day]['lat'].iloc[i])
            except:
                slons,slats=[],[]
            dpi=80
            plotit(lons,lats,slons,slats,temp,depth,time_str,path,dpi,area)
                
def read_telemetry(path):
    """read the telemetered data and fix a standard format, the return the standard data"""
    tele_df=pd.read_csv(path,sep='\s+',names=['vessel_n','esn','month','day','Hours','minates','fracyrday',\
                                          'lon','lat','dum1','dum2','depth','rangedepth','timerange','temp','stdtemp','year'])
    if len(tele_df)<6000:
        #print('check the telemetered website!')
        print('Warning: the emolt.dat file is not complete at this time.')
        #sys.exit()
        
    return tele_df


def seperate(filepathsave,ptelemetered='https://www.nefsc.noaa.gov/drifter/emolt.dat'):
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

#MAIN CODE
#def main():  
#Hardcodes#######################
area='NorthShore'
start_date='2020-02-11'
ndays=3
##################################
start_date_datetime=datetime(int(start_date[0:4]),int(start_date[5:7]),int(start_date[8:10]),0,0,0)
end_date_datetime=datetime(int(start_date[0:4]),int(start_date[5:7]),int(start_date[8:10]),0,0,0)+timedelta(days=ndays)
end_date=str(end_date_datetime.year)+'-'+str(end_date_datetime.month).zfill(2)+'-'+str(end_date_datetime.day).zfill(2)
realpath=os.path.dirname(os.path.abspath(__file__))
dpath=realpath[::-1].replace('py'[::-1],'result/GOMOFS/'[::-1],1)[::-1]  # the directory of the result
if not os.path.exists(dpath):
    os.makedirs(dpath)
dictionary=os.path.join(dpath,'dictionary_emolt.p')
gif_path=os.path.join(dpath,'gif')
map_save=os.path.join(dpath,'map')
gif_name =os.path.join(gif_path,start_date+area+'_GOMOFS.gif')
    
#############################
#run functions
seperate(filepathsave=dictionary)
gbox=getgbox(area) # gets lon and lat limits
make_images(dpath=dictionary,path=map_save,dt=start_date_datetime,interval=ndays,area=area)
make_gif(gif_name,map_save,start_time=start_date,end_time=end_date)

#if __name__ == "__main__":
#    main()
