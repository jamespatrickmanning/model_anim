# routine to plot & save model grids
# JiM in May 2018
# modified in Aug 2019 to choose particular grids
import numpy as np
from matplotlib import pyplot as plt
from netCDF4 import Dataset        # netCDF4 version


grid='GOM3'
get_input_local='yes'
save_results='no'
plot_results='no'

if get_input_local=='yes':
	lat=np.load(grid+'.lat.npy')
        lon=np.load(grid+'.lon.npy') 
else: 
	if grid=='GOM3':
		URL='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'
	elif grid=='MASSBAY':
		URL='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/NECOFS_MASS_BAY/2018/mbn_201801.nc'
	print 'getting FVCOM  Grid parameters from web'
	print URL
	ds = Dataset(URL,'r').variables   # netCDF4 version
	h=ds['h']; 
	print 'getting lon,lat'
	lat=ds['lat']; 
	lon=ds['lon']; 
	#latc=ds['latc'];
	#lonc=ds['lonc']; 

if plot_results=='yes':
	import sys
	sys.path.append("../modules")
	import basemap as bm
	plt.plot(lon,lat,'r.',markersize=1)
        #bm.basemap_detail([40.0,43.0],[-71.,-65.0],True,False)# no bathy and no lat/lon lines
        bm.basemap_region('wv')
        plt.title(grid+' Grid')
        plt.savefig(grid+'_grid.png')

if save_results=='yes':
	np.save(grid+'.h.npy',np.array(h))
	np.save(grid+'.lat.npy',np.array(lat))
	np.save(grid+'.lon.npy',np.array(lat))
        

