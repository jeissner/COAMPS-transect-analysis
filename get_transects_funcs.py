#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:56:36 2024

@author: jordaneissner
"""
import act
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import glob
from netCDF4 import Dataset
from matplotlib import ticker
#import seaborn as sns
from HaversineFunc import Haversine
from matplotlib import cbook
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm
import matplotlib as mpl
from stability_func import stability
from CAO_decouple_func import CAO_decoup_params
from matplotlib import ticker
import Ngl
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
import pandas as pd
import xarray

import sys
sys.path.insert(0, '/Users/jordaneissner/Documents/clone/')
from transect_test import pointfind2
from transect_test import get_transect

def cloudFrac(clf,cheight, lower, upper):
    
    x = np.argmin(np.abs(cheight-lower))
    y = np.argmin(np.abs(cheight-upper))    
    tclf = np.zeros([np.shape(clf)[1],np.shape(clf)[1]])
    tclf[:,:] = np.nan
    
    for i in range(0,np.shape(clf)[1]):
        for j in range(0,np.shape(clf)[2]):
            if np.any(clf[:,i,j]) == True: 
                tclf[i,j] = np.max(clf[y:x,i,j])
    return(tclf)

def cloudBounds(qc,height):
    p = qc[0,:,:]
    cbh = p.copy()
    cbh[:] = np.nan
    cth = p.copy()
    cth[:] = np.nan

    for i in range(0,np.shape(p)[0]):
        for j in range(0,np.shape(p)[0]):
                q = np.where(qc[:,i,j] > 0.00001) #.01 g/kg (Covert et al. 2022)
                if np.shape(q[0])[0] > 2 :
                    cbh[i,j] = height[q[0][-1]]
                    cth[i,j] = height[q[0][0]]
                 
    thick = (cth - cbh)/1000. 
    
    return(cbh,cth,thick)

def get_sector(xkm):
    
    sector = xkm.copy()
    sector[:] = np.nan
    
    xms = np.where(xkm <= -100) 
    x= np.where ((xkm > -100) & (xkm < 100))
    xp = np.where((xkm >= 100)) 

    sector[xms] = 0 # cold
    sector[x] = 1 # front
    sector[xp] = 2 # warm
    
    return(sector)

def map_var_plot(lon,lat,var,trans_lon,trans_lat,r,hr,v):

    fig,(ax1) = plt.subplots(1,1, figsize=(10,10), subplot_kw={'projection': ccrs.LambertConformal(central_latitude = 39, 
                            central_longitude = 332, 
                            standard_parallels = (39, 39)) })
      
    extent = [-35,-21.5,33,44]#[-46, -10.5, 25,57] #[-54,0, 10,70]#[[-35,-21.5,33,44]
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False, linewidth=0.33, color='k',alpha=0.5)
    gl.right_labels = gl.top_labels = False
    gl.ylocator = ticker.FixedLocator([25,30,35,40, 45, 50, 55,60])
    gl.xlocator = ticker.FixedLocator([-40, -30, -20, -10])
    gl.xlabel_style = {'size': 16, 'color': 'k', 'rotation':45, 'ha':'right'}
    gl.ylabel_style = {'size': 16, 'color': 'k', 'weight': 'normal'}

    ax1.set_extent(extent)
    ax1.coastlines(resolution='50m', linewidth=0.4)
    cmap0 = LinearSegmentedColormap.from_list('', ['white', *plt.cm.Blues(np.arange(255))])

    CS = ax1.pcolormesh(lon,lat,var,transform=ccrs.PlateCarree())#,cmap=cmap0,norm=LogNorm()) #vmin=0,vmax=300)

    ax1.scatter(trans_lon,trans_lat,transform=ccrs.PlateCarree(),color='k', marker='+',s=1)
    ax1.scatter([trans_lon[0], trans_lon[r], trans_lon[-1]],[trans_lat[0],trans_lat[r],trans_lat[-1]], transform=ccrs.PlateCarree(), c='r', s=45)
   
    cbar = plt.colorbar(CS, ax=ax1, orientation='vertical', cmap=cmap0, pad=0.04,extend='max')#, fraction=0.046, pad=0.04)
   # cbar.set_label('LWP (g m$^{-2}$)',size=20)
    cbar.ax.tick_params(labelsize=18) 
    
    ax1.set_title(hr +' HR ' + str(v))
    plt.show()


def line1(x):
    y = (0.7894 * x) - 212.737 #215
    return(y)
def line2(x):
    y = (0.7894 * x) - 231.4278
    return(y)


def bin_avg(x, var, nbins):
    bins = np.arange(-1450,1350,3000/(nbins-1))
    avgvar = np.zeros(len(bins)-1)
    stdvar = np.zeros(len(bins)-1)
    for b in range(0,len(bins)-1):
        print(b,bins[b])
        c,d = np.where( (x[:,:] > bins[b]) & (x[:,:] < bins[b+1]) )
        avgvar[b] = np.nanmean(var[c,d])
        stdvar[b] = np.std(var[c,d])
    return([bins,avgvar,stdvar])

def bin_med(x, var, nbins):
    bins = np.arange(-1450,1350,3000/(nbins-1))
    avgvar = np.zeros(len(bins)-1)
    stdvar = np.zeros(len(bins)-1)
    for b in range(0,len(bins)-1):
        c,d = np.where( (x[:,:] > bins[b]) & (x[:,:] < bins[b+1]) )
        avgvar[b] = np.nanmedian(var[c,d])
        stdvar[b] = np.std(var[c,d])
    return([bins,avgvar,stdvar])

def to_netcdf(DATA,filename,xrname,xrunit):
    xr = xarray.Dataset.from_dataframe(DATA)
    xr['xkm'].attrs={'units':'km', 'long_name':'km_from_coldfront'}
    xr['sector'].attrs={'units': '0=cold,1=front,2=warm'}
    for i in range(0,len(xrname)):
        xr[xrname[i]].attrs={'units': xrunit[i]}

    xr.to_netcdf(filename)
    
def get_transects(file_path,xs1,xs2, var,var2, lat, lon, r,hours, ix,nt, nbins, filename):
        
    
        ys1 = line1(xs1)
        ys2 = line2(xs2)  
        
        pn=nbins
        
        if ix > 0:
            XKM = []
            LWP = []
            CTH = []
            SEC = []
            for v in range(0,len(xs1)):
                data_prof, x_kilometers, m_grid, n_grid = get_transect(lat, lon, var2, ys1[v], xs1[v], ys2[v], xs2[v], npoints = pn, pdif = 1, norep=False)
                mgrid = [int(x) for x in m_grid]
                ngrid = [int(x) for x in n_grid]
                
                for h in range(0,len(hours)):
                    hour = hours[h]
                    print('hr', v, hour)
                    file1 = file_path+hour+'.cdf'
                    ds = act.io.read_arm_netcdf(file1)
                    lat = ds.latitude.values
                    lon = ds.longitude.values
                    lwp = ds.cloud_lwp_iwp.values
                    cth = ds.cloud_top_height.values
                    cp = ds.cloud_phase.values
                    
                    
                    x_kilometers2 = x_kilometers - x_kilometers[r[h]]

                    #filter out ice clouds
                    cpx, cpy = np.where((cp == 2) | (cp == 7)) #1 = water
                    lwp[cpx,cpy] = np.nan
                    cth[cpx,cpy] = np.nan
                    
                    secs = get_sector(x_kilometers2)
                    SEC.extend(secs)
                    LWP.extend(lwp[mgrid,ngrid])
                    CTH.extend(cth[mgrid,ngrid])
                    XKM.extend(x_kilometers2)

            d = {'xkm': XKM, 'LWP':LWP, 'CTH':CTH, 'sector': SEC}
            DATA = pd.DataFrame(data=d)
            to_netcdf(DATA, filename)
            
            return(DATA)
     
        else: 
            trans_lon = np.zeros(pn)
            trans_lat = np.zeros(pn)
     
            trans_var = np.zeros(pn)
            print('here', pn, np.shape(var2), np.shape(lon))
            data_prof, x_kilometers, m_grid, n_grid = get_transect(lat, lon, var2, ys1[18], xs1[18], ys2[18], xs2[18], npoints = pn, pdif=1, norep=False )
            mgrid = [int(x) for x in m_grid]
            ngrid = [int(x) for x in n_grid]

            for i in range(0,len(mgrid)):
                trans_lon[i] = np.median(lon[mgrid[i], ngrid[i]])
                trans_lat[i] = np.median(lat[mgrid[i], ngrid[i]])

                trans_var[i] = var[mgrid[i], ngrid[i]]

            x_kilometers = x_kilometers - x_kilometers[r] #115
            map_var_plot(lon,lat,var,lon[mgrid,ngrid],lat[mgrid,ngrid],r, hour, str(v))

            return(trans_var, x_kilometers)
    










