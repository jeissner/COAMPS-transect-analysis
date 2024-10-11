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

import sys
sys.path.insert(0, '/Users/jordaneissner/Documents/clone/')
from transect_test import pointfind2
from transect_test import get_transect
from get_transects_funcs import cloudFrac,cloudBounds, get_sector, map_var_plot, line1, line2, bin_avg, bin_med


dir_string = "/Users/jordaneissner/Documents/wintercase1_KK_nr01nc100/"
sdate = '20180123'
hours = ['012','018','024','030','036','042','048', '054']
r = [100,115,160,195,230,255,270,290]
dtype = np.dtype('>f')
m1 = '384'
n1 = int(m1)
m2 = '397'
n2 = int(m2)
kka = 45

dsigma = [40.0, 60.0, 80.0, 80.0,
          80.0, 80.0, 80.0, 80.0, 80.0,
          80.0, 80.0, 80.0, 80.0, 80.0,
          80.0, 80.0, 80.0, 80.0, 80.0,
          100.0, 100.0, 100.0, 150.0, 200.0, 
          300.0, 400.0, 600.0, 800.0, 1000.0,
          1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
          1000.0, 1000.0, 1000.0, 1000.0, 1000.0,
          1000.0, 2500.0, 4200.0, 5800.0, 7500.0]

height = np.zeros(kka)
height[0] = 20.
for i, x in enumerate(dsigma):
    height[i+1] = height[i] + x
height = height[::-1]

file1 = dir_string+"longit_sfc_000000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_00000000_fcstfld"
file2 = dir_string+"latitu_sfc_000000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_00000000_fcstfld"
f = open(r"%s" %(file1), "rb")
lon = np.fromfile(f,dtype, count=-1)
lon = np.reshape(lon,  (n2,n2))

f = open(r"%s" %(file2), "rb")
lat = np.fromfile(f,dtype, count=-1)
lat = np.reshape(lat,  (n2,n2))

# define edge points of transect width
nt = 37 # number of points to make up the 1 deg width of transect
xs1 = np.linspace(324.5,325.5,nt)
ys1 = line1(xs1)

xs2 = np.linspace(336.5,337.5,nt)
ys2 = line2(xs2) 
nbins = 31
pn = 300
    
all_lwp = np.zeros([len(r),pn])
all_eis = np.zeros([len(r),pn])
all_pbl = np.zeros([len(r),pn])
all_lts = np.zeros([len(r),pn])
all_lflx = np.zeros([len(r),pn])
all_sflx = np.zeros([len(r),pn])
all_M = np.zeros([len(r),pn])
all_D1 = np.zeros([len(r),pn])
all_D2 = np.zeros([len(r),pn])
all_D3 = np.zeros([len(r),pn])
all_precip = np.zeros([len(r),pn])
all_tclf = np.zeros([len(r),pn])
all_qcfrac = np.zeros([len(r),pn])
all_xkm = np.zeros([len(r),pn])

for p in range(0,len(r)):
    hour = hours[p]
    ii = p+2
    file1 = dir_string+"airtmp_pre_000700_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file2 = dir_string+"relhum_pre_000850_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file3 = dir_string+"wvapor_sig_032480_000010_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file4 = dir_string+"geopht_pre_000700_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file5 = dir_string+"geopht_pre_000850_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file6 = dir_string+"airtmp_pre_000850_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file23 = dir_string+"airtmp_pre_000800_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file7 = dir_string+"airtmp_pre_001000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file8 = dir_string+"slpres_msl_000000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file9 = dir_string+"longit_sfc_000000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_00000000_fcstfld"
    file10 = dir_string+"latitu_sfc_000000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_00000000_fcstfld"
    file11 = dir_string+"relhum_sig_032480_000010_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file12 = dir_string+"cldlwp_sfc_000000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file15 = dir_string+"pottmp_sig_032480_000010_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file16 = dir_string+"wwwind_sig_036230_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file17 = dir_string+"totflx_sfc_000000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file18 = dir_string+"geopht_pre_001000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file19 = dir_string+"sehflx_sfc_000000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file20 = dir_string+"lahflx_sfc_000000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file21 = dir_string+"pblzht_sfc_000000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file24 = dir_string+"cldmix_sig_032480_000010_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file25 = dir_string+"ranmix_sig_032480_000010_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file26 = dir_string+"ttlpcp_sfc_000000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file27 = dir_string+"grdtmp_sfc_000000_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    # file28 = dir_string+"uuwind_zht_000010_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    # file29 = dir_string+"vvwind_zht_000010_000000_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file30 = dir_string+"cvcldf_sig_032480_000010_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    file31 = dir_string+"stcldf_sig_032480_000010_2a0"+m2+"x0"+m2+"_"+sdate+"00_0"+hour+"0000_fcstfld"
    
  
    
    f = open(r"%s" %(file1), "rb")
    T_700 = np.fromfile(f,dtype, count=-1)
    T_700 = np.reshape(T_700,  (n2,n2))
    
    f = open(r"%s" %(file2), "rb")
    rh850 = np.fromfile(f,dtype, count=-1)
    rh850 = np.reshape(rh850,  (n2,n2))
    
    f = open(r"%s" %(file3), "rb")
    qv = np.fromfile(f,dtype, count=-1)
    qv = np.reshape(qv,  (kka,n2,n2))
    
    f = open(r"%s" %(file4), "rb")
    ht700 = np.fromfile(f,dtype, count=-1)
    ht700 = np.reshape(ht700,  (n2,n2))
    
    f = open(r"%s" %(file5), "rb")
    ht850 = np.fromfile(f,dtype, count=-1)
    ht850 = np.reshape(ht850,  (n2,n2))
    
    f = open(r"%s" %(file6), "rb")
    T_850 = np.fromfile(f,dtype, count=-1)
    T_850 = np.reshape(T_850,  (n2,n2))
    
    f = open(r"%s" %(file7), "rb")
    T_0 = np.fromfile(f,dtype, count=-1)
    T_0 = np.reshape(T_0,  (n2,n2))
    
    f = open(r"%s" %(file8), "rb")
    slp = np.fromfile(f,dtype, count=-1)
    slp = np.reshape(slp,  (n2,n2))
    
    
    f = open(r"%s" %(file11), "rb")
    rh = np.fromfile(f,dtype, count=-1)
    rh = np.reshape(rh,  (kka, n2,n2))
    
    f = open(r"%s" %(file12), "rb")
    lwp = np.fromfile(f,dtype, count=-1)
    lwp = np.reshape(lwp,  (n2,n2)) * 1000.
    
    
    f = open(r"%s" %(file15), "rb")
    theta = np.fromfile(f,dtype, count=-1)
    theta = np.reshape(theta,  (kka, n2,n2))
    
    f = open(r"%s" %(file16), "rb")
    w = np.fromfile(f,dtype, count=-1)
    w = np.reshape(w,  (kka+1, n2,n2))
    
    
    f = open(r"%s" %(file18), "rb")
    ht1000 = np.fromfile(f,dtype, count=-1)
    ht1000 = np.reshape(ht1000,  (n2,n2))
    
    f = open(r"%s" %(file19), "rb")
    sflx = np.fromfile(f,dtype, count=-1)
    sflx = np.reshape(sflx,  (n2,n2))
    
    f = open(r"%s" %(file20), "rb")
    lflx = np.fromfile(f,dtype, count=-1)
    lflx = np.reshape(lflx,  (n2,n2))
    
    f = open(r"%s" %(file21), "rb")
    pblz = np.fromfile(f,dtype, count=-1)
    pblz = np.reshape(pblz,  (n2,n2))
    
    
    f = open(r"%s" %(file23), "rb")
    T_800 = np.fromfile(f,dtype, count=-1)
    T_800 = np.reshape(T_800,  (n2,n2))
    
    f = open(r"%s" %(file24), "rb")
    qc = np.fromfile(f,dtype, count=-1)
    qc = np.reshape(qc,  (kka, n2,n2))
    
    f = open(r"%s" %(file25), "rb")
    qr = np.fromfile(f,dtype, count=-1)
    qr = np.reshape(qr,  (kka, n2,n2))
    
    f = open(r"%s" %(file26), "rb")
    precip = np.fromfile(f,dtype, count=-1)
    precip = np.reshape(precip,  (n2,n2))
    
    f = open(r"%s" %(file27), "rb")
    sst = np.fromfile(f,dtype, count=-1)
    sst = np.reshape(sst,  (n2,n2))
    
    # f = open(r"%s" %(file28), "rb")
    # u = np.fromfile(f,dtype, count=-1)
    # u = np.reshape(u,  (n2,n2))
    
    # f = open(r"%s" %(file29), "rb")
    # v = np.fromfile(f,dtype, count=-1)
    # v = np.reshape(v,  (n2,n2))
    
    f = open(r"%s" %(file30), "rb")
    cf = np.fromfile(f,dtype, count=-1)
    cf = np.reshape(cf,  (kka, n2,n2))
    
    f = open(r"%s" %(file31), "rb")
    sf = np.fromfile(f,dtype, count=-1)
    sf = np.reshape(sf,  (kka, n2,n2))
    
    
    # X = lon2[::20,::20]
    # Y = lat2[::20,::20]
    # U = u[::20,::20]
    # V = v[::20,::20]
    
    
    theta700 = (T_700)*(1000./700.)**0.286
    theta0 = (T_0)*(1000./slp)**0.286
    theta850 = (T_850)*(1000./850.)**0.286
    theta800 = (T_800)*(1000./800.)**0.286
    theta_skin = (sst)*(1000./slp)**0.286
    
    es700 = 6.11*np.exp((2.5*10.0**6./461.0)*((1./273.15)-(1./(T_700))))
    
    qv1000 = np.zeros([n2,n2])
    rh1000 = np.zeros([n2,n2])
    for i in range(0,n2):
        for j in range(0,n2):
            i1000 = (np.abs(ht1000[i,j] - height)).argmin()
            qv1000[i,j] = qv[i1000,i,j]
            rh1000[i,j] = rh[i1000,i,j]
            
    tclf = cloudFrac(sf+cf, height, 300., 2000.)
    
    EIS, LTS = stability(T_0, T_700, T_850, rh1000, ht700, qv1000)
    lwp[lwp <= 1.] = np.nan
    
    cbh, cth, thick = cloudBounds(qc,height)
    
    M, D1,D2,D3 = CAO_decoup_params(T_0, qv1000,pblz,height,theta,qc,qr,qv,cbh,theta_skin,theta800) 

    ix = 12

    trans_lon = np.zeros(pn)
    trans_lat = np.zeros(pn)
    trans_eis = np.zeros(pn)
    trans_lts = np.zeros(pn)
    trans_sflx = np.zeros(pn)
    trans_lflx = np.zeros(pn)
    trans_lwp = np.zeros(pn)
    trans_pbl = np.zeros(pn)
    trans_M = np.zeros(pn)
    trans_D1 = np.zeros(pn)
    trans_D2 = np.zeros(pn)
    trans_D3 = np.zeros(pn)
    trans_precip = np.zeros(pn)
    trans_tclf = np.zeros(pn)
    trans_qcfrac = np.zeros(pn)
    trans_qcm = np.zeros(pn)
    
    lons = np.zeros([len(xs1),pn])
    lats = np.zeros([len(xs1),pn])
    lwps = np.zeros([len(xs1),pn])
    eiss = np.zeros([len(xs1),pn])
    pbls = np.zeros([len(xs1),pn])
    sflxs = np.zeros([len(xs1),pn])
    lflxs = np.zeros([len(xs1),pn])
    d1s = np.zeros([len(xs1),pn])
    prcps = np.zeros([len(xs1),pn])
    tclfs = np.zeros([len(xs1),pn])
    xkilos = np.zeros([len(xs1),pn])
    CCs = np.zeros([len(xs1),pn])
    qcs = np.zeros([25,len(xs1),pn])
    
    if ix > 0:
        for v in range(0,len(xs1)):
            data_prof, x_kilometers, m_grid, n_grid = get_transect(lat, lon, qc, ys1[v], xs1[v], ys2[v], xs2[v], npoints = pn, pdif = 1, norep=False)
            
            mgrid = [int(x) for x in m_grid]
            ngrid = [int(x) for x in n_grid]
            x_kilometers2 = x_kilometers - x_kilometers[r[p]] 

            print(v, np.shape(lwp[mgrid, ngrid]))
            lwps[v,:] = lwp[mgrid, ngrid]
            eiss[v,:] = EIS[mgrid, ngrid]
            pbls[v,:] = pblz[mgrid, ngrid]
            sflxs[v,:] = sflx[mgrid, ngrid]
            lflxs[v,:] = lflx[mgrid, ngrid]
            d1s[v,:] = D1[mgrid, ngrid]
            prcps[v,:] = precip[mgrid, ngrid]
            lons[v,:] = lon[mgrid, ngrid]
            lats[v,:] = lat[mgrid, ngrid]
            tclfs[v,:] = tclf[mgrid, ngrid]
            xkilos[v,:] = x_kilometers2
            # qcs[:,v,:] = qc[20:,mgrid,ngrid]
            
            for e in range(0, len(m_grid)):
                if data_prof[17:,e].any() > 0.00001 :
                    CCs[v,e] = 1.
                else:
                    CCs[v,e] = 0.

    
        for i in range(0,len(m_grid)):
            trans_lon[i] = np.nanmedian(lons[:,i])
            trans_lat[i] = np.nanmedian(lats[:,i])
            trans_eis[i] = np.nanmean(eiss[:,i])
            trans_sflx[i] = np.nanmean(sflxs[:,i])
            trans_lflx[i] = np.nanmean(lflxs[:,i])
            trans_lwp[i] = np.nanmean(lwps[:,i])
            trans_pbl[i] = np.nanmean(pbls[:,i])
            trans_D1[i] = np.nanmean(d1s[:,i])
            trans_precip[i] = np.nanmean(prcps[:,i])
            trans_tclf[i] = np.nanmean(tclfs[:,i])
            x_kilometers[i] = np.nanmedian(xkilos[:,i])
            # trans_qcm[i] = np.nanmean(qcs[:,:,i])
            trans_qcfrac[i] = np.nanmean(CCs[:,i])

    else: 
        data_prof, x_kilometers, m_grid, n_grid = get_transect(lat, lon, w, ys1[18], xs1[18], ys2[18], xs2[18], npoints = pn, pdif = 1, norep=False)
        mgrid = [int(x) for x in m_grid]
        ngrid = [int(x) for x in n_grid]

        for i in range(0,len(mgrid)):
            trans_lon[i] = lon[mgrid[i], ngrid[i]]
            trans_lat[i] = lat[mgrid[i], ngrid[i]]
            trans_eis[i] = EIS[mgrid[i], ngrid[i]]
            trans_lts[i] = LTS[mgrid[i], ngrid[i]]
            trans_sflx[i] = sflx[mgrid[i], ngrid[i]]
            trans_lflx[i] = lflx[mgrid[i], ngrid[i]]
            trans_lwp[i] = lwp[mgrid[i], ngrid[i]]
            trans_pbl[i] = pblz[mgrid[i], ngrid[i]]
            trans_M[i] = M[mgrid[i], ngrid[i]]
            trans_D1[i] = D1[mgrid[i], ngrid[i]]
            trans_D2[i] = D2[mgrid[i], ngrid[i]]
            trans_D3[i] = D3[mgrid[i], ngrid[i]]
            trans_precip[i] = precip[mgrid[i], ngrid[i]]
            trans_tclf[i] = tclf[mgrid[i], ngrid[i]]
        x_kilometers = x_kilometers - x_kilometers[r[p]] 
    
    file1 = '/Users/jordaneissner/Documents/MIKE_BAUERS_MCMSV4/RUNDIR/coamps/read_coamps/fronts/fronts_coamps_20183.nc'
    ncf = Dataset(file1)

    latc = ncf.variables['lat'][:]
    lonc = ncf.variables['lon'][:]
    cfc = ncf.variables['cf'][:] 
    wfc = ncf.variables['wf'][:] 
    timec = ncf.variables['time'][:]

    ncf.close()

    fronts = cfc[ii,:,:]*-10 + wfc[ii,:,:]*10
    fronts[~((fronts == 10) | (fronts <= -10))] = np.nan

    fig,((ax1,ax2)) = plt.subplots(2,1, figsize=(15,10), subplot_kw={'projection': ccrs.LambertConformal(central_latitude = 39, 
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


    ax1.contour(lon,lat,slp, transform=ccrs.PlateCarree(),c='k')
    CS = ax1.pcolormesh(lon,lat,lwp,transform=ccrs.PlateCarree(),cmap=cmap0,norm=LogNorm(vmin=1,vmax=1000))
    ax1.pcolormesh(lonc, latc, fronts, transform=ccrs.PlateCarree(), cmap='bwr',vmin=0.1,vmax=5)
    # ax1.quiver(X, Y, U, V, color='k', transform=ccrs.PlateCarree())
    ax1.plot(xs1,ys1,transform=ccrs.PlateCarree())
    ax1.plot(xs2,ys2,transform=ccrs.PlateCarree())
    ax1.plot(trans_lon,trans_lat,transform=ccrs.PlateCarree(),color='k', linewidth=4)
    ax1.scatter(trans_lon[r[p]],trans_lat[r[p]], transform=ccrs.PlateCarree(), c='r', s=45)
    cbar = plt.colorbar(CS, ax=ax1, orientation='vertical', cmap=cmap0, pad=0.04,extend='max')#, fraction=0.046, pad=0.04)
    cbar.set_label('LWP (g m$^{-2}$)',size=20)
    cbar.ax.tick_params(labelsize=18) 
    
    ax1.set_title(hours[p]+' HR')
    
    gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, x_inline=False, y_inline=False, linewidth=0.33, color='k',alpha=0.5)
    gl.right_labels = gl.top_labels = False
    gl.ylocator = ticker.FixedLocator([25,30,35,40, 45, 50, 55,60])
    gl.xlocator = ticker.FixedLocator([-40, -30, -20, -10])
    gl.xlabel_style = {'size': 16, 'color': 'k', 'rotation':45, 'ha':'right'}
    gl.ylabel_style = {'size': 16, 'color': 'k', 'weight': 'normal'}

    ax2.set_extent(extent)
    ax2.coastlines(resolution='50m', linewidth=0.4)

    ax1.contour(lon,lat,slp, transform=ccrs.PlateCarree(),c='k')
    CS = ax2.pcolormesh(lon,lat,tclf,transform=ccrs.PlateCarree(),vmin=0,vmax=1)
    ax2.pcolormesh(lonc, latc, fronts, transform=ccrs.PlateCarree(), cmap='bwr',vmin=0.1,vmax=5)
    ax2.plot(trans_lon,trans_lat,transform=ccrs.PlateCarree(),color='k', linewidth=4)
    ax2.scatter(trans_lon[r[p]],trans_lat[r[p]], transform=ccrs.PlateCarree(), c='r', s=45)
   
    cbar = plt.colorbar(CS, ax=ax2, orientation='vertical', cmap=cmap0, pad=0.04,extend='max')#, fraction=0.046, pad=0.04)
    cbar.set_label('Cloud Cover',size=20)
    cbar.ax.tick_params(labelsize=18) 
            
    # plt.savefig('COAMPS_front_transects_LWP_'+hours[p]+'.png', format='png', bbox_inches = 'tight', dpi=300) 

    plt.show()
    
    all_lwp[p,:] = trans_lwp
    all_eis[p,:] = trans_eis
    all_pbl[p,:] = trans_pbl
    all_lts[p,:] = trans_lts
    all_lflx[p,:] = trans_lflx
    all_sflx[p,:] = trans_sflx
    all_M[p,:] = trans_M
    all_D1[p,:] = trans_D1
    all_D2[p,:] = trans_D2
    all_D3[p,:] = trans_D3
    all_precip[p,:] = trans_precip
    all_tclf[p,:] = trans_tclf
    all_qcfrac[p,:] = trans_qcfrac
    all_xkm[p,:] = x_kilometers
    

bins, avglwp, stdlwp = bin_avg(all_xkm,all_lwp,nbins)
bins, medlwp, stdlwp = bin_med(all_xkm,all_lwp,nbins)
bins, medpbl, stdpbl = bin_med(all_xkm,all_pbl,nbins)
bins, medeis, stdeis = bin_med(all_xkm,all_eis,nbins)
bins, medsflx, stdsflx = bin_med(all_xkm,all_sflx,nbins)
bins, medlflx, stdlflx = bin_med(all_xkm,all_lflx,nbins)
bins, medlts, stdlts = bin_med(all_xkm,all_lts,nbins)
bins, medM, stdM = bin_med(all_xkm,all_M,nbins)
bins, medD1, stdD1 = bin_med(all_xkm,all_D1,nbins)
bins, medD2, stdD2 = bin_med(all_xkm,all_D2,nbins)
bins, medP, stdP = bin_med(all_xkm,all_precip,nbins)
bbins, medC, stdC = bin_avg(all_xkm,all_tclf,nbins)
bbins, medCF, stdCF = bin_avg(all_xkm,all_qcfrac,nbins)

# plot middle of bins 
bins = np.arange(-1400,1300,3000/(nbins-1))

fig,((ax1),(ax2),(ax3), (ax5)) = plt.subplots(4,1, figsize=(15,22))  
for j in range(0,len(r)):
    ax1.scatter(all_xkm[j,:], all_lwp[j,:], s=3.0, label=(hours[j]+' HR'),c='silver')
    ax1.plot(bins,medlwp, color='k', linewidth=4)
    ax1.set_ylabel('LWP (g m$^{-2}$)',size=24)
    ax1.tick_params(axis='y', labelsize=22)
    ax1.tick_params(axis='x', labelsize=22)
    ax1.plot([-150,-150],[0,1100],c='lightgrey')
    ax1.plot([150,150],[0,1100],c='lightgrey')
    ax1.set_ylim(1,300)
    # ax1.legend()
    
    ax8 = ax1.twinx() 
    ax8.scatter(all_xkm[j,:], all_tclf[j,:], s=2.0, color='lightblue')
    ax8.plot(bins,medC, color='cornflowerblue', linewidth=4)
    ax8.plot(bins,medCF, color='steelblue', linewidth=4)
    ax8.set_ylabel('Cloud Cover',size=24, color='cornflowerblue') 
    ax8.tick_params(axis='y', labelsize=22,labelcolor='cornflowerblue')
    ax8.set_ylim(0,1)
    
    
    
    ax2.scatter(all_xkm[j,:], all_sflx[j,:], c = 'lightcoral', s=2.0)#, label='Sensible')
    ax2.scatter(all_xkm[j,:], all_lflx[j,:], c = 'lightgreen', s=2.0)#, label='Latent')
    ax2.plot(bins,medsflx, color='r', linewidth=4)
    ax2.plot(bins,medlflx, color='g', linewidth=4)
    ax2.plot([-150,-150],[-10,450],c='lightgrey')
    ax2.plot([150,150],[-10,450],c='lightgrey')
    ax2.set_ylim(-10,450)
    ax2.legend(['Sensible','Latent'], fontsize=22, markerscale=3.0)
    ax2.set_ylabel('Surface Heat Flux (W m$^{-2}$)',size=24)
    ax2.tick_params(axis='y', labelsize=22)
    ax2.tick_params(axis='x', labelsize=22)

    ax3.tick_params(axis='y', labelcolor='b', labelsize=22)
    ax3.scatter(all_xkm[j,:], all_pbl[j,:]/1000., c='lightblue', s=2.0,zorder=1)
    ax3.set_ylim(0,2)
    ax3.tick_params(axis='x', labelsize=22)
    ax3.set_ylabel('PBL Height (km)', color='b',size=24)
    ax3.plot([-150,-150],[0,3],c='lightgrey')
    ax3.plot([150,150],[0,3],c='lightgrey')
    ax4 = ax3.twinx()
    ax4.tick_params(axis='y', labelcolor='orange', labelsize = 22)
    ax4.set_ylim(0,15)
    ax4.set_ylabel('EIS (K)', color='orange',size=24)
    ax4.scatter(all_xkm[j,:], all_eis[j,:], c='moccasin', s=2.0,zorder=2)
    ax3.plot(bins,medpbl/1000., color='b', linewidth=4,zorder=3)
    ax4.plot(bins,medeis, color='orange',linewidth=4,zorder=4)
    
    
    
    ax5.scatter(all_xkm[j,:], all_D2[j,:], c = 'paleturquoise', s=2.0)#, label='Latent')
    ax5.plot(bins,medD2, color='teal', linewidth=4)
    ax5.set_ylabel('D (K)',size=24,color='teal')
    ax5.tick_params(axis='y', labelsize=22,labelcolor='teal')
    ax5.tick_params(axis='x', labelsize=22)
    ax5.axhline(y=0.5, c='k')
    ax5.set_ylim(-1,2)
    ax5.plot([-150,-150],[-2,3],c='lightgrey')
    ax5.plot([150,150],[-2,3],c='lightgrey')

    ax7 = ax5.twinx()
    ax7.tick_params(axis='y', labelcolor='indigo', labelsize = 22)
    ax7.set_ylim(-0.1,5)
    ax7.set_ylabel('Hourly Precipitation (mm)', color='indigo',size=24)
    ax7.scatter(all_xkm[j,:], all_precip[j,:], c='indigo', s=2.0,zorder=2)
    # ax7.plot(bins,medP, color='maroon', linewidth=4,zorder=3)
    
    # ax5.scatter(all_xkm[j,:], height, all_M[j,:], c = 'lightcoral', s=2.0)
    
    ax5.set_xlabel('Distance from cold front (km)',size=26)

# plt.savefig('../COAMPS_plots/COAMPS_front_transects_all_fine_KKNr01_1degwidth.pdf', format='pdf', bbox_inches = 'tight', dpi=300) 
plt.show()


xms,yms = np.where(all_xkm[:,:] <= -150) 
x,y= np.where ((all_xkm[:,:] > -150) & (all_xkm[:,:] < 150))
xp,yp = np.where((all_xkm[:,:] >= 150)) 

postfrontsc = [all_lwp[xms,yms],all_pbl[xms,yms]/1000., all_eis[xms,yms], all_lflx[xms,yms], all_sflx[xms,yms]]
front = [all_lwp[x,y], all_pbl[x,y]/1000.,  all_eis[x,y],all_lflx[x,y], all_sflx[x,y]]
prefront = [all_lwp[xp,yp], all_pbl[xp,yp]/1000., all_eis[xp,yp], all_lflx[xp,yp], all_sflx[xp,yp]]

        
fig, axes = plt.subplots(1, 5, figsize=[12, 8])
for i, k in enumerate(postfrontsc):

    if i==0:
        k = k[~np.isnan(k)]
        front[i] = front[i][~np.isnan(front[i])]
        prefront[i] = prefront[i][~np.isnan(prefront[i])]
        
    data = [k, front[i], prefront[i]]
    

    bp = axes[i].boxplot(data,sym='black')
    plt.setp(bp['fliers'], markersize=0.1)
    axes[0].set_ylim(-10,1500)
    axes[1].set_ylim(0,3)
    axes[2].set_ylim(-3,13)
    axes[3].set_ylim(-10,450)
    axes[4].set_ylim(-20,220)
    
    
    left = [0.65,1.65,2.65]
    right = [1.45,2.45,3.45]

plt.tight_layout()
# plt.savefig('../COAMPS_plots/COAMPS_front_transects_boxplot_fine_KKNr01.pdf', format='pdf', bbox_inches = 'tight', dpi=300) 
     
    
    
    
    