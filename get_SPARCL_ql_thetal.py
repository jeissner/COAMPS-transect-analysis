#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 11:49:07 2023

@author: jordaneissner
"""
## Calculate liquid water potential temperature and total water mixing ratio profiles from ARM data


import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from ql_thetal import liquid_vars
import math


vdate = '20180701'
shour = '2330'
dsa = xr.open_dataset('enaarsclkazr1kolliasC1.c0.'+vdate+'.000000.nc', decode_times=False)
dsl = xr.open_dataset(vdate+'stratocumulus_data.cdf', decode_times=False)
dss = xr.open_dataset('enasondewnpnC1.b1.'+vdate+'.'+shour+'00.cdf', decode_times=False)

# open ARSCL data set
time2 = dsa.time_offset.values/3600.
top = dsa.cloud_layer_top_height.values
base = dsa.cloud_base_best_estimate.values

# open SPARCL dataset 
timee = dsl.rettime.values #(time) Seconds since 1/1/1970
time = (timee - timee[0])/3600.
cLWP = dsl.pamtraCLlwp.values #(time) retrieved cloud LWP  (mm)
pLWP = dsl.pamtraDlwp.values #(time) retrieved drizzle (below+above cloud base) LWP (mm)
tLWP = dsl.pamtralwp.values #(time) retrieved total lwp (cloud + drizzle) (mm)
   
# open sounding dataset 
time3 = dss.time_offset.values[0]/3600.
pres = dss.pres.values #Pa
tdry = dss.tdry.values + 273.15 #K
rh = dss.rh.values
alt = dss.alt.values
dp = dss.dp.values
wdir = dss.deg.values


# get 30 min average LWP around sounding launch time        
y = np.abs(time - (time3-0.25)).argmin()
z = np.abs(time - (time3+0.25)).argmin()
 
LWP = np.nanmean(cLWP[y:z])
lwp = LWP/.997 #mm

# get 39 min average cloud base and top heights around sounding launch time
w = np.abs(time2 - (time3-0.25)).argmin()
u = np.abs(time2 - (time3+0.25)).argmin()
top2 = top[w:u,0]
if np.isnan(top2).any() == False:
     CTH = np.nanmean(top[w:u,0])
     CBH = np.nanmean(base[w:u])
else:
     CTH = np.nan
     CBH = np.nan

j = np.abs(alt-CBH).argmin()
k = np.abs(alt-CTH).argmin() + 1

print(vdate, time3, CBH, CTH, lwp)

# get theta_l, q_l 
thetal_lwp,thetal_ad,theta,ql_lwp,ql_ad,wvapor,LWP_ad,F,CT = liquid_vars(lwp,CBH,CTH,alt,tdry,pres,rh)

if math.isnan(F)==True:
    print('NAN',str(F))
    F = str(F)
else:
    F = str(round(F,2))
    
if F == "--":
    F = str(np.nan)
    print('NAN',F)

fig, ((ax1,ax2)) = plt.subplots(1, 2, figsize = (10,5),sharex=False, sharey=False)
plt.subplots_adjust(wspace=0.1, hspace=0.15)

ax1.plot(theta,alt/1000.,label='Theta')
ax1.plot(thetal_lwp[j:k],alt[j:k]/1000.,label='Theta l - MWR')
ax1.plot(thetal_ad[j:k],alt[j:k]/1000.,label='Theta l - AD')
ax1.set_ylim(0,3)
ax1.set_xlim(280,300)
ax1.set_xlabel('Potential Temperature (K)')
ax1.set_ylabel('Height (km)')
ax1.legend(loc='upper left')
ax1.set_title(vdate + ' ' +str(time3)+ ' UTC F = ' + F)

ax2.plot(wvapor*1000,alt/1000.,label='qv')
ax2.plot(ql_lwp*10,alt/1000.,label='10*ql - MWR')
ax2.plot(ql_ad*10,alt/1000.,label='10*ql - AD')
ax2.set_ylim(0,3)
ax2.set_xlim(0,10)
ax2.set_xlabel('Mixing ratio (g/kg)')
ax2.legend(loc='upper left')

plt.show()