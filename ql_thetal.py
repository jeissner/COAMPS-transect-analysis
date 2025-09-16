# calculate q_l,theta_l from CBH, CTH, LWP

####input#####
# LWP in g/m2
# CBH,CTH in m
# height in m
# T in K
# theta in K
# pres in hPa
# rh in fraction

####output####
# ql using adiabatic liquid water lapse rate in g/kg
# ql using liquid water lapse rate from LWP in g/kg
# theta_l using adiabatic ql in K
# theta_l using LWP ql in K



import numpy as np
import math

def liquid_vars(LWP,CBH,CTH,alt,T,pres,rh):
     theta = T * (1000./pres)**0.286 #K
     # check for liquid in layer
     if math.isnan(LWP)==True or math.isnan(CBH)==True or math.isnan(CTH)==True:
         thetal_lwp = theta.copy()
         thetal_lwp[:] = 0
         thetal_ad = theta.copy()
         thetal_ad[:] = 0
         ql_lwp = theta.copy()
         ql_lwp[:] = 0
         ql_ad = theta.copy()
         ql_ad[:] = 0
         F = np.nan
     else:    
        ###constants
         Lv = 2492000 #J/kg
         Ra = 287 #J/kg/K
         Rv = 461 #J/kg/K
         g = 9.81 #m/s/s
         cp = 1005 #J/kg K
        
         # find indices of cloud base and top heights in sounding data
         x1 = (np.abs(CBH-alt)).argmin()
         x2 = (np.abs(CTH-alt)).argmin()
        
         es=6.11*np.exp((2.5*10.0**6./461.0)*((1./273.15)-(1./(T))))
         e=(rh/100.)*es
         qs=(0.622*es)/(pres-es)
         wvapor = (rh/100.)*qs # kg/kg
         q = (rh/100.)*qs # kg/kg
         rho = (pres*100.)/(T*Ra) #kg/m3
    
         # average density in the cloud layer
         rho_cav = np.average(rho[x1:x2])
    
         ql_ad = np.zeros(len(alt))
         ql_lwp = np.zeros(len(alt))
        
         H = Ra*T/g
         gamma_w = (g/cp) * ( 1 - ( (1+((Lv*qs)/(Ra*T))) / (1+((Lv**2*qs)/(cp*Rv*T**2))) ) )
         gamma_w = np.average(gamma_w[x1:x2])
         
         # adiabatic gamma_l
         gamma_ladd = ( (((0.622+qs)*qs*Lv)/(Ra*T**2))*gamma_w ) - ( (qs*pres)/((pres-es)*H) )
         gamma_lad = np.average(gamma_ladd[x1:x2]) #1.99E-6 #np.average(gamma_ladd[x1:x2]) # kg/kg/m

         #gamma_l from LWP
         gamma_ll = (2 * LWP) / (rho_cav * ((CTH - CBH)**2))  # g/kg/m
         
        
         # ql = 0 at cloud base and increases by gamma_l to cloud top
         for k in range(x1,x2):
            ql_ad[k] = (gamma_lad * (alt[k] - CBH))*1000. # g/kg
            ql_lwp[k] = gamma_ll  * (alt[k] - CBH) # g/kg 
           
            # print(k, gamma_ll,ql_ad[k],ql_lwp[k])
            
         thetal_lwp = theta - ((Lv/cp)*(ql_lwp/1000.)) # K
         thetal_ad = theta - ((Lv/cp)*(ql_ad/1000.)) # K
        
         # F = LWP from MWR / adiabatic LWP
         LWP_ad = (1000.*0.5*rho_cav*gamma_lad*((CTH - CBH))**2)
         F = LWP/LWP_ad
         CT = np.nanmean(T[x1:x2])
         # print( 'Gamma_l = ', gamma_lad*1.0E+6, ' g kg^-1 km^-1' )
         # print( 'Max LWC from Gamma_l, CTH, and CBH = ', 
         #   gamma_lad*1.0E+3*(CTH - CBH), ' g kg^-1' )
         # print('LWP_adi = ', 1000.*0.5*rho_cav*gamma_lad*((CTH - CBH))**2, 'g m^-2')
         # print('LWP_obs', LWP, 'g m^-2')
         # print('Adiabaticity Factor', F)
         # print('avg cloud temp', np.nanmean(pres[x1:x2]))

    
     return([thetal_lwp,thetal_ad,theta,ql_lwp,ql_ad,wvapor,LWP_ad,F,CT])
 
    
 