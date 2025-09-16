import numpy as np
import metpy
from metpy.units import units
import quantities as pq

#T_0 = 1000 mb temp [K] - using 1000 mb instead of surface to avoid stable surface layers 
#rh0 = 1000 mb rh %
#T_850 = 850 mb temp [K]
#T_700 = 700 mb temp [K] 
#z700 = height of 700 mb [km]
#q1000 = wvapor at 1000 mb [kg/kg]

def stability_coamps(T_0, T_700, T_850, rh0, z_700, q1000):
    # constants
    Lv = 2492000 #J/kg
    Ra = 287 #J/kg/K
    Rv = 461 #J/kg/K
    g = 9.81 #m/s/s
    cp = 1005 #J/kg K
    p1000 = 100000.0 #Pa
    
    # lower tropospheric stability (Klein and Hartmann 1993)
    theta0 = (T_0)*(1000./p1000)**0.286
    theta700 = (T_700)*(1000./700.)**0.286
    LTS = theta700 - theta0
    
    # estimated inversion strength (Wood and Bretherton 2006)
    es = 6.11*np.exp((2.5*10.0**6./461.0)*((1./273.15)-(1./(T_850))))
    # e = (rh/100.)*es
    qs_850 = (0.622*es)/(850-es) # kg/kg
    
    ev1000 = p1000*q1000/(Ra/Rv+q1000)
    tempLCL =  2840.0/(3.5*np.log(T_0)-np.log(.01*ev1000)-4.805)+55.0
    
    #Calc LCL and height at 700 hPa, gets qv into kg/kg
    LCL = (-cp/g)*(tempLCL-T_0)
    gamma_850 = (g/cp) * ( 1 - ( (1+((Lv*qs_850)/(Ra*T_850))) / (1+((Lv**2*qs_850)/(cp*Rv*T_850**2))) ) )
    EIS = LTS - (gamma_850 * (z_700 - LCL))
    
    return([EIS,LTS])

def stability_sonde(pres, temp, rh, heightkm):
    # constants
    Lv = 2492000 #J/kg
    Ra = 287 #J/kg/K
    Rv = 461 #J/kg/K
    g = 9.81 #m/s/s
    cp = 1005 #J/kg K
    p1000 = 100000.0 #Pa

    i1000 = np.abs(pres - 1000.).argmin()
    i850 = np.abs(pres - 850.).argmin()
    i700 = np.abs(pres - 700.).argmin()
    
    T_0 = temp[i1000] + 273.15
    T_850 = temp[i850] + 273.15
    T_700 = temp[i700] + 273.15
    z_700 = heightkm[i700]

    
    theta0 = (T_0)*(1000./1000.)**0.286
    theta700 = (T_700)*(1000./700.)**0.286
    
    q1000 = metpy.calc.mixing_ratio_from_relative_humidity(pres[i1000] * units.hPa, temp[i1000] * units.degC, rh[i1000] * units.percent)
    qs_850 = metpy.calc.saturation_mixing_ratio(850 * units.hPa, T_850 * units.kelvin).to('kg/kg')

    
    LTS = theta700 - theta0
    
    ev1000 = p1000*q1000/(Ra/Rv+q1000)
    tempLCL =  2840.0/(3.5*np.log(T_0)-np.log(.01*ev1000)-4.805)+55.0
    LCL = (-cp/g)*(tempLCL-T_0)
    
    gamma_850 = (g/cp) * ( 1 - ( (1+((Lv*qs_850)/(Ra*T_850))) / (1+((Lv**2*qs_850)/(cp*Rv*T_850**2))) ) )
    EIS = LTS - (gamma_850 * (z_700 - LCL))
    
    return([EIS,LTS,LCL])





