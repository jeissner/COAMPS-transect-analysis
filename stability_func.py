import numpy as np

#T_0 = 1000 mb temp
#rh0 = 1000 mb rh
#T_850 = 850 mb temp
#T_700 = 700 mb temp
#z700 = height of 700 mb
#q1000 = wvapor at 1000 mb

def stability(T_0, T_700, T_850, rh0, z_700, q1000):
    
    theta0 = (T_0)*(1000./1000.)**0.286
   # theta850 = (T_850)*(1000./850.)**0.286
    theta700 = (T_700)*(1000./700.)**0.286
    
    es = 6.11*np.exp((2.5*10.0**6./461.0)*((1./273.15)-(1./(T_850))))
   # e = (rh/100.)*es
    qs_850 = (0.622*es)/(850-es) # kg/kg
    

    #%w=(0.622*e(ij))/(pres_met(ij)-e(ij));
    #q=(rh/100.)*qs
    #thetae = 
    Lv = 2500000 #2257000 #J/kg
    Ra = 287 #J/kg/K
    Rv = 461 #J/kg/K
    g = 9.81 #m/s/s
    cp = 1005 #J/kg K
    
    LTS = theta700 - theta0

    # zLCL = zref + cp/g * (T - 55K - ( (1/(T-55K)) - (ln(RH)/2840) )^-1 ) # Bolton 1980
    p1000 = 100000.0 #Pa
    p850  = 85000.0  #Pa
    p700  = 70000.0  #Pa
    
    # ev1000 = p1000*qv1000[i,*]/(Ra/Rv+v[i,*])
    # tempLCL[i,*] =  2840.0/(3.5*alog(temp1000[i,*])-alog(.01*ev1000[i,*])-4.805)+55.0
    
    ev1000 = p1000*q1000/(Ra/Rv+q1000)
    tempLCL =  2840.0/(3.5*np.log(T_0)-np.log(.01*ev1000)-4.805)+55.0
    
    #Calc LCL and height at 700 hPa, gets qv into kg/kg
    #LCL[i,*]  = (-cp/g)*(tempLCL[i,*]-temp1000[i,*])
    #****DAVE'S LCL CALC
    LCL = (-cp/g)*(tempLCL-T_0)
    
    
    #LCL = cp/g * (T_0 - 55 - ( (1/(T_0-55)) - (np.log(rh0/100)/2840) )**-1 )  #OLD BOLTON 1980
    gamma_850 = (g/cp) * ( 1 - ( (1+((Lv*qs_850)/(Ra*T_850))) / (1+((Lv**2*qs_850)/(cp*Rv*T_850**2))) ) )
    EIS = LTS - (gamma_850 * (z_700 - LCL))
    #EIS2 = LTS - (gamma_850 * (ht700 - LCL2))
    #print('EIS', EIS)
    
    return([EIS,LTS])