import numpy as np

#T_0 = 1000 mb temp
#rh0 = 1000 mb rh
#T_850 = 850 mb temp
#T_700 = 700 mb temp
#z700 = height of 700 mb
#slp = sea level pressure
#q1000 = wvapor at 1000 mb
#CBH,PBLZ, height in meters
#qc,qr,qv in kg/kg

# for COAMPS analysis
def CAO_decoup_params(T_0,q1000,PBLZ,height,theta,qc,qr,qv, CBH, theta_skin, theta_800):
    Lv = 2257000 #J/kg
    Ra = 287 #J/kg/K
    Rv = 461 #J/kg/K
    g = 9.81 #m/s/s
    cp = 1005 #J/kg K
    p1000 = 100000.0 #Pa
    
    ql = (qc + qr)
    thetal = theta - ((Lv/cp)*(ql))
    qT = (ql + qv)*1000 
    D1 = CBH.copy()
    D1[:] = 0
    D2 = CBH.copy()
    D2[:] = 0
    D3 = CBH.copy()
    D3[:] = 0
    
    for i in range(0,np.shape(CBH)[0]):
        for j in range(0,np.shape(CBH)[0]):
           
            x = np.where(height <= PBLZ[i,j])
            
            if PBLZ[i,j] > 20.: 
                pbl = height[x]
                pbli = x[0][0]
                up = int(np.percentile(x,25))
                low = int(np.percentile(x,75))
                
               # print(PBLZ[i,j],pbl,low, height[low], up, height[up])
                
                thetal_top = np.mean(thetal[pbli:up,i,j])
                thetal_bot = np.mean(thetal[low:,i,j])
                qT_top = np.mean(qT[pbli:up,i,j])
                qT_bot = np.mean(qT[low:,i,j])
                
                D1[i,j] = thetal_top - thetal_bot
                D2[i,j] = qT_bot - qT_top
    
    #Jones et al. (2011) decoupling #1
    #deltaq = q_bot - q_top  > 0.5 g/kg = decoupled
    #deltatheta = thetal_top - thetal_bot > 0.5 K = decoupled
    #bottom and top are mean over lower and upper 25% of the boundary lyaer below inversion
    
    
    #Jones et al. (2011) decoupling #2
    #difference between LCL and CBH
    #deltaz = z_cb - z_lcl > 150 m = decoupled 
    # if i = 2:
        
    y = np.where(CBH > 0)
    ev1000 = p1000*q1000/(Ra/Rv+q1000)
    tempLCL =  2840.0/(3.5*np.log(T_0)-np.log(.01*ev1000)-4.805)+55.0
    LCL = (-cp/g)*(tempLCL-T_0)

    D3[y] = CBH[y] - LCL[y]
    #print(CBH[y], LCL[y], D3[y])
    
    #Jones et al. (2011) decoupling #3
    # mixed layer cloud thickness: difference between inversion height and LCL
    # when the MBL is well-mixed MLCT  is equivalent to cloud thickness
    # deltaz = z_i - z_LCL 
    
    
    
    #Marine cold air outbreak parameter (Fletcher et al. 2016)
    # M > 0 = CAO absolutely unstable
    M = theta_skin - theta_800


    return(M,D1,D2,D3)


#for sounding analysis
def decouple(T_0,q1000,PBLZ,theta,ql,qv, CBH, height):
    Lv = 2257000 #J/kg
    Ra = 287 #J/kg/K
    Rv = 461 #J/kg/K
    g = 9.81 #m/s/s
    cp = 1005 #J/kg K
    p1000 = 100000.0 #Pa
    

    thetal = theta - ((Lv/cp)*(ql))
    qT = (ql + qv)*1000 
    

            
    if PBLZ > 20.: 
        x = np.where(height <= PBLZ)
        pbl = height[x]
        pbli = x[-1]
        up = int(np.percentile(x,75))
        low = int(np.percentile(x,25))
        
        # print(PBLZ,pbl,low, height[low], up, height[up])
        thetal_top = np.mean(thetal[up])
        thetal_bot = np.mean(thetal[low])
        qT_top = np.mean(qT[up])
        qT_bot = np.mean(qT[low])
        
        D1 = thetal_top - thetal_bot
        D2 = qT_bot - qT_top

    #Jones et al. (2011) decoupling #1
    #deltaq = q_bot - q_top  > 0.5 g/kg = decoupled
    #deltatheta = thetal_top - thetal_bot > 0.5 K = decoupled
    #bottom and top are mean over lower and upper 25% of the boundary lyaer below inversion
    
    
    #Jones et al. (2011) decoupling #2
    #difference between LCL and CBH
    #deltaz = z_cb - z_lcl > 150 m = decoupled 
    # if i = 2:
    D3 = np.nan
    
    if CBH > 0:
        ev1000 = p1000*q1000/(Ra/Rv+q1000)
        tempLCL =  2840.0/(3.5*np.log(T_0)-np.log(.01*ev1000)-4.805)+55.0
        LCL = (-cp/g)*(tempLCL-T_0)
    
        D3 = CBH - LCL
        # print(CBH, LCL, D3)
    
    #Jones et al. (2011) decoupling #3
    # mixed layer cloud thickness: difference between inversion height and LCL
    # when the MBL is well-mixed MLCT  is equivalent to cloud thickness
    # deltaz = z_i - z_LCL 
    
    
    
    #Marine cold air outbreak parameter (Fletcher et al. 2016)
    # M > 0 = CAO absolutely unstable
    # M = theta_skin - theta_800
    return(D1,D2,D3)