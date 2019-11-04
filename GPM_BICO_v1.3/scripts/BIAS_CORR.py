# -*- coding: utf-8 -*-
"""
Created on Wed May 22 10:14:01 2019

BIAS CORRECTOR tool
v1.3

@ author: SERVIR MEKONG 
@ correspondence M.A. LAVERDE-BARAJAS 
@ mlaverdeb@gmail.com

"""

import pandas as pd
import numpy as np
import os
import shapefile 
import scipy.stats as stats
from pykrige.ok import OrdinaryKriging
import interp_fun as Interp
import Processtools as process
from datetime import timedelta  
        
def Valid(Stations,Day,WindowsTime,RG,SRE_points_list):    
    
    ###   RainGauges Coordinates MRC
    sf1 = shapefile.Reader(os.path.join('shapes',Stations))
    Features = process.read_shapefile(sf1)
    
    HYMOS_ID = [str(name) for name in Features['HYMOS_ID'].astype(int)]
    Lat_Stations = Features['Y'].astype(float).values
    Lon_Stations = Features['X'].astype(float).values 
    Clust_group =  Features['RASTERVALU'].astype(float).values
    Coords = pd.DataFrame({'HYMOS_ID':HYMOS_ID,'Lon':Lon_Stations,'Lat':Lat_Stations,'Group':Clust_group})
    Coords = Coords.set_index('HYMOS_ID')
    
    #read time wih windows time
    Dstart = Day - timedelta(days=WindowsTime) # D2 change to D1  
    Date_win =  pd.date_range(Dstart ,Day, freq='d')

    # select valid SRE and RG
    RG_day =  pd.DataFrame()
    SRE_day = pd.DataFrame()
    for D in Date_win:
        RG_day['RG_{0}{1:02d}{2:02d}.asc'.format(D.year,D.month,D.day)] = RG.loc[D]
        SRE_day['SRE_{0}{1:02d}{2:02d}.asc'.format(D.year,D.month,D.day)] = SRE_points_list.loc[D]
        
    SRE_coord = Coords.join(SRE_day)
    RG_coord = RG_day.join(SRE_coord)
#    ColDate = 'RG_{0}{1:02d}{2:02d}.asc'.format(D.year,D.month,D.day)
#    Rain_valid = RG_coord.drop(RG_coord[~pd.notnull(RG_coord[ColDate])].index)    
    Rain_valid =RG_coord.dropna()
    
    return Rain_valid

## method 0: Distribution transformation
def BiasDefault(cfg,GridSRE,avgd):  
    
    cor_sre = GridSRE * avgd
    cor_sre[GridSRE <= 0] = GridSRE[GridSRE <= 0]
    cor_sre[cor_sre <= 0] = GridSRE[cor_sre <= 0]
        
    return cor_sre  

## method 1: Distribution transformation
def DT(cfg,Date,OBS,SAT,GridSRE,RT=0.1):  
    
#    print('GPM corrected by Distribution transformation method')
    # read default parameters    
    avgd = float(cfg.get('Bias parameters','avgd')) # Average factor if number of stations is less than min
    sdd = float(cfg.get('DT','sdd'))        # standard deviation by default
    maxavgf = float(cfg.get('DT','maxavgf'))            # maximum average factor
    maxsdf = float(cfg.get('DT','maxsdf'))          # maximum standard deviation factor      
 
    SAT_Inten = SAT.ravel()
    SAT_Inten = SAT_Inten[SAT_Inten>= RT] 
    OBS_Inten = OBS.ravel()
    OBS_Inten = OBS_Inten[OBS_Inten>= RT]
	
    # calculate statistics
    OBS_avg = np.mean(OBS_Inten)
    SAT_avg = np.mean(SAT_Inten)
    OBS_std = np.std(OBS_Inten)
    SAT_std = np.std(SAT_Inten)
    
    if (OBS_avg> 1) and (SAT_avg > 1):
        avgf = OBS_avg / SAT_avg
        if (avgf> maxavgf):
            avgf = maxavgf
        stdf = OBS_std / SAT_std
        if (stdf > maxsdf):
            stdf = maxsdf
    else:
        avgf = avgd
        stdf = sdd
    
    # distribution transformation equation    
    cor_sre = ( (GridSRE - SAT_avg ) * stdf ) + ( avgf *  SAT_std )            
    cor_sre[GridSRE <= 0] = GridSRE[GridSRE <= 0]
    cor_sre[cor_sre <= 0] = GridSRE[cor_sre <= 0]  

    return cor_sre        
    
 ## method 2: Spatial bias corrector    
def SB(cfg,Date,OBS,SAT,GridSRE,Lon,Lat,Boundaries,Interpolator,avgd):      
    
#    print('GPM corrected by Spatial bias corrector method')
    MinLon,MaxLon,MinLat, MaxLat = Boundaries
    # Read SRE Coordinates
    px = int((MaxLon - MinLon) / 0.1) #x pixel
    py = int((MaxLat - MinLat) / 0.1) #y pixel
    grid_Lat = np.linspace(MinLat,MaxLat, num=py)  # Array latitude
    grid_Lon = np.linspace(MinLon,MaxLon, num=px)  # Array longitude    
    
     # calculate statistics      
    OBS_avg = np.mean(OBS)
    SAT_avg = np.mean(SAT)            
    
    if (OBS_avg> 1) and (SAT_avg > 1):
        Loc_bias =  SAT-OBS 
        
        if Interpolator == 1:
            
            print( 'interpolating using linear radian basis function')
            grid1  = Interp.linear_rbf(Lon,Lat,Loc_bias,grid_Lon,grid_Lat)           
            bias_spa  = grid1.T 
            
        elif Interpolator == 2:
            print( 'interpolating using Ordinary Kriging function')
            ok1 = OrdinaryKriging(Lon,Lat,Loc_bias,variogram_model='spherical',nlags = 10,enable_plotting=False,verbose=False)
            Vok1,ssd1 = ok1.execute('grid',grid_Lon,grid_Lat)
            bias_spa = Vok1.data
            
        elif Interpolator >= 3:
            print( 'interpolating using IDW function')
            grid1  = Interp.iwd(Lon,Lat,Loc_bias,grid_Lon,grid_Lat,2)           
            bias_spa  = grid1.T 

        # bias correct equation 
        cor_sre = GridSRE +  bias_spa         
        cor_sre[GridSRE <= 0] = GridSRE[GridSRE <= 0]
        cor_sre[cor_sre <= 0] = GridSRE[cor_sre <= 0]              

    else:
        print('average OBS and SAT lower than 1 correct with default')
        cor_sre = BiasDefault(cfg,GridSRE,avgd)

    return(cor_sre)

# method 3: Spatiotemporal distribution transformation
def SDT(cfg,Dates,OBSw,SATw,GridSRE,Cluster,WindowsTime,Elv_Groups,minn):   
    
#    print('GPM corrected by Spatiotemporal Distribution transformation method')
    # read default parameters    
    avgd = float(cfg.get('Bias parameters','avgd')) # Average factor if number of stations is less than min
#    sdd = float(cfg.get('DT','sdd'))        # standar deviation by default
    maxavgf = float(cfg.get('DT','maxavgf'))            # maximun average factor
    maxsdf = float(cfg.get('DT','maxsdf'))          # maximum standar deviation factor
   
    flag= 0
    Zone = 1
    SubZone = Zone       

    cor_sre = np.zeros(GridSRE.shape)
    while Zone <= 3:   
        
        cl = np.isin(Cluster,SubZone)
        
        OBS_ = OBSw[cl,:]
        SAT_ = SATw[cl,:]
        
        if SAT_.shape[0] > minn:
            
            OBS_flat = OBS_.flatten()
            SAT_flat = SAT_.flatten()
            mask = ~np.isnan(OBS_flat)
            
            OBS_flat = OBS_flat[mask]           
            SAT_flat = SAT_flat[mask]
            
            # calculate statistics
            OBS_avg = np.mean(OBS_flat)
            SAT_avg = np.mean(SAT_flat)
            OBS_std = np.std(OBS_flat)
            SAT_std = np.std(SAT_flat)
                             
            if (OBS_avg> 1) and (SAT_avg > 1) :     
                
                avgf = OBS_avg / SAT_avg
                stdf = OBS_std / SAT_std
                
                # threshold mean standard deviation
                if (avgf> maxavgf):
                    avgf = maxavgf            
                if (stdf > maxsdf):
                    stdf = maxsdf
                
                Area_Zone = np.isin(Elv_Groups,SubZone)           
                cor_sre[Area_Zone] = ( (GridSRE[Area_Zone] - SAT_avg ) * stdf ) + ( avgf *  SAT_std ) 
                
                Zone = Zone +1
                SubZone = Zone 
                flag= 1
            else:
                Zone = Zone +1
                SubZone = np.append(SubZone, Zone) 
        else:            
            Zone = Zone +1
            SubZone = np.append(SubZone, Zone) 
        
        cor_sre[GridSRE <= 0] = GridSRE[GridSRE <= 0]
        cor_sre[cor_sre <= 0] = GridSRE[cor_sre <= 0]
    
    if flag == 0:
        print('average OBS and SAT lower than 1 correct with default')
        cor_sre = BiasDefault(cfg,GridSRE,avgd)
        cor_sre[GridSRE <= 0] = GridSRE[GridSRE <= 0]
        cor_sre[cor_sre <= 0] = GridSRE[cor_sre <= 0]         
 
    return(cor_sre)
 
# method 4: Empirical Quantile mapping
def EQM(cfg,OBS,SAT,GridSRE,RT=0.1):
 
#    print('GPM corrected by Empirical Quantile mapping')
   
    SAT_Inten = SAT.ravel()
    SAT_Inten = SAT_Inten[SAT_Inten >= RT] 
    OBS_Inten = OBS.ravel()
    OBS_Inten = OBS_Inten[OBS_Inten >= RT]
    
    if (np.mean(OBS_Inten) > RT) and (np.mean(SAT_Inten) > RT):
        cdfn=10            
        def map(vals,CDF_OBS,CDF_MOD,zbins):
            # calculate exact CDF values using linear interpolation
            cdf1=np.interp(vals,zbins,CDF_MOD,left=0.0,right=999.0)
            # now use interpol again to invert the obsCDF, hence reversed x,y
            corrected=np.interp(cdf1,CDF_OBS,zbins,left=0.0,right=-999.0)        
            return corrected    
        
        fit_alphaY, fit_locY, fit_betaY = stats.gamma.fit(SAT_Inten)
        fit_alphaO, fit_locO, fit_betaO = stats.gamma.fit(OBS_Inten)
        
        global_max = max(np.amax(OBS_Inten),np.amax(SAT_Inten))
        wide = global_max/cdfn
        zbins = np.arange(0.0,global_max+wide,wide)
        
        cdf_gammaMOD = stats.gamma.cdf(zbins, fit_alphaY, fit_locY, fit_betaY)
        cdf_gammaOBS = stats.gamma.cdf(zbins, fit_alphaO, fit_locO, fit_betaO)        
           
        cor_values = map(GridSRE.ravel(),cdf_gammaOBS,cdf_gammaMOD,zbins)  
          
        cor_sre = SAT.ravel() 
        cor_sre[cor_sre>= RT] = cor_values 
        cor_sre[cor_sre <= 0]  = 0 
        cor_sre = np.reshape(cor_sre,SAT.shape)
        
    else:
        avgd = float(cfg.get('Bias parameters','avgd'))
        print('average OBS and SAT lower than 1 correct with default')
        cor_sre = BiasDefault(cfg,GridSRE,avgd)
		
    return cor_sre
   
    
    