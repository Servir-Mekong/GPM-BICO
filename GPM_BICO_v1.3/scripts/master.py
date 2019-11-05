# -*- coding: utf-8 -*-
"""
BIAS CORRECTOR FOR GPM
version 1.3

@ author: SERVIR MEKONG 
@ correspondence M.A. LAVERDE-BARAJAS 
@ mlaverdeb@gmail.com
"""
import pandas as pd
import numpy as np

import os
#import traceback
from datetime import datetime
from datetime import timedelta
from osgeo import gdal,ogr
import sys
import configparser
import warnings
warnings.filterwarnings("ignore")

################################################################################# 
#%%             Master script
################################################################################# 
def GPM_correction():     
    print('=============================================================')
    print('            bias correction tool (GPM_BICO)                  ')
    print('=============================================================')
    print('     Developed by SERVIR-MEKONG                              ')
    print('     version 1.3                                             ')
    print('     last updated (01/10/2019)                               ')
    print('     contact: miguel.barajas@adpc.net                        ')
    print('                                                             ')
    print('                                       Please wait   ...     ')
    ################################################################################# 
    #%              PARAMETERS
    ################################################################################# 

    cfg =  configparser.ConfigParser()
    cfg.read('ini_conf.cfg')
    #print(os.getcwd())
    IndDir = cfg.get('Folder','Input_Folder')  
                             # read working directory
    # import functions
    sys.path.append(os.path.join(IndDir,'Scripts'))
    
    import Processtools as process
    import SERVIR_GPMextract as sat_extract
    import BIAS_CORR as biascoor    
    
    # USER DEFINED PARAMETERS
    ModelType = int(cfg.get('Parameters','ModelType') )      # save in netcdf
    saveNETCDF = int(cfg.get('Parameters','Ouput_format') )  # save in netcdf
    Zero_Val = cfg.get('Parameters','Zero_Val') 
    Hydroformat = int(cfg.get('Parameters','HYDROMET') )
    
    # SYSTEM PARAMETERS 
    SRE_type = cfg.get('Inputs','SRE_type') 
    Version = cfg.get('Inputs','Version') 
    Stations = cfg.get('Inputs','Shp_Stations') 
    HYMOS_col = int(cfg.get('Inputs','HYMOS_ID_column') )
    RainGauge = cfg.get('Inputs','RainGauge')
    ErrorMetrics = int(cfg.get('Inputs','ErrorMetrics'))
    
    MinLon = int(cfg.get('Raster','MinLon'))
    MaxLon = int(cfg.get('Raster','MaxLon'))
    MinLat = int(cfg.get('Raster','MinLat'))
    MaxLat = int(cfg.get('Raster','MaxLat'))
    Boundaries=[MinLon,MaxLon,MinLat,MaxLat]   

    ################################################################################# 
    #%%               READ INPUTS
    #################################################################################   

    # shafile  MRC database 
    Stations_ds = ogr.Open(os.path.join('shapes',Stations))
    Stations_lyr = Stations_ds.GetLayer()    

    # Read MRC daily data     
    RG_matrix = pd.read_csv(os.path.join('RainData',RainGauge),header=None)
    try:
        DateStart = datetime.strptime( RG_matrix.loc[2,1].replace('-','/') , '%d/%m/%Y') #  LOCAL TIME (GMT + 7)
    except:
        DateStart = datetime.strptime( RG_matrix.loc[2,1].replace('-','/') , '%d/%m/%y') #  LOCAL TIME (GMT + 7)
    try:
        DateEnd  = datetime.strptime(RG_matrix.loc[3,1].replace('-','/') , '%d/%m/%Y') #  LOCAL TIME (GMT + 7)
    except:
        DateEnd  = datetime.strptime(RG_matrix.loc[3,1].replace('-','/') , '%d/%m/%y') #  LOCAL TIME (GMT + 7)
    
    Headers = RG_matrix.loc[6,:]
    RG = RG_matrix.loc[8:,:]
    RG[RG == Zero_Val] = 0
    RG = RG.rename(columns = Headers)
    
    # correct lag in HYDROMET data
    if Hydroformat == 1:
        DateStart = DateStart- timedelta(days=1)
        DateEnd = DateEnd - timedelta(days=1)
        RG['DATE'] =  pd.date_range(DateStart,DateEnd, freq='d')
    else:
        RG['DATE'] =  pd.date_range(DateStart,DateEnd, freq='d')
    RG = RG.set_index('DATE')
    del RG['StaID']
    
     # convert to numeric
    RG =  RG.apply(pd.to_numeric, errors='coerce')
    
    if ModelType == 0:
        Dates =  RG.index
    else:
        Today = datetime.strptime( datetime.today().strftime('%m/%d/%Y') , '%m/%d/%Y')
        Dates = pd.date_range(Today - timedelta(days = 3),
		Today - timedelta(days=2), freq='d')	# processing the day before
        
    # BIAS CORRECTION PARAMETERS
    Bias_method = int(cfg.get('Bias parameters','Method'))
    if Bias_method == 3:
        WindowsTime = int(cfg.get('SDT','WindowsTime'))
        Elv_Group_name = cfg.get('SDT','Elv_Group_name')
        
        if len(Dates) <= 1:
            Dwin = DateStart - timedelta(days=WindowsTime) # D2 change to D1
            DatesExtract =  pd.date_range(Dwin ,DateStart, freq='d')
        else:
            Dwin = DateStart - timedelta(days=WindowsTime) # D2 change to D1
            DatesExtract =  pd.date_range(Dwin ,DateEnd, freq='d')
    else:
        WindowsTime = 0
        if len(Dates) <= 1:
            DatesExtract = Dates
        else:
            DatesExtract =  pd.date_range(DateStart ,DateEnd, freq='d')
        
    minn = int(cfg.get('Bias parameters','Min_station'))   #  min stations
    # Average factor if number of stations is less than min     
    avgd = float(cfg.get('Bias parameters','avgd')) 
    
    mask = (RG.index >= DatesExtract[0]) & (RG.index <= DatesExtract[-1])
    RG = RG[mask]
    ################################################################################# 
    ##       EXTRACT DAILY GPM (aggregated from 30 min database) 
    ##       daily information is considered from 7:00 am to 7:00 pm
    ################################################################################# 
    print('=============================================================')
    print('    Extracting DAILY GPM (aggregated from 30 min database)      ')    
    # EXTRACT GPM dataset 
    Frame_df = []
    for Day in DatesExtract:
    #    Day=DatesExtract[1]        
        # EXTRACT GPM DAY RASTER FILE
        Raster = sat_extract.extraction(IndDir,SRE_type,Version,Day,Boundaries,saveNETCDF,'day')   
        # EXTRACT POINTS PER STATION
        MRC_step = process.point_extraction(Raster,Stations_lyr,HYMOS_col,Day)  
        Frame_df.append(MRC_step)
        
    SRE_points_list = pd.concat(Frame_df, axis=0)    
    SRE_points_list.to_csv( os.path.join(IndDir,'Outputs','SRE_points_MRC.csv')) 
    
    ################################################################################# 
    #%%                      GPM BIAS correction
    #################################################################################    
    print('=============================================================')
    Type = ['Distribution transformation','Spatial bias corrector',
            'Spatiotemporal Distribution transformation','Empirical Quantile mapping']
    print('  Correcting GPM daily data by the ' + Type[Bias_method] + ' method') 

    Frame_df = []
    for Day in Dates:           
    #    Day=Dates[0]
        print('GPM corrected ' + Day.strftime('%m/%d/%Y'))
        Rain_valid = biascoor.Valid(Stations,Day,WindowsTime,RG,SRE_points_list)
    
        OBS = Rain_valid['RG_{0}{1:02d}{2:02d}.asc'.format(Day.year,Day.month,Day.day)].values
        SAT = Rain_valid['SRE_{0}{1:02d}{2:02d}.asc'.format(Day.year,Day.month,Day.day)].values 
    
        Lon = Rain_valid.Lon.values 
        Lat = Rain_valid.Lat.values
        
        Raster = 'rain_{0}{1:02d}{2:02d}000000.asc'.format(Day.year,Day.month,Day.day)
        GridSRE = gdal.Open(os.path.join(IndDir,'Outputs','UnCorr'+Version[5:],Raster)).ReadAsArray()   
#        GridSRE = gdal.Open(os.path.join(IndDir,'Outputs','UnCorr',Raster))
#        geoTrans = GridSRE.GetGeoTransform()
        
        if len(Rain_valid) > minn:
            
            if Bias_method == 0:
                Correct = biascoor.DT(cfg,Day,OBS,SAT,GridSRE,0.1)
                
            elif Bias_method == 1:
                Interpolator= int(cfg.get('SB','Interp'))
                Interpolator = 3
                Correct = biascoor.SB(cfg,Day,OBS,SAT,GridSRE,Lon,Lat,Boundaries,Interpolator,avgd)  
                
            elif Bias_method == 2:
                filter_RG = [col for col in Rain_valid if col.startswith('RG')]
                OBSw = Rain_valid[filter_RG].values
                filter_SRE = [col for col in Rain_valid if col.startswith('SRE')]
                SATw = Rain_valid[filter_SRE].values
                Elv_Groups = gdal.Open(os.path.join(IndDir,'inputmaps',Elv_Group_name)).ReadAsArray()
                Cluster = Rain_valid['Group'].values
                Correct = biascoor.SDT(cfg,Dates,OBSw,SATw,GridSRE,Cluster,WindowsTime,Elv_Groups,minn) 
                
            elif Bias_method == 3:
                Correct = biascoor.EQM(cfg,OBS,SAT,GridSRE,RT=0.1)                
        else:
            Correct = biascoor.BiasDefault(cfg,GridSRE,avgd)
            
        ## SAVE DEFAULT ASC         
        Raster_name = os.path.join(IndDir,'Outputs','BIASCor'+Version[5:],'rain_{0}{1:02d}{2:02d}240000.asc'.format(Day.year,Day.month,Day.day))   
        process.save_ascii(Raster_name,Correct,MinLon,MinLat)
        Raster_save = gdal.Open(Raster_name)   
        MRC_step = process.point_extraction(Raster_save,Stations_lyr,HYMOS_col,Day)  
        Frame_df.append(MRC_step)
        
#     Save List    
    SRE_BIASCOR_list = pd.concat(Frame_df, axis=0)    
    SRE_BIASCOR_list.to_csv( os.path.join(IndDir,'Outputs','SRE_BIASCOR_list.csv'))    
   
    ################################################################################# 
    #%%                         CALCULATE ERROR METRICS
    ################################################################################# 
    if ErrorMetrics  == 1:
        print(' Calculating Error metrics')
        Mean_SRE_points = SRE_points_list.mask(SRE_points_list.eq(0)).mean(axis = 1, skipna = True)    
        Mean_BIASCOR = SRE_BIASCOR_list.mask(SRE_BIASCOR_list.eq(0)).mean(axis = 1, skipna = True)     
        Mean_RG = RG.mask(SRE_BIASCOR_list.eq(0)).mean(axis = 1, skipna = True)          
        
        corr = Mean_BIASCOR[DateStart:DateEnd].values
        obs = Mean_RG[DateStart:DateEnd].values
        org = Mean_SRE_points[DateStart:DateEnd].values
        Rorg,RMSEorg,BIASorg = process.metrics(org,obs)
        print('ORGINAL R= {0} RMSE = {1} BIAS  = {2}'.format(Rorg,np.round(RMSEorg,4),np.round(BIASorg,4)))
        Rcorr,RMSEcorr,BIAScorr = process.metrics(corr,obs)
        print('CORRECTED R= {0} RMSE = {1} BIAS  = {2}'.format(Rcorr,np.round(RMSEcorr,4),np.round(BIAScorr,4)))
#        import matplotlib.pyplot as plt
#        plt.plot(obs,org,'o')
    else:
        print('Error metrics not calculated')
        
    # Print output  out.log     
    sys.stdout = open(os.path.join(IndDir,'Outputs','out.log'), 'a')

################################################################################# 
    #%%     RUN ERROR METRICS
#################################################################################         
if __name__ == '__main__':
    GPM_correction() 
