# -*- coding: utf-8 -*-
"""
#%%############################################################## 
#                  IMERG EXTRACTION v 1.3.5
############################################################### 

@ author: SERVIR MEKONG 
@ correspondence M.A. LAVERDE-BARAJAS 
@ mlaverdeb@gmail.com
"""
import numpy as np
import netCDF4 as nc
import os
from datetime import datetime,date,time  
from datetime import timedelta
from ftplib import FTP
from osgeo import gdal
#import sys

import Processtools as process
# extraction point per station

#Date = Day
#Boundaries=[98,110,8,23]
# time_type= 'day'
def extraction(IndDir,SRE_type,Version,Date,Boundaries,saveNETCDF,time_type,offline):      
    
    try:    
        # LOAD VRSGS database	
        ftp = FTP('203.146.112.250')
        ftp.login(user='downloader',passwd='Down0000')
        Indir = '/VRG/'+ SRE_type+ '/'+Version
        vm = 0
    except:
        ftp = FTP('216.218.240.199')
        ftp.login(user='ftpuser',passwd=  '@Smekong')
        Indir = '/home/ftpuser/wrf_precip/imerg_30m_early' 
        print('FTP 203.146.112.250 unavailable... switching to backup server mode')
        vm = 1 
    
    # create Ouputs folder
    Output = os.path.join(IndDir,'Outputs')
    Net_files = os.path.join(Output,Version)
    Un_files = os.path.join(IndDir,'Outputs','UnCorr'+Version[5:])
    Cor_files = os.path.join(IndDir,'Outputs','BIASCor'+Version[5:])
    if not os.path.exists(Output):
        os.mkdir(Output)
    if not os.path.exists(Un_files):   
        os.mkdir(Un_files)
    if not os.path.exists(Cor_files):   
        os.mkdir(Cor_files) 
    if not os.path.exists(Net_files):    
        os.mkdir(Net_files)
    
    # SET SYSTEM PARAMETERS
    # Boundaries GPM MK  DATABASE    
    # MASK LMB
    [Xmin_lbm,Xmax_lbm,Ymin_lbm,Ymax_lbm] = Boundaries # Boundary LMB    
    Lon_lmb = np.arange(Xmin_lbm,Xmax_lbm,0.1)
    Lat_lmb = np.arange(Ymin_lbm,Ymax_lbm,0.1)
    dx_c = len(Lon_lmb)
    dy_c = len(Lat_lmb)
    
    Error_file = open(os.path.join(Output,'Error_report.txt'),'a')  # error file     
    
    if time_type == 'hour':
        # Adjust daily data to local time zone (GMT+7)
        D1 = Date + timedelta(hours=7) 
        D2 = D1 + timedelta(hours=1)     
    else:
        # Adjust daily data to local time zone (GMT+7) and collection time (+7 hours)
        D1 = datetime.combine(date(Date.year,Date.month,Date.day), time(14,0))
        D2 = D1 + timedelta(hours=24) 
        
    GPM_day=np.zeros([dy_c,dx_c])
    
    while D1 < D2:         
        
        if vm == 0:    # VM VRG '203.146.112.250'
            VGR_Folder = '{0}/{1:03d}/{2:03d}/MK/'.format(Indir,Date.year,dayYear)
        else:           # VM 216.218.240.199
            VGR_Folder = '{0}/{1:03d}/{2:03d}/'.format(Indir,Date.year,dayYear)
#        ftp.cwd(VGR_Folder)
        
        if Version == '30MIN_EARLY':    
            File = 'MK_3B-HHR-E.MS.MRG.3IMERG.{0:02d}{1:02d}{2:02d}-S{3:02d}{4:02d}00'.format(D1.year,D1.month,D1.day,D1.hour,D1.minute)
        elif Version == '30MIN_LATE': 
            File = 'MK_3B-HHR-L.MS.MRG.3IMERG.{0:02d}{1:02d}{2:02d}-S{3:02d}{4:02d}00'.format(D1.year,D1.month,D1.day,D1.hour,D1.minute)
        elif Version == '30MIN_FINAL':
            File = 'MK_3B-HHR.MS.MRG.3IMERG.{0:02d}{1:02d}{2:02d}-S{3:02d}{4:02d}00'.format(D1.year,D1.month,D1.day,D1.hour,D1.minute)
        else:
            print('IMERG version not found')

        if offline == 1: 
            Files = os.listdir(Net_files)
            GPM = [f for f in Files  if f.startswith(File)]  
        else:
            Files = ftp.nlst(VGR_Folder) 
            GPM = [f for f in Files  if f.startswith(VGR_Folder + File)]           
               
        try:   
            temp=None             
            GPM_netcdf = os.path.join(Net_files, GPM[0].replace(VGR_Folder, ''))
            if not os.path.exists(GPM_netcdf):
                ftpPath =  GPM[0]
                with open(GPM_netcdf, 'wb') as outfile:
                    ftp.retrbinary("RETR " + ftpPath, outfile.write)               

            f = nc.MFDataset(GPM_netcdf) # read file
            Lat = f.variables['lat'][:]
            Lon = f.variables['lon'][:]
            Lon_index = np.asanyarray(np.where( (Lon >= Xmin_lbm ) & (Lon <= Xmax_lbm) ) )
            Lat_index = np.asanyarray(np.where( (Lat >= Ymin_lbm ) & (Lat <= Ymax_lbm) ) )
            temp = f.variables['precipitationCal'][0,Lon_index.min():Lon_index.max()+1,Lat_index.min():Lat_index.max()+1]
            temp = np.flipud(temp.T)
            f.close()          
          
            GPM_day = temp + GPM_day
            
        except Exception:                
           
            if len(GPM) == 0:
                print('no ' + str(D1))
                Error_file.write( File +' no found\n') 
            else:
                print(GPM[0] + " format invalid")
                Error_file.write( File +" format invalid\n")
            # value of error
#            D1 = D2 + timedelta(minutes=30) # D2 change to D1
            GPM_day = np.empty([dy_c,dx_c])
            GPM_day.fill(0) # fill
            pass
        
        D1 = D1 + timedelta(minutes=30)  
#        print(D1)    
    #%% ########### SAVE RASTER AS:   ###########    
                    
    if saveNETCDF == 1  :        
        # SAVE GridFIle as in a NETCDF file
        Raster_name = os.path.join(Output,'UnCorr'+Version[5:],'IMERG_MRC_{0}{1:02d}{2:02d}.nc'.format(Date.year,Date.month,Date.day))
        dataset = nc.Dataset(Raster_name,'w', format='NETCDF4') 
        
        timeo = dataset.createDimension('time',None)
        timeo = dataset.createVariable('time','f4',('time'))
        timeo.units = 'Date{0}{1:02d}{2:02d}'.format(Date.year,Date.month,Date.day)
        timeo.standard_name = 'time'
        
        latitudes = dataset.createDimension('lat', dy_c)
        latitudes = dataset.createVariable('latitudes', np.float32, ('lat',))
        latitudes.standard_name = 'latitude'
        latitudes.units = 'degree_north'
        latitudes.axis = "Y"
        latitudes[:] = Lat_lmb 
                            
        longitudes = dataset.createDimension('lon', dx_c)
        longitudes = dataset.createVariable('longitudes', np.float32,('lon',)) 
        longitudes.standard_name = 'longitude'
        longitudes.units = 'degree_east' 
        longitudes.axis = "X"
        longitudes[:] = Lon_lmb        
              
        acc_precip = dataset.createVariable('precipitationCal', 'f4', ('time', 'lat', 'lon'),
        zlib=True, complevel=9, least_significant_digit=1, fill_value=-9999)
        acc_precip[0,:,:] = GPM_day
        acc_precip.standard_name = 'precipitationCal'
        acc_precip.units = 'mm'
        acc_precip.setncattr('grid_mapping', 'spatial_ref')
    
        crs = dataset.createVariable('spatial_ref', 'i4')
        crs.spatial_ref='GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'     
         
        dataset.close()
        
    else:
        ## SAVE DEFAULT ASC 
        Raster_name = 'rain_{0}{1:02d}{2:02d}{3:02d}0000.asc'.format(Date.year,Date.month,Date.day,Date.hour)
        Raster_path = os.path.join(Output,'UnCorr'+Version[5:],Raster_name)
    #            Raster = array2raster(Raster_name,(Min_lon,Max_lat),0.1,-0.1,GPM_day)   COMPLETE            
#        Raster = process.array2raster(Raster_path,(Xmin_lbm,Ymax_lbm),0.1,-0.1,GPM_day_clip)  # LMB
        process.save_ascii(Raster_path,GPM_day,Xmin_lbm,Ymin_lbm)
        Raster = gdal.Open(Raster_path)
        print( Raster_name + ' File created')
    
    Error_file.close()   
    
    return Raster

 #
#if __name__ == '__main__':
#    GPM_extraction() 
