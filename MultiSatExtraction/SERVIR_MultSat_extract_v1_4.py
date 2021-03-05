# -*- coding: utf-8 -*-
"""
#%%############################################################## 
#                  IMERG / GSMAP EXTRACTION v 1.4
############################################################### 

@ author: SERVIR MEKONG 
@ correspondence M.A. LAVERDE-BARAJAS 
@ mlaverdeb@gmail.com
@ version  v 1.4

Variables:
    
parameters (P)

    P={}
    P['pixel_value'] = 0.1  # resolution
    P['boundary'] = [90,142,-12,35] extension IMERG #[90,118,-1,30] GSMAP #Boundary domain (MK) extended
    P['Area'] = 'MK'  Mekong basin 
    P['SRE_type'] =  SPP   IMERG* or GSMAP**
    P['Version'] =  version_type    30MIN_EARLY* 30MIN_LATE*,30MIN_FINAL*  NRT** NOW** NOWCAST*
    P['VRG'] = 0 # offline or 1 online
    P['GPM_dir'] = r'Local data folder'    
    P['time_type'] = 'hour' or 'day' (hourly agredated)
    
Date = datetime(P['year'], month, day, hour,min, 0)  

"""
import numpy as np
import netCDF4 as nc
import os
from datetime import datetime,date,time  
from datetime import timedelta
from ftplib import FTP


def extraction(P,Date):      
    
    Net_files = os.path.join(P['GPM_dir'],P['SRE_type'],P['Version']) 
    SRE_type = P['SRE_type']
    Version = P['Version']
    Boundaries= P['boundary']
    offline = P['VRG']
    time_type = P['time_type']
	
    # LOAD VRSGS database        
    Indir = '/VRG/'+ SRE_type+ '/'+Version       
    
    if not os.path.exists(Net_files):
        os.makedirs(Net_files) 
      
    # SET SYSTEM PARAMETERS
    # Boundaries GPM MK  DATABASE    
    # MASK LMB
    [Xmin_lbm,Xmax_lbm,Ymin_lbm,Ymax_lbm] = Boundaries # Boundary LMB    
    Lon_lmb = np.arange(Xmin_lbm,Xmax_lbm,0.1)
    Lat_lmb = np.arange(Ymin_lbm,Ymax_lbm,0.1)
    dx_c = len(Lon_lmb)
    dy_c = len(Lat_lmb)
    
    Error_file = open(os.path.join(Net_files,'Error_report.txt'),'a')  # error file     
    
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
        
        #%% VRG file 
        if SRE_type == 'IMERG':
            dayYear= (date(D1.year,D1.month,D1.day) - date(D1.year,1,1)).days + 1    
            VGR_Folder = '{0}/{1:03d}/{2:03d}/MK/'.format(Indir,Date.year,dayYear)
            Precip_variable = 'precipitationCal'
           
    #        ftp.cwd(VGR_Folder)            
            if Version == '30MIN_EARLY':    
                File = 'MK_3B-HHR-E.MS.MRG.3IMERG.{0:02d}{1:02d}{2:02d}-S{3:02d}{4:02d}00'.format(D1.year,D1.month,D1.day,D1.hour,D1.minute)
            elif Version == '30MIN_LATE': 
                File = 'MK_3B-HHR-L.MS.MRG.3IMERG.{0:02d}{1:02d}{2:02d}-S{3:02d}{4:02d}00'.format(D1.year,D1.month,D1.day,D1.hour,D1.minute)
            elif Version == '30MIN_FINAL':
                File = 'MK_3B-HHR.MS.MRG.3IMERG.{0:02d}{1:02d}{2:02d}-S{3:02d}{4:02d}00'.format(D1.year,D1.month,D1.day,D1.hour,D1.minute)
            else:
                File = 'IMERG version not found'
                print('IMERG version not found')
                
        if SRE_type == 'GSMAP':
            dayYear= (date(D1.year,D1.month,D1.day) - date(D1.year,1,1)).days + 1  
            VGR_Folder = '{0}/MK/{1:02d}/{2:02d}/{3:02d}/'.format(Indir,D1.year,D1.month,D1.day)           
            Precip_variable = 'precip'
            
            if Version == 'NOW':
                File = 'MK_gsmap_gauge_now.{0:02d}{1:02d}{2:02d}.{3:02d}{4:02d}'.format(D1.year,D1.month,D1.day,D1.hour,D1.minute)
            elif Version == 'NOWCAST':
                File = 'MK_gsmap_rnc_{0:02d}-{1:02d}-{2:02d}_{3:02d}z'.format(D1.year,D1.month,D1.day,D1.hour)
            elif Version == 'NRT':  
                Indir = '/VRG/'+ SRE_type+ '/RT' 
                VGR_Folder = '{0}/HOURLY/MK/{1:02d}/{2:02d}/{3:02d}/'.format(Indir,D1.year,D1.month,D1.day) 
                File = 'MK_gsmap_gauge.{0:02d}{1:02d}{2:02d}.{3:02d}{4:02d}'.format(D1.year,D1.month,D1.day,D1.hour,D1.minute)

        if offline == 1: 
            Files = os.listdir(os.path.join(Net_files,str(D1.year)))
            GPM = [f for f in Files  if f.startswith(File)]  
        else:
            if Version == 'NOWCAST':
                ftp = FTP('216.218.240.199')
                ftp.login(user='ftpuser',passwd=  '@Smekong')
                Files = ftp.nlst('gsmap_nowcast/')  

                GPM = [f for f in Files  if f.startswith('gsmap_nowcast/' + File)]    
                VGR_Folder='gsmap_nowcast/'
                               
            else:
                ftp = FTP('203.146.112.250')
                ftp.login(user='downloader',passwd='Down0000')
                Files = ftp.nlst(VGR_Folder) 
                GPM = [f for f in Files  if f.startswith(VGR_Folder + File)]           
               
        try:   
            temp=None             
            GPM_netcdf = os.path.join(Net_files,str(D1.year), GPM[0].replace(VGR_Folder, ''))
            if not os.path.exists(GPM_netcdf):
                ftpPath =  GPM[0]
                with open(GPM_netcdf, 'wb') as outfile:
                    ftp.retrbinary("RETR " + ftpPath, outfile.write)               

            f = nc.MFDataset(GPM_netcdf) # read file
            Lat = f.variables['lat'][:]
            Lon = f.variables['lon'][:]
            Lon_index = np.asanyarray(np.where( (Lon >= Xmin_lbm ) & (Lon <= Xmax_lbm) ) )
            Lat_index = np.asanyarray(np.where( (Lat >= Ymin_lbm ) & (Lat <= Ymax_lbm) ) )
            temp = f.variables[Precip_variable][0,Lon_index.min():Lon_index.max()+1,Lat_index.min():Lat_index.max()+1]
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
            temp = np.empty([dy_c,dx_c])
            temp.fill(0) # fill
            GPM_day = temp + GPM_day
            pass
        
        if Version == 'NRT':
            D1 = D1 + timedelta(hours=1)  
        else:
            D1 = D1 + timedelta(minutes=30)  
#        print(D1)    
   
    Error_file.close()     
    return GPM_day

 #
#if __name__ == '__main__':
#    GPM_extraction() 