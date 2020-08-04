# -*- coding: utf-8 -*-
"""
BIAS CORRECTOR FOR GPM

@ author: SERVIR MEKONG 
@ correspondence M.A. LAVERDE-BARAJAS 
@ mlaverdeb@gmail.com
"""
import pandas as pd
import numpy as np
from osgeo import gdal,osr
from math import sqrt
from sklearn.metrics import mean_squared_error


try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
################################################################################# 
#              PROCESSING FUNCTIONS
################################################################################# 
def point_extraction(Raster,lyr,FID,Date): 
    
    # read Raster create file    
    gt = Raster.GetGeoTransform()
    Raster_rb = Raster.GetRasterBand(1) 
    
    Name_Sta = []
    Rain = []
    for a in range(len(lyr)):
        feat=lyr[a]
        geom = feat.GetGeometryRef()
        mx,my = geom.GetX(), geom.GetY()  #coord in map units
    
        #Convert from map to pixel coordinates.
    
        px = int((mx - gt[0]) / gt[1]) #x pixel
        py = int((my - gt[3]) / gt[5]) #y pixel
        
        intval = Raster_rb.ReadAsArray()
        intval = Raster_rb.ReadAsArray(px,py,1,1)  
        Name_Var = feat.GetField(FID)
        
        if isinstance(Name_Var, str):            
            Name_Sta.append(Name_Var)
        else:
            Name_Sta.append(str(int(Name_Var)))
        Rain.append(intval[0])
#        print intval[0]
        
    df_step = pd.DataFrame(np.array(Rain).T,columns= Name_Sta)
    df_step['Date'] = Date
    df_step = df_step.set_index('Date')        
    return df_step

def array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array):

    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_UInt32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()    
    return outRaster

def save_ascii(Raster_name,Correct,MinLon,MinLat):
    f = StringIO()
    np.savetxt(f,Correct, fmt='%.3f')
    f.seek(0)
    fs = f.read().replace('-9999.000', '-9999', -1)
    f.close()
    f = open(Raster_name, 'w')
    f.write("ncols " + str(Correct.shape[1]) + "\n")
    f.write("nrows " + str(Correct.shape[0]) + "\n")
    f.write("xllcorner " + str( MinLon) + "\n")
    f.write("yllcorner " + str(MinLat) + "\n")
    f.write("cellsize " + str(0.1) + "\n")
    f.write("NODATA_value " + str(-9999) + "\n")
    f.write(fs)
    f.close() 

## read shapefile
def read_shapefile(sf):

    fields = [x[0] for x in sf.fields][1:]    
    Records=[]   
    for r in range(len(sf.shapeRecords())):
        a=sf.records()[r][:]
        Records.append(a)
    shps = [s.points for s in sf.shapes()]
    df = pd.DataFrame(columns=fields, data=Records)
    df = df.assign(coords=shps)
    return df

def metrics(pred,obs):
    try:
        pred = pred.astype(float)
        obs = obs.astype(float)
        pred[pred <=1]  = np.nan
        obs[obs <=1]  = np.nan
        mask = ~np.isnan(pred) & ~np.isnan(obs)     
        R = np.corrcoef(pred[mask],obs[mask])
        R = np.round(R[1,0],2)
        RMSE = sqrt(mean_squared_error(pred[mask],obs[mask]))    
        BIAS = 100 * ( np.sum( pred[mask] - obs[mask] ) / np.sum( obs[mask] ) )
        return  R,RMSE,BIAS
    except:
        print('the Number of values is to low to calculate error metrics')
        pass

def plot_perform(SAT,OBS,Rain_valid_BIAS,Correct,GridSRE,Boundaries,Fig_name):
    
    import matplotlib.pyplot as plt 
    zmax = np.round(np.max([SAT,OBS,Rain_valid_BIAS[:,0],Rain_valid_BIAS[:,4]])/10)*10
    
    plt.figure(num=None, figsize=(12, 4), dpi=100, facecolor='w', edgecolor='k')
    
    plt.subplot(141)
    plt.imshow(GridSRE,extent=Boundaries,vmin=0,vmax=np.max([Correct,GridSRE]))
    plt.plot(Rain_valid_BIAS[:,1],Rain_valid_BIAS[:,2],'.r')
    plt.colorbar(shrink=0.5)
    plt.title(Fig_name[-12:-4])
    
    plt.subplot(143)
    plt.imshow(Correct,extent=Boundaries,vmin=0,vmax=np.max([Correct,GridSRE]))
    plt.colorbar(shrink=0.5)
    
    plt.subplot(142)
    mask=np.logical_and(SAT>=0 , OBS>=0)    
    plt.scatter(OBS[mask],SAT[mask]) 
    plt.plot([0,zmax],[0,zmax],'r')
    plt.xlim([0,zmax])
    plt.ylim([0,zmax])
    R= np.round(np.corrcoef(OBS[mask],SAT[mask])[0,1],2)
    plt.xlabel('R: {}'.format(R))
    plt.title('original')
    plt.gca().set_aspect('equal', adjustable='box')
    
    plt.subplot(144)
    mask=np.logical_and(Rain_valid_BIAS[:,0]>=0 , Rain_valid_BIAS[:,4]>=0)    
    plt.scatter(Rain_valid_BIAS[mask,0],Rain_valid_BIAS[mask,4]) 
    R = np.round(np.corrcoef(Rain_valid_BIAS[mask,0],Rain_valid_BIAS[mask,4])[0,1],2)
    plt.xlabel('R: {}'.format(R))
    plt.plot([0,zmax], [0,zmax],'r')
    plt.xlim([0,zmax])
    plt.ylim([0,zmax])
    plt.title('corrected')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(Fig_name)
    plt.close()
       
    
    
    
    
    
    
    