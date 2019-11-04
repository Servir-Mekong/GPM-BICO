# GPM-BICO tool
Bias Correction Tool for GPM precipitation data

Operationalization of bias-corrected satellite observations via GPM IMERG products in the Mekong River Commission (MRC) Riverine Flood Forecasting System

Version 1.3


## HOW TO USE THE MODEL
### Model setup

GPM_BICO tool is a directory holding all the data needed to run the model. The Output
file is automatically created to allocate the raster files and csv files processed during
the simulation. The model contains the following folders:

● Scripts: Folder holding the scripts of the model
● RainData: Directory holding the rain gauges input from HYDROmet
● Inputmaps: Directory holdin the zonification raster file
● Shapes: Directory holding the rain station points in format *.shp
In addition, the model contains the following files

● Ini_config.cgf : The configuration setting for the GPM-BICO tool
● Runner.bat: Executable file to run the model

## Dependecies

● gdal

● numpy

● pandas

● netcdf4

● pyftpdlib

● pykrige

● scipy

● pyshp

● scikit-learn

● oi (install using pip e.g. pip install io)

This packages can be installed using anaconda via command prompt.

     conda env create -f environment.yml 
