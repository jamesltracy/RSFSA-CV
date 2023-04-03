## This code aligns and clips all environmental and spatial grids to the Background Evaluation Area (BEE)
## grid
## IMPORTANT NOTE: The geoprocessing environment must have the processing extent set to match the species background polygon shapefile (e.g., Avirbckgpoly.shp)
## in order to restrict the number or grid rows and columns to the background area, and the output coordinates, snap raster and cell size set to match the
## Background Evaluation Extent raster (e.g., avirbck)

# make arcmap python commands and spatial analyst (sa) operations available
import arcpy
from arcpy.sa import *

import os

# Overwrite pre-existing files
arcpy.env.overwriteOutput = True

################################################  


# Set input and output directories
InDirect = "E:/TexasLayersWk"
OutDirect = "E:/TexasRoadEnvLayersInt"
RoadRasterName = "txroadmsk"

# Read in Road raster 
roadraster=Raster(InDirect + "/" + RoadRasterName)


###############################################
### For processing specific list of rasters 

# set working directory
arcpy.env.workspace = OutDirect

# Read in csv file of raster names as a Python list
import csv
rastersb = [x[0] for x in csv.reader(open\
 (InDirect + "/TexasRoadLayerGridNames.csv",\
     "rb"))]
#rastersb[0]


# Print list of rasters
for raster in rastersb:
    print(raster)

# Read in csv file of raster names as a Python list
import csv
rastersc = [x[0] for x in csv.reader(open\
 (InDirect + "/TexasRoadLayerGridNamesNew.csv",\
     "rb"))]
#rastersc[0]


# Print list of rasters
for raster in rastersc:
    print(raster)

## Adjusting raster lists for interruption:
#rastersb = rastersb[29:48]
#rastersc = rastersc[29:48]

## Loop  
count=0
for raster in rastersb:
    # Set the output name to be the same as the input name, and 
    #    locate in the 'out_workspace' workspace
    #
    output = os.path.join(OutDirect, rastersc[count])
    count = count + 1

    # Add each input raster in the list to background raster
    #
    myRaster = Raster(InDirect + "/" + raster)
    out1 = Int((myRaster + roadraster) * 10000)
    out1.save(output)
    del (myRaster, out1)

######
###########################################################
## For individual raster
#
out2 = Int(roadraster + Raster(InDirect + "/seadistb"))
out2.save(OutDirect + "/seadist")
#
out2 = Int((roadraster + Raster(InDirect + "/hrtxkdeb"))*10000)
out2.save(OutDirect + "/hrtxkde")
#
out2 = Int((roadraster + Raster(InDirect + "/popden1kb"))*10000)
out2.save(OutDirect + "/popden")
#
out2 = Int((roadraster + Raster(InDirect + "/mnpopden3kra"))*10000)
out2.save(OutDirect + "/mnpopden3kr")
#
out2 = Int((roadraster + Raster(InDirect + "/mnpopden9kra"))*10000)
out2.save(OutDirect + "/mnpopden9kr")
#
out2 = Int(roadraster + Raster(InDirect + "/hiurbdista"))
out2.save(OutDirect + "/hiurbdist")
#
out2 = Int(roadraster + Raster(InDirect + "/opnwatdista"))
out2.save(OutDirect + "/opnwatdist")

