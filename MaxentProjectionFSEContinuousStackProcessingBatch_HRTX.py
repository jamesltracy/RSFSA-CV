## This code assembles selected MaxEnt models and creates a sum consensus Feature Subset Ensemble (FSE) model for
## the binary output models. Models are added to a list and then stacked
## NOTE: The geoprocessing environment must have the output coordinates, processing extent, snap raster and cell size set to match the mask
## input environmental grids used in building the MaxEnt model

# make arcmap python commands and spatial analyst (sa) operations available
import arcpy
from arcpy.sa import *

import os

# Overwrite pre-existing files
arcpy.env.overwriteOutput = True

################################################  
# set working directory
#arcpy.env.workspace = "C:/Users/james/Documents/ArcGIS/Default.gdb"
arcpy.env.scratchWorkspace = arcpy.env.workspace

# Set input and output directories
Species = "HerpRKTX10_21"
HabitatSpecies = "HerpTX10_21"
Species2 = "HRTX"
Rank = "AIC"
InDirectPre1 = "E:/MaxEntOutput/"
InDirectPre = InDirectPre1 + Species + "MaxentAll_RSFSA55Varsof55/psa_corrfilt0.7/" + Rank + "Ranked/"
InDirect1 = "E:/TexasRoadEnvLayersInt/"
OutDirectName = Species + "_Models" + Rank
OutDirect = InDirectPre1 + OutDirectName + "/"
SubsetVariableNumber = "6" # Number of variables in subset
NumProjections = 12 # Number of Maxent models to stack in loop (models numbered consecutively)
DataSetType = "Top250"
TotVars = "55"
GridDirectPre = InDirectPre + DataSetType + "_" + SubsetVariableNumber + "Vars_"
MaxentCore = Raster(OutDirect + Species2 + "mxtcore")

# Create output folder
arcpy.management.CreateFolder(InDirectPre1, OutDirectName)

LoopCountList = range(1,NumProjections+1)
#LoopCountList = [1,2,3,4,5,6,7,8,9,10,11,12]
MaxentRasterCalL = range(0,NumProjections+1) # Create list for Maxent rasters 
Count = 0

for LoopCount in LoopCountList:
    Count = LoopCount
    #Count = Count + 1
    # Read in calibrated raster
    MaxentRasterCal = Raster(GridDirectPre + str(Count) + "/" + Species + "Maxent" + DataSetType + "_" + str(Count) + "ThreshCal_" + SubsetVariableNumber + "Vars_Beta2.tif")
    MaxentRasterCalOut = Int(MaxentRasterCal)
    MaxentRasterCalOut.save(OutDirect + Species2 + SubsetVariableNumber + "of" + TotVars + "_" + str(Count))
    MaxentRasterCalL[LoopCount] = MaxentRasterCalOut

# Create stacked binary calibrated rasters
MaxentRasterMeanFSECon = Int((sum(MaxentRasterCalL)/NumProjections))

# Save raster
MaxentRasterMeanFSECon.save(OutDirect + Species2 + SubsetVariableNumber + "of" + TotVars + "FSEM")

# Mask raster by 100% consensus
MaxentRasterMeanFSEConCore = MaxentRasterMeanFSECon  * MaxentCore

# Save raster
MaxentRasterMeanFSEConCore.save(OutDirect + Species2 + SubsetVariableNumber + "of" + TotVars + "FSMC") 

# Convert above masked FSE raster to polyline shapefile
arcpy.RasterToPolyline_conversion(in_raster=OutDirect + Species2 + SubsetVariableNumber + "of" + TotVars + "FSMC",
                                  out_polyline_features=OutDirect + Species + "_MaxEnt" + Rank + SubsetVariableNumber + "of" + TotVars + "FSEMNCPolyline.shp",
                                  background_value="ZERO", minimum_dangle_length="0", simplify="NO_SIMPLIFY", raster_field="VALUE")

# Clip above polylines to species habitat (not roadkill occurrence) background evaluation extent
arcpy.Clip_analysis(in_features=OutDirect + Species + "_MaxEnt" + Rank + SubsetVariableNumber + "of" + TotVars + "FSEMNCPolyline.shp",
                    clip_features=InDirect1 + HabitatSpecies + "_BackgroundEvaluationExtentAreaNAD83Albers.shp",
                    out_feature_class=OutDirect + Species + "_MaxEnt" + Rank + SubsetVariableNumber + "of" + TotVars + "FSEMNCPolylineC.shp", cluster_tolerance="")


