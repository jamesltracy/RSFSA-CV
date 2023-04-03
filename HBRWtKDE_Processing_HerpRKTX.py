# -*- coding: cp1252 -*-
## This code spatially thins target group and HabitatSpecies data and creates 100000 random points at 10 km spacing
## over the background extent raster
## NOTE: If available, can use ArcMap project template LandCoverRiskAssessment_PointProcessing.mxd

import arcpy
from arcpy.sa import *

import os

# Overwrite pre-existing files
arcpy.env.overwriteOutput = True

################################################  
# set working directory
arcpy.env.workspace = "E:/TexasRoadEnvLayersInt"

## Set input and output directories
Species = "HerpRKTX10_21" 
SpecShort = "herprktx"
CoreRasterName = "hrtxmxtcore"
InDirect = "E:/TexasRoadEnvLayersInt/"
InDirect2 = "E:/MaxEntOutput/HerpRKTX10_21_ModelsAIC/" # Core raster directory
OutDirect = "E:/TexasLayersWk/"
ExtentRaster1 = "2010pop1kalb"
ExtentRaster2 = "devop_500mr"
CellRaster = Raster(InDirect + "devop_500mr") # raster for kde final cell size
MxtCoreRaster = Raster(InDirect2 + CoreRasterName)

################################################################################################
### Use R program "Species"..._HBRWeight.R to output human population density weighted
### occurrence points to use in KDE below
#################################################################################################
## NOTE: Calculate KDE in ArcMap Project with NAD 83 Albers Projection and Extent Raster in Project

# Set processing extent
arcpy.env.extent = InDirect + ExtentRaster1

# Set snap raster
arcpy.env.snapRaster = InDirect + ExtentRaster1

# Calculate KDE at 1 km resolution for study area
arcpy.gp.KernelDensity_sa(InDirect + Species + "_Thinned10km_PopDenHBRWtNAD83Alb.shp","PDNHBRWN", OutDirect + SpecShort + "kde1",
                          InDirect + "2010pop1kalb", "", "SQUARE_KILOMETERS", "EXPECTED_COUNTS", "GEODESIC")

# Set processing extent
arcpy.env.extent = InDirect + ExtentRaster2

# Set snap raster
arcpy.env.snapRaster = InDirect + ExtentRaster2

# Define cell size based on raster here due to bug in Arcpy
CellSize = "{0} {1}".format(arcpy.Describe(CellRaster).meanCellWidth,
                            arcpy.Describe(CellRaster).meanCellHeight)

# Resample to 30 m
arcpy.Resample_management(in_raster=OutDirect + SpecShort + "kde1",
                          out_raster=InDirect + SpecShort + "kde2", cell_size=CellSize, resampling_type="BILINEAR")

# Clip above kde raster to Texas
arcpy.Clip_management(in_raster=InDirect + SpecShort + "kde2", rectangle="",
                      out_raster=InDirect + SpecShort + "kde", in_template_dataset="TexasNAD83Albers",
                      nodata_value="-3.402823e+38", clipping_geometry="ClippingGeometry",
                      maintain_clipping_extent="NO_MAINTAIN_EXTENT")

# Set to null MaxEnt core raster values of zero
arcpy.gp.SetNull_sa(MxtCoreRaster, "1", InDirect + CoreRasterName + "1", "VALUE=0")

# Multiply above MaxEnt core raster with KDE x 10000 to produce integer output grid
arcpy.gp.RasterCalculator_sa('Int(Raster(InDirect + CoreRasterName + "1") * (Raster(InDirect + SpecShort + "kde") * 10000))',
                             InDirect + SpecShort + "kder")


# Convert above rater to polyline shapefile
arcpy.RasterToPolyline_conversion(in_raster=InDirect + SpecShort + "kder", out_polyline_features=InDirect + SpecShort + "kder_polyline.shp",
                                  background_value="ZERO", minimum_dangle_length="0", simplify="NO_SIMPLIFY", raster_field="VALUE")
