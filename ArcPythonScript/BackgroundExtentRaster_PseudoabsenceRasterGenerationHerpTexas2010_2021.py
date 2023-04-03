## For ecological niche models of a given species, this code uses the presence point data and an enviornmental raster to
## create a rasters and polygon shapefile for both the background evaluation extent and the pseudoabsence extent
## (with buffers around presence points)
## Various background extent and pseudoabsence buffer sizes must be specified using guidelines below

# These rasters can then be processed in R and converted to point shapefiles from which random
# background and pseudoabsence points can be generated and spatially thinned

# make arcmap python commands and spatial analyst (sa) operations available
import arcpy
from arcpy.sa import *

# Overwrite pre-existing files
arcpy.env.overwriteOutput = True

import os
###########################################################################################
# Set names of working directory, species occurrence shapefile (unthinned), and various buffer sizes

OutDirect = "E:/TexasRoadEnvLayersInt/" # Name of working and output directory
SpecShort= "HRTX10_21" # Identifying prefix for output files
SpecFileName = "comb_herps_adjalb"  # Name of species point shapefile (not thinned)
RegionFileName = "TexasAlbers" # Name of funnel area polygon shapefile
RoadNetworkName = "HighwaytoUnclassified_Tx3kmbufa" # Name
EnvMskName = "env30m0" # Name of an INTEGER (not floating point) final environmental raster in working directory to use for creating background and extent rasters
BiasRasterName = SpecShort1 + "bbck.tif"

######################
### IMPORTANT NOTE: In Environment Settings for ArcMap project in which running this script, set the Processing Extent
### and Snap Raster to match one of the final environmental rasters
######################

## Designate the pseudoabsence buffer as generally 20 km or 100 km for 1 km resolution data and 2 km for 30 m resolution data
## 100 km pseudoabsence buffer typically may be used for broadly ranging species across large portion of continent
## 20 km pseudoabsence buffer may be used for (1) species with more restricted range, (2) species with data availabe for
## a small potential portion of range (such as for invasive species still spreading), or (3) in the case where true absence data is available.
PseudoabsenceBuffDist = "2" # In Kilometers. Pseudoabsence buffer distance from Presence points.

## Designate buffer of evaluation extent area around convex hull polygon of shape points as generally 1,000 km or 100 km
## 1,000 km extent buffer typically may be used for broadly ranging species across large portion of continent
## 100 km extent buffer may be used for (1) species with more restricted range, or (2) species with data availabe for
## a small potential portion of range (such as for invasive species still spreading)
ExtentBuffer = "10 kilometers"

# Set work directory
arcpy.env.workspace = OutDirect

### Create an environmental raster mask called env1k
EnvRaster = Raster(OutDirect + EnvMskName)
### Set extent environment
arcpy.env.extent = EnvRaster.extent
###
##out = Int(Divide(Plus(Raster(EnvRasterName),1),Plus(Raster(EnvRasterName),1)))
##out.save(OutDirect + "env1k")

# Create convex hull polygon encompassing unthinned presence points
arcpy.MinimumBoundingGeometry_management(in_features=OutDirect + SpecFileName + ".shp",
                                         out_feature_class=OutDirect + SpecShort + "ConvexHull2",
                                         geometry_type="CONVEX_HULL", group_option="ALL",
                                         group_field="", mbg_fields_option="NO_MBG_FIELDS")

# Buffer around above convex hull polygon by ExtentBuffer
arcpy.Buffer_analysis(in_features=OutDirect + SpecShort + "ConvexHull2.shp",
                      out_feature_class=OutDirect + SpecShort + "Extent2",
                      buffer_distance_or_field=ExtentBuffer, line_side="FULL",
                      line_end_type="ROUND", dissolve_option="NONE", dissolve_field="", method="PLANAR")

# Intersect above polygon with Region shapefile to create background evaluation extent area (BEE) polygon
arcpy.Intersect_analysis(in_features=[OutDirect + SpecShort + "Extent2.shp", OutDirect + RegionFileName + ".shp"],
                         out_feature_class=OutDirect + SpecShort + "BackgroundEvaluationExtentArea.shp",
                         join_attributes="ALL", cluster_tolerance="-1 Unknown", output_type="INPUT")

# Project above shapefile to geographic
arcpy.Project_management(in_dataset=OutDirect + SpecShort + "BackgroundEvaluationExtentArea.shp",
                         out_dataset=OutDirect + SpecShort + "BackgroundEvaluationExtentAreaGeo,shp",
                         out_coor_system="GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]",
                         transform_method="WGS_1984_(ITRF00)_To_NAD_1983", in_coor_system="PROJCS['NAD_1983_Albers',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Albers'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-96.0],PARAMETER['Standard_Parallel_1',20.0],PARAMETER['Standard_Parallel_2',60.0],PARAMETER['Latitude_Of_Origin',40.0],UNIT['Meter',1.0]]",
                         preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")


# Use above NAD83 BEE polygon to clip RoadNetwork shapefile to give background polyline shapefile 
arcpy.Clip_analysis(in_features=OutDirect + RoadNetworkName + ".shp", clip_features=OutDirect + SpecShort + "BackgroundEvaluationExtentArea.shp",
                    out_feature_class=OutDirect + SpecShort + "bgrpolyline.shp", cluster_tolerance="")


# Generate random points along above background area line shapefile
arcpy.CreateRandomPoints_management(out_path=OutDirect, out_name=SpecShort + "bgr_randpnts", constraining_feature_class=OutDirect + SpecShort + "BackgroundEvaluationExtentArea.shp",
                                    constraining_extent="",
                                    number_of_points_or_field="100000", minimum_allowed_distance="", create_multipoint_output="POINT", multipoint_size="0")

### Generate background extent raster road mask of value 1
##arcpy.Clip_management(in_raster=EnvRaster + 1,
##                      out_raster=OutDirect + SpecShort + "bck1",
##                      in_template_dataset=OutDirect + SpecShort + "BackgroundEvaluationExtentArea.shp",
##                      nodata_value="0", clipping_geometry="ClippingGeometry",
##                      maintain_clipping_extent="NO_MAINTAIN_EXTENT")

###################################################################################
### In r program TM5CJT.r, created bias raster and mask to above background raster
###################################################################################

########################################################################
#
## Clip biased background raster to roadkill background evaluation extent
arcpy.Clip_management(in_raster=BiasRasterName,
                      out_raster=OutDirect + SpecShort + "bbck",
                      in_template_dataset=OutDirect + SpecShort + "BackgroundEvaluationExtentArea.shp",
                      nodata_value="0", clipping_geometry="ClippingGeometry",
                      maintain_clipping_extent="NO_MAINTAIN_EXTENT")

## Create biased background points from above bias raster masked to background raster created in r
arcpy.CreateSpatiallyBalancedPoints_ga(in_probability_raster=OutDirect + SpecShort + "bbck",
                                       number_output_points="10000", out_feature_class=OutDirect + SpecShort + "biasbckpnts.shp")


###################################################################
### Generate pseudoabsence points

### Convert background extent raster to polygon shapefile
##arcpy.RasterToPolygon_conversion(in_raster=OutDirect + SpecShort + "bck",
##                                 out_polygon_features=OutDirect + SpecShort + "bgrpoly",
##                                 simplify="NO_SIMPLIFY", raster_field="VALUE")

# Create PseudoabsenceBuffDist km buffer shapefile around unthinned presence points (may take about 10 minutes)
arcpy.Buffer_analysis(in_features=OutDirect + SpecFileName + ".shp",
                      out_feature_class=OutDirect + SpecShort + "_" + PseudoabsenceBuffDist + "kmbuf.shp",
                      buffer_distance_or_field=PseudoabsenceBuffDist + " Kilometers",
                      line_side="FULL", line_end_type="ROUND",
                      dissolve_option="ALL",
                      dissolve_field="",
                      method="GEODESIC")

# Use the above buffer shapefile as a "cut out" to remove buffered areas from evaluation area shapefile
arcpy.Erase_analysis(in_features=OutDirect + SpecShort + "BackgroundEvaluationExtentArea.shp", erase_features=SpecShort + "_" + PseudoabsenceBuffDist + "kmbuf.shp",
                     out_feature_class=OutDirect + SpecShort + "PseudoabsenceEvaluationExtentArea.shp", cluster_tolerance="")

# Generate pseudoabsence bias raster by clipping bias background raster with above shapefile
# (Clip took 20 minutes for HRTX10_21)
arcpy.Clip_management(in_raster=OutDirect + SpecShort + "bbck",
                      out_raster=OutDirect + SpecShort + "bpsa", in_template_dataset=OutDirect + SpecShort + "PseudoabsenceEvaluationExtentArea.shp",
                      nodata_value="0", clipping_geometry="ClippingGeometry", maintain_clipping_extent="NO_MAINTAIN_EXTENT")

## Create biased pseudoabsence points from above bias adjusted pseudoabsence raster created above
arcpy.CreateSpatiallyBalancedPoints_ga(in_probability_raster=OutDirect + SpecShort + "bpsa",
                                       number_output_points="10000", out_feature_class=OutDirect + SpecShort + "biaspsapnts.shp")



# Generate random points along above psaarea line shapefile
##arcpy.CreateRandomPoints_management(out_path=OutDirect, out_name=SpecShort + "psa_randpnts", constraining_feature_class=OutDirect + SpecShort + "psaarea.shp",
##                                    constraining_extent="",
##                                    number_of_points_or_field="100000", minimum_allowed_distance="", create_multipoint_output="POINT", multipoint_size="0")
##

### Convert extent with cut out polygon shapefile to raster for further processing in R where converted to points,
### randomly sampled to 20,000 and thinned to 10km
##arcpy.PolygonToRaster_conversion(in_features=OutDirect + SpecShort + "psaarea2.shp", value_field="Value",
##                                 out_rasterdataset=OutDirect + SpecShort + "psr",
##                                 cell_assignment="CELL_CENTER", priority_field="NONE")
##
