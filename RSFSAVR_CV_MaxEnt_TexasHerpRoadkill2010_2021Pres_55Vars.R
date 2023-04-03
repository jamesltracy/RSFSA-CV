##############################################################################
# This R program obtains evaluation statistics for the Envelope Score model
# by dividing the presence/background absence data into thirds and using one third
# for model building and two thirds for testing to obtain three sets of statistics
##############################################################################
#
##############################################################################
# This section loads libraries, sets working directory and loads lat/long file
##############################################################################
# Specify  directory for saved functions
FunctDirect <- "C:/Users/jamesltracy61/Documents/R/win-library/"
setwd(FunctDirect)
# NOTE: If Maxent cannot open temp directory for output due to length, create
# a shorter temp directory and put in .Renviron file such as below (make sure ba
# up copy of .Renviron file and add original additional lines back in)
#write("TMPDIR = 'C:/RTemp'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))
# Read in User Defined Functions
source("Round2_Function.R")
source("MaxentSubset_GridTrainTestEvalAIC_Calib_Function.R")
source("MaxentSubset_GridTrainTestEvalAIC_AUCbgp_Calib_Function.R")
source("MaxentSubset_TrainTestEvalK_Function.R")
source("MaxentMultiSubset_WrapperCV1TrainEvalAIC_Function.R")
source("MaxentMultiSubset_WrapperCV1TestEvalAIC_Function.R")
source("MaxentMultiSubset_WrapperCV1TrainTestEvalAIC_Function.R")
source("MaxentMultiSubset_WrapperCV2TrainTestEvalAIC_Function.R")
source("MaxentProject_Function.R")
source("MaxentProjectRAW_AICc_Function.R")
source("EvalStatVars_Summary_Function.R")
source("HeatMapDataFrame_Plot_Function.R")
source("VariableSubsetsCorrThresh_Build_Function.R")
source("VariableSubsetSizeCorrThresh_Build_Function.R")
source("GamesHowell.Test.Padj_Function.R")
#startSocketServer(port=8889)
#
# Load needed packages of raster, rgdal, dismo, rjava, and maptools (printouts not shown)
# .libPaths()
#install.packages('ENMeval_2.0.3.zip', lib='C:/Users/james.tracy/Documents/R/win-library/3.6',repos = NULL)
#install.packages(c("stingdist"))
library(sp)
library(spThin)
library(sf)
library(raster)
library(terra)
library(dismo)
library(maptools)
library(PresenceAbsence)
library(car)
library(caret)
library(ENMeval)
library(devtools)
library(MCDM)
library(plyr)
library(dplyr)
library(rJava)
library(stringdist) 
#packageVersion("ENMeval")
#??aic.maxent
#help(package="ENMeval")
# Make sure 64-bit Java is installed to match 64-bit R
## Specify prefix for cateogrical variables (or give complete name if only single
## categorical variable)
#CatVarsPrefix <- "ROADS_CAT"
# If no categorical variables being considered for maxent, need to remove CatVarsPrefix variable
## Specify base arguments for Maxent settings for evaluation run
MaxentBaseArgsIn <- c("betamultiplier=2.0", "threshold=false", "product=false", "writebackgroundpredictions=true")
#
## Specify maxent arguments for jackknife cross validation run
MaxentEvalArgs1 <- c("betamultiplier=2.0", "threshold=false", "product=false", "writebackgroundpredictions=true","outputgrids=FALSE", "replicates=3", "replicatetype=crossvalidate", "responsecurves=true", "jackknife=true")
# Specify beta regularization factor for tagging output grids
BetaMult <- 2
## Check if have categorical variable to add to maxent arguments
if(exists("CatVarsPrefix")==TRUE) {
  MaxentCatArg <- paste0("togglelayertype=", CatVarsPrefix)
  MaxentEvalArgs <- c(MaxentEvalArgs1, MaxentCatArg)
  } else {
  MaxentEvalArgs <- MaxentEvalArgs1
}
# if it exists in the active R session using the remove statement below
#remove(CatVarsPrefix)
# Specify  directory for model input data, including environmental layers,
# presence points, background points, and pseudoabsence points
InDirect <- "D:/TexasRoadEnvLayersInt/"
setwd(InDirect)
# Set input presence point projection
crs.in <- 4326 # WGS84
# Define output CRS projection to match environmental variable rasters
# Read ESRI shapefile with desired output projection
TexasNAD83AlbersShp <- st_read(paste0(InDirect, "/TexasNAD83Albers.shp"))
crs.out <- crs(TexasNAD83AlbersShp)  # ESRI NAD 83 Albers
#
# Specify directory for model outputs
Direct2 <- "D:/HerpRKMaxEntOutput/"
dir.create(paste0(Direct2))
# Specify which climatic data set is being used
DataSet <- "All"
# Specify maximum correlation allowed among variables to reduce multicollinearity
CorrThresh <- 0.7 # maximum correlation threshold
SearchLoops <- 20000 # maximum number of loops to use creating potential variable sets meeting correlation threshold
MaxModSets <- 3000 # maximum number of variable sets needed to keep below correlation threshold
NumberModSets <- 250 # maximum number of variable sets on which to calculate evaluation statistics per subset size
# Specify presence data shapefile name for unthinned data in WGS84 geographic projection
PresDataName <- "comb_herps_adj"
# Specify cutoff for presence for nested models
PresCutOff <- 1  # Default is 1 for presence model
# Cutoff label kfold partitions and output directory
if(PresCutOff > 1){
  PresCutOffLabel <- paste0("GTE", PresCutOff)
} else {
  PresCutOffLabel <- ""
}
# Specify number of observations (NObs): = number of kfold divisions for full variable subset evaluation; and
# = number of top subsets and number of random subsets compared in RSFSA evaluations
NObs <- 10  ## Number of replicates for first two exploration loops
FullNObs <- 250 ## Number of replicates for intensive search looop
FullNObs3 <- 12 #Round2(FullNObs*.2,0) # How many top versus bottom ranked of the selected FullNObs subsets to test for performance differences after RSFSA Stage II
ValidationReps <- 3 # Number of randomizations for Wrapper testing and training data in feature selection stage II (RSFSA II)
CVNum <- ValidationReps*2 # Total number of cross validation fold stats (e.g., AUCWrapperCV) to run in RSFSA Stage II, one for ranking and remaining for testing
## Specify location of temporary directory where R stores Maxent output if want
## deleted to avoid accumulation of many maxent output folders when running features selection algorithm
#TempDir <- "C:/Users/james/AppData/Local/Temp/raster/maxent"
TempDir <- ""
## Specify extent and various descriptors for maxent run
Extent <- ""
Year <- "2010_2021"
ShrtYear <- "10_21"
#ShrtYear <- substr(Year, nchar(Year)-2+1, nchar(Year))
PresLimit <- 10000 # Limited number of presence points through random sample after spatial thinning
SpecFileName <- paste0("Herp_RoadkillTX", Year)
Species <- paste0("HerpRKTX", ShrtYear)
SpecShort <- paste0("HRTX", ShrtYear)
SpeciesGen <- paste0("RoadkillTX")
Species2 <- paste0("HRTX", ShrtYear)
#
MaxBackgrndPseudoabsLimit <- 10000 # specify limit for background and pseudoabsence data
# Specify resolution of rasters in square km
ResolutionEnv <- 0.00003  ## in units of square kilometers
Climate <- "Current"
Tag <- ""
LongSpecies <- paste0("HerpRoadkillTexas", Year)
ExtentRaster <- "centralrdmsk" # name of background evaluation extent raster
EvalType <- "psa_corrfilt0.7"
AnalysisType <- paste0("RSFSA")
# Specify Feature Selection Algorithm"
FSAType <- ""
## Specify weights to criteria for ranking variable subsets by the mean Maxent
## permutation importance per variable versus the frequency of appearance of
## variables in top selected variable subsets used in Stage II
PermImpRank <- 0.6
CountRank <- 0.4
# Specify units of buffers
Units <- "km"
# Specify buffers for spatial thinning and pseudoabsence points from presence points
SpatFiltBuff <- 2 # in km
SpatFiltBuffs <- "2km"
# Specify if spatial filtering
SpatFilt <- "Filtered" # there is spatial filtering
#PsAbsBuff <- 2  # in km
#PsAbsBuffs <- "20km"
PsAbsBuff <- 2
PsAbsBuffs <- "2km"
#Conditions <- paste0(PsAbsBuff, "mPsAbsBuff_", SpatFiltBuff, "mSpatFiltBuff")
Conditions <- paste0(SpatFiltBuff, "mSpatFiltBuff")
Conditions2 <- paste("Only Native Data for Model Input")
## Create directory where to be made environmental grids masked to background area extent
## are to be stored
#dir.create(paste0(InDirect, "/", "CentralBckgrndMask"))
#
GridNamesFile <- paste0(SpecShort, "TexasRoad55GridNamesR.csv")
GridNames <- as.matrix(read.csv(GridNamesFile, na.strings=c("", " ", "NA")))
#str(BandNames)
GridNamesL <- as.list(GridNames[,1])
GridNamesL2 <- toupper(GridNamesL)
VariableNames <- data.frame(GridNamesL2, stringsAsFactors=FALSE)
colnames(VariableNames) <- "VarNames"
VariableNames$VarNames <- gsub("_NS", "", VariableNames$VarNames)
#str(VariableNames)
NumGrids <- nrow(VariableNames)  # Total number of variable grids in study
####
## Set VariableNames to Top 20 of 30 Variables
#DropList <- as.list(c("BARE_500MR", "CTI", "FOREST_500MR", "STRMDIST", "STRMLOFLODIST", "STRMMDFLODIST", "TPI500MR", "TRAFFIC_VOL", "WATER_500MR", "WETLND_500MR"))
#VariableNames <- data.frame(VariableNames[ ! VariableNames$VarNames %in% DropList, ], stringsAsFactors=FALSE)
colnames(VariableNames) <- "VarNames"
#
SetType <- 1 # use all 41 variables
ActualVars <- length(na.omit(GridNames[,SetType]))
TotVars <- nrow(VariableNames) # env + spatial variables
ModelName <- paste0("Maxent", DataSet, "_", AnalysisType)
Model <- paste(ModelName, TotVars, "Feature Subset")
###############
### Read in csv file of sets of variables for loop below
#VariableSetTypes <- read.csv(paste0("VariableSets", TotVars, ".csv"), stringsAsFactors=FALSE, na.strings = c("NA", ""))
VariableSetTypes <- GridNames
colnames(VariableSetTypes) <- c("FullSet")
#
VariableNamesSelect <- data.frame(na.omit(toupper(VariableSetTypes[,SetType])), stringsAsFactors=FALSE)
if(SetType < 10) {
  SetNumber <- paste0("0", SetType)
} else {
  SetNumber <- as.character(SetType)
}
SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes)[SetType])
SetNameF <- colnames(VariableSetTypes)[SetType]
colnames(VariableNamesSelect) <- SetName
################
#### IMPORTANT- Specify output names for files and directories
# Identifying part of output rasters
## Specify run id consisting of date
Run <- 1
Set <- Sys.Date()
output1 <- paste0(Species, PresCutOffLabel, ModelName, ActualVars, "Varsof", TotVars, Extent, Tag)
# Identifying part of output directory
# Create directory for storing envelope score related grids
dir.create(paste0(Direct2, "/", output1))
OutDirect <- paste0(Direct2, "/", output1)
output2 <- paste0(Direct2, "/", output1, "/", EvalType)
dir.create(output2)
# Specify proportion of AUC versus AIC in Multi-Object Optimization ranking
# AUCrank is 1 if doing separate AUC and AIC ranks of 1
# If doing joint ranking, then specify weights for AUC and AIC,
# such as AUC is 50% and AUCrank is 0.5
AUCrank <- 1
#
# Create subdirectories according to AUC/AIC ranking scheme
if(AUCrank==1) {                       
  # Make separate directories for AUC and AIC ranked data
  OutDirectAUC <- "AUCRanked"
  dir.create(paste0(output2, "/", OutDirectAUC))
  OutDirectAIC <- "AICRanked"
  dir.create(paste0(output2, "/", OutDirectAIC))
  RankList <- c("AUC", "AIC")
} else {
  # Make one directory specifying joint AUC/AIC ranking
  AICcrank <- 1 - AUCrank
  OutDirectAUC <- paste0("AUC", AUCrank, "_AIC", AICcrank, "_Ranked")
  dir.create(paste0(output2, "/", OutDirectAUC))
  RankList <- c(OutDirectAUC)
}
#
##############################################################################
## Read presence data from wgs84 shapefile and spatially thin to 2 km
##############################################################################
setwd(InDirect)
# Read shapefile to simple features data frame
pres.sfdf <- st_read(paste0(InDirect, "/", PresDataName, ".shp"))
#str(pres.sfdf)
head(pres.sfdf)
nrow(pres.sfdf)
pres.sfdf$Species <- Species
# Read ESRI shapefile with input projection
TexasGeoShp <- st_read(paste0(InDirect, "/TexasGeo.shp"))
# xlim sets the longitude limits and ylim the latitude limits
dev.new()
plot(st_geometry(TexasGeoShp), axes=TRUE, col='light yellow')
plot(st_geometry(pres.sfdf), add=TRUE, col='red')
#
presthin.df <- as.data.frame(thin(pres.sfdf, lat.col="NEAR_Y", long.col="NEAR_X", spec.col = "Species", thin.par=SpatFiltBuff, reps=1, write.file=FALSE, locs.thinned.list.return = TRUE))
head(presthin.df)
colnames(presthin.df) <- c("NEAR_X", "NEAR_Y")
nrow(presthin.df)
#
# Merge data with original roadkill by lat/long
presthin.df <- merge(presthin.df, pres.sfdf, by=c("NEAR_X", "NEAR_Y"))
head(presthin.df)
nrow(presthin.df)
# Remove duplicates
presthin.df2 <- presthin.df[!duplicated(presthin.df[1:2]),]
nrow(presthin.df2)
# 
# Convert back to simple features data frame
presthin.sfdf <- st_as_sf(x = presthin.df, coords = c("NEAR_X", "NEAR_Y"), crs = crs.in)
plot(st_geometry(presthin.sfdf), add=TRUE, col='blue')
# Save thinned data
st_write(presthin.sfdf, paste0(InDirect, PresDataName, "_", SpatFiltBuffs, "thin.shp"), append=FALSE)
#
# Reproject to NAD83 match environmental layers
TexasNAD83AlbersShp <- st_read(paste0(InDirect, "TexasNAD83Albers.shp"))
presthin.sfdf2 <- st_transform(presthin.sfdf, crs=crs.out)
crs(presthin.sfdf2)
# Save reprojected thinned data
st_write(presthin.sfdf2, paste0(InDirect, PresDataName, "_", SpatFiltBuffs, "thinAlb.shp"), append=FALSE)
#####################################################################
### In ArcGIS, use Pythons script BackgroundExtentRaster_PseudoabsenceGenerationHerpTX.py
### to generate road mask rasters for buffered Background Evaluation Extent (10 km beyond 
### convex hull polygon around unthinned presence data) and use as a mask for bias layer created
### in BiasBackgroundRasterCreation.r; then go back to Pythond script and convert
### bias background raster to bias pseudoabsence raster (bias backround 
### raster with holes cut 2 km from unthinned presence data); finally use python 
### script to create 10,000 bias adjusted background and pseudoabsence points from
### above rasters
########################################################################################
##
#
###################################################################
### Extract variable values for presence, absence and background
### poins for running Samples With Data (SWD) evaluations
############################################################################################
### Read in Pseudoabsence and Background points generated in ArcGIS for each species
## NOTE: If shapefiles are open in ArcGIS, a lock may give errors when reading them back in and 
## it may be necessary to close ArcGIS.
setwd(InDirect)
TexasNAD83AlbersShp <- st_read(paste0(InDirect, "TexasNAD83Albers.shp"))
crs(TexasNAD83AlbersShp)
PseudoabsencePntShp <- st_read(paste0(InDirect, SpecShort, "biaspsapnts.shp"))  
BackgroundPntShp <- st_read(paste0(InDirect, SpecShort, "biasbckpnts.shp")) 
## Read back in spatially thinned presence points
PresThinShp <- st_read(paste0(PresDataName, "_", SpatFiltBuffs, "thinAlb.shp")) 
crs(PresThinShp)
##
#str(PresThinShp)
dev.new()
plot(st_geometry(TexasNAD83AlbersShp), axes=TRUE, col='light yellow')
plot(st_geometry(BackgroundPntShp), add=TRUE, col='green')
plot(st_geometry(PseudoabsencePntShp), add=TRUE, col='blue')
plot(st_geometry(PresThinShp), add=TRUE, col='red')
crs(PresThinShp)
crs(BackgroundPntShp)
##############################################################################################
#
SubsetVariableNumber <- nrow(VariableNames)
TotVars <- nrow(VariableNames)
# Read into raster stack the selected subset of environmentalgrids
#GridNamesL <- list("elev", "tpi500mr")
#setwd(MaskGrids)
#GridNamesL[5]
#Predictors <- stack(GridNamesL[1:5], proj4string=CRS.NAAlbersEqualAreaConic)
Predictors <- stack(GridNamesL, proj4string=crs.out)
names(Predictors) <- toupper(names(Predictors))
names(Predictors) <- gsub("_NS", "", names(Predictors))
crs(Predictors[[1]])
#
#plot(Predictors[[1]], add=TRUE)
#plot(Predictors[[1]])
#writeRaster(Predictors[[1]], "test.tif", format = "GTiff", overwrite=TRUE)
#################################################################################
### Extract environmental data for presence points
# Extract Predictors for spatial points data frame of thinned roadkill points
setwd(InDirect)
## Read back in spatially thinned presence points
PresThinShp <- st_read(paste0(PresDataName, "_", SpatFiltBuffs, "thinAlb.shp"))
nrow(PresThinShp)
#
head(PresThinShp)
#
PresenceDat.df <- raster::extract(Predictors, PresThinShp)
head(PresenceDat.df)
nrow(PresenceDat.df)
#
#str(PresThinShp)
head(st_coordinates(PresThinShp))
# Re-associate data with species and x y coordinates
PresenceDat.df2 <- merge(st_coordinates(PresThinShp), PresenceDat.df, by="row.names")
head(PresenceDat.df2)
# Replace first column of row.names for PresenceDat.df with LongSpecies name
colnames(PresenceDat.df2)[1] <- c("Species")
PresenceDat.df2$Species <- LongSpecies
nrow(PresenceDat.df2)
# Omit any rows with NA values
PresenceDat.df <- na.omit(PresenceDat.df2)
any(is.na(PresenceDat.df))
nrow(PresenceDat.df)
# Rename coordinates of X/Y or long/lat as "coords.x1" and "coords.x2" for Maxent SWD runs
colnames(PresenceDat.df)[2:3] <- c("coords.x1", "coords.x2")
head(PresenceDat.df)
# Save extracted data
write.table(PresenceDat.df, sep = ",", col.names=NA, file=paste0(Species, PresCutOffLabel, "_Thinned", SpatFiltBuff, Units,"_presevaldat", DataSet, Tag, "_", TotVars, "Vars.csv"))
#####################
#####################
# Extract Predictors from pseudoabsence points (may take about 11 minutes)
setwd(InDirect)
PseudoabsencePntShp <- st_read(paste0(InDirect, SpecShort, "biaspsapnts.shp"))
head(PseudoabsencePntShp)
#
PseudoabsenceDat.df <- raster::extract(Predictors, PseudoabsencePntShp)
head(PseudoabsenceDat.df)
#
# Re-associate data with species and x y coordinates
PseudoabsenceDat.df2 <- merge(st_coordinates(PseudoabsencePntShp), PseudoabsenceDat.df, by="row.names")
head(PseudoabsenceDat.df2)
# Replace first column of row.names for PseudoabsenceDat.df with LongSpecies name
colnames(PseudoabsenceDat.df2)[1] <- c("Species")
PseudoabsenceDat.df2$Species <- "Pseudoabsence"
nrow(PseudoabsenceDat.df2)
# Omit any rows with NA values
PseudoabsenceDat.df <- na.omit(PseudoabsenceDat.df2)
any(is.na(PseudoabsenceDat.df))
nrow(PseudoabsenceDat.df)
# Rename coordinates of X/Y or long/lat as "coords.x1" and "coords.x2" for Maxent SWD runs
colnames(PseudoabsenceDat.df)[2:3] <- c("coords.x1", "coords.x2")
head(PseudoabsenceDat.df)
# Save extracted data
write.table(PseudoabsenceDat.df, sep = ",", col.names=NA, file=paste0(Species, "_", PsAbsBuff, Units,"PsAbsBuff_pseudoabsevaldat", DataSet, Tag, "_", TotVars, "Vars.csv"))
##########################
# Extract Predictors from background points (may take about 11 minutes)
setwd(InDirect)
BackgroundPntShp <- st_read(paste0(InDirect, SpecShort, "biasbckpnts.shp"))
head(BackgroundPntShp)
#
BackgroundDat.df <- raster::extract(Predictors, BackgroundPntShp)
head(BackgroundDat.df)
#
# Re-associate data with species and x y coordinates
BackgroundDat.df2 <- merge(st_coordinates(BackgroundPntShp), BackgroundDat.df, by="row.names")
head(BackgroundDat.df2)
# Replace first column of row.names for BackgroundDat.df with LongSpecies name
colnames(BackgroundDat.df2)[1] <- c("Species")
BackgroundDat.df2$Species <- "Background"
nrow(BackgroundDat.df2)
# Omit any rows with NA values
BackgroundDat.df <- na.omit(BackgroundDat.df2)
any(is.na(BackgroundDat.df))
nrow(BackgroundDat.df)
# Rename coordinates of X/Y or long/lat as "coords.x1" and "coords.x2" for Maxent SWD runs
colnames(BackgroundDat.df)[2:3] <- c("coords.x1", "coords.x2")
head(BackgroundDat.df)
# Save extracted data
setwd(InDirect)
write.table(BackgroundDat.df, sep = ",", col.names=NA, file=paste0(Species, "_backgrounddat", DataSet, Tag, "_", TotVars, "Vars.csv"))
####
#########################################################################################
# Read back in Samples With Data (SWD) formatted data if needed
setwd(InDirect)
DataSetA <- DataSet
TotVarsA <- TotVars
BackgroundDat.df <- as.data.frame(read.csv(paste0(Species, "_backgrounddat", DataSetA, Tag, "_", TotVarsA, "Vars.csv"), sep=",", row.names=1))
nrow(BackgroundDat.df)
head(BackgroundDat.df)
tail(BackgroundDat.df)
if(nrow(BackgroundDat.df)>MaxBackgrndPseudoabsLimit) {
  BackgroundDat.df <- BackgroundDat.df[1:MaxBackgrndPseudoabsLimit ,]
}
# Delete not used variables
#BackgroundDat.df <- subset(BackgroundDat.df, select=-c(BARE_500MR, CTI, FOREST_500MR, STRMDIST, STRMLOFLODIST, STRMMDFLODIST, TPI500MR, TRAFFIC_VOL, WATER_500MR, WETLND_500MR))
##
PseudoabsenceDat.df <- data.frame(read.csv(paste0(Species, "_", PsAbsBuff, Units,"PsAbsBuff_pseudoabsevaldat", DataSet, Tag, "_", TotVars, "Vars.csv"), sep=",", row.names=1, stringsAsFactors=FALSE))
head(PseudoabsenceDat.df)
tail(PseudoabsenceDat.df)
if(nrow(PseudoabsenceDat.df)>MaxBackgrndPseudoabsLimit) {
  PseudoabsenceDat.df <- PseudoabsenceDat.df[1:MaxBackgrndPseudoabsLimit,]
}
nrow(PseudoabsenceDat.df)
# Delete not used variables
#PseudoabsenceDat.df <- subset(PseudoabsenceDat.df, select=-c(BARE_500MR, CTI, FOREST_500MR, STRMDIST, STRMLOFLODIST, STRMMDFLODIST, TPI500MR, TRAFFIC_VOL, WATER_500MR, WETLND_500MR))

##
PresenceDat.df <- data.frame(read.csv(paste0(Species, PresCutOffLabel, "_Thinned", SpatFiltBuff, Units,"_presevaldat", DataSetA, Tag, "_", TotVarsA, "Vars.csv"), sep=",", row.names=1, stringsAsFactors=FALSE))
head(PresenceDat.df)
nrow(PresenceDat.df)
TotPres <- nrow(PresenceDat.df)
# Delete not used variables
#PresenceDat.df <- subset(PresenceDat.df, select=-c(BARE_500MR, CTI, FOREST_500MR, STRMDIST, STRMLOFLODIST, STRMMDFLODIST, TPI500MR, TRAFFIC_VOL, WATER_500MR, WETLND_500MR))
#####
#######################################################################################
## Create and save CVNum two-fold partition schemes
# NOTE: Only run partition scheme once for set of Presence and Pseudoabsence data
# Create kfold partitions with 2/3 data for wrapper training and 1/3 for wrapper testing
#for(i in 1:CVNum) {
#  #i=1
#  PresWrapperTrainNum <- Round2(0.6667*nrow(PresenceDat.df),0)
#  PresWrapperTestNum <- nrow(PresenceDat.df) - PresWrapperTrainNum
#  kfoldgrppl <- split(1:nrow(PresenceDat.df), sample(rep(1:2, c(PresWrapperTrainNum, PresWrapperTestNum))))
#  kfoldgrpp <- rep(1,nrow(PresenceDat.df))
#  kfoldgrpp[kfoldgrppl[[1]]]=1
#  kfoldgrpp[kfoldgrppl[[2]]]=2
#  length(kfoldgrpp)
#  occurrences <- table(unlist(kfoldgrpp))
#  # Repeat for Pseudoabsence data
#  PsaWrapperTrainNum <- Round2(0.6667*nrow(PseudoabsenceDat.df),0)
#  PsaWrapperTestNum <- nrow(PseudoabsenceDat.df) - PsaWrapperTrainNum
#  kfoldgrpal <- split(1:nrow(PseudoabsenceDat.df), sample(rep(1:2, c(PsaWrapperTrainNum, PsaWrapperTestNum))))
#  kfoldgrpa <- rep(1,nrow(PseudoabsenceDat.df))
#  kfoldgrpa[kfoldgrpal[[1]]]=1
#  kfoldgrpa[kfoldgrpal[[2]]]=2
#  length(kfoldgrpa)
#  occurrences <- table(unlist(kfoldgrpa))
#  #### Save kfold partition scheme for later testing
#  setwd(InDirect)
#  write.table(data.frame(kfoldgrpp), file=paste0(Species, PresCutOffLabel, "_PresenceDatRSFSAKfold_", i, "_.csv"), sep=",")
#  write.table(data.frame(kfoldgrpa), file=paste0(Species, "_PseudoabsenceDatRSFSAKfold_", i, "_.csv"), sep=",")
#}

##
#####################################################################################################
### Run Random Subset Features Selection Algorithm (RFSA) evaluations
####################################################################################
####################################################################################
####################################################################################
## Loop through different subset sizes of 19 variables
## Loop takes 6.4 hours for 500 sets each of two to 19 variables
## On new laptop, takes 2.3 hours for 250 sets of 3,6,8,10,12,15,20,30,40,50,60,70,85 of 86 variables
## On new laptop, takes 30 minutes with monarch roadkill data for  hours for 250 sets of 1, 2, 3, 6, 8, 10, 12, 15, 20 of 29 variables
SubsetSizes <- c(1, 2, 3, 4, 5, 6, 7, 8)
#SubsetSizes <- c(2)
#SubsetSizes <- c(15)
#SubsetSizes <- c(3)
length(SubsetSizes)
OutCount <- 0
## Loop takes 69 min for 29 variables for monarch roadkill on new laptop
####################################################################################
t1 <- Sys.time()
## Set ranking scheme for subsets
for(Rank in RankList) {
  # Rank="AUC"
  if(Rank=="AUC") {
    AUCrank <- 1
    AICcrank <- 0
  } else if (Rank=="AIC") {
    AUCrank <- 0
    AICcrank <- 1
  } else {
  }
  #
  OutDirectpsa <- paste0(Direct2, output1, "/", EvalType, "/", Rank, "Ranked")
  MaxentKeepEvalAllStatsL <- list()
  MaxentKeepEvalAllStatSummL <- list()
  OutCount <- 0
  for(SubsetLoop in SubsetSizes) {
  #SubsetLoop <- 2
    # Generate random subsets of 2 of 19 variables
    OutCount <- OutCount + 1
    VariableNamesIn <- VariableNamesSelect
    colnames(VariableNamesIn) <- "VarNames"
    TotVars <- nrow(VariableNamesIn)
    SubsetVariableNumber <- SubsetLoop
    #NumberModSets <- 250
    TestDataType <- "Wrapper"
    OutDirectIn <- OutDirectpsa
    RunType <- "RndFSAFilt"
    Run <- 1
    #
    InitVars <- nrow(VariableNamesSelect)
    #
    ## Using background values data, find variables with correlation below 0.7
    ## First calculate correlation matrix
    if(OutCount == 1) {
      # Delete first three columns to just keep environmental data
      BackgroundDatCorr.df <- BackgroundDat.df[,-1:-3]
      head(BackgroundDatCorr.df)
      #colnames(BackgroundDatCorr.df)
      # Subset data by chosen variables
      BackgroundDatCorr.df <- BackgroundDatCorr.df[, unlist(VariableNamesSelect)]
      head(BackgroundDatCorr.df)
      ### Run Spearman Rank Correlation and produce a correlation matrix on all variables
      AbsValsCor.mat <- abs(cor(BackgroundDatCorr.df, method="spearman"))
      # Find mean value of correlations in set
      LoAbsValsCor.mat <- (lower.tri(AbsValsCor.mat)*AbsValsCor.mat)
      MeanCor <- mean(LoAbsValsCor.mat)
      # Save correlation matrix
      write.table(AbsValsCor.mat, sep = ",", col.names=NA, file=paste0(Species, "_DataSet", SetName, "_CorrelationBackgroundEnvironVariables.csv"))
      ## Plot correlation matrix
      DataFrameName <- data.frame(AbsValsCor.mat)
      PlotName <- paste0(LongSpecies, "_", SetNameF, "_CorrelationBackgroundEnvironVariables.tiff")
      PlotTitle <- paste0(LongSpecies, ": VariableCorrelations ", SetNameF, " from Background\n                                                    (Mean Correlation of ", Round2(MeanCor,2), ")\n")
      ## Call Heat Map Function
      #HeatMapDataFrame.Plot(DataFrameName, PlotName, PlotTitle, OutDirectIn)
    }
    if(SubsetVariableNumber < 2) {
      VariableSubsets <- VariableNamesIn
    } else {
      ## Find subsets meeting correlation criteria
      CorrelationsIn <- AbsValsCor.mat
      SearchLoopsIn <- SearchLoops
      MaxModSetsIn <- NumberModSets # maximum number of sets needed to keep
      CorrThreshIn <- CorrThresh
      SetNameIn <- SetName
      OutDirectIn <- OutDirectpsa
      #
      system.time(VariableSubsets <- VariableSubsetSizeCorrThresh.Build(CorrelationsIn, SearchLoopsIn, SubsetVariableNumber, MaxModSetsIn, CorrThreshIn, SetNameIn, OutDirectIn))
    }
    tail(VariableSubsets)
    Sets <- nrow(VariableSubsets)
    ### Make Heat Map of above matrix for correlation values
    # Convert correlation matrix to data frame
    ## Use function to round up from .5 from http://stackoverflow.com/questions/12688717/round-up-from-5-in-r
    Round2 <- function(x, n) {
      posneg = sign(x)
      z = abs(x)*10^n
      z = z + 0.5
      z = trunc(z)
      z = z/10^n
      z*posneg
    }
    ###
    #
    #FSAType <- paste0("Wrap")
    SetName2 <- paste0(SetName, InitVars, "CorrFilt", SubsetVariableNumber)
    TestDataType <- "Wrapper"
    OutDirectIn <- OutDirectpsa
    RunType <- paste0("RndFSAFilt_CorrThresh", CorrThresh)
    Run <- 1
    #########################
    # Save row names to use in output
    SubsetSize.mat <- matrix(c("Singlets", "Doublets", "Triplets", "Quartets", "Quintets", "Sextets", "Septets", "Octets", "Nonets",
    "Dectets", "Undectets", "Duodectets","Tredectets", "Quattuordectets", "Quindectets", "Sexdectets", "Septendectets", "Octodectets", "Novemdectets",
    "Vigetets", "Unvigetets", "Duovigetets", "Trevigetets", "Quattuorvigetets", "Quinvigetets", "Sexvigetets", "Septenvigetets", "Octovigetet",
    "Novemvigetets", "Trigetets", "Untrigetets", "Duotrigetets", "Tretrigetets", "Quottuortrigetets", "Quintrigetets",
    "Sextrigetets", "Septentrigetets", "Octotrigetets", "Novemtrigetets", "Quadragetets", "Unquadragetets", "Duoquadragetets", "Trequadragetets",
    "Quattuorquadragetets", "Quinquadragetets", "Sexquadragetets", "Octoquadragetets", "Octoquadragetets", "Novemquadragetets", "Quinquagetets",
    "Unquinquagetets", "Duoquinquagetets", "Trequinguagetets", "Quattuorquinquagetets", "Quinquinquagetets",
    "Sexquinquagetets", "Septenquinquagetets", "Octoquinquagetets", "Novemquinquagetets", "Sexagetets"), ncol=1, nrow=60, byrow=TRUE, dimnames=list(c
     (seq(1:60)), c("Subset")))
    SubsetSize.df <- as.data.frame(SubsetSize.mat, stringsAsFactors=FALSE)
    Subset <- SubsetSize.df[SubsetVariableNumber,]
    if(is.na(Subset)) {
      Subset <- ""
    }
    TotVars <- nrow(VariableNamesIn)
    #TotVars <- 260
    #
    #NumberModSets <- 500
    #FSAType <- "NoWrapper"
    output2 <- paste0(NumberModSets,"Setsof", Subset, SubsetVariableNumber, "of", TotVars, "Vars")
    OutDirect2 <- paste0(OutDirectpsa, "/", output2)
    dir.create(OutDirect2)
    setwd(OutDirect2)
    OutDirectIn <- OutDirect2
    #
    ## Evaluate random subsets for NumberModSets number of subsets
    # Takes 4.48 minutes for 1000 runs
    #ncol(combn(19,5))
    #  t1a <- Sys.time()
    ###########################
    # Read in k-fold partition scheme
    setwd(InDirect)
    CVScheme <- 1
    kfoldgrpp <- as.vector(unlist(read.csv(paste0(Species, "_PresenceDatRSFSAKfold_", CVScheme, "_.csv"), stringsAsFactors=FALSE)))
    length(kfoldgrpp)
    #
    kfoldgrpa <- as.vector(unlist(read.csv(paste0(Species, "_PseudoabsenceDatRSFSAKfold_", CVScheme, "_.csv"), stringsAsFactors=FALSE)))
    length(kfoldgrpa)
    ####################
    # Designate Wrapper Training data for Wrapper test
    MaxentPresTrainData <- PresenceDat.df[kfoldgrpp==1,]
    head(MaxentPresTrainData)
    nrow(MaxentPresTrainData)
    MaxentAbsTrainData <- BackgroundDat.df
    nrow(MaxentAbsTrainData)
    head(MaxentAbsTrainData)
    WrapperTrainSWD <- rbind(MaxentPresTrainData[,4:ncol(MaxentPresTrainData)], MaxentAbsTrainData[,4:ncol(MaxentAbsTrainData)])
    head(WrapperTrainSWD)
    tail(WrapperTrainSWD)
    WrapperTrainPresID <- data.frame(rep(1,nrow(MaxentPresTrainData)))
    colnames(WrapperTrainPresID) <- "ID"
    TotPres <- nrow(WrapperTrainPresID)
    WrapperTrainAbsID <- data.frame(rep(0,nrow(MaxentAbsTrainData)))
    colnames(WrapperTrainAbsID) <- "ID"
    WrapperTrainPresAbsID <- rbind(WrapperTrainPresID, WrapperTrainAbsID)
    head(WrapperTrainPresAbsID)
    tail(WrapperTrainPresAbsID)
    VariableSubsetsIn <- VariableSubsets
    SetRunID <- ""
    #
   if(Rank!="AIC") {
      system.time(MaxentFiltSubsetEvalStats.df <- MaxentMultiSubset.WrapperCV1TestEvalAIC(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsetsIn, SetRunID, WrapperTrainSWD, WrapperTrainPresAbsID, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, OutDirectIn, MaxentArgsIn=MaxentBaseArgsIn))
      setwd(OutDirect2)
    } else {
      # If necessary, read back in evaluation stats
      ## NOTE: The next two lines are for running the loop with reranking and testing only
      OutDirect2Orig <- gsub("AIC", "AUC", OutDirect2)
      setwd(OutDirect2Orig)
    }
    ######################
    ## Rank top NObs subsets by Wrapper TSS and obtain training and final testing TSS for evaluation of each set
    ## and for ensemble set
    ##
    Run <- 1
    Sets <- min(Sets, NumberModSets) # Sets must kept the same to read output from line 653 defining MaxentFiltSubsetEvalStats.df
    TestDataType <- "WrapperCV1"
    MaxentFiltSubsetEvalStats.df <- data.frame(read.csv(paste0(Species, "MaxentResults_WrapperCV1Test_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_All", Sets, "_ForRanking", SetRunID, ".csv"), sep=",", row.names=1, stringsAsFactors=FALSE))
    head(MaxentFiltSubsetEvalStats.df)
    tail(MaxentFiltSubsetEvalStats.df)
    nrow(MaxentFiltSubsetEvalStats.df)
    # Delete any subsets with "Inf" for AICc_bg (means too many parameters in model relative to observations)
    MaxentFiltSubsetEvalStats.df2 <- subset(MaxentFiltSubsetEvalStats.df, AICc_bg!="Inf")
    nrow(MaxentFiltSubsetEvalStats.df2)
    # Delete any subsets with negative AICc_bg which appear to be spurious values
    MaxentFiltSubsetEvalStats.df2 <- subset(MaxentFiltSubsetEvalStats.df2, AICc_bg>-1)
    MaxentFiltSubsetEvalStats.df2 <- na.omit(MaxentFiltSubsetEvalStats.df2)
    head(MaxentFiltSubsetEvalStats.df2)
    nrow(MaxentFiltSubsetEvalStats.df2)
    # Delete models with duplicate EnvVarsUsed
    MaxentFiltSubsetEvalStats.df3 <-MaxentFiltSubsetEvalStats.df2[!duplicated(MaxentFiltSubsetEvalStats.df2$EnvVarsUsed), ]
    head(MaxentFiltSubsetEvalStats.df3)
    nrow(MaxentFiltSubsetEvalStats.df3)
    NewSets1 <- nrow(MaxentFiltSubsetEvalStats.df3)
    ## Use Multi-Objective Optimization from MCDM package to sort by two criteria of AUC and AICc_bg
    decision.mat <- as.matrix(MaxentFiltSubsetEvalStats.df3[,c(11,7)])
    head(decision.mat)
    # Normalize AICc_bg and AUC from zero to 1
    normalizef <- function(x) c((x+(max(x)-min(x))-max(x))/(max(x)-min(x)))
    decision.matn <- as.matrix(apply(decision.mat,2,normalizef))
    head(decision.matn)
    # Assign weights to criteria
    weightscrit <- c(AICcrank, AUCrank)
    # Assign whether cost "min", or benefit "max"
    cb <- c("min", "max")
    # Rank criteria
    MMOORA.Rank <- MMOORA(decision.matn, weightscrit, cb)
    head(MMOORA.Rank)
    MaxentFiltSubsetEvalStats.df4 <- MaxentFiltSubsetEvalStats.df3
    # Join MMOORA rank to original data
    MaxentFiltSubsetEvalStats.df4$MMOORA_Rank <- MMOORA.Rank[,8]
    head(MaxentFiltSubsetEvalStats.df4)
    # Sort data by MMOORA rank
    MaxentFiltSubsetEvalStats.df4s <- MaxentFiltSubsetEvalStats.df4[order(MaxentFiltSubsetEvalStats.df4$MMOORA_Rank),]
    head(MaxentFiltSubsetEvalStats.df4s)
    ncol(MaxentFiltSubsetEvalStats.df4s)
    # Delete rank column
    MaxentFiltSubsetEvalStats.df5 <- MaxentFiltSubsetEvalStats.df4s[,-14]
    head(MaxentFiltSubsetEvalStats.df5)
    # Keep top NObs selected subsets
    MaxentTopNObsFiltSubsetEvalStats.df <- data.frame(MaxentFiltSubsetEvalStats.df5[1:NObs,], stringsAsFactors=FALSE)
    MaxentTopNObsFiltSubsetEvalStats.df
    TopNObsSubsets <- data.frame(MaxentTopNObsFiltSubsetEvalStats.df[,1], stringsAsFactors=FALSE)
    ##############################
    # Keep random NObs random subsets
    if(Rank=="AUC") {
      MaxentRandNObsFiltSubsetEvalStats.df <- data.frame(MaxentFiltSubsetEvalStats.df3[sample(nrow(MaxentFiltSubsetEvalStats.df3), NObs), ], stringsAsFactors=FALSE)
      RandNObsSubsets <- data.frame(MaxentRandNObsFiltSubsetEvalStats.df[,1], stringsAsFactors=FALSE)
      colnames(RandNObsSubsets) <- "VarNames"
    }
    # Use difference between evaluation statistics between train and test data to evaluate overfitting
    ###############################################################################
    ## Test Top NObs Ranking Subsets and Random NObs Subsets
    ###############################################################################
    #####################################
    setwd(OutDirectpsa)
    if(Rank=="AUC") {
      SubsetRuns <- list(TopNObsSubsets, RandNObsSubsets)
      DataSetTypes <- c(paste0("Top",NObs), paste0("Random", NObs))
    } else {
      SubsetRuns <- list(TopNObsSubsets)
      DataSetTypes <- c(paste0("Top",NObs))
    }
    MaxentKeepEvalStatsL1 <- list()
    #
    for (j in 1:length(SubsetRuns)) {
      #N3Subsets <- SubsetRuns[[1]]
      #j=1
      VariableSubsets <- data.frame(SubsetRuns[[j]], stringsAsFactors=FALSE)
      DataSetType <- DataSetTypes[j]
      # Identify Wrapper data evaluation statistics from above in loop
      if(DataSetType==paste0("Top", NObs)) {
        MaxentTestSubsetEvalStats.df <- MaxentTopNObsFiltSubsetEvalStats.df
      } else {
        MaxentTestSubsetEvalStats.df <- MaxentRandNObsFiltSubsetEvalStats.df
      }
      MaxentTestSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperCV1Test")
      ##
      ## Evaluate random subsets for NumberModSets number of subsets
      system.time(MaxentTrainSubsetEvalStats.df <- MaxentMultiSubset.WrapperCV1TrainEvalAIC(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsets, WrapperTrainSWD, WrapperTrainPresAbsID, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, OutDirectIn, MaxentArgsIn=MaxentBaseArgsIn, DataSetType=DataSetType))
      ## Separate out training and testing data
      MaxentTrainSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperCV1Train")
      # Calculate difference between test and train statistics for overfitting
      MaxentDiffTestTrainSubsetEvalStats.df <- MaxentTestSubsetEvalStats.df
      head(MaxentDiffTestTrainSubsetEvalStats.df)
      MaxentDiffTestTrainSubsetEvalStats.df[,5:13] <- MaxentTrainSubsetEvalStats.df[,5:13] - MaxentTestSubsetEvalStats.df[,5:13]
      MaxentDiffTestTrainSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperCV1DiffTestTrain")
      #
      MaxentSubsetEvalStats.df <- do.call(rbind, list(MaxentTrainSubsetEvalStats.df, MaxentTestSubsetEvalStats.df, MaxentDiffTestTrainSubsetEvalStats.df))
      MaxentSubsetEvalStats.df$SetName <- SetName
      MaxentSubsetEvalStats.df$Set <- DataSetType
      MaxentSubsetEvalStats.df$TotSubsets <- Sets
      MaxentSubsetEvalStats.df$PresTrainPoints <- nrow(PresenceDat.df[kfoldgrpp==1,])
      MaxentSubsetEvalStats.df$Run <- Run
      MaxentSubsetEvalStats.df <- MaxentSubsetEvalStats.df[, colnames(MaxentSubsetEvalStats.df)[c(1:4,14:18,5:13)]]
      #
      MaxentKeepEvalStatsL1[[j]] <- MaxentSubsetEvalStats.df
    }
    ###############
    ###########################
    # Read in k-fold partition scheme for wrapper cv
    setwd(InDirect)
    CVScheme <- 2
    kfoldgrpp <- as.vector(unlist(read.csv(paste0(Species, "_PresenceDatRSFSAKfold_", CVScheme, "_.csv"), stringsAsFactors=FALSE)))
    length(kfoldgrpp)
    #
    kfoldgrpa <- as.vector(unlist(read.csv(paste0(Species, "_PseudoabsenceDatRSFSAKfold_", CVScheme, "_.csv"), stringsAsFactors=FALSE)))
    length(kfoldgrpa)
    ####################
    ## Designate Wrapper cv Training data for Wrapper cv Test
    MaxentPresTrainData <- PresenceDat.df[kfoldgrpp==1,]
    head(MaxentPresTrainData)
    nrow(MaxentPresTrainData)
    MaxentAbsTrainData <- BackgroundDat.df
    nrow(MaxentAbsTrainData)
    head(MaxentAbsTrainData)
    WrapperCVTrainSWD <- rbind(MaxentPresTrainData[,4:ncol(MaxentPresTrainData)], MaxentAbsTrainData[,4:ncol(MaxentAbsTrainData)])
    head(WrapperCVTrainSWD)
    tail(WrapperCVTrainSWD)
    WrapperCVTrainPresID <- data.frame(rep(1,nrow(MaxentPresTrainData)))
    colnames(WrapperCVTrainPresID) <- "ID"
    WrapperCVTrainAbsID <- data.frame(rep(0,nrow(MaxentAbsTrainData)))
    colnames(WrapperCVTrainAbsID) <- "ID"
    WrapperCVTrainPresAbsID <- rbind(WrapperCVTrainPresID, WrapperCVTrainAbsID)
    head(WrapperCVTrainPresAbsID)
    tail(WrapperCVTrainPresAbsID)
    #
    setwd(OutDirectpsa)
    if(Rank=="AUC") {
      SubsetRuns <- list(TopNObsSubsets, RandNObsSubsets)
      DataSetTypes <- c(paste0("Top",NObs), paste0("Random", NObs))
    } else {
      SubsetRuns <- list(TopNObsSubsets)
      DataSetTypes <- c(paste0("Top",NObs))
    }
    MaxentKeepEvalStatsL2 <- list()
    #
    for (j in 1:length(SubsetRuns)) {
      #N3Subsets <- SubsetRuns[[1]]
      #j=1
      VariableSubsets <- data.frame(SubsetRuns[[j]], stringsAsFactors=FALSE)
      DataSetType <- DataSetTypes[j]
      ## Evaluate random subsets for NumberModSets number of subsets
      system.time(MaxentTrainTestSubsetEvalStats.df <- MaxentMultiSubset.WrapperCV2TrainTestEvalAIC(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsets, WrapperCVTrainSWD, WrapperCVTrainPresAbsID, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, OutDirectIn, MaxentArgsIn=MaxentBaseArgsIn, DataSetType=DataSetType))
      ## Separate out training and testing data
      MaxentTrainSubsetEvalStats.df <- MaxentTrainTestSubsetEvalStats.df[which(MaxentTrainTestSubsetEvalStats.df$DataType=="WrapperCV2Train"),]
      MaxentTrainSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperCV2Train")
      MaxentTestSubsetEvalStats.df <- MaxentTrainTestSubsetEvalStats.df[which(MaxentTrainTestSubsetEvalStats.df$DataType=="WrapperCV2Test"),]
      MaxentTestSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperCV2Test")
      # Calculate difference between test and train statistics for overfitting
      MaxentDiffTestTrainSubsetEvalStats.df <- MaxentTestSubsetEvalStats.df
      MaxentDiffTestTrainSubsetEvalStats.df[,5:13] <- MaxentTrainSubsetEvalStats.df[,5:13] - MaxentTestSubsetEvalStats.df[,5:13]
      MaxentDiffTestTrainSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperCV2DiffTestTrain")
      #
      MaxentSubsetEvalStats.df <- do.call(rbind, list(MaxentTrainSubsetEvalStats.df, MaxentTestSubsetEvalStats.df, MaxentDiffTestTrainSubsetEvalStats.df))
      MaxentSubsetEvalStats.df$SetName <- SetName
      MaxentSubsetEvalStats.df$Set <- DataSetType
      MaxentSubsetEvalStats.df$TotSubsets <- Sets
      MaxentSubsetEvalStats.df$PresTrainPoints <- nrow(PresenceDat.df[kfoldgrpp==1,])
      MaxentSubsetEvalStats.df$Run <- Run
      MaxentSubsetEvalStats.df <- MaxentSubsetEvalStats.df[, colnames(MaxentSubsetEvalStats.df)[c(1:4,14:18,5:13)]]
      #
      MaxentKeepEvalStatsL2[[j]] <- MaxentSubsetEvalStats.df
    }
    ###############
    MaxentKeepEvalStats.df1 <- do.call(rbind, MaxentKeepEvalStatsL1)
    MaxentKeepEvalStats.df2 <- do.call(rbind, MaxentKeepEvalStatsL2)
    MaxentKeepEvalStats.df3 <- rbind(MaxentKeepEvalStats.df1, MaxentKeepEvalStats.df2)
    ## There will be duplicate EnvVarsUsed in both Random and Top data for CV2 that did not appear in CV1,
    ## but the duplicates will be left in order not to bias the results of Random or Top selected subsets
    ## producing more duplicates when variables are lost for potentially poorer less generalizable models
    ## selected in CV1
    ##
    # Sort by DataType
    MaxentKeepEvalStats.df4 <- MaxentKeepEvalStats.df3[with(MaxentKeepEvalStats.df3, order(MaxentKeepEvalStats.df3$DataType)), ]
    nrow(MaxentKeepEvalStats.df4)
    #
    MaxentKeepEvalAllStatsL[[OutCount]] <- MaxentKeepEvalStats.df4
    head(MaxentKeepEvalStats.df4)
    tail(MaxentKeepEvalStats.df4)
    ### Summarize Output
    ###Calculate Mean and Standard Deviation Values
    RunType <- "WrapperCV2"
    EvaluationStats <- MaxentKeepEvalStats.df4
    ncol(EvaluationStats)
    EvaluationStats$Model <- paste0(ModelName, Rank)
    EvaluationStats <- EvaluationStats[, colnames(EvaluationStats)[c(19,1:18)]]
    head(EvaluationStats)
    nrow(EvaluationStats)
    # Delete any subsets with "Inf" for AICc_bg (means too many parameters in model relative to observations)
    EvaluationStats <- subset(EvaluationStats, AICc_bg!="Inf")
    # Replace NAN values with zero
    EvaluationStats[is.na(EvaluationStats)] <- 0
    ##
    #FSAType <- RankStatistic
    OutName <- paste0(Species, ModelName, Rank, SetNameF, "SummaryStats_", SubsetVariableNumber, "of", TotVars, "Vars", FSAType,"_", RunType, ".csv")
    SortGroups <- c("Model", "DataType", "Set", "SetName", "SubsetVariableNumber", "TotSubsets")
    ## All statistics, and only statistics, should be at and after the column specified below
    StatVarFirstColumn <- 9
    OutDirectIn <- OutDirectpsa
    EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
    head(EvaluationStatsSummary.df)
    MaxentKeepEvalAllStatSummL[[OutCount]] <- EvaluationStatsSummary.df
    ######
    MaxentKeepEvalStatsAllSets.df <- do.call(rbind, MaxentKeepEvalAllStatSummL)
    ncol(MaxentKeepEvalStatsAllSets.df)
    ## Save data
    write.table(MaxentKeepEvalStatsAllSets.df, file=paste0(Species, ModelName, Rank, "_", SetNameF, "SummaryStats_of", min(SubsetSizes), "to", max(SubsetSizes), "Vars", FSAType,"_", RunType, ".csv"), sep=",", col.names=NA)
    ### Summarize and save raw data as loop progresses
    MaxentKeepEvalStats.df <- do.call(rbind, MaxentKeepEvalAllStatsL)
    head(MaxentKeepEvalStats.df)
    # Save data
    setwd(OutDirectpsa)
    write.table(MaxentKeepEvalStats.df , file=paste0(Species, ModelName, Rank, NObs, "_", SetNameF, "TrainVsTestStats_", min(SubsetSizes), "to", max(SubsetSizes), "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), sep=",", col.names=NA)
  }
}
#
####################################
t2 <- Sys.time()
##################
difftime(t2,t1, units = "mins")

############################################
## Process and analyze data from above run
###########################################
for(Rank in RankList) {
  #Rank <- "AIC"
  OutDirectpsa <- paste0(Direct2, output1, "/", EvalType, "/", Rank, "Ranked")
  setwd(OutDirectpsa)
  RunType <- "WrapperCV2"
  # Read data back in and sort
  MaxentKeepEvalSumStats.df1 <- data.frame(read.csv(paste0(Species, ModelName, Rank, "_", SetNameF, "SummaryStats_of", min(SubsetSizes), "to", max(SubsetSizes), "Vars", FSAType,"_", RunType, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
  head(MaxentKeepEvalSumStats.df1)
  tail(MaxentKeepEvalSumStats.df1)
  #MaxentKeepEvalSumStats.df1 <- MaxentKeepEvalSumStats.df1[grepl("Top", MaxentKeepEvalSumStats.df1$DataType), ]
  # Read in results for Random data from "AUC" run if not "AUC" rank
  if(Rank!="AUC") {
    OutDirectpsa_AUC <- gsub(Rank, "AUC", OutDirectpsa)
    setwd(OutDirectpsa_AUC)
    # Read data back in 
    MaxentKeepEvalSumStats.df1a <- data.frame(read.csv(paste0(Species, ModelName, "AUC", "_", SetNameF, "SummaryStats_of", min(SubsetSizes), "to", max(SubsetSizes), "Vars", FSAType,"_", RunType, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
    tail(MaxentKeepEvalSumStats.df1a)
    # Keep only random data
    MaxentKeepEvalSumStats.df1b <- MaxentKeepEvalSumStats.df1a[grepl("Random", MaxentKeepEvalSumStats.df1a$DataType), ]
    tail(MaxentKeepEvalSumStats.df1b)
    # Merge with ranked data
    MaxentKeepEvalSumStats.df1 <- rbind(MaxentKeepEvalSumStats.df1, MaxentKeepEvalSumStats.df1b)
  }
  setwd(OutDirectpsa)
  ## Sort data by Statistic and then by DataType for ease of plotting in excel
  MaxentKeepEvalSumStats.dfs1 <- arrange(MaxentKeepEvalSumStats.df1, Statistic, DataType)
  MaxentKeepEvalSumStats.dfs <- arrange(MaxentKeepEvalSumStats.dfs1 , DataType)
  head(MaxentKeepEvalSumStats.dfs)
  nrow(MaxentKeepEvalSumStats.dfs)
  #
  write.table(MaxentKeepEvalSumStats.dfs, file=paste0(Species, ModelName, Rank, "_", SetNameF, "SummaryStats_of", min(SubsetSizes), "to", max(SubsetSizes), "Vars", FSAType,"_", RunType, ".csv"), sep=",", col.names=NA)
  #
  #########################################################################################
  ## Run Welch's t-test to test whether Random models differ from Top models
  ## Null hypothesis is that there is no difference between Random and Top models for
  ## each individual subset size, so Bonferroni correction not needed.
  #######################################################################################
  setwd(OutDirectpsa)
  # Read in results for first FSA run to get output on the leave one out variable run with all variables
  #FSAType <- RankStatistic
  #SubsetSizes <- c(1, 2, 3, 4, 5, 6, 7, 8)
  #
  MaxentKeepEvalStats.df1 <- data.frame(read.csv(paste0(Species, ModelName, Rank, NObs, "_", SetNameF, "TrainVsTestStats_", min(SubsetSizes), "to", max(SubsetSizes), "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
  nrow(MaxentKeepEvalStats.df1)
  ncol(MaxentKeepEvalStats.df1)
  head(MaxentKeepEvalStats.df1)
  tail(MaxentKeepEvalStats.df1)
  if(Rank!="AUC") {
    OutDirectpsa_AUC <- gsub(Rank, "AUC", OutDirectpsa)
    setwd(OutDirectpsa_AUC)
    # Read data back in
    MaxentKeepEvalStats.df1a <- data.frame(read.csv(paste0(Species, ModelName, "AUC", NObs, "_", SetNameF, "TrainVsTestStats_", min(SubsetSizes), "to", max(SubsetSizes), "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
    tail(MaxentKeepEvalStats.df1a)
    # Keep only random data
    MaxentKeepEvalStats.df1b <- MaxentKeepEvalStats.df1a[grepl("Random", MaxentKeepEvalStats.df1a$DataType), ]
    tail(MaxentKeepEvalStats.df1b)
    # Merge with ranked data
    MaxentKeepEvalStats.df1 <- rbind(MaxentKeepEvalStats.df1, MaxentKeepEvalStats.df1b)
  }
  #####################################
  # Create dataframe store p values
  TResultsAll.df <- data.frame("StatType" = "AUC", "SubsetVariableNumber" = 1, "P-value" = 0.0001, stringsAsFactors=FALSE)
  #str(TResultsAll.df)
  #
  StatTypeNames <- c("AUC", "AICc_bg", "AUCdiff")
  Count <- 1
  for(StatType in StatTypeNames) {
    # StatType <- "AICc_bg"
    # Select test results according to StatType
    if(StatType!="AUCdiff") {
      MaxentKeepEvalStats.df3 <- MaxentKeepEvalStats.df1[grepl("CV2Test", MaxentKeepEvalStats.df1$DataType), ]
    } else {
      MaxentKeepEvalStats.df3 <- MaxentKeepEvalStats.df1[grepl("CV2DiffTestTrain", MaxentKeepEvalStats.df1$DataType), ]
    }
    # Delete any subsets with "Inf" for AICc_bg (means too many parameters in model relative to observations)
    if(StatType=="AICc_bg") {
      MaxentKeepEvalStats.df3 <- subset(MaxentKeepEvalStats.df3, AICc_bg!="Inf")
    }
    for(SubsetLoop in SubsetSizes) {
      # SubsetLoop <- 4
      # Count <- 1
      # Select only records for SubsetSize
      MaxentKeepEvalStats.df4 <- MaxentKeepEvalStats.df3[ which(MaxentKeepEvalStats.df3$SubsetVariableNumber==SubsetLoop), ]
      # Separate random and top (selected) model results
      TopEvalStats.df <- MaxentKeepEvalStats.df4[grepl("Top", MaxentKeepEvalStats.df4$Set), ]
      RandEvalStats.df <- MaxentKeepEvalStats.df4[grepl("Rand", MaxentKeepEvalStats.df4$Set), ]
      #  Run Welch t test
      if(StatType=="AUCdiff") {
        StatType1="AUC"
      } else {
        StatType1=StatType
      }
      TResults <- t.test(TopEvalStats.df[[StatType1]], RandEvalStats.df[[StatType1]], var.equal=FALSE)
      # Save above Welch t test output
      out1 <- paste0("Welch t test for Top vs Random model ", StatType, " for ",NObs, " of ", NumberModSets, " models with ", SubsetLoop, " of ", TotVars, " Variables")
      cat(out1, file=paste0(Species, "WelchttestforTopVsRandom", StatType, NObs, "from", NumberModSets, "Setsof", SubsetLoop, "Variables.txt"), sep="\n", append=TRUE)
      out2 <- capture.output(print(TResults))
      cat(out2, file=paste0(Species, "WelchttestforTopVsRandom", StatType, NObs, "from", NumberModSets, "Setsof", SubsetLoop, "Variables.txt"), sep="\n", append=TRUE)
      # Store results in list
      TResultsAll.df[Count,] <- c(StatType, SubsetLoop, TResults$p.value)
      Count <- Count + 1
    }
  }
  #####################
  # Save results
  write.table(TResultsAll.df, sep = ",", col.names=NA, file=paste0(Species, "WelchttestforTopVsRandom", Rank, "Ranked", NObs, "from", NumberModSets, "Setsof", min(SubsetSizes),"to",max(SubsetSizes),"Variables.csv"))
}

#####################################################################################################
### Run Random Subset Features Selection Algorithm (RFSA) evaluations for 10,000 subsets of chosen subset size
####################################################################################
## Loop takes 2.6 hours for 10,200 NumberModsets of 8 variables for flycatcher on new laptop
## Loop takes 2.8 hours for 10,200 NumberModsets of 10 variables for flycatcher on new laptop
## Loop takes 3.75 hours for 10,200 NumberModsets of 15 variables for flycatcher on new laptop
## Loop takes 4.2 hours for 10,150 NumberModsets of 15 variables for flycatcher on new laptop
## Loop takes 3.6 hours for 9,000 NumberModsets of 12 variables for flycatcher on new laptop
## Loop takes 3.07 hours for 9,000 NumberModSets of 6 variables for lo brn fire on new laptop
## Loop takes 10.6 hours for 9,000 NumberModSets of 25 variables for lo brn fire on new laptop
## Loop takes 3.55 hours for 9,000 NumberModSets of 8 variables for zizotes milkweed on new laptop
## Loop takes 2.7 hours for 9,000 NumberModSets of 6 variables for Central Funnel Roadkill on new laptop
###############################################################################
t1 <- Sys.time()
for(Rank in RankList) {
  #Rank="AIC"
  # Rank="AUC"
  if(Rank=="AUC") {
    AUCrank <- 1
    AICcrank <- 0
  } else if (Rank=="AIC") {
    AUCrank <- 0
    AICcrank <- 1
  } else {
  }
  #
  OutDirectpsa <- paste0(Direct2, output1, "/", EvalType, "/", Rank, "Ranked")
  setwd(OutDirectpsa)
  OutCount <- 0
  ##
  MaxentKeepEvalAllStatsL <- list()
  MaxentKeepEvalAllStatSummL <- list()
  TopVariablePermImpL <- list()
  TopMaxentSubsetPermImpL <- list()
  NumberModSetsIn <- 3150 # Number sets per rep: extra sets if needed due to AICc values deleted- recommend no more than 5050 sets at once in case of crash
  NumberModSets <- 3000  # Number sets per rep
  #NumberModSetsIn <- 33 # extra sets if needed due to AICc values deleted
  #NumberModSets <- 30
  SetType.Vec <- c(1)
  #SetType.Vec <- c(1)
  #SetRunsList <- c(1)
  #NumberModSetsIn <- 5200 # extra sets if needed due to AICc values deleted
  #NumberModSets <- 5000
  #
  FinalModelVariables <- 6
  #FinalModelVariables <- 3
  #choose(90,6)
  #choose(8,4)
  SubsetVariableNumber <- FinalModelVariables
  RepList <- seq(1:ValidationReps)
  #FullNObs <- 250
  #PlusNum <- 100 # extra random subsets
  #FullNObs <- NumberModSets/3
  #PlusNum <- 1
  #########################################################
  for(j in SetType.Vec) {
  #for(j in 1:2) {
    #j=1
    VariableNamesIn <- VariableNamesSelect
    colnames(VariableNamesIn) <- "VarNames"
    TotVars <- nrow(VariableNamesIn)
    #
    VariableNamesSelect <- data.frame(na.omit(toupper(VariableSetTypes[,SetType])), stringsAsFactors=FALSE)
    if(SetType < 10) {
      SetNumber <- paste0("0", SetType)
    } else {
      SetNumber <- as.character(SetType)
    }
    SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes)[SetType])
    SetNameF <- colnames(VariableSetTypes)[SetType]
    colnames(VariableNamesSelect) <- SetName
    InitVars <- nrow(VariableNamesSelect)
    #
    ## Using background values data, find variables with correlation below 0.7
    # Delete first three columns to just keep environmental data
    BackgroundDatCorr.df <- BackgroundDat.df[,-1:-3]
    head(BackgroundDatCorr.df)
    # Subset data by chosen variables
    BackgroundDatCorr.df <- BackgroundDatCorr.df[, unlist(VariableNamesSelect)]
    head(BackgroundDatCorr.df)
    ### Run Spearman Rank Correlation and produce a correlation matrix on all variables
    AbsValsCor.mat <- abs(cor(BackgroundDatCorr.df, method="spearman"))
    # Find mean value of correlations in set
    LoAbsValsCor.mat <- (lower.tri(AbsValsCor.mat)*AbsValsCor.mat)
    MeanCor <- mean(LoAbsValsCor.mat)
    # Save correlation matrix
    write.table(AbsValsCor.mat, sep = ",", col.names=NA, file=paste0(Species, "_DataSet", SetName, "_CorrelationBackgroundEnvironVariables.csv"))
    ##
    CorrelationsIn <- AbsValsCor.mat
    SearchLoopsIn <- SearchLoops
    #
    if(FinalModelVariables==2) {
      if(ncol(combn(TotVars,2))>3000) {
        NumberModSets <- 3000
        NumberModSetsIn <- 3150
        FullNObs <- 250
        PlusNum <- 100 # extra random subsets
        #FullNObs <- NumberModSets/3
        #PlusNum <- 1
      } else {
        NumberModSets <- ncol(combn(TotVars,2))
        NumberModSetsIn <- NumberModSets
        FullNObs <- Round2(0.05*NumberModSets,0)
        PlusNum <- Round2(0.025*NumberModSets, 0) # extra random subsets
        #FullNObs <- NumberModSets/3
        #PlusNum <- 1
      }
    } else {
      NumberModSetsIn <- 3150 # Number sets per rep: extra sets if needed due to AICc values deleted- recommend no more than 5050 sets at once in case of crash
      NumberModSets <- 3000  # Number sets per rep
      #NumberModSetsIn <- 33 # extra sets if needed due to AICc values deleted
      #NumberModSets <- 30
      FullNObs <- 250
      PlusNum <- 100 # extra random subsets
      #FullNObs <- NumberModSets/3
      #PlusNum <- 1
    }
    #
    ##***** For check run
  #  NumberModSets <- 21
  #  NumberModSetsIn <- 24
  #  FullNObs <- NumberModSets/3
  #  PlusNum <- 1
    ##*******
    MaxModSetsIn <- NumberModSetsIn * ValidationReps # maximum number of sets needed to keep
    CorrThreshIn <- CorrThresh
    SetNameIn <- SetName
    OutDirectIn <- OutDirectpsa
    #ncol(combn(20,3))
    #
    system.time(VariableSubsets <- VariableSubsetSizeCorrThresh.Build(CorrelationsIn, SearchLoopsIn, SubsetVariableNumber, MaxModSetsIn, CorrThreshIn, SetNameIn, OutDirectIn))
    head(VariableSubsets)
    tail(VariableSubsets)
    # Find number of variables in subsets
    VariableNamesEx <- c(unlist(VariableSubsets[1,1]))
    # Check if dash used to separate variables
    if(grepl("-", VariableNamesEx)==TRUE) {
      VarNamesEx <- unlist(strsplit(VariableNamesEx, "-"))
    } else {
      VarNamesEx <- unlist(VariableNamesEx)
    }
    # Save VariableSubsets
    write.table(VariableSubsets, file=paste0(Species, "VariableSubsetsfor", FinalModelVariables, "Vars.csv"), sep=",", row.names=FALSE)
    ##
  #
  #  ### Make Heat Map of above matrix for correlation values
  #  # Convert correlation matrix to data frame
  #  ## Use function to round up from .5 from http://stackoverflow.com/questions/12688717/round-up-from-5-in-r
  #  Round2 <- function(x, n) {
  #    posneg = sign(x)
  #    z = abs(x)*10^n
  #    z = z + 0.5
  #    z = trunc(z)
  #    z = z/10^n
  #    z*posneg
  #  }
  #  DataFrameName <- data.frame(AbsValsCor.mat)
  #  PlotName <- paste0(LongSpecies, "_", SetNameF, "_CorrelationBackgroundEnvironVariables.tiff")
  #  PlotTitle <- paste0(LongSpecies, ": VariableCorrelations ", SetNameF, " from Background\n                                                    (Mean Correlation of ", Round2(MeanCor,2), ")\n")
  #  ## Call Heat Map Function
  #  HeatMapDataFrame.Plot(DataFrameName, PlotName, PlotTitle, OutDirectIn)
    ###
    #
    #FSAType <- Rank
    SetName2 <- paste0(SetName, InitVars, "CorrFilt", SubsetVariableNumber)
    TestDataType <- "Wrapper"
    #
    RunType <- paste0("RndFSAFilt_CorrThresh", CorrThresh)
    Run <- 1
    #
    SubsetSize.mat <- matrix(c("Singlets", "Doublets", "Triplets", "Quartets", "Quintets", "Sextets", "Septets", "Octets", "Nonets",
    "Dectets", "Undectets", "Duodectets","Tredectets", "Quattuordectets", "Quindectets", "Sexdectets", "Septendectets", "Octodectets", "Novemdectets",
    "Vigetets", "Unvigetets", "Duovigetets", "Trevigetets", "Quattuorvigetets", "Quinvigetets", "Sexvigetets", "Septenvigetets", "Octovigetet",
    "Novemvigetets", "Trigetets", "Untrigetets", "Duotrigetets", "Tretrigetets", "Quottuortrigetets", "Quintrigetets",
    "Sextrigetets", "Septentrigetets", "Octotrigetets", "Novemtrigetets", "Quadragetets", "Unquadragetets", "Duoquadragetets", "Trequadragetets",
    "Quattuorquadragetets", "Quinquadragetets", "Sexquadragetets", "Octoquadragetets", "Octoquadragetets", "Novemquadragetets", "Quinquagetets",
    "Unquinquagetets", "Duoquinquagetets", "Trequinguagetets", "Quattuorquinquagetets", "Quinquinquagetets",
    "Sexquinquagetets", "Septenquinquagetets", "Octoquinquagetets", "Novemquinquagetets", "Sexagetets"), ncol=1, nrow=60, byrow=TRUE, dimnames=list(c
     (seq(1:60)), c("Subset")))
    SubsetSize.df <- as.data.frame(SubsetSize.mat, stringsAsFactors=FALSE)
    Subset <- SubsetSize.df[SubsetVariableNumber,]
    if(is.na(Subset)) {
      Subset <- ""
    }
    #
    ##########################
    for(Rep in RepList) {
      #Rep=1
      # Read in k-fold partition scheme for wrapper cv
      setwd(InDirect)
      CVScheme <- Rep
      kfoldgrpp <- as.vector(unlist(read.csv(paste0(Species, "_PresenceDatRSFSAKfold_", CVScheme, "_.csv"), stringsAsFactors=FALSE)))
      length(kfoldgrpp)
      #
      kfoldgrpa <- as.vector(unlist(read.csv(paste0(Species, "_PseudoabsenceDatRSFSAKfold_", CVScheme, "_.csv"), stringsAsFactors=FALSE)))
      length(kfoldgrpa)
      ####################
      ## Evaluate random subsets for Sets number of subsets
      # Takes 4.48 minutes for 1000 runs
      #ncol(combn(19,5))
      #  t1a <- Sys.time()
      #
      MaxentPresTrainData <- PresenceDat.df[kfoldgrpp==1,]
      head(MaxentPresTrainData)
      MaxentAbsTrainData <- BackgroundDat.df
      head(MaxentAbsTrainData)
      WrapperTrainSWD <- rbind(MaxentPresTrainData[,4:ncol(MaxentPresTrainData)], MaxentAbsTrainData[,4:ncol(MaxentAbsTrainData)])
      head(WrapperTrainSWD)
      tail(WrapperTrainSWD)
      WrapperTrainPresID <- data.frame(rep(1,nrow(MaxentPresTrainData)))
      colnames(WrapperTrainPresID) <- "ID"
      WrapperTrainAbsID <- data.frame(rep(0,nrow(MaxentAbsTrainData)))
      colnames(WrapperTrainAbsID) <- "ID"
      WrapperTrainPresAbsID <- rbind(WrapperTrainPresID, WrapperTrainAbsID)
      head(WrapperTrainPresAbsID)
      tail(WrapperTrainPresAbsID)
      ##############
      ## Partition VariableSubsets according to Rep
      RepSubsets <- nrow(VariableSubsets)
      RepSubsetMinList <- c(1, (RepSubsets/3)+1, ((RepSubsets/3)*2)+1)
      RepSubsetMaxList <- c((RepSubsets/3), ((RepSubsets/3)*2), RepSubsets)
      RepSubsetMin <- RepSubsetMinList[Rep]
      RepSubsetMax <- RepSubsetMaxList[Rep]
      VariableSubsetsRep <- as.data.frame(VariableSubsets[RepSubsetMin:RepSubsetMax,], stringsAsFactors=FALSE)
      colnames(VariableSubsetsRep) <- "VarNames"
      MaxentWrapperTrainTestSubsetEvalStatsL <- list()
      SetCount <- 0
      Sets <- nrow(VariableSubsetsRep)
      #
      output2 <- paste0(NumberModSets,"Setsof", Subset, SubsetVariableNumber, "of", TotVars, "Vars", "_Rep", Rep)
      OutDirect2 <- paste0(OutDirectpsa, "/", output2)
      dir.create(OutDirect2)
      setwd(OutDirect2)
      OutDirectIn <- OutDirect2
      #
      ### Specify three equal RunSize groups of VariableSubsetsRep of no more than around 5,050
      RunSize <- Round2(Sets/3,0)
      NumberModSetsIn <- Sets
      SetRunsList <- list(seq(1,RunSize,1), seq(RunSize+1,RunSize*2,1), seq(RunSize*2+1,Sets,1))
      # Groups of VariableSubsets
      ### NOTE: Can comment this loop if just re-ranking wrapper models by different statistic such as AICc
      ####################
      if(Rank!="AIC") {
        for(SetRun in SetRunsList) {
          #SetRun <- SetRunsList[[2]]
          VariableSubsetsIn <- as.data.frame(VariableSubsetsRep[SetRun,], stringsAsFactors=FALSE)
          colnames(VariableSubsetsIn) <- "VarNames"
          SetCount <- SetCount + 1
          SetRunID <- SetCount
          MaxentWrapperTrainTestSubsetEvalStatsList <- list()
          system.time(MaxentWrapperTrainTestSubsetEvalStatsList <- MaxentMultiSubset.WrapperCV1TrainTestEvalAIC(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsetsIn, SetRunID, WrapperTrainSWD, WrapperTrainPresAbsID, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, OutDirectIn, MaxentArgsIn=MaxentBaseArgsIn))
          MaxentWrapperTrainTestSubsetEvalStats.df1 <- MaxentWrapperTrainTestSubsetEvalStatsList[[1]]
          MaxentWrapperTrainTestSubsetPermImp.df <- MaxentWrapperTrainTestSubsetEvalStatsList[[2]]
          ### Assign ranks to variables using MMOORA joint ranking by permutation importance
          ### and number of times variables appears in models
        }
      }
      ######################
    #} # For ending loop here
    #t2 <- Sys.time()
    ######################################
    #difftime(t2,t1, units = "mins")
      # Read back in data
      if(Rank=="AIC"){
        OutDirect2Orig <- gsub("AIC", "AUC", OutDirect2)
        setwd(OutDirect2Orig)
      } else {
        setwd(OutDirect2)
      }
      SetRunIDList <- c(1:3)
      #SetRunIDList <- c(1)
      #RowMultList <- c(2)
      RowCount <- 0
      MaxentWrapperTrainTestSubsetEvalStatsL <- list()
      MaxentWrapperTrainTestSubsetPermImpL <- list()
      for(SetRunID in SetRunIDList) {
        #SetRunID<-1
        #NumberSets <- length(SetRunsList[[SetRunID]])
        # For re-reading in data
        NumberSets <- length(SetRunsList[[SetRunID]])
        MaxentWrapperTrainTestSubsetEvalStats.df1 <-  as.data.frame(read.csv(paste0(Species, "MaxentResults_WrapperCV1TrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", NumberSets, "_", SetRunID, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
        head(MaxentWrapperTrainTestSubsetEvalStats.df1)
        tail(MaxentWrapperTrainTestSubsetEvalStats.df1)
        nrow(MaxentWrapperTrainTestSubsetEvalStats.df1)
        # Convert first, second and fourth columns from factor to character
        MaxentWrapperTrainTestSubsetEvalStats.df1[,c(1,2,4)] <- sapply(MaxentWrapperTrainTestSubsetEvalStats.df1[,c(1,2,4)], function(x) as.character(x))
        Rows <- nrow(MaxentWrapperTrainTestSubsetEvalStats.df1)
        #str(MaxentWrapperTrainTestSubsetEvalStats.df1)
        #Renumber rows
        if(RowCount==0) {
          RowCount <- Rows
        } else {
          RowCount <- RowCount + Rows
          row.names(MaxentWrapperTrainTestSubsetEvalStats.df1) <- seq(RowCount-Rows+1, RowCount, 1)
        }
        MaxentWrapperTrainTestSubsetEvalStatsL[[SetRunID]] <- MaxentWrapperTrainTestSubsetEvalStats.df1
        ###
        MaxentWrapperTrainTestSubsetPermImp.df1 <- as.data.frame(read.csv(paste0(Species, "MaxentPermutationImportances_WrapperCV1TrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", NumberSets, "_", SetRunID, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
        nrow(MaxentWrapperTrainTestSubsetPermImp.df1)
        head(MaxentWrapperTrainTestSubsetPermImp.df1)
        #str(MaxentWrapperTrainTestSubsetPermImp.df1)
        # Convert first and second columns from factor to character
        MaxentWrapperTrainTestSubsetPermImp.df1[,c(1,2)] <- sapply(MaxentWrapperTrainTestSubsetPermImp.df1[,c(1,2)], function(x) as.character(x))
        MaxentWrapperTrainTestSubsetPermImpL[[SetRunID]] <- MaxentWrapperTrainTestSubsetPermImp.df1
      }
      MaxentWrapperTrainTestSubsetEvalStats.df <- do.call(rbind, MaxentWrapperTrainTestSubsetEvalStatsL)
      head(MaxentWrapperTrainTestSubsetEvalStats.df)
      nrow(MaxentWrapperTrainTestSubsetEvalStats.df)
      #MaxentWrapperTrainTestSubsetEvalStats.df[grepl("BLKDEN30CM-CONT_DFMO-IM_INDEX", MaxentWrapperTrainTestSubsetEvalStats.df$VarNames), ]
      #
      #str(MaxentWrapperTrainTestSubsetEvalStats.df)
      MaxentWrapperTrainTestSubsetPermImp.df <- do.call(rbind, MaxentWrapperTrainTestSubsetPermImpL)
      nrow(MaxentWrapperTrainTestSubsetPermImp.df)
      #MaxentWrapperTrainTestSubsetPermImp.df[grepl("BLKDEN30CM-CONT_DFMO", MaxentWrapperTrainTestSubsetPermImp.df$VarNames), ]
      ## Separate out training and testing data
      MaxentTrainSubsetEvalStats.df <- MaxentWrapperTrainTestSubsetEvalStats.df[which(MaxentWrapperTrainTestSubsetEvalStats.df$DataType=="WrapperCV1Train"),]
      MaxentTestSubsetEvalStats.df <- MaxentWrapperTrainTestSubsetEvalStats.df[which(MaxentWrapperTrainTestSubsetEvalStats.df$DataType=="WrapperCV1Test"),]
      # Calculate difference between test and train statistics for overfitting
      MaxentDiffTestTrainSubsetEvalStats.df <- MaxentTrainSubsetEvalStats.df
      MaxentDiffTestTrainSubsetEvalStats.df[,5:13] <- MaxentTrainSubsetEvalStats.df[,5:13] - MaxentTestSubsetEvalStats.df[,5:13]
      MaxentDiffTestTrainSubsetEvalStats.df$DataType <- "WrapperCV1DiffTestTrain"
      #
      # Join AUCDiff with each data set
      MaxentDiffTestTrainSubsetEvalStats.df$AUCdiff <- MaxentDiffTestTrainSubsetEvalStats.df$AUC
      MaxentTrainSubsetEvalStats.df$AUCdiff <- MaxentDiffTestTrainSubsetEvalStats.df$AUC
      MaxentTestSubsetEvalStats.df$AUCdiff <- MaxentDiffTestTrainSubsetEvalStats.df$AUC
      #
      MaxentWrapperTrainTestSubsetEvalStats.df <- do.call(rbind, list(MaxentTrainSubsetEvalStats.df, MaxentTestSubsetEvalStats.df, MaxentDiffTestTrainSubsetEvalStats.df))
      ######################
      ## Loop through different numbers of total subsets from 250 to 500 to 1,000 for ranked and random subsets
      ChosenModSetSizes <- c(Sets)
      #ChosenModSetSizes <- c(1000, 5000, 10000)
      for(ChosenModSets in ChosenModSetSizes) {
        #ChosenModSets=Sets
        ## Rank top NObs subsets by Wrapper TSS and obtain training and final testing TSS for evaluation of each set
        ## and for ensemble set
        ##
        ## NOTE: The next two lines are for running the loop with reranking and testing only
        if(Rank=="AIC"){
          OutDirect2Orig <- gsub("AIC", "AUC", OutDirect2)
          setwd(OutDirect2Orig)
        } else {
          setwd(OutDirect2)
        }
        OutCount <- OutCount + 1
        Run <- 1
        #length(unique(MaxentFiltSubsetEvalStats.df$VarNames)) # check for duplicate VarNames
        ## Separate out WrapperTest from data
        MaxentFiltSubsetEvalStats.df <- MaxentWrapperTrainTestSubsetEvalStats.df[which(MaxentWrapperTrainTestSubsetEvalStats.df$DataType=="WrapperCV1Test"),]
        # Delete any subsets with "Inf" for AICc_bg (means too many parameters in model relative to observations)
        MaxentFiltSubsetEvalStats.df2 <- subset(MaxentFiltSubsetEvalStats.df, AICc_bg!="Inf")
        NewSets1 <- nrow(MaxentFiltSubsetEvalStats.df2)
        tail(MaxentFiltSubsetEvalStats.df2)
        ## Keep only first ChosenModSets number of sets
        MaxentFiltSubsetEvalStats.df3 <- MaxentFiltSubsetEvalStats.df2[1:ChosenModSets,]
        head(MaxentFiltSubsetEvalStats.df3)
        tail(MaxentFiltSubsetEvalStats.df3)
        nrow(MaxentFiltSubsetEvalStats.df3)
        ## Use Multi-Objective Optimization from MCDM package to sort by two criteria of AUC and AICc_bg
        decision.mat <- as.matrix(MaxentFiltSubsetEvalStats.df3[,c(11,7,14)])
        head(decision.mat)
        # Normalize AICc_bg and AUC from zero to 1
        normalizef <- function(x) c((x+(max(x)-min(x))-max(x))/(max(x)-min(x)))
        decision.matn <- as.matrix(apply(decision.mat,2,normalizef))
        head(decision.matn)
        ######## Make rank schemes
        ##
        AICcrankList <- c(AICcrank)
        AUCrankList <- c(AUCrank)
        AUCdiffrankList <- c(0)
        MaxentTopNObsFiltSubsetEvalStatsL <- list()
        TopNObsSubsetsL <- list()
        for(m in 1:length(AUCrankList)) {
          # Specify ranks for AICc and AUC in Multi-Object Optimization ranking
          AICcrank <- AICcrankList[m]
          AUCrank <- AUCrankList[m]
          AUCdiffrank <- AUCdiffrankList[m]
          ##########
          # Assign weights to criteria
          weightscrit <- c(AICcrank, AUCrank, AUCdiffrank)
          # Assign whether cost "min", or benefit "max"
          cb <- c("min", "max", "min")
          # Rank criteria
          MMOORA.Rank <- MMOORA(decision.matn, weightscrit, cb)
          head(MMOORA.Rank)
          MaxentFiltSubsetEvalStats.df4 <- MaxentFiltSubsetEvalStats.df3
          # Join MMOORA rank to original data
          MaxentFiltSubsetEvalStats.df4$MMOORA_Rank <- MMOORA.Rank[,8]
          head(MaxentFiltSubsetEvalStats.df4)
          # Sort data by MMOORA rank
          MaxentFiltSubsetEvalStats.df4s <- MaxentFiltSubsetEvalStats.df4[order(MaxentFiltSubsetEvalStats.df4$MMOORA_Rank),]
          head(MaxentFiltSubsetEvalStats.df4s)
          tail(MaxentFiltSubsetEvalStats.df4s)
          nrow(MaxentFiltSubsetEvalStats.df4s)
          ncol(MaxentFiltSubsetEvalStats.df4s)
          # Delete rank column
          MaxentFiltSubsetEvalStats.df5 <- MaxentFiltSubsetEvalStats.df4s[,-15]
          head(MaxentFiltSubsetEvalStats.df5)
          # Keep top FullNObs selected subsets
          if(nrow(MaxentFiltSubsetEvalStats.df5)>FullNObs) {
            MaxentTopNObsFiltSubsetEvalStatsL[[m]] <- data.frame(MaxentFiltSubsetEvalStats.df5[1:FullNObs,], stringsAsFactors=FALSE)
            TopNObsSubsetsL[[m]] <- data.frame(MaxentFiltSubsetEvalStats.df5[1:FullNObs,1], stringsAsFactors=FALSE)
          } else {
            MaxentTopNObsFiltSubsetEvalStatsL[[m]] <- MaxentFiltSubsetEvalStats.df5
            TopNObsSubsetsL[[m]] <- MaxentFiltSubsetEvalStats.df5
          }
        }
        #######################
        MaxentTopNObsFiltSubsetEvalStats.df <- do.call(rbind, MaxentTopNObsFiltSubsetEvalStatsL)
        nrow(MaxentTopNObsFiltSubsetEvalStats.df)
        TopNObsSubsets <- do.call(rbind, TopNObsSubsetsL)
        colnames(TopNObsSubsets) <- "VarNames"
        nrow(TopNObsSubsets)
        #str(TopNObsSubsets)
        #TopNObsSubsets[grepl("BLKDEN30CM-CONT_DFMO", TopNObsSubsets$VarNames), ]
        ##############
        ### Calculate Variable Rankings for these sets
        ## Subset by TopNObsSubsets
        TopMaxentWrapperTrainTestSubsetPermImp.df1 <- subset(MaxentWrapperTrainTestSubsetPermImp.df, VarNames %in% TopNObsSubsets$VarNames)
        TopMaxentWrapperTrainTestSubsetPermImp.df1 <- MaxentWrapperTrainTestSubsetPermImp.df[MaxentWrapperTrainTestSubsetPermImp.df$VarNames %in% TopNObsSubsets$VarNames, ]
        #MaxentWrapperTrainTestSubsetPermImp.df$VarNames
        nrow(TopMaxentWrapperTrainTestSubsetPermImp.df1)
        #MaxentWrapperTrainTestSubsetPermImp.df[grepl("BLKDEN30CM-CONT_DFMO", MaxentWrapperTrainTestSubsetPermImp.df$VarNames), ]
        #TopMaxentWrapperTrainTestSubsetPermImp.df1[grepl("BLKDEN30CM-CONT_DFMO", TopMaxentWrapperTrainTestSubsetPermImp.df1$VarNames), ]
        ## Fix row numbers
        Rows <- nrow(TopMaxentWrapperTrainTestSubsetPermImp.df1)
        row.names(TopMaxentWrapperTrainTestSubsetPermImp.df1) <- seq(1, Rows, 1)
        ##
        # Drop all permutation importance values of zero since variable not used
        TopMaxentWrapperTrainTestSubsetPermImp.df1 <- TopMaxentWrapperTrainTestSubsetPermImp.df1[apply(TopMaxentWrapperTrainTestSubsetPermImp.df1!=0,1,all),]
        ## Sort by Variable
        TopMaxentWrapperTrainTestSubsetPermImp.dfs <- TopMaxentWrapperTrainTestSubsetPermImp.df1[order(TopMaxentWrapperTrainTestSubsetPermImp.df1$Variable),]
        #str(PermutationImp.dfs)
        ## Calculate mean and standard deviation per VarName
        PermutationImp.means <- aggregate(Permutation_Importance ~ Variable, TopMaxentWrapperTrainTestSubsetPermImp.dfs, mean)
        PermutationImp.sd <- aggregate(Permutation_Importance ~ Variable, TopMaxentWrapperTrainTestSubsetPermImp.dfs, sd)
        # Replace NAs for sd with zero
        PermutationImp.sd[is.na(PermutationImp.sd)] <- 0
        # Calculate count for each variable
        PermutationImp.count <- aggregate(Permutation_Importance ~ Variable, TopMaxentWrapperTrainTestSubsetPermImp.dfs, FUN = function(x){NROW(x)})
        PermutationImp <- cbind(PermutationImp.means, PermutationImp.sd[,2], PermutationImp.count[,2])
        colnames(PermutationImp) <- c("Variable", "MeanPermImp", "SDPermImp", "Count")
        ## Sort by MeanPermImp
        PermutationImp <- PermutationImp[order(-PermutationImp$MeanPermImp),]
        ## Create column with joint ranking by MeanPermImp and Count using Multi-Objective Optimization ranking
        ## Use Multi-Objective Optimization from MCDM package to sort by two criteria of MeanPermImp and Count
        decision.mat <- as.matrix(PermutationImp[,c(2,4)])
        head(decision.mat)
        # Normalize variable mean permutation importance and variable frequency in models from zero to 1
        normalizef <- function(x) c((x+(max(x)-min(x))-max(x))/(max(x)-min(x)))
        decision.matn <- as.matrix(apply(decision.mat,2,normalizef))
        head(decision.matn)
        # Assign weights to criteria
        weightscrit <- c(PermImpRank, CountRank) # PermImpRank and CountRank assigned at beginning of program
        # Assign whether cost "min", or benefit "max"
        cb <- c("max", "max")
        # Rank criteria
        MMOORA.Rank <- MMOORA(decision.matn, weightscrit, cb)
        head(MMOORA.Rank)
        # Join MMOORA rank to original data
        PermutationImp$PermImp_MMOORA_Rank <- MMOORA.Rank[,8]
        # Sort data by MMOORA rank
        PermutationImp2 <- PermutationImp[order(PermutationImp$PermImp_MMOORA_Rank),]
        head(PermutationImp2)
        # Add a column for Rep
        PermutationImp2$Rep <- Rep
        ## Save PermImp_MMOORA_Ranks for variables
        # Write Data
        write.table(PermutationImp2, file=paste0(Species, "TopSetVarsMaxentMeanPermImpRanks_WrapperCV1TrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", nrow(PermutationImp2), "VarsTotal_Rep",  Rep, ".csv"), sep=",", col.names=NA)
        TopVariablePermImpL[[Rep]] <-  PermutationImp2
        #
        # Assign PermImp_MMOORA_Rank to Variables in original file with VarNames
        TopMaxentWrapperTrainTestSubsetPermImp.dfs$PermImp_MMOORA_Rank <- PermutationImp2$PermImp_MMOORA_Rank[match(TopMaxentWrapperTrainTestSubsetPermImp.dfs$Variable, PermutationImp2$Variable)]
        # Calculate the mean PermImp_MMOORA_Rank by VarNames
        TopMaxentWrapperTrainTestSubsetPermImp.means <- aggregate(PermImp_MMOORA_Rank ~ VarNames, TopMaxentWrapperTrainTestSubsetPermImp.dfs, mean)
        # Sort by ranks
        TopMaxentWrapperTrainTestSubsetPermImp.means.dfs <- TopMaxentWrapperTrainTestSubsetPermImp.means[order(TopMaxentWrapperTrainTestSubsetPermImp.means$PermImp_MMOORA_Rank),]
        head(TopMaxentWrapperTrainTestSubsetPermImp.means.dfs)
        # Add a column for Rep
        TopMaxentWrapperTrainTestSubsetPermImp.means.dfs$Rep <- Rep
        # Write Data
        write.table(TopMaxentWrapperTrainTestSubsetPermImp.means.dfs, file=paste0(Species, "TopSetMaxentMeanPermImpRanks_WrapperCV1TrainTest_", TotVars, "TotVars_", SubsetVariableNumber, "Vars_", nrow(TopMaxentWrapperTrainTestSubsetPermImp.means.dfs), "Sets_Rep",  Rep, ".csv"), sep=",", col.names=NA)
        TopMaxentSubsetPermImpL[[Rep]] <-  TopMaxentWrapperTrainTestSubsetPermImp.means.dfs
        ##############################
        # Keep random FullNObs plus PlusNum random subsets
        # If Rank is AUC, keep Random subsets
        if(Rank=="AUC") {
          if(nrow(MaxentFiltSubsetEvalStats.df3)>FullNObs+PlusNum) {
            MaxentRandNObsFiltSubsetEvalStats.df1 <- data.frame(MaxentFiltSubsetEvalStats.df3[sample(nrow(MaxentFiltSubsetEvalStats.df3), FullNObs+PlusNum), ], stringsAsFactors=FALSE)
          } else {
            MaxentRandNObsFiltSubsetEvalStats.df1 <- MaxentFiltSubsetEvalStats.df3
          }
          ## Clean random data from aberrant AICc_bg values
          # Delete any subsets with "Inf" for AICc_bg (means too many parameters in model relative to observations)
          MaxentRandNObsFiltSubsetEvalStats.df2 <- subset(MaxentRandNObsFiltSubsetEvalStats.df1, AICc_bg!="Inf")
          # Delete any subsets with negative AICc_bg which appear to be spurious values
          MaxentRandNObsFiltSubsetEvalStats.df3 <- subset(MaxentRandNObsFiltSubsetEvalStats.df2, AICc_bg>-1)
          # Keep only FullNObs plus PlusNum/2, or around 50 sets
          if(nrow(MaxentRandNObsFiltSubsetEvalStats.df3)>FullNObs+Round2(PlusNum/2,0)) {
            MaxentRandNObsFiltSubsetEvalStats.df <- MaxentRandNObsFiltSubsetEvalStats.df3[1:(FullNObs+Round2(PlusNum/2,0)),]
          } else {
            MaxentRandNObsFiltSubsetEvalStats.df  <- MaxentRandNObsFiltSubsetEvalStats.df3
          }
          head(MaxentRandNObsFiltSubsetEvalStats.df)
          tail(MaxentRandNObsFiltSubsetEvalStats.df)
          nrow(MaxentRandNObsFiltSubsetEvalStats.df)
          ncol(MaxentRandNObsFiltSubsetEvalStats.df)
          #
          RandNObsSubsets <- data.frame(MaxentRandNObsFiltSubsetEvalStats.df[,1], stringsAsFactors=FALSE)
          colnames(RandNObsSubsets) <- "VarNames"
        }
        # Use difference between evaluation statistics between train and test data to evaluate overfitting
        ###############################################################################
        ## Test Top FullNObs Ranking Subsets and Random FullNObs Subsets
        ###############################################################################
        #####################################
        setwd(OutDirect2)
        #SubsetRuns <- list(TopNObsSubsets[1:NObs,], TopNObsSubsets[(NObs+1):(NObs*2),], TopNObsSubsets[(NObs*2+1):(NObs*3),], TopNObsSubsets[(NObs*3+1):(NObs*4),],RandNObsSubsets[1:NObs,], RandNObsSubsets[(NObs+1):(NObs*2),], RandNObsSubsets[(NObs*2+1):(NObs*3),])
        #DataSetTypes <- c(paste0("TopAUC",NObs), paste0("TopAICc",NObs),  paste0("TopAUCdiff",NObs), paste0("TopAUC_AICc_AUCdiff",NObs), paste0("Random", NObs, "A"), paste0("Random", NObs, "B"), paste0("Random", NObs, "C"))
        if(Rank=="AUC") {
          SubsetRuns <- list(TopNObsSubsets, RandNObsSubsets)
          DataSetTypes <- c(paste0("Top",FullNObs), paste0("Random", FullNObs))
        } else {
          SubsetRuns <- list(TopNObsSubsets)
          DataSetTypes <- c(paste0("Top",FullNObs))
        }
        MaxentKeepEvalStatsL1 <- list()
        #
        for (j in 1:length(DataSetTypes)) {
          #N3Subsets <- SubsetRuns[[1]]
          #j=1
          VariableSubsetsIn <- data.frame(SubsetRuns[[j]], stringsAsFactors=FALSE)
          DataSetType <- DataSetTypes[j]
          # Identify Wrapper data evaluation statistics from above in loop
          if(DataSetType==paste0("Top", FullNObs)) {
            MaxentTestSubsetEvalStats.df <- MaxentTopNObsFiltSubsetEvalStats.df
          } else {
            MaxentTestSubsetEvalStats.df <- MaxentRandNObsFiltSubsetEvalStats.df
          }
          MaxentTestSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperCV1Test")
          ##
          ## Join Wrapper Test data back with Wrapper Train and WrapperDiffTestTrain
          ## Subset original data using VarNames
          selectedRows <- (MaxentWrapperTrainTestSubsetEvalStats.df$VarNames %in% MaxentTestSubsetEvalStats.df$VarNames)
          MaxentTestSubsetEvalStats.df2 <- MaxentWrapperTrainTestSubsetEvalStats.df[selectedRows,]
          ## Take apart by DataType, reorder each by rank and reassemble
          MaxentTrainSubsetEvalStats.df3 <- MaxentTestSubsetEvalStats.df2[which(MaxentTestSubsetEvalStats.df2$DataType=="WrapperCV1Train"),]
          MaxentTrainSubsetEvalStats.df4 <- MaxentTrainSubsetEvalStats.df3[match(unlist(MaxentTestSubsetEvalStats.df$VarNames), MaxentTrainSubsetEvalStats.df3$VarNames), ]
          #
          MaxentTestSubsetEvalStats.df3 <- MaxentTestSubsetEvalStats.df2[which(MaxentTestSubsetEvalStats.df2$DataType=="WrapperCV1Test"),]
          MaxentTestSubsetEvalStats.df4 <- MaxentTestSubsetEvalStats.df3[match(unlist(MaxentTestSubsetEvalStats.df$VarNames), MaxentTestSubsetEvalStats.df3$VarNames), ]
          #
          MaxentDiffTestTrainSubsetEvalStats.df3 <- MaxentTestSubsetEvalStats.df2[which(MaxentTestSubsetEvalStats.df2$DataType=="WrapperCV1DiffTestTrain"),]
          MaxentDiffTestTrainSubsetEvalStats.df4 <- MaxentDiffTestTrainSubsetEvalStats.df3[match(unlist(MaxentTestSubsetEvalStats.df$VarNames), MaxentDiffTestTrainSubsetEvalStats.df3$VarNames), ]
          # Reassemble
          MaxentSubsetEvalStats.df <- rbind(MaxentTrainSubsetEvalStats.df4, MaxentTestSubsetEvalStats.df4, MaxentDiffTestTrainSubsetEvalStats.df4)
          MaxentSubsetEvalStats.df$DataType <- paste0(DataSetType, MaxentSubsetEvalStats.df$DataType)
          MaxentSubsetEvalStats.df$SetName <- SetName
          MaxentSubsetEvalStats.df$Set <- DataSetType
          MaxentSubsetEvalStats.df$TotSubsets <- ChosenModSets
          MaxentSubsetEvalStats.df$PresTrainPoints <- nrow(PresenceDat.df[kfoldgrpp==1,])
          MaxentSubsetEvalStats.df$Rep <- Rep
          ncol(MaxentSubsetEvalStats.df)
          MaxentSubsetEvalStats.df <- MaxentSubsetEvalStats.df[, colnames(MaxentSubsetEvalStats.df)[c(1:4,15:19,5:14)]]
          #
          MaxentKeepEvalStatsL1[[j]] <- MaxentSubsetEvalStats.df
        }
        ###############
        ###########################
        # Read in k-fold partition scheme for wrapper cv
        setwd(InDirect)
        CVScheme <- Rep+(CVNum/2)
        kfoldgrpp <- as.vector(unlist(read.csv(paste0(Species, "_PresenceDatRSFSAKfold_", CVScheme, "_.csv"), stringsAsFactors=FALSE)))
        length(kfoldgrpp)
        #
        kfoldgrpa <- as.vector(unlist(read.csv(paste0(Species, "_PseudoabsenceDatRSFSAKfold_", CVScheme, "_.csv"), stringsAsFactors=FALSE)))
        length(kfoldgrpa)
        ####################
        ## Designate Final Training data for Final Test
        MaxentPresTrainData <- PresenceDat.df[kfoldgrpp==1,]
        head(MaxentPresTrainData)
        nrow(MaxentPresTrainData)
        MaxentAbsTrainData <- BackgroundDat.df
        nrow(MaxentAbsTrainData)
        head(MaxentAbsTrainData)
        WrapperTrainSWD <- rbind(MaxentPresTrainData[,4:ncol(MaxentPresTrainData)], MaxentAbsTrainData[,4:ncol(MaxentAbsTrainData)])
        head(WrapperTrainSWD)
        tail(WrapperTrainSWD)
        WrapperTrainPresID <- data.frame(rep(1,nrow(MaxentPresTrainData)))
        colnames(WrapperTrainPresID) <- "ID"
        WrapperTrainAbsID <- data.frame(rep(0,nrow(MaxentAbsTrainData)))
        colnames(WrapperTrainAbsID) <- "ID"
        WrapperTrainPresAbsID <- rbind(WrapperTrainPresID, WrapperTrainAbsID)
        head(WrapperTrainPresAbsID)
        tail(WrapperTrainPresAbsID)
        #
        setwd(OutDirect2)
        #SubsetRuns <- list(TopNObsSubsets[1:NObs,], TopNObsSubsets[(NObs+1):(NObs*2),], TopNObsSubsets[(NObs*2+1):(NObs*3),], TopNObsSubsets[(NObs*3+1):(NObs*4),],RandNObsSubsets[1:NObs,], RandNObsSubsets[(NObs+1):(NObs*2),], RandNObsSubsets[(NObs*2+1):(NObs*3),])
        #DataSetTypes <- c(paste0("TopAUC",NObs), paste0("TopAICc",NObs),  paste0("TopAUCdiff",NObs), paste0("TopAUC_AICc_AUCdiff",NObs), paste0("Random", NObs, "A"), paste0("Random", NObs, "B"), paste0("Random", NObs, "C"))
        if(Rank=="AUC") {
          SubsetRuns <- list(TopNObsSubsets, RandNObsSubsets)
          DataSetTypes <- c(paste0("Top",FullNObs), paste0("Random", FullNObs))
        } else {
          SubsetRuns <- list(TopNObsSubsets)
          DataSetTypes <- c(paste0("Top",FullNObs))
        }
        MaxentKeepEvalStatsL2 <- list()
        #
        for (j in 1:length(DataSetTypes)) {
          #j=1
          VariableSubsetsIn <- data.frame(SubsetRuns[[j]], stringsAsFactors=FALSE)
          DataSetType <- DataSetTypes[j]
          ## Evaluate random subsets for NumberModSets number of subsets
          system.time(MaxentTrainTestSubsetEvalStats.df1 <- MaxentMultiSubset.WrapperCV2TrainTestEvalAIC(Species, VariableNamesIn, SubsetVariableNumber, TotPres, VariableSubsetsIn, WrapperTrainSWD, WrapperTrainPresAbsID, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrpp, kfoldgrpa, OutDirectIn, MaxentArgsIn=MaxentBaseArgsIn, DataSetType=DataSetType))
          # Delete any subsets with "Inf" for AICc_bg (means too many parameters in model relative to observations)
          MaxentTrainTestSubsetEvalStats.df2 <- subset(MaxentTrainTestSubsetEvalStats.df1, AICc_bg!="Inf")
          nrow(MaxentTrainTestSubsetEvalStats.df2)
          # Delete any subsets with negative values of AICc_bg which appear to be spurious values (equal to or less than -1)
          MaxentTrainTestSubsetEvalStats.df <- subset(MaxentTrainTestSubsetEvalStats.df2, AICc_bg>-1)
          nrow(MaxentTrainTestSubsetEvalStats.df)
          ## Separate out training and testing data
          MaxentTrainSubsetEvalStats.df <- MaxentTrainTestSubsetEvalStats.df[which(MaxentTrainTestSubsetEvalStats.df$DataType=="WrapperCV2Train"),]
          MaxentTrainSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperCV2Train")
          MaxentTestSubsetEvalStats.df <- MaxentTrainTestSubsetEvalStats.df[which(MaxentTrainTestSubsetEvalStats.df$DataType=="WrapperCV2Test"),]
          MaxentTestSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperCV2Test")
          # Calculate difference between test and train statistics for overfitting
          MaxentDiffTestTrainSubsetEvalStats.df <- MaxentTestSubsetEvalStats.df
          MaxentDiffTestTrainSubsetEvalStats.df[,5:13] <- MaxentTrainSubsetEvalStats.df[,5:13] - MaxentTestSubsetEvalStats.df[,5:13]
          MaxentDiffTestTrainSubsetEvalStats.df$DataType <- paste0(DataSetType, "WrapperCV2DiffTestTrain")
          #
          # Join AUCDiff with each data set
          MaxentDiffTestTrainSubsetEvalStats.df$AUCdiff <- MaxentDiffTestTrainSubsetEvalStats.df$AUC
          MaxentTrainSubsetEvalStats.df$AUCdiff <- MaxentDiffTestTrainSubsetEvalStats.df$AUC
          MaxentTestSubsetEvalStats.df$AUCdiff <- MaxentDiffTestTrainSubsetEvalStats.df$AUC
          #
          MaxentSubsetEvalStats.df <- do.call(rbind, list(MaxentTrainSubsetEvalStats.df, MaxentTestSubsetEvalStats.df, MaxentDiffTestTrainSubsetEvalStats.df))
          MaxentSubsetEvalStats.df$SetName <- SetName
          MaxentSubsetEvalStats.df$Set <- DataSetType
          MaxentSubsetEvalStats.df$TotSubsets <- ChosenModSets
          MaxentSubsetEvalStats.df$PresTrainPoints <- nrow(PresenceDat.df[kfoldgrpp==1,])
          MaxentSubsetEvalStats.df$Rep <- Rep
          MaxentSubsetEvalStats.df <- MaxentSubsetEvalStats.df[, colnames(MaxentSubsetEvalStats.df)[c(1:4,15:19,5:14)]]
          #
          MaxentKeepEvalStatsL2[[j]] <- MaxentSubsetEvalStats.df
        }
        ###############
        MaxentKeepEvalStats.df1 <- do.call(rbind, MaxentKeepEvalStatsL1)
        MaxentKeepEvalStats.df2 <- do.call(rbind, MaxentKeepEvalStatsL2)
        MaxentKeepEvalStats.df3 <- rbind(MaxentKeepEvalStats.df1, MaxentKeepEvalStats.df2)
        nrow(MaxentKeepEvalStats.df3)
        # Eliminate duplicates for EnvVarsUsed for a given DataType
        # First randomize sorting by rep to even out removals
        MaxentKeepEvalStats.df2 <- MaxentKeepEvalStats.df3[sample(nrow(MaxentKeepEvalStats.df3)), ]
        head(MaxentKeepEvalStats.df2)
        tail(MaxentKeepEvalStats.df2)
        nrow(MaxentKeepEvalStats.df2)
        # Select out the random wrapper test DataType
        DataTypeSel <- paste0("Random", FullNObs, "WrapperCV2Test")
        MaxentKeepEvalStatsRandPart.df <- MaxentKeepEvalStats.df2[MaxentKeepEvalStats.df2$DataType==DataTypeSel,]
        nrow(MaxentKeepEvalStatsRandPart.df)
        # Remove duplicates
        MaxentKeepEvalStatsRandPart.df2 <- MaxentKeepEvalStatsRandPart.df[!duplicated(MaxentKeepEvalStatsRandPart.df$EnvVarsUsed), ]
        head(MaxentKeepEvalStatsRandPart.df2)
        nrow(MaxentKeepEvalStatsRandPart.df2)
        # Save only VarNames
        RandVarNamesKeep <- MaxentKeepEvalStatsRandPart.df2$VarNames
        length(RandVarNamesKeep)
        ## Keep only above VarNames in all Random sets
        # First isolate all Random Sets
        MaxentKeepEvalStatsRand.df <- MaxentKeepEvalStats.df2[grepl("Rand", MaxentKeepEvalStats.df2$DataType), ]
        head(MaxentKeepEvalStatsRand.df)
        tail(MaxentKeepEvalStatsRand.df)
        nrow(MaxentKeepEvalStatsRand.df)
        # Only keep original data with VarNames matching that in above data
        MaxentKeepEvalStatsRand.df1 <- MaxentKeepEvalStatsRand.df[MaxentKeepEvalStatsRand.df$VarNames %in% RandVarNamesKeep, ]
        nrow(MaxentKeepEvalStatsRand.df1)
        #######  Do the same for Top data
        # Eliminate duplicates for EnvVarsUsed for a given DataType
        # Select out the top wrapper test DataType
        DataTypeSel <- paste0("Top", FullNObs, "WrapperCV2Test")
        MaxentKeepEvalStatsTopPart.df <- MaxentKeepEvalStats.df2[MaxentKeepEvalStats.df2$DataType==DataTypeSel,]
        nrow(MaxentKeepEvalStatsTopPart.df)
        # Remove duplicates
        MaxentKeepEvalStatsTopPart.df2 <- MaxentKeepEvalStatsTopPart.df[!duplicated(MaxentKeepEvalStatsTopPart.df$EnvVarsUsed), ]
        head(MaxentKeepEvalStatsTopPart.df2)
        nrow(MaxentKeepEvalStatsTopPart.df2)
        # Save only VarNames
        TopVarNamesKeep <- MaxentKeepEvalStatsTopPart.df2$VarNames
        length(TopVarNamesKeep)
        ## Keep only above VarNames in all Top sets
        # First isolate all Random Sets
        MaxentKeepEvalStatsTop.df <- MaxentKeepEvalStats.df2[grepl("Top", MaxentKeepEvalStats.df2$DataType), ]
        head(MaxentKeepEvalStatsTop.df)
        tail(MaxentKeepEvalStatsTop.df)
        nrow(MaxentKeepEvalStatsTop.df)
        # Only keep original data with VarNames matching that in above data
        MaxentKeepEvalStatsTop.df1 <- MaxentKeepEvalStatsTop.df[MaxentKeepEvalStatsTop.df$VarNames %in% TopVarNamesKeep, ]
        nrow(MaxentKeepEvalStatsTop.df1)
        ########
        ### Join back together Random and Top data cleaned of duplicates
        MaxentKeepEvalStats.df2 <- rbind(MaxentKeepEvalStatsRand.df1, MaxentKeepEvalStatsTop.df1)
        nrow(MaxentKeepEvalStats.df2)
        # Sort by DataType and Rep
        MaxentKeepEvalStats.df3 <- MaxentKeepEvalStats.df2[with(MaxentKeepEvalStats.df2, order(MaxentKeepEvalStats.df2$DataType, MaxentKeepEvalStats.df2$Rep)), ]
        nrow(MaxentKeepEvalStats.df3)
        MaxentKeepEvalAllStatsL[[OutCount]] <- MaxentKeepEvalStats.df3
        #######
        ### Summarize Output
        ###Calculate Mean and Standard Deviation Values
        setwd(OutDirectpsa)
        RunType <- "WrapperCV2"
        EvaluationStats <- MaxentKeepEvalStats.df3
        ncol(EvaluationStats)
        EvaluationStats$Model <- paste0(ModelName, Rank)
        EvaluationStats <- EvaluationStats[, colnames(EvaluationStats)[c(20,1:19)]]
        #str(EvaluationStats)
        #FSAType <- Rank
        OutName <- paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets,"SetNumber_SubsetSummaryStats_", SubsetVariableNumber, "of", TotVars, "Vars", FSAType,"_", RunType, ".csv")
        SortGroups <- c("Model", "Rep", "DataType", "Set", "SetName", "SubsetVariableNumber", "TotSubsets")
        ## All statistics, and only statistics, should be at and after the column specified below
        StatVarFirstColumn <- 9
        OutDirectIn <- OutDirectpsa
        EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
        head(EvaluationStatsSummary.df)
        MaxentKeepEvalAllStatSummL[[OutCount]] <- EvaluationStatsSummary.df
        ######
        MaxentKeepEvalStatsAllSets.df <- do.call(rbind, MaxentKeepEvalAllStatSummL)
        ncol(MaxentKeepEvalStatsAllSets.df)
        ## Save data
        write.table(MaxentKeepEvalStatsAllSets.df, file=paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets,"SetNumber_", "SummaryStats_of", FinalModelVariables, "of", TotVars, "Vars", FSAType,"_", RunType, ".csv"), sep=",", col.names=NA)
        ### Summarize and save raw data as loop progresses
        MaxentKeepEvalStats.df <- do.call(rbind, MaxentKeepEvalAllStatsL)
        nrow(MaxentKeepEvalStats.df)
        head(MaxentKeepEvalStats.df)
        # Save data
        write.table(MaxentKeepEvalStats.df , file=paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets,"SetNumber_TrainVsTestStats_of", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff", ".csv"), sep=",", col.names=NA)
        #
        # Save variable permutation importance/frequency rankings
        TopVariablePermImp.df <- do.call(rbind, TopVariablePermImpL)
        head(TopVariablePermImp.df)
        nrow(TopVariablePermImp.df)
        # Save data
        write.table(TopVariablePermImp.df, file=paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets,"SetNumber_CV1ModelTrainPermImpRanksforVariablesinTopSetsof_", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff", ".csv"), sep=",", col.names=NA)
        #
        # Save variable subset permutation importance/frequency rankings
        TopMaxentPermImp.df <- do.call(rbind, TopMaxentSubsetPermImpL)
        head(TopMaxentPermImp.df)
        nrow(TopMaxentPermImp.df)
        # Save data
        write.table(TopMaxentPermImp.df, file=paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets,"SetNumber_CV1ModelTrainPermImpRanksforTopSetsof_", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff", ".csv"), sep=",", col.names=NA)
       }
    }
  }
}
#
####################################
t2 <- Sys.time()
####################################
difftime(t2,t1, units = "mins")

############################################
## Process and analyze data from above run
###########################################
for(Rank in RankList) {
  #Rank=RankList[2]
  OutDirectpsa <- paste0(Direct2, output1, "/", EvalType, "/", Rank, "Ranked")
  setwd(OutDirectpsa)
  # Read data back in
  #ModelName <- paste0("MaxentAll_RSFSA", Rank)
  #NumberModSets<- 3000
  #ModelName <- gsub("Env", "", ModelName)
  MaxentKeepEvalStats.df1 <- data.frame(read.csv(paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets,"SetNumber_TrainVsTestStats_of", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
  #str(MaxentKeepEvalStats.df1)
  # Convert factors to character
  Factors <- sapply(MaxentKeepEvalStats.df1, is.factor)
  MaxentKeepEvalStats.df1[Factors] <- lapply(MaxentKeepEvalStats.df1[Factors], as.character)
  nrow(MaxentKeepEvalStats.df1)
  ncol(MaxentKeepEvalStats.df1)
  head(MaxentKeepEvalStats.df1)
  tail(MaxentKeepEvalStats.df1)
  # Read in results for Random data from "AUC" run if not "AUC" rank
  if(Rank!="AUC") {
    OutDirectpsa_AUC <- gsub(Rank, "AUC", OutDirectpsa)
    setwd(OutDirectpsa_AUC)
    # Read data back in 
    MaxentKeepEvalStats.df1a <- data.frame(read.csv(paste0(Species, ModelName, "AUC", FullNObs, "_", NumberModSets,"SetNumber_TrainVsTestStats_of", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
    tail(MaxentKeepEvalStats.df1a)
    # Keep only random data
    MaxentKeepEvalStats.df1b <- MaxentKeepEvalStats.df1a[grepl("Random", MaxentKeepEvalStats.df1a$DataType), ]
    tail(MaxentKeepEvalStats.df1b)
    # Merge with ranked data
    MaxentKeepEvalStats.df1 <- rbind(MaxentKeepEvalStats.df1, MaxentKeepEvalStats.df1b)
    #
    setwd(OutDirectpsa)
    write.table(MaxentKeepEvalStats.df1, file=paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets,"SetNumber_TrainVsTestStats_of", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), sep=",", col.names=NA)
  }
  setwd(OutDirectpsa)
  #
  #EvaluationStats <- MaxentKeepEvalStats.df1
  #   write.table(EvaluationStatsSummary.df, file=paste0(Species, ModelName, "_", NumberModSets,"SetNumber_", "SummaryStats_of", FinalModelVariables, "of", TotVars, "Vars", FSAType,"_", RunType, ".csv"), sep=",", col.names=NA)
  
  ######################################################################################
  ## Assemble data for bar charts and Welch t test
  ##########################################################################################
  setwd(OutDirectpsa)
  # Read in results for first FSA run to get output on the leave one out variable run with all variables
  #FSAType <- Rank
  #NumberModSets <- 10000
  #
  # Keep only DataTypes including string "WrapperCV2Test"
  MaxentKeepEvalStats.df2 <- MaxentKeepEvalStats.df1[grep("WrapperCV2Test", MaxentKeepEvalStats.df1$DataType), ]
  MaxentKeepEvalStats.df2[1:50,]
  head(MaxentKeepEvalStats.df2)
  tail(MaxentKeepEvalStats.df2)
  #
  ############################################################
  ### Keep only first FullNobs2 values of 250 for each category
  #FullNObs2 <- 10
  #FullNObs <- FullNObs2
  ## First separate Top and Random data by Rep
  #MaxentKeepEvalStats.lst <- list()
  #for(Rep in RepList) {
  #  #Rep = 1
  #  TopEvalStats.df <- MaxentKeepEvalStats.df2[which(grepl("Top", MaxentKeepEvalStats.df2$Set) & MaxentKeepEvalStats.df2$Rep==Rep), ]
  #  KeepTopEvalStats.df <- TopEvalStats.df[1:FullNObs2,]
  #  KeepTopEvalStats.df$Set <- paste0("Top", FullNObs2)
  #  RandEvalStats.df <- MaxentKeepEvalStats.df2[which(grepl("Rand", MaxentKeepEvalStats.df2$Set) & MaxentKeepEvalStats.df2$Rep==Rep), ]
  #  KeepRandEvalStats.df <- RandEvalStats.df[1:FullNObs2,]
  #  KeepRandEvalStats.df$Set <- paste0("Random", FullNObs2)
  #  # Join data back together
  #  MaxentKeepEvalStats.lst[[Rep]] <- rbind(KeepTopEvalStats.df, KeepRandEvalStats.df)
  #}
  #MaxentKeepEvalStats.df2 <- do.call(rbind, MaxentKeepEvalStats.lst)
  ###############################Calculate Mean and Standard Deviation Values
  #
  ###Calculate Mean and Standard Deviation Values
  RunType <- "WrapperCV2"
  EvaluationStats <- MaxentKeepEvalStats.df2
  OutName <- paste0(Species, Rank, "Maxent_", NumberModSets, "Setsof_", SubsetVariableNumber, "_", FullNObs, "N_StatSummary.csv")
  SortGroups <- c("Rep", "DataType", "TotSubsets", "SubsetVariableNumber")
  ## All statistics, and only statistics, should be at and after the column specified below
  StatVarFirstColumn <- 6
  OutDirectIn <- OutDirectpsa
  EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
  head(EvaluationStatsSummary.df)
  ##
  EvalStatsOut <- EvaluationStatsSummary.df
  ## Sort data by DataType for ease of plotting in excel
  EvalStatsOuts <- arrange(EvalStatsOut, DataType, Rep, Statistic, DataType)
  EvalStatsOuts <- arrange(EvalStatsOuts, Statistic)
  head(EvalStatsOuts)
  write.table(EvalStatsOuts, file=paste0(Species, Rank, "Maxent_", SetName, "_", FullNObs, "Stat", SubsetVariableNumber, "VariableRunsby", NumberModSets, "SetsSummaryTable.csv"), sep=",", col.names=NA)
  ##
  #
  ###
  #############################################
  ### Conduct statistical tests for Rep
  ##################################################################################
  FullNObs2 <- FullNObs
  # Loop through various Welch t test comparisons for each ranking scheme and the three random sets
  StatTypeNames <- c("AUC", "AICc_bg", "AUCdiff")
  RankTypeNames <- c(paste0("Top", FullNObs2))
  RandomTypeNames <- c(paste0("Random", FullNObs2))
  DataType <- "WrapperCV2Test"
  RepList <- seq(1:ValidationReps)
  ##
  for(Rep in RepList) {
    # Keep only data for Rep
    #Rep=2
    RepNum <- Rep
    MaxentKeepEvalStats.df3 <- subset(MaxentKeepEvalStats.df2, Rep==RepNum)
    for(StatType in StatTypeNames) {
      #StatType <- "AUC"
      for(RankType in RankTypeNames) {
        #RankType <- paste0("Top", FullNObs2)
        MaxentEvalStatsL <- list()
        Count <- 0
        for(RandomType in RandomTypeNames) {
          #RandomType <- paste0("Random", FullNObs2)
          ## Keep only DataType with RankType  or RandomType
          MaxentTestStats <- subset(MaxentKeepEvalStats.df3, Set %in% c(RankType, RandomType))
          Count <- Count + 1
          if(StatType==StatTypeNames[1]) {
            # Identify set as Top or Random
            MaxentTestStats$SetType <- "Top"
            head(MaxentTestStats)
            ncol(MaxentTestStats)
            MaxentTestStats <- MaxentTestStats[,c(1:6, 20, 7:19)]
            MaxentTestStats <- within(MaxentTestStats, SetType[grepl("Rand", DataType)] <- 'Random')
            # Perform summary stats for excel plotting
            ###Calculate Mean and Standard Deviation Values
            RunType <- "WrapperCV2"
            EvaluationStats <- MaxentTestStats
            OutName <- paste0(Species, "_", NumberModSets, "_", SetNameF, "_", "Setsof", SubsetVariableNumber, "Variables", RankType, "Vs", RandomType, "_", StatType, FullNObs2, "N_StatSummary", "_Rep", Rep, ".csv")
            SortGroups <- c("DataType", "SetName", "SetType", "TotSubsets", "SubsetVariableNumber")
            ## All statistics, and only statistics, should be at and after the column specified below
            StatVarFirstColumn <- 8
            OutDirectIn <- OutDirectpsa
            EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
            head(EvaluationStatsSummary.df)
            ##
            EvalStatsOut <- EvaluationStatsSummary.df
            ## Sort data by DataType for ease of plotting in excel
            EvalStatsOuts <- arrange(EvalStatsOut, DataType, Statistic, DataType)
            EvalStatsOuts <- arrange(EvalStatsOut, Statistic)
            head(EvalStatsOuts)
            #
            MaxentEvalStatsL[[Count]] <- EvalStatsOuts
          }
          #
          setwd(OutDirectpsa)
          ## Use Welch correction t test
          RankWelchttest <- t.test(as.formula(paste0(StatType, " ~ DataType")), data=MaxentTestStats)
          #str(RankWelchttest)
          ## Apply Bonferroni Correction by mulitplying the p-value by the number of t-tests per StatType
          PValueBonfCor <- length(RandomTypeNames)*RankWelchttest$p.value
          # Save above Welch t test output
          out1 <- paste0(RankType, " vs ", RandomType, ": Welch t test for ", StatType, " from ", NumberModSets, "_", SetNameF, " Sets of ", SubsetVariableNumber, " Variables")
          cat(out1, file=paste0(Species, "WelchttestStatsby", FullNObs2, "from", NumberModSets,"_", SetNameF, "Setsof", SubsetVariableNumber, "Variables", RankType, "Vs", RandomType, "_", StatType, "_Rep", Rep, ".txt"), sep="\n", append=TRUE)
          out2 <- capture.output(print(RankWelchttest))
          cat(out2, file=paste0(Species, "WelchttestStatsby", FullNObs2, "from", NumberModSets, "_", SetNameF, "Setsof", SubsetVariableNumber, "Variables", RankType, "Vs", RandomType, "_", StatType, "_Rep", Rep, ".txt"), sep="\n", append=TRUE)
          out3 <- paste0("p-value with Bonferroni Correction for ", length(RandomTypeNames), " Welch t tests: ", PValueBonfCor)
          cat(out3, file=paste0(Species, "WelchttestStatsby", FullNObs2, "from", NumberModSets, "_", SetNameF, "Setsof", SubsetVariableNumber, "Variables", RankType, "Vs", RandomType, "_", StatType, "_Rep", Rep, ".txt"), sep="\n", append=TRUE)
        }
        if(StatType==StatTypeNames[1]) {
          MaxentEvalStatsAll.df <- do.call(rbind, MaxentEvalStatsL)
          ## Sort data
          MaxentEvalStatsAll.dfs <- arrange(MaxentEvalStatsAll.df, SetType, Statistic, SetType)
          write.table(MaxentEvalStatsAll.dfs, file=paste0(Species, "WelchttestStatsby", FullNObs2, "from", NumberModSets,"_", SetNameF, "Setsof", SubsetVariableNumber, "Vars", RankType, "vsRandom_DataTable", "_Rep", Rep, ".csv"), sep=",", col.names=NA)
      }
      }
    }
  }
  ##########
  ##
  ######################################################################################
  ## Assemble data for bar charts and Welch t test for top variable sets ranked
  ## by mean value of permutation importance/frequency ranking of variables
  ##########################################################################################
  setwd(OutDirectpsa)
  # Read in results for first FSA run to get output on the leave one out variable run with all variables
  #FSAType <- Rank
  #NumberModSets <- 10000
  # FullNObs3 is how many top versus bottom ranked of the selected 250 subsets to test below
  #
  # Keep only DataTypes including string "WrapperCV2Test"
  MaxentKeepEvalStats.df2 <- MaxentKeepEvalStats.df1[grep("WrapperCV2Test", MaxentKeepEvalStats.df1$DataType), ]
  #str(MaxentKeepEvalStats.df2)
  # Keep only DataTypes including string "Top"
  MaxentKeepEvalStats.df3 <- MaxentKeepEvalStats.df2[grep("Top", MaxentKeepEvalStats.df2$DataType), ]
  MaxentKeepEvalStats.df3[1:50,]
  nrow(MaxentKeepEvalStats.df3)
  head(MaxentKeepEvalStats.df3)
  tail(MaxentKeepEvalStats.df3)
  #
  ## Merge data with CV1 variable set ranks by mean permutation importance/frequency ranks
  # Read in ranks for subsets based on permutation importance/frequency
  # Read data back in
  #NumberModSets <- 3000
  TopVariablePermImp.df1 <- data.frame(read.csv(paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets,"SetNumber_CV1ModelTrainPermImpRanksforTopSetsof_", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff", ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
  #str(TopVariablePermImp.df1)
  # Convert VarNames from factor to character
  TopVariablePermImp.df1$VarNames <- sapply(TopVariablePermImp.df1$VarNames, as.character)
  nrow(TopVariablePermImp.df1)
  ncol(TopVariablePermImp.df1)
  head(TopVariablePermImp.df1)
  tail(TopVariablePermImp.df1)
  #
  ## Match Ranking to subset by Rep
  TopMaxentKeepEvalRanksL <- list()
  BottomMaxentKeepEvalRanksL <- list()
  RandomMaxentKeepEvalRanksL <- list()
  for(i in 1:ValidationReps) {
    #i <- 1
    Rep <- i
    MaxentKeepEvalStats.df4 <- MaxentKeepEvalStats.df3[which(MaxentKeepEvalStats.df3$Rep==Rep), ]
    nrow(MaxentKeepEvalStats.df4)
    head(MaxentKeepEvalStats.df4)
    #str(MaxentKeepEvalStats.df4)
    TopVariablePermImp.df2 <- TopVariablePermImp.df1[which(TopVariablePermImp.df1$Rep==Rep), ]
    TopVariablePermImp.df2 <- TopVariablePermImp.df2[order(TopVariablePermImp.df2$VarNames),]
    #str(TopVariablePermImp.df2)
    head(TopVariablePermImp.df2)
    nrow(TopVariablePermImp.df2)
    TopVariablePermImp.df2[1:50,]
    # Add ranks to evaluation stats
    MaxentKeepEvalStats.df4$CV1MnPermImpRank <- TopVariablePermImp.df2$PermImp_MMOORA_Rank[match(MaxentKeepEvalStats.df4$VarNames, TopVariablePermImp.df2$VarNames)]
    #MaxentKeepEvalStats.df4 <- na.omit(MaxentKeepEvalStats.df4)
    head(MaxentKeepEvalStats.df4)
    nrow(MaxentKeepEvalStats.df4)
    #MissingRank.df <- MaxentKeepEvalStats.df4[rowSums(is.na(MaxentKeepEvalStats.df4))>0,]
    #nrow(MissingRank.df)
    # Delete the one or so instances where the rank may be missing
    MaxentKeepEvalStats.df4 <- na.omit(MaxentKeepEvalStats.df4)
    any(is.na(MaxentKeepEvalStats.df4))
    # Sort by rank
    MaxentKeepEvalStats.df4s <- MaxentKeepEvalStats.df4[order(MaxentKeepEvalStats.df4$CV1MnPermImpRank),]
    head(MaxentKeepEvalStats.df4s)
    tail(MaxentKeepEvalStats.df4s)
    # Identify the top 20% and bottom 20% ranked  subsets
    NumRows <- nrow(MaxentKeepEvalStats.df4s)
    # Number rows sequentially
    rownames(MaxentKeepEvalStats.df4s) <- seq(1, NumRows, 1)
    GroupRange <- min(Round2(NumRows*(FullNObs3/FullNObs),0), FullNObs3)
    FullNObs3 <- GroupRange
    TopMaxentKeepEvalStats.df4s <- MaxentKeepEvalStats.df4s[1:GroupRange,]
    head(TopMaxentKeepEvalStats.df4s)
    tail(TopMaxentKeepEvalStats.df4s)
    nrow(TopMaxentKeepEvalStats.df4s)
    BottomMaxentKeepEvalStats.df4s <- MaxentKeepEvalStats.df4s[(NumRows-GroupRange+1):NumRows,]
    tail(BottomMaxentKeepEvalStats.df4s)
    nrow(BottomMaxentKeepEvalStats.df4s)
    RandomMaxentKeepEvalStats.df4 <- MaxentKeepEvalStats.df4s[sample(nrow(MaxentKeepEvalStats.df4s),FullNObs3),]
    TopMaxentKeepEvalRanksL [[i]] <- TopMaxentKeepEvalStats.df4s
    BottomMaxentKeepEvalRanksL [[i]] <- BottomMaxentKeepEvalStats.df4s
    RandomMaxentKeepEvalRanksL [[i]] <- RandomMaxentKeepEvalStats.df4
  }
  #
  TopMaxentKeepEvalRanks.df <- do.call(rbind, TopMaxentKeepEvalRanksL)
  TopMaxentKeepEvalRanks.df$Group <- paste0("Top", FullNObs3)
  head(TopMaxentKeepEvalRanks.df)
  nrow(TopMaxentKeepEvalRanks.df)
  BottomMaxentKeepEvalRanks.df <- do.call(rbind, BottomMaxentKeepEvalRanksL)
  BottomMaxentKeepEvalRanks.df$Group <- paste0("Bottom", FullNObs3)
  head(BottomMaxentKeepEvalRanks.df)
  RandomMaxentKeepEvalRanks.df <- do.call(rbind, RandomMaxentKeepEvalRanksL)
  RandomMaxentKeepEvalRanks.df$Group <- paste0("Random", FullNObs3)
  head(RandomMaxentKeepEvalRanks.df)
  #
  MaxentKeepEvalRanks.df <- rbind(TopMaxentKeepEvalRanks.df, BottomMaxentKeepEvalRanks.df, RandomMaxentKeepEvalRanks.df)
  head(MaxentKeepEvalRanks.df)
  tail(MaxentKeepEvalRanks.df)
  ncol(MaxentKeepEvalRanks.df)
  # Move Group variable
  MaxentKeepEvalRanks.df2 <- MaxentKeepEvalRanks.df[,c(1:7,21,8:20)]
  head(MaxentKeepEvalRanks.df2)
  tail(MaxentKeepEvalRanks.df2)
  # Save data
  write.table(MaxentKeepEvalRanks.df2, file=paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets,"SetNumber_CV1ModelStatsRanksforSelected", FullNObs3, "TopVsBottomSetsof_", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff", ".csv"), sep=",", col.names=NA)
  
  ############################################################
  ### Keep only first FullNobs2 values of 250 for each category
  #FullNObs2 <- 10
  #FullNObs <- FullNObs2
  ## First separate Top and Random data by Rep
  #MaxentKeepEvalStats.lst <- list()
  #for(Rep in RepList) {
  #  #Rep = 1
  #  TopEvalStats.df <- MaxentKeepEvalStats.df2[which(grepl("Top", MaxentKeepEvalStats.df2$Set) & MaxentKeepEvalStats.df2$Rep==Rep), ]
  #  KeepTopEvalStats.df <- TopEvalStats.df[1:FullNObs2,]
  #  KeepTopEvalStats.df$Set <- paste0("Top", FullNObs2)
  #  RandEvalStats.df <- MaxentKeepEvalStats.df2[which(grepl("Rand", MaxentKeepEvalStats.df2$Set) & MaxentKeepEvalStats.df2$Rep==Rep), ]
  #  KeepRandEvalStats.df <- RandEvalStats.df[1:FullNObs2,]
  #  KeepRandEvalStats.df$Set <- paste0("Random", FullNObs2)
  #  # Join data back together
  #  MaxentKeepEvalStats.lst[[Rep]] <- rbind(KeepTopEvalStats.df, KeepRandEvalStats.df)
  #}
  #MaxentKeepEvalStats.df2 <- do.call(rbind, MaxentKeepEvalStats.lst)
  ###############################Calculate Mean and Standard Deviation Values
  #
  ###Calculate Mean and Standard Deviation Values
  RunType <- "WrapperCV2"
  EvaluationStats <- MaxentKeepEvalRanks.df2
  OutName <- paste0(Species, Rank, "Maxent_", NumberModSets, "Setsof_", SubsetVariableNumber, "_", FullNObs3, "N_SelectedTopVsBottomRank_StatSummary.csv")
  SortGroups <- c("Rep", "Group", "DataType", "TotSubsets", "SubsetVariableNumber")
  ## All statistics, and only statistics, should be at and after the column specified below
  StatVarFirstColumn <- 7
  OutDirectIn <- OutDirectpsa
  EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
  head(EvaluationStatsSummary.df)
  ##
  EvalStatsOut <- EvaluationStatsSummary.df
  ## Sort data by DataType for ease of plotting in excel
  EvalStatsOuts <- arrange(EvalStatsOut, Group, Rep, Statistic, Group)
  EvalStatsOuts <- arrange(EvalStatsOuts, Statistic)
  head(EvalStatsOuts)
  write.table(EvalStatsOuts, file=paste0(Species, Rank, "Maxent_", SetName, "_", FullNObs, "Stat", SubsetVariableNumber, "VariableRunsby", FullNObs3, "SelectedTopVsBottomRank_SetsSummaryTable.csv"), sep=",", col.names=NA)
  ##
  #
  ###
  #
  #################################################################################
  ### Conduct statistical tests for Rep Comparing Top Ranked Vs Bottom Ranked
  ##################################################################################
  # Loop through various Welch t test comparisons for each ranking scheme and the three random sets
  #FullNObs3 <- 50
  StatTypeNames <- c("AUC", "AICc_bg", "AUCdiff")
  TopTypeNames <- c(paste0("Top", FullNObs3))
  BottomTypeNames <- c(paste0("Bottom", FullNObs3))
  DataType <- "WrapperCV2Test"
  # Delete Random sets from analysis
  nrow(MaxentKeepEvalRanks.df2)
  MaxentKeepEvalRanks.df2a <- MaxentKeepEvalRanks.df2[!grepl("Random",MaxentKeepEvalRanks.df2$Group),]
  nrow(MaxentKeepEvalRanks.df2a)
  ##
  for(Rep in RepList) {
    # Keep only data for Rep
    #Rep=1
    RepNum <- Rep
    MaxentKeepEvalRanks.df3 <- subset(MaxentKeepEvalRanks.df2a, Rep==RepNum)
    head(MaxentKeepEvalRanks.df3)
    for(StatType in StatTypeNames) {
      #StatType <- "AUC"
      for(TopType in TopTypeNames) {
        #RankType <- paste0("Top", FullNObs3)
        MaxentEvalStatsL <- list()
        Count <- 0
        for(BottomType in BottomTypeNames) {
          #BottomType <- paste0("Bottom", FullNObs3)
          ## Keep only DataType with RankType  or BottomType
          MaxentTestStats <- MaxentKeepEvalRanks.df3
          Count <- Count + 1
          if(StatType==StatTypeNames[1]) {
            # Perform summary stats for excel plotting
            ###Calculate Mean and Standard Deviation Values
            RunType <- "WrapperCV2"
            EvaluationStats <- MaxentTestStats
            OutName <- paste0(Species, "_", NumberModSets, "_", SetNameF, "_", "Setsof", SubsetVariableNumber, "Variables", TopType, "Vs", BottomType, "_", StatType, FullNObs3, "N_StatSummary", "_Rep", Rep, ".csv")
            SortGroups <- c("Group", "Set", "SetName", "TotSubsets", "SubsetVariableNumber")
            ## All statistics, and only statistics, should be at and after the column specified below
            StatVarFirstColumn <- 7
            OutDirectIn <- OutDirectpsa
            EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
            head(EvaluationStatsSummary.df)
            ##
            EvalStatsOut <- EvaluationStatsSummary.df
            ## Sort data by DataType for ease of plotting in excel
              EvalStatsOuts <- arrange(EvalStatsOut, Statistic)
         EvalStatsOuts <- arrange(EvalStatsOut, Group, Statistic, Set)
             head(EvalStatsOuts)
            #
            MaxentEvalStatsL[[Count]] <- EvalStatsOuts
          }
          #
          setwd(OutDirectpsa)
          ## Use Welch correction t test
          RankWelchttest <- t.test(as.formula(paste0(StatType, " ~ Group")), data=MaxentTestStats)
          #str(RankWelchttest)
          ## Apply Bonferroni Correction by mulitplying the p-value by the number of t-tests per StatType
          PValueBonfCor <- length(BottomTypeNames)*RankWelchttest$p.value
          # Save above Welch t test output
          out1 <- paste0(TopType, " vs ", BottomType, ": Welch t test for ", StatType, " from ", NumberModSets, "_", SetNameF, " Sets of ", SubsetVariableNumber, " Variables")
          cat(out1, file=paste0(Species, "WelchttestStatsby", FullNObs3, "from", FullNObs,"Selected_", SetNameF, "Setsof", SubsetVariableNumber, "Variables", TopType, "Vs", BottomType, "_", StatType, "_Rep", Rep, ".txt"), sep="\n", append=TRUE)
          out2 <- capture.output(print(RankWelchttest))
          cat(out2, file=paste0(Species, "WelchttestStatsby", FullNObs3, "from", FullNObs,"Selected_", SetNameF, "Setsof", SubsetVariableNumber, "Variables", TopType, "Vs", BottomType, "_", StatType, "_Rep", Rep, ".txt"), sep="\n", append=TRUE)
          out3 <- paste0("p-value with Bonferroni Correction for ", length(BottomTypeNames), " Welch t tests: ", PValueBonfCor)
          cat(out3, file=paste0(Species, "WelchttestStatsby", FullNObs3, "from", FullNObs,"Selected_", SetNameF, "Setsof", SubsetVariableNumber, "Variables", TopType, "Vs", BottomType, "_", StatType, "_Rep", Rep, ".txt"), sep="\n", append=TRUE)
        }
        if(StatType==StatTypeNames[1]) {
          MaxentEvalStatsAll.df <- do.call(rbind, MaxentEvalStatsL)
          ## Sort data
          MaxentEvalStatsAll.dfs <- arrange(MaxentEvalStatsAll.df, Group, Statistic, Group)
          write.table(MaxentEvalStatsAll.dfs, file=paste0(Species, "WelchttestStatsby", FullNObs3, "from", FullNObs,"Selected_", SetNameF, "Setsof", SubsetVariableNumber, "Vars", TopType, "vs", BottomType, "_DataTable", "_Rep", Rep, ".csv"), sep=",", col.names=NA)
      }
      }
    }
  }
  #################################################################################
  ### Conduct statistical tests for Rep Comparing Top Ranked Vs Random
  ##################################################################################
  # Loop through various Welch t test comparisons for each ranking scheme and the three random sets
  #FullNObs3 <- 50
  StatTypeNames <- c("AUC", "AICc_bg", "AUCdiff")
  TopTypeNames <- c(paste0("Top", FullNObs3))
  RandomTypeNames <- c(paste0("Random", FullNObs3))
  DataType <- "WrapperCV2Test"
  # Delete Bottom sets from analysis
  nrow(MaxentKeepEvalRanks.df2)
  MaxentKeepEvalRanks.df2a <- MaxentKeepEvalRanks.df2[!grepl("Bottom",MaxentKeepEvalRanks.df2$Group),]
  nrow(MaxentKeepEvalRanks.df2a)
  ##
  for(Rep in RepList) {
    # Keep only data for Rep
    #Rep=1
    RepNum <- Rep
    MaxentKeepEvalRanks.df3 <- subset(MaxentKeepEvalRanks.df2a, Rep==RepNum)
    head(MaxentKeepEvalRanks.df3)
    for(StatType in StatTypeNames) {
      #StatType <- "AUC"
      for(TopType in TopTypeNames) {
        #RankType <- paste0("Top", FullNObs3)
        MaxentEvalStatsL <- list()
        Count <- 0
        for(RandomType in RandomTypeNames) {
          #RandomType <- paste0("Random", FullNObs3)
          ## Keep only DataType with RankType  or RandomType
          MaxentTestStats <- MaxentKeepEvalRanks.df3
          Count <- Count + 1
          if(StatType==StatTypeNames[1]) {
            # Perform summary stats for excel plotting
            ###Calculate Mean and Standard Deviation Values
            RunType <- "WrapperCV2"
            EvaluationStats <- MaxentTestStats
            OutName <- paste0(Species, "_", NumberModSets, "_", SetNameF, "_", "Setsof", SubsetVariableNumber, "Variables", TopType, "Vs", RandomType, "_", StatType, FullNObs3, "N_StatSummary", "_Rep", Rep, ".csv")
            SortGroups <- c("Group", "Set", "SetName", "TotSubsets", "SubsetVariableNumber")
            ## All statistics, and only statistics, should be at and after the column specified below
            StatVarFirstColumn <- 7
            OutDirectIn <- OutDirectpsa
            EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
            head(EvaluationStatsSummary.df)
            ##
            EvalStatsOut <- EvaluationStatsSummary.df
            ## Sort data by DataType for ease of plotting in excel
              EvalStatsOuts <- arrange(EvalStatsOut, Statistic)
         EvalStatsOuts <- arrange(EvalStatsOut, Group, Statistic, Set)
             head(EvalStatsOuts)
            #
            MaxentEvalStatsL[[Count]] <- EvalStatsOuts
          }
          #
          setwd(OutDirectpsa)
          ## Use Welch correction t test
          RankWelchttest <- t.test(as.formula(paste0(StatType, " ~ Group")), data=MaxentTestStats)
          #str(RankWelchttest)
          ## Apply Bonferroni Correction by mulitplying the p-value by the number of t-tests per StatType
          PValueBonfCor <- length(RandomTypeNames)*RankWelchttest$p.value
          # Save above Welch t test output
          out1 <- paste0(TopType, " vs ", RandomType, ": Welch t test for ", StatType, " from ", NumberModSets, "_", SetNameF, " Sets of ", SubsetVariableNumber, " Variables")
          cat(out1, file=paste0(Species, "WelchttestStatsby", FullNObs3, "from", FullNObs,"Selected_", SetNameF, "Setsof", SubsetVariableNumber, "Variables", TopType, "Vs", RandomType, "_", StatType, "_Rep", Rep, ".txt"), sep="\n", append=TRUE)
          out2 <- capture.output(print(RankWelchttest))
          cat(out2, file=paste0(Species, "WelchttestStatsby", FullNObs3, "from", FullNObs,"Selected_", SetNameF, "Setsof", SubsetVariableNumber, "Variables", TopType, "Vs", RandomType, "_", StatType, "_Rep", Rep, ".txt"), sep="\n", append=TRUE)
          out3 <- paste0("p-value with Bonferroni Correction for ", length(RandomTypeNames), " Welch t tests: ", PValueBonfCor)
          cat(out3, file=paste0(Species, "WelchttestStatsby", FullNObs3, "from", FullNObs,"Selected_", SetNameF, "Setsof", SubsetVariableNumber, "Variables", TopType, "Vs", RandomType, "_", StatType, "_Rep", Rep, ".txt"), sep="\n", append=TRUE)
        }
        if(StatType==StatTypeNames[1]) {
          MaxentEvalStatsAll.df <- do.call(rbind, MaxentEvalStatsL)
          ## Sort data
          MaxentEvalStatsAll.dfs <- arrange(MaxentEvalStatsAll.df, Group, Statistic, Group)
          write.table(MaxentEvalStatsAll.dfs, file=paste0(Species, "WelchttestStatsby", FullNObs3, "from", FullNObs,"Selected_", SetNameF, "Setsof", SubsetVariableNumber, "Vars", TopType, "vs", RandomType, "_DataTable", "_Rep", Rep, ".csv"), sep=",", col.names=NA)
      }
      }
    }
  }
  #################################################################################
  ### Compile variable rankings from all three reps, averaging values per variable
  #####################################################################################
  ##
  setwd(OutDirectpsa)
  # Read in variable rankings from all reps
  VariableRank.df <- data.frame(read.csv(paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets,"SetNumber_CV1ModelTrainPermImpRanksforVariablesinTopSetsof_", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff", ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
  head(VariableRank.df)
  #str(VariableRank.df1)
  # Convert VarNames from factor to character
  VariableRank.df$Variable <- sapply(VariableRank.df$Variable, as.character)
  nrow(VariableRank.df)
  ncol(VariableRank.df)
  head(VariableRank.df)
  tail(VariableRank.df)
  # Obtain mean rank values by variable
  VariableRank.means <- aggregate(PermImp_MMOORA_Rank ~ Variable, VariableRank.df, mean)
  # Obtain sd for each variable
  VariableRank.sd <- aggregate(PermImp_MMOORA_Rank ~ Variable, VariableRank.df, sd)
  colnames(VariableRank.sd)[2] <- "PermImp_MMOORA_Rank_SD"
  # Merge above datasets together
  VariableRank.means$PermImp_MMOORA_Rank_SD <- VariableRank.sd$PermImp_MMOORA_Rank_SD[match(VariableRank.means$Variable, VariableRank.sd$Variable)]
  nrow(VariableRank.means)
  head(VariableRank.means)
  ##
  # Obtain Mean of Mean Permutation Importance by Variable
  VariablePermImp.means <- aggregate(MeanPermImp ~ Variable, VariableRank.df, mean)
  colnames(VariablePermImp.means)[2] <- "PermImp_Mean"
  # Obtain Mean of SD Permutation Importance by Variable
  VariablePermImp.sd <- aggregate(MeanPermImp ~ Variable, VariableRank.df, sd)
  colnames(VariablePermImp.sd)[2] <- "PermImp_SD"
  # Merge above datasets together
  VariablePermImp.means$PermImp_SD <- VariablePermImp.sd$PermImp_SD[match(VariablePermImp.means$Variable, VariablePermImp.sd$Variable)]
  ##
  # Obtain Mean of Mean Permutation Importance by Variable
  VariableCount.means <- aggregate(Count ~ Variable, VariableRank.df, mean)
  colnames(VariableCount.means)[2] <- "Count_Mean"
  # Obtain SD of Mean Permutation Importance by Variable
  VariableCount.sd <- aggregate(Count ~ Variable, VariableRank.df, sd)
  colnames(VariableCount.sd)[2] <- "Count_SD"
  # Merge above datasets together
  VariableCount.means$Count_SD <- VariableCount.sd$Count_SD[match(VariableCount.means$Variable, VariableCount.sd$Variable)]
  # Merge all datasets together
  VariableRank.means2 <- cbind(VariableRank.means, VariablePermImp.means[,2:3], VariableCount.means[,2:3])
  head(VariableRank.means2)
  # Obtain N reps for each variable
  VariableRank.count <- as.data.frame(VariableRank.df %>% dplyr::count(Variable))
  colnames(VariableRank.count)[2] <- "N"
  # Merge previous datasets together
  VariableRank.means2$N <- VariableRank.count$N[match(VariableRank.means2$Variable, VariableRank.count$Variable)]
  # Sort by MeanRank
  VariableRank.means2 <- VariableRank.means2[order(VariableRank.means2$PermImp_MMOORA_Rank),]
  #VariableRank.means2$PermImp_MMOORA_Rank[VariableRank.means2$Variable=="BIO_12"]
  ##
  #######################
  ### Group variables within correlation threshold starting from topmost ranked variables
  ##
  ## Read in correlation values
  # Specify SetName
  #OutDirectpsaOrig <- gsub("AIC", "AUC", OutDirectpsa)
  #setwd(OutDirectpsaOrig)
  SetName <- paste0("01_", colnames(VariableSetTypes)[1])
  CorrelationAbsVals.df <- data.frame(read.csv(paste0(Species, "_DataSet", SetName, "_CorrelationBackgroundEnvironVariables.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
  head(CorrelationAbsVals.df)
  setwd(OutDirectpsa)
  Variables <- VariableRank.means2[,1]
  RemainingVariables <- Variables
  VariableGroups <- vector()
  VariableGroupNums <- vector()
  VariableGroupHeads <- vector()
  GroupCount <- 0
  ## First go through list and find all variables within correlation threshold of top variable and put in group
  ## and remove them from list. Then go through second variable and do the same until list is empty
  for(i in 1:length(Variables)) {
    # i=2
    # If variable within correlation threshold of top variable, put in group list
    GroupCount <- i
    Group <- vector()
    GroupNum <- vector()
    GroupHead <- vector()
    Group[1] <- RemainingVariables[1]
    if(GroupCount < 10) {
      GroupNum[1] <- paste0("0", GroupCount)
    } else {
      GroupNum[1] <- paste0(GroupCount)
    }
    GroupHead[1] <- Group[1]
    Count <- 1
    for(j in 2:length(RemainingVariables)) {
      #j=2
      if(CorrelationAbsVals.df[Group[1], RemainingVariables[j]]>=0.7) {
        Count <- Count + 1
        Group[Count]<- RemainingVariables[j]
        GroupNum[Count] <- GroupNum[1]
        GroupHead[Count] <- Group[1]
      }
    }
    # Eliminate variables in above group from remaining variables
    RemoveVariables <- Group[1:length(Group)]
    RemainingVariables <- RemainingVariables[! RemainingVariables %in% RemoveVariables]
    VariableGroups <- c(VariableGroups, Group)
    VariableGroupNums <- c(VariableGroupNums, GroupNum)
    VariableGroupHeads <- c(VariableGroupHeads, GroupHead)
    if(length(RemainingVariables)==1) {
      GroupCount <- GroupCount + 1
      Group[1] <- RemainingVariables[1]
      if(GroupCount < 10) {
        GroupNum[1] <- paste0("0", GroupCount)
      } else {
        GroupNum[1] <- paste0(GroupCount)
      }
      GroupHead[1] <- Group[1]
      VariableGroups <- c(VariableGroups, Group)
      VariableGroupNums <- c(VariableGroupNums, GroupNum)
      VariableGroupHeads <- c(VariableGroupHeads, GroupHead)
      break
    }
    if(length(RemainingVariables)==0) {
      break
    }
  }
  #######
  UniqueVariableGroupHeads <- unique(VariableGroupHeads)
  ##
  Variables <- VariableRank.means2[,1]
  VariableGroups <- vector()
  VariableGroupNums <- vector()
  VariableGroupHeads <- vector()
  VariableRanks <- vector()
  VariableRankDiffs <- vector()
  GroupCount <- 0
  ## Second, go through list and find all variables within correlation threshold of head variable and put in groups
  for(i in 1:length(UniqueVariableGroupHeads)) {
    # i=2
    # If variable within correlation threshold of top variable, put in group list
    GroupCount <- i
    Group <- vector()
    GroupNum <- vector()
    GroupHead <- vector()
    VRank <- vector()
    RankDiff <- vector()
    Group[1] <- UniqueVariableGroupHeads[i]
    if(GroupCount < 10) {
      GroupNum[1] <- paste0("0", GroupCount)
    } else {
      GroupNum[1] <- paste0(GroupCount)
    }
    GroupHead[1] <- Group[1]
    VRank[1] <- VariableRank.means2$PermImp_MMOORA_Rank[VariableRank.means2$Variable==Group[1]]
    RankDiff[1] <- 0
    Count <- 1
    for(j in 1:length(Variables)) {
      #j=2
      if(CorrelationAbsVals.df[Group[1], Variables[j]]>=0.7) {
        Count <- Count + 1
        Group[Count]<- Variables[j]
        GroupNum[Count] <- GroupNum[1]
        GroupHead[Count] <- Group[1]
        VRank[Count] <- VariableRank.means2$PermImp_MMOORA_Rank[VariableRank.means2$Variable==Group[Count]]
        RankDiff[Count] <- VRank[Count] - VRank[1]
      }
    }
    VariableGroups <- c(VariableGroups, Group)
    VariableGroupNums <- c(VariableGroupNums, GroupNum)
    VariableGroupHeads <- c(VariableGroupHeads, GroupHead)
    VariableRanks <- c(VariableRanks, VRank)
    VariableRankDiffs <- c(VariableRankDiffs, RankDiff)
  }
  #######
  VariableGroups.df <- data.frame(VariableGroups, stringsAsFactors=FALSE)
  colnames(VariableGroups.df)[1] <- "Variable"
  VariableGroupNums.df <- data.frame(VariableGroupNums, stringsAsFactors=FALSE)
  colnames(VariableGroupNums.df)[1] <- "VariableGroup"
  VariableGroupHeads.df <- data.frame(VariableGroupHeads, stringsAsFactors=FALSE)
  colnames(VariableGroupHeads.df)[1] <- "VariableGroupHeads"
  VariableRankDiffs.df <- data.frame(VariableRankDiffs, stringsAsFactors=FALSE)
  colnames(VariableRankDiffs.df)[1] <- "VariableRankDiff"
  VariableGroups2.df <- cbind(VariableGroups.df, VariableGroupNums.df, VariableGroupHeads.df, VariableRankDiffs.df)
  # Remove duplicates
  VariableGroups2.df <- unique(VariableGroups2.df)
  # Remove rows with negative values for VariablRankDiff
  VariableGroups2.df <- VariableGroups2.df[!(VariableGroups2.df$VariableRankDiff<0),]
  # For multiple values, keep one with lowest RankDiff
  # Loop through all variables, only keep row with lowest RankDiff
  VariableGroups3.df <- VariableGroups2.df[1,]
  for(i in 1:length(Variables)) {
    # i=1
    # Keep only rows with variable i
    CheckVariable.df <- VariableGroups2.df[(VariableGroups2.df$Variable==Variables[i]),]
    # Keep only row with lowest RankDiff
    VariableGroups3.df[i,] <- CheckVariable.df[which.min(CheckVariable.df$VariableRankDiff),]
  }
  nrow(VariableGroups3.df)
  ##
  # Join with other data
  VariableRank.means2$VariableGroups <- VariableGroups3.df$VariableGroup[match(VariableRank.means2$Variable, VariableGroups3.df$Variable)]
  head(VariableRank.means2)
  # Sort by VariableGroups
  VariableRank.means2 <- VariableRank.means2[order(VariableRank.means2$VariableGroups),]
  # Convert variable names to lower case
  VariableRank.means2$Variable <- tolower(VariableRank.means2$Variable)
  # Save data
  write.table(VariableRank.means2, file=paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets,"SetNumber_TrainOverallMeanPermImpRanksforVarsinTopSetsof_", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff", ".csv"), sep=",", col.names=NA)
}

#####################################################################################################
### Loop through  top 12 selected subsets to generate maxent projection maps
#####################################################################################################
## This loop takes about 80 minutes for 10 variable subsets of selected 10 and random 10 models
## This loop takes 2.2 hours for 23 models of 12 variable subsets
## This loop takes 4.1 hours for 45 models of 12 variable subsets
## This loop takes about 9.3 hours for 30 models of 15 variables subsets for fire
## This loop takes about 5.31 hours for 4 models of 6 variables subsets for fire
## This loop takes about 9.95 hours for 12 models of 8 variables for Zizotes milkweed
## This loop takes about 9.4 hours for 12 models of 6 variables for broadleaf milkweed
#
# Specify final selected ranking system
Rank <- "AIC"
OutDirectpsa <- paste0(Direct2, output1, "/", EvalType, "/", Rank, "Ranked")
#
FinalModelVariables <- 6
#FinalModelVariables <- 3
NObsJ <- 3
## Determine NumberModSets used in Stage II RSFSA
if(FinalModelVariables==2) {
  if(ncol(combn(TotVars,2))>3000) {
    NumberModSets <- 3000
    NumberModSetsIn <- 3150
    FullNObs <- 250
    PlusNum <- 100 # extra random subsets
    #FullNObs <- NumberModSets/3
    #PlusNum <- 1
  } else {
    NumberModSets <- ncol(combn(TotVars,2))
    NumberModSetsIn <- NumberModSets
    FullNObs <- 0.05*NumberModSets
    PlusNum <- 0.025*NumberModSets # extra random subsets
    #FullNObs <- NumberModSets/3
    #PlusNum <- 1
  }
} else {
  NumberModSetsIn <- 3150 # Number sets per rep: extra sets if needed due to AICc values deleted- recommend no more than 5050 sets at once in case of crash
  NumberModSets <- 3000  # Number sets per rep
  #NumberModSetsIn <- 33 # extra sets if needed due to AICc values deleted
  #NumberModSets <- 30
  FullNObs <- 250
  PlusNum <- 100 # extra random subsets
  #FullNObs <- NumberModSets/3
  #PlusNum <- 1
}
#
#NumProjections <- 9  # Specify number of models to project
#NumberModSets <- 30
NumProjections <- 12
## Read in clipped current climate and other grids
setwd(InDirect)
Predictors <- stack(GridNamesL, proj4string=crs.out)
names(Predictors) <- toupper(names(Predictors))
#names(Predictors) <- gsub("_NS", "", names(Predictors))
#
#############################
j = 1
#j=1
if(j < 10) {
  SetNumber <- paste0("0", j)
} else {
  SetNumber <- as.character(j)
}
#
SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes[j]))
SetNameF <- colnames(VariableSetTypes[j])
setwd(OutDirectpsa)
#list.files()
# Read in results for first FSA run to get output on the leave one out variable run with all variables
#FSAType <- Rank
#NumberModSets <- 33
#SubsetSizes <- c(3,8,10,12,15,20,25)
#
#################################################################################
### Assemble top 15 ranked variable subsets by AUCWrappertest from each RSFSA rep
MaxentKeepEvalStats.df1 <- data.frame(read.csv(paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets,"SetNumber_TrainVsTestStats_of", FinalModelVariables, "of", TotVars, "NumVars_", FSAType, PsAbsBuff, Units, "PsAbsBuff_", SpatFiltBuff, Units, "SpatFiltBuff.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
nrow(MaxentKeepEvalStats.df1)
ncol(MaxentKeepEvalStats.df1)
head(MaxentKeepEvalStats.df1)
tail(MaxentKeepEvalStats.df1)
#################
### Keep only first NumProjections/3 variable subsets from each Rep
## NOTE: DO NOT RESORT BY THIS AUCwrapperCV2test- THIS WILL INFLATE AUC- AUCwrapperCV1test already used to sort
#################
#
## Keep only FinalTest data
MaxentKeepEvalStats.dfa <- MaxentKeepEvalStats.df1[grep(paste0("WrapperCV2Test"), MaxentKeepEvalStats.df1$DataType), ]
nrow(MaxentKeepEvalStats.dfa)
head(MaxentKeepEvalStats.dfa)
tail(MaxentKeepEvalStats.dfa)
##################################
## Eliminate any duplicate EnvVarsUsed sets among each data type for all reps
# First randomize sorting by rep to even out removals
MaxentKeepEvalStats.df2 <- MaxentKeepEvalStats.dfa[sample(nrow(MaxentKeepEvalStats.dfa)), ]
head(MaxentKeepEvalStats.df2)
tail(MaxentKeepEvalStats.df2)
nrow(MaxentKeepEvalStats.df2)
# Select out the random wrapper test DataType
DataTypeSel <- paste0("Random", FullNObs, "WrapperCV2Test")
MaxentKeepEvalStatsRandPart.df <- MaxentKeepEvalStats.df2[MaxentKeepEvalStats.df2$DataType==DataTypeSel,]
nrow(MaxentKeepEvalStatsRandPart.df)
# Remove duplicates
MaxentKeepEvalStatsRandPart.df2 <- MaxentKeepEvalStatsRandPart.df[!duplicated(MaxentKeepEvalStatsRandPart.df$EnvVarsUsed), ]
head(MaxentKeepEvalStatsRandPart.df2)
nrow(MaxentKeepEvalStatsRandPart.df2)
# Save only VarNames
RandVarNamesKeep <- MaxentKeepEvalStatsRandPart.df2$VarNames
length(RandVarNamesKeep)
## Keep only above VarNames in all Random sets
# First isolate all Random Sets
MaxentKeepEvalStatsRand.df <- MaxentKeepEvalStats.df2[grepl("Rand", MaxentKeepEvalStats.df2$DataType), ]
head(MaxentKeepEvalStatsRand.df)
tail(MaxentKeepEvalStatsRand.df)
nrow(MaxentKeepEvalStatsRand.df)
# Only keep original data with VarNames matching that in above data
MaxentKeepEvalStatsRand.df1 <- MaxentKeepEvalStatsRand.df[MaxentKeepEvalStatsRand.df$VarNames %in% RandVarNamesKeep, ]
nrow(MaxentKeepEvalStatsRand.df1)
#######  Do the same for Top data
# Eliminate duplicates for EnvVarsUsed for a given DataType
# Select out the top wrapper test DataType
DataTypeSel <- paste0("Top", FullNObs, "WrapperCV2Test")
MaxentKeepEvalStatsTopPart.df <- MaxentKeepEvalStats.df2[MaxentKeepEvalStats.df2$DataType==DataTypeSel,]
nrow(MaxentKeepEvalStatsTopPart.df)
# Remove duplicates
MaxentKeepEvalStatsTopPart.df2 <- MaxentKeepEvalStatsTopPart.df[!duplicated(MaxentKeepEvalStatsTopPart.df$EnvVarsUsed), ]
head(MaxentKeepEvalStatsTopPart.df2)
nrow(MaxentKeepEvalStatsTopPart.df2)
# Save only VarNames
TopVarNamesKeep <- MaxentKeepEvalStatsTopPart.df2$VarNames
length(TopVarNamesKeep)
## Keep only above VarNames in all Top sets
# First isolate all Random Sets
MaxentKeepEvalStatsTop.df <- MaxentKeepEvalStats.df2[grepl("Top", MaxentKeepEvalStats.df2$DataType), ]
head(MaxentKeepEvalStatsTop.df)
tail(MaxentKeepEvalStatsTop.df)
nrow(MaxentKeepEvalStatsTop.df)
# Only keep original data with VarNames matching that in above data
MaxentKeepEvalStatsTop.df1 <- MaxentKeepEvalStatsTop.df[MaxentKeepEvalStatsTop.df$VarNames %in% TopVarNamesKeep, ]
nrow(MaxentKeepEvalStatsTop.df1)
########
### Join back together Random and Top data cleaned of duplicates
MaxentKeepEvalStats.df2 <- rbind(MaxentKeepEvalStatsRand.df1, MaxentKeepEvalStatsTop.df1)
nrow(MaxentKeepEvalStats.df2)
# Sort by DataType and Rep
MaxentKeepEvalStats.df2 <- MaxentKeepEvalStats.df2[with(MaxentKeepEvalStats.df2, order(MaxentKeepEvalStats.df2$DataType, MaxentKeepEvalStats.df2$Rep)), ]
nrow(MaxentKeepEvalStats.df2)
########################
ModelTypeNames <- c(paste0("Top",FullNObs), paste0("Random", FullNObs))
MaxentKeepEvalStatsL <- list()
Count <- 0
RepList <- seq(1:3)
##################
for(Rep in RepList) {
  #Rep=1
  for(ModelType in ModelTypeNames) {
    #ModelType <- paste0("Top",FullNObs)
    # Keep only data for Rep and ModelType
    MaxentTestKeepEvalStats.df3 <- MaxentKeepEvalStats.df2[which(MaxentKeepEvalStats.df2$Rep==Rep & MaxentKeepEvalStats.df2$DataType==paste0(ModelType, "WrapperCV2Test")),]
    # Keep only first NumProjections/3 number of rows
    Count <- Count + 1
    if(NumProjections > 6) {
      MaxentKeepEvalStatsL[[Count]] <- MaxentTestKeepEvalStats.df3[1:(NumProjections/3),]
    } else {
      MaxentKeepEvalStatsL[[Count]] <- MaxentTestKeepEvalStats.df3[1:2, ]
    }
  }
}
#############
MaxentKeepEvalStats.df4 <- do.call(rbind, MaxentKeepEvalStatsL)
head(MaxentKeepEvalStats.df4)
## Sort data by DataType
MaxentKeepEvalStats.df4s <- arrange(MaxentKeepEvalStats.df4, DataType)
nrow(MaxentKeepEvalStats.df4s)
## Save data
setwd(OutDirectpsa)
write.table(MaxentKeepEvalStats.df4s, file=paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets, "TopSetsWrapperCV2TestStats_of", FinalModelVariables, "of", TotVars, "NumVars_.csv"), sep=",", col.names=NA)
######################################
## Loop through sets 1 through FinalModelVariables
# Read back in data if needed
setwd(OutDirectpsa)
MaxentKeepEvalStats.df4s <- data.frame(read.csv(paste0(Species, ModelName, Rank, FullNObs, "_", NumberModSets, "TopSetsWrapperCV2TestStats_of", FinalModelVariables, "of", TotVars, "NumVars_.csv"), sep=",", row.names=1, stringsAsFactors=FALSE))
head(MaxentKeepEvalStats.df4s)
#
# Read Albers ESRI shapefile to obtain proper CRS
setwd(InDirect)
#AdminAlbersShp <- readOGR(".", paste0("AdminGlobalAlbers"))
#CRS.NAAlbersEqualAreaConic <- CRS(proj4string(AdminAlbersShp))
######
#ModelTypeNames <- c(paste0("Top",FullNObs), paste0("Random", FullNObs))
ModelTypeNames <- c(paste0("Top",FullNObs))
#ModelTypeNames <- c(paste0("Random", FullNObs))
MaxentKeepEvalStatsL <- list()
Count <- 0
###############################################################################
#### Create and save NObsJ usually three kfold partition scheme for multiple rounds of test with presence
#### pseudoabsence and background points using all data
#### NOTE: Only create and save partition scheme once in order to reuse
#setwd(InDirect)
#library(dismo)
##
#NObskfoldgrpp <- kfold(PresenceDat.df, NObsJ)
#NObskfoldgrpa <- kfold(PseudoabsenceDat.df, NObsJ)
### Save kfold partition scheme for later testing
#setwd(InDirect)
#write.table(data.frame(NObskfoldgrpp), file=paste0(Species, "PresenceDat", NObsJ, "Kfold", ".csv"), sep=",")
#write.table(data.frame(NObskfoldgrpa), file=paste0(Species, "_", PsAbsBuffs, "PseudoabsenceDat", NObsJ, "Kfold", ".csv"), sep=",")
########################################################################################
## Use 3-Fold run and use only one fold for presence data to derive threshold (save time): 2/3 train, 1/3 test
## But use ALL of the pseudoabsence data for testing (not the fold)
######################
#### FIRST: Find the threshold values for the Maxent model using a training grid
#### model
## Read in three-fold partition scheme
setwd(InDirect)
NObskfoldgrpp <- as.vector(unlist(read.csv(paste0(Species, "PresenceDat", NObsJ, "Kfold", ".csv"))))
length(NObskfoldgrpp)
TotPres <- length(NObskfoldgrpp)
#
NObskfoldgrpa <- as.vector(unlist(read.csv(paste0(Species, "_", PsAbsBuffs, "PseudoabsenceDat", NObsJ, "Kfold", ".csv"))))
if(nrow(PseudoabsenceDat.df)>MaxBackgrndPseudoabsLimit) {
  NObskfoldgrpa <- NObskfoldgrpa[1:MaxBackgrndPseudoabsLimit ]
}
length(NObskfoldgrpa)
# Define coordinate system
setwd(InDirect)
CRSIn <- crs.out
#################################################################################
#Count <- 21
t1 <- Sys.time()
#############################################################################
for(ModelType in ModelTypeNames) {
  #ModelType <- paste0("Top",FullNObs)
  # Keep only DataTypes including string of ModelType for Final WrapperCV2 Test
  MaxentKeepEvalStats.df5 <- MaxentKeepEvalStats.df4s[grep(paste0(ModelType, "WrapperCV2Test"), MaxentKeepEvalStats.df4s$DataType), ]
  head(MaxentKeepEvalStats.df5)
  MaxentKeepEvalStats.df5[1:20,]
  #
  ProjList <- c(1,2)
  #ProjList <- c(1,2,3,4,5,6,7,8,9,10,11,12)
  #ProjList <- c(11)
  #for(k in 1:NumProjections) {
  #for(k in 9:NumProjections) {
  for(k in ProjList) {
    #k=12
    Loop <- k
    ######################################################################################
    ## Run Threshold Calibration Projection for Top Selected Variable Model for Given Variable Set (j)
    ######################################################################################
    ##
    #################
    ## NOTE: DO NOT RESORT BY THIS AUCfinaltest- THIS WILL INFLATE AUC- AUCWrappertest already used to sort
    #################
    VariableSubset <- data.frame(as.character(MaxentKeepEvalStats.df5[k,1]), stringsAsFactors=FALSE)
    colnames(VariableSubset) <- "VarNames"
    VarNames <- unlist(strsplit(as.character(VariableSubset[,1]),"-"))
    SubsetVariableNumber <- length(VarNames)
    Rep <- MaxentKeepEvalStats.df5[k,8]
    #
    #VarNames <- VarNames[-14]
    #SubsetVariableNumber <- length(VarNames)
    #VariableSubset <- as.data.frame(paste0(unlist(VarNames), collapse="-"), stringsAsFactors=FALSE)
    #colnames(VariableSubset) <- "VarNames"
    #### FIRST: Find the threshold values for the Maxent model using a training grid
    #### model
    #
    Run=1
    ##
    VariableNamesIn <- VariableNames
    #
    ###########################################
    ##
    ##
    setwd(OutDirectpsa)
    output3 <- paste0(ModelType, "_", SubsetVariableNumber, "Vars_", Loop, SetNameF)
    #output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", Loop, "Minus", 14)
    dir.create(paste0(OutDirectpsa, "/", output3))
    OutDirectSub <- paste0(OutDirectpsa, "/", output3)
    OutDirectIn <- OutDirectSub
    ##
    kfoldgrppin <- NObskfoldgrpp
    kfoldgrpain <- NObskfoldgrpa
    #
  #  # If necessary, read back in r grids
    #setwd(MaskGrids)
  #  Predictors <- stack(GridNamesL)
  #  names(Predictors) <- toupper(names(Predictors))
  #  names(Predictors) <- gsub("_NS", "", names(Predictors))
    PredictorsIn <- subset(Predictors, VarNames)
    #PredictorsIn[[75:90]]
    #plot(PredictorsIn[[84]])
    #writeRaster(PredictorsIn[[86]], paste0("Test86"), format = "GTiff", overwrite=TRUE)
    #
    OutGridID <- paste0(ModelType, "_", Loop, "ThreshCal_", SetNameF)
    # Takes 19 minutes for 25 variable grid
    #
    SetNameIn <- k
    CRS.In <- CRSIn
    #
    Time <- system.time(MaxentKeepEvalStats.df <- MaxentSubset.GridTrainTestEvalAIC_AUCbgp_Calib(Species, VariableNamesIn, PredictorsIn, SubsetVariableNumber, TotPres, Run, SetName, VariableSubset, PresenceDat.df, PseudoabsenceDat.df, BackgroundDat.df, kfoldgrppin, kfoldgrpain, CRS.In, OutDirectIn, MaxentArgsIn=MaxentBaseArgsIn,  Output=TRUE))
    #
    head(MaxentKeepEvalStats.df)
    # Save results
    setwd(OutDirectIn)
    write.table(MaxentKeepEvalStats.df , file=paste0(Species, "Maxent", Rank, DataSet, "_TrainVsTestStatsThreshCal_of", NumberModSets, "_", SetName, SubsetVariableNumber, "_", Loop, ".csv"), sep=",", col.names=NA)
    #
    #MaxentMod <- readRDS("MaxentModel.rds") # Reads in stored Maxent model from above function
    Count <- Count + 1
    ## Add ModelType to output
    MaxentKeepEvalStats.df$ModelType <- ModelType
    ncol(MaxentKeepEvalStats.df)
    # Re-order columns
    MaxentKeepEvalStats.df <- MaxentKeepEvalStats.df[,c(1:3, 17, 3:16)]
    # Replace Run with Rep
    MaxentKeepEvalStats.df$Rep <- Rep
    ncol(MaxentKeepEvalStats.df)
    MaxentKeepEvalStats.df <- MaxentKeepEvalStats.df[,c(1:7, 18, 8:17)]
    MaxentKeepEvalStatsL[[Count]] <-  MaxentKeepEvalStats.df
  }
}
t2 <- Sys.time()
#############################################################################
difftime(t2,t1, units = "mins")
###
MaxentKeepEvalStatsAll.df <- do.call(rbind, MaxentKeepEvalStatsL)
## Sort by DataType
MaxentKeepEvalStatsAll.dfs  <- arrange(MaxentKeepEvalStatsAll.df, DataType)
head(MaxentKeepEvalStatsAll.dfs)
## Save data
setwd(OutDirectpsa)
write.table(MaxentKeepEvalStatsAll.dfs, file=paste0(Species, "Maxent", Rank, DataSet, "_TrainVsTestStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), sep=",", col.names=NA)

##################################################################################
### Analyze above data  for linear regression between AICc and AICc_bg
## Read above data back in
# Read data back in and sort
#SubsetVariableNumber <- 8
#NumberModSets <- 10000
#setwd(OutDirectpsa)
#MaxentKeepEvalGridStats.df <- data.frame(read.csv(paste0(Species, "Maxent", Rank, DataSet, "_TrainVsTestStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
###
### Keep only DataType of Test
#MaxentKeepEvalGridStats.df1 <- subset(MaxentKeepEvalGridStats.df, DataType=="Test")
## Test for significant linear regression between AICc and AICcbg
#AICc_AICcbg.lm <- lm(AICc ~ AICc_bg, data=MaxentKeepEvalGridStats.df1)
#AICc_AICcbg.lmRes <- summary(AICc_AICcbg.lm)
## Save above output
#out1 <- paste0("Linear regression F test for AICc Versus AICc_bg from", NumberModSets, SetNameF, "Setsof", SubsetVariableNumber, "Variables")
#cat(out1, file=paste0(Species, Rank, "AICcVsAICc_bg_LinearReg_by", NumProjections, "from", NumberModSets, SetNameF, "Setsof", SubsetVariableNumber, "VariablesAICcVsAICc_bg.txt"), sep="\n", append=TRUE)
#out3 <- capture.output(print(AICc_AICcbg.lmRes))
#cat(out3, file=paste0(Species, Rank, "AICcVsAICc_bg_LinearReg_by", NumProjections, "from", NumberModSets, SetNameF, "Setsof", SubsetVariableNumber, "VariablesAICcVsAICc_bg.txt"), sep="\n", append=TRUE)

#######################################################################################
### Summarize performance statistics of projected models
#######################################################################################
##
## Retrieve top selected variable from selected subset of 90 variables
# Specify variable number and row number of selected set
FinalNumberModels <- 12 # Select number of top models to keep by Regional Indices
FinalModelVariables <- 6  # Chosen variable set size based trends in AUC and AICc
SubsetVariableNumber <- FinalModelVariables
NumProjections <- 12
NumberModSets <- 3000
#
#
j = 1
#j=1
if(j < 10) {
  SetNumber <- paste0("0", j)
} else {
  SetNumber <- as.character(j)
}
#
SetName <- paste0(SetNumber, "_", colnames(VariableSetTypes[j]))
SetNameF <- ""
SetName <- paste0("01_", SetNameF)

## If above grid making loop interrupted, read back in evaluation data from grid directories and reassemble
#####################################################################################################
MaxentKeepEvalStatsL <- list()
ModelTypeNames <- c(paste0("Top",FullNObs))
count1 <- 0
####################################################
for(ModelType in ModelTypeNames) {
  #ModelType <- ModelTypeNames[1]
  count <- 0
  for(j in 1:FinalNumberModels) {
    #j=1
    setwd(OutDirectpsa)
    output3 <- paste0(ModelType, "_", SubsetVariableNumber, "Vars_", j, SetNameF)
    #output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", Loop, "Minus", 14)
    OutDirectSub <- paste0(OutDirectpsa, "/", output3)
    setwd(OutDirectSub)
    MaxentKeepEvalStats.df <- data.frame(read.csv(paste0(Species, "Maxent", Rank, DataSet, "_TrainVsTestStatsThreshCal_of", NumberModSets, "_", SetName, SubsetVariableNumber, "_", j, ".csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
    count1 <- count1 + 1
    MaxentKeepEvalStatsL[[count1]] <- MaxentKeepEvalStats.df
  }
}
###################################################
MaxentKeepEvalStatsAll.df <- do.call(rbind, MaxentKeepEvalStatsL)
## Sort by DataType
MaxentKeepEvalStatsAll.dfs  <- arrange(MaxentKeepEvalStatsAll.df, DataType)
head(MaxentKeepEvalStatsAll.dfs)
## Save data
setwd(OutDirectpsa)
write.table(MaxentKeepEvalStatsAll.dfs, file=paste0(Species, "Maxent", Rank, DataSet, "_TrainVsTestStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), sep=",", col.names=NA)
####################################
### Summarize Output
###Calculate Mean and Standard Deviation Values
RunType <- "WrapperCV2"
EvaluationStats <- MaxentKeepEvalStatsAll.dfs
ncol(EvaluationStats)
EvaluationStats$Model <- paste0(ModelName, Rank)
EvaluationStats <- EvaluationStats[, colnames(EvaluationStats)[c(19,1:18)]]
head(EvaluationStats)
nrow(EvaluationStats)
# Delete any subsets with "Inf" for AICc_bg (means too many parameters in model relative to observations)
EvaluationStats <- subset(EvaluationStats, AICc_bg!="Inf")
# Replace NAN values with zero
EvaluationStats[is.na(EvaluationStats)] <- 0
##
#FSAType <- Rank
OutName <- paste0(Species, ModelName, Rank, "_TrainVsTestStatsThreshCal_SummaryStats_", SubsetVariableNumber, "of", TotVars, "Vars", FSAType,"_", RunType, ".csv")
SortGroups <- c("Model", "DataType", "SubsetVariableNumber", "TotalPresPnts")
## All statistics, and only statistics, should be at and after the column specified below
StatVarFirstColumn <- 6
OutDirectIn <- OutDirectpsa
EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
head(EvaluationStatsSummary.df)
######

#####################################################################################################
### Run maxent 3-fold cross-valication evaluation mode with jacknife option to assess variable importance with all variables
### for top models
#####################################################################################################

## Takes about 17 minutes with 15 variables and three models
## Takes 8 minutes with 12 variables and three models
#############################################
# Read in data with rankings
SubsetVariableNumber <- 6
NumProjections <- 12
FinalNumberModels <- 12
Rank <- "AIC"
NumberModSets <- 3000
SetNameF <- ""
DataSet <- "All"
OutDirectpsa <- paste0(Direct2, output1, "/", EvalType, "/", Rank, "Ranked")
setwd(OutDirectpsa)
MaxentKeepEvalAreaStatsAll.dfs <- data.frame(read.csv(paste0(Species, "Maxent", Rank, DataSet, "_TrainVsTestStatsThreshCal_of", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_", EvalType, "_TopVsRandom", NumProjections, "Sets.csv"), header = TRUE, row.names=1, sep=','), stringsAsFactors=FALSE)
head(MaxentKeepEvalAreaStatsAll.dfs )
# Keep only rankings for model type
MaxentKeepEvalAreaStatsAll.df2 <- MaxentKeepEvalAreaStatsAll.dfs[which(MaxentKeepEvalAreaStatsAll.dfs$DataType=="Test"),]
# Keep only DataTypes including string for Top FinalTest
TopSelSubsets <- as.data.frame(MaxentKeepEvalAreaStatsAll.df2[1:FinalNumberModels,1], stringsAsFactors=FALSE)
colnames(TopSelSubsets) <- "VarNames"
##
FinalModelRowNumberList.df <- as.data.frame(c(1:NumProjections))
NumberTopModels <- 12 # Specify number of top models on which to run jacknife analysis
TopModelRowNumberList <- FinalModelRowNumberList.df[1:NumberTopModels,]
##############
t1 <- Sys.time()
for(NumberModel in 1:NumberTopModels) {
  #NumberModel <- 1
  TopModel <- FinalModelRowNumberList.df[NumberModel,]
  VariableSubset <- data.frame(as.character(TopSelSubsets[NumberModel,1]), stringsAsFactors=FALSE)
  colnames(VariableSubset) <- "VarNames"
  VarNames <- strsplit(as.character(VariableSubset[,1]),"-")
  SubsetVariableNumber <- length(unlist(VarNames))
  #########
  #
  VariableNamesIn <- VariableNames
  # Specify SubsetVariableNumber
  # Specify model training data
  MaxentPresTrainData <- PresenceDat.df
  head(MaxentPresTrainData)
  MaxentAbsTrainData <- BackgroundDat.df
  head(MaxentAbsTrainData)
  TrainSWD <- rbind(MaxentPresTrainData[,4:ncol(MaxentPresTrainData)], MaxentAbsTrainData[,4:ncol(MaxentAbsTrainData)])
  head(TrainSWD)
  TrainPresID <- data.frame(rep(1,nrow(MaxentPresTrainData)))
  colnames(TrainPresID) <- "ID"
  TrainAbsID <- data.frame(rep(0,nrow(MaxentAbsTrainData)))
  colnames(TrainAbsID) <- "ID"
  TrainPresAbsID <- rbind(TrainPresID, TrainAbsID)
  head(TrainPresAbsID)
  tail(TrainPresAbsID)
  #
  # Subset training and testing data by selected variables
  TrainSWD <- TrainSWD[,unlist(VarNames)]
  head(TrainSWD)
  ######
  # Create subdirectory for Maxent output
  output3 <- paste0("MaxentEvaluationMode_FinalSubset_", Rank, SetNameF, SubsetVariableNumber, "_Set", TopModel)
  dir.create(paste0(OutDirectpsa, "/", output3))
  OutDirectSub <- paste0(OutDirectpsa, "/", output3)
  setwd(OutDirectSub)
  #
  #system.file("java", package="dismo")
  ## Check if have categorical variable to add to maxent arguments
  if(exists("CatVarsPrefix")==TRUE) {
    MaxentCatArg <- paste0("togglelayertype=", CatVarsPrefix)
    MaxentEvalArgs <- c(MaxentEvalArgs1, MaxentCatArg)
    } else {
    MaxentEvalArgs <- MaxentEvalArgs1
  }
  # NOTE: the below jackknife run takes about 7.5 minutes with 20 variables
  # Takes 85 seconds for 8 variables
  ## Takes 42 minutes for 86 variables
  system.time(MaxentOut <- maxent(TrainSWD, TrainPresAbsID, args=MaxentEvalArgs, path=OutDirectSub))
}
t2 <- Sys.time()
###################################################################################
difftime(t2,t1, units = "mins")

