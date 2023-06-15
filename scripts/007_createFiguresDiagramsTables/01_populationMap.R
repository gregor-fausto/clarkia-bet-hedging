####
####
# Script to create map of study populations
# Produces Figure 1A
#
# Plots follow tutorial at
# http://zevross.com/blog/2016/03/15/using-the-new-r-package-feddata-to-access-federal-open-datasets-including-interactive-graphics/
# projections obtained at
# https://github.com/mousebird/proj-4/blob/master/nad/esri.extra
# https://www.neonscience.org/resources/learning-hub/tutorials/raster-data-r
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

# - Libraries ----
library(tidyverse)
library(proj4)
library(FedData)
library(rgeos)
library(RColorBrewer)
library(tmap)
library(rgdal)
library(khroma) # iridescent color scheme

source("scripts/007_createFiguresDiagramsTables/00_utilityFunctions.R")

# - Define temporary directory ----
# define temporary directory to hold shapefiles and rasters for plotting
currentDirectory <- getwd()
tmpDirectory <- "outputs/007_createFiguresDiagrams/"
figureDirectory <- "products/figures/"

# - Import data ----
# select population/site name, easting/northing coordinates (units of meters) 
data<-read.csv(file="data/siteAbioticData.csv",header=TRUE) %>% 
  dplyr::select(site,easting,northing)

# - Transform spatial data ----
# (easting, northing) -> (lat, lon)

# coordinate reference system to convert between coordinate systems
proj4string <- "+proj=utm +zone=11 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs "

# select easting and northing
locations <- data %>% dplyr::select(easting, northing) 

# transform data: inverse projection from cartographic to lat/long
pj <- proj4::project(locations, proj4string, inverse=TRUE)
latlon <- data.frame(lat=pj$y, long=pj$x)

# bind to data frame
data <- bind_cols(data,latlon)

# - Define spatial coverage/extent for raster download ----

# extent was refined manually
# define shape polygon based on lat/long coordinates
extentClarkia <- rgeos::readWKT("POLYGON((-118.25 35.45, -118.25 35.6, -118.25 35.9, -118.8 35.9, -118.8 35.45, -118.25 35.45))")
proj4string(extentClarkia) <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"

# - Get national elevation raster ----
# download National Elevation Database elevation data in extent
# default resolution is 1 arcsecond (res="1);
# to get 1/3 arcsecond (res="13)

currentDirectory<-getwd()

setwd(tmpDirectory)

# download National Elevation Database elevation data in extent
ned_kern<-FedData::get_ned(template=extentClarkia, label="ned_kern", res="1", force.redo = F)

setwd(currentDirectory)


# - Download shapefiles for lakes and rivers ----

# this requires downloading shapefiles for the rivers and water bodies
# I did this by downloading shapefiles via the system (Terminal on Mac)
# sources for shapefiles are included here

# # California hydrography: https://www.calfish.org/ProgramsData/ReferenceLayersHydrography/CaliforniaHydrography.aspx
# system(paste0("cd ",getwd(),"/hydrography/;
#        curl -O ftp://ftp.streamnet.org/pub/calfish/cdfg_100k_2003_6.zip ;
#        unzip cdfg_100k_2003_6.zip;
#        cd -"))

# water bodies of Kern County: https://geodata.lib.berkeley.edu/catalog/berkeley-s7gm4h
# system(paste0("cd ",getwd(),"/waterbodies/;
#         curl -O https://spatial.lib.berkeley.edu/public/berkeley-s7gm4h/data.zip ;
#         unzip data.zip;
#         cd -"))

# - Shapefile for Lake Isabella ----

setwd("outputs/007_createFiguresDiagrams/waterbodies")
setwd(paste0(currentDirectory,"/outputs/007_createFiguresDiagrams/waterbodies"))

bodies <- rgdal::readOGR(dsn="Kern_Waterbodies_2014.shp",layer="Kern_Waterbodies_2014")
isabella <- subset(bodies, NAME %in% c('Lake Isabella'))
isabella<-spTransform(isabella, proj4string(ned_kern))

# - Shapefile for Kern River ----

setwd(paste0(currentDirectory,"/outputs/007_createFiguresDiagrams/hydrography"))

rivers2 <- rgdal::readOGR(dsn="cdfg_100k_2003_6.shp",layer="cdfg_100k_2003_6")
kern2 <- subset(rivers2, NAME %in% c('Kern River'))
kern2<-spTransform(kern2, proj4string(ned_kern))
kern_to_plot<-raster::crop(kern2,isabella)
kernFinal=gDifference(kern2, kern_to_plot)

setwd(currentDirectory)

# - Points for population locations ----

xy <- data.frame(data$long,data$lat)
spdf <- SpatialPoints(coords = xy, 
                      proj4string = CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"))

# Reproject and change scale

# for reprojetion, follow the guidelines here
# https://datacarpentry.org/r-raster-vector-geospatial/03-raster-reproject-in-r/
# https://github.com/SpaCE-Lab-MSU/MSUGradSpatialEcology/blob/master/lab4_autocorrelation.Rmd
# https://stackoverflow.com/questions/15634882/why-the-values-of-my-raster-map-change-when-i-project-it-to-a-new-crs-projectra

# convert scale to km

cord.UTM <- sp::spTransform(spdf, CRS("+init=epsg:32629 +zone=11 +north"))

# Reprojection - followed the guidance here
# https://datacarpentry.org/r-raster-vector-geospatial/03-raster-reproject-in-r/
# https://github.com/SpaCE-Lab-MSU/MSUGradSpatialEcology/blob/master/lab4_autocorrelation.Rmd
# https://stackoverflow.com/questions/15634882/why-the-values-of-my-raster-map-change-when-i-project-it-to-a-new-crs-projectra
cord.UTM <- sp::spTransform(spdf, CRS("+init=epsg:32629 +zone=11 +north"))
cord.UTM@coords <- cord.UTM@coords/1000
cord.UTM <- data.frame(cord.UTM)

isabella.UTM <- sp::spTransform(isabella, CRS("+init=epsg:32629 +zone=11 +north"))
isabella.UTM@polygons[[1]]@Polygons[[1]]@coords <- isabella.UTM@polygons[[1]]@Polygons[[1]]@coords/1000

# split kern into lower and upper halves
kernLower = kernFinal
kernUpper = kernFinal

kernLower@lines[[1]]@Lines <- kernFinal@lines[[1]]@Lines[-2]
kernUpper@lines[[1]]@Lines <- kernFinal@lines[[1]]@Lines[-1]

kernLower <- sp::spTransform(kernLower, CRS("+init=epsg:32629 +zone=11 +north"))
kernLower@lines[[1]]@Lines[[1]]@coords <- kernLower@lines[[1]]@Lines[[1]]@coords/1000
kernLower.UTM <- SpatialLinesDataFrame(kernLower,as.data.frame("kern.lo"))

kernUpper <- sp::spTransform(kernUpper, CRS("+init=epsg:32629 +zone=11 +north"))
kernUpper@lines[[1]]@Lines[[1]]@coords <- kernUpper@lines[[1]]@Lines[[1]]@coords/1000
kernUpper.UTM <- SpatialLinesDataFrame(kernUpper,as.data.frame("kern.hi"))


shape.UTM <- raster::projectRaster(ned_kern, crs=CRS("+init=epsg:32629 +zone=11 +north +no_defs"))

# - Define elevation colors ----
# colorblind friendly

iridescent <- khroma::colour("iridescent")
plot_scheme(iridescent(23), colours = TRUE, size = 0.5)
myColIridescent=rev(iridescent(23))

# - Plot map ----
# https://github.com/mtennekes/tmap/issues/164
# https://rstudio-pubs-static.s3.amazonaws.com/289387_6d2bbaf850764a38bfea1618e78a68ef.html

# set font sizes
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12

# - +Plot map in UTM for submission ----

#   Convert rasters TO dataframes for plotting with ggplot
hdf <- raster::rasterToPoints(shape.UTM); hdf <- data.frame(hdf)
colnames(hdf) <- c("x","y","elev")
hdf <- hdf %>%
  dplyr::mutate(easting = x/1000, northing = y/1000)

dataLabel = data

# manually adjust the position of population labels
dataLabel[1,2:3] = dataLabel[1,2:3]+c(-1750,1000)
dataLabel[2,2:3] = dataLabel[2,2:3]+c(2000,-500)
dataLabel[3,2:3] = dataLabel[3,2:3]+c(-1000,1250)
dataLabel[4,2:3] = dataLabel[4,2:3]+c(2250,0)
dataLabel[5,2:3] = dataLabel[5,2:3]+c(0,-1750)
dataLabel[6,2:3] = dataLabel[6,2:3]+c(-1000,1600)
dataLabel[7,2:3] = dataLabel[7,2:3]+c(0,-1750)
dataLabel[8,2:3] = dataLabel[8,2:3]+c(-2000,500)
dataLabel[9,2:3] = dataLabel[9,2:3]+c(2350,0)
dataLabel[10,2:3] = dataLabel[10,2:3]+c(0,1750)

dataLabel[11,2:3] = dataLabel[11,2:3]+c(2250,1250)
dataLabel[12,2:3] = dataLabel[12,2:3]+c(2500,-500)
dataLabel[13,2:3] = dataLabel[13,2:3]+c(-1750,1000)
dataLabel[14,2:3] = dataLabel[14,2:3]+c(2000,-500)
dataLabel[15,2:3] = dataLabel[15,2:3]+c(2750,-250)
dataLabel[16,2:3] = dataLabel[16,2:3]+c(1250,-1750)
dataLabel[17,2:3] = dataLabel[17,2:3]+c(-2250,0)
dataLabel[18,2:3] = dataLabel[18,2:3]+c(2250,0)
dataLabel[19,2:3] = dataLabel[19,2:3]+c(2250,0)
dataLabel[20,2:3] = dataLabel[20,2:3]+c(2500,-250)

# change color of labels
colLabels=rep("white",20)

# - +Reset working directory ----

setwd(currentDirectory)


# - +Plot map  ----

tiff(filename=paste0(figureDirectory,"map-ggplot-iridescent.tif"),
     units='px',height = 5600/2, width = 5400/2,res=800,compression='lzw')

ggplot() +
  # plot elevation raster
  geom_raster(data = hdf , aes(x = easting, y = northing, fill = elev),interpolate = TRUE) + 
  # set axis limits
  scale_x_continuous(limits = c(340,385), expand = c(0, 0),
                     name = "Easting (km)", breaks = c(340, 350, 360, 370, 380)) +
  scale_y_continuous(limits = c(3925,3972.5), expand = c(0, 0),
                     name = "Northing (km)", breaks = c(3930,3940,3950,3960,3970)) +
  # use iridescent color scheme for color-blind friendliness
  scale_fill_gradientn(colours = colorRampPalette(myColIridescent)(23), na.value = NA,name="Elevation (m)",
                       breaks=c(500,1500,2500)) +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 2.5)) +
  # plot Lake Isabella
  geom_polygon(data=isabella.UTM,aes(x = long, y = lat),alpha=.5,col='lightgray',fill='lightgray') +
  # plot the lower half of the Kern River
  geom_path(data=kernLower.UTM,aes(x = long, y = lat),col='lightgray',alpha=.5) +
  # plot the upper half of the Kern River
  geom_path(data=kernUpper.UTM,aes(x = long, y = lat),col='lightgray',alpha=.5) +
  # add the population locations
  geom_point(data=data,aes(x= easting/1000, y = northing/1000),fill="white", col = 'black', 
             alpha=.75,size=3,shape=21) +
  # add the population labels
  geom_text(data=dataLabel,aes(x= easting/1000, y = northing/1000,label=site),
            col = colLabels,size=1.9) +
  # add the compass rose
  ggsn::north(x.min=340,x.max=384,y.min=3925,y.max=3971.5,symbol=19,scale=.1) +
  # add the scalebar with scalebar2, which modified the original function to be in km
  scalebar2(x.min=341,x.max=384,y.min=3925,y.max=3971,location="topleft",
            transform=FALSE,dist=4,dist_unit="m", height = 0.015,
            st.size=3,border.size=.5,st.color="white") +
  #theme with white background
  theme_bw() +
  #eliminate background, gridlines, and chart border
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.text=element_text(size=10),
    axis.title=element_text(size=12),
    legend.text=element_text(size=6),
    legend.title=element_text(size=7),
    legend.justification=c(1,1),
    legend.position=c(.30,.925),
    legend.background = element_rect(fill=alpha('white', 0.8)))

dev.off()
