
### k-means on covariates representing soil forming factors
### Author: Mercedes Roman Dobarco
### Date: 25/03/2020
### Objective: Definition of groups with similar pedogenesis (~ homogeneous soil forming factors, excluding land use and current vegetation)
### this script: Prepare covariates - Scale the numerical variables
###                                 - Sample the categorical variables
###                                 - Perform PCA on dummy variables created from categorical variables
###                                 - Predict the dim of the PCA on for the whole raster
###                                 - Scale the dims of the PCA 
###                                 - Regular sample of all covariates

library(rgdal)
library(gdalUtils)
library(raster)
library(sp)
library(sf)
library(nngeo)
library(dplyr)
library(tidyverse)
library(MASS)
library(foreach)
library(parallel)
library(doParallel)
library(snow)
library(doSNOW)
library(ggplot2) # tidyverse data visualization package
library(viridis) # color palettes
library(scales)
library(rasterVis)
library(lattice)
library(gridExtra)
library(tmap)    # for static and interactive maps
library(leaflet) # for interactive maps
library(mapview) # for interactive maps
library(shiny)   # for web applications
library(ggmap)
library(geosphere)
library(slga)


# ### A. Covariate processing and pedogenon modeling --------------------------------

# 1. Select CLORPT variables ----------------------------------------------

### what are the covariates that represent the soil forming factors?
### To make the example reproduciblem we will work with rasters from the Soil and Landscape Grid of Australiam
### available with the package slga

### Here, I define the area of interest as a rectangle, but it could be any polygon shapefile
# aoi <- st_read("aoi.shp")
AOI <- c(150,-32, 151,-31)
### Download covariates representing the soil forming factors at the moment of chosen as reference
### climate
### 'Prescott Index', 'Net Radiation [January]', 'Net Radiation [July]'
prescott <-get_lscape_data(product = 'PSIND', aoi = AOI, write_out = FALSE)
nrjan <- get_lscape_data(product = 'NRJAN' , aoi = AOI, write_out = FALSE)
nrjul <- get_lscape_data(product = 'NRJUL' , aoi = AOI, write_out = FALSE)
### Relief
### 'Slope [percent]', 'Relief [300m radius]', 'Topographic Wetness Index', 'MrVBF', 'SRTM_TopographicPositionIndex'
slope <- get_lscape_data(product = 'SLPPC' , aoi = AOI, write_out = FALSE)
rel300m <- get_lscape_data(product = 'REL3C', aoi = AOI, write_out = FALSE)
twi <- get_lscape_data(product = 'TWIND', aoi = AOI, write_out = FALSE)
mrvbf <- get_lscape_data(product = 'MRVBF', aoi = AOI, write_out = FALSE)
tpi <- get_lscape_data(product = 'TPIND', aoi = AOI, write_out = FALSE)

### To make the example easy, we are only including covariates available from the slga package, 
### but they can include any other covariate for parent material, vegetation (including land use), or time.
### Maybe also soil stable properties

### Make stack
covariates.stack <- stack(prescott,nrjan,nrjul,slope,rel300m,twi,mrvbf,tpi)
names(covariates.stack)

#  2 - Scale the numerical variables ----------------------------------------

### what are the continuous covariates?
covariates.selection.cont <- c("SLGA_PSIND", "SLGA_NRJAN", "SLGA_NRJUL", "SLGA_SLPPC", "SLGA_REL3C", "SLGA_TWIND", "SLGA_MRVBF")
### This could also be a list of tif files saved in our computer. 
### Especially if we need to mask the values outside our area of interest
### Designate a directory for storing the output raters
# HomeDir <- "C:/Covariates/" 
# dir.create(HomeDir)
# OutDir <- "C:/Covariates/Output" 
# dir.create(OutDir)

detectCores()
cl <- makeCluster(2)   ### Create cluster
registerDoParallel(cl)
getDoParWorkers()

system.time(covariates.cont.out <- foreach(i=1:length(covariates.selection.cont), .packages=c("raster", "sf"),
                                           .export = c("covariates.selection.cont")) %dopar% {
  # setwd(paste0(HomeDir)) # Set wd
  # r <- raster(covariates.selection.cont[[i]]) ### load raster
  # m <- mask(r, aoi) # Mask pixels outside the area of interest
  # s <- scale(x=m,center=TRUE, scale=TRUE) # Scale because it is a continuous variable
  s <- scale(x=covariates.stack[[covariates.selection.cont[[i]]]],center=TRUE, scale=TRUE) # Scale because it is a continuous variable
  # setwd(OutDir) # change wd
  # writeRaster(s, filename = paste0(covariates.selection.cont[[i]],"_sc.tif" ),
  #             format = "GTiff",na.rm=T, inf.rm=T, overwrite = T) # Write to file
  s # Return masked and scaled raster
})

stopCluster(cl)
covariates.cont.sc <- stack(covariates.cont.out)
plot(covariates.cont.sc)
rm(covariates.cont.out,i,s,m,r)


# 3 - Categorical variables --------------------------------------

### Mask categorical variables if your study area is bigger
#cat <- c("SLGA_TPIND")
#for (i in 1:length(cat)){
  # setwd(paste0(HomeDir)) # Set wd
  # r <- raster(cat[[i]]) ### load raster *if from file
  # setwd(OutDir) # change wd
  # m <- mask(r,aoi, # Mask pixels outside NSW
  #          filename = paste0(cat[[i]],"_mask.tif"),
  #          format = "GTiff",na.rm=T, inf.rm=T, overwrite = T) 
#}
#rm(i, cat)

### Create dummy variables
tpi.classes <- sort(unique(getValues(tpi)))
dummy.list <- list()
for(i in 1:length(tpi.classes)){
  print(i)
  # dummy.tpi <- calc(tpi, fun = function(x) {ifelse(is.na(x), NA, ifelse(x==tpi.classes[[i]],1,0))},
  #                    filename=paste0("tpi.dummy",tpi.classes[[i]],".tif"),
  #                    overwrite=TRUE)
  dummy.tpi <- calc(tpi, fun = function(x) {ifelse(is.na(x), NA, ifelse(x==tpi.classes[[i]],1,0))})
  dummy.list[[i]] <- dummy.tpi
}

dummy.stack <- stack(dummy.list)
plot(dummy.stack)

### Perform PCA
set.seed(1946)
# Regular sampling
sampleTPI <- sampleRegular(dummy.stack, size = 600000 , xy=TRUE, na.rm = TRUE, sp = FALSE)
### select only complete cases
sampleTPI <-sampleTPI[complete.cases(sampleTPI),]
sampleTPI <- as.data.frame(sampleTPI)
### Here I take a sample, but this could be done with the whole raster converted as a dataframe

### Apply PCA
tpi.pca <- prcomp(sampleTPI[,3:ncol(sampleTPI)], scale=TRUE)
# Eigenvalues
library(factoextra)
eig.val <- get_eigenvalue(tpi.pca)
eig.val ### 2 components have 100 % of the variability

### If  the study area is very large this step may be done in parallel
tpi.pcs <- predict(dummy.stack, tpi.pca, index=1:2)
plot(tpi.pcs)
### Write them to file 
# for (i in 1:nlayers(tpi.pcs)){
#   writeRaster(tpi.pcs[[i]], filename = paste0("TPI_PC",i,".tif"), overwrite=TRUE, format = "GTiff" )
# }
tpi.stack <- stack(tpi.pcs)
### You can scale them afterwards... I don't know if this step is necessary
 for (i in 1:nlayers(tpi.stack)){
   tpi.stack[[i]] <- scale(x=tpi.stack[[i]],center=TRUE, scale=TRUE)
 }
rm(tpi.pcs,tpi.pca, veg.rast,  dummy.stack, eig.val ,sampleTPI,dummy.list,tpi.classes,i,cl,dummy.tpi, covariates.selection.cont)


# 4. Regular sample of all covariates --------------------------------------

### Make stack
covariates.stack <- stack(covariates.cont.sc,tpi.stack);plot(covariates.stack)
names(covariates.stack) <- c("SLGA_PSIND", "SLGA_NRJAN", "SLGA_NRJUL", "SLGA_SLPPC",
                             "SLGA_REL3C", "SLGA_TWIND" ,"SLGA_MRVBF", "TPI_PC1","TPI_PC2" )
# Set the random number generator to reproduce the results
set.seed(1946)
# Regular sampling
sampleCLORPT<- sampleRegular(covariates.stack, size = 50000 , xy=TRUE, na.rm = TRUE, sp = TRUE)

# transform to sf
sampleCLORPT <- st_as_sf(sampleCLORPT)
### Transform into a dataframe
CLORPT.df <- st_drop_geometry(sampleCLORPT)
### select only complete cases
CLORPT.df <-CLORPT.df[complete.cases(CLORPT.df),]

## Clean 
rm(sampleCLORPT)


# 5. Apply Cholesky decomposition to sample data --------------------------

# The basic Euclidean distance treats each variable as equally important in calculating the distance.
# An alternative approach is to scale the contribution of individual variables to the distance value according
# to the variability of each variable. This approach is illustrated by the Mahalanobis distance, 
# which is a measure of the distance between each observation in a multidimensional cloud of points and
# the centroid of the cloud.

### Calculate the Mahalanobis distance, as Euclidean distance after applying the Cholesky decomposition
## Take out the coordinates
CLORPT.df.coords <- CLORPT.df[,1:2]
CLORPT.df.subset <- CLORPT.df[,3:ncol(CLORPT.df)]

## What is the average distance between pixels?
d1 <-distm(CLORPT.df.coords[200:250,], fun = distVincentyEllipsoid)
d1[1:10,1:10] ### the distance is 158 m and 238 m

# Rescale the data - Only first 22 (numerical variables)
C <- chol(var(as.matrix(CLORPT.df.subset)))
CLORPT.rs <- as.matrix(CLORPT.df.subset) %*% solve(C)

### euclidean distance of the first 10
b <- dist(CLORPT.rs[1:10,],method = "euclidean")
library(biotools)
a <- D2.dist(CLORPT.df.subset[1:10,], cov=var(CLORPT.df.subset));sqrt(a);b ### distances are the same

### Clean
rm(a,b,d1)


#############################################################################################################################################################

### 5. Optimal number of clusters

### how can we calculate the optimal number of clusters?
### there are several methods.
### A first one is to check the sum of within-cluster distances across all clusters
### Also, the package NbClust can calculate several indices to choose the optimal k

### With the package ClusterR
library(ClusterR)
set.seed(1991)
search_space <- c(2:20) ### Let's say that we want to check between 2 to 10 clusters
system.time(opt_kmeans <- Optimal_Clusters_KMeans(data = CLORPT.rs, max_clusters = search_space,
                                                  criterion = "WCSSE", num_init = 10,
                                                  max_iters = 1000, initializer = "kmeans++", plot_clusters = TRUE,
                                                  verbose = TRUE))
plot(search_space, opt_kmeans,pch=20, col="blue")
lines(search_space, opt_kmeans)


### An index that combines the within cluster similarity and the between cluster dissimilarity is the Calinski-Harbasz index

### Calculate the optimal number of k-means clusters with NbClust package
library(NbClust)
set.seed(1984)
system.time(opt.Clusters.CH <- NbClust(data = CLORPT.rs, diss=NULL, distance = "euclidean",
                                       min.nc = 2, max.nc = 20,
                                       method = "kmeans", index = "ch"))
opt.Clusters.CH
summary(opt.Clusters.CH)
opt.Clusters.CH$Best.nc
plot(search_space, opt.Clusters.CH$All.index,pch=20, col="blue")
lines(search_space, opt.Clusters.CH$All.index)
abline(v=6, col="red", lty=2)

### the silhouette method is also very popular (check here)
set.seed(1984)
system.time(opt.Clusters.silhouette <- NbClust(data = CLORPT.rs, diss=NULL, distance = "euclidean",
                                               min.nc = 2, max.nc = 20,
                                               method = "kmeans", index = "silhouette"))
opt.Clusters.silhouette
summary(opt.Clusters.silhouette)
opt.Clusters.silhouette$Best.nc
plot(search_space, opt.Clusters.silhouette$All.index,pch=20, col="blue")
lines(search_space, opt.Clusters.silhouette$All.index)
abline(v=6, col="red", lty=2)

### My parameters for running the clustering

### Apply k-means with foreach (if you want to create several maps)
# setwd(OutDir)
# detectCores()
# cl <- makeCluster(3)   ### Create cluster
# registerDoParallel(cl)
# getDoParWorkers()
#  K <- c(5,6,12)
# system.time(kmeansCLORPT <- foreach (i=1:length(K),.packages=c("ClusterR"),.export = c("CLORPT.rs","K")) %dopar% {
#   set.seed(1991+i)
#   kmeans_clorpt <- KMeans_rcpp(CLORPT.rs, clusters = K[[i]], num_init = 10, max_iters = 5000,
#                                  fuzzy = TRUE, initializer = 'kmeans++', verbose = T)
#   
#   # Save the object as we go
#   save(kmeans_clorpt, file=paste0(OutDir,"kmeans_clorpt.k",K[[i]],".RData"))
#   kmeans_clorpt # We return this
# })
# 
# stopCluster(cl)

### In this case, we choose 6 classes as optimum

# perform KMeans_rcpp clustering
my_seed <- 4587 ### Your set.seed() number
km.pedogenon.rcpp <-KMeans_rcpp(data=CLORPT.rs, clusters=6,
                                num_init = 20, max_iters = 5000,
                                initializer = "kmeans++", fuzzy = TRUE, verbose = FALSE,
                                seed = my_seed)

### But for the sake of the example we will make more classes than needed
my_seed <- 4587 ### Your set.seed() number
km.pedogenon.rcpp <-KMeans_rcpp(data=CLORPT.rs, clusters=20,
                                num_init = 20, max_iters = 5000,
                                initializer = "kmeans++", fuzzy = TRUE, verbose = FALSE,
                                seed = my_seed)

### save the cluster number in the original dataframe
CLORPT.df$Cluster.N <- as.factor(km.pedogenon.rcpp$clusters)


### Plot the clusters
### Bounding box of your study area

baseMap <- get_stamenmap(bbox = AOI, ### Here you need to change the coordinates with that of your study area
                         maptype="terrain-background", zoom=8,
                         source="stamen", crop=TRUE) 
### Plot the cluster points
library(colorspace)
ggmap(baseMap) + 
  geom_point(aes(x = x, y = y, colour = as.factor(Cluster.N)),data = CLORPT.df) +
  ggtitle("k-means of CLORPT factors, or Genosoils")+
  scale_colour_discrete_qualitative(palette = "Set 3")

### Clean and save image
rm(covariates.cont.sc, tpi.stack,opt.Clusters.silhouette,my_seed, opt_kmeans,
   opt.Clusters.CH, tpi, mrvbf, nrjan, nrjul, rel300m, prescott,slope, search_space, baseMap, tpi.stack, twi, )
#save.image(paste0(OutDir,"/1.Covariate_kmeans.RData"))


# B. Exploring the pedogenon classes ------------------------------------------------

### Load the kmeans models and raster (continue script 1.Covariates_Kmeans.R)
# InputDir <-"C:/Covariates/Output"
# setwd(InputDir)
# load("1.Covariate_kmeans.RData")

# 1. Mapping the clusters and write layers -------------------------------
### Previously, the stack with the scaled covariates were at
covariates.stack
plot(covariates.stack)

### Predict the k class with a nested foreach loop 

### Extract the index of the centroids that are na/nan/Inf
Kcent.nan <- which(apply(km.pedogenon.rcpp[["centroids"]], MARGIN = 1, FUN = function(y){any(is.na(y))}))

### Define the size of the blocks --- At each raster row do we start and finish each crop?
bs <- blockSize(covariates.stack, minblocks=2)
### I crop all raster files (across variables stacks) in a parallel process
### with a %dopar% from the foreach package
### this is thought for a large study area, but of course for the example is not needed

## My desired cluster number
K <- 20
#K <- c(6,12)
### If K is more than one instead of km.pedogenon.rcpp being a kmeans model, it would be a list of models

system.time(
  
  for(m in 1:length(K)){
    
    cl <- makeCluster(3)   
    registerDoParallel(cl)
    
    k_rast_list <- foreach(i=1:bs$n, .packages=c("raster", "ClusterR"), .export = c("C", "km.pedogenon.rcpp")) %dopar% {
      
      ### Get one tile of the raster stack
      tile <- crop(covariates.stack, extent(covariates.stack, bs$row[[i]], bs$row[[i]]+bs$nrows[[i]], 1, ncol(covariates.stack)))
      ### Transform into a dataframe
      tile.df <- as.data.frame(tile, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
      
      ### For each new pixel, I first have to rescale its values
      ## Take out the coordinates
      tile.df.coords <- tile.df[,1:2]
      tile.df <- tile.df[,3:ncol(tile.df)]
      
      # Rescale the data with CLORPT.df (sample from the stack covariates that was used to calibrate the kmeans in the first place)
      tile.df.rs <- as.matrix(tile.df) %*% solve(C)
      tile.df.rs <- as.data.frame(tile.df.rs)
      
      ### Predict cluster assignment
      
      ### Extract the index of the dataframe rows that are na/nan/Inf
      df.na <- which(apply(tile.df.rs, MARGIN = 1, FUN = function(x) {any(is.na(x))}))
      
      ### Create empty prediction column
      tile.df.rs$cluster <- NA
      
      ### If K is more than one instead of km.pedogenon.rcpp being a kmeans model, it would be a list of models
      ### km.pedogenon.rcpp[[m]]$centroids
      ### predict in those rows where there are not na
      tile.df.rs[-df.na, ]$cluster  <- predict_KMeans(data = tile.df.rs[-df.na,1:(ncol(tile.df.rs)-1)], CENTROIDS = km.pedogenon.rcpp$centroids)
      ### Assign the values to a new raster
      k.pred <- setValues(tile[[1]], tile.df.rs$cluster)
      names(k.pred) <- paste0("K",K[[m]])
      k.pred # Return this
    }
    
    stopCluster(cl)
    
    ## Assign function to mosaic
    k_rast_list$na.rm <- TRUE
    k_rast_list$fun <- min
    ## Create mosaic for whole NSW
    k.raster <- do.call(mosaic, k_rast_list)
    names(k.raster) <- paste0("K",K[[m]])
    
    ## Write to file
    #writeRaster(k.raster, filename= paste0("K",K[[m]],".tif"), na.rm=T,inf.rm=T, format="GTiff", overwrite=TRUE )
    gc()
    
  }
)



# ### 2. Visualization ----------------------------------------------------

### for this example, we use the whole study area and a smaller section. 
### Workflow
### Clip map to study area
#zone <- snmaller_area
pedogenon.zone <- crop(k.raster, extent(150.07323,150.15760,-31.36793,-31.32397))
pedogenon <- k.raster
plot(pedogenon)
plot(pedogenon.zone)
pedogenon.3857 <- projectRaster(pedogenon, crs=CRS("+init=EPSG:3857"), method = "ngb")
pedogenon.zone.3857 <- projectRaster(pedogenon.zone, crs=CRS("+init=EPSG:3857"), method = "ngb")


### Map the Whole study area
my_palette <- c("OrYel","PurpOr","TealGrn","BurgYl","Lajolla","Turku","RdPu",
                "GnBu","YlOrRd","Peach",
                "OrRd", "Greens", "Burg", "Heat 2", "Dark Mint",
                "Blues", "SunsetDark", "PuBuGn", "Viridis", "Heat")

hc.whole <- viz.map.legend.hclust(kmodel = km.pedogenon.rcpp)
plot(hc.whole) # 5 branches

### Choose number of branches
viz.branches(hc.whole, 6)
hc.6.out <-viz.map.legend.pal(kmodel = km.pedogenon.rcpp,
                              branchN =6, pal.names = my_palette,
                              kmap = pedogenon.3857,
                              need.proj = FALSE,
                              legend.name = "pedogenon.6_branches")

plot(hc.6.out$legend.plot)
hc.6.out$map.out %>%
  addScaleBar( position = "bottomright",
               options = scaleBarOptions(maxWidth = 200, metric = TRUE, imperial = FALSE, updateWhenIdle = TRUE))

## Create table that includes summary statistics on:
#  - Pedogenon class area - raster
#  - Average distance between pixels of that class - raster
#  - Mahalanobis distance across all clorpton centroids - centroid df
#  - Mahalanobis distance to the closest class - centroid df

K20.Area.df <- pedogenon.area.func(pedogenon, "K20.Area")
K20.CentrMh <- centroid.dist.func(kmodel=km.pedogenon.rcpp, k.area.df =K20.Area.df, fname = "K20.CentrMh")

# Study Area
K20.area.zone.df <- pedogenon.area.func(kmap=pedogenon.zone, fname="K20.area.zone")

### Join tables with pedogenon area present in a study area and Mahalanobis distance to its closest Genosoil (all NSW)
k20.AreaCentrMh.zone <- study.pedogenon.area.func(study.area.df= K20.area.zone.df, 
                                                  LARGE.centroid.dist.area.df =  K20.CentrMh,
                                                  fname= "k20.AreaCentrMh.zone" )
# ### Calculate the % of the large study area occupied by these pedogenons
# total.area <- sum(k6.AreaCentrMh.zone$LARGE_area_Km2, na.rm=TRUE)
# k6.AreaCentrMh.zone$AreaPerc_Large <- k6.AreaCentrMh.zone$LARGE_area_Km2/total.area*100

### Map study area with custom hierarchical cluster color ramp
### Function to calculate dendrogram  only for the genosoils present in the study area 

### Map study area with custom hierarchical cluster color ramp
### Function to calculate dendrogram  only for the genosoils present in the study area
hc.zone <- viz.map.hclust.study.area(kmodel = km.pedogenon.rcpp, study.area.map = pedogenon.zone)
plot(hc.zone) # 4 branches
### Choose number of branches
viz.branches(hc.zone, 4)
hc.zone.out <-viz.map.legend.pal.study.area(kmodel = km.pedogenon.rcpp,
                                            branchN =4, pal.names = my_palette,
                                            study.area.map = pedogenon.zone.3857,
                                            need.proj = FALSE,
                                            legend.name = "k4.zone.legend")
plot(hc.zone.out$legend.plot)
hc.zone.out$map.out %>%
  addScaleBar( position = "bottomright",
               options = scaleBarOptions(maxWidth = 200, metric = TRUE, imperial = FALSE,
                                         updateWhenIdle = TRUE))

### Map the pedogenons present in the study area across the larger area and locate rectangle around the farm
hc.zone.LARGE.out <- pedogenons.inStudy.area.func(LARGE.Geno.map = pedogenon.3857,
                                                  study.area.Geno.out = hc.zone.out,
                                                  need.proj = FALSE)


hc.zone.LARGE.out$map.out %>%
  addRectangles(color = "black", weight = 2,
                lng1=150.1576, lat1=-31.36793,
                lng2=150.0732 , lat2=-31.32397,
                fillColor = "transparent" )


### We may just be interested in dominant classes

### Map the CLORPTons present in the study area across all NSW and locate rectangle around the farm
## For example, only those larger than a certain area
hc.zone.whole.out <- pedogenons.inStudy.area.bigger.func(LARGE.Geno.map = pedogenon.3857,
                                                         study.area.Geno.out = hc.zone.out,
                                                         study.area.df = k20.AreaCentrMh.zone,
                                                         min.area = 3,
                                                         need.proj = FALSE)
hc.zone.whole.out$map.out %>%
  addRectangles(color = "black", weight = 2,
                lng1=150.1576, lat1=-31.36793,
                lng2=150.0732 , lat2=-31.32397,
                fillColor = "transparent" )

#### end of the script