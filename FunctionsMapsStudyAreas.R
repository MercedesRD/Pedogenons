####################################################################################################################################
### Explore Genosoil (CLORPTon) data in smaller study areas
### Author: Mercedes roman
### Date: 01/05/2020
### Objectives:
### 1. Table with k-prototypes for that study area (subset)
### 2. Table with Genosoil present, Area (in study area and outside the study area), Closer Genosoil, Mahalanobis distance to this Genosoil
### 3. Present dendrogram, were the present Genosoils are highlightted.
### 4. Represent with new color palette, to improve the differentiation and visibility (subset)
### 5. Overlay Phenosoil layer. Summarize which genosoil classes have remnant genosoil or phenosoils

##Load packages
library(Rtools)
library(rlang)
library(ClusterR)
library(rgdal)
library(gdalUtils)
library(raster)
library(sp)
library(sf)
library(dplyr)
library(tidyverse)
library(ggmap)
library(ggplot2) 
library(viridis) # color palettes
library(scales)
library(rasterVis)
library(lattice)
library(gridExtra)
library(tmap)    # for static and interactive maps
library(leaflet) # for interactive maps
library(mapview) # for interactive maps
library(shiny)   # for web applications
library(foreach)
library(doParallel)
library(geosphere)
library(dendsort)
library(gplots)
library(dendextend)
library(colorspace)

### Functions for examining the study areas

# ### 1. Calculate area for Genosoils -------------------------------------

### Function for calculating the area of any Genosoil map (any k) and return a summary table
### the input are
### kmap - a raster with genosoil classes
### fname - name of the file for saving the table into a csv file
### Returns a dataframe with the area of each genosoil.
genosoil.area.func <- function(kmap, fname) {
  areaNSWpixels <-area(kmap, na.rm=TRUE)
  s <- stack(kmap, areaNSWpixels)
  # K.values <- getValues(kmap)
  # K.area <- getValues(areaNSWpixels)
  k.A <-getValues(s)
  k.A <- as.data.frame(k.A)
  colnames(k.A) <- c("Genosoil", "Area_km2")
  k.A <- k.A[!is.na(k.A$Genosoil),]
  
  area_K_summary <-  k.A %>% 
    group_by(.,as.factor(Genosoil), .drop=TRUE ) %>% ## Group by Genosoils
    summarise(Genosoil_area = sum(Area_km2, na.rm=TRUE)) ### sum area by Genosoil class
  #str(area_K_summary)
  area_K_summary <- as.data.frame(area_K_summary)
  colnames(area_K_summary) <- c("Genosoil", "Area_Km2")
  write.csv(area_K_summary, file=paste0(fname,".csv")) ### Write table to csv file
  return(area_K_summary) ## and return
}

### Function that calculates the distance between entroids of the kmeans model,
### and returns a table with the genosoil, its closer genosoil, the Mahalanobis distance between genosoils (CLORPT)
### and the areas that they occupy in NSW
### Input:
### kmodel - kmeans model
### k.area.df - is the output of the genosoil.area.func function
### fname - name to export the table to csv
centroid.dist.func <- function(kmodel, k.area.df, fname){
  
  ### kmodel is a kmeans model
  ### k.area.df is the output of the genosoil.area.func function
  #kmodel <- kmeans_clorpt18NoU
  # extract the centroids
  K.centroids <- kmodel$centroids
  K.centroids <- as.data.frame(K.centroids)
  K.centroids$Genosoil <- c(1:nrow(K.centroids))
  #rownames(K.centroids) <- c(1:nrow(K.centroids))
  
  ## Is any centroid NA?
  ### Extract the index of the centroids that are na/nan/Inf
  Kcent.nan <- which(apply(K.centroids, MARGIN = 1, FUN = function(x) {any(is.na(x))}))
  
  ### Calculate distance between all centroids
  dist.centroids <- dist(x=K.centroids[,!names(K.centroids) %in% c("Genosoil")],
                         method = "euclidean")
  
  ### Create empty dataframe to store output
  outs <- data.frame(Genosoil=rep(as.integer(NA),nrow(K.centroids)),
                     ClosestG=rep(as.integer(NA),nrow(K.centroids)),
                     Distance=rep(as.double(NA), nrow(K.centroids)))
  
  outs$Genosoil <- K.centroids$Genosoil ### Assign Genosoil
  Gs <- as.numeric(as.character(outs$Genosoil))
  dist.centroids <- as.matrix(dist.centroids)
  ### Calculate distance to the closest Genosoil
  for(i in 1:nrow(outs)){
    min.dist <- sort(dist.centroids[rownames(dist.centroids)[Gs[[i]]],])[2]
    outs[i,"Distance"] <- min.dist
    outs[i,"ClosestG"] <- names(min.dist)
  }
  ### Remember that those Genosoils that don't exist are NA
  outs$ClosestG <- ifelse(outs$Genosoil %in% Kcent.nan, NA, outs$ClosestG )
  outs$Distance <- ifelse(outs$Genosoil %in% Kcent.nan, NA, outs$Distance )
  colnames(outs) <- c("Genosoil", "Closest Genosoil", "Distance")
  outs$Distance <-round(outs$Distance, digits = 3)
  
  ### Join with the Genosoil area
  outs$Genosoil <- as.character(outs$Genosoil)
  k.area.df$Genosoil <- as.character(k.area.df$Genosoil)
  outs <- left_join(outs, k.area.df, by ="Genosoil")
  
  ### Create column with area closest Genosoil
  outs$Geno2.Area <- NA
  if(length(Kcent.nan) >0) {
    G.exists <- c(1:nrow(outs))[-Kcent.nan]
  } else if(length(Kcent.nan) == 0) {
    G.exists <- c(1:nrow(outs))
  }
  
  for(i in 1:length(G.exists)){
    target.G <- outs[outs$Genosoil == G.exists[[i]], ]$`Closest Genosoil`
    target.A <- outs[outs$Genosoil == target.G, ]$Area_Km2
    outs[outs$Genosoil == G.exists[[i]], ]$Geno2.Area <- target.A
  }
  
  colnames(outs) <- c("Genosoil", "Closest.Genosoil", "MahabDist", "Area_Km2", "Closests.Geno.Area_Km2")
  
  write.csv(outs, file=paste0(fname,".csv")) ### Write table to csv file
  return(as.data.frame(outs)) ## and return
}

### function to join a table with the area per Genosoil for a particular study area, 
### which results from applying the function genosoil.area.func, with
### the output from centroid.dist.func for all NSW

### Function to join the table fot the study area and all NSW
study.geno.area.func<- function(study.area.df, NSW.centroid.dist.area.df, fname) {
  ### Change column names in study.area.df
  colnames(study.area.df) <- c("Genosoil", "Study_area_km2")
  study.area.df$Genosoil <- as.character(study.area.df$Genosoil)
  colnames(NSW.centroid.dist.area.df) <- c("Genosoil", "Closest.Genosoil","MahabDist","NSW_area_Km2","Cl.Geno.NSW_area_Km2")
  study.area.df <- left_join(study.area.df, NSW.centroid.dist.area.df, by="Genosoil")
  study.area.df <- study.area.df[,c("Genosoil","Study_area_km2","NSW_area_Km2",
                                    "Closest.Genosoil","MahabDist","Cl.Geno.NSW_area_Km2")]
  study.area.df <- arrange(study.area.df,- Study_area_km2) ### From larger to smaller Genosoil class in the study area
  #head(study.area.df)
  study.area.df$Study_area_km2 <-round(study.area.df$Study_area_km2 , digits = 2)
  study.area.df$NSW_area_Km2 <-round(study.area.df$NSW_area_Km2 , digits = 2)
  study.area.df$Cl.Geno.NSW_area_Km2 <-round(study.area.df$Cl.Geno.NSW_area_Km2 , digits = 2)
  write.csv(study.area.df, file=paste0(fname,".csv")) ### Write table to csv file
  return(study.area.df) ## and return
  
}

# ### 2. Hierarchical clustering of genosoils and color legend ------------

### First, perform the hierarchical clustering and save it to plot
### Input: kmodel - kmeans model
### Output: Hierarchical cluster (ward.D2 distance) of centroids
viz.map.legend.hclust <- function(kmodel) {
  
  ### Extract centroids from model
  centroids <- kmodel$centroids
  ### Extract the index of the centroids that are na/nan/Inf
  Kcent <- as.data.frame(centroids)
  Kcent.nan <- which(apply(Kcent, MARGIN = 1, FUN = function(x) {any(is.na(x))}))
  ### Exclude these clusters from everywhere
  if(length(Kcent.nan) >0) {
    Kcent.exist <- Kcent[-Kcent.nan,]
  } else if(length(Kcent.nan) == 0) {
    Kcent.exist <- Kcent
  }
  # Kcent.exist <- Kcent[-Kcent.nan,]
  ### Hierarchical clustering
  hc <- hclust(dist(Kcent.exist), method="ward.D2")
  plot(dendsort(hc), main="Hierarchical clustering of kmeans centroids", sub="", xlab="")
  return(hc)
  
}

### function to choose the number of branches for color ramps
### Input:
### hc.object - hierarchical cluster, output from viz.map.legend.hclust function
### branchN - number of branches
### Output: a plot with the dendrogram and colored branches
viz.branches <- function(hc.object, branchN) {
  hc.object %>% as.dendrogram(.) %>% color_branches(., k = branchN) %>%
    plot(., main = paste0("Colored ",branchN," branches"))
}

# pal.names <- c("OrYel","PurpOr", "Dark Mint", "BurgYl","Turku",
#                "YlOrRd", "RdPu", "Peach", "GnBu","Lajolla", 
#                "OrRd", "Greens", "Burg", "Heat 2","Blues", 
#                "BuPu")

my_palette <- c("OrYel","PurpOr","TealGrn","BurgYl","RdPu",
                "GnBu","YlOrRd","Peach","Turku","Lajolla",
                "OrRd", "Greens", "Burg", "Heat 2", "Dark Mint",
                "Blues", "SunsetDark", "PuBuGn", "Viridis", "Heat")

### Function to map with the selected color palettes, based on the dedrogram, the Genosoils for NSW
### input:
### kmodel - our kmeans model
### branchN - number of branches
### pal.names - selectio of palettes from colorspace
### legend.name - name for the pdf to plot the legend (dendrogram)
### kmap - raster layer with Genosoils
### Output: 
### $hc - Hierarchical cluster
### $branch.centroids.ord - Table with centroid, branch, and color
### $legend.plot - dendrogram, legend with colors.
### $binpal 
### $map.out - leaflet map
viz.map.legend.pal <- function(kmodel, branchN, pal.names, legend.name, kmap, need.proj){
  
  ### Extract centroids from model
  centroids <- kmodel$centroids
  ### Extract the index of the centroids that are na/nan/Inf
  Kcent <- as.data.frame(centroids)
  Kcent$Genosoil <- c(1:nrow(Kcent))
  Kcent.nan <- which(apply(Kcent, MARGIN = 1, FUN = function(x) {any(is.na(x))}))
  ### Exclude these clusters from everywhere
  if(length(Kcent.nan) >0) {
    Kcent.exist <- Kcent[-Kcent.nan,]
  } else if(length(Kcent.nan) == 0) {
    Kcent.exist <- Kcent
  }
  
  #Kcent.exist <- Kcent[-Kcent.nan,]
  
  ### Perform hierarchichal clustering
  hc <- hclust(dist(Kcent.exist[,!names(Kcent.exist) %in% c("Genosoil")]), method="ward.D2")
  
  ### Extract labels
  hc.labels <- hc %>% as.dendrogram(.) %>% labels %>% as.numeric()
  
  ### Extract the membership from the tree
  dend <- hc %>% as.dendrogram(.)
  Kcent.exist$branch <- dend %>% dendextend:::cutree.dendrogram(., k = branchN)
  
  # branch.centroids <- as.data.frame(cbind(c(as.numeric(as.character(rownames(Kcent.exist)))),
  #                                         as.numeric(as.character(dendextend:::cutree.dendrogram(dend,k = branchN)))))
  
  branch.centroids <- Kcent.exist[,c("Genosoil", "branch")]
  branch.centroids$Genosoil <- as.numeric(branch.centroids$Genosoil)
  branch.centroids$branch <- as.numeric(branch.centroids$branch)
  colnames(branch.centroids) <- c("Centroid", "Branch")
  
  ### sort the dataframe of branch and genosoil by the dendrogram labels
  branch.centroids.ord <- branch.centroids %>% right_join(tibble(Centroid = hc.labels), by = "Centroid")
  numbs.pal <- c((table(Kcent.exist$branch)))
  branch.count <- as.data.frame(cbind(c(1:branchN), numbs.pal))
  colnames(branch.count) <- c("Branch", "Count")
  #branch.count <- branch.count[order(- branch.count$Count),]
  branch.count <- branch.count %>% arrange(., -Count)
  
  ###Assign color to each 
  branch.centroids.ord$colors <- NA
  
  for(i in 1:length(numbs.pal)){
    ## Generate as many colors for each pallete as centroids in the branch
    branch.centroids.ord[branch.centroids.ord$Branch == branch.count[i,]$Branch,]$colors <-
      sequential_hcl(pal.names[[i]], n = branch.count[i,]$Count)
  }
  
  ### Create legend

  legend.plot <- dend %>%  set("labels_col", branch.centroids.ord$colors) %>% 
    set("branches_k_color", branch.centroids.ord$colors)
  
  pdf(file = paste0("Map_legend",legend.name,".pdf"), width = 10, height = 100 )
  plot(legend.plot,
       main = "Hierarchical histogram of genosoil centroids with the map colors",
       horiz = TRUE) # change color 
  dev.off()
  
  ### Now, reorder by Genosoil class
  branch.centroids.ord <- branch.centroids.ord %>% arrange(., Centroid)
  #branch.centroids.ord <- branch.centroids.ord[order(branch.centroids.ord$Centroid),]
  
  ### Create palette for leaflet
  #pal <- branch.centroids.ord$colors
  
  # binpal <- colorBin(palette = branch.centroids.ord$colors,
  #                    bins = c(as.numeric(as.character(rownames(Kcent.exist))),
  #                             tail(as.numeric(as.character(rownames(Kcent.exist))),1)+1),
  #                    na.color = "transparent")
  
  ### Project the map into the leaflet projection
  if (need.proj == TRUE) {
    kmap <- projectRaster(kmap, crs=CRS("+init=EPSG:3857"), method = "ngb")
  } else if (need.proj == FALSE) { 
    kmap <- kmap
    }
  
  
  binpal <- colorBin(palette = branch.centroids.ord$colors,
                     bins = c(branch.centroids.ord$Centroid,
                              tail(branch.centroids.ord$Centroid,1)+1),
                     na.color = "transparent")
  
  map.out <- leaflet() %>%
    # Base groups
    addTiles(group="OSM (default)") %>%
    addProviderTiles("Esri.WorldImagery", group = "World Imagery") %>% # , group = "World Imagery"
    addProviderTiles("OpenTopoMap", group = "Topo Map") %>%
    addRasterImage(kmap, opacity = 1, colors=binpal, project=need.proj, 
                   maxBytes = 300000000, group = "Genosoils") %>%
    fitBounds(lng1=140, lat1=-38, lng2=154, lat2=-28) %>%
    leafem::addMouseCoordinates() %>%
    addLayersControl(
      baseGroups = c("OSM (default)","World Imagery", "Topo Map"),
      overlayGroups = c("Genosoils"),
      options = layersControlOptions(collapsed = FALSE)
    )
  
  #mapshot(map.out, file = paste0(OutDir,"/Map_",legend.name,".pdf"), remove_url = FALSE)
  output <- list("hc"=hc, "branch.centroids.ord"=branch.centroids.ord,
                 "legend.plot"=legend.plot, "map.out"=map.out)
  return(output)
  
}


### Function to calculate dendrogram  only for the genosoils present in the study area
### Input: kmodel - kmeans model
### study.area.map - clip of raster genosoil only for the study area
### Output: Hierarchical cluster (ward.D2 distance) of centroids
viz.map.hclust.study.area <- function(kmodel, study.area.map) {
  ### Extract centroids from model
  K.centroids <- kmodel$centroids
  K.centroids <- as.data.frame(K.centroids)
  K.centroids$Genosoil <- c(1:nrow(K.centroids))
  #rownames(K.centroids) <- c(1:nrow(K.centroids))
  ### Extract the unique values from the genosoil maps
  Unique.Geno.sa <- unique(getValues(study.area.map))
  ## Exclude NA
  Unique.Geno.sa <- Unique.Geno.sa[!is.na(Unique.Geno.sa)]
  ## Extract the index of the centroids that are na/nan/Inf
  ### Exclude these clusters from everywhere
  Kcent.exist <- K.centroids[K.centroids$Genosoil %in% Unique.Geno.sa,]
  ### Hierarchical clustering
  hc <- hclust(dist(Kcent.exist[,!names(Kcent.exist) %in% c("Genosoil")]), method="ward.D2")
  plot(dendsort(hc), main="Hierarchical clustering of kmeans centroids", sub="", xlab="")
  return(hc)
}


### Function to map with the selected color palettes, based on the dedrogram, the Genosoils for the study area
### input:
### kmodel - our kmeans model
### branchN - number of branches
### pal.names - selectio of palettes from colorspace
### legend.name - name for the pdf to plot the legend (dendrogram)
### study.area.map - clip of raster genosoil only for the study area
### Output: 
### $hc - Hierarchical cluster
### $branch.centroids.ord - Table with centroid, branch, and color
### $legend.plot - dendrogram, legend with colors.
### $binpal 
### $map.out - leaflet map
viz.map.legend.pal.study.area <- function(kmodel, branchN, pal.names, study.area.map, legend.name, need.proj){
  
  ### Subset number of palettes
  pal.names <- pal.names[1:branchN]
  
  ### Extract centroids from model
  K.centroids <- kmodel$centroids
  K.centroids <- as.data.frame(K.centroids)
  K.centroids$Genosoil <- c(1:nrow(K.centroids))
  #rownames(K.centroids) <- c(1:nrow(K.centroids))
  ### Extract the unique values from the genosoil maps
  Unique.Geno.sa <- unique(getValues(study.area.map))
  ## Exclude NA
  Unique.Geno.sa <- Unique.Geno.sa[!is.na(Unique.Geno.sa)]
  ## Extract the index of the centroids that are na/nan/Inf
  ### Exclude these clusters from everywhere
  Kcent.exist <- K.centroids[K.centroids$Genosoil %in% Unique.Geno.sa,]
  ### Hierarchical clustering
  hc <- hclust(dist(Kcent.exist[,!names(Kcent.exist) %in% c("Genosoil")]), method="ward.D2")
  
  ### Extract labels
  hc.labels <- hc %>% as.dendrogram(.) %>% labels %>% as.numeric()
  
  ### Extract the membership from the tree
  dend <- hc %>% as.dendrogram(.)
  Kcent.exist$branch <- dend %>% dendextend:::cutree.dendrogram(., k = branchN)
  
  #branch.centroids <- as.data.frame(cbind(c(as.numeric(as.character(rownames(Kcent.exist)))),
  #                                        as.numeric(as.character(dendextend:::cutree.dendrogram(dend,k = branchN)))))
  branch.centroids <- Kcent.exist[,c("Genosoil", "branch")]
  branch.centroids$Genosoil <- as.numeric(branch.centroids$Genosoil)
  branch.centroids$branch <- as.numeric(branch.centroids$branch)
  colnames(branch.centroids) <- c("Centroid", "Branch")
  
  ### sort them by the dendrogram labels
  branch.centroids.ord <- branch.centroids %>% right_join(tibble(Centroid = hc.labels), by = "Centroid")
  numbs.pal <- c((table(Kcent.exist$branch)))
  branch.count <- as.data.frame(cbind(c(1:branchN), numbs.pal))
  colnames(branch.count) <- c("Branch", "Count")
  #branch.count <- branch.count[order(- branch.count$Count),]
  branch.count <- branch.count %>% arrange(., -Count)
  
  ###Assign color to each 
  branch.centroids.ord$colors <- NA
  
  for(i in 1:length(numbs.pal)){
    ## Generate as many colors for each pallete as centroids in the branch
    branch.centroids.ord[branch.centroids.ord$Branch == branch.count[i,]$Branch,]$colors <-
      sequential_hcl(pal.names[[i]], n = branch.count[i,]$Count)
  }
  
  legend.plot <- dend %>%  set("labels_col", branch.centroids.ord$colors) %>% 
    set("branches_k_color", branch.centroids.ord$colors)
  
  pdf(file = paste0("Map_legend",legend.name,".pdf"), width = 10, height = 40 )
  plot(legend.plot,
       main = "Hierarchical histogram of genosoil centroids with the map colors",
       horiz = TRUE) # change color 
  dev.off()
  
  ### Now, reorder by Genosoil class
  branch.centroids.ord <- branch.centroids.ord %>% arrange(., Centroid)
  #branch.centroids.ord <- branch.centroids.ord[order(branch.centroids.ord$Centroid),]
  
  ### Create palette for leaflet
  #pal <- branch.centroids.ord$colors
  ### Project the map into the leaflet projection
  if (need.proj == TRUE) {
    study.area.map <- projectRaster(study.area.map, crs=CRS("+init=EPSG:3857"), method = "ngb")
  } else if (need.proj == FALSE) {
    study.area.map <- study.area.map
    }
  
  study.area.map <- projectRaster(study.area.map, crs=CRS("+init=EPSG:3857"), method = "ngb")
  
  binpal <- colorBin(palette = branch.centroids.ord$colors,
                     bins = c(branch.centroids.ord$Centroid,
                              tail(branch.centroids.ord$Centroid,1)+1),
                     na.color = "transparent")
  
  map.out <- leaflet(options = leafletOptions(zoomControl = FALSE)) %>%
     #Base groups
    addTiles(group="OSM (default)") %>%
    addProviderTiles("Esri.WorldImagery", group = "World Imagery") %>% # , group = "World Imagery"
    addProviderTiles("OpenTopoMap", group = "Topo Map") %>%
    addRasterImage(study.area.map, opacity = 1, colors=binpal, project=FALSE, 
                   layerId = "values", maxBytes = 300000000, group="Genosoils") %>%
    #fitBounds(lng1=140, lat1=-38, lng2=154, lat2=-28) %>%
    leafem::addMouseCoordinates() %>%
    addLayersControl(
      baseGroups = c("OSM (default)","World Imagery", "Topo Map"),
      overlayGroups = c("Genosoils"),
      options = layersControlOptions(collapsed = FALSE)
    )
  
  #mapshot(map.out, file = paste0(OutDir,"/Map_",legend.name,".pdf"), remove_url = FALSE)
  output <- list("hc" = hc, "branch.centroids.ord" = branch.centroids.ord,
                 "legend.plot" =legend.plot,  "map.out" = map.out)
  
  return(output)
  
}


### Function to map the Genosoils present in a study area, their distribution across all NSW
### It works with the output from the function viz.map.legend.pal.study.area
### Input:
### study.area.Geno.out - Output from viz.map.legend.pal.study.area
### nsw.Geno.map - Map for NSW with Genosoil classes
### Output:
### $map.out - leaflet map
### $kmap - masks the genosoil classes not present in the study area
genos.inStudy.area.func <- function(nsw.Geno.map, study.area.Geno.out, need.proj) {
  
  ### Mask all Genosoils not present in the study area
  genos.present <-  study.area.Geno.out$branch.centroids.ord$Centroid ### the genosoil classes present in the study area
  #genos.present <- as.numeric(unlist(genos.present))
  kmap <- trim(calc(nsw.Geno.map, fun= function(x) {ifelse(x %in% genos.present, x, NA)}))
  
  if (need.proj == TRUE) {
    kmap <- projectRaster(kmap, crs=CRS("+init=EPSG:3857"), method = "ngb")
  } else if (need.proj == FALSE) {
    kmap <- kmap
    }
  
  #kmap <- projectRaster(kmap, crs=CRS("+init=EPSG:3857"), method = "ngb")
  #writeRaster(kmap, filename = "test.tif")
  
  ### Create palette for leaflet
  ### Reorder by Genosoil number
  study.area.Geno.out$branch.centroids.ord <- study.area.Geno.out$branch.centroids.ord %>%
     filter(., Centroid %in% genos.present) %>% arrange(., Centroid)
   
  binpal <- colorBin(palette = study.area.Geno.out$branch.centroids.ord$colors,
                     bins = c(study.area.Geno.out$branch.centroids.ord$Centroid,
                              tail(study.area.Geno.out$branch.centroids.ord$Centroid,1)+1),
                     na.color = "transparent")
  
  map.out <- leaflet() %>%
    # Base groups
    #addTiles(group="OSM (default)") %>%
    addProviderTiles("Esri.WorldImagery", group = "World Imagery") %>% # , group = "World Imagery"
    #addProviderTiles("OpenTopoMap", group = "Topo Map") %>%
   # addProviderTiles("Stamen.TonerLite", group = "Stamen.TonerLite") %>%
    addRasterImage(kmap, opacity = 1, colors=binpal, project=FALSE, 
                   maxBytes = 300000000, group="Genosoils")  #%>%
    # fitBounds(lng1=140, lat1=-38, lng2=154, lat2=-28) %>%
    # leafem::addMouseCoordinates() %>%
    # addLayersControl(
    #   baseGroups = c("OSM (default)","World Imagery", "Topo Map","Stamen.TonerLite"),
    #   overlayGroups = c("Genosoils"),
    #   options = layersControlOptions(collapsed = FALSE))

  output <- list("map.out" = map.out, "kmap" = kmap)
  return(output)
  
}


### Variation of previous function. Mapping nly thos genosoils in the study area, but the surface has to be bigger than a certain value
### It works with the output from the function viz.map.legend.pal.study.area
### Input:
### study.area.Geno.out - Output from viz.map.legend.pal.study.area
### nsw.Geno.map - Map for NSW with Genosoil classes
### study.area.df - table with the area and genosoils present in the study area, output from study.geno.area.func
### min.area - minimum area for a genosoil in order to be included in the map
### Output:
### $map.out - leaflet map
### $kmap - masks the genosoil classes not present in the study area 
### $dendro.larger.genos -  dendrogram with the larger classes in color, and the others in grey
genos.inStudy.area.bigger.func <- function(nsw.Geno.map, study.area.Geno.out, study.area.df, min.area,need.proj) {
  
  ### Subset those with an area larger than a certain value
  genos.present  <- study.area.df %>% filter(., Study_area_km2 >= min.area) %>% select(., Genosoil)
  genos.present <- as.numeric(unlist(genos.present))
  #as.vector(unlist(genos.present))
  ### Mask all Genosoils not present in the study area
  #genos.present <-  study.area.Geno.out$branch.centroids.ord$Centroid ### the genosoil classes present in the study area
  kmap <- calc(nsw.Geno.map, fun = function(x){ifelse((x %in% genos.present), yes = x, no =NA)})
  if (need.proj == TRUE) {
    kmap <- projectRaster(kmap, crs=CRS("+init=EPSG:3857"), method = "ngb")
  } else if (need.proj == FALSE) {
    kmap <- kmap
  }
  # kmap <- projectRaster(kmap, crs=CRS("+init=EPSG:3857"), method = "ngb")
  
  ### Put Color only the main (larger than min.area) Genosoils
  
  ### Extract the membership from the tree
  order.desired <- study.area.Geno.out$hc %>% as.dendrogram(.) %>% labels %>% as.numeric()

  study.area.Geno.out$branch.centroids.ord <- study.area.Geno.out$branch.centroids.ord %>%
   right_join(tibble(Centroid = order.desired), by = "Centroid")
  
  ### To put in bold
  study.area.Geno.out$branch.centroids.ord$colors <- ifelse(
    study.area.Geno.out$branch.centroids.ord$Centroid %in% genos.present,
    study.area.Geno.out$branch.centroids.ord$colors,
    "gray85")
  
  study.area.Geno.out$branch.centroids.ord$cex.label <- ifelse(
    study.area.Geno.out$branch.centroids.ord$Centroid %in% genos.present,
    1,
    0.25)
  
  dendro.larger.genos <- study.area.Geno.out$hc %>% as.dendrogram(.) %>% 
    set("labels_col", study.area.Geno.out$branch.centroids.ord$colors) %>% 
    set("branches_k_color", study.area.Geno.out$branch.centroids.ord$colors) %>%
    set("labels_cex", study.area.Geno.out$branch.centroids.ord$cex.label)
  
  plot(dendro.larger.genos)
  
  ### Create palette for leaflet
  
  ### Reorder by Genosoil number
  # centroid.leaflet.pal <- study.area.Geno.out$branch.centroids.ord %>%
  #   filter(., Centroid %in% genos.present)  %>% arrange(., Centroid)
  #pal <- centroid.leaflet.pal$colors ### Extract the colors from previous output
  # binpal <- colorBin(palette = centroid.leaflet.pal$colors,
  #                    bins = c(centroid.leaflet.pal$Centroid,
  #                             tail(centroid.leaflet.pal$Centroid,1)+1),
  #                    na.color = "transparent")
  
  study.area.Geno.out$branch.centroids.ord <- study.area.Geno.out$branch.centroids.ord %>%
    filter(., Centroid %in% genos.present) %>% arrange(., Centroid)
  
  binpal <- colorBin(palette = study.area.Geno.out$branch.centroids.ord$colors,
                     bins = c(study.area.Geno.out$branch.centroids.ord$Centroid,
                              tail(study.area.Geno.out$branch.centroids.ord$Centroid,1)+1),
                     na.color = "transparent")
  
  map.out <- leaflet() %>%
    # Base groups
    addTiles(group="OSM (default)") %>%
    addProviderTiles("Esri.WorldImagery", group = "World Imagery") %>% # , group = "World Imagery"
    addProviderTiles("OpenTopoMap", group = "Topo Map") %>%
    addRasterImage(kmap, opacity = 1,  project=FALSE, colors=binpal,
                   maxBytes = 300000000, group="Genosoils") %>%
    fitBounds(lng1=140, lat1=-38, lng2=154, lat2=-28) %>%
    leafem::addMouseCoordinates() %>%
    addLayersControl(
      baseGroups = c("OSM (default)","World Imagery", "Topo Map"),
      overlayGroups = c("Genosoils"),
      options = layersControlOptions(collapsed = FALSE)
    )
  
  output <- list("map.out" = map.out, "kmap" = kmap, 
                 "dendro.larger.genos"=dendro.larger.genos)
  return(output)
  
}


###End of script