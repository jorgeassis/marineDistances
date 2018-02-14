## ---------------------------------------------
##
## Contour : Pairwise Marine Distances
## Jorge Assis (2016)
## 
## A high-resolution polygon is converted to an infinite resistance surface. 
## Minimum distances between sites are computed with a shortest path algorithm considering the infinite resistance of land and null resistance throughout the sea. 
## The outcomes are a matrix of pairwise distances, a figure to visualize if sites are well represented in the study area and a figure depicting an example of a shortest marine distance.
## The main file with the sites should be structured as “Name Lon Lat” or “Name Lat Lon”. Coordinates must be be decimal degrees.
##
## ---------------------------------------------
## Instructions
##
## 0.   Download a high resolution polygon depicting the surface of the world (e.g., Global Self-consistent Hierarchical High-resolution Shorelines; https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html)
## 1.   Open R and set the working directory (path to)
## 2.   Load the main function "contour" into memory (bellow)
## 3.   Run the  main function "contour"
## 3.1  The main file with the locations should be text delimited
## 3.2  Provide the path of the polygon depicting the surface of the world
## 3.3  Define the delimiter and decimal character of the text file
## 3.4  Define the main file structure: 1 to "Name Lon Lat" or 2 to "Name Lat Lon"
## 3.5  Define if the text file has a header with the column names (TRUE or FALSE)
## 3.6  Define the resolution of the study area and the buffer to use around the sites. 
##      The buffer can be a simple value or a vector such as c(xmin,xmax,ymin,ymax).
## 3.7  Choose to export the results as a text delimited file (TRUE or FALSE)
##
## ---------------------------------------------
## Example
##
## contour(  global.polygon = "World/Global_CostLine_HD_Polygon.shp" ,
##           file = "example.file.txt" , 
##           file.sep = "," ,
##           file.dec = "." ,
##           file.strucutre = 2 , 
##           file.header = FALSE ,
##           resolution = 0.01 ,
##           buffer = c(5,5,5,5) ,
##           export.file = TRUE   )
## 
## ---------------------------------------------
## Citation
## 
## Assis, J., Castilho Coelho, N., Alberto, F., Valero, M., Raimondi, P., Reed, D., … Serrão, E. A. (2013). High and Distinct Range-Edge Genetic Diversity despite Local Bottlenecks. PLoS ONE, 8(7), e68646. https://doi.org/10.1371/journal.pone.0068646
## 
## ------------------------------------------------------------------------------------------------------------------

## ------------------------------------------------------------------------------
## Run the following code to load the main function "contour" into memory
## ------------------------------------------------------------------------------

global.polygon <- "/Volumes/Albacora/Dropbox/Raw Data/Shapefiles/World Present/Global_CostLine_HD_Polygon.shp"
contour( global.polygon = global.polygon , file= "Data/sampling_sites.txt" , file.sep = ";" , file.dec = "." , file.strucutre = 1 , file.header = FALSE , resolution = 0.01 , buffer = c(1,1,1,1) , export.file = TRUE )

## -----------------

contour <- function (global.polygon , file , file.sep , file.dec , file.header, resolution, buffer , file.strucutre,export.file) {
  
  site.names <- as.character(read.table(file,sep=file.sep,dec=file.dec,header = file.header)[,1])
  
  if(file.strucutre == 1) {
    site.xy <- read.table(file,sep=file.sep,dec=file.dec,header = file.header)[,c(2,3)]
  }
  if(file.strucutre == 2) {
    site.xy <- read.table(file,sep=file.sep,dec=file.dec,header = file.header)[,c(3,2)]
  }
  
  ## ---------
  
  if( !"wordcloud" %in% names(installed.packages()[,1]) ) {
    
    cat( "\r" , "\r" , "Installing Dependences.", "\r", "\r")
    suppressMessages(suppressWarnings( try(install.packages("wordcloud" , verbose=FALSE , quiet = TRUE) , silent = TRUE) ) )
    
  }
  if( !"raster" %in% names(installed.packages()[,1]) ) {
    
    cat( "\r" , "\r" , "Installing Dependences.", "\r", "\r")
    suppressMessages(suppressWarnings( try(install.packages("raster" , verbose=FALSE , quiet = TRUE) , silent = TRUE) ) )
    
  }
  if( !"rgeos" %in% names(installed.packages()[,1]) ) {
    
    cat( "\r" , "\r" , "Installing Dependences.", "\r", "\r")
    suppressMessages(suppressWarnings( try(install.packages("rgeos" , verbose=FALSE , quiet = TRUE) , silent = TRUE) ) )
    
  }
  if( !"gdistance" %in% names(installed.packages()[,1]) ) {
    
    cat( "\r" , "\r" , "Installing Dependences.", "\r", "\r")
    suppressMessages(suppressWarnings( try(install.packages("gdistance" , verbose=FALSE , quiet = TRUE) , silent = TRUE) ) )
    
  }
  cat( "\r" , "\r" ,   "                                         ", "\r", "\r")
  cat( " ", "\n")
  
  cat( " " ,   "Loading Dependences.")
  suppressMessages(suppressWarnings(try( library(raster) , silent = TRUE) ))
  suppressMessages(suppressWarnings(try( library(rgeos) , silent = TRUE) ))
  suppressMessages(suppressWarnings(try( library(gdistance) , silent = TRUE) ))
  suppressMessages(suppressWarnings(try( library(wordcloud) , silent = TRUE) ))
  
  ## ---------------------------------------------------------------------------------------------
  
  if(length(buffer) == 1) { buffer <- rep(buffer,4)}
  xmin <- min(site.xy[,1]) - buffer[1]
  xmax <- max(site.xy[,1]) + buffer[2]
  ymin <- min(site.xy[,2]) - buffer[3]
  ymax <- max(site.xy[,2]) + buffer[4]
  
  if (xmin <= -180) { xmin <- -180}
  if (xmax >= 180) { xmin <- 180}
  if (ymin <= -90) { xmin <- -90}
  if (ymax >= 90) { ymax <- 90}
  
  cat( " ", "\n")
  cat( " " , " ", "\n")
  cat( " " , ".....................................................................", "\n")
  cat( " " , "Generating region of interest.","\n")
  cat( " " , "This step may take a while depending on the resolution chosen and distance between sites.")
  
  ocean <- shapefile(global.polygon)
  ocean <- crop(ocean, extent(xmin, xmax, ymin, ymax)) 
  
  region.as.table <- matrix( NA ,nrow= ((ymax-ymin)/resolution) ,ncol= ((xmax-xmin)/resolution) )
  region.as.raster <- raster(region.as.table)
  extent(region.as.raster) <- c(xmin,xmax,ymin,ymax)
  
  ocean.r <- rasterize(ocean, region.as.raster)
  ocean.r[!is.na(ocean.r)] <- 0
  ocean.r[is.na(ocean.r)] <- 1
  
  ## ---------------------------------------------------------------------------------------------
  ## Relocate sites if needed
  
  sites.to.relocate <- extract(ocean.r,site.xy[,1:2]) == 0
  sites.to.relocate.xy <- site.xy[sites.to.relocate,]
  
  if( nrow(sites.to.relocate.xy) > 0 ) {
    
    cat(  " ", "\n")
    cat( " " , ".....................................................................", "\n")
    cat(  " ", "\n")
    
    if(nrow(sites.to.relocate.xy) == 1 ) {
      
      cat( " " , " 1 site to be relocated.")
      
    } 
    if( nrow(sites.to.relocate.xy) > 1 ) {
      
      cat( " " ,  nrow(sites.to.relocate.xy) , "sites to be relocated.")
    }  
    
    ocean.r.sites <- as.data.frame(ocean.r,xy=TRUE)[,1:2]
    ocean.r.sites <- cbind(ocean.r.sites,extract(ocean.r,ocean.r.sites))
    ocean.r.sites <- ocean.r.sites[ocean.r.sites[,3] == 1 , 1:2 ]
    
    for (i in 1:nrow(sites.to.relocate.xy)) {
      
      near.cells <- as.data.frame( sort( spDistsN1( as.matrix(ocean.r.sites), as.matrix(sites.to.relocate.xy[i,1:2]),longlat=TRUE), decreasing = FALSE,index.return = TRUE)) 
      site.xy[ which(sites.to.relocate)[i] ,  ] <- ocean.r.sites[ near.cells[1,2] , 1:2 ]
      
    }
  }
  
  ## -------------------------
  
  pdf(file="Contour _ Study Region.pdf", width=200, height=200) 
  
  plot(ocean.r , col=c("#fffafa","#87ceeb"), legend=FALSE)
  points(site.xy[,1] , site.xy[,2],  col="black",cex=10) 
  textplot(site.xy[,1] , site.xy[,2]  , site.names , cex=15, show.lines=TRUE,new=FALSE)
  
  suppressMessages( dev.off())
  
  ## -------------------------
  
  cost.surface <- ocean.r
  projection(cost.surface) <- CRS("+proj=longlat +datum=WGS84")
  
  raster_tr <- transition(cost.surface, mean, directions=8)
  raster_tr_corrected <- geoCorrection(raster_tr, type="c", multpl=FALSE)
  
  example.sites.to.plot <- sample(1:nrow(site.xy),2,replace=FALSE)
  
  pdf(file="Contour _ Example.pdf", width=200, height=200) 
  plot(ocean.r , col=c("#fffafa","#87ceeb"), legend=FALSE)
  lines( shortestPath(raster_tr_corrected, as.matrix(site.xy[example.sites.to.plot[1],]) , as.matrix(site.xy[example.sites.to.plot[2],]) , output="SpatialLines") )
  points(site.xy[example.sites.to.plot,1] , site.xy[example.sites.to.plot,2],  col="black",cex=10) 
  textplot( site.xy[example.sites.to.plot,1] , site.xy[example.sites.to.plot,2] , site.names[example.sites.to.plot] ,cex=15, show.lines=TRUE, new=FALSE)
  suppressMessages( dev.off())
  
  raster_cosDist <- costDistance(raster_tr_corrected,  as.matrix(site.xy) ,  as.matrix(site.xy) )
  distance.marine <- as.matrix(raster_cosDist) / 1000
  colnames(distance.marine) <- as.factor(site.names)
  rownames(distance.marine) <- as.factor(site.names)
  
  cat( " ", "\n")
  cat( " " , ".....................................................................", "\n")
  cat( " " , "Job is done.", "\n")
  cat( " ", "\n")
  cat( " " , "Comments or doubts?", "\n")
  cat( " " , "Please contact jorgemfa@gmail.com", "\n")
  
  if(export.file) { write.table(distance.marine, file = "Contour _ Pairwise Marine Distances.txt" , quote = FALSE , sep = file.sep, row.names = TRUE, col.names = TRUE, na = "NA", dec = file.dec) }
  
  if(!export.file) { print(distance.marine) }
}

## ------------------------------------------------------------------------------
## End of function
## ------------------------------------------------------------------------------