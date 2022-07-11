network_prep<-function(europe,diam){  
  
  library(geodist)
  library(hexbin) # bin data in hexagons  
  library(splancs) # areapl function
  library(geoscale) # plot geological time scale
  library(alphahull) # convehulls
  library(igraph) # load igraph library
  library(data.table) # rbindlist
  library(maptools) # readShapePoly
  library(RColorBrewer) # color comunities
  library(topicmodels) # perplexity
  library(dplyr) # data manipulation
  library(tidyr) # data manipulation

##############################################################################
###     Agregation of ocurrence data    ###
##############################################################################

# data: original data from PBDB

# Set agregation parameters and data
  #Hexagon_inner_diameter <- diam # hexagon inner diameter = diff(xbnds)/xbins
  mincnt <- 1 # minimun number of ocurrences per stage
  threshold.hex <- 1 # minimun number of hexagons per stage

# Agregate data
#x <- data$paleolng
#y <- data$paleolat

#europe2=europe[-which(grepl("( indet.)|( sp.)",europe[,"Genus_Species"])),]

  unique.eu=unique(europe)
  x= unique.eu$LONG
  y=unique.eu$LAT
#Use end point for long and lat of Europe Lat>35, -25<Long<40
  xbnds <- c(-10, 40)
  ybnds <- c(35, 90)
  
  wid=data.frame(long=c(-10,40),lat=c(mean(europe$LAT),mean(europe$LAT)))
  max.width=geodist(wid,measure = "geodesic")[1,2]/1000
  
  xbins <- max.width/diam
  hbin <- hexbin(x, y, xbins=xbins, shape=0.5, xbnds=xbnds,ybnds=ybnds, xlab="Longitude (dd)", ylab="Latitude (dd)", IDs=T)
  unique.eu$cid <- hbin@cID
  unique.eu$cID <- paste("H", hbin@cID, sep="")

# Select stages with a minimun number of hexagons
#stages <- stage$stage_name

# Create hexagones node set
  hex.xy <- hcell2xy(hbin)
  hex.x <- hex.xy $x
  hex.y <- hex.xy $y
  hex.ID <- paste("H", hbin@cell, sep="")
  hex.id <- hbin@cell 
  nodes <- data.frame(x=hex.x, y=hex.y, cell=hex.ID, density= as.vector(table(unique.eu$cID)))
  rownames(nodes) <- hex.ID

# Create unique spatio-temporal cells
  unique.eu$tID <- paste(unique.eu$TIMEUNIT, unique.eu$cID, sep="_")

  #create matrix with locations on the rows and species on the columns, inside the cell put the time unit number:
  M_hex=matrix(0,ncol=length(unique(unique.eu$Genus_Species)),nrow = length(unique(unique.eu$tID)))
  colnames(M_hex)=unique(unique.eu$Genus_Species)
  rownames(M_hex)=unique(unique.eu$tID)
  for(i in 1:nrow(unique.eu)){
    M_hex[rownames(M_hex)==unique.eu$tID[i],colnames(M_hex)==unique.eu$Genus_Species[i]]=unique.eu$level[i]
  }
  source("Build network.R")
  build_network(M_hex,paste0("Network_",diam,".net"))
  return(list(M_hex=M_hex, unique.eu=unique.eu))
}
