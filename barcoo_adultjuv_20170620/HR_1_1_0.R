args <- (commandArgs(TRUE))
f<-as.numeric(args[1])+1

library(adehabitatHR); library(sp)

# --------- getting mean HR area (integer)
# set dir for turtle outputs
setwd("/vlsci/VR0212/mrke/achaves/abm/results_1_1_0")
file.list<-list.files()
file.list<-file.list[grep("turtles", file.list)]
globalhr<-vector()

  turtles <- read.table(file.list[f], sep="")
  hrpath<-SpatialPointsDataFrame(turtles[1:2], turtles[3]) # creates a spatial points data frame (adehabitatHR package)
  hr<-mcp(hrpath,percent=100) # get homerange

  save(hrpath,file=paste('hrpath_',f,'.Rda',sep=''))