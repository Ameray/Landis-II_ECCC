################################################################################
## organisation ECCC
### Dominic Cyr/ Ameray Abderrahmane

#### read BDAthe outputs and calculate severity 
###################################################################
### ForCS BDA output processing
### IMPORTANT: For this script to work, disturbances must be
### simulated in the specified order.
###################################################################


# Clear environment
rm(list = ls())

require(data.table)
require(dplyr)
require(raster)
require(stringr)


dir_path <- "D:/Wena_Project/Testing/2025-04-08"
setwd(dir_path)
#### Since we're controlling the bda parameters, I've enhanced fetchBDAParamFnc so use fetchBDAParam1
##### check  fetchBDAParamFnc.R in BDA_clibaration folder
source("D:/Wena_Project/scripts/fetchBDAParamFnc.R")

### fetching outputs
a <- "mixedwood-042-51"
d <- c("bda" = 4, "wind" = 3, "harvest" = 2, "fire" = 1)  # IMPORTANT: put in the order they were simulated in LANDIS
simDir <- "D:/Wena_Project/Testing/2025-04-08"
simName <- gsub("simPkg_mixedwood-042-51_", "", basename(simDir))
simInfo <- read.csv(paste(simDir, "simInfo.csv", sep = "/"), colClasses = c("simID" = "character"))
simInfo


#  directory indexing
simIDs <- simInfo$simID
simIDs <- str_pad(simIDs, width = max(nchar(simIDs)), side = "left", pad = "0")
allDirs <- list.dirs(simDir, full.names = FALSE, recursive = FALSE)
dirIndex <- which(simIDs %in% allDirs & simInfo$areaName == a)


print(simIDs)    # Check the  simIDs
print(allDirs)   # Check the actual directory names
print(dirIndex)  # Ensure valid indices are returned

outputList <- list()  # Initialize the list for storing outputs
#### process one simulation at a time, faster, less memory consumption
for(i in dirIndex) {
  output <- list()
  ### sim variables
  simID <- simIDs[i]
  simID
  sDir <- paste(simDir, simID, sep ="/")
  areaName <- simInfo[i, "areaName"]
  scenario <- simInfo[i, "scenario"]
  mgmtScenario <- simInfo[i, "mgmt"]
  mgmtScenarioName <- mgmtScenario
  harvest <- simInfo[i, "harvest"]
  BDA<-simInfo[i,"BDA"]
  BDApara <-simInfo[i,"BDAfile"]   ## this column control the number of combinations of parameters tested 
  replicate <- simInfo[i, "replicate"]
  
  ### fetching species
  sppLvls <- read.table(paste(sDir, "species.txt", sep = "/"), skip = 1, comment.char = ">")[,1]
  if(is.factor(sppLvls)) {
    sppLvls <- levels(sppLvls)    
  }
  
  ### fetching landtypes
  landtypes <- raster(paste(sDir, "landtypes.tif", sep = "/"))
  landtypes_RAT <- read.table(paste(sDir, "landtypes.txt", sep = "/"), skip = 1, comment.char = ">")
  landtypes_RAT <- landtypes_RAT[which(landtypes_RAT[,1] %in% c("yes", "y", "Yes", "Y")),]
  
  # fetching lantypes values
  index <- which(!is.na(values(landtypes)))
  ltVal <- values(landtypes)[index]
  XY_lt <- rowColFromCell(landtypes, index)
  colnames(XY_lt) <- c("row", "column")
  XY_lt <- data.frame(XY_lt,
                      ltID = ltVal) #
  
  ### fetching mgmt areas and harvest implementation table
  if(harvest) {
    mgmtAreas <- raster(paste(sDir, "mgmt-areas.tif", sep = "/"))
    x <- paste(sDir, "biomass-harvest.txt", sep = "/")
  } else {
    mgmtAreas <- landtypes
    mgmtAreas[!is.na(landtypes)] <- 1
  }
  ## fetching management areas values
  index <- which(!is.na(values(mgmtAreas)))
  mgmtVal <- values(mgmtAreas)[index]
  XY_mgmt <- rowColFromCell(mgmtAreas, index)
  colnames(XY_mgmt) <- c("row", "column")
  XY_mgmt <- data.frame(XY_mgmt,
                        mgmtID = mgmtVal) #
  
  XY <- merge(XY_lt, XY_mgmt, all.x = T)
  
  ### Biomass processing
  agb <- fread(file = paste(sDir, "log_BiomassC.csv", sep = "/"))
  agb <- agb %>% merge(XY, all.x = F)
  #### here we used  fetchBDAParam1 
  if (BDA) {
  bdaParam <- fetchBDAParam1(paste(sDir, "Base-BDA_budworm1.txt", sep = "/"))
  bdaSpp <- bdaParam$BDASpeciesParameters[,"SpeciesName"]
  
  #### Disturbances
  d <- c( "fire" = 1, "harvest" = 2, "wind" = 3, "bda" = 4)
  ###
  
  ### Disturbance type ID: 1=fire, 2=harvest, 3=wind, 4=bda
  fluxBio <- fread(file = paste(sDir, "log_FluxBio.csv", sep = "/"))
  
  df <- fluxBio %>%
    mutate(AGBtoDOM_woody = MERCH_ToDOM + OtherWoody_ToDOM,
           FOL_ToDOM = FOL_ToDOM,
           BGBtoDOM = CrsRt_ToDOM + FRt_ToDOM,
           bioToAir = MERCH_ToAir + FOL_ToAir + OtherWoody_ToAir + CrsRt_ToAir + FRt_ToAir) %>%
    group_by(Time, row, column, ecoregion, species, Dist) %>%
    summarise(AGBtoDOM_woody = sum(AGBtoDOM_woody),
              FOL_ToDOM = sum(FOL_ToDOM),
              BGBtoDOM = sum(BGBtoDOM),
              bioToAir =  sum(bioToAir),
              bioToFPS = sum(BioToFPS))
  
  ## add a line for nonHosts for later merge
  nonHostDF <- df  %>%
    ungroup() %>%
    filter(Dist == 4) %>%
    dplyr::select(Time, row, column, ecoregion, Dist) %>%
    distinct() %>%
    mutate(species = "nonHosts", Dist = 4,
           AGBtoDOM_woody = 0, FOL_ToDOM = 0, BGBtoDOM = 0,
           bioToAir= 0,  bioToFPS = 0)
  
  df <- rbind(df, nonHostDF) %>%
    arrange(Time, row, column, species, Dist)
  
  ## fetching pre dist agb
  bio <- agb %>%
    mutate(Time = Time + 1,
           #AGB = Wood + Leaf,
           BGB = CrsRoot + FineRoot,
           species = ifelse(species %in% bdaSpp, species, "nonHosts")) %>%
    group_by(Time, row, column, ecoregion, mgmtID, species) %>%
    summarise(AGB_woody = sum(Wood),
              Leaf = sum(Leaf),
              BGB = sum(BGB))
  
  df2 <- filter(df, Dist == 4) %>%
    mutate(Dist = names(d[which(d == 4)])) %>%
    merge(bio) %>%
    dplyr::select(Time, row, column, ecoregion, mgmtID, Dist, species,
                  AGB_woody, Leaf, BGB,
                  AGBtoDOM_woody, FOL_ToDOM, BGBtoDOM, bioToAir, bioToFPS) %>%
    arrange(Time, row, column, species)
  
  df2 <-  data.frame(areaName = areaName,
                     simID = simID,
                     scenario = scenario,
                     mgmtScenario = mgmtScenario,
                     mgmtScenarioName = mgmtScenarioName,
                     BDApara=BDApara,
                     replicate =  replicate,
                     df2)
  
  output[["fluxesBDA"]] <- df2
  outputList[[i]] <- output
  }
  gc()
  
}

### Final aggregation
fluxesBDA <- do.call("rbind", lapply(outputList, function(item) item[["fluxesBDA"]]))
save(fluxesBDA, file = paste0("output_fluxesBDA_", simName, ".RData"))
fluxesBDA














