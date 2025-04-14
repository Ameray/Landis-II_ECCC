###################################################################
###################################################################
################################################################################
## organisation ECCC
### Dominic Cyr /Ameray Abderrahmane
### ForCS output processing 
rm(list = ls())

# Set working directory
wwd <- "D:/Wena_Project/Testing/Landscape/2025-04-14"
setwd(wwd)
wwd <- getwd()
wwd
#####


### fetching outputs
a <- "mixedwood-042-51"
simDir <- paste0("D:/Wena_Project/Testing/Landscape/2025-04-14")#"
simName <- gsub("simPkg_mixedwood-042-51_", "", basename(simDir))
simName
#simDir <- paste0("D:/ForCS - Test/2020-06-11")#"#Montmorency-Hereford"#"D:/ForCS - "
simInfo <- read.csv(paste(simDir, "simInfo.csv", sep = "/"),
                    colClasses=c("simID"="character"))
x <- list.dirs(simDir, full.names = F, recursive = F)

simInfo
simIDs <- simInfo$simID
###################################################################
###################################################################
###################################################################
require(data.table)
require(dplyr)
require(raster)
require(doSNOW)
require(parallel)
require(foreach)

logs <- c("summary",
          "agbAgeClasses",
          "agbTotal",
          "ageMax",
          "FPS"
          ) #"FPS" ,"summary", "agbAgeClasses", "agbAgeClasses", "agbTotal","ageMax",

# ### hereford
# mgmtLevels <- c("1" = "Intensif",
#                 "3" = "Servitude",
#                 "4" = "Nouveau zonage",
#                 "2" = "Conservation")

# ### ForMont
# mgmtLevels <- c("1" = NA,
#                 "3" = NA,
#                 "4" = NA,
#                 "2" = NA)

if("summary" %in% logs) {
    source("D:/Wena_Project/Testing/frqnt_2022-25-main/scripts/fetchHarvestImplementationFnc.R")
}

require(stringr)
clusterN <- 1
#######
cl = makeCluster(clusterN, outfile = "") ##
registerDoSNOW(cl)

file.copy(paste(simDir, "simInfo.csv", sep = "/"),
          paste0("simInfo_", simName, ".csv"), overwrite = T)

simIDs <- simInfo$simID
simIDs <- str_pad(simIDs, width = max(nchar(simIDs)),
                  side = "left", pad = "0")
dirIndex <- which(simIDs  %in% x &
                      simInfo$areaName == a)
dirIndex

outputList <- foreach(i = dirIndex)  %dopar% {

    require(dplyr)
    require(raster)
    require(reshape2)
    require(data.table)
    require(dplyr)
    require(raster)

    output <- list()

    ### sim variables
    
    simID <-  simIDs[i]  
    sDir <-    paste(simDir, simID, sep ="/")
    areaName <- simInfo[i, "areaName"]
    scenario <- simInfo[i, "scenario"]
    mgmtScenario  <- simInfo[i, "mgmt"]
    mgmtScenarioName <- mgmtScenario
    harvest <- simInfo[i,"harvest"]
    ND_scenario<-simInfo[i,"ND_scenario"]
    replicate <- simInfo[i, "replicate"]
    ###
    sppLvls <- read.table(paste(sDir, "species.txt", sep = "/"),
                          skip = 1, comment.char = ">")[,1]
    if(is.factor(sppLvls)) {
        sppLvls <- levels(sppLvls)    
    }
    

    ### fetching landtypes
    landtypes <- raster(paste(sDir, "landtypes.tif", sep = "/"))
    landtypes_RAT <- read.table(paste(sDir, "landtypes.txt", sep = "/"),
                                skip = 1, comment.char = ">")
    landtypes_RAT <- landtypes_RAT[which(landtypes_RAT[,1] %in% c("yes", "y", "Yes", "Y")),]
    
    #studyArea <- raster(paste0("../inputsLandis/studyArea_", a, ".tif"))
    #landtypes[is.na(studyArea)] <- NA
    
    ## fetching lantypes values
    index <- which(!is.na(values(landtypes)))
    ltVal <- values(landtypes)[index]
    XY_lt <- rowColFromCell(landtypes, index)
    colnames(XY_lt) <- c("row", "column")
    XY_lt <- data.frame(XY_lt,
                        ltID = ltVal) #
    
    
    if("summary" %in% logs) {
        if(harvest) {
            ### fetching mgmt areas and harvest implementation table
            mgmtAreas <- raster(paste(sDir, "mgmt-areas.tif", sep = "/"))
            x <- paste(sDir, "biomass-harvest.txt", sep = "/")
            harvImpl <- fetchHarvestImplementation(x) 
        } else {
            mgmtAreas <- landtypes
            mgmtAreas[!is.na(landtypes)] <- 1
        }

    }


    if("agbAgeClasses" %in% logs) {
        #############
        ## computing landtype area size
        areaSize_lt <- data.frame(freq(landtypes))
        areaSize_lt[,"area_ha"] <- areaSize_lt$count * (prod(res(landtypes))/10000)
        areaSize_lt <- data.frame(ltID = areaSize_lt$value,
                                  ltArea_ha = areaSize_lt$area_ha)
        areaSize_lt <- areaSize_lt[complete.cases(areaSize_lt),]

        if(!("summary" %in% logs)) {
            XY <- XY_lt
        }
    }

    if("summary" %in% logs) {
        #############
        ## computing mgmg area size
        areaSize_mgmt <- data.frame(freq(mgmtAreas))
        areaSize_mgmt[,"area_ha"] <- areaSize_mgmt$count * (prod(res(mgmtAreas))/10000)
        areaSize_mgmt <- data.frame(mgmtID = areaSize_mgmt$value,
                                    mgmtArea_ha = areaSize_mgmt$area_ha)
        areaSize_mgmt <- areaSize_mgmt[complete.cases(areaSize_mgmt),]

        #############
        ## XY df
        #############
        ## fetching management areas values
        index <- which(!is.na(values(mgmtAreas)))
        mgmtVal <- values(mgmtAreas)[index]
        XY_mgmt <- rowColFromCell(mgmtAreas, index)
        colnames(XY_mgmt) <- c("row", "column")
        XY_mgmt <- data.frame(XY_mgmt,
                              mgmtID = mgmtVal) #
        if(!("agbAgeClasses" %in% logs)) {
            XY <- XY_mgmt
        }
    }


    if("summary" %in% logs &
       "agbAgeClasses" %in% logs) {
        XY <- merge(XY_lt, XY_mgmt, all.x = T)
    }


    if("summary" %in% logs) {
        ####
        logSummary <- fread(file = paste(sDir, "log_Summary.csv", sep = "/"))
        ### subsetting outputs
        df <- merge(XY, logSummary)
        # 8. Summarize Time>=1 data by mgmtID, Time; attach simulation constants
        df2 <- df %>%
          group_by(mgmtID, Time) %>%
          summarise(
            ABio = mean(ABio),
            BBio = mean(BBio),
            TotalDOM = mean(TotalDOM),
            DelBio = mean(DelBio),
            Turnover = mean(Turnover),
            NetGrowth = mean(NetGrowth),
            NPP = mean(NPP),
            Rh = mean(Rh),
            NEP = mean(NEP),
            NBP = mean(NBP),
            .groups = "drop"
          ) %>%
          mutate(
            simID = as.character(simID),  # ensure simID is character
            areaName = areaName,
            scenario = scenario,
            mgmtScenario = mgmtScenario,
            mgmtScenarioName = mgmtScenarioName,
            ND_scenario = ND_scenario,
            replicate = replicate,
            TotalCS = ABio + BBio + TotalDOM  # total ecosystem carbon storage
          )
        
        # 9. Compute total ecosystem carbon storage at year 0 from dead and living pools
        Dead <- fread(file.path(sDir, "log_Pools.csv"))
        Live <- fread(file.path(sDir, "log_BiomassC.csv"))
        DeadSum <- Dead %>%
          filter(Time == 0) %>%
          group_by(Time, row, column) %>%
          summarise(DeadBiomass = sum(VF_A, VF_B, Fast_A, Fast_B, MED, Slow_A, Slow_B, Sng_Stem, Sng_Oth, Extra),
                    .groups = "drop")
        LiveSum <- Live %>%
          filter(Time == 0) %>%
          group_by(Time, row, column) %>%
          summarise(LiveBiomass = sum(Wood, Leaf, CrsRoot, FineRoot),
                    .groups = "drop")
        
        totalecsystem <- merge(LiveSum, DeadSum, by = c("Time", "row", "column"))
        totalecsystem <- totalecsystem %>%
          mutate(
            TotalBiomass = LiveBiomass + DeadBiomass,
            simID = simID  # use simID from our variable; already character
          )
        totalecsystem <- merge(XY, totalecsystem)
        totalecsystem_year0 <- totalecsystem %>%
          group_by(Time, simID, mgmtID) %>%
          summarise(TotalCS = mean(TotalBiomass), .groups = "drop")
        
        # 10. Combine the baseline (Time=0) rows with the Time>=1 rows (df2)
        df_combined <- bind_rows(totalecsystem_year0, df2) %>%
          arrange(simID, Time)
        
        # 11. Compute flux corrections using the baseline.
        #     Group by simID and mgmtID so that lag() applies per mgmt area.
        dfSummary <- df_combined %>%
          group_by(simID, mgmtID) %>%
          arrange(Time, .by_group = TRUE) %>%
          mutate(
            NBPc = TotalCS - lag(TotalCS, default = first(TotalCS)),  # For time=1, baseline is time=0
            D = NEP - NBP,
            NEPc = NBPc + D,
            NPPc = NEPc + Rh,
            NetGrowthc = NPPc - Turnover,
            NEPa = cumsum(NEPc),
            NBPa = cumsum(NBPc)
          ) %>%  ungroup()
        dfSummary <- dfSummary %>% filter(Time >= 1)   %>%
          group_by(simID, mgmtID) %>%
          arrange(Time, .by_group = TRUE) %>%
          mutate(
            NEPa = cumsum(NEPc),
            NBPa = cumsum(NBPc)
          ) %>%
          ungroup() %>%
          merge(areaSize_mgmt, by = "mgmtID")
        ### tidying up...
        dfSummary <- melt(dfSummary, id.vars = c("simID",
                                                 "areaName", "scenario", "mgmtScenario", "mgmtScenarioName",
                                                 "replicate","ND_scenario", "Time", "mgmtID", "mgmtArea_ha")) %>%
          arrange(simID, Time, mgmtID, variable)
        
        output[["summary"]] <- dfSummary
        
        rm("logSummary", "df")
        gc()
    }

    if("agbAgeClasses"  %in% logs |
       "agbTotal"  %in% logs |
       "ageMax"  %in% logs) {
        ####
        agb <- fread(file = paste(sDir, "log_BiomassC.csv", sep = "/"))
        #### focusing on study area
        XY_studyArea <- XY[!is.na(XY$mgmtID),]
        agb <- agb %>%
            merge(XY_studyArea, all.x = F)

        if("agbAgeClasses"  %in% logs) {
            # first reduce the size of the table (before merging the XY df)
            breaks <- c(seq(0,120, by = 20),999)
            agb[,"ageClass"] <- cut(agb$Age, breaks)


            ###
            ltVals <- landtypes_RAT$V2
            
            ###
            zeroPadDF <- expand.grid(species = unique(agb$species),
                                     ageClass = unique(agb$ageClass),
                                     landtype = ltVals,
                                     Time = unique(agb$Time),
                                     agb_tonnesTotal = NA)

            ### summarizing by landtype
            
            
            agbSummary <- agb %>%
                mutate(landtype = ltID,
                       agb_tonnes = prod(res(landtypes))/10000*
                           2*(Wood + Leaf)/100) %>%
                group_by(landtype, Time, species, ageClass) %>%
                summarise(agb_tonnesTotal = sum(agb_tonnes))


            ### summarizing AGB
            agbSummary <- agbSummary %>%
                merge(zeroPadDF,
                      by = c("species", "ageClass","landtype", "Time"),
                      all.y = T) %>%
                merge(areaSize_lt,
                      by.x = "landtype", by.y = "ltID") %>%
                mutate(agb_tonnesTotal = ifelse(is.na(agb_tonnesTotal.x), 0, agb_tonnesTotal.x),
                       agb_tonnesPerHa = round(agb_tonnesTotal / ltArea_ha, 2))


            ### tidying up
            agbSummary <- data.frame(simID = as.character(simID),
                                     areaName = areaName,
                                     scenario = scenario,
                                     mgmtScenario  = mgmtScenario,
                                     mgmtScenarioName  = mgmtScenarioName,
                                     ND_scenario=ND_scenario,
                                     replicate = replicate,
                                     agbSummary[, c("Time", "landtype", "species", "ageClass",
                                                    "agb_tonnesTotal", "ltArea_ha", "agb_tonnesPerHa")])


            colnames(agbSummary)[which(colnames(agbSummary) == "ltArea_ha")] <- "landtypeArea_ha"

            output[["agbAgeClasses"]] <- agbSummary


        }

        if("agbTotal"  %in% logs) {
            ### summarizing by landtype
            agbTotal <- agb %>%
                mutate(agb_tonnesPerHa = 2*(Wood + Leaf)/100) %>%
                group_by(row, column, Time, species) %>%
                summarise(agb_tonnesPerHa = round(sum(agb_tonnesPerHa), 2))

            agbTotal$species <- factor(agbTotal$species, levels = sppLvls)

            agbTotal <- data.frame(simID = as.character(simID),
                                   areaName = areaName,
                                   scenario = scenario,
                                   mgmtScenario  = mgmtScenario,
                                   mgmtScenarioName  = mgmtScenarioName,
                                   ND_scenario=ND_scenario,
                                   replicate = replicate,
                                   agbTotal)

            save(agbTotal, file = paste0("agbTotal_", simName, "_", simID, ".RData"))
            rm(agbTotal)
        }

        if("ageMax" %in% logs) {
            ageMax <- agb %>%
                group_by(row, column, Time) %>%
                summarise(ageMax = max(Age))

            ageMax <- data.frame(simID = as.character(simID),
                                 areaName = areaName,
                                 scenario = scenario,
                                 mgmtScenario  = mgmtScenario,
                                 ND_scenario=ND_scenario,
                                 replicate = replicate,
                                 ageMax)
            save(ageMax, file = paste0("ageMax_", simName, "_", simID, ".RData"))
            rm(ageMax)
        }

        rm(agb)
        
    }
    
    if("FPS"  %in% logs
       & harvest) {
        ### fetching targetted mgmt-areas
        ### fetching landtypes
        mgmtAreas <- raster(paste(sDir, "mgmt-areas.tif", sep = "/"))
        # if(a %in% c("Hereford", "ForMont")) {
        #     mgmtAreas <- mgmtAreas >= 10000  
        # } else {
        #     mgmtAreas[!is.na(mgmtAreas)] <- 1
        # }
        
        
        r <- mgmtAreas
        r[] <- 1:ncell(r)
        xy <- as.data.frame(zonal(mgmtAreas, r))
        xy <- cbind(xy,
                    row = rowFromCell(mgmtAreas,xy$zone),
                    column = colFromCell(mgmtAreas,xy$zone)) %>%
            filter(mean == 1)
        xy <- xy[, c("row", "column")]
        totalArea <- as.data.frame(zonal(mgmtAreas, mgmtAreas, sum))
        totalArea <- filter(totalArea, zone == 1)[,2] * prod(res(mgmtAreas))/10000
        
        ####
        FluxBio <- fread(file = paste(sDir, "log_FluxBio.csv", sep = "/"))
        ## correcting error in file format
        #cNames <- c(colnames(FluxBio)[-1], "V1")
        #colnames(FluxBio) <- cNames
        #FluxBio <- FluxBio[,c(1:17)]
        
        ## focusing on targetted area
        FluxBio <- merge(FluxBio, xy, all.y = F)
        ### summarizing FPS
        toFPS <- FluxBio %>%
            group_by(Time, species) %>%
            summarize(BioToFPS_tonnesCTotal = round(prod(res(mgmtAreas))/10000 * sum(BioToFPS)/100,2)) %>%
            mutate(areaManagedTotal_ha = totalArea)
        
        areaHarvested <-  FluxBio %>%
            filter(BioToFPS > 0) %>%
            distinct(Time, row, column) %>%
            group_by(Time) %>%
            summarise(areaHarvestedTotal_ha = prod(res(mgmtAreas))/10000 * n())

        if(nrow(areaHarvested) > 0) {
            toFPS <- merge(toFPS, areaHarvested)
        } else {
            toFPS[,"areaHarvestedTotal_ha"] <- 0
        }
        
        
        
        toFPS$species <- factor(toFPS$species, levels = sppLvls)
        
        toFPS <- data.frame(simID = as.character(simID),
                            areaName = areaName,
                            scenario = scenario,
                            mgmtScenario  = mgmtScenario,
                            mgmtScenarioName = mgmtScenarioName,
                            ND_scenario=ND_scenario,
                            replicate = replicate,
                            toFPS)
        
        output[["FPS"]] <- toFPS
        rm(toFPS)
        
    }
    
    print(paste('Done with simulation', simID))
    return(output)
    
}
stopCluster(cl)

if("summary"  %in% logs ) {
    ### summary
    outputSummary <- list()
    for(i in seq_along(outputList)) {
        outputSummary[[i]] <- outputList[[i]][["summary"]]
        
    }
    outputSummary <-do.call("rbind", outputSummary)
    save(outputSummary, file = paste0("output_summary_", simName, ".RData"))
}

if("agbAgeClasses"  %in% logs ) {
    ### agbAgeClasses
    output_agbAgeClasses <- list()
    for(i in seq_along(outputList)) {
        output_agbAgeClasses[[i]] <- outputList[[i]][["agbAgeClasses"]]
        
    }
    output_agbAgeClasses <- do.call("rbind", output_agbAgeClasses)
    save(output_agbAgeClasses, file = paste0("output_bio_", simName, ".RData"))
}
if("FPS"  %in% logs ) {
    ### summary
    outputSummary <- list()
    for(i in seq_along(outputList)) {
        outputSummary[[i]] <- outputList[[i]][["FPS"]]
        
    }
    outputSummary <-do.call("rbind", outputSummary)
    write.csv(outputSummary, file = paste0("output_BioToFPS_", simName, ".csv"), row.names = F)
    #save(outputSummary, file = paste0("output_BioToFPS_", a, ".RData"))
}










################ test script for simID=1
# ----------------------------

# 1. Load libraries
library(dplyr)
library(data.table)
library(raster)
library(reshape2)

# 2. Set working directory
setwd("D:/Wena_Project/Testing/Landscape/2025-04-08")
wwd <- getwd()
cat("Working directory:", wwd, "\n")

# 3. Read simulation information and select simID = "1"
simDir <- "D:/Wena_Project/Testing/Landscape/2025-04-08"
simInfo <- read.csv(file.path(simDir, "simInfo.csv"), colClasses = c("simID" = "character"))
# Subset the row for simID 1
simInfo_1 <- simInfo %>% filter(simID == "1")
if(nrow(simInfo_1) == 0) stop("simID 1 not found in simInfo.csv")
areaName <- simInfo_1$areaName[1]
scenario <- simInfo_1$scenario[1]
mgmtScenario <- simInfo_1$mgmt[1]
mgmtScenarioName <- mgmtScenario
harvest <- simInfo_1$harvest[1]
ND_scenario <- simInfo_1$ND_scenario[1]
replicate <- simInfo_1$replicate[1]

# For simID 1, define its simulation directory
simID <- "1"
sDir <- file.path(simDir, simID)

# 4. Read species levels (if needed)
sppLvls <- read.table(file.path(sDir, "species.txt"), skip = 1, comment.char = ">")[,1]
if (is.factor(sppLvls)) {
  sppLvls <- levels(sppLvls)
}

# 5. Fetch landtypes and derive spatial (XY) information

# Load landtypes and its RAT file
landtypes <- raster(file.path(sDir, "landtypes.tif"))
landtypes_RAT <- read.table(file.path(sDir, "landtypes.txt"), skip = 1, comment.char = ">")
landtypes_RAT <- landtypes_RAT[landtypes_RAT[,1] %in% c("yes", "y", "Yes", "Y"), ]

# Get XY coordinates from landtypes (for cells with data)
idx <- which(!is.na(values(landtypes)))
ltVal <- values(landtypes)[idx]
XY_lt <- rowColFromCell(landtypes, idx)
colnames(XY_lt) <- c("row", "column")
XY_lt <- data.frame(XY_lt, ltID = ltVal)

# 6. Get management area raster
if (harvest) {
  mgmtAreas <- raster(file.path(sDir, "mgmt-areas.tif"))
  # If needed, you could call fetchHarvestImplementation() here.
} else {
  mgmtAreas <- landtypes
  mgmtAreas[!is.na(landtypes)] <- 1
}

# Compute area size of management areas
areaSize_mgmt <- data.frame(freq(mgmtAreas))
areaSize_mgmt[,"area_ha"] <- areaSize_mgmt$count * (prod(res(mgmtAreas)) / 10000)
areaSize_mgmt <- data.frame(mgmtID = areaSize_mgmt$value, mgmtArea_ha = areaSize_mgmt$area_ha)
areaSize_mgmt <- areaSize_mgmt[complete.cases(areaSize_mgmt), ]

# Derive XY for mgmtAreas
idx_mgmt <- which(!is.na(values(mgmtAreas)))
mgmtVal <- values(mgmtAreas)[idx_mgmt]
XY_mgmt <- rowColFromCell(mgmtAreas, idx_mgmt)
colnames(XY_mgmt) <- c("row", "column")
XY_mgmt <- data.frame(XY_mgmt, mgmtID = mgmtVal)

# For summary processing, we use the management areas XY:
XY <- XY_mgmt

# 7. Read log_Summary.csv and merge with XY
logSummary <- fread(file.path(sDir, "log_Summary.csv"))
df <- merge(XY, logSummary)

# 8. Summarize Time>=1 data by mgmtID, Time; attach simulation constants
df2 <- df %>%
  group_by(mgmtID, Time) %>%
  summarise(
    ABio = mean(ABio),
    BBio = mean(BBio),
    TotalDOM = mean(TotalDOM),
    DelBio = mean(DelBio),
    Turnover = mean(Turnover),
    NetGrowth = mean(NetGrowth),
    NPP = mean(NPP),
    Rh = mean(Rh),
    NEP = mean(NEP),
    NBP = mean(NBP),
    .groups = "drop"
  ) %>%
  mutate(
    simID = as.character(simID),  # ensure simID is character
    areaName = areaName,
    scenario = scenario,
    mgmtScenario = mgmtScenario,
    mgmtScenarioName = mgmtScenarioName,
    ND_scenario = ND_scenario,
    replicate = replicate,
    TotalCS = ABio + BBio + TotalDOM  # total ecosystem carbon storage
  )

# 9. Compute total ecosystem carbon storage at year 0 from dead and living pools
Dead <- fread(file.path(sDir, "log_Pools.csv"))
Live <- fread(file.path(sDir, "log_BiomassC.csv"))
DeadSum <- Dead %>%
  filter(Time == 0) %>%
  group_by(Time, row, column) %>%
  summarise(DeadBiomass = sum(VF_A, VF_B, Fast_A, Fast_B, MED, Slow_A, Slow_B, Sng_Stem, Sng_Oth, Extra),
            .groups = "drop")
LiveSum <- Live %>%
  filter(Time == 0) %>%
  group_by(Time, row, column) %>%
  summarise(LiveBiomass = sum(Wood, Leaf, CrsRoot, FineRoot),
            .groups = "drop")

totalecsystem <- merge(LiveSum, DeadSum, by = c("Time", "row", "column"))
totalecsystem <- totalecsystem %>%
  mutate(
    TotalBiomass = LiveBiomass + DeadBiomass,
    simID = simID  # use simID from our variable; already character
  )
totalecsystem <- merge(XY, totalecsystem)
totalecsystem_year0 <- totalecsystem %>%
  group_by(Time, simID, mgmtID) %>%
  summarise(TotalCS = mean(TotalBiomass), .groups = "drop")

# 10. Combine the baseline (Time=0) rows with the Time>=1 rows (df2)
df_combined <- bind_rows(totalecsystem_year0, df2) %>%
  arrange(simID, Time)

# 11. Compute flux corrections using the baseline.
#     Group by simID and mgmtID so that lag() applies per mgmt area.
dfSummary <- df_combined %>%
  group_by(simID, mgmtID) %>%
  arrange(Time, .by_group = TRUE) %>%
  mutate(
    NBPc = TotalCS - lag(TotalCS, default = first(TotalCS)),  # For time=1, baseline is time=0
    D = NEP - NBP,
    NEPc = NBPc + D,
    NPPc = NEPc + Rh,
    NetGrowthc = NPPc - Turnover,
    NEPa = cumsum(NEPc),
    NBPa = cumsum(NBPc)
  ) %>%  ungroup()
dfSummary <- dfSummary %>% filter(Time >= 1)   %>%
    group_by(simID, mgmtID) %>%
    arrange(Time, .by_group = TRUE) %>%
    mutate(
      NEPa = cumsum(NEPc),
      NBPa = cumsum(NBPc)
    ) %>%
    ungroup() %>%
    merge(areaSize_mgmt, by = "mgmtID")
# 12. Reshape the output using melt
dfSummary <- melt(dfSummary,
                  id.vars = c("simID", "areaName", "scenario", "mgmtScenario",
                              "mgmtScenarioName", "replicate", "ND_scenario",
                              "Time", "mgmtID", "mgmtArea_ha")) %>%
  arrange(simID, Time, mgmtID, variable)

# 13. Print the output for simID = 1
print(dfSummary)










unique(dfSummary$ND_scenario)
library(ggplot2)

dfSummary_year2<-dfSummary %>% filter (Time==2)
dfSummary_year2
# Specify the flux variables you want to plot (adjust as needed)
flux_vars <- c("NEP", "NEPc")

# Filter for the relevant flux variables and ensure variable is a factor (optional, for ordering)
df_flux <- dfSummary %>%
  filter(variable %in% flux_vars) %>%
  mutate(variable = factor(variable, levels = flux_vars))

# Create the plot
ggplot(df_flux, aes(x = Time, y = value, color = variable)) +
  geom_line() +
  facet_grid(~ simID) +
  labs(
    title = "Corrected Fluxes Over Time by ND_scenario",
    x = "Time",
    y = "Flux Value",
    color = "Flux Variable"
  ) +
  theme_minimal()
