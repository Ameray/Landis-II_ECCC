################################################################################
## organisation ECCC
### Dominic Cyr/ Ameray Abderrahmane

### this script is only for calibration of BDA Model   
##### After sensitivity analysis, species parameters were identified as the most sensitive in the BDA model. 
###These parameters have been calibrated to match observed severity and  mimic the area affected by mortality.
#### we test more than 100 simulations the following set of parameters is the best one we found,
################################################################################
### we compared to the historic so the period of the outbreak was fixed  as well as the duration =7 years
# Clear the environment
rm(list = ls())

# Define the base directory
base_dir <- "D:/Wena_Project/Testing"
setwd(base_dir)

# Create a working directory for today's date
wwd <- paste(getwd(), Sys.Date(), sep = "/")
dir.create(wwd)
setwd(wwd)
wwd
###wwd Input path (LANDIS input in your laptop)
inputDir <- "D:/Wena_Project/Testing/frqnt_2022-25-main/inputsLandis"

### Load required libraries
require(stringr)
require(dplyr)
require(raster)
require(parallel)
require(doSNOW)


forCSVersion <- "3.1"
### Define simulation parameters
simDuration <- 50  
smoothAgeClasses <- TRUE
includeSnags <- FALSE

expDesign <- list(
  area = c("mixedwood-042-51"),
  scenario = c("baseline"),
  mgmt = list("mixedwood-042-51" = c("noHarvest")), # only for calibration  compare with noHarvest: Harvesting must be deactivated to match Migue's data set
  spinup = FALSE,
  cropped = list("mixedwood-042-51" = TRUE),
  
  fire = FALSE,     ### only for calibration FALSE and comparision
  BDA = TRUE,     # 
  wind = TRUE,    # Wind disturbance is  included.
  harvest = TRUE, 
  rep = 2          # Number of replicates for the simulation.
)
# Define periods for simulation only for calibration we fixed the duration  to be similar to the historic  (7 years)
### outbreak cycle is 32 years
### The outbreak cycle is 32 years. The last observed mortality in the landscape was in 1988, and we expected a return after 32 years in 2020. Our simulation starts in 2020.
### Therefore, we expect an outbreak in Year 1 (2020), lasting 7 years.
### we focus on the first outbreak to be compared to the historical :
periods <- list(c(1, 7))#,c(32,39),c(64, 71), c(96, 100))  ###  
# Base parameters for BDA files, applied to all files
base_params <- expand.grid(
  OutbreakPattern = "CyclicUniform",     # Fixed value
  MaxInterval = 0,                      # Fixed value
  MinInterval = 0,                      # Fixed value
  TimeSinceLastEpidemic = 0,            # Fixed value
  TemporalType = c("variablepulse"),             # Discrete values
  MinROS = c(0),                         # Integer values
  MaxROS = c(3),                         # Integer values linked to MinROS
  Dispersal = c("yes"),                   # Discrete values
  DispersalRate = c(20000),              # Continuous range
  EpidemicThresh = c(0.5),   # Continuous range
  InitialEpicenterNum = c(1),    # Integer range
  OutbreakEpicenterCoeff = c(0.01), # Continuous range
  OutbreakEpicenterThresh = c(1),  # Integer range
  SeedEpicenter = c("yes"),              # Discrete values
  SeedEpicenterMax = c(3),           # Integer range
  SeedEpicenterCoeff = c(0.1),      # Continuous range
  DispersalTemplate = c("MaxRadius"), # Discrete values
  NeighborFlag = c("yes"),               # Discrete values
  NeighborSpeedUp = "none",              # Fixed value
  NeighborRadius = c(1000),        # Continuous range
  NeighborShape = c("uniform"),          # Discrete values
  NeighborWeight = c(1),                 # Continuous range
  IntensityClass2_BDP = c(0.33),     # Continuous range
  IntensityClass3_BDP = c(0.67),          # Continuous range
  Age1 = 0,                              # Fixed value
  Age2 = c(30),                      # Integer range
  Age3 = c(50),                      # Integer range linked to Age2
  Age4 =c(999),  # only for PICE.MAR 
  #ABIE.BAL   prob(i,j) i for class, j for species
  SRDProb11 = c(0.25),                     # Continuous range
  SRDProb21= c(0.5),                     # Continuous range
  SRDProb31= c(1),                     # Continuous range
  VulnProb11 = c(0.2),                    # Continuous range
  VulnProb21 = c(0.42),                    # Continuous range
  VulnProb31 = c(0.85),                    # Continuous range
  #PICE.GLA   prob(i,j) i for class, j for species
  SRDProb12 = c(0.18),                     # Continuous range
  SRDProb22= c(0.36),                     # Continuous range
  SRDProb32= c(0.5),                     # Continuous range
  VulnProb12 = c(0.1),                    # Continuous range
  VulnProb22 = c(0.15),                    # Continuous range
  VulnProb32 = c(0.3),                    # Continuous range
  #PICE.MAR   prob(i,j) i for class, j for species
  SRDProb13 = c(0.07),                     # Continuous range
  SRDProb23= c(0.14),                     # Continuous range
  SRDProb33= c(0.28),                     # Continuous range
  VulnProb13 = c(0.13),                    # Continuous range
  VulnProb23 = c(0.2),                    # Continuous range
  VulnProb33 = c(1)                  # Continuous range
  
)
#### 
nrow(base_params) ### the number of combination 
# Generate BDA files for each period and each parameter combination
bda_input_files <- c()

# Generate BDA files for each period and each parameter combination
for (i in seq_along(periods)) {
  period <- periods[[i]]
  for (j in 1:nrow(base_params)) {
    param <- base_params[j, ]
    active <- i %% 1 == 0  # Active phases
    
    # Filename includes period and parameter index
    filename <- sprintf("Base-BDA_budworm%d-%d.txt", i, j)
    filePath <- file.path(wwd, filename)
    fileConn <- file(filePath, open = "w")  # open write connection
    
    # Construct species data based on phase (active/inactive)
    species_data <- if (active) {
      c(
        paste0("ABIE.BAL ", as.integer(param$Age1), " ", param$SRDProb11, " ", as.integer(param$Age2), " ", param$SRDProb21, 
               " ", as.integer(param$Age3), " ", param$SRDProb31, " ", as.integer(param$Age1), " ", param$VulnProb11, 
               " ", as.integer(param$Age2), " ", param$VulnProb21, " ", as.integer(param$Age3), " ", param$VulnProb31, " yes"),
        paste0("PICE.GLA ", as.integer(param$Age1), " ", param$SRDProb12, " ", as.integer(param$Age2), " ", param$SRDProb22, 
               " ", as.integer(param$Age3), " ", param$SRDProb32, " ", as.integer(param$Age1), " ", param$VulnProb12, 
               " ", as.integer(param$Age2), " ", param$VulnProb22, " ", as.integer(param$Age3), " ", param$VulnProb32, " yes"),
        paste0("PICE.MAR ", as.integer(param$Age1), " ", param$SRDProb13, " ", as.integer(param$Age2), " ", param$SRDProb23, 
               " ", as.integer(param$Age3), " ", param$SRDProb33, " ", as.integer(param$Age1), " ", param$VulnProb13, 
               " ", as.integer(param$Age2), " ", param$VulnProb23, " ", as.integer(param$Age4), " ", param$VulnProb33, " yes")
      )
    }
    
    # Combine all parts of the file content
    file_content <- c(
      'LandisData "BDA Agent"',
      sprintf("BDAAgentName budworm%d", i),
      "BDPCalibrator 1",
      "SRDMode mean",
      sprintf("StartYear %d", period[1]),
      sprintf("EndYear %d", period[2]),
      ">> Regional Outbreak Inputs",
      paste("OutbreakPattern", param$OutbreakPattern),
      paste("MaxInterval", param$MaxInterval),
      paste("MinInterval", param$MinInterval),
      paste("TimeSinceLastEpidemic", param$TimeSinceLastEpidemic),
      paste("TemporalType", param$TemporalType),
      paste("MinROS", param$MinROS),
      paste("MaxROS", param$MaxROS),
      ">> Dispersal Parameters",
      paste("Dispersal", param$Dispersal),
      paste("DispersalRate", param$DispersalRate),
      paste("EpidemicThresh", param$EpidemicThresh),
      paste("InitialEpicenterNum", param$InitialEpicenterNum),
      paste("OutbreakEpicenterCoeff", param$OutbreakEpicenterCoeff),
      paste("OutbreakEpicenterThresh", param$OutbreakEpicenterThresh),
      paste("SeedEpicenter", param$SeedEpicenter),
      paste("SeedEpicenterMax", param$SeedEpicenterMax),
      paste("SeedEpicenterCoeff", param$SeedEpicenterCoeff),
      paste("DispersalTemplate", param$DispersalTemplate),
      ">> Neighborhood Resource Inputs",
      paste("NeighborFlag", param$NeighborFlag),
      paste("NeighborSpeedUp", param$NeighborSpeedUp),
      paste("NeighborRadius", param$NeighborRadius),
      paste("NeighborShape", param$NeighborShape),
      paste("NeighborWeight", param$NeighborWeight),
      ">> Intensity Class Thresholds",
      paste("IntensityClass2_BDP", param$IntensityClass2_BDP),
      paste("IntensityClass3_BDP", param$IntensityClass3_BDP),
      "BDASpeciesParameters",
      ">> SpeciesName     Age1 SRDProb1  Age2 SRDProb2  Age3 SRDProb3  Age1 VulnProb1  Age2 VulnProb2  Age3 VulnProb3  fuel",
      species_data
    )
    
    writeLines(file_content, fileConn)
    close(fileConn)
    
    # Keep track of active BDA files (only once per i if needed)
    if (active && j == 1) {
      bda_input_files <- c(bda_input_files, sprintf("\tBase-BDA_budworm%d.txt", i))
    }
    ### consider replicate:
    if (file.exists(filePath)) {
      for (rep_i in seq_len(expDesign$rep)) {
        replicateName <- sprintf("Base-BDA_budworm%d-%d-rep%d.txt", i, j, rep_i)
        replicatePath <- file.path(wwd, replicateName)
        file.copy(filePath, replicatePath, overwrite = TRUE)
        
        # add them to bda_input_files
        if (active) {
          bda_input_files <- c(bda_input_files, paste0("\t", replicateName))
        }
      }
    }
  }
}


periodFiles <- character(0)

#### create base-bda.txt file based on created Base-BDA_budworm per period , simID and rep
for (i in seq_along(periods)) {
  period <- periods[[i]]
  for (j in seq_len(nrow(base_params))) {
    filenameParam <- sprintf("Base-BDA_budworm%d-%d.txt", i, j)
  }
  filePeriod <- file.path(wwd, sprintf("Base-BDA_budworm%d.txt", i))
  sourceParamFile <- file.path(wwd, sprintf("Base-BDA_budworm%d-%d.txt", i, 1))
  file.copy(sourceParamFile, filePeriod, overwrite = TRUE)
  periodFiles <- c(periodFiles, sprintf("Base-BDA_budworm%d.txt", i))
}

bdaFilePath <- file.path(wwd, "base-bda.txt")
fileConn <- file(bdaFilePath, open = "w")
# First lines:
bdaContent <- c(
  'LandisData "Base BDA"',
  "",
  "Timestep  1",
  "",
  "MapNames\t\tbda/{agentName}-severity-{timestep}.tif\t\t",
  "SRDMapNames\t\tbda/{agentName}-SRD-{timestep}.tif\t",
  "NRDMapNames\t\tbda/{agentName}-NRD-{timestep}.tif\t",
  "BDPMapNames \tbda/{agentName}-BDP-{timestep}.tif\t",
  "LogFile\t\tbda/bda-log.csv",
  ""
)

# Now you add the "BDAInputFiles" lines:
if (length(periodFiles) > 0) {
  # The first line is "BDAInputFiles" plus the first file
  bda_input_lines <- paste("BDAInputFiles", periodFiles[1])
  
  # Then subsequent lines are tab-indented
  if (length(periodFiles) > 1) {
    for (f in periodFiles[-1]) {
      bda_input_lines <- c(
        bda_input_lines,
        paste0("\t", f)
      )
    }
  }
  # Combine
  bdaContent <- c(bdaContent, bda_input_lines)
}

writeLines(bdaContent, fileConn)
close(fileConn)




###  Expand simulation info to include BDA combinations

simInfo <- expand.grid(
  areaName = names(expDesign$mgmt),
  scenario = expDesign$scenario,
  mgmt = unlist(expDesign$mgmt),
  cropped = unlist(expDesign$cropped),
  spinup = expDesign$spinup,
  includeSnags = includeSnags,
  fire = expDesign$fire,
  BDA = expDesign$BDA,
  wind = expDesign$wind,
  harvest = expDesign$harvest,
  replicate = seq_len(expDesign$rep),
  BDAfile=seq(1, nrow(base_params))
)  %>% arrange(replicate)

simInfo

# Assign simulation IDs

sID <- seq(1:nrow((simInfo)))
simInfo <- data.frame(simID = str_pad(sID, nchar(max(sID)), pad = "0"), simInfo)
simInfo$harvest <- simInfo$mgmt != "noHarvest"
row.names(simInfo) <- 1:nrow(simInfo)
simInfo
write.csv(simInfo, file = "simInfo.csv", row.names = F)
### Parallel processing to generate simulation folders and files
n <- floor(detectCores() * 0.5)
cl <- makeCluster(n, outfile = "")
registerDoSNOW(cl)

foreach(i = 1:nrow(simInfo)) %dopar% {
  require(raster)
  require(dplyr)
  
  # Extract simulation parameters
  simID <- as.character(simInfo[i, "simID"])
  areaName <- as.character(simInfo[i, "areaName"])
  scenario <- as.character(simInfo[i, "scenario"])
  mgmt <- as.character(simInfo[i, "mgmt"])
  spinup <- simInfo[i, "spinup"]
  cropped <- simInfo[i, "cropped"]
  harvest <- simInfo[i, "harvest"]
  fire <- simInfo[i, "fire"]
  wind <- simInfo[i, "wind"]
  BDA <- simInfo[i, "BDA"]
  includeSnags <- simInfo[i, "includeSnags"]
  #bdaFile <- as.character(simInfo[i, "bdaFile"])  # Dynamically generated BDA file
  BDAfileIndex <- simInfo[i, "BDAfile"]    # j
  repIndex     <- simInfo[i, "replicate"]  # replicate
  
  # Create directory for simulation
  dir.create(simID, showWarnings = FALSE)
  
  ###############################################
  ### Initial rasters and attribute files
  if (cropped) {
    # Read cropped raster
    rCrop <- raster(paste0(inputDir, "/studyArea_", areaName, "_cropped.tif"))
    
    # Initial communities raster
    r <- raster(paste0(inputDir, "/initial-communities_", areaName, ".tif"))
    r <- crop(r, rCrop)
    r[is.na(rCrop)] <- NA
    writeRaster(r, file = paste0(simID, "/initial-communities.tif"),
                datatype = 'INT4S', overwrite = TRUE, NAflag = 0)
    
    # Landtypes raster
    r <- raster(paste0(inputDir, "/landtypes_", areaName, ".tif"))
    r <- crop(r, rCrop)
    r[is.na(rCrop)] <- NA
    writeRaster(r, file = paste0(simID, "/landtypes.tif"),
                datatype = 'INT4S', overwrite = TRUE, NAflag = 0)
  } else {
    # Copy initial communities and landtypes files without cropping
    file.copy(paste0(inputDir, "/initial-communities_", areaName, ".tif"),
              paste0(simID, "/initial-communities.tif"), overwrite = TRUE)
    file.copy(paste0(inputDir, "/landtypes_", areaName, ".tif"),
              paste0(simID, "/landtypes.tif"), overwrite = TRUE)
  }
  
  ###############################################
  ### Handle Snag Files
  if (includeSnags) {
    fName <- paste0(inputDir, "/initial-snags_", areaName, ".txt")
    file.copy(fName, paste0(simID, "/initial-snags.txt"), overwrite = TRUE)
  }
  
  ###############################################
  ### Smooth Age Classes (optional)
  if (smoothAgeClasses) {
    spp <- read.table(paste0(inputDir, "/species_", areaName, ".txt"),
                      skip = 1, header = FALSE, comment.char = ">")[, 1]
    initComm <- paste0(inputDir, "/initial-communities_", areaName, ".txt")
    x <- readLines(initComm)
    tmp <- strsplit(x, " ")
    tmp <- lapply(tmp, function(x) x[which(nchar(x) > 0)])
    
    for (j in seq_along(tmp)) {
      fname <- paste0(simID, "/initial-communities.txt")
      if (j == 1) file.create(fname)
      l <- tmp[[j]]
      if (l[1] %in% spp) {
        sp <- l[1]
        cohortAge <- as.numeric(l[-1])
        cohortAge <- round(cohortAge + runif(length(cohortAge), min = -9, max = 0))
        l <- paste0(sp, "\t", paste(cohortAge, collapse = " "))
      } else {
        l <- paste0(l, collapse = " ")
      }
      write(l, file = fname, append = TRUE)
    }
  } else {
    file.copy(paste0(inputDir, "/initial-communities_", areaName, ".txt"),
              paste0(simID, "/initial-communities.txt"), overwrite = TRUE)
  }
  file.copy(paste0(inputDir, "/landtypes_",
                   areaName, ".txt"),
            paste0(simID, "/landtypes.txt"),
            overwrite = T)
  
  ###############################################
  ### Succession and Climate Files
  if (spinup) {
    file.copy(paste0(inputDir, "/forCS-input_", areaName, "_spinup.txt"),
              paste0(simID, "/forCS-input.txt"), overwrite = TRUE)
  } else {
    file.copy(paste0(inputDir, "/forCS-input_", areaName, "_", scenario, ".txt"),
              paste0(simID, "/forCS-input.txt"), overwrite = TRUE)
  }
  
  file.copy(paste0(inputDir, "/forCS-climate_", areaName, "_", scenario, ".txt"),
            paste0(simID, "/forCS-climate.txt"), overwrite = TRUE)
  
  if (as.numeric(forCSVersion) >= 3.1) {
    file.copy(paste0(inputDir, "/ForCS_DM_", areaName, "_", scenario, ".txt"),
              paste0(simID, "/ForCS_DM.txt"), overwrite = TRUE)
  }
  
  ###############################################
  ### Disturbances
  if (!spinup) {
    # Harvest Files
    if (harvest) {
      if (cropped) {
        # Stand map raster
        r <- raster(paste0(inputDir, "/stand-map_", areaName, ".tif"))
        r <- crop(r, rCrop)
        r[is.na(rCrop)] <- NA
        writeRaster(r, file = paste0(simID, "/stand-map.tif"),
                    datatype = 'INT4S', overwrite = TRUE, NAflag = 0)
        
        # Management areas raster
        r <- raster(paste0(inputDir, "/mgmt-areas_", areaName, ".tif"))
        r <- crop(r, rCrop)
        r[is.na(rCrop)] <- NA
        writeRaster(r, file = paste0(simID, "/mgmt-areas.tif"),
                    datatype = 'INT4S', overwrite = TRUE, NAflag = 0)
      } else {
        file.copy(paste0(inputDir, "/stand-map_", areaName, ".tif"),
                  paste0(simID, "/stand-map.tif"), overwrite = TRUE)
        file.copy(paste0(inputDir, "/mgmt-areas_", areaName, ".tif"),
                  paste0(simID, "/mgmt-areas.tif"), overwrite = TRUE)
      }
      file.copy(paste0(inputDir, "/biomass-harvest_", areaName, "_", mgmt, ".txt"),
                paste0(simID, "/biomass-harvest.txt"), overwrite = TRUE)
    }
    
    # Wind Files
    if (wind) {
      file.copy(paste0(inputDir, "/base-wind_", areaName, ".txt"),
                paste0(simID, "/base-wind.txt"), overwrite = TRUE)
    }
    
    # Fire Files
    if (fire) {
      file.copy(paste0(inputDir, "/base-fire_", areaName, "_", scenario, ".txt"),
                paste0(simID, "/base-fire.txt"), overwrite = TRUE)
      fireYears <- ifelse(scenario == "baseline", 0, c(10, 40, 70))
     # for (y in fireYears) {
      #  r <- raster(paste0(inputDir, "/fire-regions_", areaName, "_", scenario, "_", y, ".tif"))
       # r <- crop(r, rCrop)
        #r[is.na(rCrop)] <- NA
        #writeRaster(r, file = paste0(simID, "/fire-regions_", y, ".tif"),
        #            datatype = 'INT4S', overwrite = TRUE, NAflag = 0)
      #}
      for (y in fireYears) {
        r <- raster(paste0(inputDir, "/fire-regions_", areaName, "_", scenario, "_", y, ".tif"))
        if (cropped) {
          r <- crop(r, rCrop)
          r[is.na(rCrop)] <- NA
        }
        writeRaster(r, file = paste0(simID, "/fire-regions_", y, ".tif"),
                    datatype = 'INT4S', overwrite = TRUE, NAflag = 0)
      }
    }
    if (BDA) {
      # General BDA-related file (if exists and applicable to all simulations)
      general_bda_file <- paste0(wwd, "/base-bda.txt")
      if (file.exists(general_bda_file)) {
        file.copy(general_bda_file, paste0(simID, "/base-bda.txt"), overwrite = TRUE)
      }
      
      # Copy each base-BDA_budworm file using a simple for loop
      for (p in seq_along(periods)) {
        sourceFile  <- file.path(
          wwd,
          sprintf("Base-BDA_budworm%d-%d-rep%d.txt", p, BDAfileIndex, repIndex)
        )
        targetFile  <- file.path(simID, sprintf("Base-BDA_budworm%d.txt", p))
        if (file.exists(sourceFile)) {
          file.copy(sourceFile, targetFile, overwrite = TRUE)
        }
      }
    }
    
    
    # Adjust Scenario File
    scenFile <- paste0(inputDir, "/scenario.txt")
    x <- readLines(scenFile)
    if (!fire) x[grep("Base Fire", x)] <- paste(">>", x[grep("Base Fire", x)])
    if (!harvest) x[grep("Biomass Harvest", x)] <- paste(">>", x[grep("Biomass Harvest", x)])
    if (!wind) x[grep("Base Wind", x)] <- paste(">>", x[grep("Base Wind", x)])
    if (!BDA) x[grep("Base BDA", x)] <- paste(">>", x[grep("Base BDA", x)])
    x[grep("Duration", x)] <- paste("Duration", simDuration)
    writeLines(x, con = paste0(simID, "/scenario.txt"))
  } else {
    file.copy(paste0(inputDir, "/scenario_spinup.txt"),
              paste0(simID, "/scenario.txt"), overwrite = TRUE)
  }
  
  ###############################################
  ### Copy Species File
  file.copy(paste0(inputDir, "/species_", areaName, ".txt"),
            paste0(simID, "/species.txt"), overwrite = TRUE)
  
  ###############################################
  ### Save Simulation Information
  write.table(t(simInfo[i, ]), file = paste0(simID, "/README.txt"),
              quote = FALSE, col.names = FALSE)
}

stopCluster(cl)


