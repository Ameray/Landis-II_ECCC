## organisation ECCC
### Ameray Abderrahmane

#### this script is adapted for Local sensitivity analysis only
#### step 1 prepare inputs for local sensitivity analysis. 
#### change one para by one in BDA model input files

# Clear the environment

rm(list = ls())

# Define base directories
base_dir <- "D:/Wena_Project/Testing/test/SA"   
setwd(base_dir)

# Create a working directory for today's date
wwd <- paste(getwd(), Sys.Date(), sep = "/")
dir.create(wwd)
setwd(wwd)

# Load required libraries
require(stringr)
require(dplyr)
require(raster)
require(parallel)
require(doSNOW)
require(lhs)

# Input path (LANDIS)
inputDir <- "D:/Wena_Project/Testing/frqnt_2022-25-main/inputsLandis"

# Define simulation parameters
forCSVersion <- "3.1"

simDuration <- 32    ### run simulation during one cycle of 32 years
smoothAgeClasses <- TRUE
includeSnags <- FALSE

expDesign <- list(
  area = c("mixedwood-042-51"),
  scenario = c("baseline"),
  mgmt = list("mixedwood-042-51" = c("generic")), #"noHarvest"
  spinup = FALSE,
  cropped = list("mixedwood-042-51" = TRUE),
  fire = TRUE,
  BDA = TRUE,
  wind = TRUE,
  harvest = TRUE,
  rep = 1    ### since we fixed ndomNumberSeed  4357
)

# Reference combination values
reference_values <- list(
  OutbreakPattern = "CyclicUniform",
  MaxInterval = 32,
  MinInterval = 32,
  TimeSinceLastEpidemic = 20,
  TemporalType = "pulse",
  MinROS = 1,    # change 0 to 1 to make it continous
  MaxROS = 3,
  Dispersal = "no",
  DispersalRate = 50000,
  EpidemicThresh = 0.1,
  InitialEpicenterNum = 5,
  OutbreakEpicenterCoeff = 0.01,
  OutbreakEpicenterThresh = 0,
  SeedEpicenter = "yes",
  SeedEpicenterMax = 10,
  SeedEpicenterCoeff = 1,
  DispersalTemplate = "MaxRadius",
  NeighborFlag = "yes",
  NeighborSpeedUp = "none",
  NeighborRadius = 1000,
  NeighborShape = "uniform",
  NeighborWeight = 0.5,
  IntensityClass2_BDP = 0.33,
  IntensityClass3_BDP = 0.67,
  Age1 = 0,
  Age2 = 20,
  Age3 = 50,
  SRDProb1 = 0.1, 
  SRDProb2 = 0.3,
  SRDProb3 = 0.7,
  VulnProb1 = 0.1,
  VulnProb2 = 0.3,
  VulnProb3 = 0.7
)

# Parameter Ranges for Local Sensitivity Analysis ,test different value than the reference one
parameter_ranges <- list(
  TemporalType = setdiff(c("pulse", "variablepulse"), reference_values$TemporalType),
  MinROS = setdiff(c(1,2, 3), reference_values$MinROS),
  MaxROS = setdiff(c(1, 2,3), reference_values$MaxROS),
  Dispersal = setdiff(c("no"), reference_values$Dispersal),
  DispersalRate = setdiff(c(40000, 60000, 70000, 80000, 90000), reference_values$DispersalRate),
  EpidemicThresh = setdiff(c(0.25, 0.5, 0.75, 1.0), reference_values$EpidemicThresh),
  InitialEpicenterNum = setdiff(c(0,5, 10, 20, 30, 40), reference_values$InitialEpicenterNum),
  OutbreakEpicenterCoeff = setdiff(c(0.2, 0.4, 0.6, 0.8, 1), reference_values$OutbreakEpicenterCoeff),
  OutbreakEpicenterThresh = setdiff(c(1, 2, 3), reference_values$OutbreakEpicenterThresh),
  SeedEpicenter = setdiff(c("no"), reference_values$SeedEpicenter),
  SeedEpicenterMax = setdiff(c(0, 20, 30), reference_values$SeedEpicenterMax),
  SeedEpicenterCoeff = setdiff(c(0.2, 0.4, 0.6, 0.8), reference_values$SeedEpicenterCoeff),
  DispersalTemplate = setdiff(c("4N", "8N", "12N", "24N"), reference_values$DispersalTemplate),
  NeighborFlag = setdiff(c("no"), reference_values$NeighborFlag),
  NeighborRadius = setdiff(c(500, 1500, 2000), reference_values$NeighborRadius),
  NeighborShape = setdiff(c("gaussian", "linear"), reference_values$NeighborShape),
  NeighborWeight = setdiff(c(1, 10, 25, 50, 75, 100), reference_values$NeighborWeight),
  IntensityClass2_BDP = setdiff(c(0.2, 0.4, 0.5), reference_values$IntensityClass2_BDP),
  IntensityClass3_BDP = setdiff(c(0.5, 0.7, 0.9), reference_values$IntensityClass3_BDP),
  Age2 = setdiff(c(30, 40, 50), reference_values$Age2),
  Age3 = setdiff(c(50, 60, 70, 80), reference_values$Age3),
  SRDProb1 = setdiff(c(0.05, 0.15, 0.2), reference_values$SRDProb1),
  SRDProb2 = setdiff(c(0.2, 0.4, 0.5), reference_values$SRDProb2),
  SRDProb3 = setdiff(c(0.5, 0.75, 1.0), reference_values$SRDProb3),
  VulnProb1 = setdiff(c(0.05, 0.15, 0.2), reference_values$VulnProb1),
  VulnProb2 = setdiff(c(0.2, 0.4, 0.5), reference_values$VulnProb2),
  VulnProb3 = setdiff(c(0.5, 0.75, 1.0), reference_values$VulnProb3)
)

# Initialize list for storing combinations
local_samples <- list()
file_counter <- 0

# Number of replicates
replicates <- 1

# Add the reference combination replicates
for (rep in 1:replicates) {
  reference_entry <- reference_values
  reference_entry$fileName <- sprintf("local-BDA_budworm-%d.txt", file_counter)
  reference_entry$replicate <- rep
  local_samples <- append(local_samples, list(reference_entry))
  file_counter <- file_counter + 1
}

# Generate combinations by varying one parameter at a time with replicates
for (param in names(parameter_ranges)) {
  for (value in parameter_ranges[[param]]) {
    for (rep in 1:replicates) {  # Add replicates
      modified_values <- reference_values
      modified_values[[param]] <- value
      modified_values$fileName <- sprintf("local-BDA_budworm-%d.txt", file_counter)
      modified_values$replicate <- rep
      local_samples <- append(local_samples, list(modified_values))
      file_counter <- file_counter + 1
    }
  }
}

# Convert to a data frame
local_samples_df <- bind_rows(local_samples)

# Save results to file
write.csv(local_samples_df, "localSensitivityParameters.txt", row.names = FALSE)

# Write BDA Parameter Files including  the referece  for Local Sensitivity 
for (i in seq_len(nrow(local_samples_df))) {
  param <- local_samples_df[i, ]
  fileConn <- file(param$fileName)
  writeLines(c(
    'LandisData "BDA Agent"',
    "BDAAgentName budworm",
    "BDPCalibrator 1",
    "SRDMode mean",
    ">> Regional Outbreak Inputs",
    paste0("OutbreakPattern ", param$OutbreakPattern),
    paste0("MaxInterval ", as.integer(param$MaxInterval)),
    paste0("MinInterval ", as.integer(param$MinInterval)),
    paste0("TimeSinceLastEpidemic ", as.integer(param$TimeSinceLastEpidemic)),
    paste0("TemporalType ", param$TemporalType),
    paste0("MinROS ", as.integer(param$MinROS)),
    paste0("MaxROS ", as.integer(param$MaxROS)),
    ">> Dispersal Parameters",
    paste0("Dispersal ", param$Dispersal),
    paste0("DispersalRate ", as.integer(param$DispersalRate)),
    paste0("EpidemicThresh ", param$EpidemicThresh),
    paste0("InitialEpicenterNum ", as.integer(param$InitialEpicenterNum)),
    paste0("OutbreakEpicenterCoeff ", param$OutbreakEpicenterCoeff),
    paste0("OutbreakEpicenterThresh ", as.integer(param$OutbreakEpicenterThresh)),
    paste0("SeedEpicenter ", param$SeedEpicenter),
    paste0("SeedEpicenterMax ", as.integer(param$SeedEpicenterMax)),
    paste0("SeedEpicenterCoeff ", param$SeedEpicenterCoeff),
    paste0("DispersalTemplate ", param$DispersalTemplate),
    ">> Neighborhood Resource Inputs",
    paste0("NeighborFlag ", param$NeighborFlag),
    paste0("NeighborSpeedUp ", param$NeighborSpeedUp),
    paste0("NeighborRadius ", param$NeighborRadius),
    paste0("NeighborShape ", param$NeighborShape),
    paste0("NeighborWeight ", param$NeighborWeight),
    ">> Intensity Class Thresholds",
    paste0("IntensityClass2_BDP ", param$IntensityClass2_BDP),
    paste0("IntensityClass3_BDP ", param$IntensityClass3_BDP),
    "BDASpeciesParameters",
    ">>SpeciesName     Age1 SRDProb1  Age2 SRDProb2  Age3 SRDProb3  Age1 VulnProb1  Age2 VulnProb2  Age3 VulnProb3  fuel",
    paste0("ABIE.BAL ", as.integer(param$Age1), " ", param$SRDProb1, " ", as.integer(param$Age2), " ", param$SRDProb2, 
           " ", as.integer(param$Age3), " ", param$SRDProb3, " ", as.integer(param$Age1), " ", param$VulnProb1, 
           " ", as.integer(param$Age2), " ", param$VulnProb2, " ", as.integer(param$Age3), " ", param$VulnProb3, " yes"),
    paste0("PICE.GLA ", as.integer(param$Age1), " ", param$SRDProb1, " ", as.integer(param$Age2), " ", param$SRDProb2, 
           " ", as.integer(param$Age3), " ", param$SRDProb3, " ", as.integer(param$Age1), " ", param$VulnProb1, 
           " ", as.integer(param$Age2), " ", param$VulnProb2, " ", as.integer(param$Age3), " ", param$VulnProb3, " yes"),
    paste0("PICE.MAR ", as.integer(param$Age1), " ", param$SRDProb1, " ", as.integer(param$Age2), " ", param$SRDProb2, 
           " ", as.integer(param$Age3), " ", param$SRDProb3, " ", as.integer(param$Age1), " ", param$VulnProb1, 
           " ", as.integer(param$Age2), " ", param$VulnProb2, " ", as.integer(param$Age3), " ", param$VulnProb3, " yes")
  ), fileConn)
  close(fileConn)
}


#### step 2 prepare   all inputs of landis-ii model 
### Expand simulation info to include BDA combinations
##### add bdafileName to simInfo

simInfo <- list()
for (a in names(expDesign$mgmt)) {
  simInfo[[a]] <- expand.grid(
    areaName = a,
    scenario = expDesign$scenario,
    mgmt = expDesign$mgmt[[a]],
    cropped = expDesign$cropped[[a]],
    spinup = expDesign$spinup,
    includeSnags = includeSnags,
    fire = expDesign$fire,
    BDA = expDesign$BDA,
    wind = expDesign$wind,
    harvest = expDesign$harvest,
    rep = expDesign$rep, # Keep as 1
    bdaFile = local_samples_df$fileName
  )
}
simInfo <- do.call("rbind", simInfo) %>% arrange(replicate)

### Add replicate column for parameter variation
# Match the replicate from local_samples_df to simInfo based on bdaFile
replicate_mapping <- local_samples_df %>% dplyr::select(fileName, replicate)
simInfo <- simInfo %>%
  left_join(replicate_mapping, by = c("bdaFile" = "fileName"))

### Assign simulation IDs
sID <- ((1:nrow(simInfo)))
simInfo <- data.frame(simID = str_pad(sID, nchar(max(sID)), pad = "0"), simInfo)
simInfo$harvest <- simInfo$mgmt != "noHarvest"
row.names(simInfo) <- 1:nrow(simInfo)

### Save simulation info to CSV
write.csv(simInfo, file = "simInfo.csv", row.names = FALSE)


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
  bdaFile <- as.character(simInfo[i, "bdaFile"])  # Dynamically generated BDA file
  
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
      for (y in fireYears) {
        r <- raster(paste0(inputDir, "/fire-regions_", areaName, "_", scenario, "_", y, ".tif"))
        r <- crop(r, rCrop)
        r[is.na(rCrop)] <- NA
        writeRaster(r, file = paste0(simID, "/fire-regions_", y, ".tif"),
                    datatype = 'INT4S', overwrite = TRUE, NAflag = 0)
      }
    }
    if(BDA) {
      ### BDA - Budworm
      file.copy(paste0(inputDir, "/base-bda.txt"),
                paste0(simID),
                overwrite = T)
    }
    
    # BDA Files (Dynamic)
    if (BDA) {
      file.copy(bdaFile, paste0(simID, "/base-BDA_budworm.txt"), overwrite = TRUE)
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


###step 3: run landis model
# Clear workspace
######
rm(list = ls())

# Set working directory
wwd <- "D:/Wena_Project/Testing/Test/SA/2024-12-04"
setwd(wwd)
wwd <- getwd()
simInfo <- read.csv("simInfo.csv", colClasses = c(simID = "character"))
simDir <- simInfo$simID

# Load required libraries
require(parallel)
require(doSNOW)

# Define the number of simultaneous processes and batch size
maxProcesses <- 4  # Run up to 10 simulations per batch
numSimulations <- length(simDir)
numBatches <- ceiling(numSimulations / maxProcesses)

# Loop through batches
for (batch in 1:numBatches) {
  # Calculate start and end indices for the current batch
  startIndex <- (batch - 1) * maxProcesses + 1
  endIndex <- min(batch * maxProcesses, numSimulations)
  
  # Get the subset of simulations for this batch
  currentSimulations <- simDir[startIndex:endIndex]
  
  # Create a cluster for the current batch
  cl <- makeCluster(min(maxProcesses, length(currentSimulations)), outfile = "")
  registerDoSNOW(cl)
  
  # Run simulations in parallel for the current batch
  foreach(i = startIndex:endIndex, .packages = c("parallel")) %dopar% {
    # Set working directory for the simulation
    setwd(paste(wwd, simInfo[i, "simID"], sep = "/"))
    
    # Read the existing README file and update it
    readmeOld <- readLines("README.txt")
    readmeOld <- readmeOld[1:which(grepl(tail(colnames(simInfo), 1), readmeOld))]
    x <- as.character(shell("landis-ii-extensions.cmd list", intern = TRUE))
    
    # Append details to the README file
    sink(file = "README.txt")
    for (l in seq_along(readmeOld)) {
      cat(readmeOld[l])
      cat("\n")
    }
    sink()
    
    sink(file = "README.txt", append = TRUE)
    cat("\n")
    cat("#######################################################################\n")
    cat("########### Installed LANDIS-II extensions\n")
    cat("#######################################################################\n")
    for (l in seq_along(x)) {
      cat(x[l])
      cat("\n")
    }
    cat("\n")
    cat("#######################################################################\n")
    cat("########### System Info\n")
    cat(write.table(as.data.frame(Sys.info()),
                    quote = FALSE, col.names = FALSE))
    sink()
    
    # Run the simulation
    shell("landis-ii-7.cmd scenario.txt", wait = TRUE)
  }
  
  # Stop the cluster after the batch completes
  stopCluster(cl)
  
  # Print batch completion message
  cat(sprintf("Batch %d/%d completed.\n", batch, numBatches))
  
  # Perform garbage collection to free up memory
  gc()
}



#### step 4 read the outputs and calculate severity 

###################################################################
### ForCS BDA output processing
### IMPORTANT: For this script to work, disturbances must be
### simulated in the specified order.
###################################################################
require(data.table)
require(dplyr)
require(raster)
require(stringr)  # Removed parallel processing libraries as per request.

# Clear environment
rm(list = ls())

dir_path <- "D:/Wena_Project/Testing/Test/SA/2024-12-04"
setwd(dir_path)
#### fetchBDAParamFnc function controlling the bda parameters
source("D:/Wena_Project/Testing/Test/fetchBDAParamFnc.R")   
### fetching outputs
a <- "mixedwood-042-51"
d <- c("bda" = 4, "wind" = 3, "harvest" = 2, "fire" = 1)  # IMPORTANT: put in the order they were simulated in LANDIS
simDir <- "D:/Wena_Project/Testing/Test/SA/2024-12-04"
simName <- gsub("simPkg_mixedwood-042-51_", "", basename(simDir))
simInfo <- read.csv(paste(simDir, "simInfo.csv", sep = "/"), colClasses = c("simID" = "character"))
colnames(simInfo)


#  directory indexing
simIDs <- simInfo$simID
simIDs <- str_pad(simIDs, width = max(nchar(simIDs)), side = "left", pad = "0")
allDirs <- list.dirs(simDir, full.names = FALSE, recursive = FALSE)
dirIndex <- which(simIDs %in% allDirs & simInfo$areaName == a)


print(simIDs)    # Check the  simIDs
print(allDirs)   # Check the actual directory names
print(dirIndex)  # Ensure valid indices are returned

outputList <- list()  # Initialize the list for storing outputs

for(i in dirIndex) {
  output <- list()
  ### sim variables
  simID <- simIDs[i]
  sDir <- paste(simDir, simID, sep ="/")
  areaName <- simInfo[i, "areaName"]
  scenario <- simInfo[i, "scenario"]
  mgmtScenario <- simInfo[i, "mgmt"]
  mgmtScenarioName <- mgmtScenario
  harvest <- simInfo[i, "harvest"]
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
  
  ### fetching BDA susceptibility age classes
  
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
  
  bdaParam <- fetchBDAParam1(paste(sDir, "base-bda_budworm.txt", sep = "/"))
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
                     replicate =  replicate,
                     df2)
  
  output[["fluxesBDA"]] <- df2
  outputList[[i]] <- output
}

### Final aggregation
fluxesBDA <- do.call("rbind", lapply(outputList, function(item) item[["fluxesBDA"]]))
save(fluxesBDA, file = paste0("output_fluxesBDA_", simName, ".RData"))

fluxesBDA


#### read the final combined data and calculate severity for each simulation
# Set the directory path
dir_path <- "D:/Wena_Project/Testing/Test/SA/2024-12-04"
setwd(dir_path)

# Load necessary libraries
library(tidyverse)
library(dplyr)
# Load the dataset
fluxesBDA <- get(load("output_fluxesBDA_2024-12-04.RData"))

# Inspect column names
colnames(fluxesBDA)

# Calculate severity at the cell level
df_cell <- fluxesBDA %>%
  group_by(simID,replicate, scenario, mgmtScenario, species, Time, row, column) %>%
  summarize(
    AGB_woody_host = sum(AGB_woody),
    AGBtoDOM_woody = sum(AGBtoDOM_woody),
    severity = (AGBtoDOM_woody / AGB_woody_host)*100,
    .groups = "drop"  # Drop grouping after summarization
  )
df_cell
# Save the cell-level summary to a CSV file
write.csv(df_cell, "severity.txt", row.names = FALSE)

#########Step 5

########### Local sensitivity analysis by calculting the variation to reference values simID==1
library(dplyr)
library(ggplot2)
library(tidyr)

# Load the severity and tested parameters files
severity <- read.csv("H:/SA/2024-12-04/severity.txt")
severity
tested_parameters <- read.csv("H:/SA/2024-12-04/localSensitivityParameters.txt")
tested_parameters <- tested_parameters %>%
  mutate(simID = as.integer(sub(".*-(\\d+)\\.txt$", "\\1", fileName)) + 1,
         simID = ceiling(simID / 1)) ### divide by replicate number
head(tested_parameters)
# Remove unnecessary columns from `tested_parameters`
tested_parameters <- tested_parameters %>%
  dplyr::select(-c(OutbreakPattern, MaxInterval, MinInterval, TimeSinceLastEpidemic, fileName))
colnames(tested_parameters)
summary(tested_parameters)
# Extract the simulation number (simID) and group every 3 simID into one replicate if there is replicate if not divide by 1
severity <- severity %>%
  mutate(simID = ceiling(simID / 1))
summary(severity)
severity
# Calculate mean severity grouped by simID, replicate, and Time
severity_mean <- severity %>%
  group_by(simID, replicate) %>%
  summarise(severity = mean(severity, na.rm = TRUE), .groups = "drop")
head(severity_mean)
summary(severity_mean)
# Merge severity_mean with tested_parameters
df_merged <- severity_mean %>%
  inner_join(tested_parameters, by = c("simID", "replicate"))

write.csv(df_merged, "df_merged.csv", row.names = FALSE)

# Calculate the mean severity across replicates for each simID
df_mean <- df_merged %>%
  group_by(simID,TemporalType, MinROS, MaxROS, Dispersal, DispersalRate, EpidemicThresh, InitialEpicenterNum,
           OutbreakEpicenterCoeff, OutbreakEpicenterThresh, SeedEpicenter, SeedEpicenterMax, 
           SeedEpicenterCoeff, DispersalTemplate, NeighborFlag, NeighborSpeedUp, NeighborRadius, 
           NeighborShape, NeighborWeight, IntensityClass2_BDP, IntensityClass3_BDP, 
           Age1, Age2, Age3, SRDProb1, SRDProb2, SRDProb3, VulnProb1, VulnProb2, VulnProb3) %>%
  summarise(severity = mean(severity, na.rm = TRUE), .groups = "drop")  # Calculate mean severity
df_mean
colnames(df_mean)
# Export the mean severity data
write.csv(df_mean, "df_mean.csv", row.names = FALSE)

# Filter the reference severity for simID = 1
reference_severity <- df_mean %>%
  filter(simID == 1) %>%
  dplyr::select(-simID)

# Calculate the variation in severity relative to the reference
variation <- df_mean %>%
  filter(simID != 1) %>%
  rowwise() %>%
  mutate(variation = 100 * (severity - reference_severity$severity[1]) / reference_severity$severity[1]) %>%
  ungroup()

# Include simID and variation in the output
variation <- variation %>%
  dplyr::select(simID, variation, everything())
variation
# Reshape variation data for plotting
variation_long <- variation %>%
  mutate(across(everything(), as.character)) %>%  # Convert all columns to character
  pivot_longer(cols = -c(simID, variation), names_to = "Parameter", values_to = "Value")

# Reshape reference_severity data for comparison
reference_long <- reference_severity %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "ReferenceValue")

# Identify varied parameters
varied_parameters <- variation_long %>%
  left_join(reference_long, by = "Parameter") %>%
  filter(Value != ReferenceValue)

# Define parameter groups similar to the BDA model
varied_parameters <- varied_parameters %>%
  mutate(
    ParameterGroup = case_when(
      Parameter %in% c(
        "OutbreakPattern", "MaxInterval", "MinInterval", 
        "TimeSinceLastEpidemic", "TemporalType", "MinROS", "MaxROS"
      ) ~ "Regional Outbreak Inputs",
      Parameter %in% c(
        "Dispersal", "DispersalRate", "EpidemicThresh", 
        "InitialEpicenterNum", "OutbreakEpicenterCoeff", 
        "OutbreakEpicenterThresh", "SeedEpicenter", 
        "SeedEpicenterMax", "SeedEpicenterCoeff", 
        "DispersalTemplate"
      ) ~ "Dispersal Parameters",
      Parameter %in% c(
        "NeighborFlag", "NeighborSpeedUp", "NeighborRadius", 
        "NeighborShape", "NeighborWeight"
      ) ~ "Neighborhood Resource Inputs",
      Parameter %in% c(
        "IntensityClass2_BDP", "IntensityClass3_BDP"
      ) ~ "Intensity Class Thresholds",
      Parameter %in% c(
        "Age1", "SRDProb1", "Age2", "SRDProb2", "Age3", 
        "SRDProb3", "VulnProb1", "VulnProb2", "VulnProb3"
      ) ~ "BDASpeciesParameters",
      TRUE ~ "severity"
    )
  )

varied_parameters<-varied_parameters%>% filter(Parameter!="severity")
varied_parameters
a<-varied_parameters %>% group_by(ParameterGroup) %>% summarise(n = n())
a
# Plotting function for each parameter group
plot_group <- function(group_name, data) {
  ggplot(data, aes(x = paste(Parameter, Value, sep = " = "), y = as.numeric(variation), fill = Parameter)) +
    geom_bar(stat = "identity", position = "dodge") +  # Bar plot with dodge position
    labs(
      title = paste("Variation Across", group_name),
      x = "Parameter and Value",
      y = "Severity Variation compared to the reference  (%)",
      fill = "Parameters"
    ) +
    theme_minimal() +
    scale_fill_viridis_d(option = "cividis") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "right"
    )+
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "top"
    )
}

# Filter data for each group and create plots
grouped_plots <- varied_parameters %>%
  split(.$ParameterGroup) %>%
  lapply(function(data) plot_group(unique(data$ParameterGroup), data))

# Display plots for each group
print(grouped_plots[["Regional Outbreak Inputs"]])
print(grouped_plots[["Dispersal Parameters"]])
print(grouped_plots[["Neighborhood Resource Inputs"]])
print(grouped_plots[["Intensity Class Thresholds"]])
print(grouped_plots[["BDASpeciesParameters"]])


# Function to save plots
output_directory <- "D:/NewVersion/Local_sensitivity_analysis/plots_LSA_severity/"
save_plots <- function(plots) {
  for (group_name in names(plots)) {
    plot_file <- paste0(output_directory, group_name, ".png") # Define the file path and name
    ggsave(plot_file, plot = plots[[group_name]], width = 10, height = 8, dpi = 300)
  }
}
# Call the function to save all plots
save_plots(grouped_plots)


##### plot all group in th same plot
# Create heat map data
heatmap_data1 <- varied_parameters %>%
  dplyr::select(-c(simID)) %>%  # Exclude simID from heatmap analysis
  mutate(variation = as.numeric(variation)) %>%
  group_by(ParameterGroup, Parameter) %>%
  summarise(
    MeanVariation = mean(variation, na.rm = TRUE),
    .groups = "drop"
  )

# Create a heat map
heatmap_plot <- ggplot(heatmap_data1, aes(x = Parameter, y = ParameterGroup, fill = MeanVariation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-10, 10), name = "Variation (%)"
  ) +
  labs(
    title = "Heat Map of Parameter Sensitivity",
    x = "Parameter",
    y = "Parameter Group"
  ) +
  theme_minimal() +
  theme_minimal() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
         axis.text.y = element_text(size = 10),
         panel.background = element_rect(fill = "white", color = NA),
         plot.background = element_rect(fill = "white", color = NA),
         legend.position = "top"
  )

# Print the heat map
print(heatmap_plot)

# Save the heat map to a file
ggsave("D:/NewVersion/Local_sensitivity_analysis/plots_LSA_severity/All_parameter_sensitivity.png", plot = heatmap_plot, width = 12, height = 8)

