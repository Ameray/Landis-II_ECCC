################################################################################
################################################################################
## organisation ECCC
### Preparation of Forcs simulation file packages  to calibrate biomass harvest
### Dominic Cyr /Ameray Abderrahmane
################################################################################

#############
# Clear the environment 
rm(list = ls())

# Define the base directory

dir_path <- "D:/JaneProject"
setwd(dir_path)
# Create a working directory for today's date
wwd <- paste(getwd(), Sys.Date(), sep = "/")
dir.create(wwd)
setwd(wwd)
### Input path (LANDIS input in your laptop)
inputDir <- "D:/Wena_Project/Testing/frqnt_2022-25-main/inputsLandis"
#inputPathScripts <- "../scripts"

### Load required libraries
require(stringr)
require(dplyr)
require(raster)
require(parallel)
require(doSNOW)

simDuration <- 20  #overridden if spinup == T
t0 <- 2020 
forCSVersion <- "3.1"
smoothAgeClasses <- T
includeSnags <- F ## require snag tables, also will occupy some cohorts (no living trees allowed)
expDesign <- list(area = c("temperate-2a-3b"),#temperate-2a-3b", "boreal-5a", "mixedwood-042-51"),#", "Hereford"
                  scenario = c("baseline"),#c("baseline", "RCP45", "RCP85")
                  mgmt = list(#Hereford = "1"),#c("1", "2", "3", "4", "noHarvest")),
                    # "mixedwood-042-51" =  c("generic", "noHarvest")),
                    "temperate-2a-3b" =  c("generic")), #noHarvest
                  spinup = F,
                  cropped  = list("temperate-2a-3b" = F),
                  rep = 1,
                  #ND natural disturbances 
                  ND = data.frame(
                    wind = c( TRUE), 
                    BDA = c(TRUE),
                    fire = c(TRUE),
                    ND_scenario=c("Wind_Sbw_Fire")
                   )
                  )

### BDA params

bdaParams <- list(ongoingAtT0 = F, ### HAS NOT BEEN TESTED WHEN TRUE ###
                  cycleMean = 35, ###NRCan entre 30 et 40 ans
                  cycleSD = 2.5,
                  durMortMean = 7,
                  durMortSD = 3)

#— PREPARE three biomass‐harvest files with adjusted rates
base_biomass  <- file.path(dir_path, "biomass-harvest.txt")
harvestrate   <- c(0.75,1,1.25,1.5)  ### 1.25%  to 1.3% is good rate for temperate study area 
# these are the % contributions in your original HarvestImplementations table
base_contrib  <- c(
  CPRS                 = 31.9,
  SeedTree             =  0.4,
  CommercialThinning   =  2.9,
  PC1                  = 10.5,   ## Regular Shelterwood
  PC2                  = 22.4,    ## Irregular Shelterwood (Slow Regeneration)
  PC3                  = 31.9     ## Shelterwood with Permanent Canopy
) / 100   # convert to fraction

# read once
template_lines <- readLines(base_biomass)

for (j in seq_along(harvestrate)) {
  hr <- harvestrate[j]
  out <- template_lines
  
  # loop by treatment name
  for (tr in names(base_contrib)) {
    # find the exact line starting with "1 <treatment>"
    idx <- grep(paste0("^1[[:space:]]+", tr, "\\b"), out)
    if (length(idx) != 1) {
      stop("Did not find exactly one line for treatment: ", tr)
    }
    # compute new rate in percent
    new_pct <- base_contrib[tr] * hr
    # rebuild that one line
    out[idx] <- sprintf("1 %s %.4f%%", tr, new_pct)
  }
  
  # write out
 writeLines(out, file.path(dir_path, paste0("biomass-harvest-", j, ".txt")))
}


simInfo <- list()
for (a in names(expDesign$mgmt)) {
  for (i in seq_along(expDesign$ND$wind)) {
    simInfo[[paste(a, i, sep = "_")]] <- expand.grid(
      areaName = a,
      scenario = expDesign$scenario,
      mgmt = expDesign$mgmt[[a]],
      cropped = expDesign$cropped[[a]],
      spinup = expDesign$spinup,
      includeSnags = includeSnags,
      wind = expDesign$ND$wind[i],
      BDA = expDesign$ND$BDA[i],
      fire = expDesign$ND$fire[i],
      ND_scenario=expDesign$ND$ND_scenario[i],
      #harvest = expDesign$disturbances$harvest[i],
      #dist = i - 1, #index of disturbances scenarios
      replicate = seq_len(expDesign$rep),
      harvestrate=harvestrate
    )
  }
}

simInfo <- do.call("rbind", simInfo) %>%
  arrange(replicate)



sID <- ((1:nrow(simInfo))-1)
simInfo <- data.frame(simID = str_pad(sID, nchar(max(sID)),
                                      pad = "0"),
                      simInfo)
## add harvest column
simInfo[,"harvest"] <- simInfo[,"mgmt"]!= "noHarvest"
write.csv(simInfo, file = "simInfo.csv", row.names = F)
simInfo
max_simID <- as.integer(max(simInfo$simID))
max_simID
### 
row.names(simInfo) <- 1:nrow(simInfo)
n <- floor(detectCores() * 0.5)
# #######
cl = makeCluster(n, outfile = "") ## 
registerDoSNOW(cl)

foreach(i = 1:nrow(simInfo)) %dopar% { 
    require(raster)
    simID <- as.character(simInfo[i,"simID"])
    areaName <- as.character(simInfo[i,"areaName"])
    scenario <- as.character(simInfo[i,"scenario"])
    mgmt <- as.character(simInfo[i,"mgmt"])
    spinup <- simInfo[i,"spinup"]
    replicate <- as.character(simInfo[i,"replicate"])
    cropped <- simInfo[i,"cropped"]
    harvest <- simInfo[i,"harvest"]
    fire <- simInfo[i,"fire"]
    wind <- simInfo[i,"wind"]
    BDA <- simInfo[i,"BDA"]
    includeSnags <- simInfo[i,"includeSnags"]
    dir.create(simID)
    
    ###############################################
    ### initial rasters and attribute files
    

    if(cropped) {
      rCrop <- raster(paste0(inputDir, "/studyArea_",
                             areaName, "_cropped.tif"))
      # initial communities
      r <- raster(paste0(inputDir, "/initial-communities_", areaName, ".tif"))
      # r <- crop(r, e)
      r <- crop(r, rCrop)
      r[is.na(rCrop)] <- NA
      writeRaster(r, file = paste0(simID, "/initial-communities.tif"),
                   datatype='INT4S', overwrite=TRUE, NAflag = 0)
      # landtypes
      r <- raster(paste0(inputDir, "/landtypes_", areaName, ".tif"))
      # r <- crop(r, e)
      r <- crop(r, rCrop)
      r[is.na(rCrop)] <- NA
      writeRaster(r, file = paste0(simID, "/landtypes.tif"),
                  datatype='INT4S', overwrite=TRUE, NAflag = 0)
      
    } else {
      # initial communities
      file.copy(paste0(inputDir, "/initial-communities_", areaName, ".tif"),
                paste0(simID, "/initial-communities.tif"),
                overwrite = T)
      # landtypes
      file.copy(paste0(inputDir, "/landtypes_",
                       areaName, ".tif"),
                paste0(simID, "/landtypes.tif"),
                overwrite = T)
    }
    
    if(includeSnags) {
      ### storing info to exclude adjust living cohort appropriately
      fName <- paste0(inputDir, "/initial-snags_", areaName, ".txt")
      snags <- read.table(fName, skip = 3,
                        header = FALSE, comment.char = ">")
      file.copy(paste0(fName),
                paste0(simID, "/initial-snags.txt"),
                overwrite = T)
    }

    #############################################
    #############################################
    ### adjust living cohort to snag cohorts
    ## 
    
    #############################################
    #############################################
    
    if(smoothAgeClasses) {
      spp <- read.table(paste0(inputDir, "/species_",areaName, ".txt"), skip = 1,
                         header = FALSE, comment.char = ">")[,1]
      
      initComm <- paste0(inputDir, "/initial-communities_",
                         areaName, ".txt")
      x <- readLines(initComm)
      tmp <- strsplit(x, " ")
      tmp <- lapply(tmp, function(x) x[which(nchar(x)>0)])
      
     
      for (j in seq_along(tmp)) {
        fname <- paste0(simID, "/initial-communities.txt")
        if(j == 1) {
          file.create(fname)
        }
        l <- tmp[[j]]
  
        if(l[1] %in% spp) {
          sp <- l[1]
          cohortAge <- as.numeric(l[-1])
          cohortAge <- round(cohortAge + runif(length(cohortAge), min = -9, max = 0))
          l <- paste0(sp, "\t", paste(cohortAge, collapse = " "))
    
        } else {
          l <- paste0(l, collapse = " ")
        }
        write(l, file = fname,
              append = T)
      }
     
      
    } else {
      file.copy(paste0(inputDir, "/initial-communities_",
                       areaName, ".txt"),
                paste0(simID, "/initial-communities.txt"),
                overwrite = T)
    }
    
    file.copy(paste0(inputDir, "/landtypes_",
                     areaName, ".txt"),
              paste0(simID, "/landtypes.txt"),
              overwrite = T)
    
    ###############################################
    ### Succession extension
    
    # ForC-succession
    
    if(spinup) {
        file.copy(paste0(inputDir, "/forCS-input_",
                         areaName, "_spinup.txt"),
                  paste0(simID, "/forCS-input.txt"),
                  overwrite = T)
    } else {
        file.copy(paste0(inputDir, "/forCS-input_",
                         areaName,"_", scenario, ".txt"),
                  paste0(simID, "/forCS-input.txt"),
                  overwrite = T)
    }
    
    
    # Climate inputs
    file.copy(paste0(inputDir, "/forCS-climate_",
                     areaName, "_", scenario, ".txt"),
              paste0(simID, "/forCS-climate.txt"),
              overwrite = T)
    
    # ForC_DM 
    if(as.numeric(forCSVersion) >=3.1){
      file.copy(paste0(inputDir, "/ForCS_DM_",
                       areaName, "_", scenario, ".txt"),
                paste0(simID, "//ForCS_DM.txt"),
                overwrite = T)
    }
      
      

    
    ###############################################
    ### Disturbances
    if(!spinup) {
      if(harvest) {
        if(cropped) {
          
          ### Harvesting
          # stand map
          r <- raster(paste0(inputDir, "/stand-map_",
                             areaName, ".tif"))
          # r <- crop(r, e)
          r <- crop(r, rCrop)
          r[is.na(rCrop)] <- NA
          writeRaster(r, file = paste0(simID, "/stand-map.tif"),
                      datatype='INT4S', overwrite=TRUE, NAflag = 0)
          # management areas
          r <- raster(paste0(inputDir, "/mgmt-areas_",
                             areaName, ".tif"))
          r <- crop(r, rCrop)
          r[is.na(rCrop)] <- NA
          writeRaster(r, file = paste0(simID, "/mgmt-areas.tif"),
                      datatype='INT4S', overwrite=TRUE, NAflag = 0)
          
        } else {
          # stand map
          file.copy(paste0(inputDir, "/stand-map_",
                           areaName, ".tif"),
                    paste0(simID, "/stand-map.tif"),
                    overwrite = T)
          
          # management areas
          file.copy(paste0(inputDir, "/mgmt-areas_",
                           areaName, ".tif"),
                    paste0(simID, "/mgmt-areas.tif"),
                    overwrite = T)
          
          
        }
        # ### Harvesting
        
        # base-harvest.txt
        ## add a loop here to copy past each biomass havest-i to each corresponding folder simID
        #file.copy(paste0(inputDir, "/biomass-harvest_",
         #                areaName, "_", mgmt, ".txt"),
        #          paste0(simID, "/biomass-harvest.txt"),
         #         overwrite = T)
        for (j in 1:(max_simID+1)) {
          # source file
          src <- file.path(dir_path, paste0("biomass-harvest-", j, ".txt"))
          # target simID folder is j-1 (since simID starts at 0)
          sim_folder <- file.path(getwd(), as.character(j - 1))
          if (!dir.exists(sim_folder)) {
            warning("Skipping j=", j, ": folder ", sim_folder, " does not exist")
            next
          }
          # destination path
          dest <- file.path(sim_folder, "biomass-harvest.txt")
          file.copy(src, dest, overwrite = TRUE)
        }
      }
    
      if(wind) {
        ### Wind
        file.copy(paste0(inputDir, "/base-wind_",
                         areaName, ".txt"),
                  paste0(simID, "/base-wind.txt"),
                  overwrite = T)
      }
 
      if(fire) {
        ### Fire
        file.copy(paste0(inputDir, "/base-fire_",
                         areaName, "_", scenario, ".txt"),
                  paste0(simID, "/base-fire.txt"),
                  overwrite = T)
        if(scenario == "baseline") {
          yrs <- 0
        } else {
          yrs <- c(10, 40, 70)
        }
        for (y in yrs) {
          if(cropped) {
            r <-  raster(paste0(inputDir, "/fire-regions_",
                                areaName, "_", scenario, "_", y, ".tif"))
            r <- crop(r, rCrop)
            r[is.na(rCrop)] <- NA
            writeRaster(r, file = paste0(simID, "/fire-regions_",
                                         y, ".tif"),
                        datatype='INT4S', overwrite=TRUE, NAflag = 0)
          } else {
            file.copy(paste0(inputDir, "/fire-regions_",
                             areaName, "_", scenario, "_", y, ".tif"),
                      paste0(simID, "/fire-regions_",
                             y, ".tif"),
                      overwrite = T)
          }
          
        }
      }
      
      if (BDA) {
        # original file, we need to import all other parameters of BDA model to be the same in the new create files
        original_content <- readLines(file.path(inputDir, "Base-BDA_budworm.txt"))
        # The base-Bda file well need to be updated each time base on the number created of Base-BDA_budworm
        base_bda_template <- readLines(file.path(inputDir, "base-bda.txt"))
        bda_input_files_section_index <- grep("BDAInputFiles", base_bda_template) 
        
        
        # Function to update Start and End year in the BDA content
        update_years_in_content <- function(content, start_year, duration) {
          end_year <- start_year + duration
          content <- sub("StartYear \\d+", paste("StartYear", start_year), content)
          content <- sub("EndYear \\d+", paste("EndYear", end_year), content)
          return(content)
        }
        
        
        tsle <- t0 - 2006 ### here we should indicate the year where 'mortality' really started, not only defoliation
        
        # Initialize parameters for the BDA file 
        start_year <- round(rnorm(1,
                                  mean = bdaParams$cycleMean,
                                  sd = bdaParams$cycleSD)) - tsle
        
        
        BDAduration <- round(rnorm(1,
                                   mean = bdaParams$durMortMean,
                                   sd = bdaParams$durMortSD))
        
        current_year <- start_year + BDAduration  # Initialize current_year based on first outbreak
        j <- 1
        
        # List to store generated BDA file names
        bda_file_names <- character()
        
        # Generate BDA files until we reach the end of simulation  
        while(start_year <= simDuration) {
          
          if(j >1) {
            BDAduration <- round(rnorm(1,
                                       mean = bdaParams$durMortMean,
                                       sd = bdaParams$durMortSD))
          }
          updated_content <- update_years_in_content(original_content, start_year, BDAduration)
          new_bda_file_name <- sprintf("Base-BDA_budworm-%d.txt", j)
          new_bda_file_path <- file.path(simID, new_bda_file_name)
          writeLines(updated_content, new_bda_file_path)
          
          # Append the new BDA file name to the list
          bda_file_names <- append(bda_file_names, new_bda_file_name)
          
          # Prepare for the next iteration
          start_year <- start_year + round(rnorm(1,
                                                 mean = bdaParams$cycleMean,
                                                 sd = bdaParams$cycleSD))
          
          
           if (start_year > simDuration) {
             break
           }
          
          j <- j + 1
        }
        
        # Update the base BDA template to include all BDA file names
        bda_input_files_lines <- paste0(bda_file_names, collapse = "\n\t")
        base_bda_template[bda_input_files_section_index] <- paste0("BDAInputFiles\t",
                                                                   bda_input_files_lines)
        
        # Save the updated base BDA configuration file
        writeLines(base_bda_template, file.path(simID, "base-bda.txt"))
      }
      
      ### scenarios file
      scenFile <- paste0(inputDir, "/scenario.txt")
      x <- readLines(scenFile)
      if(!fire) {
        index <- grep("Base Fire", x)
        x[index] <- paste(">>", x[index])
      }
      if(!harvest) {
        index <- grep("Biomass Harvest", x)
        x[index] <- paste(">>", x[index])
      }
      if(!wind) {
        index <- grep("Base Wind", x)
        x[index] <- paste(">>", x[index])
      }
      if(!BDA) {
        index <- grep("Base BDA", x)
        x[index] <- paste(">>", x[index])
      }
      
      index <- grep("Duration", x)
      x[index] <- paste("Duration", simDuration)
      ###############################################
      ### scenario.txt
      sink(paste(simID, "scenario.txt", sep = "/"))
      writeLines(x)
      sink()
    } else {
        ###############################################
        ### scenario.txt
        file.copy(paste0(inputDir, "/scenario_spinup.txt"),
                  paste0(simID, "/scenario.txt"),
                  overwrite = T)
    }
    
    ### species.txt
    file.copy(paste0(inputDir, "/species_",
                     areaName, ".txt"),
              paste0(simID, "/species.txt"),
              overwrite = T)
    
    write.table(t(simInfo[i,]), file = paste0(simID, "/README.txt"),
                quote = F, col.names = F)
    
}

stopCluster(cl)
