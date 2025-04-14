### Clear environment
################################################################################
## organisation ECCC
### Ameray Abderrahmane

rm(list = ls())

require(data.table)
require(dplyr)
require(raster)
require(stringr)
library(tidyr)
library(ggplot2)

### directory 
dir_path <- "D:/Wena_Project/Testing//Landscape/2025-04-10"
setwd(dir_path)
source("D:/Wena_Project/scripts/fetchBDAParamFnc.R")

a <- "mixedwood-042-51"
d <- c("BDA" = 4, "wind" = 3, "harvest" = 2, "fire" = 1)  # all  disturbance types
simDir <- dir_path
simName <- gsub("simPkg_mixedwood-042-51_", "", basename(simDir))
simInfo <- read.csv(file.path(simDir, "simInfo.csv"), 
                    colClasses = c("simID" = "character"))
simInfo
# 
simIDs <- str_pad(simInfo$simID,
                  width = max(nchar(simInfo$simID)),
                  side = "left",
                  pad = "0")

# Identify relevant directories
allDirs <- list.dirs(simDir, full.names = FALSE, recursive = FALSE)
dirIndex <- which(simIDs %in% allDirs & simInfo$areaName == a)
dirIndex

list_dist <- list()

# 2) Loop over each disturbance name in 'd'
for (distName in names(d)) {
  
  outputList <- list()
  
  # Dist code 
  distCode <- d[[distName]]
  
  for (i in dirIndex) {
    
    simID <- simIDs[i]
    sDir  <- file.path(simDir, simID)
    
    scenario     <- simInfo[i, "scenario"]
    mgmtScenario <- simInfo[i, "mgmt"]
    replicate    <- simInfo[i, "replicate"]
    areaName     <- simInfo[i, "areaName"]
    mgmtScenarioName <- mgmtScenario
    
    #Read fluxBio for this simID
    fluxBio_raw <- fread(file.path(sDir, "log_FluxBio.csv"))
    
    # If there are no rows, skip
    if (nrow(fluxBio_raw) == 0) next
    
    # Aggregate across species so that each (Time, row, column, Dist) 
    fluxBio <- fluxBio_raw %>%
      group_by(Time, row, column, Dist) %>%
      summarise(
        ToDOM = sum(MERCH_ToDOM + OtherWoody_ToDOM),
        TOAir = sum(MERCH_ToAir + FOL_ToAir + OtherWoody_ToAir + CrsRt_ToAir + FRt_ToAir),
        ToFPS = sum(BioToFPS),
        .groups = "drop"
      ) %>%
      mutate(
        cell = paste(row, column, sep = "-"),
        isHarvest = (Dist == 2 & ToFPS > 0),
        isBDA     = (Dist == 4 & ToDOM  > 0),
        isWind    = (Dist == 3 & ToDOM  > 0),
        isFire    = (Dist == 1 & TOAir  > 0)
      )
    
    # 2C) If the simulation does NOT actually include the disturbance we're focusing on, skip
    if (!simInfo[i, distName]) {
      next
    }
    
    # 3) Summarize how many times each cell was disturbed by each type
    distSummary <- fluxBio %>%
      group_by(cell) %>%
      summarise(
        nHarvest = sum(isHarvest),
        nBDA     = sum(isBDA),
        nWind    = sum(isWind),
        nFire    = sum(isFire),
        .groups  = "drop"
      )
    
    # 4) According to the chosen 'distName', we filter exactly 1 event of that type 
    #    and 0 events of the others.
    if (distName == "harvest") {
      distSummary_oneEventOnly <- distSummary %>%
        filter(nHarvest == 1, nBDA == 0, nWind == 0, nFire == 0)
    } else if (distName == "BDA") {
      distSummary_oneEventOnly <- distSummary %>%
        filter(nBDA == 1, nHarvest == 0, nWind == 0, nFire == 0)
    } else if (distName == "wind") {
      distSummary_oneEventOnly <- distSummary %>%
        filter(nWind == 1, nHarvest == 0, nBDA == 0, nFire == 0)
    } else if (distName == "fire") {
      distSummary_oneEventOnly <- distSummary %>%
        filter(nFire == 1, nHarvest == 0, nBDA == 0, nWind == 0)
    }
    
    # If no cell matches, skip this sim
    if (nrow(distSummary_oneEventOnly) == 0) next
    
    cells_oneEvent <- distSummary_oneEventOnly$cell
    
    # 5) Extract the single DistYear
    distCol <- switch(
      distName,
      "harvest" = "isHarvest",
      "BDA"     = "isBDA",
      "wind"    = "isWind",
      "fire"    = "isFire"
    )
    oneEvent <- fluxBio %>%
      filter(cell %in% cells_oneEvent, .data[[distCol]] == TRUE) %>%
      dplyr::select(cell, Time) %>%
      rename(DistYear = Time)
    
    # 6) Read log_summary
    poolfluxes_raw <- fread(file.path(sDir, "log_summary.csv"))
    
    # Add cell column
    poolfluxes_raw <- poolfluxes_raw %>%
      mutate(cell = paste(row, column, sep = "-"))
    
    # Keep only the one-event cells
    poolfluxes_oneDist <- poolfluxes_raw %>%
      filter(cell %in% cells_oneEvent)
    
    # 7) Join DistYear
    df <- poolfluxes_oneDist %>%
      left_join(oneEvent, by = "cell")
    
    # 8) Compute TimeSinceDist & DistPeriod
    df_final <- df %>%
      mutate(
        TimeSinceDist = Time - DistYear,
        DistPeriod = case_when(
          Time < DistYear ~ "Before",
          Time == DistYear ~ "During",
          Time > DistYear ~ "After",
          TRUE ~ NA_character_
        )
      )
    
    # 9) Add scenario, replicate, mgmt info
    df_all <- df_final %>%
      mutate(
        areaName         = areaName,
        simID            = simID,
        scenario         = scenario,
        mgmtScenario     = mgmtScenario,
        mgmtScenarioName = mgmtScenarioName,
        replicate        = replicate,
        DistType         = distName
      )
    
    # 10) Store in outputList for this simID
    outputList[[i]] <- df_all
  }
  
  # Combine all simIDs for this disturbance
  df_dist <- bind_rows(outputList)
  
  # Save as RData file (unchanged from your script)
  save(df_dist, file = paste0("DistTable_", distName, "_", simName, ".RData"))
  
  # Also store in a list if you want to keep them in memory
  list_dist[[distName]] <- df_dist
}


###############################################################################
## Combine all DistTypes into a single file
###############################################################################
df_allDist <- bind_rows(list_dist)
save(df_allDist, file = paste0("DistTable_allDist_", simName, ".RData"))
###############################################################################


######### visualisation 

### 0. Clean environment & load packages
rm(list = ls())
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

### 1. Define the disturbances
d <- c("harvest", "Fire", "BDA", "wind")

### 2. Loop over each disturbance (UNCHANGED)
for (distName in d) {
  
  # Build the file name, e.g. "DistTable_harvest_2025-04-03.RData"
  distFile <- paste0("DistTable_", distName, "_2025-04-10.RData")
  
  # Load the data (assume it creates an object called df_dist)
  load(distFile)
  
  # We'll assign the loaded object to df_final for clarity.
  df_final <- df_dist
  
  ### 2a) Pivot 'df_final' to long format
  df_long <- df_final %>%
    pivot_longer(
      cols = c("ABio", "BBio", "TotalDOM", "NPP", "Rh", "NEP", "NBP"),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      group = if_else(parameter %in% c("ABio", "BBio", "TotalDOM"),
                      "pool",
                      "flux")
    )
  
  ### 3) Summarize
  df_summary <- df_long %>%
    group_by(simID, scenario, mgmtScenario, TimeSinceDist, parameter, group) %>%
    summarise(
      n = n(),
      mean_val = mean(value/100, na.rm = TRUE),
      sd_val   = sd(value/100, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(
      se_val = sd_val / sqrt(n),
      lower  = mean_val - 1.96 * se_val,
      upper  = mean_val + 1.96 * se_val
    )
  
  ### 4) Split into flux & pool
  df_summary_flux <- df_summary %>% filter(TimeSinceDist>-10, group == "flux")
  df_summary_pool <- df_summary %>% filter(TimeSinceDist>-10, group == "pool")
  
  ### 4a) Build the fluxes plot
  p_fluxes <- ggplot(df_summary_flux, aes(x = TimeSinceDist, y = mean_val, color = parameter)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = parameter),
                alpha = 0.2, color = NA) +
    facet_grid( ~ simID, scales = "free_y") +
    labs(
      x = "Time Since Disturbance",
      y = "tC/ha/yr",
      title = paste("Fluxes vs. Time Since Disturbance -", distName)
    ) +
    scale_fill_discrete(name = "Parameter") +
    scale_color_discrete(name = "Parameter")
  
  outFile_fluxes <- paste0("carbondynamic_fluxes_", distName, ".png")
  ggsave(filename = outFile_fluxes, plot = p_fluxes, width = 8, height = 6, dpi = 300)
  message("Saved fluxes plot for ", distName, " to file: ", outFile_fluxes)
  
  ### 4b) Build the pools plot
  p_pools <- ggplot(df_summary_pool, aes(x = TimeSinceDist, y = mean_val, color = parameter)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = parameter),
                alpha = 0.2, color = NA) +
    facet_grid( ~ simID, scales = "free_y") +
    labs(
      x = "Time Since Disturbance",
      y = "tC/ha",
      title = paste("Pools vs. Time Since Disturbance -", distName)
    ) +
    scale_fill_discrete(name = "Parameter") +
    scale_color_discrete(name = "Parameter")
  
  outFile_pools <- paste0("carbondynamic_pools_", distName, ".png")
  ggsave(filename = outFile_pools, plot = p_pools, width = 8, height = 6, dpi = 300)
  message("Saved pools plot for ", distName, " to file: ", outFile_pools, "\n")
}


# test 
distFile <- paste0("DistTable_harvest_2025-04-08.RData")
harvest <- get(load(distFile))
harvest
unique(harvest$cell)
