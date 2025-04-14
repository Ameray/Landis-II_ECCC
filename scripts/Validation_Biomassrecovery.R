################################################################################
## organisation ECCC
### Ameray Abderrahmane
#####   
### Clear environment
rm(list = ls())

require(data.table)
require(dplyr)
require(raster)
require(stringr)
require(tidyr)
require(ggplot2)

# Set the working directory and source the function for BDA parameters
dir_path <- "D:/Wena_Project/Testing/2025-04-08/"  
setwd(dir_path)
source("D:/Wena_Project/scripts/fetchBDAParamFnc.R") 

### Fetching outputs
a <- "mixedwood-042-51"
d <- c("bda" = 4, "wind" = 3, "harvest" = 2, "fire" = 1)  
simDir <- "D:/Wena_Project/Testing/2025-04-08"
simName <- gsub("simPkg_mixedwood-042-51_", "", basename(simDir))

# Read simulation info
simInfo <- read.csv(paste(simDir, "simInfo.csv", sep = "/"), colClasses = c("simID" = "character"))
simIDs <- simInfo$simID
simIDs <- str_pad(simIDs, width = max(nchar(simIDs)), side = "left", pad = "0")
simIDs
# subdirs
allDirs <- list.dirs(simDir, full.names = FALSE, recursive = FALSE)
dirIndex <- which(simIDs %in% allDirs & simInfo$areaName == a)

outputList <- list()

for (i in dirIndex) {
  # Pull scenario attributes
  simID <- simIDs[i]
  sDir <- file.path(simDir, simID)
  areaName <- simInfo[i, "areaName"]
  scenario <- simInfo[i, "scenario"]
  mgmtScenario <- simInfo[i, "mgmt"]
  mgmtScenarioName <- mgmtScenario
  wind <- simInfo[i,"wind"]
  fire <- simInfo[i,"fire"]
  BDA <- simInfo[i,"BDA"]
  harvest <- simInfo[i, "harvest"]
  BDApara <-simInfo[i,"BDAfile"]   ## this column control the number of combinations of parameters tested 
  replicate <- simInfo[i, "replicate"]
  
  # Only proceed if BDA == TRUE and harvest == FALSE
  if (!harvest && BDA) {
    # Fetch parameters if BDA is used
    bdaParam <- fetchBDAParam1(file.path(sDir, "Base-BDA_budworm1.txt"))
    bdaSpp <- bdaParam$BDASpeciesParameters$SpeciesName
    
    # Species levels (read from the LANDIS species.txt)
    sppLvls <- read.table(file.path(sDir, "species.txt"), skip = 1, comment.char = ">")[,1]
    if (is.factor(sppLvls)) sppLvls <- levels(sppLvls)
    
    # Read the flux biomass file
    fluxBio <- fread(file.path(sDir, "log_FluxBio.csv"))
    
    # Filter fluxBio to BDA disturbance only (Dist==4) and host species
    dfBDA <- fluxBio %>%
      filter(Dist == 4, species %in% bdaSpp) %>%
      group_by(Time, row, column, ecoregion, species) %>%
      mutate(
        AGBtoDOM_woody = MERCH_ToDOM + OtherWoody_ToDOM,
        FOL_ToDOM      = FOL_ToDOM,
        BGBtoDOM       = CrsRt_ToDOM + FRt_ToDOM,
        bioToAir       = MERCH_ToAir + FOL_ToAir + OtherWoody_ToAir + CrsRt_ToAir + FRt_ToAir
      ) %>%
      summarise(
        AGBtoDOM_woody = sum(AGBtoDOM_woody),
        FOL_ToDOM      = sum(FOL_ToDOM),
        BGBtoDOM       = sum(BGBtoDOM),
        bioToAir       = sum(bioToAir),
        bioToFPS       = sum(BioToFPS)
      ) %>%
      ungroup()
    
    dfBDA$cell <- paste(dfBDA$row, dfBDA$column, sep = "-")
    
    # Read total biomass
    Biomass <- fread(file.path(sDir, "log_BiomassC.csv"))
    Biomass$cell <- paste(Biomass$row, Biomass$column, sep = "-")
    
    # Filter biomass to only the same cells (which have BDA disturbance)
    BiomassAffected <- Biomass %>%
      filter(cell %in% dfBDA$cell) %>%
      group_by(Time, row, column, cell, ecoregion, species) %>%
      mutate(
        Time = Time + 1,
        AGB  = Wood + Leaf,
        BGB  = CrsRoot + FineRoot
      ) %>%
      summarise(
        AGB_woody = sum(Wood),
        Leaf      = sum(Leaf),
        BGB       = sum(BGB)
      ) %>%
      ungroup()
    
    # Merge flux-based mortality data with biomass data to compute severity
    df_allspecies <- merge(
      BiomassAffected,
      dfBDA,
      by = c("Time", "row", "column", "cell", "ecoregion", "species"),
      all.x = TRUE
    )
    
    # Identify host vs non-host groups (if needed)
    df_allspecies <- df_allspecies %>%
      mutate(group = ifelse(species %in% bdaSpp, "host", "no-host"))
    
    # Summarize severity by cell
    df_cell <- df_allspecies %>%
      group_by(Time, row, column, cell, ecoregion) %>%
      summarise(
        deadwood = sum(AGBtoDOM_woody, na.rm = TRUE),
        AGB      = sum(AGB_woody, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      # Only cells/timesteps where Dist occurs (deadwood not NA)
      filter(!is.na(deadwood)) %>%
      group_by(Time, row, column, cell, ecoregion) %>%
      summarise(severity = 100 * sum(deadwood) / sum(AGB)) %>%
      ungroup()
    
    # Average severity across all outbreak timesteps for each cell
    df_cell_mean <- df_cell %>%
      group_by(row, column, cell, ecoregion) %>%
      summarise(severity2 = mean(severity[severity != 0], na.rm = TRUE)) %>%
      ungroup()
    
    # Merge back to get a final table for all time steps / species
    dfall <- merge(
      BiomassAffected,
      df_cell_mean,
      by = c("row", "column", "cell", "ecoregion"),
      all.x = TRUE
    ) %>%
      arrange(Time)
    
    ## ----------------------------------------------------------------------
    ## CLASSIFY SEVERITY USING THE 20TH AND 80TH PERCENTILES
    ## (lowest 20% = low, highest 20% = high, remainder = moderate)
    ##the severity classes were determined using the quantiles of the severity index value distribution. 
    #The Low class represents the lowest 20% of severity indices, 
    ###the High class includes the highest 20%,
    #### and the Moderate class, it's the 20% of either side of the median.
    ## ----------------------------------------------------------------------
    quantiles_sim <- quantile(dfall$severity2, probs = c(0.2, 0.8), na.rm = TRUE)
    qs20 <- quantiles_sim[1]
    qs80 <- quantiles_sim[2]
    
    dfall2 <- dfall %>%
      group_by(row, column, cell, ecoregion) %>%
      mutate(
        class = case_when(
          is.na(severity2)            ~ NA,  # e.g., no outbreak or no data
          severity2 <= qs20           ~ "low",
          severity2 >= qs80           ~ "high",
          TRUE                        ~ "moderate"
        )
      ) %>%
      ungroup()
    
    # Add scenario metadata
    dfall2 <- dfall2 %>%
      mutate(
        areaName         = areaName,
        simID            = simID,
        scenario         = scenario,
        mgmtScenario     = mgmtScenario,
        mgmtScenarioName = mgmtScenarioName,
        BDApara          =BDApara,
        replicate        = replicate
      )
    
    # Save to output list
    out_list <- list()
    out_list[["dfall"]] <- dfall2
    outputList[[i]] <- out_list
  }
}

all_sims_combined <- do.call(rbind, lapply(outputList, `[[`, "dfall"))
head(all_sims_combined)
save(all_sims_combined, file = paste0("Biomassall", simName, ".RData"))
## distribution by class:
ggplot(all_sims_combined, aes(x = class)) +
  geom_bar() +
  facet_wrap(~simID, scales = "free_y") +
  labs(title = "Distribution of severity classes across simIDs",
       x = "Severity class", y = "Count of cell-time records")

all_sims_combined
unique(all_sims_combined$class)
summary(all_sims_combined)
## test
filter(all_sims_combined,Time==1, row==100,column==    103)
###dfall<-get(load(paste0(dir_path,"/Biomassall", simName, ".RData")))
#dfall
#bdaSpp<-c("ABIE.BAL" ,"PICE.GLA", "PICE.MAR")
dfall2 <-all_sims_combined %>%  
  mutate(group = ifelse(species %in% bdaSpp, "host", "no-host"))

filter(dfall2,Time==1, cell=="100-103")

dfall3 <-dfall2 %>%
  group_by(simID ,scenario, mgmtScenario,BDApara,replicate, cell,ecoregion, Time, class,   group) %>%
  summarise(AGB=sum(AGB_woody+Leaf),TotalBiomass=sum(AGB_woody+Leaf+BGB),BGB=sum(BGB))  # sum of all species in the cell
head(dfall3)
summary(dfall3)
filter(dfall3,Time==1, cell=="100-103")

# sum the biomass for 'host' and 'no-host' for each time and class
aggregated_data <- dfall3 %>%
  group_by(simID ,scenario, mgmtScenario,BDApara,replicate,cell,ecoregion, Time, class) %>%
  summarize(AGB = sum(AGB),
            BGB=sum(BGB),
            TotalBiomass=sum(TotalBiomass), .groups = 'drop') %>%
  mutate(group = "allspecies")
aggregated_data
filter(aggregated_data,Time==1, cell=="100-103")
# Bind this new data back to the original data frame
alldata<- bind_rows(dfall3, aggregated_data)
alldata
summary(alldata)
### calculate the means for replicates
alldatamean<- alldata %>% 
  group_by(scenario, mgmtScenario,BDApara,cell,ecoregion, Time, class,group) %>%
  summarize(AGB = mean(AGB),
            BGB=mean (BGB),
            TotalBiomass= mean(TotalBiomass),
            .groups = 'drop')
alldatamean
# Calculate mean of all cells, standard deviation, and confidence intervals  only the 1st outbreaks
meandata <- alldatamean %>%  ##
  group_by(scenario, mgmtScenario,BDApara, Time,group, class) %>%
  summarise(
    meanAGB = mean(AGB),
    sdAGB = sd(AGB),
    n = n(),  #
    semAGB = sdAGB / sqrt(n),  # 
    lowerCI = meanAGB - 1.96 * semAGB,  
    upperCI = meanAGB + 1.96 * semAGB,  
    .groups = 'drop'
  )
meandata
###  

### filter best comb founded in claibration_Severity.R script

meandata2 <- meandata  %>% filter(BDApara==1)
meandata2
unique(meandata2$BDApara)
#plot
#  geom_area
meandata2$class <- factor(meandata2$class,
                          levels = c("low", "moderate", "high"))
BDAplot2 <- ggplot(meandata2, aes(x = Time, y = meanAGB/100, group = group, fill = as.factor(group))) +
  geom_area(position = "stack", alpha = 0.5) +
  facet_grid( ~ class) +
  labs(x = "Year", y = "AGB (tC/ha/yr)") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 12)
  ) +
  scale_fill_manual(values = c("allspecies" = "darkgreen", "no-host" = "gray", "host" = "red"),
                    labels=c("All species", "Host", "No Host"))
BDAplot2
# Save the plot
ggsave(filename = "Abundance.png", 
       path = "D:/dom_plots/", 
       plot = BDAplot2, 
       width = 7, height = 2.5, dpi = 300)
# Function to find continuous periods in a vector
find_periods <- function(time_vector) {
  diffs <- diff(time_vector)
  breaks <- which(diffs != 1)
  starts <- c(time_vector[1], time_vector[breaks + 1])
  ends <- c(time_vector[breaks], tail(time_vector, 1))
  data.frame(start = starts, end = ends)
}
meandata2
## select best bdapara  for validation
BDApara <- 1
simIDs_sub <- simInfo$simID[ simInfo$BDAfile == BDApara ]
simID_to_process <- simIDs_sub[1]  ## 

periods_list <- list()
for (simID in simID_to_process) {
  sDir <- file.path(simDir, simID, "BDA")
  bda_log_file <- file.path(sDir, "bda-log.csv")
  
  if (!file.exists(bda_log_file)) {
    next
  }
  
  BDAlog <- read.csv(bda_log_file)
  periods <- find_periods(BDAlog$Time)  # your find_periods() function
  
  period_label <- paste("period", seq_along(periods$start), sep = "_")
  
  period_data <- data.frame(
    simID = rep(simID, length(periods$start)),
    period = period_label,
    start = periods$start,
    end   = periods$end
  )
  
  periods_list[[as.character(simID)]] <- period_data
}

# Combine them
all_periods_data <- do.call(rbind, periods_list)
all_periods_data


### use patchwork libraries to create plots
library(patchwork)

# Initialize an empty list to store the plots
plot_list <- list()

# Loop through each simID and generate the plot
for (id in BDApara) {
  period_info <- filter(all_periods_data, BDApara == as.character(id), period == "period_1")
  if (nrow(period_info) > 0) {
    start_of_period1 <- period_info$start[1]
    end_of_period1 <- period_info$end[1]
    
    meandata_filtered <- filter(meandata2,  BDApara == as.character(id)) %>%
      mutate(Time2 = Time - start_of_period1)
    
    plot <- ggplot(meandata_filtered, aes(x = Time2, y = meanAGB/100, group = class, color = as.factor(class))) +
      geom_line(size = 0.5, alpha = 1) +
      geom_ribbon(aes(ymin = lowerCI/100, ymax = upperCI/100, fill = as.factor(class)), color = NA, alpha = 0.1) +
      facet_grid(BDApara~ group, scales = "free") +
      labs(x = "Years Since Outbreak", y = "tC/ha/yr") +
      theme_minimal() +
      theme(legend.position = "none") +  # Remove legend from individual plots
      scale_color_manual(values = c("low" = "darkgreen", "moderate" = "brown", "high" = "red")) +
      scale_fill_manual(values = c("low" = "darkgreen", "moderate" = "brown", "high" = "red")) +
      annotate("rect", xmin = 0, xmax = end_of_period1 - start_of_period1 + 1, ymin = -Inf, ymax = Inf, alpha = 0.5, fill = "gray")
    
    # Add the plot to the list
    plot_list[[as.character(id)]] <- plot
  }
}

# Combine all plots 
combined_plot <- wrap_plots(plot_list, ncol = 1) +
  plot_layout(guides = 'collect') &  
  theme(legend.position = "bottom")  

print(combined_plot)
meandata_filtered
unique(meandata_filtered$BDApara)


####################### filter plot to compare with Migue plots for biomass recovery  

meandata_filtered1 <-meandata_filtered %>%
                    filter(BDApara=="1",group=="allspecies")   ### 

meandata_filtered1
common_ylim<-c(0,80)
# Simulated plot with fixed y-axis range
plot_simulated <- ggplot(meandata_filtered1, aes(x = Time2, y = meanAGB/100, group = class, color = as.factor(class))) +
  geom_line(size = 0.5, alpha = 1) +
  geom_ribbon(aes(ymin = lowerCI/100, ymax = upperCI/100, fill = as.factor(class)), color = NA, alpha = 0.1) +
  labs(x = "Time Since Outbreak", y = "tC/ha", title = "Simulated") +
  theme_minimal() +
  coord_cartesian(ylim = common_ylim) +  # Set common y-axis
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "none"
  ) +
  scale_color_manual(values = c("low" = "darkgreen", "moderate" = "brown", "high" = "red")) +
  scale_fill_manual(values = c("low" = "darkgreen", "moderate" = "brown", "high" = "red")) +
  annotate("rect", xmin = 0, xmax = 9,
           ymin = -Inf, ymax = Inf, alpha = 0.5, fill = "gray")
plot_simulated
################## emperical obseved biomass recovery from migue dataset
BiomassR_emp<-read.csv("D:/dom_plots/Biomassrecovery_emp.csv")
head(BiomassR_emp)
BiomassR_emp1<-BiomassR_emp %>% filter (BIO_ZONE=="Balsam_fir-white_birch")
BiomassR_emp1
# Filter data living carbon
live_data <- BiomassR_emp1 %>% filter(Carbon_Type == "Live_Carbon")
live_data
# Main plot

# Empirical plot with same y-axis range
plot_empirical <- ggplot() +
  geom_smooth(
    data = live_data,
    aes(x = TSO_yrs, y = Mean_Carbon, color = SEVERITY_cl, fill = SEVERITY_cl),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    size = 1,
    alpha = 0.1
  ) +
  coord_cartesian(ylim = common_ylim) +  # Set common y-axis
  scale_color_manual(values = c("Low" = "darkgreen", "Moderate" = "brown", "High" = "red")) +
  scale_fill_manual(values = c("Low" = "darkgreen", "Moderate" = "brown", "High" = "red")) +
  labs(x = "Time since Outbreak", y = "tC/ha", title = "Empirical") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right"
  ) +
  annotate("rect", xmin = 0, xmax = 9,
           ymin = -Inf, ymax = Inf, alpha = 0.5, fill = "gray")
plot_empirical
# Combine
plotCompar<-plot_simulated + plot_empirical
plotCompar
# Save the plot
ggsave(filename = "biomassrecovery.png", 
       path = "D:/dom_plots/", 
       plot = plotCompar, 
       width = 7, height = 2.5, dpi = 300)

