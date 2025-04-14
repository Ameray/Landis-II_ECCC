################################################################################
## organisation ECCC
### Ameray Abderrahmane

######################################################################
# SBW Severity Analysis at the entire  landscape only.  turn    cropped=FALSE
######################################################################

# Clear workspace
rm(list = ls())

library(raster)
library(ggplot2)
library(dplyr)
library(gridExtra)  # 
library(rasterVis)  # 
library(readr)
### Integrate simID by setting a dynamic working directory
wwd <- "D:/Wena_Project/Testing/2025-04-03"
wwd
# Reset working directory
setwd(wwd)
#### simID 
simInfo <- read.csv("simInfo.csv", colClasses = c(simID = "character"))
simIDs <- as.integer(simInfo$simID)
simIDs
# Function to process rasters for cumulative severity analysis
process_budworm_severity <- function(years, simID) {
  subdir <- file.path(wwd, simID)
  if (!dir.exists(subdir)) {
    stop(paste("Directory does not exist:", subdir))
  }
  setwd(subdir)  # 
  
  cumulative_raster <- NULL  
  cumulative_affected_area <- numeric(length(years))
  yearly_rasters <- list()   
  
  for (year in years) {
    # Load severity raster
    file_path <- sprintf("bda/budworm1-severity-%d.tif", year)
    if (!file.exists(file_path)) {
      stop(paste("Severity raster for year", year, "not found."))
    }
    
    r <- raster(file_path)
    
    # 
    if (is.null(proj4string(r))) {
      crs(r) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    }
    
    # Classify raster values
    r[r >= 2] <- 2  # Affected forest
    r[r == 1] <- 1  # Forest not affected
    r[r < 1] <- 0   # No forest
    
    # Store raster for current year
    yearly_rasters[[year]] <- r
    
    # Update cumulative raster
    if (is.null(cumulative_raster)) {
      cumulative_raster <- r
    } else {
      cumulative_raster <- overlay(cumulative_raster, r, fun=max)
    }
    
    # Calculate cumulative affected area for this year (values >= 2)
    cumulative_affected_area[year - min(years) + 1] <- sum(values(cumulative_raster) >= 2, na.rm = TRUE)
  }
  
  # Save cumulative affected area data to CSV
  write.csv(data.frame(Year = years, Cumulative_Affected_Area = cumulative_affected_area),
            sprintf("%s/Cumulative_Affected_Area_SimID_%d.csv", wwd, simID),
            row.names = FALSE)
  
  # Return results
  return(list(yearly_rasters = yearly_rasters, cumulative_raster = cumulative_raster, cumulative_affected_area = cumulative_affected_area))
}

# Loop over each simID and process rasters
results <- lapply(simIDs, function(simID) {
  process_budworm_severity(1:7, simID)
})
results
# Reset working directory
setwd(wwd)



# Define base path and SimIDs
basepath <- "D:/Wena_Project/Testing/2025-04-03"
simIDs <- 1:2

# Read and combine CSVs
combined_data <- do.call(rbind, lapply(simIDs, function(simID) {
  filename <- file.path(basepath, sprintf("Cumulative_Affected_Area_SimID_%s.csv", simID))
  df <- read_csv(filename)
  df$SimID <- simID
  return(df)
}))
combined_data1<-combined_data %>% filter(SimID==1)
# Ensure SimID is treated as a factor (for plot legend order)
combined_data$SimID <- factor(combined_data$SimID, levels = simIDs)

# Plot
accArea<-ggplot(combined_data, aes(x = Year, y = 100*Cumulative_Affected_Area/ (9561), color = SimID)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(title = "Cumulative Affected Area over Time by (best simID)",
       x = "Year",
       y = "Affected Area (%",
       color = "Simulation ID") +
  theme_minimal()
accArea
# Save the plot
ggsave(filename = "accArea.png", 
       path = "D:/dom_plots/", 
       plot = accArea, 
       width = 7, height = 4, dpi = 300)



## Plot results for each simID
lapply(seq_along(results), function(i) {
  sim_result <- results[[i]]
  # Define color scale and breaks
  color_scale <- c("white", "green", "red")
  breaks <- c(0, 0.5, 1.5, 3)  # Ensure breaks cover exact integer ranges for 0, 1, and >=2
  
  # Prepare plots for each year with titles
  plots <- lapply(1:length(sim_result$yearly_rasters), function(idx) {
    rasterVis::levelplot(sim_result$yearly_rasters[[idx]], col.regions=color_scale,
                         at=breaks,
                         colorkey=FALSE,  # Disable colorkey for individual plots
                         scales=list(draw=TRUE),  # Enable scales to show axes
                         margin=FALSE,
                         main=paste("Year", idx))  # Title for each year
  })
  
  # Prepare the cumulative plot with its specific title
  cumulative_plot <- rasterVis::levelplot(sim_result$cumulative_raster, col.regions=color_scale,
                                          at=breaks,
                                          colorkey=FALSE,  # Only this plot shows the color key
                                          margin=FALSE,
                                          main="Accumulated")  # Title for cumulative plot
  
  # Add cumulative plot to the list of plots
  plots <- c(plots, list(cumulative_plot))
  
  # Combine plots for each simID into one figure with equal sizing
  plot_grid <- gridExtra::grid.arrange(grobs=plots, ncol=4, widths=rep(1, 4), heights=rep(1, 2))
  # Save the combined plot grid
  ggsave(sprintf("%s/Plots_All_Years_SimID_%d.png", wwd, simIDs[i]), 
         plot_grid,  # Save directly the plot grid
         width = 20, height = 10)
})



