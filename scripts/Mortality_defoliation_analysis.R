## organisation ECCC
### Ameray Abderrahmane

######  mortality analysis from East (Location 1) to West (location 7)
######################################################################
######################################################################
# Clear workspace
rm(list = ls())

library(raster)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set working directory to rasters dir
setwd("D:/wannaproject/mortalityRasters")

# Function to process rasters for a single location
process_location_mortality <- function(location) {
  years <- 1967:1992
  propArea <- numeric()
  propAreaCumul <- numeric()
  outbreak_start <- NA
  outbreak_end <- NA
  first_positive_year_found <- FALSE
  
  # Load ForestLand raster
  forestland_path <- sprintf("D:/wannaproject/mortalityRasters/ForestLand_%d.tif", location)
  if (!file.exists(forestland_path)) {
    stop(paste("ForestLand raster for location", location, "not found."))
  }
  
  reference_raster <- raster(forestland_path)
  total_forest_cells <- sum(values(reference_raster) == 1, na.rm = TRUE)  # Total forest cells
  
  if (total_forest_cells == 0) {
    stop(paste("No forest cells found in ForestLand raster for location", location))
  }
  
  cumulative_raster <- NULL  # Initialize the cumulative raster for max calculations
  for (year in years) {
    # Load mortality raster
    file_path <- sprintf("D:/wannaproject/mortalityRasters/mortality_%d_%d.tif", year, location)
    if (file.exists(file_path)) {
      r <- raster(file_path)
      if (!compareRaster(r, reference_raster, stopiffalse = FALSE)) {
        r <- resample(r, reference_raster, method = "ngb")  # Align raster
      }
    } else {
      # Create empty raster with the same alignment as ForestLand raster
      r <- raster(nrows = nrow(reference_raster), ncols = ncol(reference_raster), 
                  ext = extent(reference_raster), crs = crs(reference_raster))
      r[] <- 0
      warning(sprintf("Mortality raster for year %d location %d not found. Using empty raster.", year, location))
    }
    
    # Consider only mortality levels 11 to 15 and apply proportions
    mortality_proportions <- c("11" = 0.125, "12" = 0.375, "13" = 0.625, "14" = 0.875, "15" = 1)
    r[] <- mortality_proportions[as.character(values(r))]
    r[!(values(r) %in% c(0.125, 0.375, 0.625, 0.875, 1))] <- 0
    r[is.na(r)] <- 0  # Replace NA with 0
    
    # Calculate weighted area affected for this year
    weighted_affected_area <- sum(values(r * reference_raster), na.rm = TRUE)
    yearly_prop <- weighted_affected_area / total_forest_cells
    propArea <- append(propArea, yearly_prop)
    
    # Update cumulative raster and calculate cumulative weighted proportion
    if (is.null(cumulative_raster)) {
      cumulative_raster <- r
    } else {
      cumulative_raster <- calc(stack(cumulative_raster, r), fun=max)
    }
    cumulative_weighted_area <- sum(values(cumulative_raster), na.rm = TRUE)
    propAreaCumul <- append(propAreaCumul, cumulative_weighted_area / total_forest_cells)
    
    if (yearly_prop > 0 && !first_positive_year_found) {
      outbreak_start <- year
      first_positive_year_found <- TRUE
    }
    if (yearly_prop > 0) {
      outbreak_end <- year
    }
  }
  
  outbreak_duration <- if (is.na(outbreak_start)) 0 else outbreak_end - outbreak_start + 1
  
  return(data.frame(
    Year = years,
    propArea = propArea * 100,
    propAreaCumul = propAreaCumul * 100,
    Location = location,
    OutbreakDuration = rep(outbreak_duration, length(years))
  ))
}

# Main execution: Process all locations and collect the results
all_locations_mortality <- bind_rows(lapply(1:7, function(location) {
  cat("Processing location", location, "...\n")
  return(process_location_mortality(location))
}))

# Calculate statistics and reshape data for the table
summary_data <- all_locations_mortality %>% 
  group_by(Location) %>% 
  summarise(
    Min = min(propArea),
    Max = max(propArea),
    Mean = mean(propArea),
    MinCumul = min(propAreaCumul),
    MaxCumul = max(propAreaCumul),
    MeanCumul = mean(propAreaCumul),
    OutbreakDuration = unique(OutbreakDuration)
  )

# Print summary of data
print(summary_data)

# Calculate overall statistics across all locations
overall_stats <- summary_data %>% 
  summarise(
    OverallMin = min(Mean),
    OverallMax = max(Mean),
    verallmean = mean(Mean),
    OverallMinCumul = min(MaxCumul),
    OverallMaxCumul = max(MaxCumul),
    OverallMeanCumul = mean(MaxCumul),
    OverallMinDuration = min(OutbreakDuration),
    OverallMaxDuration = max(OutbreakDuration),
    OverallMeanDuration = mean(OutbreakDuration),
  )

# Print overall summary of data
print(overall_stats)



# Convert mortality data to long format and add data type
all_locations_mortality_long <- all_locations_mortality %>%
  pivot_longer(cols = c("propArea", "propAreaCumul"), names_to = "Variable", values_to = "Value") %>%
  mutate(datatype = "mortality")

combined_plot <- ggplot(all_locations_mortality_long, aes(x = Year, y = Value, color = Variable)) +
    geom_line() +
    facet_wrap(~ Location + Variable, scales = "free_y", ncol = 2) +
    labs(
      title = "Mortality Statistics Across Locations",
      x = "Year", y = "Proportion (%)",
      color = "Statistic"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "top"
    )
print(combined_plot)
  # Save the plot
ggsave("D:/wannaproject/mortalityRasters/plots/all_mortality_plot.png", plot = combined_plot, width = 8, height = 12)
print(combined_plot)






######  defoliation analysis fro m East (Location 1) to West (location 7)
######################################################################
#######################################################################

library(raster)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set working directory
setwd("D:/wannaproject/mortalityRasters")

# Function to process rasters for a single location
process_location_defoliation <- function(location) {
  years <- 1969:1987
  propArea <- numeric()
  propAreaCumul <- numeric()
  
  # Load ForestLand raster
  forestland_path <- sprintf("D:/wannaproject/mortalityRasters/ForestLand_%d.tif", location)
  if (!file.exists(forestland_path)) {
    stop(paste("ForestLand raster for location", location, "not found."))
  }
  
  reference_raster <- raster(forestland_path)
  if(is.null(reference_raster) || length(values(reference_raster)) == 0) {
    stop("Reference raster is NULL or empty.")
  }
  
  total_forest_cells <- sum(values(reference_raster) == 1, na.rm = TRUE)
  if (total_forest_cells == 0) {
    stop(paste("No forest cells found in ForestLand raster for location", location))
  }
  
  cumulative_raster <- NULL  # Initialize the cumulative raster for max calculations
  for (year in years) {
    file_path <- sprintf("D:/wannaproject/mortalityRasters/defoliation_%d_%d.tif", year, location)
    if (file.exists(file_path)) {
      r <- raster(file_path)
    } else {
      # Create empty raster with the same alignment as ForestLand raster
      r <- raster(nrows = nrow(reference_raster), ncols = ncol(reference_raster), 
                  ext = extent(reference_raster), crs = crs(reference_raster))
      values(r) <- 0
      warning(sprintf("Defoliation raster for year %d location %d not found. Using empty raster.", year, location))
    }
    
    if(!is.null(r) && !compareRaster(reference_raster, r)) {
      r <- resample(r, reference_raster, method = "ngb")
    }
    
    # Apply only values within defoliation classes
    r[!(values(r) %in% 1:3)] <- 0
    r[is.na(r)] <- 0
    
    # Calculate the number of affected cells
    affected_area <- sum(values(r) %in% c(1, 2, 3), na.rm = TRUE)
    yearly_prop <- (affected_area / total_forest_cells)*100
    propArea <- append(propArea, yearly_prop)
    
    # Update cumulative raster and calculate cumulative proportion
    if (is.null(cumulative_raster)) {
      cumulative_raster <- r
    } else {
      cumulative_raster <- calc(stack(cumulative_raster, r), fun=max)
    }
    
    cumulative_affected_area <- sum(values(cumulative_raster) %in% c(1, 2, 3), na.rm = TRUE)
    propAreaCumul <- append(propAreaCumul, (cumulative_affected_area / total_forest_cells)*100)
  }
  
  return(data.frame(
    Year = years,
    propArea = propArea,
    propAreaCumul = propAreaCumul,
    Location = location
  ))
}


# Main execution: Process all locations and collect the results
all_locations_defoliation <- bind_rows(lapply(1:7, function(location) {
  cat("Processing location", location, "...\n")
  return(process_location_defoliation(location))
}))

# Calculate statistics and reshape data for the table
summary_data <- all_locations_defoliation %>%
  group_by(Location) %>%
  summarise(
    Min = min(propArea),
    Max = max(propArea),
    Mean = mean(propArea),
    MinCumul = min(propAreaCumul),
    MaxCumul = max(propAreaCumul),
    MeanCumul = mean(propAreaCumul)
  )

# Print summary of data
print(summary_data)
# Plot defoliation data

# Convert defoliation data to long format and add data type
all_locations_defoliation_long <- all_locations_defoliation %>%
  pivot_longer(cols = c("propArea", "propAreaCumul"), names_to = "Variable", values_to = "Value") %>%
  mutate(datatype = "defoliation")
all_locations_defoliation_long
combined_plot2 <- ggplot(all_locations_defoliation_long, aes(x = Year, y = Value, color = Variable)) +
    geom_line() +
    facet_wrap(~ Location + Variable, scales = "free_y", ncol = 2) +
    labs(
      title = "Defoliation Statistics Across Locations",
      x = "Year", y = "Proportion (%)",
      color = "Statistic"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "top"
    )
print(combined_plot2) 
  # Save the plot
ggsave("D:/wannaproject/mortalityRasters/plots/all_defoliation_plot.png", plot = combined_plot2, width = 8, height = 12)





# Calculate yearly mean and standard deviation across all locations for mortality and defoliation
datasetMD_combined <- bind_rows(all_locations_mortality_long, all_locations_defoliation_long)
datasetMD_combined
yearly_stats_combined <- datasetMD_combined %>%
  group_by(Year, datatype, Variable) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    .groups = "drop"
  )
yearly_stats_combined

# Plot with confidence intervals and separate by data type and variable
plot3<-ggplot(yearly_stats_combined, aes(x = Year, y = mean, group = Variable)) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = Variable), alpha = 0.1) +
  geom_line(aes(color = Variable),linewidth=0.5) +
  facet_wrap(Variable~ datatype , scales = "free_y", ncol = 2) +
  labs(title = "Yearly Statistics for Mortality and Defoliation",
       x = "Year",
       y = "Value (%)",
       color = "Variable",
       fill = "Variable") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "top"
  )
plot3
# Save the plot
ggsave("D:/wannaproject/mortalityRasters/plots/mortaltiy_defoliation.png", plot = plot3, width = 6, height = 6)
