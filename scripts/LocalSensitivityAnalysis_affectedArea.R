## organisation ECCC
### Ameray Abderrahmane

##### this script is applied on the same output of  LSA for severity 
#### this script is adapted for Local sensitivity analysis  of affected area by mortality
#### step 1 prepare inputs for local sensitivity analysis. 
#### change one para by one in BDA model input files


# Clear the environment
rm(list = ls())

#####  use the same output from LocalSensitivityAnalysis_severity.R
#### step 4 read the outputs and calculate affected area  

# Clear workspace
rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Set the base directory path to the same output in step0 for severity 
base_dir <- "D:/Wena_Project/Testing/Test/SA/2024-12-04"

# List all simulation IDs (subdirectories)
simIDs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)

# Initialize a list to store the data
all_data <- list()

# Define the corrected SBW periodic outbreak periods
periods <- list(c(0, 12), c(13, 32))  # #cycle of 32 years  

# Process each simulation independently
for (simID in simIDs) {
  # Define the path to the bda-log.csv file
  bda_log_path <- file.path(base_dir, simID, "bda", "bda-log.csv")
  
  # Check if the file exists
  if (!file.exists(bda_log_path)) {
    message(paste("File not found for simID:", simID))
    next
  }
  
  # Read the bda-log.csv file, ignoring unnamed columns
  df <- read_csv(bda_log_path, show_col_types = FALSE) %>%
    dplyr::select(-starts_with("...")) # Drop unnamed columns like in bda file
  
  # Create a complete sequence of years (Time)
  full_years <- data.frame(Time = seq(1, 32)) # Ensure all years from 1 to n are included
  
  # Merge the full sequence with the original data, filling missing values with 0
  df <- full_years %>%
    left_join(df, by = "Time") %>%
    mutate(DamagedSites = ifelse(is.na(DamagedSites), 0, DamagedSites))
  
  # Initialize an empty data frame for the current simID
  sim_data <- data.frame()
  
  # Calculate cumulative areas for each period
  for (period in periods) {
    start <- period[1]
    end <- period[2]
    
    # Filter data for the given period
    period_df <- df %>%
      filter(Time >= start & Time <= end) %>%
      mutate(CumulativeArea = cumsum(DamagedSites * 6.25)) # Reset cumulative area per period
    
    # Add period and simID information
    period_df$Period <- paste0("Period_", start, "-", end)
    period_df$simID <- simID
    
    # Combine the results for the current period
    sim_data <- bind_rows(sim_plot_data, period_df)
  }
  
  # Store the processed data for the current simID
  all_data[[simID]] <- sim_plot_data
}

# Combine all simulation data into a single data frame
affectedarea <- bind_rows(all_plot_data)
affectedarea
# Save  to a CSV file
write.csv(affectedarea, "D:/Wena_Project/Testing/Test/SA/affectedarea.txt", row.names = FALSE)


#########Step 5

########### Local sensitivity analysis by calculting the variation to reference values simID==1
library(dplyr)
library(ggplot2)
library(tidyr)

# Load the severity and tested parameters files
affectedarea <- read.csv("H:/SA/affectedarea.txt")
affectedarea
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

# Calculate mean severity grouped by simID, replicate, and Time
affectedarea_mean <- affectedarea %>%
  group_by(simID) %>%
  summarise(affectedarea_mean = mean(DamagedSites*6.25, na.rm = TRUE), .groups = "drop")
head(affectedarea_mean)
summary(affectedarea_mean)
# Merge severity_mean with tested_parameters
df_merged <- affectedarea_mean %>%
  inner_join(tested_parameters, by = c("simID"))
colnames(df_merged)
write.csv(df_merged, "D:/Wena_Project/Testing/Test/SA/df_merged.csv", row.names = FALSE)

# Calculate the mean severity across replicates for each simID
df_mean <- df_merged %>%
  group_by(simID,TemporalType, MinROS, MaxROS, Dispersal, DispersalRate, EpidemicThresh, InitialEpicenterNum,
           OutbreakEpicenterCoeff, OutbreakEpicenterThresh, SeedEpicenter, SeedEpicenterMax, 
           SeedEpicenterCoeff, DispersalTemplate, NeighborFlag, NeighborSpeedUp, NeighborRadius, 
           NeighborShape, NeighborWeight, IntensityClass2_BDP, IntensityClass3_BDP, 
           Age1, Age2, Age3, SRDProb1, SRDProb2, SRDProb3, VulnProb1, VulnProb2, VulnProb3) %>%
  summarise(affectedarea= mean(affectedarea_mean, na.rm = TRUE), .groups = "drop")  # Calculate mean severity
df_mean
colnames(df_mean)
# Export the mean severity data
write.csv(df_mean, "df_mean.csv", row.names = FALSE)

# Filter the reference severity for simID = 1
reference_affectedarea<- df_mean %>%
  filter(simID == 1) %>%
  dplyr::select(-simID)

# Calculate the variation in severity relative to the reference
variation <- df_mean %>%
  filter(simID != 1) %>%
  rowwise() %>%
  mutate(variation = 100 * (affectedarea - reference_affectedarea$affectedarea[1]) / reference_affectedarea$affectedarea[1]) %>%
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
reference_long <- reference_affectedarea %>%
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

varied_parameters<-varied_parameters%>% filter(Parameter!="affectedarea")
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
      y = "Affected Area Variation (%)",
      fill = "Parameters"
    ) +
    theme_minimal() +
    scale_fill_viridis_d(option = "cividis") +
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
output_directory <- "D:/NewVersion/Local_sensitivity_analysis/plots_LSA_Area/"
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
heatmap_data <- varied_parameters %>%
  dplyr::select(-c(simID)) %>%  # Exclude simID from heatmap analysis
  mutate(variation = as.numeric(variation)) %>%
  group_by(ParameterGroup, Parameter) %>%
  summarise(
    MeanVariation = mean(variation, na.rm = TRUE),
    .groups = "drop"
  )

# Create a heat map
heatmap_plot <- ggplot(heatmap_data, aes(x = Parameter, y = ParameterGroup, fill = MeanVariation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-100, 100), name = "Variation (%)"
  ) +
  labs(
    title = "Heat Map of Parameter Sensitivity",
    x = "Parameter",
    y = "Parameter Group"
  ) +
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
ggsave("D:/NewVersion/Local_sensitivity_analysis/plots_LSA_Area/All_parameter_sensitivity.png", plot = heatmap_plot, width = 12, height = 8)
