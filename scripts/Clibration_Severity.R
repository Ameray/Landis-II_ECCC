################################################################################
#### organisation ECCC
### Dominic Cyr/ Ameray Abderrahmane 


######## severity calculation and validation 
#### read the final combined data and calculate severity for each simulation
# Clear workspace and load necessary libraries
rm(list = ls())
library(tidyverse)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

### ---------------------------------------------

# Set working directory and load data
dir_path <- "D:/Wena_Project/Testing/2025-04-08/"
setwd(dir_path)

fluxesBDA <- get(load("output_fluxesBDA_2025-04-08.RData"))
unique(fluxesBDA$simID)
head(fluxesBDA)
# ---------------------------------------------
# CALCULATE SEVERITY PER SPECIES PER CELL
# ---------------------------------------------
df_cell <- fluxesBDA %>%
  group_by(simID, scenario, mgmtScenario, species, Time, row, column, ecoregion,replicate , Dist,BDApara ) %>%
  summarize(
    AGB_woody_host = sum(AGB_woody, na.rm = TRUE),
    AGBtoDOM_woody_all = sum(AGBtoDOM_woody, na.rm = TRUE),
    severity = (AGBtoDOM_woody_all / AGB_woody_host) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    budworm = case_when(Time %in% 1:7 ~ "budworm1")
  )
df_cell
# ---------------------------------------------
# CALCULATE TOTAL AGB, totalAGBtoDOM AND SEVERITY2 PER CELL
# ---------------------------------------------
total_AGB_per_cell <- df_cell %>%
  group_by(simID, scenario, mgmtScenario, Time, row, column, ecoregion,BDApara ,replicate, Dist) %>%
  summarize(
    total_AGB = sum(AGB_woody_host, na.rm = TRUE),
    totalAGBtoDOM = sum(AGBtoDOM_woody_all, na.rm = TRUE),
    severity2 = 100*totalAGBtoDOM / (total_AGB),
    .groups = "drop"
  ) 

# Merge with df_cell
mean_severity <- df_cell %>%
  left_join(total_AGB_per_cell, 
            by = c("simID", "scenario", "mgmtScenario", "Time", "row", "column", "ecoregion","BDApara" ,"replicate", "Dist"))
mean_severity

# mean_severity---------------------------------------------
# AGGREGATE SIMULATED DATA PER CELL
# ---------------------------------------------
simulated_data <- mean_severity %>%
  filter(budworm == "budworm1") %>%
  group_by(scenario, mgmtScenario, Time, row, column, ecoregion,budworm ,BDApara) %>%
  summarize(severity = mean(severity2), .groups = "drop")
simulated_data
# Calculate mean severity and its confidence interval
#######
# Summary for Quantile Plots
dfSummary_by_sp <- simulated_data %>%
  group_by(BDApara , Time) %>%
  summarize(
    severityP01 = quantile(severity, 0.01, na.rm = TRUE),
    severityP10 = quantile(severity, 0.1, na.rm = TRUE),
    severityP25 = quantile(severity, 0.25, na.rm = TRUE),
    severityP50 = quantile(severity, 0.50, na.rm = TRUE),
    severityP75 = quantile(severity, 0.75, na.rm = TRUE),
    severityP90 = quantile(severity, 0.9, na.rm = TRUE),
    severityP99 = quantile(severity, 0.99, na.rm = TRUE),
    severityMean = mean(severity, na.rm = TRUE),
    .groups = "drop"
  )
dfSummary_by_sp
# Quantile plots with legends for each quantile line
quantile_gg <- ggplot(dfSummary_by_sp, aes(x = Time)) +
  geom_line(aes(y = severityP01, color = "P01"), size = 1) +
  geom_line(aes(y = severityP10, color = "P10"), size = 1) +
  geom_line(aes(y = severityP25, color = "P25"), size = 1) +
  geom_line(aes(y = severityP50, color = "P50"), size = 1) +
  geom_line(aes(y = severityP75, color = "P75"), size = 1) +
  geom_line(aes(y = severityP90, color = "P90"), size = 1) +
  geom_line(aes(y = severityP99, color = "P99"), size = 1) +
  scale_color_manual(values = c(
    P01 = "blue",
    P10 = "cyan",
    P25 = "green",
    P50 = "yellow",
    P75 = "orange",
    P90 = "red",
    P99 = "purple"
  ), name = "Quantile", labels = c("1st", "10th", "25th", "50th", "75th", "90th", "99th")) +
  theme_minimal() +
  labs(
    title = "Quantiles of Severity Over Time by Species",
    x = "Time",
    y = "Severity Quantiles",
    color = "Quantile"
  ) +
  facet_wrap(~ BDApara , scales = "free_y")

# Print the quantile plot
print(quantile_gg)

ggplot(dfSummary_by_sp, aes(x = Time)) +
  # Interquartile Range Ribbon
  geom_ribbon(aes(ymin = severityP25, ymax = severityP75), fill = "skyblue", alpha = 0.4) +
  
  # Optional: 10–90% range ribbon
  geom_ribbon(aes(ymin = severityP10, ymax = severityP90), fill = "lightgray", alpha = 0.2) +
  
  # Median line (P50)
  geom_line(aes(y = severityP50), color = "blue", size = 1.2, linetype = "solid") +
  
  # Mean line
  geom_line(aes(y = severityMean), color = "red", size = 1.2, linetype = "dashed") +
  
  labs(
    title = "Severity Over Time: Median vs. Mean",
    x = "Time",
    y = "Severity",
    caption = "Shaded areas show variability (IQR and 10–90%)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ BDApara )


####### compare over mean using cross validation methode 
# ---------------------------------------------
# LOAD EMPIRICAL DATA
# ---------------------------------------------
empirical_data <- read_csv("D:/Wena_Project/Testing/Test/Empirical_severity.csv")
empirical_data$severity <- as.numeric(empirical_data$severity)
empirical_data
# ---------------------------------------------
# CROSS-VALIDATION FUNCTION
# ---------------------------------------------
cross_validate_severity <- function(empirical_data, simulated_data, n_subsamples = 1000, n_iterations = 1000) {
  empirical_means <- replicate(n_iterations, mean(sample(empirical_data$severity, n_subsamples, replace = TRUE)))
  simulated_means <- replicate(n_iterations, mean(sample(simulated_data$severity, n_subsamples, replace = TRUE)))
  
  list(
    empirical_means = empirical_means,
    simulated_means = simulated_means,
    empirical_ci = quantile(empirical_means, probs = c(0.025, 0.975)),
    simulated_ci = quantile(simulated_means, probs = c(0.025, 0.975))
  )
}

# ---------------------------------------------
# VALIDATION LOOP AND PLOTTING DATA PREP
# ---------------------------------------------
validation_results <- list()
combined_plot_data <- list()
density_plot_data <- list()

for (budworm in unique(simulated_data$budworm)) {
  for (id in unique(simulated_data$BDApara )) {
    sim_data <- simulated_data %>% filter(budworm == budworm, BDApara  == id)
    results <- cross_validate_severity(empirical_data, sim_data)
    validation_results[[paste(budworm, id)]] <- results
    
    # Boxplot data
    combined_plot_data[[paste(budworm, id)]] <- tibble(
      DataCategory = c(rep("Empirical", length(results$empirical_means)),
                       rep(paste("Simulated", budworm), length(results$simulated_means))),
      MeanSeverity = c(results$empirical_means, results$simulated_means),
      CI_Lower = c(rep(results$empirical_ci[1], length(results$empirical_means)),
                   rep(results$simulated_ci[1], length(results$simulated_means))),
      CI_Upper = c(rep(results$empirical_ci[2], length(results$empirical_means)),
                   rep(results$simulated_ci[2], length(results$simulated_means))),
      BDApara  = rep(id, length(results$empirical_means) + length(results$simulated_means))
    )
    
    # Density plot data
    density_plot_data[[paste(budworm, id)]] <- tibble(
      MeanSeverity = c(results$empirical_means, results$simulated_means),
      Dataset = rep(c("Empirical", paste("Simulated", budworm)), each = length(results$empirical_means)),
      BDApara  = rep(id, 2 * length(results$empirical_means))
    )
    
    # Print CI info
    cat("Budworm Period:", budworm, "Simulation ID:", id, "\n")
    cat("Empirical severity 95% CI:", paste(round(results$empirical_ci, 2), collapse = " - "), "\n")
    cat("Simulated severity 95% CI:", paste(round(results$simulated_ci, 2), collapse = " - "), "\n\n")
  }
}

# ---------------------------------------------
# BIND DATA FOR PLOTS
# ---------------------------------------------
combined_plot_data <- bind_rows(combined_plot_data)
density_plot_data <- bind_rows(density_plot_data)

# ---------------------------------------------
# PLOT: BOXPLOT BY SIM ID
# ---------------------------------------------
ggplot(combined_plot_data, aes(x = DataCategory, y = MeanSeverity, fill = DataCategory)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  labs(
    title = "Comparison of Empirical and Simulated Severity Data by Simulation ID",
     y = "Mean Severity"
  ) +
  theme_minimal() +
  facet_wrap(~ BDApara , scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---------------------------------------------
# PLOT: DENSITY PLOT
# ---------------------------------------------
ggplot(density_plot_data, aes(x = MeanSeverity, fill = Dataset)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of Severity by Dataset and Simulation ID",
       x = "Mean Severity", y = "Density") +
  scale_fill_manual(values = c("Empirical" = "red", "Simulated budworm1" = "chartreuse")) +
  theme_minimal() +
  xlim(25, 55) +
  facet_wrap(~ BDApara , scales = "free") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "top"
  )



#### compared severity with emeprical ones only by budwroms  average accross all simulation 
#### when you find the appropriate BDApara , you can replicate the simulation and calculate the average, but deactivating RandomNumberSeed
#########################################################################
#########################################################################
##################### Validation method: Sub-sampling for cross-validation
# Read and prepare empirical data
empirical_data <- read.csv("D:/Wena_Project/Testing/Test/Empirical_severity.csv")
empirical_data$severity <- as.numeric(empirical_data$severity)
summary(empirical_data)
mean_severity
# 
simulated_data <- mean_severity %>%
  filter(budworm == "budworm1",BDApara==2) %>%
  group_by(scenario, mgmtScenario, Time, row, column, ecoregion,budworm ,BDApara) %>%
  summarize(severity = mean(severity2), .groups = "drop")
simulated_data
# Function for cross-validation with fixed subsamples
cross_validate_severity <- function(empirical_data, simulated_data, n_subsamples = 1000, n_iterations = 1000) {
  empirical_means <- replicate(n_iterations, mean(sample(empirical_data$severity, n_subsamples, replace = TRUE)))
  simulated_means <- replicate(n_iterations, mean(sample(simulated_data$severity, n_subsamples, replace = TRUE)))
  
  empirical_ci <- quantile(empirical_means, probs = c(0.025, 0.975))
  simulated_ci <- quantile(simulated_means, probs = c(0.025, 0.975))
  
  list(empirical_means = empirical_means, simulated_means = simulated_means, empirical_ci = empirical_ci, simulated_ci = simulated_ci)
}

# Perform cross-validation for each budworm period
validation_results <- list()
for (budworm in unique(simulated_data$budworm)) {
  sim_data <- simulated_data %>% filter(budworm == budworm)
  results <- cross_validate_severity(empirical_data, sim_data)
  validation_results[[budworm]] <- results
  
  cat("Budworm Period:", budworm, "\n")
  cat("Empirical severity 95% CI:", paste(results$empirical_ci, collapse = ", "), "\n")
  cat("Simulated severity 95% CI:", paste(results$simulated_ci, collapse = ", "), "\n\n")
}

# Prepare data for plotting (boxplot and density plot)
plot_data <- do.call(rbind, lapply(names(validation_results), function(budworm) {
  results <- validation_results[[budworm]]
  data.frame(
    DataCategory = rep(c("Empirical", paste("Simulated", budworm)), each = 1000),
    MeanSeverity = c(results$empirical_means, results$simulated_means),
    CI_Lower = c(rep(results$empirical_ci[1], 1000), rep(results$simulated_ci[1], 1000)),
    CI_Upper = c(rep(results$empirical_ci[2], 1000), rep(results$simulated_ci[2], 1000))
  )
}))
###prepare data for the density plot
density_plot_data <- do.call(rbind, lapply(names(validation_results), function(budworm) {
  results <- validation_results[[budworm]]
  data.frame(
    MeanSeverity = c(results$empirical_means, results$simulated_means),
    Dataset = rep(c("Empirical", paste("Simulated", budworm)), each = length(results$empirical_means))
  )
}))


head(plot_data)
head(density_plot_data)
library(patchwork)

# Rename 'Simulated budworm1' to 'simulated' in both datasets
plot_data$FillGroup <- gsub("Simulated budworm1", "simulated", plot_data$DataCategory)
density_plot_data$FillGroup <- gsub("Simulated budworm1", "simulated", density_plot_data$Dataset)

# Define consistent color mapping with new label
fill_colors <- c("Empirical" = "red", "simulated" = "chartreuse")

# Updated boxplot
plot1_best_BDApara  <- ggplot(plot_data, aes(x = FillGroup, y = MeanSeverity, fill = FillGroup)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  labs(x = NULL, y = "Severity") +
  scale_fill_manual(values = fill_colors, name = "Dataset") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )
plot1_best_BDApara 
# Updated density plot
plot2_best_BDApara  <- ggplot(density_plot_data, aes(x = MeanSeverity, fill = FillGroup)) +
  geom_density(alpha = 0.5) +
  labs(
    x = "Mean Severity",
    y = "Density"
  ) +
  scale_fill_manual(values = fill_colors, name = "Dataset") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )
plot2_best_BDApara
# Combine with shared legend
combined_plot <- plot1_best_BDApara  + plot2_best_BDApara  + plot_layout(guides = "collect")
combined_plot




ggsave("D:/dom_plots/combined_plot.png", plot = combined_plot, width = 8, height = 3)






