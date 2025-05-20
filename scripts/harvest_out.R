### Ameray Abderrahmane 
#### organisation :ECCC
##### 06-05-2025
###  this script allows to compare the total ACC simulated and empirical  

# Clear the environment 
rm(list = ls())
#########

library(tidyverse)
library(dplyr)
#Set working directory 
setwd("D:/JaneProject/2025-05-16")

simInfo <- read_csv("simInfo.csv", show_col_types = FALSE)
simIDs   <- simInfo$simID
simIDs
df_all <- map_dfr(simIDs, function(sim) {
  summary_path <- file.path(sim, "harvest", "summarylog.csv")
  if (!file.exists(summary_path)) {
    warning("Missing summarylog.csv for simID = ", sim)
    return(NULL)
  }
  read_csv(summary_path, show_col_types = FALSE) %>%
    group_by(Time) %>%
    summarize(
      total_Mg = sum(TotalBiomassHarvested, na.rm = TRUE) / 1e6,   ### convert to Mt
      .groups = "drop"
    ) %>%
    mutate(simID = sim)
})
df_all<-df_all %>% filter()

# Compute meanAAC per simID
meanAAC_simID <- df_all %>%
  group_by(simID) %>%
  summarise(meanAAC = mean(total_Mg), .groups = "drop")

# Plot with  emprical and simulated meanAAC lines
p <- ggplot(df_all, aes(x = Time, y = total_Mg)) +
  geom_line(color = "steelblue", size = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +  # Empirical AAC
  geom_hline(data = meanAAC_simID, aes(yintercept = meanAAC), 
             color = "black", linetype = "dotted", linewidth = 1) +   # Simulated meanAAC
  facet_grid(~ simID) +
  labs(
    title    = "Total Biomass Harvested (Mt) Over Time per simID  in temperate-2a-3b ",
    subtitle = "Red dashed: Empirical AAC | Black dotted: Simulated Mean AAC",
    x        = "Time",
    y        = "Total Biomass (Mt)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

# Display the plot
print(p)
# Save the plot as a JPG file
ggsave(
  filename = "D:/JaneProject/biomass_harvest_plot.jpg",  # File name
  plot = p,                               # The plot object
  width = 8,                             # Width in inches
  height = 3,                             # Height in inches
  dpi = 300                               # Resolution (dots per inch)
)




# -----------------------------------------------------------------------------
#  species‐level harvest across all sims
# -----------------------------------------------------------------------------
df_sp_all <- map_dfr(simIDs, function(sim) {
  path <- file.path(sim, "harvest", "summarylog.csv")
  if (!file.exists(path)) {
    warning("Missing summarylog.csv for simID = ", sim)
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE) %>%
    dplyr::select(Time, Prescription, starts_with("BiomassHarvestedMg_")) %>%
    pivot_longer(
      cols      = starts_with("BiomassHarvestedMg_"),
      names_to  = "species",
      values_to = "harvested_Mg"
    ) %>%
    mutate(
      species = str_remove(species, "^BiomassHarvestedMg_"),
      simID   = sim
    )
}) %>%
  group_by(simID, Prescription, Time, species) %>%
  summarise(
    harvested_Mg = sum(harvested_Mg, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# species‐level time series, faceted by simID × Prescription
# -----------------------------------------------------------------------------
p1 <- ggplot(df_sp_all, aes(x = Time, y = harvested_Mg, color = species)) +
  geom_line(size = 1) +
  facet_grid(simID ~ Prescription) +
  labs(
    title = "Harvest by Species over Time",
    x     = "Time",
    y     = "Harvested Biomass (Mg)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))
p1
# -----------------------------------------------------------------------------
#  species‐level area plot
# -----------------------------------------------------------------------------
p2 <- ggplot(df_sp_all, aes(x = Time, y = harvested_Mg/10e6, fill =species   )) +
  geom_area(alpha = 0.7, position = "stack") +
  facet_grid(simID ~ Prescription,scales="free_y") +
  labs(
    title = "Harvest by Species (Stacked Area)",
    x     = "Time",
    y     = "Harvested Biomass (Mg)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))
p2
