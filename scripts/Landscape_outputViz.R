###################################################################
###################################################################
################################################################################
## organisation ECCC
### Dominic Cyr
### ForCS output processing  : updated with corrected fluxes
###########################################################
###################################################################
### ForCS output processing 
rm(list = ls())

# Set working directory
wwd <- "D:/Wena_Project/Testing/Landscape/2025-04-14"
setwd(wwd)
wwd <- getwd()
wwd
#####
require(dplyr)

initYear <- 2020
unitConvFact <- 0.01 ### from gC /m2 to tonnes per ha
simName <- "2025-04-14"
a <- ifelse(grepl("mixedwood-042-51", simName), "mixedwood-042-51",
            ifelse(grepl("boreal-085-51", simName), "boreal-085-51",
                   ifelse(grepl("temperate-2a-3b", simName), "temperate-2a-3b")))

areaName <- a
# areaName <- ifelse(a == "ForMont", "For?t Montmorency",
#                    ifelse(a == "Hereford", "For?t Hereford", a))

require(ggplot2)
require(dplyr)
require(tidyr)
require(scales)

##########################################################
##### pools
variableLvl <- c("TotalEcosys", "TotalDOM", "ABio", "BBio") ## ordering levels for plotting


treatLevels <- list("boreal-085-51" = c("generic" = "Generic",
                                        #"CPI-CP" = "Couvert permanent",
                                        "noHarvest" = "Conservation"),
                    "mixedwood-042-51" = c("No_ND" = "No natural disturbances",
                                           "Wind" = "Wind",
                                           "Wind_Sbw" = "Wind and Spruce Budworm",
                                           "Wind_Sbw_Fire" = "Wind, Spruce Budworm and wildfires"),
                    #c("generic" = "Generic",
                    #"CPI-CP" = "Couvert permanent",
                    #"noHarvest" = "Conservation"),
                    "temperate-2a-3b" = c("generic" = "Generic",
                                          #"CPI-CP" = "Couvert permanent",
                                          "noHarvest" = "Conservation"))

mgmtLevels <- list("boreal-085-51" = c("generic" = "Generic",
                                       #"CPI-CP" = "Couvert permanent",
                                       "noHarvest" = "Conservation"),
                   "mixedwood-042-51" = c("generic" = "Generic",
                                          #"CPI-CP" = "Couvert permanent",
                                          "noHarvest" = "Conservation"),
                   #c("generic" = "Generic",
                   #"CPI-CP" = "Couvert permanent",
                   #"noHarvest" = "Conservation"),
                   "temperate-2a-3b" = c("generic" = "Generic",
                                         #"CPI-CP" = "Couvert permanent",
                                         "noHarvest" = "Conservation"))

mgmtLvls <- c("CPRS", "Couvert permanent", "Conservation")




cols <- list("boreal-085-51" = c("Generic" = "black",
                                 #"Couvert permanent" = "dodgerblue3",
                                 "Conservation" = "goldenrod3"),
             "mixedwood-042-51" = c("No_ND" = "black",
                                    "Wind" = "lightgreen",#"Couvert permanent" = "dodgerblue3",
                                    "Wind_Sbw" = "goldenrod3",
                                    "Wind_Sbw_Fire" = "indianred"),
             "temperate-2a-3b" = c("Generic" = "black",
                                   #"Couvert permanent" = "dodgerblue3",
                                   "Conservation" = "goldenrod3"))

# scenRef <-  list("boreal-085-51" = "generic",
#                  "mixedwood-042-51" = "Wind_Sbw_Fire",
#                  "temperate-2a-3b" = "generic")
################################################################################
################################################################################
################################################################################
outputSummary <- get(load(paste0(wwd,"/output_summary_", simName, ".RData")))
fps <- read.csv(paste0(wwd,"/output_BioToFPS_", simName, ".csv"))#read.csv("output_BioToFPS_Hereford.csv")#
AGB <- get(load(paste0(wwd,"/output_bio_", simName, ".RData")))

################################################################################
### renaming stuff
outputSummary <- outputSummary %>%
  filter(variable != "mgmtScenarioName") %>%
  mutate(value = as.numeric(value)) %>%
  mutate(ND_scenarioName  = factor(treatLevels[[a]][match(as.character(ND_scenario), names(treatLevels[[a]]))],
                                   levels = treatLevels[[a]]),
         mgmtScenarioName  = factor(mgmtLevels[[a]][match(as.character(mgmtScenarioName), names(mgmtLevels[[a]]))],
                                    levels = mgmtLevels[[a]]))

outputSummary <- droplevels(outputSummary)

AGB <- AGB %>%
  mutate(ND_scenarioName  = factor(treatLevels[[a]][match(as.character(ND_scenario), names(treatLevels[[a]]))],
                                   levels = treatLevels[[a]]),
         mgmtScenarioName  = factor(mgmtLevels[[a]][match(as.character(mgmtScenarioName), names(mgmtLevels[[a]]))],
                                    levels = mgmtLevels[[a]]))

fps <- fps %>%
  mutate(ND_scenarioName  = factor(treatLevels[[a]][match(as.character(ND_scenario), names(treatLevels[[a]]))],
                                   levels = treatLevels[[a]]),
         mgmtScenarioName  = factor(mgmtLevels[[a]][match(as.character(mgmtScenarioName), names(mgmtLevels[[a]]))],
                                    levels = mgmtLevels[[a]]))



require(ggplot2)
require(dplyr)
require(tidyr)


### pools
df <- outputSummary %>%
  filter(Time >=1,variable %in% variableLvl) %>%
  group_by(areaName, scenario, 
           ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName, 
           ##################
           ##plantedSp,
           
           # mgmtID,
           Time, variable) %>%
  summarise(value = mean(value),
            mgmtArea_ha = unique(mgmtArea_ha)) %>%
  group_by(areaName, scenario, 
           ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName, 
           ##################
           ##plantedSp,
           
           # mgmtID,
           Time, variable) %>%
  summarise(valueTotal = sum(value*mgmtArea_ha),
            mgmtArea_ha = sum(mgmtArea_ha)) %>%
  mutate(value = valueTotal/mgmtArea_ha) %>%
  as.data.frame()


df <- df %>%
  group_by(areaName, scenario,
           ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName, 
           ##################
           #plantedSp, mgmt,
           Time) %>%
  summarise(valueTotal = sum(valueTotal),
            mgmtArea_ha = unique(mgmtArea_ha)) %>%
  mutate(value = valueTotal / mgmtArea_ha,
         variable = "TotalEcosys") %>%
  as.data.frame()%>%
  rbind(df) %>%
  mutate(variable = factor(variable, levels = variableLvl)) 



# x <- df %>%
#     filter(areaName == a,
#            ND_scenarioName == scenRef[[a]])
# cNames <- colnames(df)
# cNames <- cNames[-which(cNames %in% c("treatment", "mgmtScenario", "mgmtName", "ndName", "valueTotal", "value", "mgmtArea_ha"))]

# dfRef <- merge(df, x,
#              by = cNames,
#              suffixes = c("",".ref"))
# 
# cNames <- c(colnames(df), "value.ref")
# dfRef <- dfRef[,cNames]
# 
# ref <- df


p <- ggplot(df, aes(x = initYear+Time, y = value*unitConvFact,
                    colour = ND_scenario, linetype = mgmtScenarioName  )) +
  facet_grid(variable ~ scenario, scale = "free") +
  theme_bw() +
  scale_color_manual(name = "Natural disturbance\nscenario",
                     values = cols[[a]],
                     labels = treatLevels[[a]]) +
  scale_linetype_manual(name = "Management\nscenario",
                        values = c("Generic" = "solid",
                                   "Conservation" = "dotted")) +
  
  geom_line(#linetype = mgmtName,
    linewidth = 0.5) +
  theme(plot.caption = element_text(size = rel(.5), hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Evolution of carbon density",
       subtitle = paste(areaName, simName),
       x = "",
       y = expression(paste("tonnes C"," ha"^"-1","\n")),
       caption = paste0("Allometric equations from: Ung et al. 2008. Canadian national biomass equations: new parameter estimates that include British Columbia data. Can. J. For. Res 38:1123-2232.\n",
                        "Taper equations from: Honer et al. 1983. Metric timber tables for the commercial tree species of Central and Eastern Canada."))



png(filename= paste0("pools_", simName, ".png"),
    width = 9, height = 6, units = "in", res = 600, pointsize=10)

print(p)

dev.off()



### fluxes


variableLvl <- c("DelBio",  "Turnover",
                 "NetGrowthc",  "NPPc",
                 "Rh",  "NEPc",
                 "NBPc")

caption <- c("DelBio: Annual change in biomass stocks",
             "Turnover: Annual transfer of biomass (above-and belowground) to dead organic matter and soil pools before disturbances occur",
             "NetGrowthc: corrected Change in biomass from growth alone: the difference between the biomass at the beginning and the end of the growth routine in the timestep. This value could be negative\n         as the stand ages and mortality outpaces growth. DelBio and NetGrowth will be the same when there are no losses caused by disturbances.",
             "NPPc :corrected Net Primary Production (includes above and belowground). This includes growth and replacement of litterfall and annual turnover, i.e., the sum of NetGrow and turnover.",
             "Rh : Heterotrophic respiration. This is the sum of the 'To Air' fluxes through decomposition, not disturbance.",
             "NEPc: corrected Net Ecosystem Productivity. NPP minus Rh.",
             "NBPc :corrected Net Biome Productivity. NEP minus losses from the ecosystem due to disturbances. (Both emissions to air from combustion and losses to the forest products sector.)")


df <- outputSummary %>%
  filter(Time >=1,
         mgmtID >= 10000 |
           mgmtID == 1,
         variable %in% variableLvl) %>%
  ##################
group_by(areaName, scenario,
         ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName, 
         replicate, Time, variable) %>%
  summarize(totalArea = sum(mgmtArea_ha),
            value = weighted.mean(value, mgmtArea_ha)) %>%
  ungroup() %>%
  group_by(areaName, scenario,
           ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName, 
           Time, variable) %>%
  summarise(value = mean(value))


p <- ggplot(df, aes(x = initYear+Time, y = value*unitConvFact, #group = simID,
                    colour = ND_scenario, linetype = mgmtScenarioName)) +
  facet_grid(variable ~ scenario, scale = "free") +#facet_wrap( ~ scenario) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 1, color = "grey25", size = 0.35) +
  geom_line() +
  theme(plot.caption = element_text(size = rel(.5), hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(name = "Natural disturbance\nscenario",
                     values = cols[[a]],
                     labels = treatLevels[[a]]) +
  scale_linetype_manual(name = "Management\nscenario",
                        values = c("Generic" = "solid",
                                   "Conservation" = "dotted")) +
  theme(plot.caption = element_text(size = rel(.5), hjust = 0)) +
  labs(title = paste("Carbon dynamics"),
       subtitle = paste(areaName, simName),
       x = "",
       y = expression(paste("tonnes C"," ha"^"-1", " yr"^"-1", "\n")),
       caption = paste(caption, collapse = "\n")) 


png(filename= paste0("fluxes_", simName, ".png"),
    width = 7, height = 8, units = "in", res = 600, pointsize=10)

print(p +  theme(legend.text=element_text(size=rel(.75))))
dev.off()

################################################################################
################################################################################
### to FPS
df <- fps %>%
  group_by(areaName, scenario,
           ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName,
           species, Time) %>%
  summarise(BioToFPS_tonnesCTotal = mean(BioToFPS_tonnesCTotal),
            areaManagedTotal_ha = unique(areaManagedTotal_ha),
            areaHarvestedTotal_ha = mean(areaHarvestedTotal_ha))

dfTotal <- df %>%
  group_by(areaName, scenario,
           ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName,
           Time) %>%
  summarise(BioToFPS_tonnesCTotal = sum(BioToFPS_tonnesCTotal),
            areaManagedTotal_ha = unique(areaManagedTotal_ha),
            areaHarvestedTotal_ha = mean(areaHarvestedTotal_ha))

labdf <- df %>%
  group_by(areaName, scenario,
           ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName) %>%
  summarise(areaManagedTotal_ha = unique(areaManagedTotal_ha),
            areaHarvestedTotal_ha = mean(areaHarvestedTotal_ha)) 

yMax <- df %>%
  group_by(areaName, scenario,
           ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName,
           Time) %>%
  summarise(BioToFPS_tonnesCTotal = sum(BioToFPS_tonnesCTotal/areaHarvestedTotal_ha)) %>%
  group_by() %>%
  summarise(yMax = max(BioToFPS_tonnesCTotal))
yMax <- as.numeric(yMax)


require(RColorBrewer)
colourCount = length(unique(df$species))
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

### stacked (per species)
pWidth <- 4*length(unique(df$ndName))+2
png(filename= paste0("fps_spp_", simName, ".png"),
    width = pWidth, height = 6, units = "in", res = 600, pointsize=10)

#ggplot(df, aes(x = 2010+Time, y = BioToFPS_tonnesCTotal/areaHarvestedTotal_ha)) + 
ggplot(df, aes(x = initYear+Time, y = BioToFPS_tonnesCTotal)) + 
  stat_summary(aes(fill = species), fun.y="sum", geom="area", position = "stack") +
  facet_grid(scenario ~ ND_scenarioName  ) +
  scale_fill_manual(values = getPalette(colourCount)) +
  scale_y_continuous(labels = label_number(suffix = "kt C", scale = 1e-3)) +
  theme_dark() +
  theme(plot.caption = element_text(size = rel(.5), hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Transfers to harvested wood products",
       subtitle = paste(areaName, simName),
       x = "",
       y = expression(paste("harvested wood",))) +
  geom_text(data = labdf, aes(label = paste("Managed area:", areaManagedTotal_ha, "ha"),
                              y = yMax, x = initYear),
            hjust = 0, vjust = 1, size = 2)

dev.off()



### total

p <- ggplot(dfTotal, aes(x = initYear+Time, y = BioToFPS_tonnesCTotal,
                         colour = ND_scenario)) +
  geom_line() +
  facet_wrap( ~ scenario) +
  scale_color_manual(name = "Disturbance\nscenario",
                     values = cols[[a]],
                     labels = treatLevels[[a]]) +
  scale_y_continuous(labels = label_number(suffix = "kt C", scale = 1e-3)) +
  theme_dark() +
  theme(plot.caption = element_text(size = rel(.5), hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  #lims(y = c(0,20000)) +
  labs(title = "Transfers to harvested wood products",
       subtitle = paste(areaName, simName),
       x = "",
       y = expression(paste("harvested wood",))) 
#y = expression(paste("tonnes C"," ha"^"-1", "r?colt?","\n")))


png(filename= paste0("fps_total_", simName, ".png"),
    width = 8, height = 4, units = "in", res = 600, pointsize=10)
print(p)

dev.off()




################################################################################
################################################################################
################################################################################
### Aboveground biomass
################################################################################
require(raster)
# 
# ### must correct this
# # 
# if(a == "ForMont") {
#     r <- raster("../inputsLandis/landtypes_ForMont_cropped.tif")
# }
# if(a == "Hereford") {
#     r <- raster("../inputsLandis/landtypes_Hereford_cropped.tif")
# }
# 
# # totalArea <- AGB %>%
# #     distinct(areaName, landtype, landtypeArea_ha)
# # totalArea <- sum(totalArea$landtypeArea_ha)
# 
# # 
# # foo <- filter(AGB,
# #               scenario == "baseline",
# #               mgmtScenario == "0",
# #               replicate == 1, Time == 10, species == "ABIE.BAL")
# 
# # sum(foo$agb_tonnesTotal)/totalArea


##without the age classes
df <- AGB %>%
  
  group_by(areaName, scenario,
           ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName,
           Time, replicate,
           species) %>%
  summarise(agb_tonnesTotal = sum(agb_tonnesTotal),
            areaTotal_ha = sum(unique(landtypeArea_ha))) %>%
  group_by(areaName, scenario,
           ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName,
           Time, species) %>%
  summarise(agb_tonnesTotal = mean(agb_tonnesTotal),
            areaTotal_ha = unique(areaTotal_ha)) %>%
  as.data.frame()




require(RColorBrewer)

### stacked (total)
colourCount = length(unique(df$species))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
pHeight <- 2 * length(unique(df$ND_scenario))
png(filename= paste0("agb_sppStack_", simName, ".png"),
    width = 8, height = pHeight, units = "in", res = 600, pointsize=10)

ggplot(df, aes(x = 2020 + Time, y = agb_tonnesTotal/areaTotal_ha)) + 
  stat_summary(aes(fill = species), fun.y="sum", geom="area", position = "stack") +
  facet_grid(mgmtScenarioName ~ ND_scenarioName) +
  scale_fill_manual("",
                    values = getPalette(colourCount)) +
  theme_bw() +
  theme(plot.caption = element_text(size = rel(.5), hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Evolution of landscape composition  - Aboveground biomass*",
       subtitle = paste(areaName, simName),
       x = "",
       y = expression(paste("tonnes"," ha"^"-1")),
       caption = "*Values are expressed as total dry mass (biomass), not carbon mass.")


dev.off()


### line (total)
colourCount = length(unique(df$species))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

png(filename= paste0("agb_sppLine_", simName, ".png"),
    width = 8, height = pHeight, units = "in", res = 600, pointsize=10)

ggplot(df, aes(x = 2020+Time, y = agb_tonnesTotal/areaTotal_ha,
               colour = species,
               group = species)) + 
  geom_line() +
  #stat_summary(aes(fill = species), fun.y="sum", geom="area", position = "stack") +
  facet_grid(mgmtScenarioName ~ ND_scenarioName) +
  scale_colour_manual("",
                      values = getPalette(colourCount)) +
  #scale_fill_manual(values = getPalette(colourCount)) +
  theme_bw() +
  theme(plot.caption = element_text(size = rel(.5), hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Evolution of landscape composition  - Aboveground biomass*",
       subtitle = paste(areaName, simName),
       x = "",
       y = expression(paste("tonnes"," ha"^"-1")),
       caption = "*Values are expressed as total dry mass (biomass), not carbon mass.")


dev.off()
# 



################################################################################
################################################################################
################################################################################
### Aboveground biomass (by age classes)
################################################################################
####
cols = brewer.pal(n = 9, name = 'Greens')[3:9]
##########################




df <- AGB %>% mutate(mgmt =  mgmtScenario) %>% 
  group_by(areaName, scenario,
           ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName,
           Time, replicate,
           ageClass,
           species) %>%
  summarise(agb_tonnesTotal = sum(agb_tonnesTotal),
            areaTotal_ha = sum(unique(landtypeArea_ha))) %>%
  group_by(areaName, scenario,
           ND_scenario, ND_scenarioName, mgmtScenario, mgmtScenarioName,
           ageClass,
           Time, species) %>%
  summarise(agb_tonnesTotal = mean(agb_tonnesTotal),
            areaTotal_ha = unique(areaTotal_ha)) %>%
  as.data.frame()


#require(RColorBrewer)


################################################################################
### stacked (age classes, global)
cols = rev(brewer.pal(n = 9, name = 'Greens')[3:9])

# ################################################################################


acLvls <- levels(df$ageClass)
acLvls <- rev(acLvls)
df[,"ageClass"] <- factor(df$ageClass, levels = acLvls)


#### mgmt ~ scenario

p <- ggplot(df, aes(x = 2020+Time, y = agb_tonnesTotal/areaTotal_ha)) + 
  stat_summary(aes(fill = ageClass), fun.y="sum", geom="area", position = "stack") +
  facet_grid(mgmtScenarioName ~ ND_scenarioName, scales = "fixed") +#scales = "free_y") +
  scale_fill_manual("Age classes",
                    values = cols)+
  theme_bw() +
  theme(plot.caption = element_text(size = rel(1), hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste0("Evolution of age structure - ", areaName, " - ", simName),
       subtitle = paste0("Aboveground biomass broken down by age classes"),
       x = "",
       y = expression(paste("tonnes"," ha"^"-1")),
       caption = "*Values are expressed as total dry mass (biomass), not carbon mass.")

png(filename= paste0("agb_AgeClassStacked_",simName, ".png"),
    width = 10, height = 6, units = "in", res = 600, pointsize=10)
print(p)
dev.off()



