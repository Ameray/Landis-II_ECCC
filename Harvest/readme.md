# 🌲 LANDIS-II Biomass Harvest Prescriptions

- This repository provides a full implementation of **Biomass Harvest** prescriptions for use with the LANDIS-II model. The strategy includes both **even-aged** and **uneven-aged silvicultural systems**, targeting species based on **shade tolerance**, **maturity age**, **longevity**, and **ecological role** in boreal and mixedwood forests.
- To run `harvest_calibration_simPkgs.R`, copy the initial file named `/harvest/biomass-harvest.txt` into your working directory (`dir_path`). This file contains the harvest prescription details. Then run simulation using `script/Landscape_SimPilot.R`.
- `mgmt-areas_temperate-2a-3b.tif` defines the new managed area excluding protected areas (Mapcode 0), and includes management units (Mapcodes 1 to 5) as well as private forest (Mapcode 6).

---

## 🎯 Objectives

These prescriptions aim to:

- Simulate realistic forest management aligned with ecological forestry principles.
- Reflect current provincial practices such as CPRS, seed-tree retention, commercial thinning, and shelterwood systems.
- Target species and stand ages based on silvics and economic importance.

---

## 📋 Prescription Summary

| Treatment ID | Name                                | System Type                              | Target Stand Age | Harvest Intensity         | Frequency       |
|--------------|-------------------------------------|------------------------------------------|------------------|---------------------------|-----------------|
| CPRS         | Clearcut with protection            | Even-aged                                | ≥ 40 yrs         | 95% of mature cohorts     | One-time        |
| SeedTree     | Seed tree retention                 | Even-aged                                | ≥ 40 yrs         | 100% then 0% of old forest| One-time        |
| CommThin     | Commercial thinning                 | Thinning                                 | ≥ 40 yrs         | 25% of mature cohorts     | One-time        |
| PC1          | Regular shelterwood                 | Partial cut (even-aged regen)            | 40–100 yrs       | 70%                       | One-time        |
| PC2          | Irregular shelterwood (slow regen)  | Partial cut (tolerant regen)             | 41–100 yrs       | 50%                       | Repeat (60 yrs) |
| PC3          | Shelterwood with permanent canopy   | Multi-aged retention                     | ≥ 70 yrs         | 50% of older cohorts      | Repeat (60 yrs) |

---

## 🌿 Species and Age Logic

Each treatment targets species based on their **shade tolerance**, **regeneration strategy**, and **economic importance**:

| Species Code     | Common Name         | Shade Tolerance  | Maturity (yrs)  | Longevity (yrs)  | Assigned Prescriptions                 |
|------------------|---------------------|------------------|------------------|------------------|----------------------------------------|
| ABIE.BAL         | Balsam fir          | Very tolerant    | 30               | 150              | CPRS, SeedTree, CommThin, PC2, PC3     |
| ACER.RUB         | Red maple           | Intermediate     | 10               | 150              | CPRS, SeedTree, PC1, PC3               |
| ACER.SAH         | Sugar maple         | Very tolerant    | 40               | 300              | PC1, PC2, PC3                          |
| BETU.ALL         | Yellow birch        | Mod. tolerant    | 40               | 300              | CPRS, SeedTree, CommThin, PC1          |
| BETU.PAP         | White birch         | Intermediate     | 20               | 150              | CPRS, SeedTree, CommThin, PC1          |
| FAGU.GRA         | Beech               | Very tolerant    | 40               | 250              | PC1, PC2, PC3                          |
| LARI.LAR         | Larch               | Intolerant       | 40               | 150              | CPRS, SeedTree                         |
| PICE.GLA         | White spruce        | Tolerant         | 30               | 200              | PC1, PC2                               |
| PICE.MAR         | Black spruce        | Intermediate     | 30               | 200              | CPRS, SeedTree, CommThin               |
| PICE.RUB         | Red spruce          | Mod. tolerant    | 30               | 300              | CPRS, SeedTree, PC1, PC2, PC3          |
| PINU.BAN         | Jack pine           | Intolerant       | 20               | 150              | CPRS, SeedTree, CommThin               |
| PINU.RES         | Red pine            | Intolerant       | 40               | 200              | CPRS, SeedTree, CommThin, PC2          |
| PINU.STR         | White pine          | Intermediate     | 20               | 300              | CPRS, SeedTree, CommThin, PC2          |
| POPU.TRE         | Trembling aspen     | Intolerant       | 20               | 150              | CPRS, SeedTree, PC1                    |
| QUER.RUB         | Red oak             | Mod. tolerant    | 30               | 250              | PC1, PC2, PC3                          |
| THUJ.SPP.ALL     | Eastern cedar       | Very tolerant    | 30               | 300              | CPRS, SeedTree, PC2, PC3               |
| TSUG.CAN         | Hemlock             | Very tolerant    | 60               | 300              | CPRS, SeedTree, PC1, PC2, PC3          |

---

## 🔧 Harvest Rules by Treatment

### 🌲 **CPRS – Clearcut with Regeneration Protection**
- **Goal:** Remove most merchantable biomass of shade-intolerant and intermediate species.
- **Cohorts Removed:**  
  - < Maturity: 0%  
  - Maturity to Longevity: 95%

---

### 🌱 **SeedTree – Seed Tree Retention**
- **Goal:** Retain seed trees to facilitate natural regeneration.
- **Cohorts Removed:**  
  - All except oldest: 100%  
  - Oldest: 0%

---

### 🌳 **Commercial Thinning**
- **Goal:** Reduce density and promote residual tree growth.
- **Cohorts Removed:**  
  - < Maturity: 0%  
  - ≥ Maturity: 25%

---

### 🌲 **PC1 – Regular Shelterwood**
- **Goal:** Prepare for even-aged regeneration with moderate canopy removal.
- **Cohorts Removed:**  
  - 40–100 yrs: 70%  
  - Others: 0%

---

### 🌲 **PC2 – Irregular Shelterwood (Slow Regeneration)**
- **Goal:** Regenerate tolerant species through repeated light cuts.
- **Cohorts Removed:**  
  - 41–100 yrs: 50%  
  - Others: 0%

---

### 🌲 **PC3 – Permanent Canopy Shelterwood**
- **Goal:** Maintain canopy while removing older, declining cohorts.
- **Cohorts Removed:**  
  - ≥ 70 yrs: 50%  
  - < 70 yrs: 0%

---

## 🗺️ Harvest Maps

Output files generated:

- `prescripts-{timestep}.tif`: Spatial raster of treatment codes  
- `biomass-removed-{timestep}.tif`: Raster of removed biomass  
- `harvest/log.csv`: Cohort-level harvest log  
- `harvest/summarylog.csv`: Stand-level harvest summary  

---

## 📈 Harvest Rates by Management Area

| Management Area | Treatment            | Rate (%)   |
|------------------|----------------------|------------|
| 1                | CPRS                 | 0.004788   |
| 1                | PC1                  | 0.001050   |
| 1                | PC2                  | 0.017556   |
| 1                | PC3                  | 0.018606   |
| 2                | CPRS                 | 0.044880   |
| 2                | CommercialThinning   | 0.002860   |
| 2                | PC1                  | 0.003410   |
| 2                | PC2                  | 0.032230   |
| 2                | PC3                  | 0.026510   |
| 3                | CPRS                 | 0.003570   |
| 3                | SeedTree             | 0.001820   |
| 3                | CommercialThinning   | 0.000175   |
| 3                | PC1                  | 0.000175   |
| 3                | PC2                  | 0.008435   |
| 3                | PC3                  | 0.015680   |
| 4                | CPRS                 | 0.005760   |
| 4                | PC2                  | 0.034624   |
| 4                | PC3                  | 0.023616   |
| 5                | CPRS                 | 0.099944   |
| 5                | CommercialThinning   | 0.011408   |
| 5                | PC1                  | 0.047616   |
| 5                | PC2                  | 0.017608   |
| 5                | PC3                  | 0.071424   |
| 6                | CPRS                 | 0.159500   |
| 6                | SeedTree             | 0.002000   |
| 6                | CommercialThinning   | 0.014500   |
| 6                | PC1                  | 0.052500   |
| 6                | PC2                  | 0.112000   |
| 6                | PC3                  | 0.159500   |

--- 
