
##🌲 LANDIS-II Biomass Harvest Prescriptions

This repository provides a full implementation of **Biomass Harvest** prescriptions for use with the LANDIS-II model. The strategy includes both **even-aged** and **uneven-aged silvicultural systems**, targeting species based on **shade tolerance**, **maturity age**, **longevity**, and **ecological role** in boreal and mixedwood forests.

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
| CPRS         | Clearcut with protection            | Even-aged                                | ≥ 60 yrs         | 95% of mature cohorts     | One-time        |
| SeedTree     | Seed tree retention                 | Even-aged                                | ≥ 60 yrs         | 100% then 0% of old forest| One-time        |
| CommThin     | Commercial thinning                 | Thinning                                 | ≥ 40 yrs         | 25% of mature cohorts     | One-time        |
| PC1          | Regular shelterwood                 | Partial cut (even-aged regen)            | 40–100 yrs       | 70%                       | One-time        |
| PC2          | Irregular shelterwood (slow regen)  | Partial cut (tolerant regen)             | 41–100 yrs       | 50%                       | Repeat (60 yrs) |
| PC3          | Shelterwood with permanent canopy   | Multi-aged retention                     | ≥ 70 yrs         | 50% of older cohorts      | Repeat (60 yrs) |

---

## 🌿 Species and Age Logic

Each treatment targets species based on their **shade tolerance**, **regeneration strategy**, and **economic importance**:

| Species Code   | Common Name         | Tolerance | Maturity | Longevity | Primary Treatments                     |
|----------------|---------------------|-----------|----------|-----------|----------------------------------------|
| ABIE.BAL       | Balsam fir          | Very tol. | 30       | 150       | CPRS, SeedTree, PC2, PC3               |
| ACER.RUB       | Red maple           | Int.      | 10       | 150       | All treatments                         |
| ACER.SAH       | Sugar maple         | Very tol. | 40       | 300       | CommThin, PC1, PC2, PC3                |
| BETU.ALL       | Yellow birch        | Mod tol.  | 40       | 300       | CPRS, CommThin, PC1                    |
| BETU.PAP       | White birch         | Int.      | 20       | 150       | CPRS, CommThin, PC1                    |
| FAGU.GRA       | Beech               | Very tol. | 40       | 250       | CommThin, PC1, PC2, PC3                |
| LARI.LAR       | Larch               | Intol.    | 40       | 150       | CPRS, SeedTree                         |
| PICE.GLA       | White spruce        | Tol.      | 30       | 200       | CommThin, PC1, PC2                     |
| PICE.MAR       | Black spruce        | Int.      | 30       | 200       | CPRS, SeedTree                         |
| PICE.RUB       | Red spruce          | Mod tol.  | 30       | 300       | All treatments                         |
| PINU.BAN       | Jack pine           | Intol.    | 20       | 150       | CPRS, SeedTree                         |
| PINU.RES       | Red pine            | Intol.    | 40       | 200       | CPRS, SeedTree                         |
| PINU.STR       | White pine          | Int.      | 20       | 300       | CPRS, SeedTree                         |
| POPU.TRE       | Trembling aspen     | Intol.    | 20       | 150       | CPRS, CommThin, PC1                    |
| QUER.RUB       | Red oak             | Mod tol.  | 30       | 250       | CommThin, PC1                          |
| THUJ.SPP.ALL   | Eastern cedar       | Very tol. | 30       | 300       | CPRS, SeedTree, PC2, PC3               |
| TSUG.CAN       | Hemlock             | Very tol. | 60       | 300       | All treatments                         |

---

## 🔧 Harvest Rules by Treatment

This section outlines the harvest rules applied to each prescription in the Biomass Harvest extension for LANDIS-II. Species were selected based on silvicultural characteristics such as shade tolerance, regeneration strategy, and economic value.

---

### 🌲 **CPRS – Clearcut with Regeneration Protection**
- **Goal:** Remove most merchantable biomass of shade-intolerant and intermediate species.
- **Rationale:** These species require full sun for regeneration and do not tolerate competition.

#### 🪓 Cohorts Removed:
- Ages < Maturity: 0%  
- Ages Maturity–Longevity: 95%

#### 🌿 Target Species:
- **Intolerant conifers:** PINU.BAN, PICE.MAR, PINU.STR, PINU.RES  
- **Pioneer hardwoods:** BETU.PAP, POPU.TRE, ACER.RUB  
- **Other intermediates:** LARI.LAR, THUJ.SPP.ALL, TSUG.CAN, PICE.RUB, BETU.ALL, ABIE.BAL

---

### 🌱 **SeedTree – Two-Stage Shelterwood**
- **Goal:** Retain scattered seed trees to promote natural regeneration.
- **Rationale:** Seed-bearing individuals are left in place while harvesting younger cohorts.

#### 🪓 Cohorts Removed:
- All cohorts except oldest: 100%  
- Oldest (approaching longevity): 0%

#### 🌿 Target Species:
- Same as CPRS — species that naturally regenerate via seed, such as PICE.MAR, PINU.BAN, POPU.TRE

---

### ✂️ **Commercial Thinning**
- **Goal:** Reduce stand density and enhance growth of residual trees.
- **Rationale:** Especially effective for mid-tolerant and tolerant species that benefit from reduced competition.

#### 🪓 Cohorts Removed:
- Ages < Maturity: 0%  
- Ages ≥ Maturity: 25%

#### 🌿 Target Species:
- **Hardwoods:** ACER.SAH, FAGU.GRA, ACER.RUB, QUER.RUB, BETU.PAP, BETU.ALL  
- **Conifers:** TSUG.CAN, PICE.GLA, PICE.RUB, POPU.TRE

---

### 🌳 **PC1 – Regular Shelterwood**
- **Goal:** Promote even-aged regeneration of moderately shade-tolerant species.
- **Rationale:** Mid-tolerant species respond well to phased overstory removal and light availability.

#### 🪓 Cohorts Removed:
- < 40 yrs: 0%  
- 40–100 yrs: 70%  
- > 100 yrs: 0%

#### 🌿 Target Species:
- **Hardwoods:** ACER.SAH, ACER.RUB, FAGU.GRA, QUER.RUB, BETU.PAP, BETU.ALL, POPU.TRE  
- **Conifers:** PICE.RUB, PICE.GLA, TSUG.CAN

---

### 🌲 **PC2 – Irregular Shelterwood (Slow Regeneration)**
- **Goal:** Encourage slow, natural regeneration of very shade-tolerant species with moderate canopy removal.
- **Rationale:** These species regenerate well in shaded environments but require reduced competition.

#### 🪓 Cohorts Removed:
- ≤ 40 yrs: 0%  
- 41–100 yrs: 50%  
- > 100 yrs: 0%

#### 🌿 Target Species:
- **Conifers:** TSUG.CAN, THUJ.SPP.ALL, ABIE.BAL, PICE.RUB, PICE.GLA  
- **Hardwoods:** ACER.SAH, FAGU.GRA

---

### 🌲 **PC3 – Shelterwood with Permanent Canopy**
- **Goal:** Maintain continuous canopy while selectively removing senescent trees.
- **Rationale:** Mimics multi-aged forest dynamics and is ideal for very shade-tolerant species.

#### 🪓 Cohorts Removed:
- < 70 yrs: 0%  
- ≥ 70 yrs: 50%

#### 🌿 Target Species:
- **Conifers:** TSUG.CAN, THUJ.SPP.ALL, ABIE.BAL, PICE.RUB  
- **Hardwoods:** ACER.SAH, ACER.RUB, FAGU.GRA

---

## 🗺️ Harvest Maps

Output files generated by the LANDIS-II Biomass Harvest extension:

- `prescripts-{timestep}.tif`: Raster of harvest prescriptions  
- `biomass-removed-{timestep}.tif`: Biomass removal raster  
- `harvest/log.csv`: Log of harvest events  
- `harvest/summarylog.csv`: Summary of cohorts & biomass harvested  

---

## 🧪 Harvest Rates Applied

The table below summarizes initial harvest rates. These will be calibrated to match sustainable AAC (allowable annual cut) targets:

| Management Area | Treatment            | Rate (%)   |
|------------------|----------------------|------------|
| 1                | CPRS                 | 0.22648%   |
| 1                | SeedTree             | 0.00310%   |
| 1                | CommercialThinning   | 0.02076%   |
| 1                | PC1                  | 0.07432%   |
| 1                | PC2                  | 0.15954%   |
| 1                | PC3                  | 0.22646%   |

---
