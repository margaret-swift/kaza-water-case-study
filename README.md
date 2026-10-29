# kaza-water-case-study

Data and scripts used in the elephant water-use case study for: 

Swift, M.E., Songhurst, A., McCulloch, G., Beytell, P., Naidoo, R. 2026.
  "Mapping small-scale ephemeral surface water to inform transfrontier conservation 
  planning in southern Africa". Ecosphere.

  ECS26-0338

This repository holds code for two analyses run for this paper: water_analysis and elephant_analysis. 

## 01_data

Contains all data for analysis. Elephant GPS collar data are protected and not included in this repository.

- elevation - elevation data
- fences - fence shapes
- kaza - KAZA boundary shapefile
- ng13 - area of interest
- rainfall - rainfall data
- surfacewater - ESW

## water_analysis

Includes code to run analysis in Methods Sec. 2.7 and Results Sec. 3.3 (Seasonal surface water fill analysis) 

- 00_validation.R - takes water validation data and turns it into tables
- 01_createRasterStats.R - creates datasets for further water models
- 02_createPlots.R - creates plots Fig. 4 and Fig. 5
- 03_runWaterModels.R - runs GLMMs for seasonal water availability analysis

## elephant_analysis

Includes code to run analysis in Methods Sec. 2.8	and Results Sec. 3.4 (Savanna elephant water use analysis).

- 00_cleanEle.R - cleans up elephant GPS collar data, removes outliers and filters to only hourly GPS collar data from female elephants.
- 01_cleanGeographicData.R - cleans geographic data for use in further analyses
- 02_addWaterDistsToEle.R - takes elephant GPS collar data and calculates distance-from-water for several water datasets: MNDWI, NDWI2, WI2015, ESWr, ESWr15, GSW, RSW, KAZA05, and KAZA21.
- 03_runEleCaseStudies.R - runs GLMMs for elephant water use analysis and creates Fig. 7
