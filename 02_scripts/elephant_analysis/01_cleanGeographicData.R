# cleanGeographicData.R
# Created 26 Sept 2023
# Margaret Swift <margaret.swift@cornell.edu>
# 
# Purpose: To clean up fence, water, and AOI data into a usable format

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
source(here::here('02_scripts/utilities.R'))
i_am('02_scripts/elephant_analysis/01_cleanGeographicData.R')

ng13 <- st_read(here('01_data', 'ng13', 'raw', 'NG13_AOI.shp'))
fences <- st_read(here('01_data', 'fences', 'raw', 'fences.shp'))

path_water <- here('01_data', 'surfacewater')
esw_raw <- terra::rast(here(path_water, 'raw', 'NG13_ESW.tif'))
gsw_raw <- terra::rast(here(path_water, 'raw', 'NG13_GSW.tif'))


# ******************************************************************************
#                                 TRANSFORM DATA
# ******************************************************************************
fences_utm <- fences %>% st_transform(crs=CRS_utm)
save(fences_utm, file=here('01_data', 'fences', 'processed', 'fences_utm.rdata'))
ng13_utm <- ng13 %>% st_transform(crs=CRS_utm)
save(ng13_utm, file=here('01_data', 'ng13', 'processed', 'ng13_utm.rdata'))

esw_utm = terra::project(esw_raw, CRS_utm)
gsw_utm = terra::project(gsw_raw, fences_utm)
# terra::writeRaster(esw.utm, here(path_water, 'processed', 'NG13_ESW_utm.tif'))
# terra::writeRaster(esw.utm, here(path_water, 'processed', 'NG13_GSW_utm.tif'))

# rasters to points
esw.utm.df = esw_utm %>% 
  as.data.frame(xy = TRUE) %>%
  na.omit() %>%
  st_as_sf(coords=c('x', 'y'), crs=CRS_utm, remove=FALSE) %>% 
  filter(recurrence > 0)
gsw.utm.df = gsw_utm %>% 
  as.data.frame(xy = TRUE) %>%
  na.omit() %>%
  st_as_sf(coords=c('x', 'y'), crs=CRS_utm, remove=FALSE) %>% 
  filter(recurrence > 0)
save(esw.utm.df, file=here(path_water, 'processed', 'NG13_ESW_utm.rdata'))
save(gsw.utm.df, file=here(path_water, 'processed', 'NG13_GSW_utm.rdata'))



# EOF