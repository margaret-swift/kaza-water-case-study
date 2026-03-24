# cleanEle.R
# Created 26 Sept 2023
# Margaret Swift <margaret.swift@cornell.edu>
# 
# Purpose: To clean up elephant GPS data into a usable format

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
source(here::here('02_scripts/utilities.R'))
i_am('02_scripts/elephant_analysis/00_cleanEle.R')

# elephant data
rawdatadir <- here('01_data', 'elephant', 'raw')
procdatadir<- here('01_data', 'elephant', 'processed')
ele.nam <- read.csv(here(rawdatadir, 'nam.eles.fences.csv'))
ele.bots <- read.csv(here(rawdatadir, 'ecoexist.border.fence.pts.csv'))
ng13 <- st_read(here('01_data', 'ng13', 'raw', 'NG13_AOI.shp'))

# set CRS parameters
CRS_data <- 32734
CRS_utm <- 20935

# ******************************************************************************
#                       Pull together initial data cleaning
# ******************************************************************************
message('  initial data cleaning ...')
ele <- rbind(ele.nam, ele.bots) %>% 
  rename_all(toupper) %>% 
  arrange(ID, DATE.TIME) %>% 
  mutate( 
    # add index to rows
    INX =  1:nrow(.),
    
    # Date-time coordination
    DATE.TIME = as.POSIXlt(DATE.TIME, tz="UTC"),
    DATE = date(DATE.TIME),
    MONTH = month(DATE.TIME),
    YEAR = year(DATE.TIME),
    DTS.TRUE = DATE.TIME - lag(DATE.TIME),
    DTS.TRUE = ifelse(ID != lag(ID), NA, DTS.TRUE),
    DTS = DT, 
    DTM = DTS / 60,
    MPS = DIST / DTS,
    YEAR = year(DATE.TIME), 
    SEASON = ifelse(MONTH %in% 5:10, 'DRY', 'WET'),
    
    # reassigning IDs and fixing sex data
    ANIMAL_ID = ID,
    SEX = toupper(SEX)
  ) %>%   
  dplyr::select(INX, ID, ANIMAL_ID, SEX, 
                DATE.TIME, SEASON, 
                X, Y, DTS.TRUE, DTS, DIST, DTM, REL.ANGLE, SPEED, 
                COUNTRY, PROVIDER)

# anonymize IDs
ele$ID <- match(ele$ANIMAL_ID, unique(ele$ID))

# Assign bursts to data. Change burst whenever ID changes or if there is 
# a 400min or larger gap in the data
ele <- ele %>% 
  mutate(
    START.COUNT = ifelse(is.na(DTS.TRUE), TRUE, DTS.TRUE > 400),
    BURST = ID*1000 + cumsum(START.COUNT)) %>% 
  relocate(BURST, START.COUNT, .after="ID")

# transform to SF type
ele.sf <- ele %>% 
  st_as_sf(coords=c('X', 'Y'), crs=CRS_data, remove=FALSE) %>% 
  st_transform(crs=CRS_utm)

# add latlon
latlon <- st_coordinates(ele.sf); colnames(latlon) <- c('LON', "LAT")
ele.sf <- cbind(ele.sf, latlon) %>% relocate(LON, LAT, .after=Y)
ele.df <- cbind(ele, latlon) %>% relocate(LON, LAT, .after=Y)

# ******************************************************************************
#                                       STS
# ******************************************************************************

message('  saving rdata file ...')
save(ele.df, ele.sf, file=here(procdatadir, 'elephant.rdata'))
message("FINISHED cleanEle.R")
