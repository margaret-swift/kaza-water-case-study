# createRasterStats.R
# Margaret Swift <margaret.swift@cornell.edu>
# 
# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
source(here::here('02_scripts/utilities.R'))
i_am('02_scripts/water_analysis/01_createRasterStats.R')
datadir <- here('01_data', 'surfacewater', 'raw')
outdir <- here('03_output', 'water_analysis')

##############################################################
#                       RAINFALL DATA
##############################################################

# load shapefiles
hydro_shp = st_read(here('01_data', 'elevation', 'raw','all_regions_dem.shp')) %>% 
  mutate(ID = ifelse(is.na(ID), HYBAS_ID, ID),
         ID = as.character(ID)) %>% 
  dplyr::select(ID, mean_slope, min_slope, max_slope, 
                mean_SRTM, min_SRTM, max_SRTM)

# load rainfall data for hydrosheds
mwf_stats <- read.csv(here(datadir,'HydroshedMWFStatsSlim.csv')) %>% 
  rename(ID=HYBAS_ID) %>% mutate(ID = as.character(ID))
hydro_mar <- read.csv(here('01_data', 'rainfall', 'raw','Hydroshed_MAR.csv')) %>% 
  relocate('HYBAS_ID') %>% select(-system.index, -.geo)
hydro_mar_m = reshape2::melt(hydro_mar, id.var="HYBAS_ID", 
                             variable.name="YEAR", value.name="MAR") %>% 
  mutate(YEAR = as.numeric(gsub("X", "", YEAR))) %>% 
  rename(ID=HYBAS_ID)

# load rainfall data for SIZ catchments
siz_catch_mar <- read.csv(here('01_data','rainfall','raw','SIZ_Catchments_MAR.csv')) %>% 
  relocate('name') %>% select(-system.index, -.geo)

# Chobe and Linyanti SIZ also use Zambezi catchment; 
#   Makgadikgadi gets water from both Nata and Okavango
zambezi = siz_catch_mar %>% filter(name=="Zambezi")
okavango = siz_catch_mar %>% filter(name=="Okavango")
nata = siz_catch_mar %>% filter(name=="Nata")
makgadikgadi = siz_catch_mar %>% 
  filter(name != "Zambezi") %>% 
  select(-name) %>% 
  colMeans() %>% t() %>% 
  as.data.frame() %>% 
  mutate(name="Makgadikgadi") %>% 
  relocate(name)
siz_mar_m = rbind(
  okavango,
  zambezi, 
  zambezi %>% mutate(name="Chobe"),
  zambezi %>% mutate(name="Linyanti"),
  makgadikgadi 
) %>% 
  reshape2::melt(id.var="name", 
                variable.name="YEAR", value.name="MAR") %>% 
  mutate(YEAR = as.numeric(gsub("mm", "", YEAR))) %>% 
  rename(ID=name)

#combine them all
mar_m = rbind(siz_mar_m, hydro_mar_m)


##############################################################
#                       ESW FILL DATA
##############################################################
# load hydroshed data
hdata = NULL
files <- list.files(here(datadir, 'ESW_sizes_HY'), full.names=TRUE)
for (f in files) {
  data <- read.csv(f) %>% 
    group_by(ID, START, END, YEAR) %>% 
    summarize(MWF_SIZE = mean(MWF_SIZE),
              ESW_SIZE = sum(ESW_SIZE)) %>% 
    mutate(TYPE='HYDROSHED',
           PERIOD = paste(START, END, sep="_"),
           START = ymd(paste0(YEAR, '-', START, "-01")),
           END = ymd(paste0(YEAR, '-', END, "-01")), 
           ID = as.character(ID)) %>% 
    select(ID, TYPE, PERIOD, YEAR, START, END,
           ESW_SIZE, MWF_SIZE) %>% 
    ungroup()
  if (is.null(hdata)) hdata <- data
  else hdata <- rbind(hdata, data)
}

# define periods
pds = c('NOV_FEB', 'FEB_MAY', 'MAY_JUL', 'JUL_SEP', 'SEP_NOV')
pds.df <- data.frame(MONTH=1:12, 
                     PERIOD=c('NOV_FEB', 'FEB_MAY', 'FEB_MAY', 'FEB_MAY', 'MAY_JUL', 'MAY_JUL',
                              'JUL_SEP', 'JUL_SEP', 'SEP_NOV', 'SEP_NOV', 'NOV_FEB', 'NOV_FEB'))

# load SIZ data
files <- list.files(here('data', 'ESW_sizes_SIZ'), full.names=TRUE)
sizdata = NULL
for (f in files) {
  name = gsub('.*\\/|_ESW_sizes.csv.*', '', f)
  data <- read.csv(f) %>% 
    mutate(TYPE='SIZ', MWF_SIZE=NA,
           START=as.POSIXct(as.Date(START, "%m/%d/%y")),
           END=as.POSIXct(as.Date(END, "%m/%d/%y")),
           ) %>% 
    left_join(pds.df, by="MONTH") %>% 
    select_at(names(hdata))
  data$MWF_SIZE = max(data$ESW_SIZE)
  if (is.null(sizdata)) sizdata <- data
  else sizdata <- rbind(sizdata, data)
}

alldata = rbind(hdata, sizdata) 

##############################################################
#                       DOUBLE CHECKING
##############################################################

# massage data
myAnomaly = function(dat) {
  mu = mean(dat, na.rm=TRUE)
  chg = (dat-mu) / mu
  return(chg)
}
stats <- alldata %>% 
  group_by(ID) %>% 
  left_join(mar_m, by=c("ID", "YEAR")) %>% 
  left_join(mwf_stats, by=c("ID")) %>%
  mutate(
    MWF_SIZE = max(ESW_SIZE),
    ID = factor(ID),
    PERIOD = factor(PERIOD, levels=pds),
    COVER_P = (ESW_SIZE) / MWF_SIZE,
    COVER_D = myAnomaly(ESW_SIZE),
    MAR_D = myAnomaly(MAR),
    YEAR = factor(YEAR, levels=2019:2025),
    ISFILL = COVER_P > 0,
    AVG_ESW_SIZE = mean(ESW_SIZE), 
    MONTH=month(START)
  ) %>% 
  # remove hydrosheds with max fill less than 1000 pixels
  filter(MWF_SIZE > 1000)

fill_stats = stats %>% 
  group_by(ID) %>% 
  summarize(mean_fill = mean(COVER_P, na.rm=TRUE), 
            sd_fill = sd(COVER_P, na.rm=TRUE)) %>% 
  mutate(INX = 1:n())

sdat = stats %>% st_drop_geometry() %>% 
  group_by(TYPE, PERIOD, YEAR) %>% 
  summarize(cover_mu = mean(COVER_P, na.rm=TRUE), 
            cover_sd = sd(COVER_P, na.rm=TRUE), 
            mar_mu = mean(MAR),
            mar_sd = sd(MAR))
stats = stats %>% 
  left_join(fill_stats, by="ID") %>% 
  left_join(hydro_shp, by="ID") %>% 
  st_as_sf()

logidata = stats %>% st_drop_geometry() %>% as.data.frame()
betadata = logidata %>% filter(ISFILL)
save(stats, fill_stats, sdat, hdata, sizdata, alldata,
     siz_catch_mar, mar_m, 
     logidata, betadata, file=here(outdir, "hydrostats.rdata"))
