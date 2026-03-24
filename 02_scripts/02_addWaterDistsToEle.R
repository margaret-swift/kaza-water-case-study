# addWaterDistsToEle.R
# Created 24 Aug 2025
# Margaret Swift <margaret.swift@cornell.edu>
# 
#  Combines the water and elephant datasets

# ******************************************************************************
#                                LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/03_addWaterDistsToEle.R')
pacman::p_load(sf, nngeo) # nngeo for nearest-neighbor analyses
CRS_utm = 20935

path_ele <- here('01_data', 'elephant')
load(here(path_ele, 'processed', 'elephant.rdata'))
load(here('01_data', 'fences', 'processed', 'fences_utm.rdata'))
load(here('01_data', 'ng13', 'processed', 'ng13_utm.rdata'))
load(here('01_data', 'surfacewater', 'processed', 'NG13_ESW_utm.rdata'))
load(here('01_data', 'surfacewater', 'processed', 'NG13_GSW_utm.rdata'))

# ******************************************************************************
#                               LOADING OTHER DATA
# ******************************************************************************


# We need to create a buffer around the AOI to catch water just outside
aoi <- st_bbox(esw.utm.df) %>% st_as_sfc() %>% st_as_sf()
aoi.buff = aoi %>% st_buffer(1000)

#  Let's split the AOI in two along the fence
aoi.split <- lwgeom::st_split(aoi.buff, fences_utm) %>% 
  st_collection_extract(c("POLYGON"))
fences.aoi <- fences_utm %>% st_crop(aoi.buff)

# load elephant data and crop to buffer
ele.dat <- ele.sf %>%
  st_crop(aoi.buff) %>%
  mutate(INX = 1:n())

# create random data frames for points inside our buffer zone
rsw_1800_1 = st_sample(aoi.buff, size=1800)
rsw_1800_2 = st_sample(aoi.buff, size=1800)
rsw_1800_3 = st_sample(aoi.buff, size=1800)
rsw_2500_1 = st_sample(aoi.buff, size=2500)
rsw_2500_2 = st_sample(aoi.buff, size=2500)
rsw_2500_3 = st_sample(aoi.buff, size=2500)
rsw_5000_1 = st_sample(aoi.buff, size=5000)
rsw_5000_2 = st_sample(aoi.buff, size=5000)
rsw_5000_3 = st_sample(aoi.buff, size=5000)

# now run dist to water on both sides of the fence
rsw_names <- c('RSW_DIST_1800_1', 'RSW_DIST_1800_2', 'RSW_DIST_1800_3', 
               'RSW_DIST_2500_1', 'RSW_DIST_2500_2', 'RSW_DIST_2500_3', 
               'RSW_DIST_5000_1', 'RSW_DIST_5000_2', 'RSW_DIST_5000_3')
esw_names = c('ESW_NN_ID', 'ESW_NN_DIST', 'ESW_NN_RECUR')
eswp_names = c('ESWP_NN_ID', 'ESWP_NN_DIST', 'ESWP_NN_RECUR')
gsw_names = c('GSW_NN_ID', 'GSW_NN_DIST', 'GSW_NN_RECUR')
frq_names = "MFREQ"
colnames <- c(esw_names, eswp_names, gsw_names, frq_names, rsw_names)
nn.dat = matrix(NA, nrow=nrow(ele.dat), ncol=length(colnames)) %>% as.data.frame()
names(nn.dat) = colnames

# function to find nearest neighbors
nnVectors <- function(epts, wpts) {
  nn = st_nn(epts, wpts, returnDist=TRUE)
  ids = unlist(nn$nn)
  freq = wpts$recurrence[ids]
  dist = unlist(nn$dist)
  return(data.frame(NN_ID = ids, NN_DIST = dist, NN_FREQ = freq))
}

for (i in 1:nrow(aoi.split)) {
  roi <- aoi.split[i,]
  
  # crop ele data to AOI
  nn.roi <- ele.dat %>% st_crop(roi)
  inx = nn.roi$INX
  
  # all esw
  water.roi = esw.utm.df %>% st_crop(roi) %>% rename(recurrence=isWater)
  nn_esw = nnVectors(nn.roi, water.roi)
  nn.dat[inx,esw_names] <- nn_esw
  
  # esw permanent
  water.roi.perm <- water.roi %>% filter(recurrence > 15)
  nn_eswp = nnVectors(nn.roi, water.roi.perm)
  nn.dat[inx,eswp_names] <- nn_eswp
  
  # gsw
  gsw.roi = gsw.utm.df %>% st_crop(roi)
  nn_gsw = nnVectors(nn.roi, gsw.roi)
  nn.dat[inx,gsw_names] <- nn_gsw
  
  # randoms
  # function for nearest neighbor random points
  nnRand = function(rsw_rand) {
    rsw_roi = rsw_rand %>% st_crop(roi)
    nn_rsw = st_nn(nn.roi, rsw_roi, returnDist=TRUE)
    return(nn_rsw$dist %>% unlist())
  }
  
  nn_rsw = cbind( nnRand(rsw_1800_1),  nnRand(rsw_1800_2),  nnRand(rsw_1800_3),
                  nnRand(rsw_2500_1),  nnRand(rsw_2500_2),  nnRand(rsw_2500_3),
                  nnRand(rsw_5000_1),  nnRand(rsw_5000_2),  nnRand(rsw_5000_3))
  nn.dat[inx,rsw_names] <- nn_rsw
  
  # freqs
  nn = st_nn(nn.roi, water.roi, k=10, returnDist=TRUE)
  nn_mfreq = unlist(lapply(nn$nn, function(e) mean(water.roi$recurrence[e], na.rm=TRUE)))
  nn.dat[inx,frq_names] <- nn_mfreq
}

# Add NN data and crop it back to the original non-buffer AOI.
ele.dat <- cbind(ele.dat, nn.dat) %>% st_crop(aoi)

# ******************************************************************************
# combine ele data from inside and outside the fence, then 
# remove bursts with too little data or edge cases.
to.rm <- c(7212:7242, 8246, 480059)
b.df <- ele.dat %>% st_drop_geometry() %>% 
  group_by(BURST) %>% summarize(n=n()) %>% 
  filter(n > 10, !(BURST %in% to.rm))
bursts <- b.df$BURST

ele.dat = ele.dat %>% 
  filter(BURST %in% bursts) %>% 
  mutate(ID=factor(ID),
         ANIMAL_ID = factor(ANIMAL_ID),
         DATE = date(DATE.TIME),
         YEAR = year(DATE.TIME),
         SEASON = factor(SEASON))

# save
save(ele.dat, file=here(path_ele, 'processed', 'NG13_elephant.rdata'))
