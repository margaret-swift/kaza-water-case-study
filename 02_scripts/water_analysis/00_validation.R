# validation.R
# Margaret Swift <margaret.swift@cornell.edu>
# 
# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
source(here::here('02_scripts/utilities.R'))
i_am('02_scripts/water_analysis/00_validation.R')
pacman::p_load(tidyverse, reshape2)

datadir <- here('01_data', 'surfacewater', 'raw', 'validation')
regions = c('A_Kalala', 'B_Kwando', 'C_Okavango', 'D_Chobe', 
            'E_Zambezi', 'F_Ngamiland', 'G_Makgadikgadi')

##############################################################
#                       VALIDATION DATA 
##############################################################

massageData = function(name) {
  fname = paste0("Validation_Data_", name, ".csv")
  message('massaging data from ', fname, '...')
  data = read.csv(here(datadir, fname))
  data_m = data %>% rename_all(toupper) %>% 
    mutate(SEASON = factor(SEASON, levels=c('DRY', 'WET')),
           TYPE = factor(TYPE, levels=c('WATER', 'NOTWATER')), 
           P_TRUE = round(P_TRUE*100, 0)) %>% 
    filter(!grepl('mountains', REGION)) %>% 
    dcast(REGION ~ SEASON + TYPE, 
          value.var = "P_TRUE")
  data_m$REGION = regions
  message('DONE!')
  return(data_m)
}
esw_df = massageData('ESW')
gsw_df = massageData('GSW')

esw_df[,2:5] - gsw_df[,2:5]
