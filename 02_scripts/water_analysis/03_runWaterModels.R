# runWaterModels.R
# Margaret Swift <margaret.swift@cornell.edu>
# # run some models on waterhole raster fill level statistics

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
source(here::here('02_scripts', 'utilities.R'))
i_am('02_scripts/water_analysis/03_runWaterModels.R')
outdir <- here('03_output', 'water_analysis')
load(here(outdir, "hydrostats.rdata"))
pacman::p_load(LaplacesDemon,emmeans,
               brms, parameters, tidybayes, tinytable, betareg)
# ******************************************************************************

# ESW report

stats %>% 
  st_drop_geometry() %>% 
  group_by(TYPE, PERIOD) %>% 
  summarize(m=mean(ESW_SIZE_M2/10000), 
            sd=sd(ESW_SIZE_M2/10000), 
  )
stats %>% 
  filter(TYPE == "HYDROSHED") %>% 
  st_drop_geometry() %>% 
  group_by(PERIOD, YEAR) %>% 
  summarize(sum(ESW_SIZE_M2/10000)) %>% 
  View()


# ******************************************************************************

## HOW DIFFERENT are water fill levels by pentad? Do this by SIZ and hydroshed 
# separately

## How does cover percentage change over period for SIZ?
logidata.siz <- logidata %>% filter(TYPE=="SIZ")
aov.beta.siz <- aov(COVER_P ~ PERIOD, data = logidata.siz)

## How does cover percentage change over period for HYDROSHEDS?
aov.beta.hy <- aov(COVER_P ~ PERIOD, data = logidata.hy)

# need to run a Kruskal-Wallis test to see if any groups are different
logidata.hy <- logidata %>% filter(TYPE!="SIZ")

kruskal.test(COVER_P ~ PERIOD, data = logidata.hy)
# ... and a pairwise Wilcox test to see which ones
pairwise.wilcox.test(logidata.hy$COVER_P, logidata.hy$PERIOD,
                     p.adjust.method = "bonf")


## HOW DOES MAR INFLUENCE COVER PERCENTAGE?
# run this model with SIZ and Hydrosheds together

# data are zero-inflated, so let's first run a logistic model to see if there
## is a significant effect of MAR on whether there is any water at all.
logidata.hy <- logidata %>% filter(TYPE=="HYDROSHED")
m.logi <- glm(ISFILL ~ MAR_D * PERIOD, family=binomial(), data=logidata.hy)
summary(m.logi)

# next run a beta model just on (0,1) data:
betadata.hy <- betadata %>% filter(TYPE=="HYDROSHED") %>% filter(COVER_P<1)
m.beta <- betareg(COVER_P ~ MAR_D * PERIOD, 
                  data = betadata.hy)
summary(m.beta)
expcoefs = exp(coef(m.beta))

# marginal effect of MAR in general
cMAR = expcoefs[2]
cMAR

# marginal effect of MAR during FEB-MAY time period
cFM = expcoefs[3]
cMARFM = expcoefs[7]
change_in_OR_FEBMAY = cFM*cMARFM
change_in_OR_FEBMAY
change_in_OR_FEBMAY * cMAR

# marginal effect of MAR during MAY-JULY time period
cMJ = expcoefs[4]
cMARMJ = expcoefs[8]
change_in_OR_MAYJULY =cMJ*cMARMJ
change_in_OR_MAYJULY
change_in_OR_MAYJULY * cMAR
