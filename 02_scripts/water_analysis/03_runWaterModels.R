# runWaterModels.R
# Margaret Swift <margaret.swift@cornell.edu>
# # run some models on waterhole raster fill level statistics

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
source(here::here('scripts', 'utilities.R'))
i_am('02_scripts/water_analysis/03_runWaterModels.R')
outdir <- here('03_output', 'water_analysis')
load(here(outdir, "hydrostats.rdata"))
pacman::p_load(LaplacesDemon,emmeans,
               brms, parameters, tidybayes, tinytable, betareg)

# ******************************************************************************

## HOW DIFFERENT are water fill levels by pentad? Do this by SIZ and hydroshed 
# separately

## How does cover percentage change over period for SIZ?
logidata.siz <- logidata %>% filter(TYPE=="SIZ")
aov.beta.siz <- aov(COVER_P ~ PERIOD, data = logidata.siz)
TukeyHSD(aov.beta.siz)
GGTukey(TukeyHSD(aov.beta.siz))

## How does cover percentage change over period for HYDROSHEDS?
logidata.hy <- logidata %>% filter(TYPE!="SIZ")
aov.beta.hy <- aov(COVER_P ~ PERIOD, data = logidata.hy)
TukeyHSD(aov.beta.hy)
GGTukey(TukeyHSD(aov.beta.hy))


## HOW DOES MAR INFLUENCE COVER PERCENTAGE?
# run this model with SIZ and Hydrosheds together

# data are zero-inflated, so let's first run a logistic model to see if there
## is a significant effect of MAR on whether there is any water at all.
m.logi <- glm(ISFILL ~ MAR_D * PERIOD * TYPE, family=binomial(), data=logidata)
summary(m.logi)

# next run a beta model just on (0,1) data:
m.beta <- betareg(COVER_P ~ MAR_D * PERIOD * TYPE, data = betadata)
summary(m.beta)
m.beta2 <- betareg(COVER_P ~ MAR_D * PERIOD, data = betadata)
summary(m.beta2)
expcoefs = exp(coef(m.beta2))

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
