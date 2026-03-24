# eleCaseStudies.R
# Created 24 Aug 2025
# Margaret Swift <margaret.swift@cornell.edu>
# 
# This script runs all the analyses given in the Elephant Case Studies section 
# of the manuscript.

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/03_runEleCaseStudies.R')
pacman::p_load(terra, lme4)
CRS_utm = 20935

path_ele <- here('01_data', 'elephant')
load(here(path_ele, 'processed', 'NG13_elephant.rdata'))
load(here('01_data', 'fences', 'processed', 'fences_utm.rdata'))
load(here('01_data', 'ng13', 'processed', 'ng13_utm.rdata'))
load(here('01_data', 'surfacewater', 'processed', 'NG13_ESW_utm.rdata'))
load(here('01_data', 'surfacewater', 'processed', 'NG13_GSW_utm.rdata'))

# my thresholds
time_T = 48 #(hours)
dist_T = 50 #(meters)

# ******************************************************************************
#                                     FUNCTIONS
# ******************************************************************************

mymean <- function(x) mean(x[!is.infinite(x)], na.rm=TRUE)
createEle13 <- function(dist_T) {
  # Creates elephant data summaries - how many points are near water,
  # adds thirst counter, etc. 
  ele.13 <- ele.dat %>% 
    mutate(
       # fixrate is 1 hour; added in case of later analyses with other fix rates
       FIX = 1,
       # Flag - is current point within dist_T of water?
       NEAR_GSW = GSW_NN_DIST < dist_T,
       NEAR_ESW = ESW_NN_DIST < dist_T,
       NEAR_ESWP = ESWP_NN_DIST < dist_T,
       NEAR_RSW_1800_1 = RSW_DIST_1800_1 < dist_T,
       NEAR_RSW_1800_2 = RSW_DIST_1800_2 < dist_T,
       NEAR_RSW_1800_3 = RSW_DIST_1800_3 < dist_T,
       NEAR_RSW_2500_1 = RSW_DIST_2500_1 < dist_T,
       NEAR_RSW_2500_2 = RSW_DIST_2500_2 < dist_T,
       NEAR_RSW_2500_3 = RSW_DIST_2500_3 < dist_T,
       NEAR_RSW_5000_1 = RSW_DIST_5000_1 < dist_T,
       NEAR_RSW_5000_2 = RSW_DIST_5000_2 < dist_T,
       NEAR_RSW_5000_3 = RSW_DIST_5000_3 < dist_T,
       # This is the elephant water counter - count up (cumsum) how many 
       # steps it has been since NEAR_water was TRUE
       GP_ESW = cumsum(START.COUNT | NEAR_ESW),
       GP_ESWP = cumsum(START.COUNT | NEAR_ESWP),
       GP_GSW = cumsum(START.COUNT | NEAR_GSW),
       GP_RSW_1800_1 = cumsum(START.COUNT | NEAR_RSW_1800_1),
       GP_RSW_1800_2 = cumsum(START.COUNT | NEAR_RSW_1800_2),
       GP_RSW_1800_3 = cumsum(START.COUNT | NEAR_RSW_1800_3),
       GP_RSW_2500_1 = cumsum(START.COUNT | NEAR_RSW_2500_1),
       GP_RSW_2500_2 = cumsum(START.COUNT | NEAR_RSW_2500_2),
       GP_RSW_2500_3 = cumsum(START.COUNT | NEAR_RSW_2500_3),
       GP_RSW_5000_1 = cumsum(START.COUNT | NEAR_RSW_5000_1),
       GP_RSW_5000_2 = cumsum(START.COUNT | NEAR_RSW_5000_2),
       GP_RSW_5000_3 = cumsum(START.COUNT | NEAR_RSW_5000_3),
    ) %>% 
    # now for each thirst counter grouping, find out the max # hours of thirst
    # the elephant reaches.
    ungroup() %>% group_by(GP_ESW) %>% mutate(THIRST_ESW = cumsum(!NEAR_ESW) * FIX) %>% 
    ungroup() %>% group_by(GP_ESWP) %>% mutate(THIRST_ESWP= cumsum(!NEAR_ESWP) * FIX) %>% 
    ungroup() %>% group_by(GP_GSW) %>% mutate(THIRST_GSW = cumsum(!NEAR_GSW) * FIX) %>% 
    ungroup() %>% group_by(GP_RSW_1800_1) %>% mutate(THIRST_RSW_1800_1 = cumsum(!NEAR_RSW_1800_1) * FIX) %>% 
    ungroup() %>% group_by(GP_RSW_1800_2) %>% mutate(THIRST_RSW_1800_2 = cumsum(!NEAR_RSW_1800_2) * FIX) %>% 
    ungroup() %>% group_by(GP_RSW_1800_3) %>% mutate(THIRST_RSW_1800_3 = cumsum(!NEAR_RSW_1800_3) * FIX) %>% 
    ungroup() %>% group_by(GP_RSW_2500_1) %>% mutate(THIRST_RSW_2500_1 = cumsum(!NEAR_RSW_2500_1) * FIX) %>% 
    ungroup() %>% group_by(GP_RSW_2500_2) %>% mutate(THIRST_RSW_2500_2 = cumsum(!NEAR_RSW_2500_2) * FIX) %>% 
    ungroup() %>% group_by(GP_RSW_2500_3) %>% mutate(THIRST_RSW_2500_3 = cumsum(!NEAR_RSW_2500_3) * FIX) %>% 
    ungroup() %>% group_by(GP_RSW_5000_1) %>% mutate(THIRST_RSW_5000_1 = cumsum(!NEAR_RSW_5000_1) * FIX) %>% 
    ungroup() %>% group_by(GP_RSW_5000_2) %>% mutate(THIRST_RSW_5000_2 = cumsum(!NEAR_RSW_5000_2) * FIX) %>% 
    ungroup() %>% group_by(GP_RSW_5000_3) %>% mutate(THIRST_RSW_5000_3 = cumsum(!NEAR_RSW_5000_3) * FIX) %>% 
    ungroup()
  return(ele.13)
}
makeAOVDat <- function(df) {
  # Makes AOV data for Tukey HSD analyses
  typelist = c('ESW', 'GSW', 'RSW_1800_1', 'RSW_2500_1', 'RSW_5000_1')
  inx = which(grepl("THIRST|^ID$|DATE.TIME", names(df)))
  thirst.df = df[,inx] %>% st_drop_geometry() %>% 
    mutate(DATE.TIME = as.POSIXct(DATE.TIME)) %>% 
    reshape2::melt(id.vars=c('ID', 'DATE.TIME'), 
                   variable.name="watertype", value.name="thirst") %>% 
    arrange(ID, DATE.TIME) %>% 
    mutate(watertype = gsub("THIRST_", "", watertype)) %>% 
    filter(watertype %in% typelist) %>% 
    mutate(watertype = gsub("_1", "", watertype))
  return(thirst.df)
}
GGTukey <- function(Tukey){
  # a pretty version of Tukey HSD plot
  A<-require("tidyverse")
  if(A==TRUE){
    library(tidyverse)
  } else {
    install.packages("tidyverse")
    library(tidyverse)
  }
  B<-as.data.frame(Tukey[1])
  colnames(B)[2:3]<-c("min",
                      "max")
  C<-data.frame(id=row.names(B),
                min=B$min,
                max=B$max)
  D<-C%>%
    ggplot(aes(id))+
    geom_errorbar(aes(ymin=min,
                      ymax=max),
                  width = 0.2)+
    geom_hline(yintercept=0,
               color="red")+
    labs(x=NULL)+
    coord_flip()+
    theme(text=element_text(family="TimesNewRoman"),
          title=element_text(color="black",size=15),
          axis.text = element_text(color="black",size=10),
          axis.title = element_text(color="black",size=10),
          panel.grid=element_line(color="grey75"),
          axis.line=element_blank(),
          plot.background=element_rect(fill="white",color="white"),
          panel.background=element_rect(fill="white"),
          panel.border = element_rect(colour = "black", fill = NA,size=0.59),
          legend.key= element_rect(color="white",fill="white")
    )
  return(D)
}
makeStats <- function(dist_T) {
  # A table of visitation rates within 48 hours
  df = createEle13(dist_T)
  stats_all=df %>% st_drop_geometry() %>% 
    summarize(total = n(), 
              ESW_n48 = sum(THIRST_ESW <= time_T),
              GSW_n48 = sum(THIRST_GSW <= time_T),
              RSW_1800_1 = sum(THIRST_RSW_1800_1 <= time_T),
              RSW_1800_2 = sum(THIRST_RSW_1800_2 <= time_T),
              RSW_1800_3 = sum(THIRST_RSW_1800_3 <= time_T),
              RSW_2500_1 = sum(THIRST_RSW_2500_1 <= time_T),
              RSW_2500_2 = sum(THIRST_RSW_2500_2 <= time_T),
              RSW_2500_3 = sum(THIRST_RSW_2500_3 <= time_T),
              RSW_5000_1 = sum(THIRST_RSW_5000_1 <= time_T),
              RSW_5000_2 = sum(THIRST_RSW_5000_2 <= time_T),
              RSW_5000_3 = sum(THIRST_RSW_5000_3 <= time_T),
    ) %>% 
    mutate(ESW_p48 = ESW_n48 / total,
           GSW_p48 = GSW_n48 / total,
           RSW_1800_1_p48 = RSW_1800_1 / total,
           RSW_1800_2_p48 = RSW_1800_2 / total,
           RSW_1800_3_p48 = RSW_1800_3 / total,
           RSW_2500_1_p48 = RSW_2500_1 / total,
           RSW_2500_2_p48 = RSW_2500_2 / total,
           RSW_2500_3_p48 = RSW_2500_3 / total,
           RSW_5000_1_p48 = RSW_5000_1 / total,
           RSW_5000_2_p48 = RSW_5000_2 / total,
           RSW_5000_3_p48 = RSW_5000_3 / total,
    ) %>% 
    mutate(RSW_1800_p48 = mean(c(RSW_1800_1_p48, RSW_1800_2_p48, RSW_1800_3_p48)),
           RSW_2500_p48 = mean(c(RSW_2500_1_p48, RSW_2500_2_p48, RSW_2500_3_p48)),
           RSW_5000_p48 = mean(c(RSW_5000_1_p48, RSW_5000_2_p48, RSW_5000_3_p48)),
           THRESH=dist_T) %>% 
    dplyr::select(THRESH, ESW_p48, GSW_p48, RSW_1800_p48, RSW_2500_p48, RSW_5000_p48)
  return(stats_all)
}
runAllAnalyses <- function(dist_T) {
  # Create data, get stats, and run Tukey HSD analysis for threshold
  message("running analyses for dist_T = ", dist_T)
  ele.13.nt = createEle13(dist_T)
  thirst.nt = makeAOVDat(ele.13.nt)
  lm.nt =  lm(thirst ~ watertype, data=thirst.nt)
  avg.thirst = thirst.nt %>% group_by(watertype) %>% 
    summarize(mu=mean(thirst, na.rm=TRUE) %>% round(2), 
            sd=sd(thirst, na.rm=TRUE)%>% round(2)) %>% 
    mutate(val = paste0(mu, " + ", sd)) %>% 
    dplyr::select(watertype, val)
  print(avg.thirst)
  aov.nt =  aov(lm.nt)
  tuk.nt  = TukeyHSD(aov.nt)
  print(GGTukey(tuk.nt))
}

# ******************************************************************************
#                                   ANALYSES!
# ******************************************************************************

# Run Tukey HSD and AOV on water visitation data
timeSteps = c(50, 100, 200, 322, 644)
for (t in timeSteps) {
  runAllAnalyses(t)
}

# Create plot (Fig 5) for % of data under threshold for each T
mdat = do.call(rbind, lapply(timeSteps, makeStats)) %>% as.data.frame() %>% 
  rename_all(function(x) gsub('_p48','', x)) %>% 
  reshape2::melt(id.var="THRESH", variable.name="watertype")
ggplot(mdat, aes(x=THRESH, y=value*100, color=watertype)) + 
  geom_line(linewidth=1) + 
  geom_line(data=mdat %>% filter(watertype=="ESW"), color='red', linewidth=1.5) +
  geom_line(data=mdat %>% filter(watertype=="GSW"), color='blue', linewidth=1.5) +
  geom_point(color='black', size=3) + 
  xlab('distance (m) from water classified as visitation') + 
  ylab('% points within 48hours of water') + 
  scale_color_manual(values=c('red', 'blue', 'gray40', 'gray80', '#e6e6e6')) + 
  theme_minimal() + theme(text=element_text(size=14)) + 
  scale_x_continuous(breaks=timeSteps) + ylim(0, 100)

#visitation rates by sex
ele.13 = createEle13(322)
ele.13 %>% st_drop_geometry() %>% 
  group_by(SEX) %>% 
  summarize(total = n(), 
            ESW_n48 = sum(THIRST_ESW <= t),
            GSW_n48 = sum(THIRST_GSW <= t),
  ) %>% 
  mutate(ESW_p48 = ESW_n48 / total,
         GSW_p48 = GSW_n48 / total)

# visitation rates by season
ele.13 %>% st_drop_geometry() %>% 
  group_by(SEX, SEASON) %>% 
  summarize(total = n(), 
            ESW_n48 = sum(THIRST_ESW <= t),
            ESWP_n48 = sum(THIRST_ESWP <= t),
            GSW_n48 = sum(THIRST_GSW <= t),
            
  ) %>% 
  mutate(ESWP_p48 = ESWP_n48 / total,
         ESWE_p48 = 1 - ESWP_p48,
         ESW_p48 = ESW_n48 / total,
         GSW_p48 = GSW_n48 / total)

#  H2A: % of ele GPS points that drink every 48hr is higher for ESW than GSW
sdat2 <- ele.13 %>% group_by(SEX) %>% st_drop_geometry() %>% 
  summarize(n=n(),
            mhours_GSW = round(mean(THIRST_GSW)),
            mhours_esw = round(mean(THIRST_ESW)),
            n2days_GSW = round(sum(THIRST_GSW < t) / n * 100),
            n2days_esw = round(sum(THIRST_ESW < t) / n * 100)
  )
mGSW = round(mean(ele.13$THIRST_GSW))
mesw = round(mean(ele.13$THIRST_ESW))
ggplot(ele.13 %>% mutate(is_close = THIRST_GSW < t+2)) + 
  geom_histogram(aes(x=THIRST_GSW, fill=is_close), binwidth=19.8, color='white')  +
  geom_vline(xintercept=t, linetype='dashed', linewidth=0.7, color='#ED008D') +
  scale_x_continuous(n.breaks=15, limits=c(0, 1200)) +
  scale_fill_manual(values=c('black', 'purple')) +
  theme_classic() +
  guides(fill="none") + xlab('') + ylab('') + 
  theme(text=element_text(size=16))

ggplot(ele.13 %>% mutate(is_close = THIRST_ESW < t)) + 
  geom_histogram(aes(x=THIRST_ESW, fill=is_close), bins=100, color='white')  +
  scale_fill_manual(values=c('black', 'purple')) +
  geom_vline(xintercept=t, linetype='dashed', linewidth=0.7, color='#ED008D') +
  scale_x_continuous(n.breaks=15, limits=c(0, 100)) +
  scale_y_continuous(limits=c(0, 17000)) +
  theme_classic() +
  guides(fill="none") + xlab('') + ylab('') + 
  theme(text=element_text(size=16))



# ******************************************************************************
# We also tested whether, given our surface water maps, elephant rely on 
# seasonal surface water in the wet season and permanent surface water in the 
# dry season. To do so, we calculated the mean recurrence of each elephant 
# GPS point’s ten nearest neighbors (found using ESWr) and used a Generalized 
# Linear Mixed-effects Model (GLMM, using ‘glmmTMB’ with family = beta_family() 
# from the ‘glmer’ package in R) to model the response of the beta-distributed 
# mean recurrence to a binary season (“wet season” as the base category) and 
# elephant ID as the random effect. Our hypothesis was supported if the nearest
# mean recurrence increased significantly in the dry season (i.e., the 
# coefficient for “dry season”, beta_dry, was significantly positive).
# ******************************************************************************

ez <- ele.13 %>% st_drop_geometry() %>% mutate(FREQ = MFREQ / 30, 
                                  FREQ = ifelse(FREQ == 1, 0.99, FREQ))
m1 = lm(ESWP_NN_DIST ~ SEASON, data=ez)
summary(m1)
m2 = glmmTMB::glmmTMB(formula = FREQ ~ SEASON,
                 data = ez, 
                 family = glmmTMB::beta_family() )
summary(m2)
ez %>% group_by(SEASON) %>% summarize(mean(MFREQ), sd(MFREQ), mean(ESWP_NN_DIST), sd(ESWP_NN_DIST))


# distances from water
p1 = ggplot(ez) + geom_histogram(aes(x=ESW_NN_DIST),
                                 color='black',
                                 bins=200) + 
  scale_x_continuous(n.breaks=20) + 
  facet_wrap(~SEASON, ncol=1)
p2 = ggplot(ez) + geom_histogram(aes(x=ESWP_NN_DIST),
                                 color='black',
                                 bins=200) + 
  scale_x_continuous(n.breaks=20) + 
  facet_wrap(~SEASON, ncol=1)
p3 = ggplot(ez) + geom_histogram(aes(x=GSW_NN_DIST),
                                 color='black',
                                 bins=200) +
  scale_x_continuous(n.breaks=20) + 
  facet_wrap(~SEASON, ncol=1)
library(patchwork)
p1 + p2 + p3

ggplot(ez) + geom_histogram(aes(x=DIST_MIN_ESW),
                            color='black',
                            bins=100) + 
  facet_wrap(~ID)



# ******************************************************************************
#                                 MODELING DISTANCE FREQUENCIES
# ******************************************************************************
# 
m1 = mymean(ele.13$ESWP_NN_DIST)
m2 = mymean(ele.13$GSW_NN_DIST)
ele.norm <- ele.13 %>% st_drop_geometry() %>%
  mutate(TOWARDS_PERM = ifelse(START.COUNT, NA, ESWP_NN_DIST < lag(ESWP_NN_DIST)),
         TOWARDS_WATER = ifelse(START.COUNT, NA, ESW_NN_DIST < lag(ESW_NN_DIST)),
         TOWARDS_GSW = ifelse(START.COUNT, NA, GSW_NN_DIST < lag(GSW_NN_DIST)),
         PERM_DIST_NORM = (ESWP_NN_DIST - m1) / m1 ,
         GSW_DIST_NORM = (GSW_NN_DIST - m2) / m2)

mod1 = lme4::glmer(SEASON ~ PERM_DIST_NORM + MFREQ + (1|ID),
    family = binomial(link="logit"),
    data = ele.norm)
mod2 = lme4::glmer(formula = SEASON ~ MFREQ + GSW_FREQ + (1|ANIMAL_ID),
                   family = binomial(link="logit"),
                   data = ele.norm)
mod3 = lme4::glmer(TOWARDS_WATER ~ SEASON * THIRST_ESW + (1|ID),
                   family = binomial(link="logit"),
                   data = ele.norm)
summary(mod1)
summary(mod2)
summary(mod3)
ggplot(ele.norm) +
  geom_histogram(aes(x=LOCAL_WATER_FREQ), bins=40) +
  facet_wrap(~YEAR)

# ******************************************************************************
#                               SEASONAL DISTANCES FROM WATER
# ******************************************************************************

sdat = ele.13 %>% st_drop_geometry() %>% 
  group_by(DATE, SEASON, SEX) %>% 
  summarize(mdist =  mymean(ESW_NN_DIST),
            mpdist = mean(ESWP_NN_DIST, na.rm=TRUE),
            mjdist = mean(GSW_NN_DIST, na.rm=TRUE),
            n=n(), 
            MFREQ=mean(MFREQ))

p1=ggplot(sdat) + 
  geom_bar(aes(x=DATE, y=MFREQ, fill=SEASON), stat='identity') + 
  scale_fill_brewer(palette="Dark2", direction=-1)+ 
  theme(axis.title.x=element_blank())
p2 = ggplot(sdat) + 
  geom_bar(aes(x=DATE, y=mdist, fill=SEASON), stat='identity') + 
  scale_fill_brewer(palette="Dark2", direction=-1)+ 
  theme(axis.title.x=element_blank())
p3 = ggplot(sdat) + 
  geom_bar(aes(x=DATE, y=mpdist, fill=SEASON), stat='identity') + 
  scale_fill_brewer(palette="Dark2", direction=-1)+ 
  theme(axis.title.x=element_blank())
p4 = ggplot(sdat) + 
  geom_bar(aes(x=DATE, y=mjdist, fill=SEASON), stat='identity') + 
  scale_fill_brewer(palette="Dark2", direction=-1)+ 
  theme(axis.title.x=element_blank())
p1 / p2 / p3 / p4


