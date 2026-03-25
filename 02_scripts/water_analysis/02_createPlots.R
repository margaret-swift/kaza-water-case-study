# createPlots.R
# Margaret Swift <margaret.swift@cornell.edu>
# 
# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
source(here::here('02_scripts', 'utilities.R'))
i_am('02_scripts/water_analysis/02_createPlots.R')
outdir <- here('03_output', 'water_analysis')
load(here(outdir, "hydrostats.rdata"))
pacman::p_load(LaplacesDemon) # logit and invlogit

############################################################
############################################################
# SANDBOX
x=stats %>% st_drop_geometry() %>% 
  filter(TYPE == "HYDROSHED") %>% 
  group_by(PERIOD, YEAR) %>% 
  summarize(fill = sum(ESW_SIZE))
hist(x$fill, 30)

############################################################

# Mean Annual Rainfall (MAR) for KAZA aggregated regionally
mardata = mar_m %>% group_by(YEAR) %>% 
  summarize(mMAR=mean(MAR), sdMAR=sd(MAR))
p.rainfall = ggplot(mardata, aes(x=factor(YEAR), y=mMAR)) + 
  geom_bar(stat='identity') + 
  geom_segment(aes(y=mMAR-sdMAR, yend=mMAR+sdMAR), color='darkgray') +
  theme_classic() + ylab('mm MAR') + xlab('year') +
  theme(text=element_text(size=15))
p.rainfall
ggsave(p.rainfall, 
       file=here(outdir, 'plots', 'KAZA_MAR.png'),
       width=10, height=3)

############################################################
# Plotting total water fill over time 

# Hydroshed water fill levels
colors = RColorBrewer::brewer.pal(6, 'BrBG')[c(5,6,3,2,1)]
p.fill.hy = ggplot() + 
  geom_bar(data=stats %>% filter(TYPE == "HYDROSHED"),
           aes(x=PERIOD, y=ESW_SIZE_M2, fill=PERIOD),
           position="stack", stat="identity") +
  facet_wrap(~YEAR, nrow=1, scales="free_x") +
  theme_classic() + ylab('ESW fill (km2)') + 
  xlab('time') +
  theme(text=element_text(size=15)) +
  scale_fill_manual(values=colors)
ggsave(p.fill.hy, 
       file=here(outdir, 'plots', 'KAZA_total_fill_HY.png'),
       width=10, height=6)

# SIZ water fill levels
p.fill.siz = ggplot(data=stats %>% filter(TYPE == "SIZ"),
       aes(x=START, y=ESW_SIZE_M2, 
           color=ID)) + 
  geom_point() + geom_line() +
  theme_classic() + ylab('ESW fill (km2)') + 
  xlab('time') +
  theme(text=element_text(size=15)) + 
  scale_x_date(date_breaks="3 months", 
               date_labels='%b')
ggsave(p.fill.siz, 
       file=here(outdir, 'plots', 'KAZA_total_fill_SIZ.png'),
       width=10, height=6)

############################################################
# Checking whether data are normally distributed
ggplot(betadata) + 
  geom_histogram(aes(x=COVER_P, 
                     y=..density..), bins=40,color='darkgray', fill='lightgray') +
  geom_density(aes(x=COVER_P), color='black', linewidth=1.2) + 
  xlab('proportion of maximum fill') +
  theme_classic() +
  facet_wrap(~TYPE) +
  theme(text=element_text(size=15))


############################################################
### PLOTTING SEASONALITY BY GROUP ###
newstats = stats %>%   
  mutate(FGROUP = case_when( mean_fill > 0.8 ~ 1,
                            mean_fill > 0.25 ~ 2, 
                            .default = 3))
ggplot(newstats %>% filter(TYPE == "HYDROSHED", sd_fill<0.2) %>% 
         mutate(ID = as.numeric(ID)),
       aes(x=START, y=COVER_P, group=ID, color=ID)) + 
  geom_point() + geom_line() + 
  facet_wrap(~FGROUP, ncol=1)

# histograms of % fill by hydroshed type
ggplot(newstats) + 
  geom_histogram(aes(x=COVER_P, y=..density..), bins=20, color='white') + 
  geom_density(aes(x=COVER_P), color='black', linewidth=1.2) + 
  facet_wrap(~FGROUP, ncol=1)

# plotting KAZA regions by different coefficients
mar_avg = mar_m %>% 
  group_by(ID) %>% 
  summarize(MMAR = mean(MAR))
fill_sf = fill_stats %>% 
  left_join(hydro_shp, by="ID") %>% 
  left_join(mar_avg, by="ID") %>% 
  st_as_sf()
ggplot(fill_sf) + geom_sf(aes(fill=MMAR * mean_slope)) +
  scale_fill_distiller(palette="BrBG", direction=-1)
ggplot(fill_sf) + geom_sf(aes(fill=mean_fill)) +
  scale_fill_distiller(palette="BrBG", direction=-1)
ggplot(fill_sf) + geom_sf(aes(fill=mean_slope)) + 
  scale_fill_distiller(palette="BrBG", direction=-1)

# Figure 4a - fill anomaly vs period
p4a = ggplot(stats, aes(x=PERIOD, group=YEAR, color=YEAR)) + 
  geom_hline(yintercept=0, linetype='dashed') +
  geom_jitter(aes(y=COVER_P), width=0.215, shape=1, alpha=0.7) +
  geom_line(data=sdat, aes(y=cover_mu)) +
  geom_point(data=sdat, aes(y=cover_mu), size=2.5, color="black") +
  geom_point(data=sdat, aes(y=cover_mu), width=0.2) +
  xlab("period") + ylab('% of maximum fill') +
  facet_wrap(~TYPE, nrow=2) +
  theme_classic() +
  theme(text=element_text(size=15)) + 
  guides(color="none") +
  theme( strip.text.x = element_blank() )

# Figure 4b - Beta model of fill vs precip by period
p4b = ggplot(data=betadata, aes(x=MAR_D, y=logit(COVER_P), color=TYPE)) + 
  geom_point(size=0.3) + 
  geom_smooth(method="lm") +
  geom_vline(xintercept=0, linetype='dashed') +
  facet_wrap(~PERIOD, nrow=1) + 
  theme_classic() +
  xlab('normalized MAR') + 
  ylab('logit (pMWF)') +
  theme(text=element_text(size=15)) +
  theme( strip.text.x = element_blank() ) + 
  guides(color=guide_legend(title='')) +
  scale_color_manual(values=c('darkgray', 'blue'))
p4b
p4a + p4b

