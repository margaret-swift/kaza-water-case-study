# createPlots.R
# Margaret Swift <margaret.swift@cornell.edu>
# 
# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************
source(here::here('scripts', 'utilities.R'))
i_am('02_scripts/water_analysis/02_createPlots.R')
outdir <- here('03_output', 'water_analysis')
load(here(outdir, "hydrostats.rdata"))

############################################################

mardata = mar_m %>% group_by(YEAR) %>% 
  summarize(m=mean(MAR), sd=sd(MAR))
ggplot(mardata, aes(x=factor(YEAR), y=m)) + 
  geom_bar(stat='identity') + 
  geom_segment(aes(y=m-sd, yend=m+sd), color='darkgray') +
  theme_classic() + ylab('mm MAR') + xlab('year') +
  theme(text=element_text(size=15))

hist(hdata$ESW_SIZE)

x=stats %>% st_drop_geometry() %>% 
  filter(TYPE == "HYDROSHED") %>% 
  group_by(PERIOD, YEAR) %>% 
  summarize(fill = sum(ESW_SIZE))
hist(x$fill, 30)
TukeyHSD(aov(lm(fill~PERIOD, data=x)))

# diagnostic plots - cover size
colors = RColorBrewer::brewer.pal(6, 'BrBG')[c(5,6,3,2,1)]
ggplot() + 
  geom_bar(data=stats %>% filter(TYPE == "HYDROSHED"),
           aes(x=PERIOD, y=ESW_SIZE / 10000, fill=PERIOD),
           position="stack", stat="identity") +
  facet_wrap(~YEAR, nrow=1, scales="free_x") +
  theme_classic() + ylab('ESW fill (km2)') + 
  xlab('time') +
  theme(text=element_text(size=15)) +
  scale_fill_manual(values=colors)

# SIZ water fill levels
ggplot(data=stats %>% filter(TYPE == "SIZ"),
       aes(x=START, y=ESW_SIZE / 10000, 
           color=ID)) + 
  geom_point() + geom_line() +
  theme_classic() + ylab('ESW fill (ha)') + 
  xlab('time') +
  theme(text=element_text(size=15)) + 
  scale_x_date(date_breaks="3 months", 
               date_labels='%b')


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
newstats %>% ggplot() + 
  geom_histogram(aes(x=(COVER_P)), bins=20, 
                 color='white') + 
  facet_wrap(~FGROUP, ncol=1)

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

## BY YEAR ##
ggplot(stats%>% filter(TYPE=="HYDROSHED"), 
       aes(x=PERIOD, y=ESW_SIZE, 
           group=YEAR, fill=YEAR)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(text=element_text(size=15))

# plotting SIZ catchment rainfall over time
siz_catch_mar %>% 
  reshape2::melt(id.var='name') %>% 
  ggplot() + geom_bar(aes(x=variable, y=value, fill=name), 
                      stat='identity', position='dodge') +
  scale_fill_manual(values=c('magenta', 'red', 'orange')) + 
  theme_classic()

# plotting SIZ fill over time
begin = "-10-15"; end="-04-30";
years = 2019:2025
starts = as.Date(paste0(years-1, begin))
ends = as.Date(paste0(years, end))
date.df = data.frame(START=starts, END=ends, YEAR=years)
my.df = reshape2::melt(siz_catch_mar, id.var="name", 
               variable.name="YEAR", value.name="MAR") %>% 
  mutate(YEAR = as.numeric(gsub("mm", "", YEAR))) %>% 
  rename(ID=name) %>% 
  left_join(date.df, by="YEAR")

hdata %>% filter(TYPE == "SIZ") %>% mutate(START=as.POSIXct(START)) %>% 
  ggplot() + 
  geom_rect(data = my.df, 
            aes(xmin = START, xmax = END, 
                ymin = -Inf, ymax = MAR*50000), 
            fill='blue', alpha = 0.2)  +
  geom_point(aes(x=START, y=ESW_SIZE)) + 
  geom_line(aes(x=START, y=ESW_SIZE)) + 
  facet_wrap(~ID, ncol=1, scales="free_y")+
  scale_x_datetime(date_breaks="4 months") + 
  scale_color_manual(values=c('red', 'orange')) + 
  theme_classic()

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

ggplot(betadata) + 
  geom_histogram(aes(x=COVER_P, 
                     y=..density..), bins=40,color='darkgray', fill='lightgray') +
  geom_density(aes(x=COVER_P), color='black', linewidth=1.2) + 
  xlab('proportion of maximum fill') +
  theme_classic() +
  facet_wrap(~TYPE) +
  theme(text=element_text(size=15))
  

