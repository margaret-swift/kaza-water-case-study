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
#                           THEMES
############################################################

# https://emilhvitfeldt.github.io/r-color-palettes/discrete/awtools/a_palette/
colors.ptd = c( "#019875FF", "#315c4d","#99B898FF", "#FECEA8FF", "#FF847CFF")
colors.ys = c('#fa3628', '#ff7a70', '#ffaca6',
             '#c4a8f0', "#a071eb", '#5d27b3', '#4c2a80', '#2d026e')
colors.siz = c('#6988cf', '#c0cef0', 'black', '#346beb',
               '#24304f', '#b8c0d9')
plot.theme = theme_minimal() + theme(text=element_text(size=15))

############################################################
#                           FIGURE 4
############################################################

# Fig. 4a: Mean Annual Rainfall (MAR) for KAZA aggregated regionally
p.rainfall = ggplot(mardata,
                    aes(x=as.character(YEAR), y=mMAR, fill=YEAR)) + 
  geom_hline(yintercept=mean(mardata$mMAR), color='#35323b', linetype="dashed") +
  geom_bar(stat='identity') + 
  geom_point(aes(x=YEAR, y=mMAR-sdMAR)) +
  geom_point(aes(x=YEAR, y=mMAR+sdMAR)) +
  geom_segment(aes(y=mMAR-sdMAR, yend=mMAR+sdMAR), 
               color='#35323b') +
  ylab('mm MAR') + xlab('year')  +
  plot.theme + guides(fill="none") + 
  scale_fill_manual(values=colors.ys)
p.rainfall
ggsave(p.rainfall,
       file=here(outdir, 'plots', 'KAZA_MAR.png'),
       width=12, height=3)

# Fig. 4b: Hydroshed water fill levels
p.fill.hy = ggplot(data=stats %>% filter(TYPE == "HYDROSHED"),
                   aes(x=PERIOD, y=ESW_SIZE_M2/10000, fill=PERIOD)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~YEAR, nrow=1, scales="free_x") +
  ylab('ESW fill (km2)') + xlab('') +
  scale_fill_manual(values=colors) +
  plot.theme + guides(fill="none") + 
  theme(axis.text.x = element_blank())
p.fill.hy
ggsave(p.fill.hy,
       file=here(outdir, 'plots', 'KAZA_total_fill_HY.png'),
       width=10, height=4)

# Fig. 4c: SIZ water fill levels
buildPlot <- function(ids, colors) {
  p = ggplot(data=stats %>% filter(id %in% ids), 
             aes(x=START, y=ESW_SIZE_M2/10000, color=ID)) + 
    geom_point(size=0.5) + geom_line() +
    ylab('ESW fill (km2)') + xlab('') +
    plot.theme +
    scale_color_manual(values=colors) +
    scale_x_date(date_breaks="4 months", 
                 date_labels='%b') + 
    guides(color="none")
  p
}
p1 = buildPlot(c("Kariba", "Zambezi", "Okavango", "Makgadikgadi"), colors.siz[1:4])
p2 = buildPlot(c("Chobe", "Linyanti"), colors.siz[5:6])
p.fill.siz = (p1 + theme(axis.text.x = element_blank())) / 
              p2 + plot_layout(ncol=1, heights=c(3,1))
p.fill.siz
ggsave(p.fill.siz, 
       file=here(outdir, 'plots', 'KAZA_total_fill_SIZ.png'),
       width=10, height=5)


############################################################
#                           FIGURE 5
############################################################

# Figure 5a - fill anomaly vs period
sdat.hy = sdat%>% filter(TYPE == "HYDROSHED") %>% 
  mutate(YEAR = factor(YEAR, levels=mardata$YEAR))
stats.hy = stats %>% filter(COVER_P<1, TYPE == "HYDROSHED")%>% 
  mutate(YEAR = factor(YEAR, levels=mardata$YEAR))
p5a = ggplot(stats.hy, aes(x=PERIOD, group=YEAR, color=YEAR)) + 
  geom_hline(yintercept=0, linetype='dashed') +
  geom_jitter(aes(y=COVER_P), width=0.215, shape=1, alpha=0.7) +
  geom_line(data=sdat.hy, aes(y=cover_mu), linewidth=1) +
  geom_point(data=sdat.hy, aes(y=cover_mu), size=2.5, color="black") +
  geom_point(data=sdat.hy, aes(y=cover_mu), width=0.2) +
  xlab("pentad") + ylab('pMWF') +
  guides(color="none") +
  plot.theme +
  scale_color_manual(values=colors.ys)
p5a
ggsave(p5a,
       file=here(outdir, 'plots', 'Figure5_a.png'),
       width=5, height=5)

# Figure 5b - Beta model of fill vs precip by period
p5b = ggplot(data=stats %>% filter(COVER_P<1, COVER_P>0),
             aes(x=MAR_D, y=COVER_P, color=PERIOD)) + 
  geom_point(size=0.3) + 
  geom_smooth(method="lm") +
  geom_vline(xintercept=0, linetype='dashed') +
  xlab('dMAR') + 
  ylab('logit (pMWF)') +
  plot.theme +
  theme( strip.text.x = element_blank() ) + 
  guides(color=guide_legend(title='')) + 
  scale_color_manual(values=colors.ptd) + 
  guides(color="none")
p5b
ggsave(p5b,
       file=here(outdir, 'plots', 'Figure5_b.png'),
       width=5, height=5)

############################################################
#                      EXTRA FIGURES
############################################################

# Line graph of fill levels over time
ggplot(groups, aes(x=START, y=COVER_D, group=ID)) + 
  geom_point(size=0.5) + geom_line(alpha=0.3) +
  facet_wrap(~GROUP, ncol=1, scales="free_y") + 
  theme_minimal() + theme(text=element_text(size=15))
groups %>% 
  filter(START==as.Date('2019-02-01')) %>% 
  ggplot() + 
  geom_point(aes(y=mean_fill, x=sd_fill, color=GROUP)) + 
  geom_abline(slope=1, intercept=0)+ 
  theme_minimal() + theme(text=element_text(size=15))
