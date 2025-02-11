# READING -- Libraries ----------------------------------------------------
#Data mungering
library(tidyverse)
library(data.table)
library(broom)

select <- dplyr::select

genTheme <- function(x) {
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 25),
    legend.key = element_rect(fill = "transparent"),
    legend.key.size = unit(3, 'line'),
    legend.box= "vertical",
    plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
    # plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title = element_text(size = 30),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
    ) # bg of the plot
}

# READING -- data ---------------------------------------------------------
raw <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/kappaMortalitySimulations.csv')

phiRaw <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/phiData.csv')

grazingRaw<- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/grazingData.csv')

mumbyMortality <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/mumbyIncreasingMortality.csv')

ccaAdditions <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/ccaAdditionsFreeSpace.csv')
# biyearlyAdditions <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/biyearlyAdditions.csv')
yearlyAdditions <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/yearlyAdditions.csv')
yearlySubtractions <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/macroMinus.csv')
yearlyAddandSubs <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/yearlyCCAMacro.csv')

minCCAAdds <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/MinimumYearlyInterventions/yearlyAdditions.csv')
minMacroMinus <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/MinimumYearlyInterventions/macroMinus.csv')

minCCACoral <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/MinimumYearlyInterventions/yearlyCCACoral.csv')
minMacroCoral <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/MinimumYearlyInterventions/yearlyMacroCoral.csv')
minCoral <- read_csv('~/Documents/GitHub/ccaReefModeling/data/raw/MinimumYearlyInterventions/yearlyCoral.csv')

# Organism colors  --------------------------------------------------------
orgColors <-  c('#8C1AF5', '#A88353', '#6FA86C')

# VISUALIZATIONS -- Max coral cover ----------------------------------------------------------------
maxVis <- phiRaw%>%
  rename(mortality = xParam, kappa = yParam)%>%
  group_by(mortality, kappa, phi, g0)%>%
  summarize_if(is.numeric, max)%>%
  mutate(phi = round(phi, 2))%>%
  filter(phi == 0 | phi == 0.32 | phi == 0.63 | phi == 1)

coralDiff <- maxVis%>%
  ungroup()%>%
  mutate(kappa = round(kappa, 2))%>%
  filter(phi == 1,
         kappa == 0.05 | kappa == 1)%>%
  select(mortality, kappa, cEnd)%>%
  pivot_wider(names_from = 'kappa', values_from = 'cEnd')%>%
  mutate(difference = abs(`0.05` - `1`))
  

maxVisGrazing <- grazingRaw%>%
  rename(mortality = xParam, kappa = yParam)%>%
  group_by(mortality, kappa, phi, g0)%>%
  summarize_if(is.numeric, max)%>%
  mutate(g0 = round(g0,2))%>%
  filter(g0 == 0.10 | g0 == 0.18 | g0 == 0.29 | g0 == 0.31)

pdf("~/Documents/GitHub/ccaReefModeling/data/plots/kapMortMaxCoralSize.pdf", width = 20, height = 10)
maxVis%>%
  filter(g0 == 0.2,
         phi == 1)%>%
  mutate(type = ifelse(aEnd <= 0.5, 'square', 'circle'))%>%
  ggplot(aes(mortality, kappa, fill = cEnd)) +
  geom_raster() +
  geom_point(aes(color = aEnd, size = 3, shape = type)) +
  # scale_color_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
  scale_color_gradient(low = 'white', high = '#8C1AF5') +
  # scale_color_viridis_b() +
  scale_fill_gradient2(low = 'white', mid = 'grey', high = 'black', midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(color = 'Max CCA community size',
       fill = 'Max coral community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'right')

maxVis%>%
  filter(g0 == 0.2,
         phi == 1)%>%
  mutate(type = ifelse(aEnd <= 0.5, 'square', 'circle'))%>%
  ggplot(aes(mortality, kappa, fill = cEnd)) +
  geom_raster() +
  geom_point(aes(color = aEnd, size = 3)) +
  # scale_color_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
  scale_color_gradient(low = 'white', high = '#8C1AF5') +
  scale_fill_gradient2(low = 'white', mid = 'grey', high = 'black', midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(color = 'Max CCA community size',
       fill = 'Max coral community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 'right')

maxVis%>%
  filter(g0 == 0.2,
         phi == 1)%>%
  ggplot(aes(mortality, kappa, fill = aEnd)) +
  geom_raster() +
  scale_fill_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(fill = 'Max coral community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

dev.off()

## PLotting difference between 0 and 1 kappa
pdf('~/Documents/GitHub/ccaReefModeling/data/plots/kappaDifference.pdf', width = 15, height = 10)
coralDiff%>%
  ggplot(aes(mortality, difference, color = 'coral')) +
  geom_line(size = 3) +
  scale_color_manual(values = '#A88353') +
  labs(x = 'Coral Mortality (µ)', y = 'coral cover difference') +
  gen_theme()
dev.off()

pdf('~/Documents/GitHub/ccaReefModeling/data/plots/Phig0Shifts.pdf', width = 30, height = 10)
maxVis%>%
  filter(g0 == 0.2)%>%
  ggplot(aes(mortality, kappa, fill = cEnd)) +
  geom_raster() +
  scale_fill_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(fill = 'Max coral community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~reorder(phi, -phi), nrow = 1) +
  labs(title = 'shifting inhibition')
  

# maxVis%>%
#   filter(g0 == 0.2)%>%
#   ggplot(aes(mortality, kappa, fill = aEnd)) +
#   geom_raster() +
#   scale_fill_gradient2(midpoint = 0.5) +
#   # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
#   labs(fill = 'Max CCA community size',
#        y = 'Coral settlement preference (kappa)',
#        x = 'Coral Mortality') +
#   genTheme() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   facet_wrap(~phi) +
#   labs(title = 'shifting inhibition')

# maxVisGrazing%>%
#   filter(phi == 1)%>%
#   ggplot(aes(mortality, kappa, fill = cEnd)) +
#   geom_raster() +
#   scale_fill_gradient2(midpoint = 0.5) +
#   # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
#   labs(fill = 'Max coral community size',
#        y = 'Coral settlement preference (kappa)',
#        x = 'Coral Mortality') +
#   genTheme() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   facet_wrap(~round(g0, 2)) +
#   labs(title = 'shifting mortality')

maxVisGrazing%>%
  filter(phi == 1)%>%
  ggplot(aes(mortality, kappa, fill = aEnd)) +
  geom_raster() +
  scale_fill_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(fill = 'Max CCA community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~g0, nrow = 1) +
  labs(title = 'shifting mortality')

maxVisGrazing%>%
  filter(phi == 1)%>%
  ggplot(aes(mortality, kappa, fill = mEnd)) +
  geom_raster() +
  scale_fill_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(fill = 'Max Macroalgal community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~g0, nrow = 1) +
  labs(title = 'shifting mortality')

dev.off()
  


# VISUALIZATIONS -- Mumby mortality ---------------------------------------
maxMumby <- mumbyMortality%>%
  # rename(mortality = xParam, kappa = yParam)%>%
  group_by(mu)%>%
  summarize_if(is.numeric, max)

minMumby <- mumbyMortality%>%
  # rename(mortality = xParam, kappa = yParam)%>%
  group_by(mu)%>%
  summarize_if(is.numeric, min)

pdf('~/Documents/GitHub/ccaReefModeling/data/plots/mumbyMortality.pdf', width = 13, height = 10)
maxMumby%>%
  ggplot(aes(mu)) +
  geom_line(aes(y = cEnd), color = '#A88353', size = 5) +
  geom_line(data = minMumby, aes(y = mEnd),color = '#6FA86C', size = 5) +
  genTheme() 
dev.off()

# VIZUALIZATIONS -- settlement inhibition shift ---------------------------
shiftPhi <- phiRaw%>%
  rename(mortality = xParam, kappa = yParam)%>%
  mutate(kappa = round(kappa, 2))%>%
  filter(mortality == 0.11,
         kappa == 0.05)%>%
  group_by(phi)%>%
  summarize_if(is.numeric, max)%>%
  ungroup()%>%
  select(phi, cEnd, mEnd, aEnd)%>%
  gather(organism, abundance, 2:4)

pdf('~/Documents/GitHub/ccaReefModeling/data/plots/shiftingPhi.pdf', width = 15, height = 10)
shiftPhi%>%
  ggplot(aes(phi, abundance, color = organism)) +
  geom_line(aes(group = organism, linetype = organism), size = 3) +
  genTheme() +
  scale_color_manual(values = orgColors) +
  scale_linetype_manual(values = c('dashed', 'dashed', 'solid')) +
  scale_x_reverse() +
  labs(y = 'Relative abundance at equilibrium', x = 'Settlement inhibition')
dev.off()  

pdf('~/Documents/GitHub/ccaReefModeling/data/plots/shiftingPhiHeatmap.pdf', width = 15, height = 30)
maxVis%>%
  filter(g0 == 0.2,
         mortality > 0.07)%>%
  ggplot(aes(mortality, kappa, fill = cEnd)) +
  geom_raster() +
  geom_point(aes(color = aEnd, size = 3)) +
  # scale_color_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
  scale_color_gradient(low = 'white', high = '#8C1AF5') +
  scale_fill_gradient2(low = 'white', mid = 'grey', high = 'black', midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(color = 'Max CCA community size',
       fill = 'Max coral community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  facet_wrap(~-phi, nrow = 4) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(hjust = -0.25, angle = -45),
        legend.position = 'None')
dev.off()

pdf('~/Documents/GitHub/ccaReefModeling/data/plots/shiftingPhiHeatmapSupplement.pdf', width = 30, height = 10)
maxVis%>%
  filter(g0 == 0.2)%>%
         # mortality > 0.07)%>%
  ggplot(aes(mortality, kappa, fill = cEnd)) +
  geom_raster() +
  geom_point(aes(color = aEnd, size = 3)) +
  # scale_color_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
  scale_color_gradient(low = 'white', high = '#8C1AF5') +
  scale_fill_gradient2(low = 'white', mid = 'grey', high = 'black', midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(color = 'Max CCA community size',
       fill = 'Max coral community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  facet_wrap(~-phi, nrow = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(hjust = -0.25, angle = -45),
        legend.position = 'None')
dev.off()


# CCA addition  -------------------------------------------------------
ccaMax <- ccaAdditions%>%
  rename(mortality = xParam, kappa = yParam)%>%
  group_by(mortality, kappa, additions)%>%
  summarize_if(is.numeric, max)

ccaMax%>%
  # filter(g0 == 0.2)%>%
  ggplot(aes(mortality, kappa, fill = cEnd)) +
  geom_raster() +
  scale_fill_gradient2(midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(fill = 'Max coral community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~additions) +
  labs(title = 'shifting inhibition')

ccaMax%>%
  # filter(g0 == 0.2)%>%
  ggplot(aes(mortality, kappa, fill = mEnd)) +
  geom_raster() +
  scale_fill_gradient2(midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(fill = 'Max CCA community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~additions) +
  labs(title = 'shifting inhibition')


# Timed addition ----------------------------------------------------------
ccaTimedMax <- yearlyAdditions%>%
  rename(mortality = xParam, kappa = yParam)%>%
  group_by(mortality, kappa, additions)%>%
  summarize_if(is.numeric, max)%>%
  mutate(additionPercent = str_c(round(additions*100, 1), '%'))

facet_labels <- ccaTimedMax$additionPercent
names(facet_labels) <- ccaTimedMax$additions

yearlyAddPhasePlanes <- yearlyAdditions%>%
  rename(mortality = xParam, kappa = yParam)%>%
  mutate(kappa = round(kappa, 2),
         mortality = round(mortality, 2))%>%
  mutate(equilibrium = case_when(cEnd < 0.00001 ~ 'algal',
                                 TRUE ~ 'calcifier'))%>%
  group_by(additions, mortality, kappa, equilibrium)%>%
  summarize_if(is.numeric, max)%>%
  ungroup()%>%
  select(additions, mortality, equilibrium, kappa, cEnd, mEnd, aEnd)%>%
  gather(organism, abundance, 5:7)%>%
  filter(equilibrium == 'algal')%>%
  unite(eOrg, c(organism, equilibrium), sep = '_', remove = FALSE)

pdf('~/Documents/GitHub/ccaReefModeling/data/plots/yearlyAdditionsOrgs.pdf', width = 10, height = 10)
yearlyAddPhasePlanes%>%
  filter(mortality == 0.02,
         kappa == 0.38)%>%
  ggplot(aes((additions*100), (abundance*100), color = organism)) +
  geom_line(aes(group = eOrg, linetype = equilibrium), size = 5) +
  genTheme() +
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_color_manual(values = orgColors) +
  facet_wrap(~equilibrium, nrow = 2) +
  labs(y = 'Abundance at equilibrium (%)', x = 'Freespace converted to CCA each year (%)')


yearlyAddPhasePlanes%>%
  filter(mortality == 0.08,
         kappa == 0.38)%>%
  ggplot(aes((additions*100), (abundance*100), color = organism)) +
  geom_line(aes(group = eOrg, linetype = equilibrium), size = 5) +
  genTheme() +
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_color_manual(values = orgColors) +
  facet_wrap(~equilibrium, nrow = 2) +
  labs(y = 'Abundance at equilibrium (%)', x = 'Freespace converted to CCA each year (%)')


yearlyAddPhasePlanes%>%
  filter(mortality == 0.11,
         kappa == 0.38)%>%
  ggplot(aes((additions*100), (abundance*100), color = organism)) +
  geom_line(aes(group = eOrg, linetype = equilibrium), size = 5) +
  genTheme() +
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_color_manual(values = orgColors) +
  facet_wrap(~equilibrium, nrow = 2) +
  labs(y = 'Abundance at equilibrium (%)', x = 'Freespace converted to CCA each year (%)')


yearlyAddPhasePlanes%>%
  filter(mortality == 0.16,
         kappa == 0.38)%>%
  ggplot(aes((additions*100), (abundance*100), color = organism)) +
  geom_line(aes(group = eOrg, linetype = equilibrium), size = 5) +
  genTheme() +
  theme(legend.position = 'None',
       strip.background = element_blank(),
       strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_color_manual(values = orgColors) +
  facet_wrap(~equilibrium, nrow = 2) +
  labs(y = 'Abundance at equilibrium (%)', x = 'Freespace converted to CCA each year (%)')

dev.off()

# pdf('~/Documents/GitHub/ccaReefModeling/data/plots/yearlyAdditions.pdf', width = 15, height = 10)
# ccaTimedMax%>%
#   # filter(g0 == 0.2)%>%
#   ggplot(aes(mortality, kappa, fill = cEnd)) +
#   geom_raster() +
#   scale_fill_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
#   # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
#   labs(fill = 'Max coral community size',
#        y = 'Coral settlement preference (kappa)',
#        x = 'Coral Mortality') +
#   genTheme() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   facet_wrap(~additions, labeller = labeller(additions = facet_labels)) +
#   labs(title = 'shifting inhibition')
# 
# ccaTimedMax%>%
#   # filter(g0 == 0.2)%>%
#   ggplot(aes(mortality, kappa, fill = aEnd)) +
#   geom_raster() +
#   scale_fill_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
#   # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
#   labs(fill = 'Max CCA community size',
#        y = 'Coral settlement preference (kappa)',
#        x = 'Coral Mortality') +
#   genTheme() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   facet_wrap(~additions, labeller = labeller(additions = facet_labels)) +
#   labs(title = 'shifting inhibition')
# 
# 
# pdf('~/Downloads/ccaAddedyearly.pdf', width = 15, height = 10)
# ccaTimedMax%>%
#   # filter(g0 == 0.2)%>%
#   ggplot(aes(mortality, kappa, fill = mEnd)) +
#   geom_raster() +
#   scale_fill_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
#   # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
#   labs(fill = 'Max macroalgae community size',
#        y = 'Coral settlement preference (kappa)',
#        x = 'Coral Mortality') +
#   genTheme() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   facet_wrap(~additions, labeller = labeller(additions = facet_labels)) +
#   labs(title = 'Changing Amount of CCA added')
# dev.off()


# yearly subtractions -----------------------------------------------------
yearlySubPhasePlanes <- yearlySubtractions%>%
  rename(mortality = xParam, kappa = yParam)%>%
  mutate(kappa = round(kappa, 2),
         mortality = round(mortality, 2))%>%
  mutate(equilibrium = case_when(cEnd < 0.00001 ~ 'algal',
                                 TRUE ~ 'calcifier'))%>%
  group_by(additions, mortality, kappa, equilibrium)%>%
  summarize_if(is.numeric, max)%>%
  ungroup()%>%
  select(additions, mortality, equilibrium, kappa, cEnd, mEnd, aEnd)%>%
  gather(organism, abundance, 5:7)%>%
  filter(equilibrium == 'algal')%>%
  unite(eOrg, c(organism, equilibrium), sep = '_', remove = FALSE)

pdf('~/Documents/GitHub/ccaReefModeling/data/plots/yearlySubtractions.pdf', width = 10, height = 10)
yearlySubPhasePlanes%>%
  filter(mortality == 0.02,
         kappa == 0.38)%>%
  ggplot(aes((additions*100), (abundance*100), color = organism)) +
  geom_line(aes(group = eOrg, linetype = equilibrium), size = 5) +
  genTheme() +
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_color_manual(values = orgColors) +
  facet_wrap(~equilibrium, nrow = 2) +
  # scale_linetype_manual(values = c('solid', 'dashed')) +
  labs(y = 'Abundance at equilibrium (%)', x = 'Macroalgae removed each year (%)')

yearlySubPhasePlanes%>%
  filter(mortality == 0.08,
         kappa == 0.38)%>%
  ggplot(aes((additions*100), (abundance*100), color = organism)) +
  geom_line(aes(group = eOrg, linetype = equilibrium), size = 5) +
  genTheme() +
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_color_manual(values = orgColors) +
  facet_wrap(~equilibrium, nrow = 2) +
  # scale_linetype_manual(values = c('solid', 'dashed')) +
  labs(y = 'Abundance at equilibrium (%)', x = 'Macroalgae removed each year (%)')

yearlySubPhasePlanes%>%
  filter(mortality == 0.11,
         kappa == 0.38)%>%
  ggplot(aes((additions*100), (abundance*100), color = organism)) +
  geom_line(aes(group = eOrg, linetype = equilibrium), size = 5) +
  genTheme() +
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_color_manual(values = orgColors) +
  facet_wrap(~equilibrium, nrow = 2) +
  # scale_linetype_manual(values = c('solid', 'dashed')) +
  labs(y = 'Abundance at equilibrium (%)', x = 'Macroalgae removed each year (%)')

yearlySubPhasePlanes%>%
  filter(mortality == 0.16,
         kappa == 0.38)%>%
  ggplot(aes((additions*100), (abundance*100), color = organism)) +
  geom_line(aes(group = eOrg, linetype = equilibrium), size = 5) +
  genTheme() +
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_color_manual(values = orgColors) +
  facet_wrap(~equilibrium, nrow = 2) +
  # scale_linetype_manual(values = c('solid', 'dashed')) +
  labs(y = 'Abundance at equilibrium (%)', x = 'Macroalgae removed each year (%)')

dev.off()

subtractionMax <- yearlySubtractions%>%
  rename(mortality = xParam, kappa = yParam)%>%
  group_by(mortality, kappa, additions)%>%
  summarize_if(is.numeric, max)%>%
  mutate(additionPercent = str_c(round(additions*100, 1), '%'))

subfacet_labels <- subtractionMax$additionPercent
names(subfacet_labels) <- subtractionMax$additions

subtractionMax%>%
  # filter(g0 == 0.2)%>%
  ggplot(aes(mortality, kappa, fill = cEnd)) +
  geom_raster() +
  scale_fill_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(fill = 'Max coral community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~additions, labeller = labeller(additions = subfacet_labels)) +
  labs(title = 'shifting inhibition')

subtractionMax%>%
  # filter(g0 == 0.2)%>%
  ggplot(aes(mortality, kappa, fill = mEnd)) +
  geom_raster() +
  scale_fill_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(fill = 'Max coral community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~additions, labeller = labeller(additions = subfacet_labels)) +
  labs(title = 'shifting inhibition')


subtractionMax%>%
  # filter(g0 == 0.2)%>%
  ggplot(aes(mortality, kappa, fill = aEnd)) +
  geom_raster() +
  scale_fill_gradient2(low = '#78B7C5', mid = 'white', high = '#FD6467', midpoint = 0.5) +
  # scale_color_gradientn(colors= c('red', 'white', 'blue')) +
  labs(fill = 'Max coral community size',
       y = 'Coral settlement preference (kappa)',
       x = 'Coral Mortality') +
  genTheme() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~additions, labeller = labeller(additions = subfacet_labels)) +
  labs(title = 'shifting inhibition')


# yearlySubs and Adds -----------------------------------------------------
yearlySubAddPhase <- yearlyAddandSubs%>%
  rename(mortality = xParam, kappa = yParam)%>%
  mutate(kappa = round(kappa, 2),
         mortality = round(mortality, 2))%>%
  mutate(equilibrium = case_when(cEnd < 0.00001 ~ 'algal',
                                 TRUE ~ 'calcifier'))%>%
  group_by(additions, mortality, kappa, equilibrium)%>%
  summarize_if(is.numeric, max)%>%
  ungroup()%>%
  select(additions, mortality, equilibrium, kappa, cEnd, mEnd, aEnd)%>%
  gather(organism, abundance, 5:7)%>%
  filter(equilibrium == 'algal')%>%
  unite(eOrg, c(organism, equilibrium), sep = '_', remove = FALSE)

pdf('~/Documents/GitHub/ccaReefModeling/data/plots/yearlySubtractionsAdditions.pdf', width = 10, height = 10)
yearlySubAddPhase%>%
  filter(mortality == 0.02,
         kappa == 0.38)%>%
  ggplot(aes((additions*100), (abundance*100), color = organism)) +
  geom_line(aes(group = eOrg, linetype = equilibrium), size = 5) +
  genTheme() +
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_color_manual(values = orgColors) +
  facet_wrap(~equilibrium, nrow = 2) +
  # scale_linetype_manual(values = c('solid', 'dashed')) +
  labs(y = 'Abundance at equilibrium (%)', x = 'Macroalgae removed and CCA added each year (%)')

yearlySubAddPhase%>%
  filter(mortality == 0.08,
         kappa == 0.38)%>%
  ggplot(aes((additions*100), (abundance*100), color = organism)) +
  geom_line(aes(group = eOrg, linetype = equilibrium), size = 5) +
  genTheme() +
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_color_manual(values = orgColors) +
  facet_wrap(~equilibrium, nrow = 2) +
  # scale_linetype_manual(values = c('solid', 'dashed')) +
  labs(y = 'Abundance at equilibrium (%)', x = 'Macroalgae removed and CCA added each year (%)')

yearlySubAddPhase%>%
  filter(mortality == 0.11,
         kappa == 0.38)%>%
  ggplot(aes((additions*100), (abundance*100), color = organism)) +
  geom_line(aes(group = eOrg, linetype = equilibrium), size = 5) +
  genTheme() +
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_color_manual(values = orgColors) +
  facet_wrap(~equilibrium, nrow = 2) +
  # scale_linetype_manual(values = c('solid', 'dashed')) +
  labs(y = 'Abundance at equilibrium (%)', x = 'Macroalgae removed and CCA added each year (%)')

yearlySubAddPhase%>%
  filter(mortality == 0.16,
         kappa == 0.38)%>%
  ggplot(aes((additions*100), (abundance*100), color = organism)) +
  geom_line(aes(group = eOrg, linetype = equilibrium), size = 5) +
  genTheme() +
  theme(legend.position = 'None',
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_color_manual(values = orgColors) +
  facet_wrap(~equilibrium, nrow = 2) +
  # scale_linetype_manual(values = c('solid', 'dashed')) +
  labs(y = 'Abundance at equilibrium (%)', x = 'Macroalgae removed and CCA added each year (%)')

dev.off()




# minimum Additions for no bistability ------------------------------------
MaxMinCCA <- minCCAAdds%>%
  group_by(additions)%>%
  summarize_if(is.numeric, max)

MaxMinMacro <- minMacroMinus%>%
  group_by(additions)%>%
  summarize_if(is.numeric, max)

MaxMinCCACoral <- minCCACoral%>%
  group_by(additions, ccaAdds)%>%
  summarize_if(is.numeric, max)

MaxMinMacroCoral <- minMacroCoral%>%
  group_by(additions, ccaAdds)%>%
  summarize_if(is.numeric, max)

MaxMinCoral <- minCoral%>%
  group_by(ccaAdds)%>%
  summarize_if(is.numeric, max)


pdf('~/Documents/GitHub/ccaReefModeling/data/plots/InterventionTimeLines.pdf', width = 15, height = 10)
MaxMinCCA%>%
  ggplot(aes(additions*100, mEnd)) +
  geom_line(size = 3, color = '#6FA86C') +
  labs(y = 'Benthic cover of macroalage',
       x = '% of free space converted to CCA per year') +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

MaxMinMacro%>%
  ggplot(aes(additions*100, mEnd)) +
  geom_line(size = 3, color = '#6FA86C') +
  labs(y = 'Benthic cover of macroalage',
       x = '% of macroalgae removed per year') +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

MaxMinCoral%>%
  ggplot(aes(ccaAdds*100, mEnd)) +
  geom_line(size = 3, color = '#6FA86C') +
  labs(y = 'Benthic cover of macroalage',
       x = '% of free space converted to coral per year') +
  scale_x_continuous(limits = c(0,2.5)) +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

MaxMinCCACoral%>%
  ggplot(aes(ccaAdds*100, additions*100, fill = mEnd)) +
  geom_raster() +
  scale_fill_gradient(low = 'white', high = '#6FA86C') +
  labs(y = '% of free space converted to CCA per year',
       x = '% of free space converted to coral per year') +
  scale_x_continuous(limits = c(0,1)) +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

MaxMinMacroCoral%>%
  ggplot(aes(ccaAdds*100, additions*100, fill = mEnd)) +
  geom_raster() +
  scale_fill_gradient(low = 'white', high = '#6FA86C') +
  labs(y = '% of macroalgae removed per year',
       x = '% of free space converted to coral per year') +
  scale_x_continuous(limits = c(0,1)) +
  genTheme() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

