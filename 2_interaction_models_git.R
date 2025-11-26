# Libraries ------
library(lme4)
library(lmerTest)
library(effects)
library(performance)
library(ggplot2)
library(extrafont)
library(patchwork)
library(viridis)

# 1. Data ----
zfst <- readRDS('zfst_long_se.rds') # from GDAR_script

# 2. Analysis -----
## Gene diversity -----
imod_gd <- lmer(log_gene_diversity ~ log_area + global_fst + log_area*global_fst + (1|filename), data = zfst)
summary(imod_gd)

amod_gd <- lmer(log_gene_diversity ~ log_area + (1|filename), data = zfst)
summary(amod_gd)

## Allelic richness ------

imod_ar <- lmer(log_allelic_richness ~ log_area + global_fst + log_area*global_fst + (1|filename), data = zfst)
summary(imod_ar)

amod_ar <- lmer(log_allelic_richness ~ log_area + (1|filename), data = zfst)
summary(amod_ar)

## Allele count -----
imod_ac <- lmer(log_allele_count ~ log_area + global_fst + log_area*global_fst + (1|filename), data = zfst)
summary(imod_ac)

amod_ac <- lmer(log_allele_count ~ log_area + (1|filename), data = zfst)
summary(amod_ac)

# 3. Plots ------
## GD ------
ee_gd <- Effect(c('log_area', 'global_fst'), imod_gd)

eedf_gd <- as.data.frame(ee_gd)
eedf_gd$global_fst_fac <- as.factor(eedf_gd$global_fst)

blues <- c('#392C57', '#36688B', '#229BA6', '#27A69B', '#3ABF9D')

int_gd <- ggplot(eedf_gd,
                 aes(log_area, fit, color = global_fst_fac, fill = global_fst_fac))+
  geom_line() +
  geom_ribbon(colour=NA, 
              alpha=0.1,
              aes(ymin=lower,ymax=upper)) +
  scale_color_manual(values=blues) +
  scale_fill_manual(values=blues) +
  labs(color = expression(global~F[ST]), fill = expression(global~F[ST]), x = 'log area', title = 'Gene diversity') +
  theme_minimal() +
  theme(text=element_text(family='Roboto Medium'))

## AR -----
ee_ar <- Effect(c('log_area', 'global_fst'), imod_ar)

eedf_ar <- as.data.frame(ee_ar)
eedf_ar$global_fst_fac <- as.factor(eedf_ar$global_fst)

int_ar <- ggplot(eedf_ar,
       aes(log_area, fit, color = global_fst_fac, fill = global_fst_fac))+
  geom_line() +
  geom_ribbon(colour=NA,
              alpha=0.1,
              aes(ymin=lower,ymax=upper)) +
  scale_color_manual(values=blues) +
  scale_fill_manual(values=blues) +
  labs(color = expression(global~F[ST]), fill = expression(global~F[ST]), x = 'log area', title = 'Allelic richness') +
  theme_minimal() +
  theme(text=element_text(family='Roboto Medium'))

## AC ------
ee_ac <- Effect(c('log_area', 'global_fst'), imod_ac)

eedf_ac <- as.data.frame(ee_ac)
eedf_ac$global_fst_fac <- as.factor(eedf_ac$global_fst)

int_ac <- ggplot(eedf_ac,
                 aes(log_area, fit, color = global_fst_fac, fill = global_fst_fac))+
  geom_line() +
  geom_ribbon(colour=NA,
              alpha=0.1,
              aes(ymin=lower,ymax=upper)) +
  scale_color_manual(values=blues) +
  scale_fill_manual(values=blues) +
  labs(color = expression(global~F[ST]), fill = expression(global~F[ST]), x = 'log area', title = 'Allele count') +
  theme_minimal() +
  theme(text=element_text(family='Roboto Medium'))


## Fig. 3 ------
intplot <- int_ac + int_ar + int_gd + plot_layout(guides = 'collect')

zplotl <- plot_common(legend = T)

Fig3 <- zplotl/intplot
