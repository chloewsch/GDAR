# Libraries -----
library(tidyr)
library(dplyr)

# Plots
library(ggplot2)
library(extrafont)
library(patchwork)


# 1. Data ------
zfst <- readRDS('zfst_long_se.rds') # from GDAR_script

# 2. Observed vs Predicted loss -----
# for the largest area loss (maximum area vs minimum area per nested sample): 
# Predict from GDAR 

# Split into a list of dataframes based on file
zfst_list <- zfst %>% 
  group_by(filename) %>% 
  group_split()

# Create a grouping factor that identifies east nested sample site within a study
zfst3 <- lapply(zfst_list, function(x){
  n <- nrow(x)/6
  x$group <- rep(rep(c(1:n), each = 6))
  x$group <- paste0('g', x$group)
  return(x)
})

# Combine back into 1 dataframe
zfst3 <- do.call('rbind', zfst3)

# Keep only maximum and minimum area samples
zhi <- zfst3 %>% 
  group_by(filename, group) %>% 
  slice_max(area_m2) %>% 
  mutate(value = 'max')


zlo <- zfst3 %>% 
  group_by(filename, group) %>% 
  slice_min(area_m2) %>% 
  mutate(value = 'min')

zhilo <- rbind(zhi, zlo)

# Mean z values for each metric:
mean_zgd <- 0.02
mean_zar <- 0.03
mean_zac <- 0.13

# Calculate observed and predicted genetic diversity for each metric
zhilo2 <- zhilo %>% 
  group_by(filename, group) %>% 
  summarise(area = range(area_m2),
            gd = range(gene_diversity),
            ar = range(allelic_richness),
            ac = range(allele_count),
            zgd = unique(zGD),
            zar = unique(zAR),
            zac = unique(zAC)) %>%
  # Within a nested site group use the minimum or the max value for area and genetic diversity
  mutate(a1_a0 = min(area)/max(area), # Note: the ratio of min to max area is always 0.01 because we sampled up to 99% of the area
  # Observed proportion of diversity remaining = minimum/maximum
         obs_g1_g0 = min(gd)/max(gd), 
         obs_ar1_ar0 = min(ar)/max(ar),
         obs_ac1_ac0 = min(ac)/max(ac),
         obs_perc_lossgd = (1-obs_g1_g0)*100, # Observed percentage of genetic diversity lost = (1 - proportion remaining)*100
         obs_perc_lossar = (1-obs_ar1_ar0)*100,
         obs_perc_lossac = (1-obs_ac1_ac0)*100) %>% 
  # Predicted proportion of diversity remaining with STUDY z = (min area/max area)^z (estimated z for each study)
  mutate(predict_ratio_gd = a1_a0^zgd, 
         predict_ratio_ar = a1_a0^zar,
         predict_ratio_ac = a1_a0^zac,
         pred_perc_lossgd = (1-predict_ratio_gd)*100, # Predicted percentage diversity lost = (1-prop remaining)*100
         pred_perc_lossar = (1-predict_ratio_ar)*100,
         pred_perc_lossac = (1-predict_ratio_ac)*100) %>% 
  # Predicted proportion of diversity remaining with INTERMEDIATE z = (min area/max area)^z (mean z for each metric)
  # Note because the proportion of area loss is always the same, the predicted losses will be the same for all studies
  mutate(predict_ratio_gd_INT = a1_a0^mean_zgd, # mean zGD
         predict_ratio_ar_INT = a1_a0^mean_zar, # mean zAR
         predict_ratio_ac_INT = a1_a0^mean_zac, # mean zAC
         pred_perc_lossgd_INT = (1-predict_ratio_gd_INT)*100,
         pred_perc_lossar_INT = (1-predict_ratio_ar_INT)*100,
         pred_perc_lossac_INT = (1-predict_ratio_ac_INT)*100)


# 3. Mean observed vs predicted loss ---------
# Estimate mean observed loss for each study (average % loss across nested sample sites)
lossdf <- zhilo2 %>% group_by(filename) %>% 
  summarise(# Mean observed diversity lost
            mean_perc_lost_ogd = mean(obs_perc_lossgd, na.omit = TRUE),
            mean_perc_lost_oar = mean(obs_perc_lossar, na.omit = TRUE),
            mean_perc_lost_oac = mean(obs_perc_lossac, na.omit = TRUE),
            # Min and max observed loss
            min_perc_lost_ogd = min(obs_perc_lossgd, na.omit = TRUE),
            min_perc_lost_oar = min(obs_perc_lossar, na.omit = TRUE),
            min_perc_lost_oac = min(obs_perc_lossac, na.omit = TRUE),
            max_perc_lost_ogd = max(obs_perc_lossgd, na.omit = TRUE),
            max_perc_lost_oar = max(obs_perc_lossar, na.omit = TRUE),
            max_perc_lost_oac = max(obs_perc_lossac, na.omit = TRUE),
            # Predicted with study z
            perc_lost_pred_gd = unique(pred_perc_lossgd),
            perc_lost_pred_ar = unique(pred_perc_lossar),
            perc_lost_pred_ac = unique(pred_perc_lossac),
            # Predicted with intermediate z
            perc_lost_pred_gdi = unique(pred_perc_lossgd_INT),
            perc_lost_pred_ari = unique(pred_perc_lossar_INT),
            perc_lost_pred_aci = unique(pred_perc_lossac_INT))

# Create long form dataframe to plot all metrics on 1 plot
# Observed
loss_longo <- lossdf %>% 
  select(filename, mean_perc_lost_ogd, mean_perc_lost_oar, mean_perc_lost_oac) %>% 
  pivot_longer(c(mean_perc_lost_ogd, mean_perc_lost_oar, mean_perc_lost_oac), names_to = 'obsvariable',
               values_to = 'observed_loss') 

# Predicted (study z)
loss_longp <- lossdf %>% 
  select(filename, perc_lost_pred_gd, perc_lost_pred_ar, perc_lost_pred_ac) %>%  
  pivot_longer(c(perc_lost_pred_gd, perc_lost_pred_ar, perc_lost_pred_ac), names_to = 'predvariable',
               values_to = 'predicted_loss')

# Predicted (mean z)
loss_longpi <- lossdf %>% 
  select(filename, perc_lost_pred_gdi, perc_lost_pred_ari, perc_lost_pred_aci) %>%  
  pivot_longer(c(perc_lost_pred_gdi, perc_lost_pred_ari, perc_lost_pred_aci), names_to = 'predvariablei',
               values_to = 'predicted_loss_int')

# Put back into 1 dataframe
loss_long <- cbind.data.frame(loss_longo, loss_longp[,-1])
loss_long <- cbind.data.frame(loss_long, loss_longpi[,-1])

# Add a column for metric
loss_long$variable <- NA
loss_long$variable[grep('gd', loss_long$obsvariable)] <-'gene diversity'
loss_long$variable[grep('ar', loss_long$obsvariable)] <-'allelic richness'
loss_long$variable[grep('ac', loss_long$obsvariable)] <-'allele count'

# 4. Plots -------

## 4A: Observed vs predicted with mean z ----
Fig4a <- ggplot(loss_long, aes(shape = variable)) +
  geom_ribbon(aes(x = predicted_loss_int, y = observed_loss, xmin = observed_loss*0.5, xmax = observed_loss*1.5), fill = NA, color = 'grey65', linetype=3) +
  geom_ribbon(aes(x = predicted_loss_int, y = observed_loss, xmin = observed_loss*0.75, xmax = observed_loss*1.25), fill = NA, color = 'grey40', linetype=3) +
  geom_ribbon(aes(x = predicted_loss_int, y = observed_loss, xmin = observed_loss*0.9, xmax = observed_loss*1.1), fill = "grey85", color = NA) +
  geom_abline(slope = 1, intercept = 0, color = 'red', lty = 'dashed') +
  ggdist::stat_interval(aes(x = round(predicted_loss_int,2), y = observed_loss)) +
  geom_point(aes(x = predicted_loss_int, y = observed_loss), alpha =0.6, size=0.7) +
  scale_color_brewer() +
  labs(x = 'GDAR predicted % lost', y = 'mean observed % lost', shape = '', color = '', 
       title = 'A    Mean z') +
  theme_minimal() +
  theme(text=element_text(family="Roboto Medium"))


## 4B: study level z -----
gd_study_level <- ggplot(lossdf) +
  geom_ribbon(aes(x = perc_lost_pred_gd, y = mean_perc_lost_ogd, xmin = mean_perc_lost_ogd*0.5, xmax = mean_perc_lost_ogd*1.5), fill = NA, color = 'grey65', linetype=3) +
  geom_ribbon(aes(x = perc_lost_pred_gd, y = mean_perc_lost_ogd, xmin = mean_perc_lost_ogd*0.75, xmax = mean_perc_lost_ogd*1.25), fill = NA, color = 'grey40', linetype=3) +
  geom_ribbon(aes(x = perc_lost_pred_gd, y = mean_perc_lost_ogd, xmin = mean_perc_lost_ogd*0.9, xmax = mean_perc_lost_ogd*1.1), fill = "grey85", color = NA) +
  geom_abline(slope = 1, intercept = 0, color = 'red', lty = 'dashed') +
  geom_linerange(aes(x = perc_lost_pred_gd, ymin = min_perc_lost_ogd, ymax = max_perc_lost_ogd), color = '#45C4B0', alpha = 0.25, size = 1) +
  geom_point(aes(x = perc_lost_pred_gd, y = mean_perc_lost_ogd), color = '#45C4B0') +
  guides(color="none") +
  labs(x = 'GDAR predicted % lost', y = 'observed % loss') +
  theme_minimal() +
  theme(text=element_text(family="Roboto Medium"))

ar_study_level <- ggplot(lossdf) +
  geom_ribbon(aes(x = perc_lost_pred_ar, y = mean_perc_lost_oar, xmin = mean_perc_lost_oar*0.5, xmax = mean_perc_lost_oar*1.5), fill = NA, color = 'grey65', linetype=3) +
  geom_ribbon(aes(x = perc_lost_pred_ar, y = mean_perc_lost_oar, xmin = mean_perc_lost_oar*0.75, xmax = mean_perc_lost_oar*1.25), fill = NA, color = 'grey40', linetype=3) +
  geom_ribbon(aes(x = perc_lost_pred_ar, y = mean_perc_lost_oar, xmin = mean_perc_lost_oar*0.9, xmax = mean_perc_lost_oar*1.1), fill = "grey85", color = NA) +
  geom_abline(slope = 1, intercept = 0, color = 'red', lty = 'dashed') +
  geom_linerange(aes(x = perc_lost_pred_ar, ymin = min_perc_lost_oar, ymax = max_perc_lost_oar), color = '#13678A', alpha = 0.2, size = 1) +
  geom_point(aes(x = perc_lost_pred_ar, y = mean_perc_lost_oar), color = '#13678A') +
  guides(color="none") +
  labs(x = 'GDAR predicted % lost', y = 'observed % loss') +
  theme_minimal() +
  theme(text=element_text(family="Roboto Medium"))

ac_study_level <- ggplot(lossdf) +
  geom_ribbon(aes(x = perc_lost_pred_ac, y = mean_perc_lost_oac, xmin = mean_perc_lost_oac*0.5, xmax = mean_perc_lost_oac*1.5), fill = NA, color = 'grey65', linetype=3) +
  geom_ribbon(aes(x = perc_lost_pred_ac, y = mean_perc_lost_oac, xmin = mean_perc_lost_oac*0.75, xmax = mean_perc_lost_oac*1.25), fill = NA, color = 'grey40', linetype=3) +
  geom_ribbon(aes(x = perc_lost_pred_ac, y = mean_perc_lost_oac, xmin = mean_perc_lost_oac*0.9, xmax = mean_perc_lost_oac*1.1), fill = "grey85", color = NA) +
  geom_abline(slope = 1, intercept = 0, color = 'red', lty = 'dashed') +
  geom_linerange(aes(x = perc_lost_pred_ac, ymin = min_perc_lost_oac, ymax = max_perc_lost_oac), color = '#012030', alpha = 0.2, size = 1) +
  geom_point(aes(x = perc_lost_pred_ac, y = mean_perc_lost_oac), color = '#012030') +
  guides(color="none") +
  labs(x = 'GDAR predicted % lost', y = 'observed % loss', title = 'B    Study-specific z') +
  theme_minimal() +
  theme(text=element_text(family="Roboto Medium"))



## 4C: Prediction residuals ------
resi <- lossdf %>% 
        mutate(#  study z error
               res_gd = mean_perc_lost_ogd - perc_lost_pred_gd,
               res_ar = mean_perc_lost_oar - perc_lost_pred_ar,
               res_ac = mean_perc_lost_oac - perc_lost_pred_ac,
               # intermediate z error
               res_gdi = mean_perc_lost_ogd - perc_lost_pred_gdi,
               res_ari = mean_perc_lost_oar - perc_lost_pred_ari,
               res_aci = mean_perc_lost_oac - perc_lost_pred_aci,
               # over vs predict 
               under_predgd = ifelse(res_gd>0, 'under', 'over'),
               under_predar = ifelse(res_ar>0, 'under', 'over'),
               under_predac = ifelse(res_ac>0, 'under', 'over'),
               under_predgdi = ifelse(res_gdi>0, 'under', 'over'),
               under_predari = ifelse(res_ari>0, 'under', 'over'),
               under_predaci = ifelse(res_aci>0, 'under', 'over'))

table(resi$under_predac)/sum(table(resi$under_predac))
table(resi$under_predar)/sum(table(resi$under_predar))
table(resi$under_predgd)/sum(table(resi$under_predgd))

table(resi$under_predaci)/sum(table(resi$under_predaci))
table(resi$under_predari)/sum(table(resi$under_predari))
table(resi$under_predgdi)/sum(table(resi$under_predgdi))

## plot
resgd <- ggplot(data= resi) +
  geom_histogram(aes(x = res_gdi), bins = 10, alpha = 0.35, fill = '#2b786c', color = NA) +
  geom_histogram(aes(x = res_gd), bins = 10, alpha = 0.75, fill = '#45C4B0') +
  geom_vline(xintercept = 0, lty = 'dashed', lwd = 0.05) +
  labs(x = 'gene diversity') +
  theme_minimal() +
  theme(text=element_text(family="Roboto Medium"))

resar <- ggplot(data= resi) +
  geom_histogram(aes(x = res_ari), bins = 10, alpha = 0.35, fill = '#0c3e52', color = NA) +
  geom_histogram(aes(x = res_ar), bins = 10, alpha = 0.75, fill = '#13678A') +
  geom_vline(xintercept = 0, lty = 'dashed', lwd = 0.05) +
  labs(x = 'allelic richness') +
  theme_minimal() +
  theme(text=element_text(family="Roboto Medium"))

resac <- ggplot(data= resi) +
  geom_histogram(aes(x = res_aci), bins = 10, alpha = 0.35, fill = '#000608', color = NA) +
  geom_histogram(aes(x = res_ac), bins = 10, alpha = 0.75, fill = '#012030') +
  geom_vline(xintercept = 0, lty = 'dashed', lwd = 0.05) +
  labs(title ='C    Prediction error (observed - predicted % loss)' , x = 'allele count') +
  theme_minimal() +
  theme(text=element_text(family="Roboto Medium"))

## Fig 4 ----
Fig4 <- mean_loss_pred_ob_ploti/(ac_study_level+ar_study_level+gd_study_level)/(resac+resar+resgd)