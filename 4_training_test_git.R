# Libraries ------
# tidy
library(tidyr)
library(dplyr)

# Models
library(caret)
library(betareg)

# Plots
library(ggplot2)
library(ggrepel)
library(viridis)
library(patchwork)
library(extrafont)
source('functions_git.R')

# 1. Data -----
mam_dat <- read.csv('mammal_traits.csv', header = TRUE)
mam_dat <- mam_dat %>%
  mutate(log_mass = log10(mass_g),
         log_range = log10(range_area_km2),
         log_HR = log10(mean_HR_km2),
         log_area = log10(area_km2))

mpg_test_data <- read.csv('mpg_test_data.csv', header = TRUE) %>% 
  filter(!species_ref=='Peromyscus_leucopus_R870')

# 2. Analysis -----
## Original data ----
### FST ----
mod <- betareg(global_fst ~ log_HR + log_range + log_mass + log_area, data = mam_dat)
summary(mod)

### z ----
modAC <- lm(zAC ~ log_HR + log_range + log_mass + log_area, data = mam_dat)
modAR <- lm(zAR ~ log_HR + log_range + log_mass + log_area, data = mam_dat)
modGD <- lm(zGD ~ log_HR + log_range + log_mass + log_area, data = mam_dat)

summary(modAC)
summary(modAR)
summary(modGD)

## Train vs test data -----
### FST -----
trainingsamples <-createDataPartition(mam_dat$global_fst, times= 1000, p=0.8)

train.data <- lapply(trainingsamples, function(x) mam_dat[x,])
test.data <- lapply(trainingsamples, function(x) mam_dat[-x,])

train_models <- lapply(train.data, function(x) betareg(global_fst ~ log_HR + log_range + log_mass + log_area, data = x))

model_predictions <- Map(function(modellist, testlist) predict(modellist, testlist), 
                         modellist = train_models, testlist = test.data)

error_rates <- Map(function(predictions, testdat) RMSE(predictions, testdat$global_fst)/mean(testdat$global_fst), # CV
                   predictions = model_predictions, testdat = test.data)

# mean CV
mean(unlist(error_rates))

fstcv <- data.frame(z = 'FST', 
                    cv = unlist(error_rates))

### z ----
# zGD
train_modelszgd <- lapply(train.data, function(x) lm(zGD ~ log_HR + log_range + log_mass + log_area, data = x))

model_predictionszgd <- Map(function(modellist, testlist) predict(modellist, testlist), 
                          modellist = train_modelszgd, testlist = test.data)

error_rateszgd <- Map(function(predictions, testdat) RMSE(predictions, testdat$zGD)/mean(testdat$zGD), 
                    predictions = model_predictionszgd, testdat = test.data)

hist(unlist(error_rateszgd)) 
abline(v = mean(unlist(error_rateszgd)), col = 'red', lty='dashed', lwd = 2)

# Mean CV
mean(unlist(error_rateszgd))

zgdcv <- data.frame(z = 'GD', 
                    cv = unlist(error_rateszgd))

# zAR
train_modelszar <- lapply(train.data, function(x) lm(zAR ~ log_HR + log_range + log_mass + log_area, data = x))

model_predictionszar <- Map(function(modellist, testlist) predict(modellist, testlist), 
                          modellist = train_modelszar, testlist = test.data)

error_rateszar <- Map(function(predictions, testdat) RMSE(predictions, testdat$zAR)/mean(testdat$zAR), 
                    predictions = model_predictionszar, testdat = test.data)

hist(unlist(error_rateszar)) 
abline(v = mean(unlist(error_rateszar)), col = 'red', lty='dashed', lwd = 2)

# Mean CV
mean(unlist(error_rateszar))

zarcv <- data.frame(z = 'AR', 
                    cv = unlist(error_rateszar))

# zAC
train_modelszac <- lapply(train.data, function(x) lm(zAC ~ log_HR + log_range + log_mass + log_area, data = x))

model_predictionszac <- Map(function(modellist, testlist) predict(modellist, testlist), 
                          modellist = train_modelszac, testlist = test.data)

error_rateszac <- Map(function(predictions, testdat) RMSE(predictions, testdat$zAC)/mean(testdat$zAC), 
                    predictions = model_predictionszac, testdat = test.data)

hist(unlist(error_rateszac)) 
abline(v = mean(unlist(error_rateszac)), col = 'red', lty='dashed', lwd = 2)

# Mean CV
mean(unlist(error_rateszac))

zaccv <- data.frame(z = 'AC', 
                    cv = unlist(error_rateszac))

## Testing with MacroPopGen ------
# Model
moder <- betareg(global_fst ~ log_HR + log_range + log_mass + log_area, data = mam_dat)

# Coefficient of variation
# predicted values
mpg_pred <- predict(moder, mpg_test_data)
(mpg_mean_error <- RMSE(mpg_pred, mpg_test_data$global_fst)/mean(mpg_test_data$global_fst))

# Predict mean FST of species:
mpg_test_data_means <- mpg_test_data %>% 
  group_by(species) %>% 
  mutate(mean_FST = mean(global_fst))

mpg_test_data_means$global_fst <- mpg_test_data_means$mean_FST
mpg_pred2 <- predict(moder, mpg_test_data_means)
RMSE(mpg_pred2, mpg_test_data_means$global_fst)/mean(mpg_test_data_means$global_fst)



# 3. Plots ------
## Fig. 4 ------
pal <- c('#CCC354', '#012030', '#13678A', '#45C4B0') 

cvplotlist <- Map(cvplot, cvdf=list(fstcv, zaccv, zarcv, zgdcv), colo=pal, xaxislab=list(F,F,T,T))
cvplotlist[[1]] <- cvplotlist[[1]] + geom_vline(xintercept = mpg_mean_error*100, color = 'black', lty = 'dashed')

cv_fig4 <- (cvplotlist[[1]] + cvplotlist[[2]])/ (cvplotlist[[3]] + cvplotlist[[4]]) +plot_annotation(tag_levels = 'A')

## S3: Compare FST across datasets ------
fst_compare <- merge(mpg_test_data %>% select(species, global_fst), mam_dat %>% select(species, global_fst), by = 'species', all=FALSE)
fst_compare$species <- gsub('_', ' ', fst_compare$species)

fst_compare1 <- fst_compare %>% 
  group_by(species) %>% 
  summarise(mpg_mean = mean(global_fst.x),
            mmpop_mean = mean(global_fst.y)) %>% 
  ungroup()

figs3 <- ggplot(data = fst_compare1, aes(y = mpg_mean, x = mmpop_mean, label = species)) +
  geom_abline(slope=1, intercept = 0, color = 'red', lty = 'dashed', lwd = 0.5) +
  geom_point(size = 2) +
  geom_text(hjust = "inward", vjust = 1.75, size = 3, fontface = "italic") +
  labs(x = expression("F"[ST]~"(this paper)"), y=expression("F"[ST]~"(MacroPopGen)")) +
  theme_minimal() +
  theme(text=element_text(family='Roboto Medium'))


## S4: MPG FST predictions -----
pointsize=3
mpg_predA <- ggplot(data = mpg_test_data, aes(y = global_fst, x = mpg_pred)) +
  geom_abline(slope=1, intercept = 0, color = 'red', lty = 'dashed', lwd = 1) +
  geom_point(aes(color=species), size = pointsize) +
  labs(y = expression('observed F'[ST]), x=expression('predicted F'[ST]),
       title = 'Study-level predictions') +
  scale_color_viridis(discrete = T) +
  guides(color="none") +
  theme_minimal() +
  theme(text=element_text(family='Roboto Medium'))

mpg_predB <- ggplot(data = mpg_test_data_means, aes(y = global_fst, x = mpg_pred2)) +
  geom_abline(slope=1, intercept = 0, color = 'red', lty = 'dashed', lwd = 1) +
  geom_point(size = pointsize) +
  labs(y = expression('observed F'[ST]), x=expression('predicted F'[ST]),
       title = 'Species mean predictions') +
  theme_minimal() +
  theme(text=element_text(family='Roboto Medium'))

figs4 <- mpg_predA + mpg_predB
