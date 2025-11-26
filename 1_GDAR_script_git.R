# Libraries ------
library(adegenet)
library(dplyr)
library(tidyr)
library(sf)
library(hierfstat)
library(lme4)
library(ggplot2)
library(extrafont)
library(patchwork)
source('functions_git.R')

# 1. Genetic data -------------
meta <- read.csv('gdar_metadata.csv', h=T) %>% 
  drop_na(lon, lat)
ogwd <- getwd()


## Read in genotypes ------
# Set wd to folder with files
str <- 'str'
setwd(str)  

file_list<- as.list(list.files())

geno_list <- lapply(file_list, read_genot)

setwd(ogwd)

## Data summary: count species and datasets -----
genopops <- unlist(lapply(geno_list, function(x) unique(pop(x))))

popfilt <- meta %>% 
  filter(pop %in% c(genopops))

# Number of datasets
length(unique(popfilt$data_doi))

# Number of species
length(unique(popfilt$species))

# Median number of individuals and sites in each dataset
pf <- popfilt %>% 
  group_by(data_doi) %>% 
  summarise(total_inds = sum(num_individuals),
            total_sites = n())

median(pf$total_inds)
range(pf$total_inds)

median(pf$total_sites)
range(pf$total_sites)

## Estimate global FST and mean Hs-----
fst <- Map(function(geno, file){
  spmeta <- meta %>% 
    dplyr::filter(pop %in% pop(geno))
  splitty <- seppop(geno)
  popsinmetadata <- repool(splitty[spmeta$pop])
  print(file)
  bs <- basic.stats(popsinmetadata)
  fstdf <- data.frame(filename = file,
                      fst = bs$overall['Fst'],
                      hs = bs$overall['Hs'],
                      ht = bs$overall['Ht'])
  rownames(fstdf) <- NULL
  return(fstdf)}, 
  geno = geno_list, file = unlist(file_list))

fstdf <- do.call('rbind', fst)


# 2. GDARs -----
## Nested area sampling --------
mollweide <- 'ESRI:54009'
empty_df <- data.frame(allele_count = NA,
                       allelic_richness = NA,
                       He = NA,
                       individuals = NA,
                       num_sites = NA,
                       area_m2 = NA)

area_dfl <- list()
for(i in 1:length(geno_list)){
  rads <- get_radii(geno_list[[i]], meta, c(0.1, 0.25, 0.5, 0.75, 0.9, 0.99), mollweide)
  if(sum(rads) != 0){
    try({area_samples <- ctr_nested_site_samp(microsat = geno_list[[i]], 
                                              radius = rads, 
                                              data = meta, eecrs = mollweide)[[1]]
    
    
    area_dfl[[i]] <- area_samples
    })
  } else{
    area_dfl[[i]] <- empty_df
  }
  print(i)
}

# save file
#saveRDS(area_dfl, 'area_dfl.rds')

# attach global FST
area_dffi <- Map(function(x, y){
  if(is.null(x)){
    x <- empty_df
  }
  x$filename <- y
  fstt <- fstdf %>% filter(filename == y)
  x$global_fst <- fstt$fst
  x$global_hs <- fstt$hs
  x$global_ht <- fstt$ht
  x$log_allelic_richness <- log10(x$allelic_richness)
  x$log_gene_diversity <- log10(x$He)
  x$log_allele_count <- log10(x$allele_count)
  x$log_area <- log10(x$area_m2)
  x <- drop_na(x, He)
  x <- x %>% rename(gene_diversity = He)
  return(x)
}, x = area_dfl, y = file_list)

## Estimate z -----
# Z for allele counts
zAC <- lapply(area_dffi, function(x) get.z2(area_data = x, variable = x$log_allele_count))

# Z for allelic richness
zAR <- lapply(area_dffi, function(x) get.z2(area_data = x, variable = x$log_allelic_richness))

# Z for gene diversity
zGD <- lapply(area_dffi, function(x) get.z2(area_data = x, variable = x$log_gene_diversity))


z.attach <- Map(function(zac, zar, zgd, area_df){
  area_df$zAC <- zac$zGDAR
  area_df$zAR <- zar$zGDAR
  area_df$zGD <- zgd$zGDAR
  area_df$zACse <- zac$z_se
  area_df$zARse <- zar$z_se
  area_df$zGDse <- zgd$z_se
  area_df$r2_AC <- zac$rsq
  area_df$r2_AR <- zar$rsq
  area_df$r2_GD <- zgd$rsq
  area_df$c_AC <- zac$int_c
  area_df$c_AR <- zar$int_c
  area_df$c_GD <- zgd$int_c
  return(area_df)
}, zac=zAC, zar=zAR, zgd=zGD, area_df=area_dffi)

zfst <- do.call('rbind', z.attach)

# save file for interaction models
#saveRDS(zfst, 'zfst_long_se.rds')

# count species with NA z values:
#zsub<-zfst %>% 
#  distinct(filename, .keep_all = TRUE)
#
#nrow(zsub[which(is.na(zsub$zAC)),]) #0
#nrow(zsub[which(is.na(zsub$zGD)),]) #6
#nrow(zsub[which(is.na(zsub$zAR)),]) #4
#rm(zsub)

##
zfst <- zfst %>% 
  drop_na(zAC, zAR, zGD) %>% 
  group_by(filename) %>% 
  distinct(global_fst, zAC, zAR, zGD, .keep_all = T) %>% 
  ungroup() 

zfst1 <- merge(zfst, meta %>% select(species, class, filename) %>% distinct(), by = 'filename', all.x=TRUE, all.y = F)

# save file
#saveRDS(zfst1, 'zfst1_df.rds')

## z summary ------
z.df <- zfst1 %>% 
  select(class, global_fst, zAC, zAR, zGD) %>% 
  pivot_longer(c(zAC, zAR, zGD), names_to = 'variable',
               values_to = 'z')
z.dfse <- zfst1 %>% 
  select(class, global_fst, zACse, zARse, zGDse) %>% 
  pivot_longer(c(zACse, zARse, zGDse), names_to = 'variable',
               values_to = 'z_se')
z.df$z_se <- z.dfse$z_se
z.df$genetic_metric <- as.factor(z.df$variable)
levels(z.df$genetic_metric) <- c('allele count', 'allelic richness', 'gene diversity') # reorder for plot

z.df %>% 
  group_by(genetic_metric) %>% 
  summarise(mean_z = mean(z),
            sd_z =sd(z))
z.df %>% 
  group_by(genetic_metric, class) %>% 
  summarise(mean_z = mean(z),
            sd_z =sd(z))

r2df %>% 
  group_by(variable) %>% 
  summarise(mean_r2 = mean(r2),
            sd_r2 = sd(r2))
# 3. Plots -------

## Fig 1 example species -----
m_luci <- alldata %>% 
  filter(filename == 'Johnson_Myotis_lucifugus.gen')
h_platy <- alldata %>% 
  filter(filename=='Rovito_Hydromantes_platycephalus.gen')
Fig1.twosp<- ggplot() +
  geom_jitter(data = m_luci, aes(y = log_allele_count, x = log_area), color = '#229BA6', alpha = 0.5) +
  geom_smooth(data = m_luci, aes(y = log_allele_count, x = log_area), method = 'lm', color = '#229BA6', fill = '#229BA6') +
  geom_jitter(data = h_platy, aes(y = log_allele_count, x = log_area), color = '#392C57', alpha = 0.5) +
  geom_smooth(data = h_platy, aes(y = log_allele_count, x = log_area), method = 'lm', color = '#392C57', fill='#392C57') +
  labs(y = 'log genetic diversity', x = 'log area') +
  theme_classic(base_size = 20) +
  theme(text=element_text(family='Roboto Medium'),
        legend.position = 'none',
        axis.ticks=element_blank())

## Fig 2 -------
### bottom: Z across groups ----
pointsize <- 0.25
fatten <- 3

zac_class <- ggplot(z.df %>% filter(genetic_metric=='allele count'), aes(y = z, x = class)) +
  geom_violin(color = NA, alpha = 0.4, fill = '#012030') +
  geom_pointrange(aes(ymin = z-z_se, ymax = z+z_se), fatten = fatten, size =pointsize, position = position_jitter(width = 0.1)) +
  labs(y = expression(z[AC]), x = '') +
  theme_minimal(base_size = 11) +
  theme(text=element_text(family='Roboto Medium'))

zar_class <- ggplot(z.df %>% filter(genetic_metric=='allelic richness'), aes(y = z, x = class)) +
  geom_violin(color = NA, alpha = 0.4, fill = '#13678A') +
  geom_pointrange(aes(ymin = z-z_se, ymax = z+z_se), fatten = fatten, size =pointsize, position = position_jitter(width = 0.1)) +
  labs(y = expression(z[AR]), x = '') +
  theme_minimal(base_size = 11) +
  theme(text=element_text(family='Roboto Medium'))

zgd_class <- ggplot(z.df %>% filter(genetic_metric=='gene diversity'), aes(y = z, x = class)) +
  geom_violin(color = NA, alpha = 0.4, fill = '#45C4B0') +
  geom_pointrange(aes(ymin = z-z_se, ymax = z+z_se), fatten = fatten, size =pointsize, position = position_jitter(width = 0.1)) +
  labs(y = expression(z[GD]), x = '') +
  theme_minimal(base_size = 11) +
  theme(text=element_text(family='Roboto Medium'))

z_by_class <- zac_class + zar_class + zgd_class

### top:  log GD vs log area ------
z.attach1 <- lapply(z.attach, function(x){
  
  areass <- unique(x$log_area)
  area_perc <- c(0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  x <- x %>% 
    mutate(perc_area = case_when(log_area == areass[1] ~ area_perc[1],
                                 log_area == areass[2] ~ area_perc[2],
                                 log_area == areass[3] ~ area_perc[3],
                                 log_area == areass[4] ~ area_perc[4],
                                 log_area == areass[5] ~ area_perc[5],
                                 log_area == areass[6] ~ area_perc[6]))
}
)

percac_plots <- lapply(z.attach1, perc_loglog_plot, 'log_allele_count', 'zAC')


alldata <- do.call('rbind', z.attach1)

## overall effect
modac <- lmer(log_allele_count ~ log_area + (log_area|filename), data=alldata)
modar <- lmer(log_allelic_richness ~ log_area + (log_area|filename), data=alldata)
modgd <- lmer(log_gene_diversity ~ log_area + (log_area||filename), data=alldata)


alldata$pred_ac <- predict(modac,re.form=NA)  ## population level
alldata$pred_ar <- predict(modar,re.form=NA)  ## population level
alldata$pred_gd <- predict(modgd,re.form=NA)  ## population level


lines1 <- ggplot() +
  geom_smooth(data = alldata, aes(y = log_allele_count, x = log_area, shape = filename), 
              color = '#acc2ce', method = 'lm', fill = NA, lwd = 0.4) +
  geom_smooth(data = alldata, aes(y = pred_ac, x = log_area), method = 'lm', color = '#012030', fill='#012030') +
  labs(y = 'log allele count', x = 'log area') +
  theme_minimal() +
  theme(text=element_text(family='Roboto Medium'),
        legend.position = 'none')

lines2 <- ggplot() +
  geom_smooth(data = alldata, aes(y = log_allelic_richness, x = log_area, shape = filename), 
              color = '#baddec', method = 'lm', fill = NA, lwd = 0.4) +
  geom_smooth(data = alldata, aes(y = pred_ar, x = log_area), method = 'lm', color = '#13678A', fill='#13678A') +
  labs(y = 'log allelic richness', x = 'log area') +
  theme_minimal() +
  theme(text=element_text(family='Roboto Medium'),
        legend.position = 'none')

lines3 <- ggplot() +
  geom_smooth(data = alldata, aes(y = log_gene_diversity, x = log_area, shape = filename), 
              color = '#b0f4ea', method = 'lm', fill = NA, lwd = 0.4) +
  geom_smooth(data = alldata, aes(y = pred_gd, x = log_area), method = 'lm', color = '#45C4B0', fill='#45C4B0') +
  labs(y = 'log gene diversity', x = 'log area') +
  theme_minimal() +
  theme(text=element_text(family='Roboto Medium'),
        legend.position = 'none')

linee <- lines1+lines2+lines3

Fig2 <- linee/z_by_class +  plot_annotation(tag_levels = 'a')


## Fig S1: Plot R2 vs. fst -------- 
r2df <- zfst1 %>% 
  select(species, class, r2_AC, r2_AR, r2_GD) %>% 
  pivot_longer(c(r2_AC, r2_AR, r2_GD), names_to = 'variable',
               values_to = 'r2')
r2df$genetic_metric <- as.factor(r2df$variable)
levels(r2df$genetic_metric) <- c('allele count', 'allelic richness', 'gene diversity')

ora<- c('#012030', '#13678A', '#45C4B0')

(r2plot <- ggplot(r2df, aes(y = r2, x = genetic_metric, color = genetic_metric)) +
    geom_boxplot(notch = TRUE, outlier.alpha = 0, lwd = 1) + 
    geom_jitter(width = 0.1) +
    scale_color_manual(values = ora) +
    labs(y = expression(R^2), x = '') +
    theme_minimal(base_size = 14) +
    theme(text=element_text(family="Roboto Medium"), 
          legend.position = 'none')
)

## Fig. S2 (log log plots) -----

llgdac_plots <- lapply(z.attach, loglog_plot, 'log_allele_count', 'zAC')
llgdar_plots <- lapply(z.attach, loglog_plot, 'log_allelic_richness', 'zAR')
llgdgd_plots <- lapply(z.attach, loglog_plot, 'log_gene_diversity', 'zGD')

## Fig. S5 -------
# Number of loci & alleles
nallo_df <- data.frame(filename = unlist(file_list),
                       n_loci = unlist(lapply(geno_list, nLoc)),
                       n_alleles = unlist(lapply(geno_list, function(b) sum(b@loc.n.all))))

zfst2 <- merge(zfst1, nallo_df, by = 'filename', all=FALSE)

p1<-ggplot(data=zfst2, aes(y=r2_AC, x = n_alleles)) +
  geom_point() +
  geom_smooth(method='lm', color = '#012030', fill = '#012030', alpha=0.25) +
  labs(title = "", y = expression(R^2), x = '# alleles') +
  theme_classic()

p2<-ggplot(data=zfst2, aes(y=r2_AR, x = n_alleles)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_smooth(method='lm', color = '#13678A', fill = '#13678A', alpha=0.25) +
  labs(title = "", y = expression(R^2), x = '# alleles') +
  theme_classic()

p3<-ggplot(data=zfst2, aes(y=r2_GD, x = n_alleles)) +
  geom_point() +
  geom_smooth(method='lm', color = '#45C4B0', fill = '#45C4B0', alpha=0.25) +
  labs(title = "", y = expression(R^2), x = '# alleles') +
  theme_classic()

p4<-ggplot(data=zfst2, aes(y=r2_AC, x = n_loci)) +
  geom_point() +
  geom_smooth(method='lm', color = '#012030', fill = '#012030', alpha=0.25) +
  labs(title = "Allele count", y = expression(R^2), x = '# loci') +
  theme_classic()

p5<-ggplot(data=zfst2, aes(y=r2_AR, x = n_loci)) +
  geom_point() +
  geom_smooth(method='lm') +
  geom_smooth(method='lm', color = '#13678A', fill = '#13678A', alpha=0.25) +
  labs(title = "Allelic richness", y = expression(R^2), x = '# loci') +
  theme_classic()

p6<-ggplot(data=zfst2, aes(y=r2_GD, x = n_loci)) +
  geom_point() +
  geom_smooth(method='lm', color = '#45C4B0', fill = '#45C4B0', alpha=0.25) +
  labs(title = "Gene diversity", y = expression(R^2), x = '# loci') +
  theme_classic()

FigS5 <- (p4/p5/p6) | (p1/p2/p3)