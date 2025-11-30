# Libraries ------
# Traits
library(traitdata)
library(HomeRange)

# Taxonomy
library(taxize)

# tidy
library(tidyr)
library(dplyr)

# spatial
library(terra)

# Plots
library(ggplot2)
library(ggrepel)
library(viridis)
library(patchwork)
source('functions_git.R')

load('GDAR_env_figs.Rdata')

# 1. Compile traits ------
## Mass -----
data("pantheria")

traits <- pantheria %>% 
          select(Genus, Species, AdultBodyMass_g) %>% 
          mutate(Genus_species = paste(Genus, Species, sep = '_'))

# Mammal species names in our data
mammals <- zfst1 %>% 
  filter(class == 'mammal')

mammal_sp <- unique(mammals$species)

# Find missing species or mismatched names 
setdiff(mammal_sp, unique(traits$Genus_species))

# Update fisher scientific name
traits$Genus_species[which(traits$Genus_species=='Martes_pennanti')] <- 'Pekania_pennanti'

# Add mass information for Escalera's bat
traits <- rbind(traits, c('Myotis', 'escalerai', 7.25, 'Myotis_escalerai'))

# Filter and merge
massg <- traits %>% 
            #drop_na(AdultBodyMass_g) %>% 
            filter(Genus_species %in% mammal_sp) %>% 
            select(-c(Genus, Species)) %>% 
            rename(species = Genus_species)

# Missing data for Myotis septentrionalis and Sorex antinorii
# https://www.gbif.org/species/2432436
# https://www.gbif.org/species/2435998

massg$AdultBodyMass_g[which(is.na(massg$AdultBodyMass_g))] <- c(7.5, 9)

## Range area --------
## IUCN global terrestrial mammal ranges: https://www.iucnredlist.org/resources/spatial-data-download
## METADATA: 
# Presence: 1 = Extant; 2 = Probably extant; 3 = Possibly extant; 4 = Possibly extinct; 5 = Extinct; 6 = Uncertain
# Origin: 1 = Native; 2 = Reintroduced; 3 = Introduced; 4 = Vagrant; 5 = Uncertain; 6 = Assisted colonization
# Seasonal: 1 = Resident; 2 = Breeding Season; 3 = Non-breeding season; 4 = Passage; 5 = Uncertain
mam <- read_sf('TERRESTRIAL_MAMMALS.shp')

## Filter range data

## Create column to match with genetic data
# (Replace spaces with underscores)
mam$species <- gsub(' ', '_', mam$binomial)

# Match species across datasets:
species <- match(mammal_sp, mam$species)

# Find species that didn't match range data:
NAsp <- which(is.na(species))
mammal_sp[NAsp]

## Replace Martes pennanti with Pekania pennanti in both datasets:
mam$species <- gsub('Martes_pennanti', 'Pekania_pennanti', mam$species)

## Check again for NAs:
#which(is.na(match(mammal_sp, mam$species)))

species <- match(mammal_sp, mam$species)

ranges <- mam %>% 
  filter(species %in% c(mammal_sp)) %>% 
  filter(presence == 1,
         origin == 1, 
         seasonal == 1)

# area
range_vect <- vect(ranges)

rareas <- expanse(range_vect, unit = 'km', transform = TRUE)

ranges$area_km <- rareas

sp_areas <- ranges %>% 
  st_drop_geometry() %>% 
  group_by(species) %>%
  mutate(range_area_km2 = sum(area_km)) %>% 
  select(-area_km) %>%
  select(species, range_area_km2) %>% 
  distinct()


## Home ranges -----
homerange <- GetHomeRangeData()

homerange <- homerange %>% 
              filter(HR_Level == 'Individual', Life_Stage == 'Adult', Context == 'Wild') %>% 
              group_by(Species) %>% 
              summarise(mean_HR = mean(Home_Range_km2, na.rm = TRUE))

homerange$species <- gsub(' ', '_', homerange$Species)

# Check Pekania pennanti
homerange[which(homerange$species=='Martes_pennanti'),]
homerange$species[which(homerange$species=='Martes_pennanti')] <- 'Pekania_pennanti'

homerange$genus <- unlist(lapply(strsplit(homerange$Species, ' '), function(x) x[1]))

# Find missing species
homerange1 <- homerange %>% 
                filter(species %in% mammal_sp) %>% 
                select(genus, species, mean_HR)
missing <- setdiff(mammal_sp, homerange1$species)

homerange[which(homerange$genus=='Myotis'), ]
homerange[which(homerange$genus=='Tamiasciurus'), ]
homerange[which(homerange$genus=='Microtus'), ]
homerange[which(homerange$genus=='Rousettus'), ] # None
homerange[which(homerange$genus=='Vicugna'), ] # None
homerange[which(homerange$genus=='Lemmus'), ]

missing_genus <- unlist(lapply(strsplit(setdiff(mammal_sp, homerange1$species), '_'), function(x)x[1]))

missing_g_s <- data.frame(species = missing, 
                                genus = missing_genus)

missing_gdf <- homerange %>% 
                  filter(genus %in% missing_genus) %>% 
                  group_by(genus) %>% 
                  summarise(mean_mean_HR = mean(mean_HR))


missing_repl <- merge(missing_g_s, missing_gdf, by = 'genus') %>% 
                rename(mean_HR = mean_mean_HR)

# Vicugna: Franklin, W. 1974. The social behavior of the vicuna. The Behaviour of Ungulates and its relation to management, Vol. 1: 477-487 # per individual
# R aegyptiacus: https://cdnsciencepub.com/doi/full/10.1139/z11-013?casa_token=ekrKJUH9gZAAAAAA%3AWFq_v5QJ1MpQ_rSPSX4ggo7yNDxDtSYR7YyPQkxB78d0YpkKVsBvtvNnWfq9CHkn1ngi5JVkNC7uwbE

missing2 <- data.frame(genus = c('Vicugna', 'Rousettus'), species = c('Vicugna_vicugna', 'Rousettus_aegyptiacus'), mean_HR = c(0.034, 5.16))

homerange2 <- rbind(homerange1, missing_repl)

homerange3 <- rbind(homerange2, missing2)

setdiff(mammal_sp, homerange3$species)

HR_km2 <- homerange3 %>% 
                  select(species, mean_HR)
## Taxonomy ------
# Get NCBI unique identifier (UID) for each species:
uids <- lapply(mammal_sp, function(x) get_uid(x, messages = FALSE)[1])

uid_strip <- unlist(unique(uids))

tax_df <- data.frame(uid = uid_strip,
                     species = mammal_sp)


# Check for missing names:
nrow((tax_df[is.na(tax_df$uid),])) # none
#missingtax <- (tax_df[is.na(tax_df$uid),])


# Assign classifications:
taxize_class <- classification(uids, db = "ncbi")

pulltax <- lapply(taxize_class, function(x) as.data.frame(x))

pulltax1 <- tibble(uid = names(pulltax), pulltax) %>% 
  unnest(cols = c(pulltax)) %>% 
  filter(rank %in% c("order","family","genus", "species")) %>% 
  select(-id) %>% 
  distinct() %>%
  spread(rank, name) %>% 
  mutate(species = gsub(" ", "_", species))


## Merge -----
mam_zfst <- Reduce(function(x,y) merge(x,y, by="species", all=TRUE, incomparables = "NA"), 
                        list(mammals, pulltax1, massg, sp_areas, HR_km2))

nrow(mammals) == nrow(mam_zfst)

mam_zfst <- mam_zfst %>% 
            rename(mass_g = AdultBodyMass_g,
                   mean_HR_km2 = mean_HR) %>% 
            mutate(mass_g = as.numeric(mass_g),
                   log_mass = log10(mass_g),
                   log_range = log10(range_area_km2),
                   log_HR = log10(mean_HR_km2),
                   )

# Check for NAs
mam_zfst[which(is.na(mam_zfst$mass_g)),]
mam_zfst[which(is.na(mam_zfst$mean_HR_km2)),]
mam_zfst[which(is.na(mam_zfst$range_area_km2)),]

## Area of sample extent --------
bbox_areas <- list()

fns <- unique(zfst$filename)
for(i in 1:length(fns)){
  species_coord <- meta %>% 
    dplyr::filter(filename == fns[[i]])
  
  sites <- st_as_sf(species_coord, coords = c("lon", "lat"), crs = 4326)
  
  sites_bbox <- st_bbox(sites)
  bboxpoly <- makepolygon(sites_bbox, fns[i])
  bboxt <- vect(bboxpoly)
  area <- expanse(bboxt, unit = 'km', transform = TRUE)
  bbox_areas[[i]] <- area
}

gdar_area <- data.frame(filename = fns,
                       area_km2 = unlist(bbox_areas))

mam_dat <- merge(mam_zfst, gdar_area, by = 'filename')

# Fix number of sites and individuals
sites_indivs <- meta %>% 
                  select(filename, num_individuals) %>% 
                  filter(filename %in% mam_dat$filename) %>% 
                  group_by(filename) %>% 
                  summarise(individuals = sum(num_individuals),
                            num_sites = n())

colnm <- names(mam_dat)
mam_dat2 <- merge(mam_dat %>% select(-c(individuals, num_sites)), sites_indivs, by = 'filename')
mam_dat <- mam_dat2 %>% select(c(colnm))

## Save file -----
#write.csv(mam_dat, 'mammal_traits.csv', row.names = F, quote = F)

# 2. Validation with MacroPopGen (Lawrence et al. 2019)-----
# MacroPopGen database: https://figshare.com/articles/dataset/MacroPopGen_Database_Geo-referenced_population-specific_microsatellite_data_across_the_American_continents/7207514/3?file=23392853
## Data -----
mpg <- read.csv('mpg_sheet2.csv', h=T)
# MacroPopGen database with updated taxonomy: https://figshare.com/articles/dataset/Data_from_Population_demography_maintains_biogeographic_boundaries/19879864?file=35296222
mpg_bg <- read.csv('GD_dist_figshare.csv', header = TRUE)

mpg_no_tax <- mpg %>% select(PopId, RefID, n, FST, GlobalFST, Min.Distance, Max.Distance, Mean.Distance, N.Population, Lat, Long)

mpg2 <- merge(mpg_bg %>% select(species, family, order, class, PopId, RefID), mpg_no_tax, all.x=TRUE, all.y = FALSE)


mpg3 <- mpg2 %>% 
        filter(n>=5, N.Population>=10) %>% # Remove sites with fewer than 5 individuals and fewer than 10 sites
        drop_na(Lat, Long) %>%             # Remove sites with no location data
        select(species, family, order, GlobalFST, FST, PopId, RefID, Lat, Long) %>% # Select columns
        mutate(genus = unlist(lapply(strsplit(species, '_'), function(x) x[1]))) %>%
        group_by(species, RefID) %>% 
        mutate(n_pops = n()) %>% 
        ungroup() %>% 
        filter(n_pops >=10) %>% 
        select(order, family, order, genus, species, RefID, PopId, GlobalFST, FST, Lat, Long) %>% 
        rename(lon = Long, lat = Lat)


# Global FST data:
mpg_glob <- mpg3 %>% 
  filter(GlobalFST=='Global') %>% 
  drop_na(FST) %>% 
  rename(gFST = FST)

# Population FST data:
mpg_pop <- mpg3 %>% 
  filter(GlobalFST=='Population') %>% 
  drop_na(FST) %>% 
  group_by(RefID, species) %>% 
  mutate(gFST = mean(FST))

# 
mpg.f <- rbind(mpg_glob %>% select(-GlobalFST), mpg_pop %>% select(order, family, genus, species, RefID, PopId, gFST, lat, lon))
mpg.f <- mpg.f %>% 
  rename(global_fst = gFST)

mpg.f$species_ref <- paste(mpg.f$species, mpg.f$RefID, sep='_')

## Traits -----
### Area of sample extent --------
sp_ref <- unique(mpg.f$species_ref)

bbox_area <- list()

for(i in 1:length(sp_ref)){
  species_coord <- mpg.f %>% 
    dplyr::filter(species_ref %in% sp_ref[i])
  
  sites <- st_as_sf(species_coord, coords = c("lon", "lat"), crs = 4326)
  
  sites_bbox <- st_bbox(sites)
  bboxpoly <- makepolygon(sites_bbox, sp_ref[i])
  bboxt <- vect(bboxpoly)
  area <- expanse(bboxt, unit = 'km', transform = TRUE)
  bbox_area[[i]] <- area
}

mpg_area <- data.frame(species_ref = sp_ref,
                       area_km2 = unlist(bbox_area))
### Mass -----
traits <- pantheria %>% 
  select(Genus, Species, AdultBodyMass_g) %>% 
  mutate(Genus_species = paste(Genus, Species, sep = '_'))

# Update Pekania
traits$Genus_species[which(traits$Genus_species=='Martes_pennanti')] <- 'Pekania_pennanti'
mpg.f$species[which(mpg.f$species=='Martes_pennanti')] <- 'Pekania_pennanti'
mpg_sp <- unique(mpg.f$species)

# Find missing species or mismatched names 
setdiff(mpg_sp, unique(traits$Genus_species))

traits$Genus_species[which(traits$Genus_species=='Spermophilus_brunneus')] <- 'Urocitellus_brunneus'

# Add mass information
traits <- rbind(traits, c('Lama', 'guanicoe', 115000, 'Lama_guanicoe'))

# Filter and merge
mpg_massg <- traits %>% 
  #drop_na(AdultBodyMass_g) %>% 
  filter(Genus_species %in% mpg_sp) %>% 
  select(-c(Genus, Species)) %>% 
  rename(species = Genus_species)

# Check NAs
mpg_massg[which(is.na(mpg_massg$AdultBodyMass_g)),]

# Missing data for Tamias striatus
mpg_massg$AdultBodyMass_g[which(is.na(mpg_massg$AdultBodyMass_g))] <- 90.5



### Range area -------

# Match species across datasets:
mpg_sp <- unique(mpg.f$species) # Update P. pennanti in species list
species <- match(mpg_sp, mam$species)

# Find species that didn't match range data:
NAsp <- which(is.na(species))
mpg_sp[NAsp]

## Edit other names:
mam$species <- gsub('Xerospermophilus_mohavensis', 'Spermophilus_mohavensis', mam$species)
mam$species <- gsub('Neotamias_ruficaudus', 'Tamias_ruficaudus', mam$species)

## Check again for NAs:
#which(is.na(match(mpg_sp, mam$species)))

species <- match(mpg_sp, mam$species)

ranges <- mam %>% 
  filter(species %in% c(mpg_sp)) %>% 
  filter(presence == 1,
         origin == 1, 
         seasonal == 1)

# area
range_vect <- vect(ranges)

rareas <- expanse(range_vect, unit = 'km', transform = TRUE)

ranges$area_km <- rareas

mpg_sp_areas <- ranges %>% 
  st_drop_geometry() %>% 
  group_by(species) %>%
  mutate(range_area_km2 = sum(area_km)) %>% 
  select(-area_km) %>%
  select(species, range_area_km2) %>% 
  distinct()
mpg_sp_areas

### Home range -----

# Find missing species
mpghomerange1 <- homerange %>% 
  filter(species %in% mpg_sp) %>% 
  select(genus, species, mean_HR)
(missing <- setdiff(mpg_sp, mpghomerange1$species))

homerange[which(homerange$genus=='Cynomys'), ]
homerange[which(homerange$genus=='Lama'), ] # None
homerange[which(homerange$genus=='Ovis'), ]
homerange[which(homerange$genus=='Vicugna'), ] # None
homerange[which(homerange$genus=='Zapus'), ] # None
homerange[which(homerange$genus=='Spermophilus'), ] # None
homerange[which(homerange$genus=='Tamias'), ] # None
homerange[which(homerange$genus=='Urocitellus'), ] # None

missing_genus <- unlist(lapply(strsplit(setdiff(mpg_sp, mpghomerange1$species), '_'), function(x)x[1]))

missing_g_s <- data.frame(species = missing, 
                          genus = missing_genus)

missing_gdf <- homerange %>% 
  filter(genus %in% missing_genus) %>% 
  group_by(genus) %>% 
  summarise(mean_mean_HR = mean(mean_HR))


missing_repl <- merge(missing_g_s, missing_gdf, by = 'genus') %>% 
  rename(mean_HR = mean_mean_HR)

missing2 <- data.frame(genus = c('Vicugna'), species = c('Vicugna_vicugna'), mean_HR = 0.034)

mpghomerange2 <- rbind(mpghomerange1, missing_repl)

mpghomerange3 <- rbind(mpghomerange2, missing2)

setdiff(mpg_sp, mpghomerange3$species)

mpgHR_km2 <- mpghomerange3 %>% 
  select(species, mean_HR)

### Merge -----
mpg_testset <- Reduce(function(x,y) merge(x,y, by="species", all=TRUE, incomparables = "NA"), 
                   list(mpg.f, mpg_massg, mpg_sp_areas, mpgHR_km2))

nrow(mpg_testset) == nrow(mpg.f)

mpg_testset1 <- merge(mpg_testset, mpg_area, by = 'species_ref', all=TRUE)

mpg_testset2 <- mpg_testset1 %>% 
  rename(mass_g = AdultBodyMass_g,
         mean_HR_km2 = mean_HR) %>% 
  mutate(mass_g = as.numeric(mass_g),
         log_mass = log10(mass_g),
         log_range = log10(range_area_km2),
         log_HR = log10(mean_HR_km2),
         log_area = log10(area_km2)) %>% 
  select(-c(PopId, lat, lon)) %>% 
  filter(species != 'Ctenomys_sp')



# Check for NAs
mpg_testset2[which(is.na(mpg_testset2$mass_g)),]
mpg_testset2[which(is.na(mpg_testset2$mean_HR_km2)),]
mpg_testset2[which(is.na(mpg_testset2$range_area_km2)),]
mpg_testset2[which(is.na(mpg_testset2$area_km2)),]


mpg_test_data <- distinct(mpg_testset2 %>% drop_na(mean_HR_km2, range_area_km2))

### save file -----
#write.csv(mpg_test_data, 'mpg_test_data.csv', row.names = F, quote = F)