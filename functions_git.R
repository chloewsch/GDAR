# Functions

# Genetic data ---------
read_genot <- function(fi){
  # check genetic file exists:
  if(!fi %in% list.files()){
    print(fi) 
    stop('Genotype file is not here (did you set the right working directory?)')
  }
  # read STR or GENEPOP files:
  if(grepl(".str", fi)){
    junk <- read.table(fi)
    s <- read.structure(fi, n.ind = (nrow(junk)), n.loc = ((ncol(junk)-1)/2),
                        onerowperind = TRUE, col.lab = 0, col.pop = 1, row.marknames = 0, ask = FALSE)}
  else {
    k <- try(read.genepop(fi, ncode = 2))
    s <- if(inherits(k, "try-error")){read.genepop(fi, ncode = 3)}
    else {
      s <- (read.genepop(fi, ncode = 2))}
  }
  return(s)
}

# GDAR Sampling --------------

# radii
get_radii <- function(microsat, metadata, percents, eecrs){
species_coord <- metadata %>% 
  dplyr::filter(pop %in% pop(microsat))

sites <- st_as_sf(species_coord, coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform(crs = st_crs(eecrs))

sites_bbox <- st_bbox(sites)
lon_range <- sites_bbox[3] - sites_bbox[1]
lat_range <- sites_bbox[4] - sites_bbox[2]
widest_dim <- max(c(lon_range, lat_range))
print(c(lon_range, lat_range))

radii <- (widest_dim/2)*percents # radius is half the length of the widest side of the box
return(radii)
}


ctr_nested_site_samp <- function(microsat, radius, data, eecrs){
  species_coord <- data %>% 
    dplyr::filter(pop %in% pop(microsat))
  
  sites_cea <- st_as_sf(species_coord, coords = c("lon", "lat"), crs = 4326) %>% 
    st_transform(crs = st_crs(eecrs))
  
    site_bufs <- list()
  for(i in 1:nrow(sites_cea)){
    #print(paste0('Site sample #: ', i))
    centrs <- sites_cea[i,]
    site_nest <- list()
    for(j in 1:length(radius)){
      #print(paste0('Radius: ', j))
      q <- st_as_sfc(st_bbox(st_buffer(centrs, dist = radius[j]))) # radius
      site_nest[[j]] <- q
      }
    site_bufs[[i]] <- site_nest
  }
  
  ## Identify sites within each boundary:
    sites_in_bounds <- lapply(site_bufs, function(p){
    lapply(p, function(h) st_intersection(sites_cea, h))
    }
  )
  
    pop_names_inb <- lapply(sites_in_bounds, function(p){
    lapply(p, function(h) h$pop)
  }
  )
  
  pop_names_inb_no0 <- lapply(pop_names_inb, function(x) x[lapply(x,length)>0]) # drop squares with no sites inside
  
  pools <- lapply(pop_names_inb_no0, function(y) lapply(y, function(g){
    if(length(g) == 1){
      separate_pops <- seppop(microsat, drop = TRUE)
      separate_pops[[g]]
      #print('only 1 pop!')
    }
    else({
      separate_pops <- seppop(microsat, drop = FALSE)
      repool(separate_pops[g])})
  })
  )
  
  
    ## Regroup all individuals in each grid as 1 population
  new_pops <- lapply(pop_names_inb_no0, function(u)
    paste0("pool_", c(1:length(u))))
  
  pools_newpop <- pools
  
  for(i in 1:length(pools_newpop)){
    for(j in 1:length(pools_newpop[[i]])){
      pop(pools_newpop[[i]][[j]]) <- rep(new_pops[[i]][[j]], length(pop(pools_newpop[[i]][[j]])))
    }
  }
  
  ### summarize
  allele_count <- lapply(pools_newpop, function(a){
    unlist(lapply(a, function(b) sum(b@loc.n.all)))
    }
  )
  allelic_richness <- lapply(pools_newpop, function(x){
    unlist(lapply(x, function(y) {
      mean(allelic.richness(y, min.n = min(species_coord$num_individuals)*2, diploid = TRUE)$Ar[,1], na.rm = TRUE)}
    ))
  }
      )
  
 
  He <- lapply(pools_newpop, function(a){
    unlist(lapply(a, Hs))
    }
  )
  He <- lapply(He, function(x) x[!is.nan(x)]) # remove dummy populations with NaN values for He
  individuals <- lapply(pools_newpop, function(a){ 
    unlist(lapply(a, function(b) summary(b)$n.by.pop))
    }
  )
  
  
  num_sites <- lapply(pop_names_inb_no0, function(a){
    unlist(lapply(a, length))
    })
  
  
  areaas <- unlist(lapply(site_bufs[[1]], st_area)) # gives the sampled area size (the same across all random sites, so just take from the first)
  
  # Attach areas
  area_repeat_len <- lapply(pop_names_inb, function(x){
    replen <- unlist(lapply(x, length))
    rl <- ifelse(replen>0, 1, 0)
    }
  )
  
  area_m2 <- lapply(area_repeat_len, function(x){
    rep(areaas, x)
  })
  
  
  # Combine all the data
  dat0 <- Map(function(allele_count, allelic_richness, He,
                       individuals, num_sites, area_m2){
    cbind.data.frame(allele_count, allelic_richness, He,
                     individuals, num_sites, area_m2)
  }, allele_count = allele_count, allelic_richness = allelic_richness, He = He, 
  individuals = individuals, num_sites = num_sites, area_m2 = area_m2)
  
  dat <- do.call('rbind', dat0)
  return(list(dat, site_bufs))
}





# estimate z
get.z2 <- function(area_data, variable){
  if(nrow(area_data)>1){
  mod <- try(lm(variable~ log_area, data = area_data))
  if(summary(mod)$coefficients['log_area', 'Pr(>|t|)']<0.05){
    intc <- coef(mod)[1]
    zgd <- coef(mod)[2]
    zse <- summary(mod)$coefficients[2,2]
    r2 <- summary(mod)$r.squared} else{
      intc <- NA
      zgd <- NA
      zse <- NA
      r2 <- NA
    }} else{
      intc <- NA
      zgd <- NA
      zse <- NA
      r2 <- NA
    }
df <- data.frame(zGDAR = zgd,
                 z_se = zse,
                 int_c = intc,
                 rsq = r2)
return(df)
}


# Plots ---------------
plot_common <- function(classt=NULL, legend=TRUE){
  if(is.null(classt)) {
    dat <- z.df
    font = 'Roboto Medium'
    fontsize = 14
  }
  
  else{
    dat <- z.df %>% filter(class==classt)
    font = "Roboto"
    fontsize = 11
  }
  
  p <- ggplot(dat, aes(y = z, x = global_fst, color = genetic_metric, fill = genetic_metric)) +
    geom_pointrange(aes(ymin = z-z_se, ymax = z+z_se), shape=21, fatten = 0.5, size =0.25) +
    geom_smooth(method='lm', alpha = 0.2) +
    scale_color_manual(values = ora) +
    scale_fill_manual(values = ora) +
    labs(y = 'z', x = expression(paste("global ", F[ST])), fill = '', color = '') +
    theme_minimal(base_size = fontsize) +
    theme(text=element_text(family=font))
  
  if(legend == TRUE){
    return(p)
  } else return(p+ theme(legend.position = 'none'))
  
}

oranges <- c('#6c0730', '#a61354', '#dc2b00', '#f08700')

ora<- c('#012030', '#13678A', '#45C4B0')

plot_gdar <- function(random_site_output, variable, zvar){
  ylabel <- gsub("_", ' ', variable)
  speciesname <- strsplit(random_site_output$filename[1], "\\_|\\.str|\\.gen")[[1]][c(2,3)]
  speciesname <- paste(speciesname[1], speciesname[2], sep=' ')
  title <- paste0(speciesname, 
                 '; z=', 
                 round(random_site_output[,zvar], 3), 
                 '; Fst=',
                 round(random_site_output$global_fst, 3))

  p <- ggplot(data = random_site_output, aes(x = area_m2, y = random_site_output[,variable], 
                               color = individuals)) +
  geom_jitter(width = ggplot2:::resolution(random_site_output$area_m2, FALSE) * 0.2,
              height = ggplot2:::resolution(random_site_output[,variable], FALSE) * 0.2) +
  geom_smooth(method='lm', formula = y~log10(x), color = '#000000', fill = '#ff0000', alpha = 0.25) +
  scale_colour_gradientn(
    limits  = range(random_site_output$individuals),
    colours = oranges[c(1, seq_along(oranges), length(oranges))],
    guide="none") +
  labs(title = title,
       y = ylabel, x = 'area (m2)',
       color = "# individuals") +
  theme_minimal(base_size = 11) #+
    if(variable == 'gene_diversity'){
    plot(p+ geom_hline(yintercept = random_site_output$global_ht, lty = 'dashed'))
  } else(plot(p))
}

# log-log GDAR
loglog_plot <- function(random_site_output, variable, zvar){
  ylabel <- gsub("_", ' ', variable)
  speciesname <- strsplit(random_site_output$filename[1], "\\_|\\.str|\\.gen")[[1]][c(2,3)]
  speciesname <- paste(speciesname[1], speciesname[2], sep=' ')
  title <- paste0(speciesname, 
                  '; z=', 
                  round(random_site_output[,zvar], 3))
  
  ggplot(data = random_site_output, aes(x = log_area, y = random_site_output[,variable], 
                                  color = individuals)) +
  geom_jitter(width = ggplot2:::resolution(random_site_output$log_area, FALSE) * 0.2,
                height = ggplot2:::resolution(random_site_output[,variable], FALSE) * 0.2) +
  geom_smooth(method = 'lm', color = 'black', fill = '#ff0000', alpha = 0.25) +
  scale_colour_gradientn(
      limits  = range(random_site_output$individuals),
      colours = oranges[c(1, seq_along(oranges), length(oranges))],
      guide="none") +
  labs(title = title, 
       y = ylabel, x = 'log area',
       color = "# individuals") +
    theme_minimal(base_size = 11)
}

perc_loglog_plot <- function(random_site_output, variable, zvar){
  ylabel <- gsub("_", ' ', variable)
  speciesname <- strsplit(random_site_output$filename[1], "\\_|\\.str|\\.gen")[[1]][c(2,3)]
  speciesname <- paste(speciesname[1], speciesname[2], sep=' ')
  title <- paste0(speciesname, 
                  '; z=', 
                  round(random_site_output[,zvar], 3))
  
  ggplot(data = random_site_output, aes(x = perc_area, y = random_site_output[,variable], 
                                  color = individuals)) +#, 
                                  #size = individuals)) +
  geom_jitter(width = ggplot2:::resolution(random_site_output$log_area, FALSE) * 0.2,
                height = ggplot2:::resolution(random_site_output[,variable], FALSE) * 0.2) +
  geom_smooth(method = 'lm', color = 'black', fill = '#ff0000', alpha = 0.25) +
  scale_colour_gradientn(
      limits  = range(random_site_output$individuals),
      colours = oranges[c(1, seq_along(oranges), length(oranges))],
      guide="none") +
  labs(title = title, 
       y = ylabel, x = 'log area',
       color = "# individuals") +
    theme_minimal(base_size = 11) 
    
}

cvplot <- function(cvdf, colo, xaxislab=TRUE){
  if(xaxislab==TRUE){
    histo<-ggplot() +
      geom_histogram(data = cvdf, aes(x = cv*100), fill = colo, alpha = 0.7, bins = 20) +
      labs(y = 'frequency', x = 'CV of prediction error (%)') +
      scale_x_continuous(limits=c(0, 1.4*100)) +
      theme_minimal(base_size = 11) +
      theme(text=element_text(family='Roboto Medium'))
  } else
  {
    histo<- ggplot() +
      geom_histogram(data = cvdf, aes(x = cv*100), fill = colo, alpha = 0.7, bins = 20) +
      scale_x_continuous(limits=c(0, 1.4*100)) +
      labs(y = 'frequency', x = '') +
      theme_minimal(base_size = 11) +
      theme(text=element_text(family='Roboto Medium'))}
  return(histo)
}

# Bounding box to polygon -----
makepolygon <- function(bbox, spp){ # turn a bounding box into a polygon
  bbox_2_df <- data.frame(binomial = spp,
                          lon = c(bbox[1], bbox[1], bbox[3], bbox[3]), # box XY coordinates
                          lat = c(bbox[2], bbox[4], bbox[4], bbox[2]))
  bbox_2_df$lon <- ifelse(bbox_2_df$lon >=180, 179.999, bbox_2_df$lon)
  bbox_polygon <- bbox_2_df %>% # make coordinates spatial; convert points -> polygon
    st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>%
    group_by(binomial) %>% 
    summarise(geometry = st_combine(geometry)) %>%
    st_cast('POLYGON')
  crs <- st_crs(bbox_polygon)
  # Next steps are to avoid warping the box when using unprojected coordinates
  bbox_polygon <- st_set_crs(bbox_polygon, NA)
  bbox_seg <- bbox_polygon %>% sf::st_segmentize(1) %>% st_set_crs(crs) # segmentize box edges
}

cvplot <- function(cvdf, colo, xaxislab=TRUE){
    if(xaxislab==TRUE){
    histo<-ggplot() +
    geom_histogram(data = cvdf, aes(x = cv*100), fill = colo, alpha = 0.7, bins = 20) +
    labs(y = 'frequency', x = 'CV of prediction error (%)') +
    scale_x_continuous(limits=c(0, 1.4*100)) +
    theme_minimal(base_size = 11) +
    theme(text=element_text(family='Roboto Medium'))
  
  } else{
    histo<- ggplot() +
    geom_histogram(data = cvdf, aes(x = cv*100), fill = colo, alpha = 0.7, bins = 20) +
    scale_x_continuous(limits=c(0, 1.4*100)) +
    labs(y = 'frequency', x = '') +
    theme_minimal(base_size = 11) +
    theme(text=element_text(family='Roboto Medium'))}
  
  return(histo)
  
}
