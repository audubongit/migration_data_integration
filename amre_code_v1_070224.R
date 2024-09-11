# setup ------------------------------------------------------------------------
# time zone
Sys.setenv(TZ="America/New_York") # should be ET

# load libs
library(sf)
library(raster)
library(ebirdst)
library(ggplot2)
library(concaveman)
library(purrr)
library(gdistance)
library(doParallel)
library(stringr)
library(tidyr)
library(dplyr)
library(terra)
library(stats)
library(foreach)

# load ebird key and then restart
# set_ebirdst_access_key("XXXXXX", overwrite=T)

# settings
extract <- raster::extract 
select <- dplyr::select
options(scipen=9999999)

# directories
setwd("D:/Users/tmeehan/Box/-.Tim.Meehan individual/amre_mmm/") # set to code directory

# crs
mbi_analysis <- "+proj=aeqd +lat_0=15 +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
ebird_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# get cores helper function
get_cores <- function(ras=as_breed,
                      pols=breeding_mcrs,
                      keeps=0.99) { 
  cr <- crs(ras) 
  rs <- res(ras)
  ex <- ext(ras)
  outs <- c()
  nams <- unique(pols$mcr_class)
  for(i in 1:length(nams)){
    test1 <- mask(ras, vect(pols[i,]))
    test2 <- raster(test1)
    test3 <- rasterToPoints(test2)
    test4 <- as.data.frame(test3)
    bp_i <- test4 %>%
      rename(z=3) %>%
      arrange(desc(z)) %>%
      mutate(z2=z/sum(z)) %>%
      mutate(cs=cumsum(z2)) %>%
      mutate(keep=cs / max(cs)) %>%
      filter(keep <= keeps) %>%
      select(1:3) %>%
      mutate(mcr_class=nams[i],
             core=1)
    outs <- rbind(outs, bp_i)
    print(paste("finished", i, "of", length(nams)))
  }
 return(outs)
}

# get path helper function
get_path <- function(wp = all_win_points, sp = all_sum_points, 
                   cond = spring_conduct, j = i, seas="prebreeding",
                   theta_i = 0.003){
  require("raster")
  require("gdistance")
  require("sf")
  o <- as(wp[j,], "Spatial")
  g <- as(sp[j,], "Spatial")
  t <- passage(x=cond, 
               origin=o,
               goal=g, 
               theta=theta_i, # dial
               output="RasterLayer")
  names(t) <- paste(seas, paste("start", o$mcr_class, sep="."),
        paste("stop", g$mcr_class, sep="."), sep=".")
  return(t)
}

# from grid to df for ggplot helper function
gplot_data <- function(x, maxpixels = 50000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  dat <- tibble::as_tibble(data.frame(coords, dat))
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  dat
}

# rescale helper functions
rescale_it <- function(x){
  m <- cellStats(x, max, na.rm=T)
  y <- x / m
  return(y)
}
rescale_01 <- function(x){
  minx <- min(x, na.rm=T)
  maxx <- max(x, na.rm=T)
  y <- (x-minx) / (maxx - minx)
  return(y)
}
rescale_ras_01 <- function(x, new.min = 0, new.max = 1) {
  x.min = quantile(x, probs=c(0.01), na.rm=T)
  x.max = quantile(x, probs=c(0.99), na.rm=T)
  x[x<x.min] <- x.min
  x[x>x.max] <- x.max
  y <- new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
  return(y)
}
raster01 = function(r){
  minmax_r = range(values(r), na.rm=TRUE) 
  return( (r-minmax_r[1]) / (diff(minmax_r)))
}

# baseline helper function
add_baselines <- function(x, lmask=land_mask, 
                             base_land=0.01,
                             base_ocean=0.01){
    y <- x
    y[(is.na(x) & lmask==1) | x <= base_land] <- base_land
    y[is.na(x)] <- base_ocean
    return(y)
}

# view simple sf dataframe
vu <- function(x) {x %>% st_drop_geometry() %>% View()}

# ggplot themes helper function
theme_map <- function(base_size = 9, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.title = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0, "lines"),
          plot.background = element_blank(),
          legend.background=element_rect(fill=NA, colour=NA),
          legend.direction="vertical",
          legend.key=element_rect(fill=NA, colour="white"),
          legend.text.align=1,
          legend.text = element_text(size=9),
          legend.title=element_text(hjust=0, size=11),
          legend.justification=c(0, 0.5),
          plot.title = element_text(size=14, hjust = 0.7))
}
theme_timeseries <- function (base_size = 11, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.9), angle = 0),
          axis.text.y = element_text(size = rel(0.9), angle = 0),
          strip.background = element_rect(fill = "grey50"),
          legend.key = element_rect(fill = "white", colour = NA),
          plot.title = element_text(size=14, hjust = 0.5,
                                    margin=margin(t=5, b=10)),
          legend.position="right",
          complete = TRUE)
}; theme_set(theme_timeseries())
  
# get hemisphere map
west_hem <- read_sf("./input/west_hemi_ebird.shp") %>% 
  st_transform(mbi_analysis) %>%
  st_simplify(dTolerance = 5000) %>% 
  dplyr::select(1) %>% 
  mutate(mask=1) %>% 
  select(mask)
# plot(west_hem$geometry)
# ------------------------------------------------------------------------------



# get species table and number of points ---------------------------------------
spp_table <- read.csv("./input/spp_tab.csv")
spec_list <- spp_table %>% pull(species) %>% as.character()
overwater_pre_lst <- spp_table %>% pull(overwater_pre)
overland_pre_lst <- spp_table %>% pull(overland_pre)
overwater_post_lst <- spp_table %>% pull(overwater_post)
overland_post_lst <- spp_table %>% pull(overland_post)
theta_pre_lst <- spp_table %>% pull(theta_pre)
theta_post_lst <- spp_table %>% pull(theta_post)
# spec_list %in% ebirdst_runs$common_name

# set the number of points/paths per mcr
n_connect_points <- 5 # use around 5 for fdevelopment and 50 for good results
# ------------------------------------------------------------------------------




# start for --------------------------------------------------------------------
# for(s in 1:length(spec_list)){
s <- 1
# ------------------------------------------------------------------------------




# pick species------------------------------------------------------------------
spp1 <- spec_list[s]
spp2 <- ebirdst_runs$species_code[ebirdst_runs$common_name==spp1]
overwater_pre <- overwater_pre_lst[s]
overland_pre <- overland_pre_lst[s]
overwater_post <- overwater_post_lst[s]
overland_post <- overland_post_lst[s]
theta_pre <- theta_pre_lst[s]
theta_post <- theta_post_lst[s]
print(paste("Starting", spp1))

# save stuff
out_list <- list()
out_list$species_code <- spp2
out_list$species_name <- spp1
out_list$hemisphere_map_sf <- west_hem
# ------------------------------------------------------------------------------



  
# get migration occurrence surfaces --------------------------------------------
dl_path <- ebirdst_data_dir()
ebirdst_download_status(species = spp2, 
                        path=dl_path,
                        download_abundance = T,
                        download_occurrence = T,
                        pattern="seasonal")
occs <- load_raster(species = spp2, 
                    product = "occurrence", 
                    period = "seasonal", 
                    metric = "max", 
                    resolution = "27km") # 27km2

# view spp stats
ebirdst_runs %>%
  filter(species_code == spp2) %>%
  glimpse()
  
# combine, trim, and reproject migration occ stack
occs_stack <- rast(list(occs$prebreeding_migration, 
                        occs$postbreeding_migration)) %>% 
  project(res=res(occs), y=mbi_analysis) %>%
  crop(west_hem, mask=T) %>% 
  trim(.) 
names(occs_stack) <- c("prebreeding_migration_occurrence", 
                        "postbreeding_migration_occurrence")

# rescale occs 0 to 1 after # added
occs_stack[[1]] <- raster01(occs_stack[[1]])
occs_stack[[2]] <- raster01(occs_stack[[2]])
# plot(occs_stack)

# save it
out_list$occurrence_stack <- occs_stack
# ------------------------------------------------------------------------------
    



# get stationary abundance surfaces --------------------------------------------
# stack abund stems
abd <- load_raster(species = spp2, 
                    product = "abundance", 
                    period = "seasonal", 
                    metric = "mean", 
                    resolution = "27km") # 27km2

# combine, trim, and reproject migration occ stack
abund_stack <- rast(list(abd$breeding, abd$nonbreeding)) %>% 
  project(res=res(abd), y=mbi_analysis) %>% 
  crop(west_hem, mask=T) %>% 
  trim(.) 
names(abund_stack) <- c("breeding_abundance", "nonbreeding_abundance")
# plot(abund_stack)
# ------------------------------------------------------------------------------
  



# create land mask for baselining ----------------------------------------------
land_mask <- rasterize(west_hem, occs_stack)
  
# make migration stack and add baseline conductance values set above
mig_stack <- rast(list(occs_stack$prebreeding_migration_occurrence,
                    occs_stack$postbreeding_migration_occurrence))
mig_stack$prebreeding_migration_occurrence <- 
  add_baselines(x=mig_stack$prebreeding_migration_occurrence,
                lmask=land_mask, 
                base_land=overland_pre,
                base_ocean=overwater_pre)  
mig_stack$postbreeding_migration_occurrence <- 
  add_baselines(x=mig_stack$postbreeding_migration_occurrence, 
                lmask=land_mask, 
                base_land=overland_post,
                base_ocean=overwater_post)
# plot(mig_stack)
# ------------------------------------------------------------------------------
  
  


# get mcrs and select breeding abundance cores per mcr -------------------------
mcrs <- read_sf(paste0("./input/", spp2, "_bird_migration_regions.shp")) %>%
  st_transform(mbi_analysis) %>%
  st_simplify(dTolerance = 5000) 
breeding_mcrs <- mcrs %>% filter(season=="breeding") 
nonbreeding_mcrs <- mcrs %>% filter(season=="nonbreeding") 
# plot(mcrs$geometry)
print("got all maps")
# ------------------------------------------------------------------------------
  



# select summer and winter core areas per mcr ----------------------------------
# get rid of zero abundance for speed
as_breed <- abund_stack$breeding_abundance
as_nonbreed <- abund_stack$nonbreeding_abundance
as_breed[as_breed<=0] <- NA
as_nonbreed[as_nonbreed<=0] <- NA
  
# run get_cores function
breed_cores <- get_cores(ras=as_breed, pols=breeding_mcrs, keeps=0.99)
nonbreed_cores <- get_cores(ras=as_nonbreed, pols=nonbreeding_mcrs, 
                            keeps=0.99)
  
# make core rasters
breed_core_ras <- rast(breed_cores[, c(1,2,3)], type="xyz", crs=crs(as_breed))
nonbreed_core_ras <- rast(nonbreed_cores[, c(1,2,3)], type="xyz",
                                    crs=crs(as_breed))
# plot(breed_core_ras) 
# plot(west_hem$geometry, add=TRUE) 
# plot(nonbreed_core_ras)
# plot(west_hem$geometry, add=TRUE) 
# ------------------------------------------------------------------------------
  



# select breeding points within cores and match with winter points -------------
breeding_points <- breed_cores %>% group_by(mcr_class) %>% 
  sample_n(n_connect_points, replace=T) %>% dplyr::select(1,2,4) %>% ungroup()
# plot(st_geometry(breeding_mcrs), col="gray")
# plot(st_as_sf(breeding_points, coords=c("x","y")), add=TRUE)
  
# list breeding mcrs for species
breeding_mcr_list <- sort(mcrs$mcr_class[mcrs$season=="breeding"])
nonbreeding_mcr_list <- sort(mcrs$mcr_class[mcrs$season=="nonbreeding"])
  
# pull in the connectivity matrix if there is one. if not, then script assumes
# equal connectivity across mcrs
cm1 <- paste0("./input/", spp2, "_MigConnect_RANGEfiltered_Results_Jan2023.csv")
if(file.exists(cm1)) { 
  mcmat <- read.csv(cm1) %>% select(proportion=mean, mcr_class_breeding=brdg,
                                    mcr_class_nonbreeding=nonbrdg)
} else mcmat <- expand.grid(proportion=1/length(unique(nonbreeding_mcr_list)),
            mcr_class_breeding=breeding_mcr_list, 
            mcr_class_nonbreeding=nonbreeding_mcr_list)
# head(mcmat)

# match winter and summer points per breeding mcr based on connectivity matrix
all_sum_points <- c()
all_win_points <- c()
for(m in 1:length(breeding_mcr_list)){
  print(paste("m=",m))
  # subset breeding points
  breeding_points_i <- breeding_points %>% 
    filter(mcr_class==breeding_mcr_list[m]) 
  n_pts_i <- nrow(breeding_points_i)
  # subset mcmat and get sample sizes per winter poly
  mcmat_i <- mcmat %>% filter(mcr_class_breeding==breeding_mcr_list[m]) %>%
    select(mcr_class=mcr_class_nonbreeding, proportion) %>%
    mutate(samp_n=floor(n_pts_i*proportion))
  samp1 <- sample(1:nrow(mcmat_i), n_pts_i-sum(mcmat_i$samp_n), 
                prob=mcmat_i$proportion, replace=T)
  row1 <- as.integer(names(table(samp1)))
  add1 <- as.integer(table(samp1))
  mcmat_i[row1,3] <- mcmat_i[row1,3] + add1
  winter_mcrs_i <- as.data.frame(mcmat_i[,c(1,3)])
  # randomly locate nonbreeding sites
  sum_points <- c()
  win_points <- c()
  for(t in 1:nrow(winter_mcrs_i)){
    print(paste("t=",t))
    target_feature <- winter_mcrs_i[t, "mcr_class"]
    target_size <- winter_mcrs_i[t, "samp_n"]
    sum_samps <- breeding_points_i %>% 
      sample_n(size=target_size, replace=T) %>%
      arrange(x, y)
    win_samps <- nonbreed_cores %>% filter(mcr_class==target_feature) %>%
      sample_n(size=target_size, replace=T) %>% select(1,2,4) %>%
      arrange(x, y)
    sum_points <- rbind(sum_points, sum_samps)
    win_points <- rbind(win_points, win_samps)
  }
    all_sum_points <- rbind(all_sum_points, sum_points)
    all_win_points <- rbind(all_win_points, win_points)
}
  
# make points into sf objects
all_sum_points <- all_sum_points %>% 
  st_as_sf(coords=c("x", "y"), crs=mbi_analysis)
all_win_points <- all_win_points %>% 
  st_as_sf(coords=c("x", "y"), crs=mbi_analysis)
# plot(west_hem$geometry)
# plot(all_sum_points$geometry, add=T, col="red")
# plot(all_win_points$geometry, add=T, col="blue") 
#-------------------------------------------------------------------------------
    



# make conductance surfaces from migration occs surfaces -----------------------
envelope <- read_sf(paste0("./input/", spp2, "_envelope.shp")) %>%
  st_transform(mbi_analysis) %>%
  st_simplify(dTolerance = 5000) %>%
  st_buffer(dist=100000)
out_list$envelope <- envelope
# plot(west_hem$geometry)
# plot(envelope$geometry, add=TRUE)

# spring
spr_con <- mask(crop(mig_stack$prebreeding_migration_occurrence, 
                   vect(envelope)), vect(envelope))
spr_con <- raster(spr_con)
spring_conduct <- transition(spr_con, transitionFunction=mean, directions=8) %>%
  geoCorrection(scl=T)
# plot(raster(spring_conduct)) 
# plot(all_win_points$geometry, add=T, pch=1, col="blue")
# plot(all_sum_points$geometry, add=T, pch=1, col="red")

# fall
fal_con <- mask(crop(mig_stack$postbreeding_migration_occurrence, 
                   vect(envelope)), vect(envelope))
fal_con <- raster(fal_con)
fall_conduct <- transition(fal_con, transitionFunction=mean, directions=8) %>%
  geoCorrection(scl=T)  
# plot(raster(fall_conduct))
# plot(all_win_points$geometry, add=T, pch=1, col="blue")
# plot(all_sum_points$geometry, add=T, pch=1, col="red")
# ------------------------------------------------------------------------------
  



################################################################################
################################################################################
# get probability of spring passage --------------------------------------------
print("getting spring corridors")
  
# make cluster
cl <- makeCluster(5)
registerDoParallel(cl)
prob_stack <- foreach(i=1:nrow(all_sum_points), 
                    .packages=c("raster", "gdistance", "sf")) %dopar% {
  get_path(wp=all_win_points, sp=all_sum_points, 
            cond=spring_conduct, j=i, seas="prebreeding",
           theta_i = theta_pre)
}
stopCluster(cl)
  
# make stack from parallel output list and add names
prob_stack <- stack(prob_stack)
# plot(prob_stack[[c(1:3, 7:9)]])

# summary layer
# prebreeding_migration_mean <- rescale_it((mean(prob_stack))) # choose one
prebreeding_migration_mean <- rescale_it(sqrt(mean(prob_stack))) # choose one
names(prebreeding_migration_mean) <- "prebreeding.migration.mean"
# plot(((prebreeding_migration_mean)))
# plot(st_geometry(west_hem), add=T)
# ------------------------------------------------------------------------------
  
  
# get probability of fall passage ----------------------------------------------
print("getting fall corridors")
cl <- makeCluster(5)
registerDoParallel(cl)
prob_stack2 <- foreach(i=1:nrow(all_sum_points), 
                    .packages=c("raster", "gdistance",
                                                    "sf")) %dopar% {
  get_path(wp=all_sum_points, sp=all_win_points, 
            cond=fall_conduct, j=i, seas="postbreeding",
           theta_i = theta_post)
}
stopCluster(cl)
  
# make stack from parallel list output list and add names
prob_stack2 <- stack(prob_stack2)
# plot(prob_stack2[[c(1:3, 7:9)]])
  
# summary layer
# postbreeding_migration_mean <- rescale_it(((mean(prob_stack2)))) # choose one
postbreeding_migration_mean <- rescale_it(sqrt((mean(prob_stack2)))) # choose one
names(postbreeding_migration_mean) <- "postbreeding.migration.mean"
# plot(((postbreeding_migration_mean)))
# plot(st_geometry(west_hem), add=T)
# ------------------------------------------------------------------------------


# stack and save all lcps ------------------------------------------------------
lcp_stack <- stack(prebreeding_migration_mean,
                    postbreeding_migration_mean)

# save lcp stack
dir.create(paste0("./output/"))
terra::writeRaster(lcp_stack, paste0("./output/",
                        gsub(" ", "_", spp1), "_",
                                "lcp_stack_27km.tif"), overwrite=T)

# stash it
out_list$lcp_stack <- lcp_stack
# ------------------------------------------------------------------------------
################################################################################
################################################################################




# reopen lcps ------------------------------------------------------------------
lcp_stack <- stack(paste0("./output/", 
                        gsub(" ", "_", spp1), "_",
                                "lcp_stack_27km.tif"))
names(lcp_stack) <- c("prebreeding.migration.mean", "postbreeding.migration.mean")
# ------------------------------------------------------------------------------




# set parameters for integration and validation --------------------------------
species_code <- spec1 <- s1 <- spp2
species_name <- spec2 <- s2 <- gsub(" ", "_", spp1)

# number of background points to compare with presence-only points
n_background <- 10000

# thin_per_day
thin_per_day <- FALSE

# keep only one tracking point per 27 km grid cell?
thin_per_cell <- TRUE

# keep only this many band encounters
max_bands <- 5000
  
# load more libs for integration and validation
library(cowplot)
library(dismo)
library(lubridate)
library(mgcViz)
library(mgcv)
library(scam)
library(pROC)
  
# more options
options(stringsAsFactors = FALSE)
# ------------------------------------------------------------------------------
  



# trim integration area -------------------------------------------------------- 
hem_poly <- st_convex_hull(west_hem)
envelope <- st_intersection(envelope, hem_poly)
out_list$study_region_map_sf <- envelope

# trim lcp maps to integration/validation area
out_list$mbi_crs <- mbi_analysis
crs(lcp_stack) <- mbi_analysis
lcp_map <- lcp_stack[[1:2]]
lcp_map <- raster::crop(lcp_map, envelope) %>% 
  raster::mask(envelope)
names(lcp_map) <- c("prebreeding_migration", "postbreeding_migration")
  
# trim eBird S&T maps to integration/validation area
stem_map <- stack(occs_stack[[1:2]])
stem_map <- raster::resample(stem_map, lcp_map)
stem_map[is.na(stem_map)] <- 0
stem_map <- raster::crop(stem_map, envelope) %>% 
  raster::mask(envelope)
names(stem_map) <- c("prebreeding_migration", "postbreeding_migration")
print("got trimmed maps")

# stash it
out_list$lcp_map_rs <- lcp_map
out_list$stem_map_rs <- stem_map
# ------------------------------------------------------------------------------
  


  
# get, clean, and sample movement data -----------------------------------------
move_dat <- read.csv(paste0("./input/", spec1, "_clean_data.csv")) %>%
  mutate(season=ifelse(season=="summer", "breeding", season)) %>%
  mutate(season=ifelse(season=="winter", "nonbreeding", season)) %>%
  mutate(season=ifelse(season=="spring", "prebreeding_migration", season)) %>%
  mutate(season=ifelse(season=="fall", "postbreeding_migration", season)) %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  drop_na(season) %>%
  mutate(date_time=str_sub(date_time, 1, 10)) %>% 
  mutate(date_time=lubridate::ymd(date_time)) %>% 
  mutate(date_time=lubridate::date(date_time)) %>%
  dplyr::filter(season!="breeding") %>%
  dplyr::filter(season!="nonbreeding") %>%
  dplyr::filter(tech_type %in% c("bnd", "ptt", "gps", "llg", "rtl")) %>% 
  st_as_sf(coords=c("x_coord","y_coord"), crs = 4236) %>%
  st_transform(crs=mbi_analysis) %>% 
  st_join(envelope, join = st_within)
out_list$initial_move_dat_ss <- nrow(move_dat)
  
# remove geolocator data near equinox
move_dat <- move_dat %>% filter(!((move_dat$week==10 |
        move_dat$week==11 |
        move_dat$week==12 |
        move_dat$week==37 |
        move_dat$week==38 |
        move_dat$week==39) & move_dat$tech_type=="llg")) 
  
# thin out banding data
spring_bands <- move_dat[move_dat$season=="prebreeding_migration" &
           move_dat$tech_type=="bnd", ]
if(nrow(spring_bands) > max_bands) spring_bands <- 
  sample_n(spring_bands, max_bands)
fall_bands <- move_dat[move_dat$season=="postbreeding_migration" &
            move_dat$tech_type=="bnd", ]
if(nrow(fall_bands) > max_bands) fall_bands <- sample_n(fall_bands, max_bands)
md0 <- filter(move_dat, tech_type!="bnd")
move_dat <- rbind(md0, fall_bands, spring_bands)
out_list$thinned_bnd_move_dat_ss <- nrow(move_dat)
  
# sample 1 track point per day
if(thin_per_day==TRUE){
  move_dat <- move_dat %>%
    dplyr::group_by(bird_id, date_time, season, tech_type) %>%
    dplyr::sample_n(size = 1) %>%
    dplyr::ungroup()
}
out_list$after_temporal_move_dat_ss <- nrow(move_dat)
  
# sample 1 track point per cell
if(thin_per_cell==TRUE){
  r <- lcp_map$prebreeding_migration
  values(r) <- 1:ncell(r)
  move_dat$grid_id <- raster::extract(r, move_dat)
  # plot(move_dat["grid_id"])
  move_dat <- move_dat %>% 
    dplyr::ungroup() %>%
    dplyr::group_by(bird_id, season, grid_id, tech_type) %>%
    dplyr::sample_n(size = 1)
}
out_list$after_spatial_move_dat_ss <- nrow(move_dat) 
  
# get sample sizes and data holders
out_list$move_dat_ss_after_time_space_filter <- nrow(move_dat)
out_list$tech_type_sample_size <- move_dat %>% as.data.frame() %>% 
  group_by(season) %>% 
  count(tech_type)
out_list$bird_id_sample_size <- move_dat %>% as.data.frame() %>% 
  distinct(bird_id, season) %>%
  count(season)
out_list$data_owners <- move_dat %>% ungroup() %>%
    st_drop_geometry() %>%
    group_by(dataset_key, data_source, study_code,
      tech_type, species_code) %>%
    group_by(tech_type) %>%
    summarise(individual_bird_ids=length(unique(bird_id)),
              records=n()) %>% 
    arrange(desc(records))
out_list$clean_move_dat <- move_dat
  
# get fall data for analysis and add random effects
dat1 <- move_dat %>% ungroup() %>%
  dplyr::filter(season=="postbreeding_migration") %>% 
  dplyr::select(bird_id, date_time, season, tech_type) %>%
  dplyr::mutate(presence=1)
dat1$bird_id <- ifelse(dat1$tech_type=="bnd", "bnd", dat1$bird_id) 
  
# add zeros
dat2 <- st_sample(envelope, size = n_background) %>%
  st_sf() %>% 
  dplyr::mutate(tech_type="ran",
          season="postbreeding_migration",
          date_time=mean(dat1$date_time),
          bird_id=sample(dat1$bird_id, n_background, replace=T)) %>%
  dplyr::select(bird_id, date_time, season, tech_type) %>% 
  dplyr::mutate(presence=0) 

# combine presence and background
dat3 <- rbind(dat1, dat2) 
dat3$lcp_prob <- extract(lcp_map$postbreeding_migration, 
                          dat3)
dat3$stem_prob <- extract(stem_map$postbreeding_migration, 
                          dat3)
dat3 <- dat3[complete.cases(dat3$lcp_prob, dat3$stem_prob),]
out_list$fall_df <- dat3
  
# spring data for analysis
dat1 <- move_dat %>% ungroup() %>%
  dplyr::filter(season=="prebreeding_migration") %>% 
  dplyr::select(bird_id, date_time, season, tech_type) %>% 
  dplyr::mutate(presence=1)
dat1$bird_id <- ifelse(dat1$tech_type=="bnd", "bnd", dat1$bird_id) 

# add zeros
dat2 <- st_sample(envelope, size = n_background) %>%
  st_sf() %>% 
  dplyr::mutate(tech_type="ran",
          season="prebreeding_migration",
          date_time=mean(dat1$date_time),
          bird_id=sample(dat1$bird_id, n_background, replace=T)) %>%
  dplyr::select(bird_id, date_time, season, tech_type) %>% 
  dplyr::mutate(presence=0) 
  
# combine presence and background
dat3 <- rbind(dat1, dat2) 
dat3$lcp_prob <- extract(lcp_map$prebreeding_migration, 
                          (dat3))
dat3$stem_prob <- extract(stem_map$prebreeding_migration, 
                          (dat3))
dat3 <- dat3[complete.cases(dat3$lcp_prob, dat3$stem_prob),]
out_list$spring_df <- dat3
  
# combine all
all_movement_df <- rbind(out_list$spring_df, out_list$fall_df) 
out_list$all_movement_df <- all_movement_df
  
# save integration and validation data for the species
print("got available movement data")
#-------------------------------------------------------------------------------
  
  


# put integration and validation data back into env and make plots -------------
#list2env(out_list, environment())
# plot it
p1 <- all_movement_df %>%
  dplyr::filter(season=="prebreeding_migration", tech_type!="ran") %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "bnd", "Band")) %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "ptt", "PTT")) %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "gps", "GPS")) %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "llg", "LLG")) %>%
  ggplot() +
  geom_tile(data=filter(gplot_data(lcp_map[["prebreeding_migration"]],
                                    maxpixels=ncell(lcp_map)),
                        !is.na(value), value>=0.01), aes(x=x, y=y, fill=value)) +
  geom_sf(data=west_hem, fill=NA) +
  geom_sf(aes(col=tech_type), size=0.1, shape=16) +
  scale_fill_gradientn("Spring\nLCP\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  scale_color_brewer("Data\ntype", palette="Set1") +
  coord_sf(xlim=st_bbox(envelope)[c(1,3)],
            ylim=st_bbox(envelope)[c(2,4)]) +
  labs(x="Longitude", y="Latitude")
p2 <- all_movement_df %>% dplyr::filter(season=="prebreeding_migration",
                                  tech_type!="ran") %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "bnd", "Band")) %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "ptt", "PTT")) %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "gps", "GPS")) %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "llg", "LLG")) %>%
  ggplot() +
  geom_tile(data=filter(gplot_data((stem_map[["prebreeding_migration"]]),
                                    maxpixels=ncell(lcp_map)),
                        !is.na(value), value>=0.01), aes(x=x, y=y, fill=value)) +
  geom_sf(data=west_hem, fill=NA) +
  geom_sf(aes(col=tech_type), size=0.1, shape=16) +
  scale_fill_gradientn("Spring\nSTEM\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  scale_color_brewer("Data\ntype", palette="Set1") +
  coord_sf(xlim=st_bbox(envelope)[c(1,3)],
            ylim=st_bbox(envelope)[c(2,4)]) +
  labs(x="Longitude", y="Latitude")
p3 <- all_movement_df %>% dplyr::filter(season=="postbreeding_migration",
                                  tech_type!="ran") %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "bnd", "Band")) %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "ptt", "PTT")) %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "gps", "GPS")) %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "llg", "LLG")) %>%
  ggplot() +
  geom_tile(data=filter(gplot_data((lcp_map[["postbreeding_migration"]]),
                                    maxpixels=ncell(lcp_map)),
                        !is.na(value), value>=0.01), aes(x=x, y=y, fill=value)) +
  geom_sf(data=west_hem, fill=NA) +
  geom_sf(aes(col=tech_type), size=0.1, shape=16) +
  scale_fill_gradientn("Fall\nLCP\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  scale_color_brewer("Data\ntype", palette="Set1") +
  coord_sf(xlim=st_bbox(envelope)[c(1,3)],
            ylim=st_bbox(envelope)[c(2,4)]) +
  labs(x="Longitude", y="Latitude")
p4 <- all_movement_df %>% dplyr::filter(season=="postbreeding_migration",
                                  tech_type!="ran") %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "bnd", "Band")) %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "ptt", "PTT")) %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "gps", "GPS")) %>%
  dplyr::mutate(tech_type=str_replace(tech_type, "llg", "LLG")) %>%
  ggplot() +
  geom_tile(data=filter(gplot_data((stem_map[["postbreeding_migration"]]),
                                    maxpixels=ncell(lcp_map)),
                        !is.na(value), value>=0.01), aes(x=x, y=y, fill=value)) +
  geom_sf(data=west_hem, fill=NA) +
  geom_sf(aes(col=tech_type), size=0.1, shape=16) +
  scale_fill_gradientn("Fall\nSTEM\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  scale_color_brewer("Data\ntype", palette="Set1") +
  coord_sf(xlim=st_bbox(envelope)[c(1,3)],
            ylim=st_bbox(envelope)[c(2,4)]) +
  labs(x="Longitude", y="Latitude")
pg1 <- plot_grid(p1, p2, p3, p4, nrow = 2); pg1
ggsave2(filename = paste0("./output/", species_name,
                  "_", "lcp_stem_and_movement_map.pdf"),  width=9.6, height=9.6)

# stash it
out_list$lcp_stem_and_movement_map <- pg1
# ------------------------------------------------------------------------------
  



# analyze spring data ----------------------------------------------------------
# choose spring STEM and LCP
focal_season <- "prebreeding_migration"
lcp_map <- out_list$lcp_map_rs[[focal_season]]; names(lcp_map) <- "lcp_prob"
stem_map <- out_list$stem_map_rs[[focal_season]]; names(stem_map) <- "stem_prob"
  
# extract and prepare data for that season
dat1_train <- all_movement_df %>% 
  dplyr::filter(season==focal_season) %>% 
  st_as_sf(sf_column_name="geometry", crs=mbi_analysis) %>%
  dplyr::mutate(bird_id=factor(as.numeric(factor(bird_id)))) %>%
  dplyr::arrange(tech_type, bird_id, date_time)
  
# add a max predictor value
dat1_train$max_prob <-
  apply(dat1_train[,c("lcp_prob", "stem_prob")] %>%
                                st_drop_geometry() %>% as.matrix(), 1,
                              function(x) max(x, na.rm=T))

# models
f0 <- presence ~ 1 
r0 <- scam(f0, data=dat1_train, family="binomial")
f3 <- presence ~ 1 + s(stem_prob, k=5, bs="mpi") + s(lcp_prob, k=5, bs="mpi") 
r3 <- scam(f3, data=dat1_train, family="binomial")
f4 <- presence ~ 1 + s(max_prob, k=5, bs="mpi") 
r4 <- scam(f4, data=dat1_train, family="binomial")

# make a model fit summaries table
aic_tab <- AIC(r0, r3, r4) %>%
  mutate(Model=c("GAM_Null", "GAM_STEM_LCP", "GAM_Max"),
          D_explained=c(summary(r0)$dev.expl,
          summary(r3)$dev.expl,
          summary(r4)$dev.expl)) %>%
  arrange(AIC) %>% mutate(Delta_AIC=c(0, diff(AIC))) %>%
  mutate_if(is.numeric, round, 2) %>%
  select(Model, D_explained, AIC, Delta_AIC)

# partial plot for full GAm model with STEM and LCP
p1 <- function() {
  par(mar = c(3, 3, 1, 1), mgp = c(2, 1, 0), mfrow = c(1,2))
  plot(r3, select=1, shade=TRUE, xlab="STEM Probability", 
       ylab="Standardized effect")
  plot(r3, select=2, shade=TRUE, xlab="LCP Probability", 
       ylab="Standardized effect")
}
p_par <- ggdraw(p1)

# get model fitted values for more fit stats
dat1_train$r0_fit <- predict(r0, dat1_train, type="response", exclude="s(bird_id)")
dat1_train$r3_fit <- predict(r3, dat1_train, type="response", exclude="s(bird_id)")
dat1_train$r4_fit <- predict(r4, dat1_train, type="response", exclude="s(bird_id)")

# get auc and make ROC plot
roc0 <- roc(presence~r0_fit, data=dat1_train) 
roc3 <- roc(presence~r3_fit, data=dat1_train)
roc4 <- roc(presence~r4_fit, data=dat1_train)
p_roc <- ggroc(list(roc3, roc4), legacy.axes = T) + 
    geom_abline(intercept=0, slope=1, lty=2, col="gray60") +
    labs(x="1 - Specificity", y="Sensitivity", color="Model") +
    scale_color_brewer(labels=c("STEM+LCP\nmodel", "Max\nmodel"), palette="Set1") +
  annotate("text", x = 1, y = 0.1, hjust = 1,
            label = paste("STEM+LCP model AUC =",
                          str_pad(round(roc3$auc, 2), width=4,
                                  side ="right", pad=0))) +
  annotate("text", x = 1, y = 0, hjust = 1,
            label = paste("Max model AUC =", 
                          str_pad(round(roc4$auc, 2), width=4, 
                                  side ="right", pad=0)))
  
# add auc to fit table and save
aic_tab$AUC <- NA
aic_tab[aic_tab$Model=="GAM_STEM_LCP", "AUC"] <- round(roc3$auc, 2)
aic_tab[aic_tab$Model=="GAM_Max", "AUC"] <- round(roc4$auc, 2)
row.names(aic_tab) <- NULL
write.csv(aic_tab, paste0("./output/", species_name,
                  "_", focal_season, "_", "aic_table.csv"), row.names=F, na="")
  
# make predictor map stack to feed to prediction function
predictors <- stack((stem_map), (lcp_map))
predictors$max_prob <- max(predictors[[1:2]], na.rm=T) 
predictors$bird_id <- predictors$stem_prob
values(predictors$bird_id) <- dat1_train$bird_id[1]
# plot(predictors)
  
# make prediction map stack
r3_pred <- predict(predictors, r3, exclude="s(bird_id)")
r4_pred <- predict(predictors, r4, exclude="s(bird_id)")
pred_stack <- stack(r3_pred, r4_pred)
names(pred_stack) <- c("r3_pred", "r4_pred")
for(s in 1:nlayers(pred_stack)){
  ra <- calc(pred_stack[[s]], plogis)
  rb <- rescale_ras_01(ra)
  pred_stack[[s]] <- rb
}
names(pred_stack) <- c("r3_pred", "r4_pred")
crs(pred_stack) <- mbi_analysis
terra::writeRaster(pred_stack$r3_pred, 
                   paste0("./output/", species_name, "_", 
                          focal_season, "_", 
                          "stem_lcp_integrated_surface.tif"), 
                   overwrite=T)
terra::writeRaster(pred_stack$r4_pred, paste0("./output/", 
                                              species_name, "_", focal_season, 
                                              "_", 
                                              "max_integrated_surface.tif"), 
                   overwrite=T)
# plot(pred_stack)
  
# make df for ggplots
df_lcp <- gplot_data((lcp_map), maxpixels=ncell((lcp_map)))
df_stem <- gplot_data((stem_map), maxpixels=ncell((stem_map)))
df_max <- gplot_data(predictors$max_prob, maxpixels=ncell(predictors$max_prob))
df_r3 <- gplot_data(pred_stack$r3_pred, 
                      maxpixels=ncell(pred_stack$r3_pred))
df_r4 <- gplot_data(pred_stack$r4_pred, 
                      maxpixels=ncell(pred_stack$r4_pred))
  
# make maps
p_lcp <- ggplot() +
  geom_tile(data=filter(df_lcp, !is.na(value), value>=0.01),
            aes(x=x, y=y, fill=value)) +
  scale_fill_gradientn("Spring\nLCP\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  guides(colour = "none") +
  geom_sf(data=west_hem, fill=NA) +
  ylim(range(st_coordinates(envelope)[,"Y"])) +
  xlim(range(st_coordinates(envelope)[,"X"])) +
  theme_map()
p_stem <- ggplot() +
  geom_tile(data=filter(df_stem, !is.na(value), value>=0.01),
            aes(x=x, y=y, fill=value)) +
  scale_fill_gradientn("Spring\nSTEM\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  guides(colour = "none") +
  geom_sf(data=west_hem, fill=NA) +
  ylim(range(st_coordinates(envelope)[,"Y"])) +
  xlim(range(st_coordinates(envelope)[,"X"])) +
  theme_map()
p_max <- ggplot() +
  geom_tile(data=filter(df_max, !is.na(value), value>=0.01),
            aes(x=x, y=y, fill=value)) +
  scale_fill_gradientn("Spring\nMax\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  guides(colour = "none") +
  geom_sf(data=west_hem, fill=NA) +
  ylim(range(st_coordinates(envelope)[,"Y"])) +
  xlim(range(st_coordinates(envelope)[,"X"])) +
  theme_map()
p_r3 <- ggplot() +
  geom_tile(data=filter(df_r3, !is.na(value), value>=0.01), aes(x=x, y=y, 
                                                              fill=value)) +
  scale_fill_gradientn("Spring\nmodel\npredicted\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  scale_color_brewer("Data\ntype", palette="Set1") +
  geom_sf(data=west_hem, fill=NA) +
  ylim(range(st_coordinates(envelope)[,"Y"])) + 
  xlim(range(st_coordinates(envelope)[,"X"])) +
  theme_map()
p_r4 <- ggplot() +
  geom_tile(data=filter(df_r4, !is.na(value), value>=0.01), aes(x=x, y=y, 
                                                              fill=value)) +
  scale_fill_gradientn("Max\nmodel\npredicted\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  scale_color_brewer("Data\ntype", palette="Set1") +
  geom_sf(data=west_hem, fill=NA) +
  ylim(range(st_coordinates(envelope)[,"Y"])) + 
  xlim(range(st_coordinates(envelope)[,"X"])) +
  theme_map()
p_pts <- ggplot() +
  geom_sf(data=west_hem, fill="gray60") +
  geom_sf(data=filter(dat1_train, presence==1) %>%
            dplyr::mutate(tech_type=str_replace(tech_type, "bnd", "Band")) %>%
            dplyr::mutate(tech_type=str_replace(tech_type, "ptt", "PTT")) %>%
            dplyr::mutate(tech_type=str_replace(tech_type, "gps", "GPS")) %>%
            dplyr::mutate(tech_type=str_replace(tech_type, "llg", "LLG")),
          size=0.5, shape=16, aes(col=factor(tech_type))) +
  scale_color_brewer("Spring\nmovement\ndata", palette="Set1") +
  ylim(range(st_coordinates(envelope)[,"Y"])) + 
  xlim(range(st_coordinates(envelope)[,"X"])) +
  theme_map()
pg2 <- plot_grid(p_stem, p_lcp, p_max, p_pts, p_par, p_r3, nrow=2); pg2
ggsave2(filename = paste0("./output/", species_name,
                  "_", focal_season, "_", "prediction_maps.pdf"),  
        width=12, height=7)

# stash it
out_list$prediction_maps_spring <- pg2
# ------------------------------------------------------------------------------




# analyze fall data ------------------------------------------------------------
# choose fall STEM and LCP
focal_season <- "postbreeding_migration"
lcp_map <- out_list$lcp_map_rs[[focal_season]]; names(lcp_map) <- "lcp_prob"
stem_map <- out_list$stem_map_rs[[focal_season]]; names(stem_map) <- "stem_prob"

# extract and prepare data for that season
dat1_train <- all_movement_df %>% 
  dplyr::filter(season==focal_season) %>% 
  st_as_sf(sf_column_name="geometry", crs=mbi_analysis) %>%
  dplyr::mutate(bird_id=factor(as.numeric(factor(bird_id)))) %>%
  dplyr::arrange(tech_type, bird_id, date_time)

# add a max predictor value
dat1_train$max_prob <-
  apply(dat1_train[,c("lcp_prob", "stem_prob")] %>%
          st_drop_geometry() %>% as.matrix(), 1,
        function(x) max(x, na.rm=T))

# models
f0 <- presence ~ 1 
r0 <- scam(f0, data=dat1_train, family="binomial")
f3 <- presence ~ 1 + s(stem_prob, k=5, bs="mpi") + s(lcp_prob, k=5, bs="mpi") 
r3 <- scam(f3, data=dat1_train, family="binomial")
f4 <- presence ~ 1 + s(max_prob, k=5, bs="mpi") 
r4 <- scam(f4, data=dat1_train, family="binomial")

# make a model fit summaries table
aic_tab <- AIC(r0, r3, r4) %>%
  mutate(Model=c("GAM_Null", "GAM_STEM_LCP", "GAM_Max"),
          D_explained=c(summary(r0)$dev.expl,
          summary(r3)$dev.expl,
          summary(r4)$dev.expl)) %>%
  arrange(AIC) %>% mutate(Delta_AIC=c(0, diff(AIC))) %>%
  mutate_if(is.numeric, round, 2) %>%
  select(Model, D_explained, AIC, Delta_AIC)

# partial plot for full GAm model with STEM and LCP
p1 <- function() {
  par(mar = c(3, 3, 1, 1), mgp = c(2, 1, 0), mfrow = c(1,2))
  plot(r3, select=1, shade=TRUE, xlab="STEM Probability", 
       ylab="Standardized effect")
  plot(r3, select=2, shade=TRUE, xlab="LCP Probability", 
       ylab="Standardized effect")
}
p_par <- ggdraw(p1)

# get model fitted values for more fit stats
dat1_train$r0_fit <- predict(r0, dat1_train, type="response", 
                             exclude="s(bird_id)")
dat1_train$r3_fit <- predict(r3, dat1_train, type="response", 
                             exclude="s(bird_id)")
dat1_train$r4_fit <- predict(r4, dat1_train, type="response", 
                             exclude="s(bird_id)")

# get auc and make ROC plot
roc0 <- roc(presence~r0_fit, data=dat1_train) 
roc3 <- roc(presence~r3_fit, data=dat1_train)
roc4 <- roc(presence~r4_fit, data=dat1_train)
p_roc <- ggroc(list(roc3, roc4), legacy.axes = T) + 
    geom_abline(intercept=0, slope=1, lty=2, col="gray60") +
    labs(x="1 - Specificity", y="Sensitivity", color="Model") +
    scale_color_brewer(labels=c("STEM+LCP\nmodel", "Max\nmodel"), 
                       palette="Set1") +
  annotate("text", x = 1, y = 0.1, hjust = 1,
            label = paste("STEM+LCP model AUC =",
                          str_pad(round(roc3$auc, 2), width=4,
                                  side ="right", pad=0))) +
  annotate("text", x = 1, y = 0, hjust = 1,
            label = paste("Max model AUC =", 
                          str_pad(round(roc4$auc, 2), width=4, 
                                  side ="right", pad=0)))
# p_roc
  
# add auc to fit table and save
aic_tab$AUC <- NA
aic_tab[aic_tab$Model=="GAM_STEM_LCP", "AUC"] <- round(roc3$auc, 2)
aic_tab[aic_tab$Model=="GAM_Max", "AUC"] <- round(roc4$auc, 2)
row.names(aic_tab) <- NULL
write.csv(aic_tab, paste0("./output/", species_name,
                  "_", focal_season, "_", "aic_table.csv"), row.names=F, na="")
  
# make predictor map stack to feed to prediction function
predictors <- stack(stem_map, lcp_map)
predictors$max_prob <- max(predictors[[1:2]], na.rm=T) # product to fall back on, if necessary
predictors$bird_id <- predictors$stem_prob
values(predictors$bird_id) <- dat1_train$bird_id[1]
# plot(predictors)
  
# make prediction map stack
r3_pred <- predict(predictors, r3, exclude="s(bird_id)")
r4_pred <- predict(predictors, r4, exclude="s(bird_id)")
pred_stack <- stack(r3_pred, r4_pred)
names(pred_stack) <- c("r3_pred", "r4_pred")
for(s in 1:nlayers(pred_stack)){
  ra <- calc(pred_stack[[s]], plogis)
  rb <- rescale_ras_01(ra)
  pred_stack[[s]] <- rb
}
names(pred_stack) <- c("r3_pred", "r4_pred")
crs(pred_stack) <- mbi_analysis
terra::writeRaster(pred_stack$r3_pred, paste0("./output/", species_name,
                  "_", focal_season, "_", "stem_lcp_integrated_surface.tif"), overwrite=T)
terra::writeRaster(pred_stack$r4_pred, paste0("./output/", species_name,
                  "_", focal_season, "_", "max_integrated_surface.tif"), overwrite=T)
# plot(pred_stack)
  
# make df for ggplots
df_lcp <- gplot_data((lcp_map), maxpixels=ncell((lcp_map)))
df_stem <- gplot_data((stem_map), maxpixels=ncell((stem_map)))
df_max <- gplot_data(predictors$max_prob, maxpixels=ncell(predictors$max_prob))
df_r3 <- gplot_data(pred_stack$r3_pred, 
                      maxpixels=ncell(pred_stack$r3_pred))
df_r4 <- gplot_data(pred_stack$r4_pred, 
                      maxpixels=ncell(pred_stack$r4_pred))
  
# make maps
p_lcp <- ggplot() +
  geom_tile(data=filter(df_lcp, !is.na(value), value>=0.01),
            aes(x=x, y=y, fill=value)) +
  scale_fill_gradientn("Fall\nLCP\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  guides(colour = "none") +
  geom_sf(data=west_hem, fill=NA) +
  ylim(range(st_coordinates(envelope)[,"Y"])) +
  xlim(range(st_coordinates(envelope)[,"X"])) +
  theme_map()
p_stem <- ggplot() +
  geom_tile(data=filter(df_stem, !is.na(value), value>=0.01),
            aes(x=x, y=y, fill=value)) +
  scale_fill_gradientn("Fall\nSTEM\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  guides(colour = "none") +
  geom_sf(data=west_hem, fill=NA) +
  ylim(range(st_coordinates(envelope)[,"Y"])) +
  xlim(range(st_coordinates(envelope)[,"X"])) +
  theme_map()
p_max <- ggplot() +
  geom_tile(data=filter(df_max, !is.na(value), value>=0.01),
            aes(x=x, y=y, fill=value)) +
  scale_fill_gradientn("Fall\nMax\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  guides(colour = "none") +
  geom_sf(data=west_hem, fill=NA) +
  ylim(range(st_coordinates(envelope)[,"Y"])) +
  xlim(range(st_coordinates(envelope)[,"X"])) +
  theme_map()
p_r3 <- ggplot() +
  geom_tile(data=filter(df_r3, !is.na(value), value>=0.01), aes(x=x, y=y, 
                                                              fill=value)) +
  scale_fill_gradientn("Fall\nmodel\npredicted\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  scale_color_brewer("Data\ntype", palette="Set1") +
  geom_sf(data=west_hem, fill=NA) +
  ylim(range(st_coordinates(envelope)[,"Y"])) +
  xlim(range(st_coordinates(envelope)[,"X"])) +
  theme_map()
p_r4 <- ggplot() +
  geom_tile(data=filter(df_r4, !is.na(value), value>=0.01), aes(x=x, y=y, 
                                                              fill=value)) +
  scale_fill_gradientn("Max\nmodel\npredicted\nindex", 
                        colours = rev(terrain.colors(10)),
                        limits=c(0,1)) +
  scale_color_brewer("Data\ntype", palette="Set1") +
  geom_sf(data=west_hem, fill=NA) +
  ylim(range(st_coordinates(envelope)[,"Y"])) +
  xlim(range(st_coordinates(envelope)[,"X"])) +
  theme_map()
p_pts <- ggplot() +
  geom_sf(data=west_hem, fill="gray60") +
  geom_sf(data=filter(dat1_train, presence==1) %>%
            dplyr::mutate(tech_type=str_replace(tech_type, "bnd", "Band")) %>%
            dplyr::mutate(tech_type=str_replace(tech_type, "ptt", "PTT")) %>%
            dplyr::mutate(tech_type=str_replace(tech_type, "gps", "GPS")) %>%
            dplyr::mutate(tech_type=str_replace(tech_type, "llg", "LLG")),
          size=0.5, shape=16, aes(col=factor(tech_type))) +
  scale_color_brewer("Fall\nmovement\ndata", palette="Set1") +
  ylim(range(st_coordinates(envelope)[,"Y"])) + 
  xlim(range(st_coordinates(envelope)[,"X"])) +
  theme_map()
pg3 <- plot_grid(p_stem, p_lcp, p_max, p_pts, p_par, p_r3, nrow=2); pg3
ggsave2(filename = paste0("./output/", species_name,
                  "_", focal_season, "_", "prediction_maps.pdf"),  
        width=12, height=7)

# stash it
out_list$prediction_maps_fall <- pg3
save(out_list, file=paste0(out_path, species_name, "/", species_name,
                    "_", "analysis_data.RData"))

# finish species
print(paste("Done with", spp1))
# ------------------------------------------------------------------------------




# end for loop -----------------------------------------------------------------  
# }
# ------------------------------------------------------------------------------







