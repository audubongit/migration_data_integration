# LCP code
# Tim Meehan
# 28 Sep 2021


# setup ------------------------------------------------------------------------
# time zone
Sys.setenv(TZ="America/Denver") 

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

# ebird key
set_ebirdst_access_key("pquolo3jdeur", overwrite=T)

# settings
extract <- raster::extract
select <- dplyr::select
options(scipen=9999999)

# directories
setwd(getwd())
out_path <- "./outputs/"
movedat_path <- "./inputs/"
mcr_path <- "./inputs/"
mcmat_path <- "./inputs/"

# crs
mbi_analysis <- "+proj=aeqd +lat_0=15 +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
ebird_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# helper functions
# get cores function
get_cores <- function(ras=as_breed,
                      pols=breeding_mcrs,
                      keeps=0.30) {
  cr <- crs(ras)
  rs <- res(ras)
  ex <- extent(ras)
  outs <- c()
  nams <- unique(pols$mcr_class)
  for(i in 1:length(nams)){
    bp_i <- as.data.frame(rasterToPoints(mask(ras,
                                    pols[i,]))) %>%
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

# get path function
get_path <- function(wp = all_win_points, sp = all_sum_points, 
                   cond = spring_conduct, j = i, seas="prebreeding"){
  require("raster")
  require("gdistance")
  require("sf")
  o <- as(wp[j,], "Spatial")
  g <- as(sp[j,], "Spatial")
  t <- passage(x=cond, 
               origin=o,
               goal=g, 
               theta=cellStats(raster(spring_conduct), min, na.rm=T)/2,
               output="RasterLayer")
  names(t) <- paste(seas, paste("start", o$mcr_class, sep="."),
        paste("stop", g$mcr_class, sep="."), sep=".")
  return(t)
}

# from grid to df for ggplot
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

# rescale functions
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
# baseline function
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

# ggplot themes
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
east_hem <- read_sf("./inputs/east_hemi_ebird.shp") %>% 
  st_transform(mbi_analysis) %>%
  st_simplify(dTolerance = 5000) %>% 
  dplyr::select(1) %>% 
  mutate(mask=1) %>% 
  select(mask)
# check it out: plot(east_hem$geometry)
# ------------------------------------------------------------------------------



# load species list and specify number of lcp point pairs per mcr --------------
spp_table <- read.csv("./inputs/species_table.csv")
spec_list <- spp_table %>% pull(species) %>% as.character()
overwater_list <- spp_table %>% pull(overwater)
overland_list <- spp_table %>% pull(overland)
model_int <- spp_table %>% pull(model_mig_integr_done)
max_int <- spp_table %>% pull(max_mig_integr_done)

# make sure eBird ST surfaces exist for your species
spec_list %in% ebirdst_runs$common_name

# set points
n_connect_points <- 5
# ------------------------------------------------------------------------------



# loop through species, in this example a single species -----------------------
# it is set up as a loop so a multispecies analysis can be run as a paralleled 
# loop on multiple cores
s <- 1
for(s in 1:length(1)){
  # try catch
  tryCatch({
  
    # set species pars
  spp1 <- spec_list[s]
  overwater <- overwater_list[s]
  overland <- overland_list[s]
  spp2 <- ebirdst_runs$species_code[ebirdst_runs$common_name==spp1]
  model_int_s <- model_int[s]
  max_int_s <- max_int[s]
  print(paste("Starting", spp1))
  
  
  # get migration occurrence surfaces ------------------------------------------
  # get raster
  dl_path <- ebirdst_download(species = spp2, force=T)
  occs <- load_raster("occurrence", path = dl_path, resolution="lr")

  # get occ seasons
  sp_dates <- dplyr::filter(ebirdst_runs, common_name == spp1) %>%
    dplyr::select(setdiff(matches("(start)|(end)"), matches("year_round"))) %>%
    gather("label", "date") %>%
    separate(label, c("season", "start_end"), "_(?=s|e)") %>%
    spread(start_end, date) %>%
    dplyr::select(season, start_dt=start, end_dt=end) %>%
    filter(season != "resident")
  
  # make sure seasons pass qaqc
  run_review <- dplyr::filter(ebirdst_runs, common_name == spp1)
  quality_rating <- run_review[paste0(sp_dates$season, "_quality")]
  sp_dates$pass <- as.integer(quality_rating) > 1
  sp_dates <- mutate(sp_dates, pass = !(is.na(start_dt) | is.na(end_dt)))
  
  # get occs weeks
  weeks <- parse_raster_dates(label_raster_stack(occs))
  weeks_season <- rep(NA_character_, length(weeks))
  for (i in seq_len(nrow(sp_dates))) {
    s <- sp_dates[i, ]
    if (!s$pass) {next()}
    if (s$start_dt <= s$end_dt) {
      in_season <- weeks >= s$start_dt & weeks <= s$end_dt
    } else {
      in_season <- weeks >= s$start_dt | weeks <= s$end_dt
    }
    weeks_season[in_season] <- s$season
  }
  # check it out: table(weeks_season)

  # drop weeks not assigned to season
  week_pass1 <- !is.na(weeks_season)
  occs <- occs[[which(week_pass1)]]
  weeks1 <- weeks[week_pass1]
  weeks_season1 <- weeks_season[week_pass1]

  # save spring layers
  week_pass4 <- weeks_season1=="prebreeding_migration"
  occ_spring <- occs[[which(week_pass4)]]
  weeks4 <- weeks1[week_pass4]
  weeks_season4 <- weeks_season1[week_pass4]
  
  # save fall layers
  week_pass5 <- weeks_season1=="postbreeding_migration"
  occ_fall <- occs[[which(week_pass5)]]
  weeks5 <- weeks1[week_pass5]
  weeks_season5 <- weeks_season1[week_pass5]
  
  # combine, trim, average, and reproject migration occ stack
  occs_stack <- stack(mean(occ_spring, na.rm=T), mean(occ_fall, na.rm=T)) %>% 
    projectRaster(crs=mbi_analysis) %>%
    trim(.)
  names(occs_stack) <- c("prebreeding_migration_occurrence", 
                         "postbreeding_migration_occurrence")
  # check it out: plot(occs_stack)
  # ----------------------------------------------------------------------------
    
  
  # get stationary abundance surfaces ------------------------------------------
  # stack abund stems
  abd <- load_raster("abundance_seasonal", path = dl_path, resolution="lr")
  abund_stack <- stack(abd[["breeding"]], abd[["nonbreeding"]])
  abund_stack <- projectRaster(abund_stack, crs=mbi_analysis)
  names(abund_stack) <- c("breeding_abundance", "nonbreeding_abundance")
  # check it out: plot(abund_stack)
  
  # backup occ and abund stacks
  dir.create(paste0(out_path, gsub(" ", "_", spp1)))
  writeRaster(occs_stack, paste0(out_path, gsub(" ", "_", spp1), "/", 
                          gsub(" ", "_", spp1), 
                          "_migration_occurrence_stack_27km.grd"),
              format="raster", overwrite=T)
  writeRaster(abund_stack, paste0(out_path, gsub(" ", "_", spp1), "/", 
                          gsub(" ", "_", spp1),
                          "_stationary_abundance_stack_27km.grd"),
              format="raster", overwrite=T)
  # ----------------------------------------------------------------------------
  
  
  # create land mask for baselining --------------------------------------------
  land_mask <- rasterize(east_hem, occs_stack)
  
  # make migration stack and add baseline conductance values set above
  mig_stack <- stack(occs_stack$prebreeding_migration_occurrence,
                     occs_stack$postbreeding_migration_occurrence)
  mig_stack$prebreeding_migration_occurrence <- 
    add_baselines(x=mig_stack$prebreeding_migration_occurrence,
                  lmask=land_mask, 
                  base_land=overland,
                  base_ocean=overwater)  
  mig_stack$postbreeding_migration_occurrence <- 
    add_baselines(x=mig_stack$postbreeding_migration_occurrence, 
                  lmask=land_mask, 
                  base_land=overland,
                  base_ocean=overwater)
  # check it out: plot(mig_stack)
  # ----------------------------------------------------------------------------
  
  
  # get mcrs and select breeding abundance cores per mcr -----------------------
  mcrs <- read_sf(paste0(mcr_path, spp2, "_bird_migration_regions.shp")) %>%
    st_transform(mbi_analysis) %>%
    st_simplify(dTolerance = 5000) 
  breeding_mcrs <- mcrs %>% filter(season=="breeding") 
  nonbreeding_mcrs <- mcrs %>% filter(season=="nonbreeding") 
  # check it out: plot(mcrs$geometry)
  
  print("got all maps")
  # ----------------------------------------------------------------------------
  

  # select summer and winter core areas per mcr --------------------------------
  # get rid of zero abundance for speed
  as_breed <- abund_stack$breeding_abundance
  as_nonbreed <- abund_stack$nonbreeding_abundance
  as_breed[as_breed<=0] <- NA
  as_nonbreed[as_nonbreed<=0] <- NA
  
  # run get_cores function
  breed_cores <- get_cores(ras=as_breed, pols=breeding_mcrs, keeps=0.30)
  nonbreed_cores <- get_cores(ras=as_nonbreed, pols=nonbreeding_mcrs, 
                              keeps=0.30)
  
  # make core rasters
  breed_core_ras <- rasterFromXYZ(breed_cores[, c(1,2,3)], crs=crs(as_breed), 
                                  res=res(as_breed))
  nonbreed_core_ras <- rasterFromXYZ(nonbreed_cores[, c(1,2,3)], 
                                     crs=crs(as_breed), res=res(as_breed))
  # check it out: plot(breed_core_ras); plot(nonbreed_core_ras)
  # ----------------------------------------------------------------------------
  
  
  # select breeding points within cores and match with winter points -----------
  breeding_points <- breed_cores %>% group_by(mcr_class) %>% 
    sample_n(n_connect_points, replace=T) %>% dplyr::select(1,2,4) %>% ungroup()
  # check it out: plot(st_geometry(breeding_mcrs), col="gray")
  # check it out: plot(st_as_sf(breeding_points, coords=c("x","y")), add=T)
  
  # list breeding mcrs for species
  breeding_mcr_list <- sort(mcrs$mcr_class[mcrs$season=="breeding"])
  nonbreeding_mcr_list <- sort(mcrs$mcr_class[mcrs$season=="nonbreeding"])
  
  # pull in the connectivity matrix if there is one. if not, then script assumes
  # equal connectivity across mcrs
  cm1 <- paste0(mcmat_path, spp2, "_connect_matrix.csv")
  if(length(cm1) > 0) { 
    mcmat <- read.csv(cm1) %>% select(proportion=1,mcr_class_breeding=7,
                                      mcr_class_nonbreeding=8)
  } else mcmat <- expand.grid(proportion=1/length(unique(nonbreeding_mcr_list)),
              mcr_class_breeding=breeding_mcr_list, 
              mcr_class_nonbreeding=nonbreeding_mcr_list)

  # match winter and summer points per breeding mcr based on connectivity matrix
  all_sum_points <- c()
  all_win_points <- c()
  for(m in 1:length(breeding_mcr_list)){
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
  #-----------------------------------------------------------------------------
    
  
  # make conductance surfaces from migration occs surfaces ---------------------
  spring_conduct <- transition(mig_stack$prebreeding_migration_occurrence, 
                             transitionFunction=mean, directions=8) %>%
    geoCorrection(scl=T)
  fall_conduct <- transition(mig_stack$postbreeding_migration_occurrence, 
                             transitionFunction=mean, directions=8) %>%
    geoCorrection(scl=T)  
  # check it out: plot(raster(fall_conduct))
  # check it out: plot(all_win_points$geometry, add=T, pch=16, col="blue")
  # check it out: plot(all_sum_points$geometry, add=T, pch=16, col="red")
  # ----------------------------------------------------------------------------
  
  
  # get probability of spring passage ------------------------------------------
  print("getting spring corridors")
  
  # make cluster
  ##########  WARNING, ADJUST AS NECESSARY TO NOT CRASH YOUR COMPUTER ##########
  ########## If looping through species on a parallel loop, maybe set ##########
  ########## this to a value = 1.                                     ##########
  ########## If doing this one species at a time, maybe set as high   ##########
  ########## as your machine can handle (number cores - 2)            ##########
  cl <- makeCluster(5)
  registerDoParallel(cl)
  prob_stack <- foreach(i=1:nrow(all_sum_points), 
                      .packages=c("raster", "gdistance",
                                                      "sf")) %dopar% {
    get_path(wp=all_win_points, sp=all_sum_points, 
             cond=spring_conduct, j=i, seas="prebreeding")
  }
  stopCluster(cl)
  
  # make stack from parallel list output list and add names
  prob_stack <- stack(prob_stack)
  # check out individual paths: plot(prob_stack[[c(1:3, 7:9)]])

  # summary layer
  prebreeding_migration_mean <- rescale_it(sqrt(mean(prob_stack)))
  names(prebreeding_migration_mean) <- "prebreeding.migration.mean"
  # check out squashed paths: plot(((prebreeding_migration_mean)))
  # check out squashed paths: plot(st_geometry(east_hem), add=T)
  # ----------------------------------------------------------------------------
  
  
  # get probability of fall passage --------------------------------------------
  print("getting fall corridors")
  cl <- makeCluster(5)
  registerDoParallel(cl)
  prob_stack2 <- foreach(i=1:nrow(all_sum_points), 
                      .packages=c("raster", "gdistance",
                                                      "sf")) %dopar% {
    get_path(wp=all_sum_points, sp=all_win_points, 
             cond=fall_conduct, j=i, seas="postbreeding")
  }
  stopCluster(cl)
  
  # make stack from parallel list output list and add names
  prob_stack2 <- stack(prob_stack2)
  # check out individual paths: plot(prob_stack2[[c(1:3, 7:9)]])
  
  # summary layer
  postbreeding_migration_mean <- rescale_it(sqrt(mean(prob_stack2)))
  names(postbreeding_migration_mean) <- "postbreeding.migration.mean"
  # check out squashed paths: plot(((postbreeding_migration_mean)))
  # check out squashed paths: plot(st_geometry(east_hem), add=T)
  # ----------------------------------------------------------------------------
  
  
  # get paths per spring migration starting mcr --------------------------------
  prob_stack_targets <- unique(unlist(lapply(str_split(names(prob_stack), 
                        pattern="[.]", n=4, simplify=F),
         function(x) {paste(x[1], x[2], x[3], sep=".")})))
  
  # get mean for first starter mcr
  means_prebreeding <- rescale_it(sqrt(mean(prob_stack[[grep(prob_stack_targets[1], 
                            names(prob_stack), value = T)]], na.rm=T)))
  names(means_prebreeding) <- prob_stack_targets[1]

  # add means for other starter mcrs
  for(pt in 2:length(prob_stack_targets)){
    mean_i <- rescale_it(sqrt(mean(prob_stack[[grep(prob_stack_targets[pt], 
                            names(prob_stack), value = T)]], na.rm=T)))
    names(mean_i) <- prob_stack_targets[pt]
    means_prebreeding <- addLayer(means_prebreeding, mean_i)
  }
  # check it out: plot(means_prebreeding)
  # ----------------------------------------------------------------------------

  
  # get paths per fall migration starting mcr ----------------------------------
  prob_stack2_targets <- unique(unlist(lapply(str_split(names(prob_stack2), 
                        pattern="[.]", n=4, simplify=F),
         function(x) {paste(x[1], x[2], x[3], sep=".")})))
  
  # get mean for first starter mcr
  means_postbreeding <- rescale_it(sqrt(mean(prob_stack2[[grep(prob_stack2_targets[1], 
                            names(prob_stack2), value = T)]], na.rm=T)))
  names(means_postbreeding) <- prob_stack2_targets[1]
  
  # add means for other starter mcrs
  for(pt in 2:length(prob_stack2_targets)){
    mean_i <- rescale_it(sqrt(mean(prob_stack2[[grep(prob_stack2_targets[pt], 
                            names(prob_stack2), value = T)]], na.rm=T)))
    names(mean_i) <- prob_stack2_targets[pt]
    means_postbreeding <- addLayer(means_postbreeding, mean_i)
  }
  # check it out: plot(means_postbreeding)
  # ----------------------------------------------------------------------------
  
    
  # stack and save all lcps ----------------------------------------------------
  lcp_stack <- stack(prebreeding_migration_mean,
                     postbreeding_migration_mean,
                     means_prebreeding,
                     means_postbreeding)
  # save lcp stack
  writeRaster(lcp_stack, paste0(out_path, gsub(" ", "_", spp1), "/", 
                          gsub(" ", "_", spp1), "_",
                                 "lcp_stack_27km.grd"), format="raster", 
              overwrite=T)
  # ----------------------------------------------------------------------------

  
  # clean up and dump garbage --------------------------------------------------
  rm(prob_stack2, prob_stack,
     means_postbreeding, means_prebreeding,
     spring_conduct, fall_conduct,
     as_nonbreed, as_breed)
  gc()
  # ----------------------------------------------------------------------------
  
  
  # set parameters for integration and validation ------------------------------
  species_code <- s1 <- spp2
  species_name <- s2 <- gsub(" ", "_", spp1)
  mod_int_i <- model_int_s 
  max_int_i <- max_int_s
  # number of background points to compare with presence-only points
  n_background <- 10000
  # keep only one tracking point per day?
  thin_per_day <- F
  # keep only one tracking point per 27 km grid cell?
  thin_per_cell <- T
  # keep only this many band encounters
  max_bands <- 5000
  
  # load more libs for integration and validation
  library(cowplot)
  library(dismo)
  library(lubridate)
  library(mgcViz)
  library(gammit)
  library(mgcv)
  library(pROC)
  
  # more options
  options(stringsAsFactors = FALSE)
  # ----------------------------------------------------------------------------
  
  
  # get and trim lcp and stem layers -------------------------------------------
  spec1 <- species_code
  spec2 <- species_name
  
  # save stuff to output list
  out_list <- list()
  out_list$species_code <- species_code
  out_list$species_name <- species_name
  
  # save hemisphere map
  out_list$hemisphere_map_sf <- east_hem
  hem_poly <- st_convex_hull(east_hem)
  out_list$hemisphere_poly_sf <- hem_poly
  
  # get extent of integration/validation area, where background points 
  # could be located
  envelope <- read_sf(paste0(movedat_path, "/", 
                     list.files(movedat_path, 
                                pattern="envelope.shp"))) %>%
    st_transform(mbi_analysis) %>%
    st_simplify(dTolerance = 5000)
  envelope <- st_intersection(envelope, hem_poly)
  out_list$study_region_map_sf <- envelope
  # check out extent: plot(envelope$geometry)
  print("got extent maps")
  
  # trim lcp maps
  lcp_map <- lcp_stack[[1:2]]
  lcp_map <- raster::crop(lcp_map, envelope) %>% 
    raster::mask(envelope)
  names(lcp_map) <- c("prebreeding_migration", "postbreeding_migration")
  out_list$lcp_map_rs <- lcp_map
  # check out trimmed lcp maps: plot(lcp_map)
  print("got trimmed lcp maps")
  
  # get and trim eBird S&T maps
  stem_map <- occs_stack[[1:2]]
  stem_map <- raster::resample(stem_map, lcp_map)
  stem_map[is.na(stem_map)] <- 0
  stem_map <- raster::crop(stem_map, envelope) %>% 
    raster::mask(envelope)
  names(stem_map) <- c("prebreeding_migration", "postbreeding_migration")
  out_list$stem_map_rs <- stem_map
  # check out trimmed s&t maps: plot(stem_map)
  print("got trimmed eBird S&T occurrence maps")
  # ----------------------------------------------------------------------------
  
  
  # get, clean, and sample movement data ---------------------------------------
  move_dat <- read.csv(paste0(movedat_path, spec1, "_movement_data.csv")) %>%
    mutate(season=ifelse(season=="summer", "breeding", season)) %>%
    mutate(season=ifelse(season=="winter", "nonbreeding", season)) %>%
    mutate(season=ifelse(season=="spring", "prebreeding_migration", season)) %>%
    mutate(season=ifelse(season=="fall", "postbreeding_migration", season)) %>%
    dplyr::filter(!grepl("rng", outlier),
           !grepl("ang", outlier),
           !grepl("dst", outlier),
           !grepl("mrk", outlier),
           !grepl("ddd", outlier),
           !grepl("dup", outlier),
           !grepl("org", outlier),
           tech_type %in% c("bnd", "ptt", "gps", "llg"),
           season!="breeding",
           season!="nonbreeding") %>%
    dplyr::mutate(date_time=lubridate::ymd_hms(date_time)) %>%
    dplyr::mutate(date_time=lubridate::date(date_time)) %>%
    st_as_sf(coords=c("x_coord","y_coord"), crs = 4236) %>%
    st_transform(crs=mbi_analysis) %>% 
    st_join(envelope, join = st_within)
  out_list$initial_move_dat_ss <- nrow(move_dat)
  # check it out: plot(envelope$geometry); plot(move_dat["tech_type"], add=T)
  
  # remove geolocator data near equinox
  move_dat <- move_dat %>% filter(!((week(move_dat$date_time)==10 | 
         week(move_dat$date_time)==11 |
         week(move_dat$date_time)==12 |
         week(move_dat$date_time)==37 |
         week(move_dat$date_time)==38 |
         week(move_dat$date_time)==39) & move_dat$tech_type=="llg"))
  
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
  if(thin_per_day==T){
    move_dat <- move_dat %>%
      dplyr::group_by(bird_id, date_time, season, tech_type) %>%
      dplyr::sample_n(size = 1) %>%
      dplyr::ungroup()
  }
  out_list$after_temporal_move_dat_ss <- nrow(move_dat)
  
  # sample 1 track point per cell
  if(thin_per_cell==T){
    r <- lcp_map$prebreeding_migration
    values(r) <- 1:ncell(r)
    move_dat$grid_id <- raster::extract(r, move_dat)
    # plot(move_dat["grid_id"])
    move_dat <- move_dat %>% 
      dplyr::ungroup() %>%
      dplyr::group_by(bird_id, season, grid_id, tech_type) %>%
      dplyr::sample_n(size = 1) %>%
      dplyr::select(-FID)
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
      summarise(individual_bird_ids=length(unique(bird_id)),
                records=n()) %>% 
      arrange(desc(records))
  
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
  dat3$lcp_prob <- raster::extract(lcp_map$postbreeding_migration, 
                           as(dat3, "Spatial"))
  dat3$stem_prob <- raster::extract(stem_map$postbreeding_migration, 
                           as(dat3, "Spatial"))
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
  dat3$lcp_prob <- raster::extract(lcp_map$prebreeding_migration, 
                           as(dat3, "Spatial"))
  dat3$stem_prob <- raster::extract(stem_map$prebreeding_migration, 
                           as(dat3, "Spatial"))
  dat3 <- dat3[complete.cases(dat3$lcp_prob, dat3$stem_prob),]
  out_list$spring_df <- dat3
  
  # combine all
  all_movement_df <- rbind(out_list$spring_df, out_list$fall_df)
  out_list$all_movement_df <- all_movement_df
  View(all_movement_df)
  
  # save integration and validation data for the species
  save(out_list, file=paste0(out_path, "/", species_name, "/", species_name,
                     "_", "analysis_data.RData"))
  print("got available movement data")
  #-----------------------------------------------------------------------------
  
  
  # put integration and validation data back into env and make plots -----------
  list2env(out_list, environment())
  
  # plot it
  p1 <- all_movement_df %>%
    dplyr::filter(season=="prebreeding_migration", tech_type!="ran") %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "bnd", "Band")) %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "ptt", "PTT")) %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "gps", "GPS")) %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "llg", "LLG")) %>%
    ggplot() +
    geom_tile(data=filter(gplot_data(lcp_map_rs[["prebreeding_migration"]],
                                     maxpixels=ncell(lcp_map_rs)),
                          !is.na(value), value>=0.01), aes(x=x, y=y, fill=value)) +
    geom_sf(data=hemisphere_map_sf, fill=NA) +
    geom_sf(aes(col=tech_type), size=0.1, shape=16) +
    scale_fill_gradientn("Spring\nLCP\nindex", 
                         colours = rev(terrain.colors(10)),
                         limits=c(0,1)) +
    scale_color_brewer("Data\ntype", palette="Set1") +
    coord_sf(xlim=st_bbox(study_region_map_sf)[c(1,3)],
             ylim=st_bbox(study_region_map_sf)[c(2,4)]) +
    labs(x="Longitude", y="Latitude")
  p2 <- all_movement_df %>% dplyr::filter(season=="prebreeding_migration",
                                   tech_type!="ran") %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "bnd", "Band")) %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "ptt", "PTT")) %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "gps", "GPS")) %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "llg", "LLG")) %>%
    ggplot() +
    geom_tile(data=filter(gplot_data(stem_map_rs[["prebreeding_migration"]],
                                     maxpixels=ncell(lcp_map_rs)),
                          !is.na(value), value>=0.01), aes(x=x, y=y, fill=value)) +
    geom_sf(data=hemisphere_map_sf, fill=NA) +
    geom_sf(aes(col=tech_type), size=0.1, shape=16) +
    scale_fill_gradientn("Spring\nSTEM\nindex", 
                         colours = rev(terrain.colors(10)),
                         limits=c(0,1)) +
    scale_color_brewer("Data\ntype", palette="Set1") +
    coord_sf(xlim=st_bbox(study_region_map_sf)[c(1,3)],
             ylim=st_bbox(study_region_map_sf)[c(2,4)]) +
    labs(x="Longitude", y="Latitude")
  p3 <- all_movement_df %>% dplyr::filter(season=="postbreeding_migration",
                                   tech_type!="ran") %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "bnd", "Band")) %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "ptt", "PTT")) %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "gps", "GPS")) %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "llg", "LLG")) %>%
    ggplot() +
    geom_tile(data=filter(gplot_data(lcp_map_rs[["postbreeding_migration"]],
                                     maxpixels=ncell(lcp_map_rs)),
                          !is.na(value), value>=0.01), aes(x=x, y=y, fill=value)) +
    geom_sf(data=hemisphere_map_sf, fill=NA) +
    geom_sf(aes(col=tech_type), size=0.1, shape=16) +
    scale_fill_gradientn("Fall\nLCP\nindex", 
                         colours = rev(terrain.colors(10)),
                         limits=c(0,1)) +
    scale_color_brewer("Data\ntype", palette="Set1") +
    coord_sf(xlim=st_bbox(study_region_map_sf)[c(1,3)],
             ylim=st_bbox(study_region_map_sf)[c(2,4)]) +
    labs(x="Longitude", y="Latitude")
  p4 <- all_movement_df %>% dplyr::filter(season=="postbreeding_migration",
                                   tech_type!="ran") %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "bnd", "Band")) %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "ptt", "PTT")) %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "gps", "GPS")) %>%
    dplyr::mutate(tech_type=str_replace(tech_type, "llg", "LLG")) %>%
    ggplot() +
    geom_tile(data=filter(gplot_data(stem_map_rs[["postbreeding_migration"]],
                                     maxpixels=ncell(lcp_map_rs)),
                          !is.na(value), value>=0.01), aes(x=x, y=y, fill=value)) +
    geom_sf(data=hemisphere_map_sf, fill=NA) +
    geom_sf(aes(col=tech_type), size=0.1, shape=16) +
    scale_fill_gradientn("Fall\nSTEM\nindex", 
                         colours = rev(terrain.colors(10)),
                         limits=c(0,1)) +
    scale_color_brewer("Data\ntype", palette="Set1") +
    coord_sf(xlim=st_bbox(study_region_map_sf)[c(1,3)],
             ylim=st_bbox(study_region_map_sf)[c(2,4)]) +
    labs(x="Longitude", y="Latitude")
  plot_grid(p1, p2, p3, p4, nrow = 2)
  ggsave2(filename = paste0(out_path, "/", species_name, "/", species_name,
                    "_", "lcp_stem_and_movement_map.pdf"),  width=9.6, height=9.6)
  # ----------------------------------------------------------------------------
  
  
  # analyze spring data  -------------------------------------------------------
  # choose spring STEM and LCP
  focal_season <- "prebreeding_migration"
  lcp_map <- lcp_map_rs[[focal_season]]; names(lcp_map) <- "lcp_prob"
  stem_map <- stem_map_rs[[focal_season]]; names(stem_map) <- "stem_prob"
  
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

  # make models depending on type of data, where models with tracking data
  # have random effects for individual birds
  if(paste(unique(dat1_train$tech_type), collapse="_") == "bnd_ran"){
    f0 <- presence ~ 1
    r0 <- gam(f0, data=dat1_train, family="binomial")
    f1 <- presence ~ 1 + s(stem_prob, k=3)
    r1 <- gam(f1, data=dat1_train, family="binomial")
    f3 <- presence ~ 1 + s(stem_prob, lcp_prob, k=6)
    r3 <- gam(f3, data=dat1_train, family="binomial")
    f4 <- presence ~ 1 + s(max_prob, k=3)
    r4 <- gam(f4, data=dat1_train, family="binomial")
  }
  if(paste(unique(dat1_train$tech_type), collapse="_") != "bnd_ran"){
    f0 <- presence ~ 1 + s(bird_id, bs="re")
    r0 <- gam(f0, data=dat1_train, family="binomial")
    f1 <- presence ~ 1 + s(stem_prob, k=3) + s(bird_id, bs="re")
    r1 <- gam(f1, data=dat1_train, family="binomial")
    f3 <- presence ~ 1 + s(stem_prob, lcp_prob, k=6) + s(bird_id, bs="re")
    r3 <- gam(f3, data=dat1_train, family="binomial")
    f4 <- presence ~ 1 + s(max_prob, k=3) + s(bird_id, bs="re")
    r4 <- gam(f4, data=dat1_train, family="binomial")
  }
  # check it out: summary(r1); summary(r3); summary(r4)
  
  # make a model fit summaries table
  aic_tab <- AIC(r0, r1, r3, r4) %>% 
    mutate(Model=c("GAM_Null", "GAM_STEM", "GAM_STEM_LCP", "GAM_Max_STEM_LCP"),
            D_explained=c(summary(r0)$dev.expl,
            summary(r1)$dev.expl,
            summary(r3)$dev.expl,
            summary(r4)$dev.expl)) %>% 
    arrange(AIC) %>% mutate(Delta_AIC=c(0, diff(AIC))) %>%
    mutate_if(is.numeric, round, 2) %>%
    select(Model, D_explained, AIC, Delta_AIC)

  # partial plot for full GAm model with STEM and LCP
  b <- getViz(r3)
  p_par <- plot(sm(b, 1), trans=plogis, n=100) + l_fitRaster() + 
    l_fitContour(colour=1) + 
    labs(title = NULL, x="STEM index", y="LCP index") + 
    scale_fill_viridis_c("Full\nmodel\npredicted\nindex", limits=c(0,1))
  # check it out: p_par
  
  # get model fitted values for more fit stats
  dat1_train$bird_id <- dat1_train$bird_id[1]
  dat1_train$r1_fit <- predict_gamm(model=r1, newdata=dat1_train, re_form = NA, 
                                type="response", exclude="s(bird_id)")[,1]
  dat1_train$r3_fit <- predict_gamm(model=r3, newdata=dat1_train, re_form = NA, 
                                type="response", exclude="s(bird_id)")[,1]
  dat1_train$r4_fit <- predict_gamm(model=r4, newdata=dat1_train, re_form = NA, 
                                type="response", exclude="s(bird_id)")[,1]
  
  # get auc and make ROC plot
  roc1 <- roc(presence~r1_fit, data=dat1_train)
  roc3 <- roc(presence~r3_fit, data=dat1_train)
  roc4 <- roc(presence~r4_fit, data=dat1_train)
  p_roc <- ggroc(list(roc1, roc3), legacy.axes = T) + 
      geom_abline(intercept=0, slope=1, lty=2, col="gray60") +
      labs(x="1 - Specificity", y="Sensitivity", color="Model") +
      scale_color_brewer(labels=c("STEM\nmodel", "STEM + LCP\nmodel"), palette="Set1") +
    labs(title=NULL, x="False positive rate", y="True positive rate") +
    annotate("text", x = 1, y = 0.1, hjust = 1, 
             label = paste("STEM model AUC =", 
                           str_pad(round(roc1$auc, 2), width=4, 
                                   side ="right", pad=0))) +
    annotate("text", x = 1, y = 0, hjust = 1,
             label = paste("STEM + LCP model AUC =", 
                           str_pad(round(roc3$auc, 2), width=4, 
                                   side ="right", pad=0)))
  # check it out: p_roc
  
  # add auc to fit table and save
  aic_tab$AUC <- NA
  aic_tab[aic_tab$Model=="GAM_Max_STEM_LCP", "AUC"] <- round(roc4$auc, 2)
  aic_tab[aic_tab$Model=="GAM_STEM_LCP", "AUC"] <- round(roc3$auc, 2)
  aic_tab[aic_tab$Model=="GAM_STEM", "AUC"] <- round(roc1$auc, 2)
  row.names(aic_tab) <- NULL
  write.csv(aic_tab, paste0(out_path, "/", species_name, "/", species_name,
                    "_", focal_season, "_", "aic_table.csv"), row.names=F, na="")
  # check it out: aic_tab
  
  # make predictor map stack to feed to prediction function
  predictors <- stack(stem_map, lcp_map)
  predictors$max_prob <- max(predictors[[1:2]], na.rm=T)
  predictors$bird_id <- predictors$stem_prob
  values(predictors$bird_id) <- dat1_train$bird_id[1]
  # check it out: plot(predictors)
  
  # make prediction map stack
  r1_pred <- predict(predictors, r1, exclude="s(bird_id)")
  r3_pred <- predict(predictors, r3, exclude="s(bird_id)")
  r4_pred <- predict(predictors, r4, exclude="s(bird_id)")
  pred_stack <- stack(r1_pred, r3_pred, r4_pred)
  names(pred_stack) <- c("r1_pred", "r3_pred", "r4_pred")
  for(s in 1:nlayers(pred_stack)){
    ra <- calc(pred_stack[[s]], plogis)
    rb <- rescale_ras_01(ra)
    pred_stack[[s]] <- rb
  }
  names(pred_stack) <- c("r1_pred", "r3_pred", "r4_pred")
  writeRaster(pred_stack$r3_pred, paste0(out_path, "/", species_name, "/", species_name,
                    "_", focal_season, "_", "integrated_surface.tif"), overwrite=T)
  writeRaster(pred_stack$r4_pred, paste0(out_path, "/", species_name, "/", species_name,
                    "_", focal_season, "_", "max_surface.tif"), overwrite=T)
  # check it out: plot(pred_stack)
  
  # make df for ggplots
  df_lcp <- gplot_data(lcp_map, maxpixels=ncell(lcp_map))
  df_stem <- gplot_data(stem_map, maxpixels=ncell(stem_map))
  df_r1 <- gplot_data(pred_stack$r1_pred, 
                       maxpixels=ncell(pred_stack$r1_pred))
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
    geom_sf(data=hemisphere_map_sf, fill=NA) +
    ylim(range(st_coordinates(study_region_map_sf)[,"Y"])) +
    xlim(range(st_coordinates(study_region_map_sf)[,"X"])) +
    theme_map()
  p_stem <- ggplot() +
    geom_tile(data=filter(df_stem, !is.na(value), value>=0.01),
              aes(x=x, y=y, fill=value)) +
   scale_fill_gradientn("Spring\nSTEM\nindex", 
                         colours = rev(terrain.colors(10)),
                         limits=c(0,1)) +
    guides(colour = "none") +
    geom_sf(data=hemisphere_map_sf, fill=NA) +
    ylim(range(st_coordinates(study_region_map_sf)[,"Y"])) +
    xlim(range(st_coordinates(study_region_map_sf)[,"X"])) +
    theme_map()
  p_r3 <- ggplot() +
    geom_tile(data=filter(df_r3, !is.na(value), value>=0.01), aes(x=x, y=y, 
                                                               fill=value)) +
    scale_fill_gradientn("Spring\nmodel\npredicted\nindex", 
                         colours = rev(terrain.colors(10)),
                         limits=c(0,1)) +
    scale_color_brewer("Data\ntype", palette="Set1") +
    geom_sf(data=hemisphere_map_sf, fill=NA) +
    ylim(range(st_coordinates(study_region_map_sf)[,"Y"])) + 
    xlim(range(st_coordinates(study_region_map_sf)[,"X"])) +
    theme_map()
  p_r4 <- ggplot() +
    geom_tile(data=filter(df_r4, !is.na(value), value>=0.01), aes(x=x, y=y, 
                                                               fill=value)) +
    scale_fill_gradientn("Max\nmodel\npredicted\nindex", 
                         colours = rev(terrain.colors(10)),
                         limits=c(0,1)) +
    scale_color_brewer("Data\ntype", palette="Set1") +
    geom_sf(data=hemisphere_map_sf, fill=NA) +
    ylim(range(st_coordinates(study_region_map_sf)[,"Y"])) + 
    xlim(range(st_coordinates(study_region_map_sf)[,"X"])) +
    theme_map()
  p_pts <- ggplot() +
    geom_sf(data=hemisphere_map_sf, fill="gray60") +
    geom_sf(data=filter(dat1, presence==1) %>%
              dplyr::mutate(tech_type=str_replace(tech_type, "bnd", "Band")) %>%
              dplyr::mutate(tech_type=str_replace(tech_type, "ptt", "PTT")) %>%
              dplyr::mutate(tech_type=str_replace(tech_type, "gps", "GPS")) %>%
              dplyr::mutate(tech_type=str_replace(tech_type, "llg", "LLG")),
            size=0.5, shape=16, aes(col=factor(tech_type))) +
    scale_color_brewer("Spring\nmovement\ndata", palette="Set1") +
    ylim(range(st_coordinates(study_region_map_sf)[,"Y"])) + 
    xlim(range(st_coordinates(study_region_map_sf)[,"X"])) +
    theme_map()
  plot_grid(p_stem, p_lcp, p_pts, p_par$ggObj, p_roc, p_r3, nrow=2)
  ggsave2(filename = paste0(out_path, "/", species_name, "/", species_name,
                    "_", focal_season, "_", "prediction_maps.pdf"),  
          width=12, height=7)
  # ----------------------------------------------------------------------------
  
  
  # analyze fall data  ---------------------------------------------------------
  # choose fall STEM and LCP
  focal_season <- "postbreeding_migration"
  lcp_map <- lcp_map_rs[[focal_season]]; names(lcp_map) <- "lcp_prob"
  stem_map <- stem_map_rs[[focal_season]]; names(stem_map) <- "stem_prob"
  
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

  # make models
  if(paste(unique(dat1_train$tech_type), collapse="_") == "bnd_ran"){
    f0 <- presence ~ 1
    r0 <- gam(f0, data=dat1_train, family="binomial")
    f1 <- presence ~ 1 + s(stem_prob, k=3)
    r1 <- gam(f1, data=dat1_train, family="binomial")
    f3 <- presence ~ 1 + s(stem_prob, lcp_prob, k=6)
    r3 <- gam(f3, data=dat1_train, family="binomial")
    f4 <- presence ~ 1 + s(max_prob, k=3)
    r4 <- gam(f4, data=dat1_train, family="binomial")
  }
  if(paste(unique(dat1_train$tech_type), collapse="_") != "bnd_ran"){
    f0 <- presence ~ 1 + s(bird_id, bs="re")
    r0 <- gam(f0, data=dat1_train, family="binomial")
    f1 <- presence ~ 1 + s(stem_prob, k=3) + s(bird_id, bs="re")
    r1 <- gam(f1, data=dat1_train, family="binomial")
    f3 <- presence ~ 1 + s(stem_prob, lcp_prob, k=6) + s(bird_id, bs="re")
    r3 <- gam(f3, data=dat1_train, family="binomial")
    f4 <- presence ~ 1 + s(max_prob, k=3) + s(bird_id, bs="re")
    r4 <- gam(f4, data=dat1_train, family="binomial")
  }
  # check it out: summary(r1); summary(r3); summary(r4)
  
  # model fit summaries
  aic_tab <- AIC(r0, r1, r3, r4) %>% 
    mutate(Model=c("GAM_Null", "GAM_STEM", "GAM_STEM_LCP", "GAM_Max_STEM_LCP"),
            D_explained=c(summary(r0)$dev.expl,
            summary(r1)$dev.expl,
            summary(r3)$dev.expl,
            summary(r4)$dev.expl)) %>% 
    arrange(AIC) %>% mutate(Delta_AIC=c(0, diff(AIC))) %>%
    mutate_if(is.numeric, round, 2) %>%
    select(Model, D_explained, AIC, Delta_AIC)

  # partial plot for full model
  b <- getViz(r3)
  p_par <- plot(sm(b, 1), trans=plogis, n=100) + l_fitRaster() + 
    l_fitContour(colour=1) + 
    labs(title = NULL, x="STEM index", y="LCP index") + 
    scale_fill_viridis_c("Full\nmodel\npredicted\nindex", limits=c(0,1))
  # check it out: p_par
  
  # get model fitted values
  dat1_train$bird_id <- dat1_train$bird_id[1]
  dat1_train$r1_fit <- predict_gamm(model=r1, newdata=dat1_train, re_form = NA, 
                                type="response", exclude="s(bird_id)")[,1]
  dat1_train$r3_fit <- predict_gamm(model=r3, newdata=dat1_train, re_form = NA, 
                                type="response", exclude="s(bird_id)")[,1]
  dat1_train$r4_fit <- predict_gamm(model=r4, newdata=dat1_train, re_form = NA, 
                                type="response", exclude="s(bird_id)")[,1]
  
  # get auc 
  roc1 <- roc(presence~r1_fit, data=dat1_train)
  roc3 <- roc(presence~r3_fit, data=dat1_train)
  roc4 <- roc(presence~r4_fit, data=dat1_train)
  p_roc <- ggroc(list(roc1, roc3), legacy.axes = T) + 
      geom_abline(intercept=0, slope=1, lty=2, col="gray60") +
      labs(x="1 - Specificity", y="Sensitivity", color="Model") +
      scale_color_brewer(labels=c("STEM\nmodel", "STEM + LCP\nmodel"), palette="Set1") +
    labs(title=NULL, x="False positive rate", y="True positive rate") +
    annotate("text", x = 1, y = 0.1, hjust = 1, 
             label = paste("STEM model AUC =", 
                           str_pad(round(roc1$auc, 2), width=4, 
                                   side ="right", pad=0))) +
    annotate("text", x = 1, y = 0, hjust = 1,
             label = paste("STEM + LCP model AUC =", 
                           str_pad(round(roc3$auc, 2), width=4, 
                                   side ="right", pad=0)))
  # check it out: p_roc
  
  # add auc to fit table and save
  aic_tab$AUC <- NA
  aic_tab[aic_tab$Model=="GAM_Max_STEM_LCP", "AUC"] <- round(roc4$auc, 2)
  aic_tab[aic_tab$Model=="GAM_STEM_LCP", "AUC"] <- round(roc3$auc, 2)
  aic_tab[aic_tab$Model=="GAM_STEM", "AUC"] <- round(roc1$auc, 2)
  write.csv(aic_tab, paste0(out_path, "/", species_name, "/", species_name,
                    "_", focal_season, "_", "aic_table.csv"), row.names=F, na="")
  row.names(aic_tab) <- NULL
  # check it out: aic_tab
  
  # make predictor stack
  predictors <- stack(stem_map, lcp_map)
  predictors$max_prob <- max(predictors[[1:2]], na.rm=T)
  predictors$bird_id <- predictors$stem_prob
  values(predictors$bird_id) <- dat1_train$bird_id[1]
  # check it out: plot(predictors)
  
  # make predictions stack
  r1_pred <- predict(predictors, r1, exclude="s(bird_id)")
  r3_pred <- predict(predictors, r3, exclude="s(bird_id)")
  r4_pred <- predict(predictors, r4, exclude="s(bird_id)")
  pred_stack <- stack(r1_pred, r3_pred, r4_pred)
  names(pred_stack) <- c("r1_pred", "r3_pred", "r4_pred")
  for(s in 1:nlayers(pred_stack)){
    ra <- calc(pred_stack[[s]], plogis)
    rb <- rescale_ras_01(ra)
    pred_stack[[s]] <- rb
  }
  names(pred_stack) <- c("r1_pred", "r3_pred", "r4_pred")
  writeRaster(pred_stack$r3_pred, paste0(out_path, "/", species_name, "/", species_name,
                    "_", focal_season, "_", "integrated_surface.tif"), overwrite=T)
  writeRaster(pred_stack$r4_pred, paste0(out_path, "/", species_name, "/", species_name,
                    "_", focal_season, "_", "max_surface.tif"), overwrite=T)
  # check it out: plot(pred_stack)
  
  # make df for plots
  df_lcp <- gplot_data(lcp_map, maxpixels=ncell(lcp_map))
  df_stem <- gplot_data(stem_map, maxpixels=ncell(stem_map))
  df_r1 <- gplot_data(pred_stack$r1_pred, 
                       maxpixels=ncell(pred_stack$r1_pred))
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
    geom_sf(data=hemisphere_map_sf, fill=NA) +
    ylim(range(st_coordinates(study_region_map_sf)[,"Y"])) +
    xlim(range(st_coordinates(study_region_map_sf)[,"X"])) +
    theme_map()
  p_stem <- ggplot() +
    geom_tile(data=filter(df_stem, !is.na(value), value>=0.01),
              aes(x=x, y=y, fill=value)) +
   scale_fill_gradientn("Fall\nSTEM\nindex", 
                         colours = rev(terrain.colors(10)),
                         limits=c(0,1)) +
    guides(colour = "none") +
    geom_sf(data=hemisphere_map_sf, fill=NA) +
    ylim(range(st_coordinates(study_region_map_sf)[,"Y"])) +
    xlim(range(st_coordinates(study_region_map_sf)[,"X"])) +
    theme_map()
  p_r3 <- ggplot() +
    geom_tile(data=filter(df_r3, !is.na(value), value>=0.01), aes(x=x, y=y, 
                                                               fill=value)) +
    scale_fill_gradientn("Full\nmodel\npredicted\nindex", 
                         colours = rev(terrain.colors(10)),
                         limits=c(0,1)) +
    scale_color_brewer("Data\ntype", palette="Set1") +
    geom_sf(data=hemisphere_map_sf, fill=NA) +
    ylim(range(st_coordinates(study_region_map_sf)[,"Y"])) + 
    xlim(range(st_coordinates(study_region_map_sf)[,"X"])) +
    theme_map()
  p_r4 <- ggplot() +
    geom_tile(data=filter(df_r4, !is.na(value), value>=0.01), aes(x=x, y=y, 
                                                               fill=value)) +
    scale_fill_gradientn("Max\nmodel\npredicted\nindex", 
                         colours = rev(terrain.colors(10)),
                         limits=c(0,1)) +
    scale_color_brewer("Data\ntype", palette="Set1") +
    geom_sf(data=hemisphere_map_sf, fill=NA) +
    ylim(range(st_coordinates(study_region_map_sf)[,"Y"])) + 
    xlim(range(st_coordinates(study_region_map_sf)[,"X"])) +
    theme_map()
  p_pts <- ggplot() +
    geom_sf(data=hemisphere_map_sf, fill="gray60") +
    geom_sf(data=filter(dat1, presence==1) %>%
              dplyr::mutate(tech_type=str_replace(tech_type, "bnd", "Band")) %>%
              dplyr::mutate(tech_type=str_replace(tech_type, "ptt", "PTT")) %>%
              dplyr::mutate(tech_type=str_replace(tech_type, "gps", "GPS")) %>%
              dplyr::mutate(tech_type=str_replace(tech_type, "llg", "LLG")),
            size=0.5, shape=16, aes(col=factor(tech_type))) +
    scale_color_brewer("Fall\nmovement\ndata", palette="Set1") +
    ylim(range(st_coordinates(study_region_map_sf)[,"Y"])) + 
    xlim(range(st_coordinates(study_region_map_sf)[,"X"])) +
    theme_map()
  plot_grid(p_stem, p_lcp, p_pts, p_par$ggObj, p_roc, p_r3, nrow=2)
  ggsave2(filename = paste0(out_path, "/", species_name, "/", species_name,
                    "_", focal_season, "_", "prediction_maps.pdf"),  
          width=12, height=7)
  # ----------------------------------------------------------------------------

  
  # finish species
  print(paste("Done with", spp1))
  }, error=function(e){})
} # end spp loop
# ------------------------------------------------------------------------------



