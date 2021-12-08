###########################################################################################################
# Estimating migratory connectivity (proportions of indivs moving between breeding and wintering areas)
# First part of script pulls tracking and banding data for a given species, then creates m-arrays
# Second part of script uses these inputs in a model adapted from Korner-Nivergelt et al. 2017
# Purpose: Meehan et al. Integrating data types to estimate spatial patterns of avian migration
# across the Western Hemisphere.
# Outputs: data frame of mean proportions per breeding-wintering connection and associated error; 
# Fig illustrating those quantities
# Code developed by S. Saunders and L. Taylor in 2020 - 2021.
#########################################################################################################

library(data.table)
library(tibble)
library(dplyr)
library(tidyr)
library(sf)
library(mapview)
library(rgdal)
library(sp)
library(jagsUI)
library(ggplot2)

modefx <- function(x) {  x=x[!is.na(x)]; return(unique(x)[which.max(tabulate(match(x, unique(x))))]) } #use mode to sort individuals' locations into polys

setwd('Z:/Migratory_Bird_Initiative/') #set working directory

# projections
crs_sinu <- '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
crs_laea <- '+proj=laea +lat_0=15 +lon_0=-90 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs'
crs_wgs84 <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

# inputs
start_yr <- 1930

sp.list <- c('amwpel','bkpwar','brwhaw','graspa','greegr',
             'osprey','ovenbi1','prawar','prowar','swahaw','treswa',
             'turvul') #specify species list

for (sp in sp.list){ #loop through all species
  species_code <- sp 
  
  # create output folder for each species
  dir.create(path = paste0("Z:/Migratory_Bird_Initiative/Migratory Connectivity/", species_code, '_new'), recursive = TRUE)
  
  #load range-filtered BMRs
  clust_path <- paste0('LeastCostPaths/make_polys/output/', species_code, '_bird_migration_regions.shp')
  clust <- st_read(clust_path) %>% st_transform(crs_laea) %>%
    st_cast('POLYGON') %>% dplyr::select(season=season, region=mcr_class)
  
  # load banding and tracking data and filter outliers
  data_path <- paste0('Migration_Data/Species/', species_code, '/', species_code, '_clean_data.csv')
  points <- fread(data_path) %>% dplyr::filter(year >= start_yr) %>%  
    dplyr::mutate(season=case_when(season=='spring' ~ 'prebreeding_migration',
                                   season=='summer' ~ 'breeding',
                                   season=='fall' ~ 'postbreeding_migration',
                                   season=='winter' ~ 'nonbreeding',
                                   TRUE ~ '')) %>%
    dplyr::filter(!grepl("rng",outlier),
                  !grepl("ang",outlier),
                  !grepl("dst",outlier),
                  !grepl("mrk",outlier),
                  !grepl("dup",outlier),
                  !grepl("ddd",outlier)) %>%
    st_as_sf(coords=c('x_coord', 'y_coord'), remove=FALSE, crs=crs_wgs84) %>%
    st_transform(crs_laea)
  
  # season-specific inputs
  points_full <- NULL
  for (s in c('breeding', 'nonbreeding')) {
    
    # find nearest density cluster
    clust_s <- dplyr::filter(clust, season==s)
    buff_s <- st_buffer(clust_s, 250000) %>% group_by(season) %>% summarise() #250 km max distance
    points_s <- dplyr::filter(points, season==s) %>%
      st_join(clust_s, join=st_nearest_feature) %>% # sort to nearest feature 
      st_intersection(buff_s) %>% 
      st_set_geometry(NULL)
    
    # create capture history
    if (s=='breeding') {
      chistories <- dplyr::filter(points_s) %>% mutate(dummy=1) %>%    
        full_join(tibble(year=start_yr:year(Sys.Date())), by='year') %>%
        distinct(data_source, study_code, tech_type, species_code,
                 bird_id, region, year, dummy) %>%
        arrange(year) %>% spread(year, dummy, fill=0)
    }
    
    # combine seasonal data
    points_full <- rbind(points_full, points_s)
    rm(clust_s,points_s) # clean up objects
    
  }
  
  # count totals for regions and tech type
  clust_counts <- points_full %>%
    dplyr::select(tech_type,species_code,study_code,data_source,bird_id,season=season.x,x_coord,y_coord,grid_id,outlier,region) %>%
    mutate(tech_type=ifelse(tech_type=='bnd', 'bnd', 'trk')) %>%
    group_by(data_source, study_code, tech_type, species_code, bird_id, season) %>%
    summarize(region=modefx(region)) %>% spread(season, region) %>%
    dplyr::filter(!is.na(breeding) & !is.na(nonbreeding)) %>%
    group_by(tech_type, breeding, nonbreeding) %>%
    count() %>% ungroup()
  
  # regions to include
  clust_b <- dplyr::filter(clust, season=='breeding')
  br_regions <- unique(sort(clust_b$region))
  clust_nb <- dplyr::filter(clust, season=='nonbreeding')
  nb_regions <- tibble(id=unique(sort(clust_nb$region)),
                       breeding=0, dummy=0) %>% spread(id, dummy)
  
  # banding matrix
  m_banding <- dplyr::filter(clust_counts, tech_type=='bnd') %>% dplyr::select(-tech_type) %>%
    spread(nonbreeding, n, fill=0) %>% full_join(nb_regions) %>% 
    right_join(tibble(breeding=br_regions), by='breeding') %>%
    mutate_all(~replace_na(., 0)) %>% arrange(breeding) %>% column_to_rownames(var='breeding')
  
  # tracking matrix
  m_tracking <- dplyr::filter(clust_counts, tech_type=='trk') %>% dplyr::select(-tech_type) %>%
    spread(nonbreeding, n, fill=0) %>% full_join(nb_regions) %>% 
    right_join(tibble(breeding=br_regions), by='breeding') %>%
    mutate_all(~replace_na(., 0)) %>% arrange(breeding) %>% column_to_rownames(var='breeding')
  
  # make sure order of columns matches between m_banding and m_tracking
  m_tracking<-m_tracking[,names(m_banding)]
  
  # save sample sizes for reference
  bnd_ind <- sum(m_banding)
  trk_ind <- sum(m_tracking)
  
  ssizes <- c(bnd_ind,trk_ind)
  names(ssizes) <- c('BndResight','Tracker')
  ssizes <- as.data.frame(ssizes)
  write.csv(ssizes, file=paste0("Z:/Migratory_Bird_Initiative/Migratory Connectivity/", species_code, '_new/',species_code,"_MigConnect_SampleSizes.csv"))
  
  # function to create a m-array based on capture-recapture data (CH)
  marray <- function(CH){
    nind <- dim(CH)[1]
    n.occasions <- dim(CH)[2]
    m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
    
    # Calculate the number of released individuals at each time period
    for (t in 1:n.occasions){
      m.array[t,1] <- sum(CH[,t])
    }
    for (i in 1:nind){
      pos <- which(CH[i,]!=0)
      g <- length(pos)
      for (z in 1:(g-1)){
        m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
      } #z
    } #i
    
    # Calculate the number of individuals that is never recaptured
    for (t in 1:n.occasions){
      m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
    }
    out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
    return(out)
  }
  
  # function to subset capture histories by region, then create marray
  func <- function(reg){
    chistories.reg <- subset(chistories,region==reg)
    CH.Reg <- chistories.reg[,as.character(start_yr:year(Sys.Date()))]
    CH.AReg <- data.matrix(CH.Reg)
    CH.A.marray <- marray(CH.AReg)
    return(CH.A.marray)
  }
  
  # supplement capture histories with dummy histories for missing breeding regions
  br_regions_ch <- unique(sort(chistories$region))
  toadd <- setdiff(br_regions, br_regions_ch)
  
  if (length(toadd) > 0) {
    dumhist <- c(rep(NA,6),rep(0,91),1) 
    addmat <- matrix(dumhist,nrow=length(toadd),ncol=98, byrow = TRUE) 
    addmat[,6] <- toadd
    colnames(addmat) <- colnames(chistories)
    # add back to chistories
    chistories <- rbind(addmat,chistories)
  } else if (length(toadd==0)){
    chistories <- chistories
  }
  
  # convert to array where third dimension is number of unique regions
  n_yrs <- year(Sys.Date()) - start_yr
  multi.marray <- array(data=NA, dim=c(n_yrs,n_yrs+1,length(br_regions)))
  for (i in 1:length(br_regions)){
    multi.marray[,,i] <- func(reg=br_regions[i])
  } 
  
  # Specify data for model --------------------------------------------------------------------------------------
  
  #Number of cols in m-arrays (1930 - 2021)
  J <- ncol(multi.marray)
  
  #Total number of 'marked and released' individuals (sum of rows in each m-array)
  n <- colSums(multi.marray[,ncol(multi.marray),]) 
  
  #Number of populations (number of m-arrays)
  npop <- nrow(m_banding) #breeding pops
  nonpop <- ncol(m_banding) #non breeding pops
  
  # rename non-breeding reencounter data
  wirec <- m_banding
  
  #rename geolocator data
  geolocs <- m_tracking
  
  #Total number of indivs with geoloc data (row sums of above)
  ngeolocs <- rowSums(geolocs)
  
  ## Model ##########################
  
  sink("migconnectivity.gen")
  cat("
model{
  # priors  for models to estimate number of ringed  (see Korner-Nievergelt et al. 2012)
  phi~dunif(0,1)           # yearly survival rate
  for(j in 1:npop){
    p[j]~dunif(0,1)        # recapture probability per breeding area
  
  # multinomial likelihood for model to estimate numbers of ringed birds
  for(t in 1:(J-1)){
    multi.marray[t,1:J,j]~dmulti(pr[t,,j], rel[t,j])
}
  # cell probabilities
    q[j] <- 1-p[j]
  
  for(t in 1:(J-1)){      # main diagonal
    pr[t,t,j] <- phi*p[j]
    for(k in (t+1):(J-1)){        # above main diagonal
      pr[t,k,j] <- pow(phi, (k-t+1))*pow(q[j], (k-t))*p[j]
    }
}
  for(t in 2:(J-1)){      # below main diagonal
    for(k in 1:(t-1)){
      pr[t,k,j] <- 0
    }
  }
  # last column: probability of non recapture
  for(t in 1:(J-1)){
    pr[t,J,j] <- 1-sum(pr[t,1:(J-1),j])
  }
  
  # number of ringed in the population
    N[j] <- n[j]/(1-(1-phi)/(1-phi*(1-p[j])))  
  
  # migration rate model based on ringing data
    for(a in 1:nonpop){ #number of non-breeding pops
      r.korr[j,a]<-r[a]*m[j,a]
      muW[j,a] <- r.korr[j,a] * N[j]
      wirec[j,a]~dpois(muW[j,a])
    }
  
  # migration rate model based on tracking data
    geolocs[j,1:nonpop] ~ dmulti(m[j,], ngeolocs[j]) 
  }
  
  # priors for migration rate model   
  for(a in 1:nonpop){  # based on nonbrdg pops
    r[a]~dunif(0,1)         # re-encounter probabilities
  }
  
  for(j in 1:npop){
    for(a in 1:nonpop){
      m0[j,a]~dunif(0,1)
      m[j,a] <- m0[j,a]/sum(m0[j,])
    }
  }
}
    ",fill = TRUE)
  sink()
  
  # Bundle data
  #make release matrix of number of years by number of breeding regions
  rel <- c()
  for (j in 1:npop){
    dat <- rowSums(multi.marray[,,j])
    rel <- cbind(rel,dat)
  }
  
  jags.data <- list(multi.marray = multi.marray, J = J, n = n, 
                    rel = rel, wirec = wirec, npop = npop, nonpop = nonpop, geolocs = geolocs, ngeolocs = ngeolocs) 
  
  # Initial values
  inits <- function(){list(phi = runif(1, 0.1, 0.9), p = runif(npop, 0.1, 0.9), r = runif(nonpop, 0.1, 0.9))}   
  
  # Parameters monitored
  parameters <- c("phi", "p", "m", "r", "N", "muW") 
  
  # MCMC settings
  nt <- 3
  nb <- 50000
  nc <- 3
  
  # Call JAGS from R
  #use autojags within species loop to run each model to convergence
  migconnect <- autojags(jags.data, inits, parameters, "migconnectivity.gen", n.chains = nc, n.thin = nt, n.burnin = nb, 
                         parallel = TRUE, Rhat.limit=1.1, iter.increment=30000, max.iter=200000, verbose = FALSE) 
  print(migconnect, digits=3)
  
  #---------------------------------------------------------
  #figures and summary stats
  #---------------------------------------------------------
  nbrd <- npop
  nnon <- nonpop
  
  quants.low <- matrix(NA,nrow=nbrd,ncol=nnon)
  quants.low85 <- matrix(NA,nrow=nbrd,ncol=nnon)
  quants.high <- matrix(NA,nrow=nbrd,ncol=nnon)
  quants.high85 <- matrix(NA,nrow=nbrd,ncol=nnon)
  medians <- matrix(NA,nrow=nbrd,ncol=nnon)
  for (i in 1:nbrd){
    for (j in 1:nnon){
      quants.low[i,j] <- quantile(migconnect$sims.list$m[,i,j], probs=0.025)
      quants.low85[i,j] <- quantile(migconnect$sims.list$m[,i,j], probs=0.075)
      quants.high[i,j] <- quantile(migconnect$sims.list$m[,i,j], probs=0.975)
      quants.high85[i,j] <- quantile(migconnect$sims.list$m[,i,j], probs=0.925)
      medians[i,j] <- quantile(migconnect$sims.list$m[,i,j], probs=0.5)
    }
  }
  
  means <- c(migconnect$mean$m)
  
  #create data frame to save summary results
  df <- data.frame(
    mean = means,
    median = c(medians),
    low95 = c(quants.low),
    up95 = c(quants.high),
    low85 = c(quants.low85),
    up85 = c(quants.high85),
    brdg = rep(c(rownames(m_tracking)),times=nnon),
    nonbrdg = rep(c(colnames(m_tracking)),each=nbrd)
  )
  
  #save data frame
  write.csv(df, file=paste0("Z:/Migratory_Bird_Initiative/Migratory Connectivity/", species_code, '_new/',species_code,"_MigConnect_RANGEfiltered_Results.csv"),row.names = FALSE)
  
  #violin plot with full posterior ----------------------------------------------------------------
  
  fullpost <- c()
  for (i in 1:nbrd){
    for (j in 1:nnon){
      dat = migconnect$sims.list$m[,i,j]
      fullpost = cbind(fullpost,dat)
    }
  }
  
  colnames(fullpost) <- seq(1,ncol(fullpost),by=1)
  fullpost <- as.data.frame(fullpost)
  full.stack <- fullpost %>%
    gather() #stack output
  full.means <- full.stack[,2] #save means
  
  dframe <- data.frame(full.means,breed=c(rep(rownames(m_tracking),each=nnon*migconnect$mcmc.info$n.samples)),
                       nonbreed=c(rep(colnames(m_tracking),each=migconnect$mcmc.info$n.samples,times=nbrd)))
  
  dodge <- position_dodge(width = 0.7)
  p <- ggplot(dframe, aes(x=nonbreed, y=full.means, color=breed)) + 
    geom_violin(position=dodge, lwd=0.75) +
    stat_summary(fun=base::mean, geom="point", size=3, position=dodge)
  p1 <- p +
    theme_bw()+
    xlab("Wintering Area")+
    ylab("Proportion")+
    theme(axis.title=element_text(size=14),
          axis.text=element_text(size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black"),
          panel.background = element_rect(fill = 'gray97', colour = 'gray97'))+
    theme(legend.position = "none")
  
  plot_final <- p1 + facet_wrap(.~breed,ncol=round(npop/4,0)+1)
  
  #save output (fig)
  png(filename=paste0("Z:/Migratory_Bird_Initiative/Migratory Connectivity/", species_code, '_new/',species_code,"_MigConnect_RANGEfiltered_Props.png"), width=2.5*(round(npop/4,0)+1), height=10.5, units="in", res=300)
  print(plot_final)
  dev.off()
  
} #end species loop
