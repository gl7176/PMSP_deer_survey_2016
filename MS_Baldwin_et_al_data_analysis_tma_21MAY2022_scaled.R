# Comparison of Analysis Methods for White-tailed Deer Population Surveys
# Baldwin et al. 2019
# Prepared November 2019

# Reanalysis of PMSP deer abundance data
# 07 Jan 2022

# updated 15 May 2022 to model by day using FUN = max and an area
# of pi*(318)^2 * 22 sites = 6.99 km2

# name of weather station:
# PILOT MOUNTAIN, NC. CBTN7 (RAWS)

# Note from Jared 4/9/2022:
# The reason we had 22 sites is because Jacobson et al. showed using marked 
# animals that a baited camera site for every 100 acres would capture ~90% of 
# the deer in the area. We used a 100 acre grid system for camera placement, so 
# distance between camera sites would be 636.15 meters or a radius of 318 meters
# per camera site.

# clear workspace
  rm(list=ls())

# options
  options("digits" = 8, "scipen" = 12)
  
# setwd
  setwd("E:/ANDERSON_BACKUP_26MAR2022/_WORKING_PROJECTS/BALDWIN_PMSP_DEER_PROJECT/REANALYSIS_JAN2022")

# load libraries
  library(car)
  library(unmarked)
  library(AICcmodavg)
  library(reshape)
  library(circular)

# ------------------ TOTAL SURVEY AREA ----------------------------------#    
  
# 10.24 comes from total area/22 sites  
# PMSP is 3872 acres or 15.67 km2
  
# calculate camera area per trap:    
  cam.area.m <- pi*(636.15/2)^2 # in m
  cam.area.km <- pi*(0.63615/2)^2 # in km
  
  cam.area.m/(1000*1000) * 22
  survey.area <- cam.area.km * 22 # 6.9924863 km2  

# ----------------- FLIR VALIDATION DENSITY ESTIMATION ----------------- #

# read-in observer data
  flir.dat <- read.csv("flir_obs.csv")
  transect.dat <- read.csv("transect_information.csv")
  
# set seed for reproducibility (applies for all future analyses as well)
  set.seed(4)

# merge dataframes for anaysis
  flir.dat <- merge(flir.dat, transect.dat, by.x = "TRANSECT", by.y = "TRANSECT")

# convert counts to densities
  flir.dat$DENSITY <- (flir.dat$COUNT / flir.dat$AREA)

# subset data by flight
  f.1 <- subset(flir.dat, FLIGHT == 1, drop = T)
  f.2 <- subset(flir.dat, FLIGHT == 2, drop = T)
  f.3 <- subset(flir.dat, FLIGHT == 3, drop = T)
  f.4 <- subset(flir.dat, FLIGHT == 4, drop = T)
  f.5 <- subset(flir.dat, FLIGHT == 5, drop = T)

# aggregate data by transect and observer
  f1trans <- aggregate(f.1$DENSITY, by = list(f.1$TRANSECT, f.1$OBSERVER), FUN = "sum")
  names(f1trans) <- c("TRANSECT","OBSERVER", "DENSITY")
  f1trans <- aggregate(f1trans$DENSITY, by = list(f1trans$TRANSECT), FUN = 'mean')
  names(f1trans) <- c("TRANSECT","DENSITY")

  f2trans <- aggregate(f.2$DENSITY, by = list(f.2$TRANSECT, f.2$OBSERVER), FUN = "sum")
  names(f2trans) <- c("TRANSECT","OBSERVER", "DENSITY")
  f2trans <- aggregate(f2trans$DENSITY, by = list(f2trans$TRANSECT), FUN = 'mean')
  names(f2trans) <- c("TRANSECT","DENSITY")
  
  f3trans <- aggregate(f.3$DENSITY, by = list(f.3$TRANSECT, f.3$OBSERVER), FUN = "sum")
  names(f3trans) <- c("TRANSECT","OBSERVER", "DENSITY")
  f3trans <- aggregate(f3trans$DENSITY, by = list(f3trans$TRANSECT), FUN = 'mean')
  names(f3trans) <- c("TRANSECT","DENSITY")
  
  f4trans <- aggregate(f.4$DENSITY, by = list(f.4$TRANSECT, f.4$OBSERVER), FUN = "sum")
  names(f4trans) <- c("TRANSECT","OBSERVER", "DENSITY")
  f4trans <- aggregate(f4trans$DENSITY, by = list(f4trans$TRANSECT), FUN = 'mean')
  names(f4trans) <- c("TRANSECT","DENSITY")

  f5trans <- aggregate(f.5$DENSITY, by = list(f.5$TRANSECT, f.5$OBSERVER), FUN = "sum")
  names(f5trans) <- c("TRANSECT","OBSERVER", "DENSITY")
  f5trans <- aggregate(f5trans$DENSITY, by = list(f5trans$TRANSECT), FUN = 'mean')
  names(f5trans) <- c("TRANSECT","DENSITY")

# join data to calculate average density between flights
  density.dat <- rbind(f1trans, f2trans, f3trans, f4trans, f5trans)
  mean(density.dat$DENSITY)
  density <- c(mean(f1trans$DENSITY), mean(f2trans$DENSITY), mean(f3trans$DENSITY), mean(f4trans$DENSITY), mean(f5trans$DENSITY))

# create dataframe for modelling
  fl.agg <- data.frame(Flight = c(rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 6), rep(5, 6)), 
                       Transect = c(f1trans$TRANSECT, f2trans$TRANSECT, f3trans$TRANSECT, f4trans$DENSITY, f5trans$DENSITY), 
                       Density = c(f1trans$DENSITY, f2trans$DENSITY, f3trans$DENSITY, f4trans$DENSITY, f5trans$DENSITY))
  
# bootstrapped CI for density
  bstrap <- NULL
  for (i in 1:10000){
  bstrap <- c(bstrap, mean(sample(density, 5, replace=T)))}
  fl.ci <- rbind(c(quantile(bstrap,.025), quantile(bstrap, .975)))

  fl.mn <- mean(density)
  fl.sd <- sd(density)
  fl.df <- data.frame(TYPE='UAI', MEAN=fl.mn, SD=fl.sd, LOWER=fl.ci[1], UPPER=fl.ci[2])
  row.names(fl.df) <- 1:length(fl.df$MEAN)  
  
# run ANOVA on flight estimates
  m1 <- lm(fl.agg$Density ~ fl.agg$Flight)
  Anova(m1, typ = '2')

# ----------------- MR BOOTSTRAPPING FOR PRECISON ESTIMATION ----------------- #

# remove old objects
  #rm(list = ls())
  
# read in data (UM stands for number of unique males at each site)
  mr.dat <- read.csv("mr_error.csv")
  head(mr.dat)
  str(mr.dat)
  
# simulate survey 10,000 times 
  set.seed(4)
  
  N <- NULL
  sp <- NULL
  for(i in 1:10000){
    sp <- data.frame(SITE = sample(c(1:22), 22, replace = T))
    mg <- merge(sp, mr.dat, by.y = 'SITE')
    correction <- sum(mg$MALE) / sum(mg$UM)
    FC <- sum(mg$FEMALE)/correction
    YC <- sum(mg$FAWN)/correction
    N[i] <- sum(c(FC, YC, sum(mg$UM)))
  }
  
  mr.l <- quantile(N, 0.025) / 10.24
  mr.u <- quantile(N, 0.975) / 10.24
  mr.mn <- mean(N) / 10.24
  # mr.l <- quantile(N, 0.025) / survey.area
  # mr.u <- quantile(N, 0.975) / survey.area
  # mr.mn <- mean(N) / survey.area
  mr.sd <- sd(N)
  
  mr.df <- data.frame(TYPE='MR', MEAN=mr.mn, SD=mr.sd, LOWER=mr.l, UPPER=mr.u)
  row.names(mr.df) <- 1:length(mr.df$MEAN)

# ----------------- N-MIXTURE MODELLING DEER ABUNDANCE SUMMER 2016 - SPRING 2018 ----------------- #

# remove old objects
  #rm(list=ls())

# load camera covariate data
  cam.dat <- read.csv("camera_covariates.csv")  

# convert aspect data to something more meaningful:  
  abs.rad <- abs(cam.dat$aspect * pi/180) # gives absolute value of radians
  cos(abs.rad) # "northness", a value from 1 to -1, where positive values are north-facing
  
  # NOTE: in the Northern Hemisphere, north-facing slopes in latitudes from about 30 to 55 degrees 
  # receive less direct sunlight than south-facing slopes. The lack of direct sunlight throughout 
  # the day, whether in winter or summer, results in north-facing slopes being cooler than south-facing 
  # slopes. During winter months, portions of north-facing slopes may remain shaded throughout the day 
  # due to the low angle of the sun. This causes snow on north-facing slopes to melt slower than on 
  # south-facing ones.
  
  #sin(x) # "eastness", a value also from 1 to -1, where positive values are east-facing 
  
# create site covariate dataframe and scale data; ROAD is not used
  deer.covars <- data.frame(ELE=cam.dat$site.ele, ASPECT = cos(abs(cam.dat$aspect*pi/180)), SLOPE = cam.dat$slope, EDGE = cam.dat$EDGE, ROAD = cam.dat$ROAD)
  deer.covars.scaled <- scale(deer.covars)
  deer.covars.scaled <- data.frame(deer.covars.scaled[1:22, ])
  
# look at scaed variables
  str(deer.covars.scaled)
  head(deer.covars.scaled)

# load camera trap data
  deer.dat <- read.csv("deer_obs_full_period.csv")

# create list of ordered levels for sorting
  out <- NULL
  for(i in 1:22) {
  tmp <- noquote(paste("PM",i,sep=""))
  out <- append(out, tmp)
  } 
  deer.dat$SITE <- ordered(deer.dat$SITE, levels=out)

# convert DATE.TIME to POSIX
  deer.dat$DATE.TIME <- as.POSIXct(strptime(deer.dat$DATE.TIME, format = "%m/%d/%Y %H:%M"))

# order observations by site and date.time
  deer.dat <- deer.dat[order(deer.dat[,"SITE"], deer.dat[,"DATE.TIME"]),]

### HERE IS WHERE WE AGGREGATE BY HOUR - THIS MAKES SENSE  
  
# aggregate by max count by hour prior to merge
  deer.dat$DATE.TIME.HOUR <- round.POSIXt(deer.dat$DATE.TIME, "hours")
  count.by.hour <- aggregate(deer.dat$COUNT, by = list(deer.dat$SITE, as.character(deer.dat$DATE.TIME.HOUR)), FUN = max)
  names(count.by.hour) <- c("SITE", "DATE.HOUR", "COUNT")
  count.by.hour <- count.by.hour[order(count.by.hour$SITE, count.by.hour$DATE.HOUR),]

# make continuous calender of all observation periods for merge
  all.sites <- unique(deer.dat$SITE)

# convert date to POSIXct
  deer.dat$DATE.TIME <- as.POSIXct(deer.dat$DATE.TIME, format = '%m/%d/%Y %H:%M')

# order observations by site and time
  deer.dat <- deer.dat[order(deer.dat$SITE, deer.dat$DATE.TIME),]

# create a sequence of date-times from the start to the finish of the survey
  min.date.time <- as.POSIXct(strptime("2016-06-09 01:00:00", "%Y-%m-%d %H:%M:%S"))
  max.date.time <- as.POSIXct(strptime("2018-12-31 23:00:00", "%Y-%m-%d %H:%M:%S"))  
  
# create backbone of all possible hour, date and site combinations
  seq.date.time <- seq.POSIXt(from = min.date.time, to = max.date.time, by = "hour")
  rep.date.time <- rep(seq.date.time, length(all.sites))
  length(rep.date.time) # 494208
  
# create site reps for new backbone data frame
  rep.sites <- sort(rep(all.sites, length(seq.date.time)))
  length(rep.sites) # 494208
  
# create backbone dataframe for merging with data 
  survey.sites.df <- data.frame(SITE = rep.sites, DATE.HOUR = rep.date.time)
  head(survey.sites.df)
  
# look at count.by.hour  
  min(count.by.hour$DATE.HOUR) # "2015-09-12 17:00:00"
  max(count.by.hour$DATE.HOUR) # "2019-01-29 12:00:00"
  
# merge site by survey backbone dataframe with aggregated count data
  deer.df <- merge(x = survey.sites.df, y = count.by.hour, by.x = c("SITE","DATE.HOUR"), by.y = c("SITE","DATE.HOUR"), all.x = TRUE, all.y = FALSE)
  deer.df[is.na(deer.df)] <- 0
  head(deer.df, 100)
  tail(deer.df, 100)
  summary(deer.df)
  min(deer.df$DATE.HOUR) # "2016-06-09 01:00:00 EDT"
  max(deer.df$DATE.HOUR) # "2019-12-31 23:00:00 EST"
  
# prepare backbone for aggregation by hour  
  (days.on <- as.numeric(round(max.date.time - min.date.time, 0))) # 936
  temp <- deer.df[order(deer.df$DATE.HOUR, deer.df$SITE),]
  rownames(temp) <- 1:nrow(temp)
  temp$DATE <- as.Date(temp$DATE.HOUR, tz = '')
  length(temp$DATE) # 494208
  head(temp)
  tail(temp)
  #deer.df <- temp

### INSPECT THIS STEP - WHY AGGREGATE BY DAY???  
### on 17 MAY 2022 I changed this step and reran the analyses on the sum of the one hour data 
### models failed to converge

### on 18 MAY 2022 I changed this step and reran the analyses on the max of a six-hour window 
  # 6-hour: estimate too high
  # 12-hour: estimate too high 
  # one-day:

# # 6-hour
#   temp.6 <- temp[order(temp$SITE, temp$DATE.HOUR), ] # sort by site and time
#   temp.6.days.on <- as.numeric(max(temp.6$DATE)-min(temp.6$DATE))+1
#   temp.6.all <- temp.6.days.on*22*4
#   tmp6h <- sort(rep(1:temp.6.all, 6)); length(tmp6h) # 22 sites * 935 days * 4 x 6-hour periods, repeated for 6 hours = 494208
#   temp.6$AGG <- tmp6h[1:length(temp.6$SITE)]
#   temp.6agg <- aggregate(COUNT ~ SITE + AGG, data = temp.6, FUN = "max")
#   temp.6agg <- temp.6agg[order(temp.6agg$SITE, temp.6agg$AGG), ]
#   row.names(temp.6agg) <- 1:length(temp.6agg$SITE); length(temp.6agg$SITE)
#   temp.6.dates <- rep(sort(rep(unique(temp.6$DATE), 4)), 22); length(temp.6.dates)
#   temp.6.sites <- sort(rep(unique(temp.6$SITE), 4*temp.6.days.on)); length(temp.6.sites)
#   temp.6.sites.dates <- data.frame(SITE=temp.6.sites, DATE=temp.6.dates)
#   temp.6agg <- data.frame(temp.6.sites.dates, temp.6agg)
#   temp.6agg <- temp.6agg[, c('SITE', 'DATE', 'AGG', 'COUNT')]
#   temp.6agg$SIX.HOUR <- rep(1:(temp.6.days.on*4), 22)
#   deer.df <- temp.6agg

# # 12-hour
#   temp.12 <- temp[order(temp$SITE, temp$DATE.HOUR), ] # sort by site and time
#   temp.12.days.on <- as.numeric(max(temp.12$DATE)-min(temp.12$DATE))+1
#   temp.12.all <- temp.12.days.on*22*2
#   tmp12h <- sort(rep(1:temp.12.all, 12)); length(tmp12h) # 22 sites * 935 days * 4 x 6-hour periods, repeated for 6 hours = 494208
#   temp.12$AGG <- tmp12h[1:length(temp.12$SITE)]
#   temp.12agg <- aggregate(COUNT ~ SITE + AGG, data = temp.12, FUN = "max")
#   temp.12agg <- temp.12agg[order(temp.12agg$SITE, temp.12agg$AGG), ]
#   row.names(temp.12agg) <- 1:length(temp.12agg$SITE); length(temp.12agg$SITE)
#   temp.12.dates <- rep(sort(rep(unique(temp.12$DATE), 2)), 22); length(temp.12.dates)
#   temp.12.sites <- sort(rep(unique(temp.12$SITE), 2*temp.12.days.on)); length(temp.12.sites)
#   temp.12.sites.dates <- data.frame(SITE=temp.12.sites, DATE=temp.12.dates)
#   temp.12agg <- data.frame(temp.12.sites.dates, temp.12agg)
#   temp.12agg <- temp.12agg[, c('SITE', 'DATE', 'AGG', 'COUNT')]
#   temp.12agg$TWELVE.HOUR <- rep(1:(temp.12.days.on*2), 22)
#   deer.df <- temp.12agg

# one-day  
  temp.agg <- aggregate(COUNT ~ SITE + DATE, data = temp, FUN = "max")
  temp.agg <- temp.agg[order(temp.agg$SITE, temp.agg$DATE),]
  deer.df <- temp.agg
  
# subset count dataframe by month for change in density analysis
  JUN16 <- subset(deer.df, format.Date(DATE, "%m")=="06" & format.Date(DATE, "%Y")=="2016")
  JUL16 <- subset(deer.df, format.Date(DATE, "%m")=="07" & format.Date(DATE, "%Y")=="2016")
  AUG16 <- subset(deer.df, format.Date(DATE, "%m")=="08"  & format.Date(DATE, "%Y")=="2016")
  OCT16 <- subset(deer.df, format.Date(DATE, "%m")=="10"& format.Date(DATE, "%Y")=="2016" )
  SEP16 <- subset(deer.df, format.Date(DATE, "%m")=="09" & format.Date(DATE, "%Y")=="2016") 
  NOV16 <- subset(deer.df, format.Date(DATE, "%m")=="11" & format.Date(DATE, "%Y")=="2016")
  DEC16 <- subset(deer.df, format.Date(DATE, "%m")=="12" & format.Date(DATE, "%Y")=="2016")
  JAN17 <- subset(deer.df, format.Date(DATE, "%m")=="01" & format.Date(DATE, "%Y")=="2017")
  FEB17 <- subset(deer.df, format.Date(DATE, "%m")=="02" & format.Date(DATE, "%Y")=="2017")
  MAR17 <- subset(deer.df, format.Date(DATE, "%m")=="03" & format.Date(DATE, "%Y")=="2017")
  APR17 <- subset(deer.df, format.Date(DATE, "%m")=="04" & format.Date(DATE, "%Y")=="2017")
  MAY17 <- subset(deer.df, format.Date(DATE, "%m")=="05" & format.Date(DATE, "%Y")=="2017")
  JUN17 <- subset(deer.df, format.Date(DATE, "%m")=="06" & format.Date(DATE, "%Y")=="2017")
  JUL17 <- subset(deer.df, format.Date(DATE, "%m")=="07" & format.Date(DATE, "%Y")=="2017")
  AUG17 <- subset(deer.df, format.Date(DATE, "%m")=="08" & format.Date(DATE, "%Y")=="2017")
  SEP17 <- subset(deer.df, format.Date(DATE, "%m")=="09" & format.Date(DATE, "%Y")=="2017")
  OCT17 <- subset(deer.df, format.Date(DATE, "%m")=="10" & format.Date(DATE, "%Y")=="2017")
  NOV17 <- subset(deer.df, format.Date(DATE, "%m")=="11" & format.Date(DATE, "%Y")=="2017")
  DEC17 <- subset(deer.df, format.Date(DATE, "%m")=="12" & format.Date(DATE, "%Y")=="2017") 
  JAN18 <- subset(deer.df, format.Date(DATE, "%m")=="01" & format.Date(DATE, "%Y")=="2018")
  FEB18 <- subset(deer.df, format.Date(DATE, "%m")=="02" & format.Date(DATE, "%Y")=="2018")
  MAR18 <- subset(deer.df, format.Date(DATE, "%m")=="03" & format.Date(DATE, "%Y")=="2018")
  APR18 <- subset(deer.df, format.Date(DATE, "%m")=="04" & format.Date(DATE, "%Y")=="2018")
  MAY18 <- subset(deer.df, format.Date(DATE, "%m")=="05" & format.Date(DATE, "%Y")=="2018")
  JUN18 <- subset(deer.df, format.Date(DATE, "%m")=="06" & format.Date(DATE, "%Y")=="2018")
  JUL18 <- subset(deer.df, format.Date(DATE, "%m")=="07" & format.Date(DATE, "%Y")=="2018")
  AUG18 <- subset(deer.df, format.Date(DATE, "%m")=="08" & format.Date(DATE, "%Y")=="2018")

# rbind into seasons
  summer16 <- rbind(JUN16, JUL16, AUG16)
  fall16 <- rbind(SEP16, OCT16, NOV16)
  winter16 <- rbind(DEC16, JAN17, FEB17)
  spring17 <- rbind(MAR17, APR17, MAY17)
  summer17 <- rbind(JUN17,JUL17, AUG17)
  fall17 <- rbind(SEP17, OCT17, NOV17)
  winter17 <- rbind(DEC17, JAN18, FEB18)
  spring18 <- rbind(MAR18, APR18, MAY18)
  summer18 <- rbind(JUN18,JUL18, AUG18)

# change months in winter17 to test sensitivity of timing, model, to the UAI estimates  
  #winter17 <- rbind(NOV17, DEC18, JAN18) 
  #winter17 <- rbind(JAN18, FEB18, MAR18)

# create a list of seasons to run functions on simultaneously
  season.list <- list(summer16, fall16, winter16, spring17, summer17, fall17, winter17, spring18, summer18)
  
# create simplified df of counts
  # season.list <- lapply(season.list, function(x) x <- x[, c("SITE","DATE", "SIX.HOUR", "COUNT")])
  # season.list <- lapply(season.list, function(x) x <- x[, c("SITE","DATE", "TWELVE.HOUR", "COUNT")])
  season.list <- lapply(season.list, function(x) x <- x[, c("SITE","DATE", "COUNT")])
  
# cast for unmarked  
  # season.list <- lapply(season.list, function (x) x <- cast(data = x, formula = SITE~SIX.HOUR, mean))      #timemeans <- cast(mdata, time~variable, mean) 
  # season.list <- lapply(season.list, function (x) x <- cast(data = x, formula = SITE~TWELVE.HOUR, mean))   #timemeans <- cast(mdata, time~variable, mean) 
  season.list <- lapply(season.list, function (x) x <- cast(data = x, formula = SITE~DATE, mean))          #timemeans <- cast(mdata, time~variable, mean) 
  season.list <- lapply(season.list, function (x) x <- as.matrix.cast_df(x[,-1])) 

# create unmarked frames with data 
  season.list.scaled <- lapply(season.list, function(x) unmarkedFramePCount(y = x, siteCovs = deer.covars.scaled))
  
# model of Royle 2004 based on counts
# general formula is: model <- pcount(~detection_formula ~occupancy_formula, dataframe, K=100, se=TRUE)

# run source command to execute pcount models in a separate script
  
# *** this will take a very long time ***  
#  source(file = "MS_Baldwin_et_al_run_models_21MAY2022.R", local = FALSE, echo = TRUE, )

# now pull out and combine all the results by season
  # setwd('./model_results_scaled_6hr')
  # setwd('./model_results_scaled_12hr')
  setwd('./model_results_scaled_1day')
  file.list <- list.files()
  
# import all results files
  for(i in 1:length(file.list)){
    f <- read.csv(file.list[i])
    n <- regmatches(file.list[i], regexpr( "\\d+", file.list[i]))
    name <- paste("mod", n, sep = '')
    assign(name, f)
  }
  
# subset each mod file into the nine seasons
  season_1 <- NULL
  season_2 <- NULL
  season_3 <- NULL
  season_4 <- NULL
  season_5 <- NULL
  season_6 <- NULL
  season_7 <- NULL
  season_8 <- NULL
  season_9 <- NULL

  for(i in 1:length(file.list)){
    f <- read.csv(file.list[i])
    n <- regmatches(file.list[i], regexpr( "\\d+", file.list[i]))
    m <- paste("mod", n, sep = '')
    f$MODEL <- rep(m, length(f$X))
    
      for(j in 1:9){
        seas <- paste("S_", j, sep = '')
        tmp <- data.frame(f[f$SEASON== seas, ])
        out <- rbind(get(paste("season_", j, sep = "")), tmp)
        assign(x = paste("season_", j, sep = ""), value = out)
      }
    
  }
  
  summer.2016 <- season_1
  fall.2016 <- season_2
  winter.2016 <- season_3
  spring.2017 <- season_4
  summer.2017 <- season_5 
  fall.2017 <- season_6
  winter.2017 <- season_7
  spring.2018 <- season_8
  summer.2018 <- season_9
  
  options(digits = 3, scipen = 3)
  
  summer.2016.df <- data.frame(summer.2016[summer.2016$AIC == min(unique(summer.2016$AIC)), ])  # mod 28
  fall.2016.df <- data.frame(fall.2016[fall.2016$AIC == min(unique(fall.2016$AIC)), ])          # mod 43
  winter.2016.df <- data.frame(winter.2016[winter.2016$AIC == min(unique(winter.2016$AIC)), ])  # mod 41
  
  spring.2017.df <- data.frame(spring.2017[spring.2017$AIC == min(unique(spring.2017$AIC)), ])  # mod 41
  summer.2017.df <- data.frame(summer.2017[summer.2017$AIC == min(unique(summer.2017$AIC)), ])  # mod 45
  fall.2017.df <- data.frame(fall.2017[fall.2017$AIC == min(unique(fall.2017$AIC)), ])          # mod 45
  winter.2017.df <- data.frame(winter.2017[winter.2017$AIC == min(unique(winter.2017$AIC)), ])  # mod 40
  
  spring.2018.df <- data.frame(spring.2018[spring.2018$AIC == min(unique(spring.2018$AIC)), ])  # mod 47
  summer.2018.df <- data.frame(summer.2018[summer.2018$AIC == min(unique(summer.2018$AIC)), ])  # mod 47

# make data frames and rbind
  options(scipen = 4, digits = 4)
  all.df <- rbind(summer.2016.df, fall.2016.df, winter.2016.df,
        spring.2017.df, summer.2017.df, fall.2017.df, winter.2017.df,
        spring.2018.df, summer.2018.df)

  names(all.df) <- c("X","Estimate","SE","z","P","AIC","SEASON","MODEL")
  setwd("E:/ANDERSON_BACKUP_26MAR2022/_WORKING_PROJECTS/BALDWIN_PMSP_DEER_PROJECT/REANALYSIS_JAN2022")
  #write.csv(all.df, "All PMSP occupancy results 22MAY2022_1day.csv")
  #all.df <- read.csv("All PMSP occupancy results 22MAY2022_1day.csv")
  
# AIC differences - check all with 1.5 units
  
# pcount(~detection_formula ~occupancy_formula
# 
# INT1 = abundance  
# VAR1 = abundance
# INT2 = detection
# VAR2 = detection
  
### 2016 ###  
  sort(unique(summer.2016$AIC) - min(unique(summer.2016$AIC)))    ### 
  summer.2016[order(summer.2016$AIC), ]                           ### mod33
  summer.2016[summer.2016$AIC == min(unique(summer.2016$AIC)), ]
  
  sort(unique(fall.2016$AIC) - min(unique(fall.2016$AIC)))        ### 
  fall.2016[order(fall.2016$AIC), ]                               ### mod13
  fall.2016[fall.2016$AIC == min(unique(fall.2016$AIC)), ] 
  
  sort(unique(winter.2016$AIC) - min(unique(winter.2016$AIC)))    ###
  winter.2016[order(winter.2016$AIC), ]                           ### mod25
  winter.2016[winter.2016$AIC == min(unique(winter.2016$AIC)), ] 
  
### 2017 ###  
  sort(unique(spring.2017$AIC) - min(unique(spring.2017$AIC)))    ###
  spring.2017[order(spring.2017$AIC), ]                           ### mod26
  spring.2017[spring.2017$AIC == min(unique(spring.2017$AIC)), ] 

  sort(unique(summer.2017$AIC) - min(unique(summer.2017$AIC)))    ###
  summer.2017[order(summer.2017$AIC), ]                           ### mod25
  summer.2017[summer.2017$AIC == min(unique(summer.2017$AIC)), ] 

  sort(unique(fall.2017$AIC) - min(unique(fall.2017$AIC)))        ###
  fall.2017[order(fall.2017$AIC), ]                               ### mod38
  fall.2017[fall.2017$AIC == min(unique(fall.2017$AIC)), ] 

  sort(unique(winter.2017$AIC) - min(unique(winter.2017$AIC)))    ###
  winter.2017[order(winter.2017$AIC), ]                           ### mod22
  winter.2017[winter.2017$AIC == min(unique(winter.2017$AIC)), ] 

### 2018 ###  
  sort(unique(spring.2018$AIC) - min(unique(spring.2018$AIC)))    ###
  spring.2018[order(spring.2018$AIC), ]                           ### mod8
  spring.2018[spring.2018$AIC == min(unique(spring.2018$AIC)), ] 

  sort(unique(summer.2018$AIC) - min(unique(summer.2018$AIC)))    ###
  summer.2018[order(summer.2018$AIC), ]                           ### mod16
  summer.2018[summer.2018$AIC == min(unique(summer.2018$AIC)), ] 

# example of parameter extractions
  #model.list38[[1]]
  #str(model.list38[[1]])
  #model.list38[[1]]@AIC
  #model.list38[[1]]@opt$par
  #model.list38[[1]]@estimates@estimates$state
  #model.list38[[1]]@estimates@estimates$det

# run the source command to run best models in each season from another R scirpt 
  #source(file = , local = FALSE, echo = TRUE)
  
### SUMMER 2016  
  # summer16.df <- summer16[, c('SITE', 'DATE', 'TWELVE.HOUR', 'COUNT')]
  summer16.df <- summer16[, c('SITE', 'DATE', 'COUNT')]
  
  # cast for unmarked  
  # summer16.df <- cast(data = summer16.df, formula = SITE~TWELVE.HOUR, mean)
  summer16.df <- cast(data = summer16.df, formula = SITE~DATE, mean)
  summer16.df <- as.matrix.cast_df(summer16.df[,-1])
  summer16.scaled <- unmarkedFramePCount(y = summer16.df, siteCovs = deer.covars.scaled)
  
  # m33
  summer16.model33 <- pcount(~ELE + EDGE ~SLOPE, data=summer16.scaled, K = 150, se = T,
                             starts = c(0, 0, 0, 0, 0),
                             control = list(trace=T, REPORT = 1))
  
  # final estimates for model
  re.sum16 <- ranef(summer16.model33)
  summer16.EBUP <- bup(re.sum16, stat = "mean")
  summer16.CI <- confint(re.sum16, level = 0.95)
  summer16.CI <- as.data.frame(summer16.CI)
  summer16.EST <- data.frame(Estimate = sum(summer16.EBUP), Lower = sum(summer16.CI$`2.5%`), Upper = sum(summer16.CI$`97.5%`))

### FALL 2016  
  # fall16.df <- fall16[, c('SITE', 'DATE', 'TWELVE.HOUR', 'COUNT')]
  fall16.df <- fall16[, c('SITE', 'DATE', 'COUNT')]
  
  # cast for unmarked  
  fall16.df <- cast(data = fall16.df, formula = SITE~DATE, mean)
  fall16.df <- as.matrix.cast_df(fall16.df[,-1])
  fall16.scaled <- unmarkedFramePCount(y = fall16.df, siteCovs = deer.covars.scaled)
  
  # m47
  fall16.model13 <- pcount(~ASPECT ~SLOPE + ELE, data=fall16.scaled, K = 150, se = T,
                           starts = c(0, 0, 0, 0, 0),
                           control = list(trace=T, REPORT = 1))
  
  # final estimates for model
  re.fall16 <- ranef(fall16.model13)
  fall16.EBUP <- bup(re.fall16, stat = "mean")
  fall16.CI <- confint(re.fall16, level = 0.95)
  fall16.CI <- as.data.frame(fall16.CI)
  fall16.EST <- data.frame(Estimate = sum(fall16.EBUP), Lower = sum(fall16.CI$`2.5%`), Upper = sum(fall16.CI$`97.5%`))
  
### WINTER 2016  
  # winter16.df <- winter16[, c('SITE', 'DATE', 'TWELVE.HOUR', 'COUNT')]
  winter16.df <- winter16[, c('SITE', 'DATE', 'COUNT')]
  
  # cast for unmarked  
  winter16.df <- cast(data = winter16.df, formula = SITE~DATE, mean)
  winter16.df <- as.matrix.cast_df(winter16.df[,-1])
  winter16.scaled <- unmarkedFramePCount(y = winter16.df, siteCovs = deer.covars.scaled)
  
  # m25
  winter16.model25 <- pcount(~SLOPE + ELE ~ASPECT, data=winter16.scaled, K = 150, se = T,
                           starts = c(0, 0, 0, 0, 0),
                           control = list(trace=T, REPORT = 1))
  
  # final estimates for model
  re.winter16 <- ranef(winter16.model25)
  winter16.EBUP <- bup(re.winter16, stat = "mean")
  winter16.CI <- confint(re.winter16, level = 0.95)
  winter16.CI <- as.data.frame(winter16.CI)
  winter16.EST <- data.frame(Estimate = sum(winter16.EBUP), Lower = sum(winter16.CI$`2.5%`), Upper = sum(winter16.CI$`97.5%`))
  
### SPRING 2017 
  # spring17.df <- spring17[, c('SITE', 'DATE', 'TWELVE.HOUR', 'COUNT')]
  spring17.df <- spring17[, c('SITE', 'DATE', 'COUNT')]
  
  # cast for unmarked  
  spring17.df <- cast(data = spring17.df, formula = SITE~DATE, mean)
  spring17.df <- as.matrix.cast_df(spring17.df[,-1])
  spring17.scaled <- unmarkedFramePCount(y = spring17.df, siteCovs = deer.covars.scaled)
  
  # m26
  spring17.model26 <- pcount(~SLOPE + ELE ~EDGE, data=spring17.scaled, K = 150, se = T,
                             starts = c(0, 0, 0, 0, 0),
                             control = list(trace=T, REPORT = 1))
  
  # final estimates for model
  re.spring17 <- ranef(spring17.model26)
  spring17.EBUP <- bup(re.spring17, stat = "mean")
  spring17.CI <- confint(re.spring17, level = 0.95)
  spring17.CI <- as.data.frame(spring17.CI)
  spring17.EST <- data.frame(Estimate = sum(spring17.EBUP), Lower = sum(spring17.CI$`2.5%`), Upper = sum(spring17.CI$`97.5%`))
  
### SUMMER 2017 
  # summer17.df <- summer17[, c('SITE', 'DATE', 'TWELVE.HOUR', 'COUNT')]
  summer17.df <- summer17[, c('SITE', 'DATE', 'COUNT')]
  
  # cast for unmarked  
  summer17.df <- cast(data = summer17.df, formula = SITE~DATE, mean)
  summer17.df <- as.matrix.cast_df(summer17.df[,-1])
  summer17.scaled <- unmarkedFramePCount(y = summer17.df, siteCovs = deer.covars.scaled)
  
  # m25
  summer17.model25 <- pcount(~SLOPE + ELE ~ASPECT, data=summer17.scaled, K = 150, se = T,
                             starts = c(0, 0, 0, 0, 0),
                             control = list(trace=T, REPORT = 1))
  
  # final estimates for model
  re.summer17 <- ranef(summer17.model25)
  summer17.EBUP <- bup(re.summer17, stat = "mean")
  summer17.CI <- confint(re.summer17, level = 0.95)
  summer17.CI <- as.data.frame(summer17.CI)
  summer17.EST <- data.frame(Estimate = sum(summer17.EBUP), Lower = sum(summer17.CI$`2.5%`), Upper = sum(summer17.CI$`97.5%`))
  
### FALL 2017 
  # fall17.df <- fall17[, c('SITE', 'DATE', 'TWELVE.HOUR', 'COUNT')]
  fall17.df <- fall17[, c('SITE', 'DATE', 'COUNT')]
  
  # cast for unmarked  
  fall17.df <- cast(data = fall17.df, formula = SITE~DATE, mean)
  fall17.df <- as.matrix.cast_df(fall17.df[,-1])
  fall17.scaled <- unmarkedFramePCount(y = fall17.df, siteCovs = deer.covars.scaled)
  
  # m38
  fall17.model38 <- pcount(~ASPECT + ELE ~SLOPE + EDGE, data=fall17.scaled, K = 150, se = T,
                             starts = c(0, 0, 0, 0, 0, 0),
                             control = list(trace=T, REPORT = 1))
  
  # final estimates for model
  re.fall17 <- ranef(fall17.model38)
  fall17.EBUP <- bup(re.fall17, stat = "mean")
  fall17.CI <- confint(re.fall17, level = 0.95)
  fall17.CI <- as.data.frame(fall17.CI)
  fall17.EST <- data.frame(Estimate = sum(fall17.EBUP), Lower = sum(fall17.CI$`2.5%`), Upper = sum(fall17.CI$`97.5%`))
  
### WINTER 2017
# flir flight was Feb 2018 - so matches with winter 2017
  # winter17.df <- winter17[, c('SITE', 'DATE', 'TWELVE.HOUR', 'COUNT')]
  winter17.df <- winter17[, c('SITE', 'DATE', 'COUNT')]
  max(winter17.df$DATE)
  min(winter17.df$DATE)
  
  # cast for unmarked  
  # winter17.df <- cast(data = winter17.df, formula = SITE~TWELVE.HOUR, mean)
  winter17.df <- cast(data = winter17.df, formula = SITE~DATE, mean)
  winter17.df <- as.matrix.cast_df(winter17.df[,-1])
  winter17.scaled <- unmarkedFramePCount(y = winter17.df, siteCovs = deer.covars.scaled)
  
  # m22 = pcount(formula = ~SLOPE ~ASPECT + ELE, data = deer.pcount, 
  winter17.model22 <- pcount(~SLOPE ~ASPECT + ELE, data=winter17.scaled, K = 150, se = T,
                             starts = c(0, 0, 0, 0, 0),
                             control = list(trace=T, REPORT = 1))

  # final estimates for model
  re.winter17 <- ranef(winter17.model22)
  winter17.EBUP <- bup(re.winter17, stat = "mean")
  winter17.CI <- confint(re.winter17, level = 0.95)
  winter17.CI <- as.data.frame(winter17.CI)
  winter17.EST <- data.frame(Estimate = sum(winter17.EBUP), Lower = sum(winter17.CI$`2.5%`), Upper = sum(winter17.CI$`97.5%`))
  
  Nmixture.estimate.win17 <- rbind(c(Estimate = sum(winter17.EBUP), colSums(winter17.CI)))
  # Nmixture.estimate.win17 / 10.24
  Nmixture.estimate.win17 / 7
  
  # Plots Winter 2017
  # plot(re.winter17, subset=site %in% c(1:22), layout=c(8, 3), xlim=c(-1,40))  
  
### SPRING 2018
  # spring18.df <- spring18[, c('SITE', 'DATE', 'TWELVE.HOUR', 'COUNT')]
  spring18.df <- spring18[, c('SITE', 'DATE', 'COUNT')]
  
  # cast for unmarked  
  spring18.df <- cast(data = spring18.df, formula = SITE~DATE, mean)
  spring18.df <- as.matrix.cast_df(spring18.df[,-1])
  spring18.scaled <- unmarkedFramePCount(y = spring18.df, siteCovs = deer.covars.scaled)
  
  # m08
  spring18.model08 <- pcount(~ELE ~SLOPE, data=spring18.scaled, K = 150, se = T,
                           starts = c(0, 0, 0, 0),
                           control = list(trace=T, REPORT = 1))
  
  # final estimates for model
  re.spring18 <- ranef(spring18.model08)
  spring18.EBUP <- bup(re.spring18, stat = "mean")
  spring18.CI <- confint(re.spring18, level = 0.95)
  spring18.CI <- as.data.frame(spring18.CI)
  spring18.EST <- data.frame(Estimate = sum(spring18.EBUP), Lower = sum(spring18.CI$`2.5%`), Upper = sum(spring18.CI$`97.5%`))
  
### SUMMER 2018
  # summer18.df <- summer18[, c('SITE', 'DATE', 'TWELVE.HOUR', 'COUNT')]
  summer18.df <- summer18[, c('SITE', 'DATE', 'COUNT')]
  
  # cast for unmarked  
  summer18.df <- cast(data = summer18.df, formula = SITE~DATE, mean)
  summer18.df <- as.matrix.cast_df(summer18.df[,-1])
  summer18.scaled <- unmarkedFramePCount(y = summer18.df, siteCovs = deer.covars.scaled)
  
  # m47
  summer18.model16 <- pcount(~ASPECT ~SLOPE + ELE, data=summer18.scaled, K = 150, se = T,
                             starts = c(0, 0, 0, 0, 0),
                             control = list(trace=T, REPORT = 1))
  
  # final estimates for model
  re.summer18 <- ranef(summer18.model16)
  summer18.EBUP <- bup(re.summer18, stat = "mean")
  summer18.CI <- confint(re.summer18, level = 0.95)
  summer18.CI <- as.data.frame(summer18.CI)
  summer18.EST <- data.frame(Estimate = sum(summer18.EBUP), Lower = sum(summer18.CI$`2.5%`), Upper = sum(summer18.CI$`97.5%`))
  
# bind estimates together
  month.est <- rbind(summer16.EST, fall16.EST, winter16.EST, spring17.EST, summer17.EST, fall17.EST, winter17.EST, spring18.EST, summer18.EST)
  # month.est <- month.est/10.24
  month.est <- month.est/7
  names(month.est) <- c("Estimate", "Lower", "Upper")
  month.est$SEASON <- c('summer16', 'fall16', 'winter16', 'spring17', 'summer17', 'fall17', 'winter17', 'spring18', 'summer18')
  names(month.est) <- c('ESTIMATE', 'LOWER', 'UPPER', 'SEASON')
  month.est <- month.est[, c('SEASON', 'ESTIMATE', 'LOWER', 'UPPER')]
  
# write month.est to file as .csv  
  write.csv(month.est, "Monthly deer density estimates 22MAY2022_1day.csv")
  month.est <- read.csv('Monthly deer density estimates 22MAY2022_1day.csv')
  
# read-in ndvi season data
  setwd("E:/ANDERSON_BACKUP_26MAR2022/_WORKING_PROJECTS/BALDWIN_PMSP_DEER_PROJECT/REANALYSIS_JAN2022")
  ndvi.dat <- read.csv("ndvi_full_period.csv")
  
# convert date to date
  ndvi.dat$DATE <- as.Date(ndvi.dat$DATE, format = "%m/%d/%Y")

# subset ndvi by season
  JUN16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="06" & format.Date(DATE, "%Y")=="2016")
  JUL16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="07" & format.Date(DATE, "%Y")=="2016")
  AUG16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="08"  & format.Date(DATE, "%Y")=="2016")
  OCT16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="10"& format.Date(DATE, "%Y")=="2016")
  SEP16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="09" & format.Date(DATE, "%Y")=="2016")
  NOV16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="11" & format.Date(DATE, "%Y")=="2016")
  DEC16 <- subset(ndvi.dat, format.Date(DATE, "%m")=="12" & format.Date(DATE, "%Y")=="2016")
  JAN17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="01" & format.Date(DATE, "%Y")=="2017")
  FEB17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="02" & format.Date(DATE, "%Y")=="2017")
  MAR17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="03" & format.Date(DATE, "%Y")=="2017")
  APR17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="04" & format.Date(DATE, "%Y")=="2017")
  MAY17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="05" & format.Date(DATE, "%Y")=="2017")
  JUN17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="06" & format.Date(DATE, "%Y")=="2017")
  JUL17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="07" & format.Date(DATE, "%Y")=="2017")
  AUG17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="08" & format.Date(DATE, "%Y")=="2017")
  SEP17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="09" & format.Date(DATE, "%Y")=="2017")
  OCT17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="10" & format.Date(DATE, "%Y")=="2017")
  NOV17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="11" & format.Date(DATE, "%Y")=="2017")
  DEC17 <- subset(ndvi.dat, format.Date(DATE, "%m")=="12" & format.Date(DATE, "%Y")=="2017") 
  JAN18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="01" & format.Date(DATE, "%Y")=="2018")
  FEB18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="02" & format.Date(DATE, "%Y")=="2018")
  MAR18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="03" & format.Date(DATE, "%Y")=="2018")
  APR18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="04" & format.Date(DATE, "%Y")=="2018")
  MAY18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="05" & format.Date(DATE, "%Y")=="2018")
  JUN18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="06" & format.Date(DATE, "%Y")=="2018")
  JUL18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="07" & format.Date(DATE, "%Y")=="2018")
  AUG18 <- subset(ndvi.dat, format.Date(DATE, "%m")=="08" & format.Date(DATE, "%Y")=="2018")
  
# rbind into seasons
  summer16.ndvi <- rbind(JUN16, JUL16, AUG16)
  fall16.ndvi <- rbind(SEP16, OCT16, NOV16)
  winter16.ndvi <- rbind(DEC16, JAN17, FEB17)
  spring17.ndvi <- rbind(MAR17, APR17, MAY17)
  summer17.ndvi <- rbind(JUN17,JUL17, AUG17)
  fall17.ndvi <- rbind(SEP17, OCT17, NOV17)
  winter17.ndvi <- rbind(DEC17, JAN18, FEB18)
  spring18.ndvi <- rbind(MAR18, APR18, MAY18)
  summer18.ndvi <- rbind(JUN18,JUL18, AUG18)

# add factor for season
  summer16.ndvi$SEASON <- rep(1, nrow(summer16.ndvi))
  sum16mean <- mean(summer16.ndvi$NDVI)
  bstrap <- NULL
  for (i in 1:10000){
    bstrap <- c(bstrap, mean(sample(summer16.ndvi$NDVI, nrow(summer16.ndvi),replace=T)))}
  ndvisum16 <- data.frame(Mean = sum16mean, Lower = quantile(bstrap,.025), Upper =quantile(bstrap, .975))
  
  fall16.ndvi$SEASON <- rep(2, nrow(fall16.ndvi))
  fstrap <- NULL
  for (i in 1:10000){
    fstrap <- c(fstrap, mean(sample(fall16.ndvi$NDVI, nrow(fall16.ndvi),replace=T)))}
  ndvifall16 <- data.frame(Mean = mean(fall16.ndvi$NDVI), Lower = quantile(fstrap,.025), Upper =quantile(fstrap, .975))
  
  winter16.ndvi$SEASON <- rep(3, nrow(winter16.ndvi))
  wstrap <- NULL
  for (i in 1:10000){
    wstrap <- c(wstrap, mean(sample(winter16.ndvi$NDVI, nrow(winter16.ndvi),replace=T)))}
  ndviwinter16 <- data.frame(Mean = mean(winter16.ndvi$NDVI), Lower = quantile(wstrap,.025), Upper =quantile(wstrap, .975))
  
  spring17.ndvi$SEASON <-rep(4, nrow(spring17.ndvi))
  sprstrap <- NULL
  for (i in 1:10000){
    sprstrap <- c(sprstrap, mean(sample(spring17.ndvi$NDVI, nrow(spring17.ndvi),replace=T)))}
  ndvispring17 <- data.frame(Mean = mean(spring17.ndvi$NDVI), Lower = quantile(sprstrap,.025), Upper =quantile(sprstrap, .975))
  
  summer17.ndvi$SEASON <- rep(5, nrow(summer17.ndvi))
  sumstrap <- NULL
  for (i in 1:10000){
    sumstrap <- c(sumstrap, mean(sample(summer17.ndvi$NDVI, nrow(summer17.ndvi),replace=T)))}
  ndvisummer17 <- data.frame(Mean = mean(summer17.ndvi$NDVI), Lower = quantile(sumstrap,.025), Upper =quantile(sumstrap, .975))
  
  fall17.ndvi$SEASON <- rep(6, nrow(fall17.ndvi))
  fastrap <- NULL
  for (i in 1:10000){
    fastrap <- c(fastrap, mean(sample(fall17.ndvi$NDVI, nrow(fall17.ndvi),replace=T)))}
  ndvifall17 <- data.frame(Mean = mean(fall17.ndvi$NDVI), Lower = quantile(fastrap,.025), Upper =quantile(fastrap, .975))
  
  winter17.ndvi$SEASON <- rep(7, nrow(winter17.ndvi))
  wstrap <- NULL
  for (i in 1:10000){
    wstrap <- c(wstrap, mean(sample(winter17.ndvi$NDVI, nrow(winter17.ndvi),replace=T)))}
  ndviwinter17 <- data.frame(Mean = mean(winter17.ndvi$NDVI), Lower = quantile(wstrap,.025), Upper =quantile(wstrap, .975))
  
  spring18.ndvi$SEASON <-rep(8, nrow(spring18.ndvi))
  sprstrap <- NULL
  for (i in 1:10000){
    sprstrap <- c(sprstrap, mean(sample(spring18.ndvi$NDVI, nrow(spring18.ndvi),replace=T)))}
  ndvispring18 <- data.frame(Mean = mean(spring18.ndvi$NDVI), Lower = quantile(sprstrap,.025), Upper =quantile(sprstrap, .975))
  
  summer18.ndvi$SEASON <- rep(9, nrow(summer18.ndvi))
  sumstrap <- NULL
  for (i in 1:10000){
    sumstrap <- c(sumstrap, mean(sample(summer18.ndvi$NDVI, nrow(summer18.ndvi),replace=T)))}
  ndvisummer18 <- data.frame(Mean = mean(summer18.ndvi$NDVI), Lower = quantile(sumstrap,.025), Upper =quantile(sumstrap, .975))
  
# rebind
  ndvi.dat <- rbind(summer16.ndvi, fall16.ndvi, winter16.ndvi, spring17.ndvi, summer17.ndvi, fall17.ndvi, winter17.ndvi, spring18.ndvi, summer18.ndvi)
  ndvi.plot <- rbind(ndvisum16, ndvifall16, ndviwinter16, ndvispring17, ndvisummer17, ndvifall17, ndviwinter17, ndvispring18, ndvisummer18)

#  ----------------------------------------- N-MIXTURE MODELLING DEER ABUNDANCE DEC 2017 - FEB 2018 --------------------------------------------------- #

# try this just for the month of February, thus the date:
  # deer.df.feb18 <- rbind(JAN18, FEB18)
  deer.df.feb18 <- deer.df[deer.df$DATE > as.Date("2018-01-15") & deer.df$DATE < as.Date("2018-03-15"), ]
  #deer.df.feb18 <- deer.df[deer.df$DATE > as.Date("2017-12-01") & deer.df$DATE < as.Date("2018-02-28"), ]
  
# cast df
  deer.df.feb18 <- cast(data = deer.df.feb18, formula = SITE~DATE, mean)
  deer.df.feb18 <- as.matrix.cast_df(deer.df.feb18[,-1])
  
# create unmarked frames with data 
  deer.df.feb18 <- unmarkedFramePCount(y = deer.df.feb18, siteCovs = deer.covars.scaled)
  
# run model
  feb18.mod <- pcount(~SLOPE ~ASPECT + ELE, data=deer.df.feb18, K = 150, se = T,
         starts = c(0, 0, 0, 0, 0),
         control = list(trace=T, REPORT = 1))
  
# extract EBUP  
  re.feb18 <- ranef(feb18.mod)
  feb18.EBUP <- bup(re.feb18, stat = "mean")
  feb18.CI <- confint(re.feb18, level = 0.95)
  feb18.CI <- as.data.frame(feb18.CI)
  feb18.EST <- data.frame(Estimate = sum(feb18.EBUP), Lower = sum(feb18.CI$`2.5%`), Upper = sum(feb18.CI$`97.5%`))
  FEB18.est <- feb18.EST/7
  
# run model on Winter 2017 data with best model at different aggregations 
# as a sensitivity analysis for survey length
  
# the one-day agg should equal this data set:  
  #winter17.scaled
  
# the winter 2017 survey began on Dec 1st and ended 28th Feb, thus the date:
  min.date.time <- as.POSIXct(strptime("2017-12-01 00:00:00", "%Y-%m-%d %H:%M:%S"))
  max.date.time <- as.POSIXct(strptime("2018-02-28 23:00:00", "%Y-%m-%d %H:%M:%S"))
  seq.date.time <- seq.POSIXt(from = min.date.time, to = max.date.time, by = "hour")
  rep.date.time <- rep(seq.date.time, length(all.sites))

# create site reps for new data frame
  rep.sites <- sort(rep(all.sites, length(seq.date.time)))
  
# create data frame for merge
  survey.sites.df <- data.frame(SITE = rep.sites, DATE.HOUR = rep.date.time)
  head(survey.sites.df)
  tail(survey.sites.df)
  
# merge site by survey df with aggregated count data keeping ALL the survey sites
  deer.df.win17 <- merge(x = survey.sites.df, y = count.by.hour, by.x = c("SITE","DATE.HOUR"), by.y = c("SITE","DATE.HOUR"), all.x = T, all.y = F)
  deer.df.win17[is.na(deer.df.win17)] <- 0
  head(deer.df.win17)
  tail(deer.df.win17)
 
# trim dataset so it is exactly 90 days
  (days.on <- as.numeric(round(max.date.time - min.date.time, 0))) # 90
  temp <- deer.df.win17[order(deer.df.win17$SITE, deer.df.win17$DATE.HOUR),]
  length(temp$DATE.HOUR) # 47520
  temp$DATE <- as.Date.POSIXct(temp$DATE.HOUR, tz = '', format = "%Y-%m-%d")
  temp$YEAR.MONTH <- paste(format(temp$DATE, "%Y"),format(temp$DATE, "%m"), sep='.')
  tail(temp, 50)
  str(temp)
 
# make duplicate dataframes for sensitivity analysis of survey length
  temp.1 <- temp
  temp.2 <- temp
  temp.3 <- temp
  temp.4 <- temp
  temp.6 <- temp
  temp.8 <- temp
  temp.12 <- temp
  temp.day <- temp
  temp.48 <- temp
  temp.5day <- temp
  temp.week <- temp

  # add new column to aggregate data for modelling
  # length(sort(rep(1:blocks, 22*24*5))) # rep(1:number of blocks, sites * hours * days)
  #12 weeks up to Feb 22nd, then one 6-day window up until Feb 28
  
  mean.x <- function(x){mean(x, na.rm=TRUE)}

# 7-day  
  tmp1 <- sort(rep(seq.Date(as.Date(min.date.time), as.Date("2018-02-22 23:00:00 EST"), by = 'day'), 22))
  tmp2 <- rep(all.sites, length(tmp1)/22)
  tmp3 <- data.frame(SITE=tmp2, DATE=tmp1)
  tmp3 <- tmp3[order(tmp3$SITE, tmp3$DATE),]
  row.names(tmp3) <- 1:length(tmp3$SITE)
  tmp3$AGG <- sort(rep(1:264, 7))
  
  tmp4 <- sort(rep(seq.Date(as.Date("2018-02-23"), as.Date("2018-02-28"), by = 'day'), 22))
  tmp5 <- rep(all.sites, length(tmp4)/22)
  tmp6 <- data.frame(SITE=tmp5, DATE=tmp4)
  tmp6 <- tmp6[order(tmp6$SITE, tmp6$DATE),]
  row.names(tmp6) <- 1:length(tmp6$SITE)
  tmp6$AGG <- sort(rep(265:286, 6))
  tmp7 <- rbind(tmp3, tmp6)
  tmp7 <- tmp7[order(tmp7$SITE, tmp7$DATE),]
  row.names(tmp7) <- 1:length(tmp7$SITE)
  
  temp.week <- merge(temp.week, tmp7, by = c('SITE', 'DATE'), all = TRUE)
  temp.weeknew <- aggregate(x = temp.week$COUNT, by = list(temp.week$SITE, temp.week$AGG), FUN = max)
  names(temp.weeknew) <- c('SITE', 'AGG', 'COUNT')
  temp.weeknew <- temp.weeknew[order(temp.weeknew$SITE, temp.weeknew$AGG),]
  temp.weeknew$WEEK <- rep(1:13, 22)
  row.names(temp.weeknew) <- 1:length(temp.weeknew$SITE)
  
  count.week <- cast(data = temp.weeknew, formula = SITE~WEEK, fun.aggregate = mean.x, value = 'COUNT')
  
# 5-day 
  tmp1 <- sort(rep(seq.Date(as.Date(min.date.time), as.Date("2018-02-28 23:00:00 EST"), by = 'day'), 22))
  tmp2 <- rep(all.sites, length(tmp1)/22)
  tmp3 <- data.frame(SITE=tmp2, DATE=tmp1)
  tmp3 <- tmp3[order(tmp3$SITE, tmp3$DATE),]
  row.names(tmp3) <- 1:length(tmp3$SITE)
  length(tmp3$SITE)
  tmp3$AGG <- sort(rep(1:396, 5))
  
  temp.5day <- merge(temp.5day, tmp3, by = c('SITE', 'DATE'), all = TRUE)
  temp.5daynew <- aggregate(x = temp.5day$COUNT, by = list(temp.5day$SITE, temp.5day$AGG), FUN = max)
  names(temp.5daynew) <- c('SITE', 'AGG', 'COUNT')
  temp.5daynew <- temp.5daynew[order(temp.5daynew$SITE, temp.5daynew$AGG),]
  temp.5daynew$FIVE.DAY <- rep(1:18, 22)
  row.names(temp.5daynew) <- 1:length(temp.5daynew$SITE)
  
  count.5day <- cast(data = temp.5daynew, formula = SITE~FIVE.DAY, fun.aggregate = mean.x, value = 'COUNT')
  
# 2-day 
  tmp1 <- sort(rep(seq.Date(as.Date(min.date.time), as.Date("2018-02-28 23:00:00 EST"), by = 'day'), 22))
  tmp2 <- rep(all.sites, length(tmp1)/22)
  tmp3 <- data.frame(SITE=tmp2, DATE=tmp1)
  tmp3 <- tmp3[order(tmp3$SITE, tmp3$DATE),]
  row.names(tmp3) <- 1:length(tmp3$SITE)
  length(tmp3$SITE)
  tmp3$AGG <- sort(rep(1:990, 2))
  
  temp.48 <- merge(temp.48, tmp3, by = c('SITE', 'DATE'), all = TRUE)
  temp.48new <- aggregate(x = temp.48$COUNT, by = list(temp.48$SITE, temp.48$AGG), FUN = max)
  names(temp.48new) <- c('SITE', 'AGG', 'COUNT')
  temp.48new <- temp.48new[order(temp.48new$SITE, temp.48new$AGG),]
  temp.48new$TWO.DAY <- rep(1:45, 22)
  row.names(temp.48new) <- 1:length(temp.48new$SITE)
  
  count.48 <- cast(data = temp.48new, formula = SITE~TWO.DAY, fun.aggregate = mean.x, value = 'COUNT')

# 1-day  
  temp.day$AGG <- sort(rep(1:1980, 24))   # 22 sites * 90 days, repeated for 24 hours
  temp.dayagg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.day, FUN = "max")
  count.day <- cast(data = temp.dayagg, formula = SITE~DATE, value = 'COUNT')
  
# 12-hour  
  tmp12h <- sort(rep(1:3960, 12)); length(tmp12h) # 22 sites * 90 days * 2 half days, repeated for 12 hours
  temp.12$AGG <- tmp12h[1:length(temp.12$SITE)]
  temp.12agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.12, FUN = "max")
  temp.12agg <- temp.12agg[order(temp.12agg$SITE, temp.12agg$AGG),]
  row.names(temp.12agg) <- 1:length(temp.12agg$SITE)
  temp.12agg$TWELVE.HOUR <- rep(1:180, 22)
  count.12 <- cast(data = temp.12agg, formula = SITE~TWELVE.HOUR, value = 'COUNT')
  
# 8-hour  
  tmp8h <- sort(rep(1:5940, 8)); length(tmp8h) # 22 sites * 90 * 3 8-hour periods, repeated for 8 hours
  temp.8$AGG <- tmp8h[1:length(temp.8$SITE)]
  temp.8agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.8, FUN = "max")
  temp.8agg <- temp.8agg[order(temp.8agg$SITE, temp.8agg$AGG),]
  row.names(temp.8agg) <- 1:length(temp.8agg$SITE)
  temp.8agg$EIGHT.HOUR <- rep(1:270, 22)
  count.8 <- cast(data = temp.8agg, formula = SITE~EIGHT.HOUR, value = 'COUNT')
  
# 6-hour  
  tmp6h <- sort(rep(1:7920, 6)); length(tmp6h) # 22 sites * 90 * 4 6-hour periods, repeated for 6 hours
  temp.6$AGG <- tmp6h[1:length(temp.6$SITE)]
  temp.6agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.6, FUN = "max")
  temp.6agg <- temp.6agg[order(temp.6agg$SITE, temp.6agg$AGG),]
  row.names(temp.6agg) <- 1:length(temp.6agg$SITE)
  temp.6agg$SIX.HOUR <- rep(1:360, 22)
  count.6 <- cast(data = temp.6agg, formula = SITE~SIX.HOUR, value = 'COUNT')
  
# 4-hour  
  tmp4h <- sort(rep(1:11880, 4)); length(tmp4h) # 22 sites * 90 * 6 4-hour periods, repeated for 4 hours
  temp.4$AGG <- tmp4h[1:length(temp.4$SITE)]
  temp.4agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.4, FUN = "max")
  temp.4agg <- temp.4agg[order(temp.4agg$SITE, temp.4agg$AGG),]  
  row.names(temp.4agg) <- 1:length(temp.4agg$SITE)
  temp.4agg$FOUR.HOUR <- rep(1:540, 22)
  count.4 <- cast(data = temp.4agg, formula = SITE~FOUR.HOUR, value = 'COUNT')
  
# 3-hour  
  tmp3h <- sort(rep(1:15840, 3)); length(tmp3h) # 22 sites * 90 * 8 3-hour periods, repeated for 3 hours
  temp.3$AGG <- tmp3h[1:length(temp.3$SITE)]
  temp.3agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.3, FUN = "max")
  temp.3agg <- temp.3agg[order(temp.3agg$SITE, temp.3agg$AGG),]  
  row.names(temp.3agg) <- 1:length(temp.3agg$SITE)
  temp.3agg$THREE.HOUR <- rep(1:720, 22)
  count.3 <- cast(data = temp.3agg, formula = SITE~THREE.HOUR, value = 'COUNT')
  
# 2-hour  
  tmp2h <- sort(rep(1:23760, 2)); length(tmp2h) # 22 sites * 90 * 12 2-hour periods, repeated for 2 hours
  temp.2$AGG <- tmp2h[1:length(temp.2$SITE)]
  temp.2agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.2, FUN = "max")
  temp.2agg <- temp.2agg[order(temp.2agg$SITE, temp.2agg$AGG),]  
  row.names(temp.2agg) <- 1:length(temp.2agg$SITE)
  temp.2agg$TWO.HOUR <- rep(1:1080, 22)
  count.2 <- cast(data = temp.2agg, formula = SITE~TWO.HOUR, value = 'COUNT')
  
# 1-hour  
  tmp1h <- sort(rep(1:47520, 1)) # 24, 1-hour blocks per day
  temp.1$AGG <- tmp1h[1:length(temp.1$SITE)]
  temp.1agg <- aggregate(COUNT ~ SITE + AGG + DATE, data = temp.1, FUN = "max")
  temp.1agg <- temp.1agg[order(temp.1agg$SITE, temp.1agg$AGG),]  
  row.names(temp.1agg) <- 1:length(temp.1agg$SITE)
  temp.1agg$ONE.HOUR <- rep(1:2160, 22)
  count.1 <- cast(data = temp.1agg, formula = SITE~ONE.HOUR, value = 'COUNT')
 
# convert to matrix and remove first column
  count.1 <- as.matrix.cast_df(count.1,-1)
  count.2 <- as.matrix.cast_df(count.2,-1)
  count.3 <- as.matrix.cast_df(count.3,-1)
  count.4 <- as.matrix.cast_df(count.4,-1)
  count.6 <- as.matrix.cast_df(count.6,-1)
  count.8 <- as.matrix.cast_df(count.8,-1)
  count.12 <- as.matrix.cast_df(count.12,-1)
  count.day <- as.matrix.cast_df(count.day,-1)
  count.48 <- as.matrix.cast_df(count.48,-1)
  count.5day <- as.matrix.cast_df(count.5day,-1)
  count.week <- as.matrix.cast_df(count.week,-1)

# clear memory 
  gc()
  memory.limit()
  memory.limit(size=200000)
  
# create unmarked frame with data
  pcount.list <- list(count.1, count.2, count.3, count.4, count.6, count.8, count.12, count.day, count.48, count.5day, count.week)
  pcount.list <- lapply(pcount.list, function(x) unmarkedFramePCount(y = x, siteCovs = deer.covars.scaled))
  #deer.pcount <- unmarkedFramePCount(y = count.df, siteCovs = deer.covars.scaled)
  #deer.day <- unmarkedFramePCount(y = count.day, siteCovs = deer.covars.scaled)
  
  # model of Royal 2004b based on counts
  # general formula is: model <- pcount(~detection_formula ~occupancy_formula, dataframe, K=100, se=TRUE)
  # starting values generated by iteratively running models
  
  # run null model on different aggregations as sensitivity analysis for survey length
  # this takes a long time...
  model.list <- lapply(pcount.list, function(x) pcount(~SLOPE ~ASPECT + ELE, x, K = 150, se = T, 
                                                       starts = c(2.6, -0.2, 0.2, -4.7, -0.3), 
                                                       control = list(trace=T, REPORT = 1)) )  
  
# repeat process on the sensitivity analysis models
  re.list <- lapply(model.list, function(x) ranef(x))
  EBUP.list <- lapply(re.list, function (x) bup(x, stat = 'mean'))
  CI.list <- lapply(re.list, function(x) confint(x, level = 0.95))
  
  hour1.EBUP <- unlist(EBUP.list[[1]])
  hour1.CI <- unlist(CI.list[[1]])
  hour1.CI <- as.data.frame(hour1.CI)
  hour1.EST <- data.frame(Estimate = sum(hour1.EBUP), Lower = sum(hour1.CI$`2.5%`), Upper = sum(hour1.CI$`97.5%`))
  
  hour2.EBUP <- unlist(EBUP.list[[2]])
  hour2.CI <- unlist(CI.list[[2]])
  hour2.CI <- as.data.frame(hour2.CI)
  hour2.EST <- data.frame(Estimate = sum(hour2.EBUP), Lower = sum(hour2.CI$`2.5%`), Upper = sum(hour2.CI$`97.5%`))
  
  hour3.EBUP <- unlist(EBUP.list[[3]])
  hour3.CI <- unlist(CI.list[[3]])
  hour3.CI <- as.data.frame(hour3.CI)
  hour3.EST <- data.frame(Estimate = sum(hour3.EBUP), Lower = sum(hour3.CI$`2.5%`), Upper = sum(hour3.CI$`97.5%`))
  
  hour4.EBUP <- unlist(EBUP.list[[4]])
  hour4.CI <- unlist(CI.list[[4]])
  hour4.CI <- as.data.frame(hour4.CI)
  hour4.EST <- data.frame(Estimate = sum(hour4.EBUP), Lower = sum(hour4.CI$`2.5%`), Upper = sum(hour4.CI$`97.5%`))
  
  hour6.EBUP <- unlist(EBUP.list[[5]])
  hour6.CI <- unlist(CI.list[[5]])
  hour6.CI <- as.data.frame(hour6.CI)
  hour6.EST <- data.frame(Estimate = sum(hour6.EBUP), Lower = sum(hour6.CI$`2.5%`), Upper = sum(hour6.CI$`97.5%`))
  
  hour8.EBUP <- unlist(EBUP.list[[6]])
  hour8.CI <- unlist(CI.list[[6]])
  hour8.CI <- as.data.frame(hour8.CI)
  hour8.EST <- data.frame(Estimate = sum(hour8.EBUP), Lower = sum(hour8.CI$`2.5%`), Upper = sum(hour8.CI$`97.5%`))
  
  hour12.EBUP <- unlist(EBUP.list[[7]])
  hour12.CI <- unlist(CI.list[[7]])
  hour12.CI <- as.data.frame(hour12.CI)
  hour12.EST <- data.frame(Estimate = sum(hour12.EBUP), Lower = sum(hour12.CI$`2.5%`), Upper = sum(hour12.CI$`97.5%`))
  
  day.EBUP <- unlist(EBUP.list[[8]])
  day.CI <- unlist(CI.list[[8]])
  day.CI <- as.data.frame(day.CI)
  day.EST <- data.frame(Estimate = sum(day.EBUP), Lower = sum(day.CI$`2.5%`), Upper = sum(day.CI$`97.5%`))
  
  hour48.EBUP <- unlist(EBUP.list[[9]])
  hour48.CI <- unlist(CI.list[[9]])
  hour48.CI <- as.data.frame(hour48.CI)
  hour48.EST <- data.frame(Estimate = sum(hour48.EBUP), Lower = sum(hour48.CI$`2.5%`), Upper = sum(hour48.CI$`97.5%`))
  
  day5.EBUP <- unlist(EBUP.list[[10]])
  day5.CI <- unlist(CI.list[[10]])
  day5.CI <- as.data.frame(day5.CI)
  day5.EST <- data.frame(Estimate = sum(day5.EBUP), Lower = sum(day5.CI$`2.5%`), Upper = sum(day5.CI$`97.5%`))
  
  week.EBUP <- unlist(EBUP.list[[11]])
  week.CI <- unlist(CI.list[[11]])
  week.CI <- as.data.frame(week.CI)
  week.EST <- data.frame(Estimate = sum(week.EBUP), Lower = sum(week.CI$`2.5%`), Upper = sum(week.CI$`97.5%`))
  
  window.est <- rbind(hour1.EST, hour2.EST, hour3.EST, hour4.EST, hour6.EST, hour8.EST, hour12.EST, day.EST, hour48.EST, day5.EST, week.EST)
  window.est$Window <- c('1 hour','2 hour','3 hour', '4 hour', '6 hour', '8 hour', '12 hour', '24 hour', '48hour', '5day', '1 week')
  
  plot(window.est$Estimate/10.24 ~ 1)
  
  #extract coefficients
  state.list <- lapply(model.list, function(x) coef(x, type = 'state'))
  det.list <- lapply(model.list, function(x) coef(x, type = 'det'))
  
  # ext state
  hour1.state <- unlist(state.list[[1]])
  hour2.state <- unlist(state.list[[2]])
  hour3.state <- unlist(state.list[[3]])
  hour4.state <- unlist(state.list[[4]])
  hour6.state <- unlist(state.list[[5]])
  hour8.state <- unlist(state.list[[6]])  
  hour12.state <- unlist(state.list[[7]])
  day.state <- unlist(state.list[[8]])
  hour48.state <- unlist(state.list[[9]])
  day5.state <- unlist(state.list[[10]])
  week.state <- unlist(state.list[[11]])
  state.df <- rbind(hour1.state, hour2.state, hour3.state, hour4.state, hour6.state, hour8.state,
                    hour12.state, day.state, hour48.state, day5.state, week.state)
  
  # extract det
  hour1.det <- unlist(det.list[[1]])
  hour2.det <- unlist(det.list[[2]])
  hour3.det <- unlist(det.list[[3]])
  hour4.det <- unlist(det.list[[4]])
  hour6.det <- unlist(det.list[[5]])
  hour8.det <- unlist(det.list[[6]])  
  hour12.det <- unlist(det.list[[7]])
  day.det <- unlist(det.list[[8]])
  hour48.det <- unlist(det.list[[9]])
  day5.det <- unlist(det.list[[10]])
  week.det <- unlist(det.list[[11]])
  det.df <- rbind(hour1.det, hour2.det, hour3.det, hour4.det, hour6.det, hour8.det,
                  hour12.det, day.det, hour48.det, day5.det, week.det)
  
  window.est <- data.frame(window.est, state.df, det.df)
  window.est$TIME.BLOCK <- row.names(window.est)
  window.est$HOURS <- c(1,2,3,4,6,8,12,24,(2*24),(5*24),(7*24))
  row.names(window.est) <- seq(1:length(window.est$Estimate))
  
  names(window.est) <- c("ESTIMATE","LOWER","UPPER","WINDOW","OCC.INT","ASPECT",
                         "ELE","DET.INT","SLOPE","TIME.BLOCK","HOURS")
  
  window.est$DENSITY <- window.est$ESTIMATE/10.24
  window.est$DENSITY.L <- window.est$LOWER/10.24
  window.est$DENSITY.U <- window.est$UPPER/10.24
  
  getwd()
  #write.csv(window.est, "Variable.time.windows.occupancy.models.csv")
  #window.est <- read.csv("Variable.time.windows.occupancy.models.csv")
  window.est <- read.csv("Variable.time.windows.occupancy.models.csv")
  summary(lm(window.est$DENSITY ~ log(window.est$HOURS)))

# ------------------------------------- MAKE FIGURES ---------------------------- #  

# make figure 2: deer density vs method
  z <- seq(0, 50, 10)
  yl1 <- "Deer Density Estimate" 
  yl2 <- expression("(Deer" ~ km^2 ~ ")")
  
  # with three points
  nmm.df <- data.frame(TYPE = 'NMM', MEAN = FEB18.est$Estimate, SD = NA, LOWER = FEB18.est$Lower, UPPER = FEB18.est$Upper)
  fig2.df <- rbind(fl.df, mr.df, nmm.df)
  fig2.df$TYPE <- factor(fig2.df$TYPE, levels=c('MR', 'NMM', 'UAI'))
  fig2.df$ORDER <- as.numeric(fig2.df$TYPE)
  fig2.df <- fig2.df[order(fig2.df$ORDER),]
  
  # with four points
  # win17.df <- month.est[month.est$SEASON=='winter17', ]
  # spr18.df <- month.est[month.est$SEASON=='spring18', ]
  # nmm.df.17 <- data.frame(TYPE = 'NMM-Win17', MEAN = win17.df$ESTIMATE, SD = NA, LOWER = win17.df$LOWER, UPPER = win17.df$UPPER)
  # nmm.df.18 <- data.frame(TYPE = 'NMM-Spr18', MEAN = spr18.df$ESTIMATE, SD = NA, LOWER = spr18.df$LOWER, UPPER = spr18.df$UPPER)
  # fig2.df <- rbind(fl.df, mr.df, nmm.df.17, nmm.df.18)
  # fig2.df$TYPE <- factor(fig2.df$TYPE, levels=c('MR', 'NMM-Win17', 'NMM-Spr18', 'UAI'))
  # fig2.df$ORDER <- as.numeric(fig2.df$TYPE)
  # fig2.df <- fig2.df[order(fig2.df$ORDER),]
  
  png(file="Deer density vs method 19MAY2022_1day.png",
      width=1100, height=850)
  
  par(mar=c(10,10,5,5)) #par(mar = c(bottom, left, top, right)
  
  plot(fig2.df$ORDER, fig2.df$MEAN, type = 'p', pch = 21, col = 'black', bg = c('grey','black', 'black'), 
       axes = FALSE, cex = 4, ylim = c(0, 50), xlim = c(0.8, 3.2), ylab = '', xlab = '')
  #plot(fig2.df$ORDER, fig2.df$MEAN, type = 'p', pch = 21, col = 'black', bg = c('grey','black', 'black', 'black'), 
  #     axes = FALSE, cex = 4, ylim = c(0,50), xlim = c(0.8,4.2), ylab = '', xlab = '')
  
  mtext(1, text = 'Method', cex = 3, line = 8)
  mtext(2, text = yl1, cex = 3, line = 8)
  mtext(2, text = yl2, cex = 3, line = 4)
  
  axis(side = 1, at = as.numeric(fig2.df$TYPE), labels = fig2.df$TYPE, cex.axis=3, line = 0, padj = 1)
  axis(side = 2, at = z, cex.axis = 3)
  
  arrows(1, fig2.df$MEAN[1], 1, fig2.df$UPPER[1], length = 0.05, angle = 90, lwd = 4, lty=1, col='grey')
  arrows(1, fig2.df$MEAN[1], 1, fig2.df$LOWER[1], length = 0.05, angle = 90, lwd = 4, lty=1, col='grey')
  points(fig2.df$ORDER, fig2.df$MEAN, type='p', pch = 21, col = 'black', bg = c('black', 'black', 'black'), cex = 4)
  #points(fig2.df$ORDER, fig2.df$MEAN, type='p', pch = 21, col = 'black', bg = c('black', 'black', 'black', 'black'), cex = 4)
  
  arrows(2, fig2.df$MEAN[2], 2, fig2.df$UPPER[2], length = 0.05, angle = 90, lwd = 4)
  arrows(2, fig2.df$MEAN[2], 2, fig2.df$LOWER[2], length = 0.05, angle = 90, lwd = 4)
  
  arrows(3, fig2.df$MEAN[3], 3, fig2.df$UPPER[3], length = 0.05, angle = 90, lwd = 4)
  arrows(3, fig2.df$MEAN[3], 3, fig2.df$LOWER[3], length = 0.05, angle = 90, lwd = 4)
  
  #arrows(4, fig2.df$MEAN[4], 4, fig2.df$UPPER[4], length = 0.05, angle = 90, lwd = 4)
  #arrows(4, fig2.df$MEAN[4], 4, fig2.df$LOWER[4], length = 0.05, angle = 90, lwd = 4)
  
  dev.off()

# make figure 3: deer density vs survey length 
  
  yl1 <- "Deer Density Estimate" 
  yl2 <- expression("(Deer" ~ km^2 ~ ")")
  
  png(file="Deer density vs survey length 18MAY2022.png",
      width=1100, height=850)

  par(mar=c(10,10,5,5)) #par(mar = c(bottom, left, top, right)
  
  plot(window.est$DENSITY ~ log(window.est$HOURS), ylim = c(10,45), type = 'b', pch = 19, cex = 4, axes = FALSE, xlab = '', ylab = '', lwd = 3)
  
  axis(1, at = 0:5, labels = 0:5, cex.axis = 3, line = 0, padj = 1)
  axis(2, at = seq(10, 45, 5), cex.axis = 3, labels = seq(10, 45, 5), padj = -0.5)

  mtext(1, text = 'Survey Time (log[hours])', cex = 3, line = 6)
  mtext(2, text = yl1, cex = 3, line = 8)
  mtext(2, text = yl2, cex = 3, line = 4.5)
  
  rect(-0.19, 23.7, 5.32, 40, border = NA, col = 'grey86', lwd=0)
  
  abline(h = 31, lty = 2, col = 'grey50', lwd=3)
  
  points(window.est$DENSITY ~ log(window.est$HOURS), ylim=c(20,45), type = 'b', pch=19, cex = 4, lwd=2)
  
  arrows(x0 = log(window.est$HOURS), x1 = log(window.est$HOURS), y0 = window.est$DENSITY, y1 = window.est$DENSITY.U, angle = 90, length = 0.05, lwd = 3)
  arrows(x0 = log(window.est$HOURS), x1 = log(window.est$HOURS), y0 = window.est$DENSITY, y1 = window.est$DENSITY.L, angle = 90, length = 0.05, lwd = 3)

  box()
  
  dev.off()
  
# make figure 4: deer density vs NDVI
  
  # define labels for axis
  labs <- c('Summer 2016', 'Fall 2016', 'Winter 2016', 'Spring 2017', 'Summer 2017', 'Fall 2017', 'Winter 2017', 'Spring 2018', 'Summer 2018')
  
  # define x and y axis text
  yl1 <- "Deer Density Estimate" 
  yl2 <- expression("(Deer" ~ km^2 ~ ")")
  
  # get month est of NDVI
  y <- month.est$ESTIMATE
  ndvi <- ndvi.plot$Mean / 10000
  ndvi <- ndvi*40
  ndvi <- matrix(ndvi,1,9)
  
  x.ticks <- seq(0.5, 8.5, 1)
  
  png(file="Deer density vs NDVI 22MAY2022_1day.png",
      width=1100, height=850)
  
  par(mar=c(10,10,5,8)) #par(mar = c(bottom, left, top, right)
  
  # barplot
  df.bar <- barplot(ndvi, xlab = '', ylab = '', ylim = c(0, 60), col = 'darkseagreen1', space = 0, width = 1, cex.axis = 3)
  
  # axis labels
  axis(side = 1, at = x.ticks, tick = TRUE, labels = FALSE)
  axis(side = 4, cex.axis = 3, col.axis = 'forest green', labels = c(0, 0.25, 0.5, 0.75, 1), at = c(0, 12.5, 25, 37.5, 50), line = 0, padj = 1)
  
  mtext("NDVI", side = 4, cex = 3, line = 6)
  mtext(1, text = 'Season', cex = 3, line = 8)
  mtext(2, text = yl1, cex = 3, line = 8)
  mtext(2, text = yl2, cex = 3, line = 4)
  
  # points
  points(y ~ df.bar,
         pch = 19,
         ylim = c(0,40),
         xaxt = 'n',
         cex = 4)
  
  text(x = x.ticks,
       y = par("usr")[3] - 2,
       labels = labs,
       xpd = NA,
       srt = 15,
       cex = 1.5)
  
  arrows(x0 = df.bar, x1 = df.bar, y0 = month.est$ESTIMATE, y1 = month.est$UPPER, angle = 90, length = 0.05, lwd = 3)
  arrows(x0 = df.bar, x1 = df.bar, y0 = month.est$ESTIMATE, y1 = month.est$LOWER, angle = 90, length = 0.05, lwd = 3)
  
  arrows(x0 = 7, x1 = 7, y0 = 31, y1 = 40, angle = 90, length = 0.05, lwd = 3)  
  arrows(x0 = 7, x1 = 7, y0 = 31, y1 = 23.7, angle = 90, length = 0.05, lwd = 3) 
  
  lines(df.bar, y, col='black', lwd=3)
  
  points(7, 31,
         pch = 19,
         cex = 4,
         col = 'dodgerblue')
  
  dev.off()

# make additional figure of detection vs survey window length 
  
  yl3 <- "Detection Intercept" 
  
  png(file="Detection intercept vs survey length 23MAR2022.png",
      width=1100, height=850)
  
  par(mar=c(10,10,5,5)) #par(mar = c(bottom, left, top, right)
  
  plot(window.est$DET.INT ~ log(window.est$HOURS), ylim = c(-7,0), type = 'p', pch = 19, cex = 4, axes = FALSE, xlab = '', ylab = '', lwd = 3)
  
  axis(1, at = 0:5, labels = 0:5, cex.axis = 3, line = 0, padj = 1)
  axis(2, at = seq(-7, 0, 0.5), cex.axis = 3, labels = seq(-7, 0, 0.5), padj = -0.5)
  
  mtext(1, text = 'Survey Time (log[hours])', cex = 3, line = 6)
  mtext(2, text = yl3, cex = 3, line = 6)

  points(window.est$DET.INT ~ log(window.est$HOURS), ylim=c(-7,0), type = 'p', pch=19, cex = 4, lwd=2)

  box()
  
  dev.off()
  
  # # add precip data  
  # prcp.dat <- read.csv('PMSP precip data 2016-2018.csv')
  # prcp.dat$Ordernew <- prcp.dat$Order-0.5
  # points(prcp.dat$PRCP_cm ~ prcp.dat$Ordernew, pch=19, type='b', col='dodgerblue')
  # 
  # plot(prcp.dat$PRCP_cm ~ month.est$ESTIMATE, pch=19, col='black', cex=1.5)
  
### EOF