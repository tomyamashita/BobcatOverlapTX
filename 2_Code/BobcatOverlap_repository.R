# Bobcat Spatial Overlap

rm(list = ls()); gc()

## Install required packages
devtools::install_github("ctmm-initiative/ctmm")


#############################################################################################################################
################################################# Part 1: Prepping the data ################################################# 
#############################################################################################################################

rm(list = ls()); gc()

# Step 1: Loading the data ####
## Load the telemetry data
ds1 <- read.csv(file.path("Data", "Bobcats_cleaned_v4_20241021.csv"))
## Load the metadata file
meta1 <- read.csv(file.path("Data", "Bobcats_Metadata_20241021.csv"))
with(meta1, table(Location, CollarType, Useable))
with(meta1, table(Sex, CollarType, Useable))
### 71 GPS collared cats (30 female, 41 male), 70 useable, 1 not useable (1 male)
### 66 VHF collared cats (25 female, 35 male, 6 unknown), 55 useable, 11 not useable (4 female, 5 male, 2 unknown)

times <- difftime(lubridate::mdy_hm(meta1$EndDate), lubridate::mdy_hm(meta1$StartDate), units = "days")
summary(as.numeric(times))

# Step 2: Prep the data for use in CTMM ####
## Identify VHF data to exclude from DOP check
ds1$checkDOP <- ifelse(ds1$CollarType == "VHF" | ds1$Location == "KingRanch", F, T)

## Remove points with high HDOP
ds1$goodDOP <- ifelse(ds1$checkDOP == F, T, ifelse(!is.na(ds1$HDOP) & ds1$HDOP <= 10.0, T, F))
ds2 <- ds1[ds1$goodDOP,]

## Check for and remove duplicates and make sure that the points are in order
#ds2_split <- split(ds2, ds2$ID)
#ds2_split2 <- lapply(ds2_split, function(x){
#  x$dup <- duplicated(x$DateTime)
#  x$badfix <- x$UTM_14N_X == x$UTM_14N_Y
#  x1 <- x[x$dup == FALSE & x$badfix == FALSE,]
#  x2 <- x1[order(x1$DateTime),]
#  return(x2)
#})
#ds2 <- do.call(rbind, ds2_split2)
#row.names(ds2) <- NULL

#write.csv(ds2[,-c(17,18,19,20)], "Bobcats_cleaned_v3_20240122.csv", row.names = FALSE)

## Load the cleaned data
#ds2 <- read.csv(file.path("Data", "Bobcats_cleaned_v3_20240122.csv"))

## Subset only the required information for creating the telemetry object
ds3 <- ds2[,c("ID", "DateTime", "WGS84_Long", "WGS84_Lat", "HDOP")]
colnames(ds3) <- c("individual.local.identifier", "timestamp", "location.long", "location.lat", "GPS.HDOP")


# Step 3: Create the telemetry object for each individual ####
tAll1 <- ctmm::as.telemetry(ds3, projection = "epsg:26914")
all(names(tAll1) == meta1$ID)

## Calculate net squared displacements
NSD <- do.call(rbind, lapply(1:length(tAll1), function(i){
  x <- tAll1[[i]]
  x1 <- amt::make_track(x, .x = longitude, .y = latitude, .t = timestamp, crs = "epsg:4326")
  x1 <- amt::transform_coords(x1, crs_to = "epsg:32614")
  x1$nsd_m2 <- amt::nsd(x1)
  x1$nsd_km2 <- x1$nsd_m2/1000000
  name <- names(tAll1)[i]
  x1$ID <- name
  x1$t2 <- 1:nrow(x1)
  #x2 <- data.frame(ID = name, av_nsd = mean(x1$nsd))
  return(x1)
}))

NSD2 <- merge.data.frame(NSD, meta1[,1:4], by = "ID")

library(ggplot2)

ggplot(data = NSD2, mapping = aes(x = t2, y = nsd_km2, color = ID)) + 
  geom_line() + 
  facet_grid(rows = vars(Location), cols = vars(CollarType), scales = "free") + 
  theme_bw()

nsd.plot <- ggplot(data = NSD2, mapping = aes(x = t2, y = nsd_km2)) + 
  geom_line() + 
  facet_wrap(facets = vars(ID), scales = "free") + 
  theme_bw()
ggsave(filename = paste("nsd_all_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = ""), plot = nsd.plot, device = "tiff", width = 16, height = 13, units = "in", dpi = 300, compression = "lzw")

NSD_mean <- dplyr::summarise(dplyr::group_by(NSD2, Location, ID, Sex, CollarType), av_NSD = mean(nsd_km2))
summary(NSD_mean$av_NSD)
quantile(NSD_mean$av_NSD, probs = c(0.75, 0.80, 0.90, 0.95, 0.99))

meta2 <- merge.data.frame(meta1, NSD_mean, by = "ID")

ggplot(data = NSD_mean, aes(x = av_NSD)) + 
  geom_boxplot()

NSD_summary <- dplyr::summarise(dplyr::group_by(NSD2, ID), 
                                av_NSD = mean(nsd_km2), med_NSD = median(nsd_km2), 
                                q90 = quantile(nsd_km2, 0.90), q95 = quantile(nsd_km2, 0.95), q99 = quantile(nsd_km2, 0.99))
meta2 <- merge.data.frame(meta1, NSD_summary, by = "ID")

## Remove cats that are not usable
tAll2 <- tAll1[meta1$ID[meta1$Useable == "YES"]]

## Separate them by the fix rate
tVHF <- tAll1[which(meta1$AnalysisGroup == "VHF" & meta1$Useable == "YES")]
tGPS <- tAll1[which(meta1$AnalysisGroup == "Constant" & meta1$Useable == "YES")]
tHix <- tAll1[which(meta1$AnalysisGroup == "Hixon" & meta1$Useable == "YES")]
tMoon <- tAll1[which(meta1$AnalysisGroup == "Moon" & meta1$Useable == "YES")]
tUnk <- tAll1[which(meta1$AnalysisGroup == "Unk" & meta1$Useable == "YES")]

## Plot the data
par(mfrow = c(1,1))
plot(tVHF, col = rainbow(length(tVHF)))
plot(tGPS, col = rainbow(length(tGPS)))
plot(tHix, col = rainbow(length(tHix)))
plot(tMoon, col = rainbow(length(tMoon)))
plot(tUnk, col = rainbow(length(tUnk)))

# Step 4: Check for autocorrelation to determine which telemetry objects can be pooled for variogram estimation ####
## The way this is done in ctmm is by using a variogram to check for autocorrelation among locations
## However, when the sampling rate varies, you need to account for this in the variogram calculation

## Creating the variograms for different sampling methods
library(ctmm)
varVHF <- lapply(tVHF, ctmm::variogram)
varGPS <- lapply(tGPS, ctmm::variogram)
varHix <- lapply(tHix, ctmm::variogram, dt = c(1, 2) %#% "hour")
varMoon <- lapply(tMoon, ctmm::variogram, dt = c(12, 0.5) %#% "hour")
varUnk <- lapply(tUnk, ctmm::variogram, dt = c(12, 0.5, 5/60) %#% "hour")

## For calculation of the fitted home range, we need to account for autocorrelation among data points

### Plot the variograms for each cat to determine if any can be pooled to get a more robust estimate of the population-level variability in the animal movements
#### VHF Cats
#par(mfrow = c(5,5))
#lapply(1:length(varVHF), function(i){plot(varVHF[[i]], main = names(varVHF)[i], xlim = c(0,60 %#% "day"), ylim = c(0,4 %#% "km^2"))})
all(names(tVHF) == names(varVHF))
names(varVHF)
## Pool EB13M, LB28M, LB071M, LB146M, LB152M, LB154M, LB163M, WB12M, WB39F
tVHF_1 <- tVHF[names(tVHF) %in% c("EB13M", "LB028M", "LB071M", "LB146M", "LB152M", "LB154M", "LB163M", "WB12M", "WB39F")]
varVHF_1 <- varVHF[names(varVHF) %in% c("EB13M", "LB028M", "LB071M", "LB146M", "LB152M", "LB154M", "LB163M", "WB12M", "WB39F")]
varVHF_1p <- mean(varVHF_1)
## Pool all others
tVHF_2 <- tVHF[!(names(tVHF) %in% c("EB13M", "LB028M", "LB071M", "LB146M", "LB152M", "LB154M", "LB163M", "WB12M", "WB39F"))]
varVHF_2 <- varVHF[!(names(varVHF) %in% c("EB13M", "LB028M", "LB071M", "LB146M", "LB152M", "LB154M", "LB163M", "WB12M", "WB39F"))]
varVHF_2p <- mean(varVHF_2)

#### Hixon Cats
#par(mfrow = c(2,5))
#lapply(1:length(varHix), function(i){plot(varHix[[i]], main = names(varHix)[i], xlim = c(0,300 %#% "day"), ylim = c(0,8 %#% "km^2"))})
all(names(tHix) == names(varHix))
names(varHix)
## HR01F on its own (dispersal event)
tHix_1 <- tHix[["HR01F"]]
varHix_1 <- varHix[["HR01F"]]
## HR03M on its own (high variability)
tHix_2 <- tHix[["HR03M"]]
varHix_2 <- varHix[["HR03M"]]
## HR09F on its own (not working pooled)
tHix_3 <- tHix[["HR09F"]]
varHix_3 <- varHix[["HR09F"]]
## HR10M on its own (exploratory movement)
tHix_4 <- tHix[["HR10M"]]
varHix_4 <- varHix[["HR10M"]]
## Pool all others
tHix_5 <- tHix[!(names(tHix) %in% c("HR01F", "HR03M", "HR09F", "HR10M"))]
varHix_5 <- varHix[!(names(varHix) %in% c("HR01F", "HR03M", "HR09F", "HR10M"))]
varHix_5p <- mean(varHix_5)

#### El Sauz, Yturria, King GPS
#par(mfrow = c(4,5))
#lapply(1:length(varGPS), function(i){plot(varGPS[[i]], main = names(varGPS)[i], xlim = c(0,300 %#% "day"), ylim = c(0,4 %#% "km^2"))})
all(names(tGPS) == names(varGPS))
names(varGPS)
## EB39M on its own (exploratory movement)
tGPS_1 <- tGPS[["EB39M"]]
varGPS_1 <- varGPS[["EB39M"]]
## Pool EB22F, EB31M, EB32F, EB41F, EB67M, EB72M, EB73F, EB76M (high variability)
tGPS_2 <- tGPS[names(tGPS) %in% c("EB22F", "EB31M", "EB32F", "EB41F", "EB67M", "EB72M", "EB73F", "EB76M")]
varGPS_2 <- varGPS[names(varGPS) %in% c("EB22F", "EB31M", "EB32F", "EB41F", "EB67M", "EB72M", "EB73F", "EB76M")]
varGPS_2p <- mean(varGPS_2)
## King Ranch cats on their own (no HDOP information)
tGPS_3 <- tGPS[names(tGPS) %in% c("KRB1M", "KRB2M", "KRB3F", "KRB4F", "KRB5F", "KRB6M", "KRB8M", "KRB9M")]
varGPS_3 <- varGPS[names(varGPS) %in% c("KRB1M", "KRB2M", "KRB3F", "KRB4F", "KRB5F", "KRB6M", "KRB8M", "KRB9M")]
varGPS_3p <- mean(varGPS_3)
## HB729232F on its own (dispersal)
tGPS_4 <- tGPS[["HB729232F"]]
varGPS_4 <- varGPS[["HB729232F"]]
## Pool all others
tGPS_5 <- tGPS[!(names(tGPS) %in% c("EB39M", "EB22F", "EB31M", "EB32F", "EB41F", "EB67M", "EB72M", "EB73F", "EB76M", "KRB1M", "KRB2M", "KRB3F", "KRB4F", "KRB5F", "KRB6M", "KRB8M", "KRB9M", "HB729232F"))]
varGPS_5 <- varGPS[!(names(varGPS) %in% c("EB39M", "EB22F", "EB31M", "EB32F", "EB41F", "EB67M", "EB72M", "EB73F", "EB76M", "KRB1M", "KRB2M", "KRB3F", "KRB4F", "KRB5F", "KRB6M", "KRB8M", "KRB9M", "HB729232F"))]
varGPS_5p <- mean(varGPS_5)
## Try EB67M by itself
tGPS_6 <- tGPS[["EB67M"]]
varGPS_6 <- varGPS[["EB67M"]]
## Try EB83F by itself
tGPS_7 <- tGPS[["EB83F"]]
varGPS_7 <- varGPS[["EB83F"]]

#### El Sauz Moon
#par(mfrow = c(1,3))
#lapply(1:length(varMoon), function(i){plot(varMoon[[i]], main = names(varMoon)[i], xlim = c(0, 120 %#% "day"), ylim = c(0,12 %#% "km^2"))})
all(names(tMoon) == names(varMoon))
names(varMoon)
## EB16M on its own
tMoon_1 <- tMoon[["EB16M"]]
varMoon_1 <- varMoon[["EB16M"]]
## Pool all others
tMoon_2 <- tMoon[!(names(tMoon) %in% c("EB16M"))]
varMoon_2 <- varMoon[!(names(varMoon) %in% c("EB16M"))]
varMoon_2p <- mean(varMoon_2)

#### El Sauz WTF
#par(mfrow = c(1,1))
#lapply(1:length(varUnk), function(i){plot(varUnk[[i]], main = names(varUnk)[i])})
## There is only 1 cat with this schedule
tUnk_1 <- tUnk[["EB18M"]]
varUnk_1 <- varUnk[["EB18M"]]


# Step 5: Calculate the pooled variograms ####
## Calculate Pooled Variograms and estimates of sigma and tau (REPLACE ALL XX with actual numbers)
par(mfrow = c(1,1))
### VHF cats
#### Assume IID anisotropic (only estimate sigma, not tau)
ctmm::variogram.fit(varVHF_1p)
ctmmVHF_1 <- ctmm::ctmm(sigma = 2.62951 %#% "km^2")
ctmm::variogram.fit(varVHF_2p)
ctmmVHF_2 <- ctmm::ctmm(sigma = 40.66743 %#% "hm^2")

### Hixon cats
#### Assume OUF anisotropic (estimate sigma and tau)
ctmm::variogram.fit(varHix_1)
ctmmHix_1 <- ctmm::ctmm(sigma = 193.31112 %#% "km^2", tau = c(6.60762 %#% "mon", 61.64582 %#% "min"), error = T)
ctmm::variogram.fit(varHix_2)
ctmmHix_2 <- ctmm::ctmm(sigma = 0.77091 %#% "km^2", tau = c(0.93033 %#% "day", 48.2501 %#% "min"), error = T)
ctmm::variogram.fit(varHix_3)
ctmmHix_3 <- ctmm::ctmm(sigma = 11.27514 %#% "km^2", tau = c(10.92535 %#% "day", 131.70327 %#% "min"), error = T)
ctmm::variogram.fit(varHix_4)
ctmmHix_4 <- ctmm::ctmm(sigma = 29.18778 %#% "hm^2", tau = c(0.54726 %#% "day", 65.66866 %#% "min"), error = T)
ctmm::variogram.fit(varHix_5p)
ctmmHix_5 <- ctmm::ctmm(sigma = 43.34432 %#% "hm^2", tau = c(13.31281 %#% "hr", 63.86975 %#% "min"), error = T)

### El Sauz, Yturria, Hacienda, Laguna, King GPS
#### Assume OUF anisotropic (estimate sigma and tau)
ctmm::variogram.fit(varGPS_1)
ctmmGPS_1 <- ctmm::ctmm(sigma = 15.32283 %#% "km^2", tau = c(36.01369 %#% "hour", 4.64151 %#% "hour"), error = T)
ctmm::variogram.fit(varGPS_2p)
ctmmGPS_2 <- ctmm::ctmm(sigma = 2.75563 %#% "km^2", tau = c(4.72170 %#% "day", 22.88803 %#% "min"), error = T)
ctmm::variogram.fit(varGPS_3p)
ctmmGPS_3 <- ctmm::ctmm(sigma = 1.00244 %#% "km^2", tau = c(34.10727 %#% "hr", 38.30668 %#% "min"), error = F)
ctmm::variogram.fit(varGPS_4)
ctmmGPS_4 <- ctmm::ctmm(sigma = 17.20793 %#% "km^2", tau = c(22.05273 %#% "day", 149.33779 %#% "min"), error = T)
ctmm::variogram.fit(varGPS_5p)
ctmmGPS_5 <- ctmm::ctmm(sigma = 26.09349 %#% "hm^2", tau = c(12.32579 %#% "hr", 19.82107 %#% "min"), error = T)
ctmm::variogram.fit(varGPS_6)
ctmmGPS_6 <- ctmm::ctmm(sigma = 4.76586 %#% "km^2", tau = c(6.52887 %#% "day", 126.41529 %#% "min"), error = T)
ctmm::variogram.fit(varGPS_7)
ctmmGPS_7 <- ctmm::ctmm(sigma = 2.05901 %#% "km^2", tau = c(9.70996 %#% "day", 7.21100 %#% "min"), error = T)

### El Sauz Moon
#### Assume OUF anisotropic (estimate sigma and tau)
ctmm::variogram.fit(varMoon_1)
ctmmMoon_1 <- ctmm::ctmm(sigma = 5.68432 %#% "km^2", tau = c(8.54387 %#% "day", 58.82682 %#% "min"), error = T)
ctmm::variogram.fit(varMoon_2p)
ctmmMoon_2 <- ctmm::ctmm(sigma = 25.45342 %#% "hm^2", tau = c(8.24854 %#% "hr", 41.04689 %#% "min"), error = T)

### El Sauz WTF
#### Assume OUF anisotropic (estimate sigma and tau)
ctmm::variogram.fit(varUnk_1)
ctmmUnk_1 <- ctmm::ctmm(sigma = 44.11173 %#% "hm^2", tau = c(8.78027 %#% "hr", 8.06504 %#% "min"), error = T)

# Step 6: Estimating the error distributions for each subset of cats ####

## Calculate ctmm.select objects to determine the best variogram structure for our data
### ctmm.select() starts with the most complex model and runs subsequently simpler models and compares them with AICc
## Calcualte maximum number of cores that can be used for this (leave 2 available for other use)
#parallel::detectCores()

## Running ctmm.select on each group of cats
### Create a function for running these ctmm.select() objects using parallel processing
ctmmSelectFun <- function(x, pp, cores.left){
  #x <- selHix_start
  #pp <- TRUE
  #cores.left <- 2
  
  print(paste("This function started at ", Sys.time(), sep = ""))
  arg1 <- list(verbose = T, trace = 2, cores = 1)
  
  if (isTRUE(pp)) {
    message("Parallel processing enabled")
    if (is.null(cores.left)) {
      cores.left <- 2
    }
    else {
      cores.left <- tryCatch(as.numeric(cores.left), error = function(e) {
        message("There was an error coercing cores.left to a number. The default of 2 cores are not utilized")
        return(2)
      }, warning = function(w) {
        message("Could not coerce cores.left to a number. The default of 2 cores are not utilized")
        return(2)
      })
    }
    cl1 <- parallel::makeCluster(parallel::detectCores() - cores.left, outfile = "out.txt")
    parallel::clusterExport(cl1, varlist = c("x", "arg1"), envir = environment())
  }
  else {
    cl1 <- NULL
  }
  
  out <- pbapply::pblapply(1:length(x), cl = cl1, function(i){
    print(paste("Started ", names(x)[i], " at ", Sys.time(), sep = ""))
    x1 <- c(x[[i]], arg1)
    x2 <- do.call(ctmm::ctmm.select, x1)
    saveRDS(x2, file = paste("ctmmSelect_", names(x)[i], ".RDS", sep = ""))
    print(paste("Finished with ", names(x)[i], " at ", Sys.time(), sep = ""))
    return(x2)
  })
  names(out) <- names(x)
  
  if (isTRUE(pp)) {
    parallel::stopCluster(cl1)
  }
  
  print(paste("This function finished at ", Sys.time(), sep = ""))
  return(out)
  rm(arg1, cl1, out)
}

### If you have to run each type of error separately
ctmmfitFun <- function(x, ID, sigma, tau1, tau2, selections = c("OUF_aniso", "OUF_iso", "OU_aniso", "OU_iso"), pp, cores.left){
  #x = selGPS_start
  #ID = "EB68M"
  #sigma = 26.09349 %#% "hm^2"
  #tau1 = 12.32579 %#% "hr"
  #tau2 = 19.82107 %#% "min"
  #selections = c("OUF_aniso", "OUf_iso", "OU_aniso")
  #pp = FALSE
  #cores.left = 17
  
  print(paste("This function started at ", Sys.time(), sep = ""))
  
  # Extract the correct individual
  ds <- x[[ID]]
  
  # Define the different error structures
  start1 <- list("OUF_aniso"  = list(ds$data, ctmm::ctmm(sigma = sigma, tau = c(tau1, tau2), isotropic = F, error = T)), 
                 "OUF_iso"    = list(ds$data, ctmm::ctmm(sigma = sigma, tau = c(tau1, tau2), isotropic = T, error = T)), 
                 "OU_aniso"   = list(ds$data, ctmm::ctmm(sigma = sigma, tau = c(tau1), isotropic = F, error = T)), 
                 "OU_iso"     = list(ds$data, ctmm::ctmm(sigma = sigma, tau = c(tau1), isotropic = T, error = T)), 
                 "OUf__aniso" = list(ds$data, ctmm::ctmm(sigma = sigma, tau = c(tau1, tau1), isotropic = F, error = T)), 
                 "OUf__iso"   = list(ds$data, ctmm::ctmm(sigma = sigma, tau = c(tau1, tau1), isotropic = T, error = T)))
  start2 <- start1[selections]
  
  # Enable parallel processing
  if (isTRUE(pp)) {
    message("Parallel processing enabled")
    if (is.null(cores.left)) {
      cores.left <- 2
    }
    else {
      cores.left <- tryCatch(as.numeric(cores.left), error = function(e) {
        message("There was an error coercing cores.left to a number. The default of 2 cores are not utilized")
        return(2)
      }, warning = function(w) {
        message("Could not coerce cores.left to a number. The default of 2 cores are not utilized")
        return(2)
      })
    }
    cl1 <- parallel::makeCluster(parallel::detectCores() - cores.left, outfile = "out_fit.txt")
    parallel::clusterExport(cl1, varlist = c("start2"), envir = environment())
  }else {
    cl1 <- NULL
  }
  
  # Run each item
  out <- pbapply::pblapply(1:length(start2), cl = cl1, function(i){
    x1 <- start2[[i]]
    name <- names(start2)[i]
    names(x1) <- c("data", "CTMM")
    x1 <- c(x1, list(trace = 1))
    x2 <- do.call(ctmm::ctmm.fit, x1)
    saveRDS(x2, paste(ID, "_", name, "_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
    print(paste("Finished with ", name, " at ", Sys.time(), sep = ""))
    return(x2)
    #rm(x1, name, x2)
  })
  
  if (isTRUE(pp)) {
    parallel::stopCluster(cl1)
  }
  
  names.ds <- data.frame(short = c("OUF_aniso", "OUF_iso", 
                                   "OU_aniso", "OU_iso", 
                                   "OUf__aniso", "OUf__iso"), 
                         full = c("OUF anisotropic error", "OUF error", 
                                  "OU anisotropic error", "OU error", 
                                  "OUf anisotropic error", "OUf error"))
  
  names(out) <- names.ds$full[names.ds$short %in% names(out)]
  
  out1 <- out[rownames(summary(out))]
  
  print(paste("This function finished at ", Sys.time(), sep = ""))
  return(out1)
}
names.ds <- data.frame(short = c("OUF_aniso", "OUF_iso", "OU_aniso", "OU_iso", "OUf__aniso", "OUf__iso"), full = c("OUF anisotropic error", "OUF error", "OU anisotropic error", "OU error", "OUf anisotropic error", "OUf error"))

### VHF cats (COMPLETE)
#### Set up the data for running in the function
selVHF_start <- c(lapply(tVHF_1, function(x){list(data = x, CTMM = ctmmVHF_1)}), 
                  lapply(tVHF_2, function(x){list(data = x, CTMM = ctmmVHF_2)}))
selVHF_start <- selVHF_start[names(tVHF)]
#### Run the function on all the cats in the group
#selVHF_end <- ctmmSelectFun(selVHF_start, pp = T, cores.left = 2)  # NEED TO UPDATE
#### Resort the cats into their original order
#selVHF <- selVHF_end[names(tVHF)]
#### Save the output as an RDS file for easy loading later
#saveRDS(selVHF, file = paste("ctmm.select_VHF_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
#### Load the completed output file
selVHF <- readRDS(file.path("Data", "ctmm.select_VHF_20240122.RDS"))   # This one is correct
#### Select only the best model for each cat
selVHF_best <- lapply(selVHF, function(x){x[[1]]})

### Hixon cats (COMPLETE)
#### Same steps as above
selHix_start <- c(list(HR01F = list(data = tHix_1, CTMM = ctmmHix_1)), 
                  list(HR03M = list(data = tHix_2, CTMM = ctmmHix_2)), 
                  list(HR10M = list(data = tHix_3, CTMM = ctmmHix_3)), 
                  list(HR09F = list(data = tHix_4, CTMM = ctmmHix_4)),
                  lapply(tHix_5, function(x){list(data = x, CTMM = ctmmHix_5)}))
selHix_start <- selHix_start[names(tHix)]
#selHix_end <- ctmmSelectFun(selHix_start, pp = T, cores.left = 4)
#### HR09F models are getting hung up on the OUf anisotropic model so need to do manual selection of all other models
#HR09F_end1 <- ctmmfitFun(x = selHix_start, ID = "HR09F", sigma = 29.18778 %#% "hm^2", tau1 = 0.54726 %#% "day", tau2 = 64.66866 %#% "min", selections = names.ds[c(1:4,6),1], pp = TRUE, cores.left = 6)
#summary(HR09F_end1)
#HR09F_end <- HR09F_end1[rownames(summary(HR09F_end1))]
#temp_Hix1 <- list.files(path = "Individual_Results_20241022", pattern = "ctmmSelect")
#temp_Hix2 <- temp_Hix1[sub("ctmmSelect_", "", fs::path_ext_remove(temp_Hix1)) %in% names(selHix_start)]
#selHix_end2 <- lapply(file.path("Individual_Results_20241022", temp_Hix2), readRDS)
#names(selHix_end2) <- sub("ctmmSelect_", "", fs::path_ext_remove(temp_Hix2))
#rm(temp_Hix1, temp_Hix2)
#### Bring it all back together
#selHix_end <- c(list(HR09F = HR09F_end2), selHixend2)
#selHix <- selHix_end[names(tHixu)]
#saveRDS(selHix, file = paste("ctmm.select_Hix_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
selHix <- readRDS(file.path("Data", "ctmm.select_Hix_20230619.RDS"))   # Adjust these as they finish
selHix_best <- lapply(selHix, function(x){x[[1]]})

### GPS cats
selGPS_start <- c(list(EB39M = list(data = tGPS_1, CTMM = ctmmGPS_1)), 
                  lapply(tGPS_2, function(x){list(data = x, CTMM = ctmmGPS_2)}), 
                  lapply(tGPS_3, function(x){list(data = x, CTMM = ctmmGPS_3)}), 
                  list(HB729232F = list(data = tGPS_4, CTMM = ctmmGPS_4)), 
                  lapply(tGPS_5, function(x){list(data = x, CTMM = ctmmGPS_5)}), 
                  list(EB67M = list(data = tGPS_6, CTMM = ctmmGPS_6)), 
                  list(EB83F = list(data = tGPS_7, CTMM = ctmmGPS_7)))
selGPS_start <- selGPS_start[names(tGPS)]
#selGPS_end <- ctmmSelectFun(selGPS_start[57:58], pp = T, cores.left = 16)
#### Because we had to rerun some, but not all, we need to only run those that we need to
#selGPS_end1 <- ctmmSelectFun(selGPS_start["EB40F"], pp = F, cores.left = 4)
#temp_GPS1 <- list.files(path = "Individual_Results_20241022", pattern = "ctmmSelect")
#temp_GPS2 <- temp_GPS1[sub("ctmmSelect_", "", fs::path_ext_remove(temp_GPS1)) %in% names(selGPS_start)]
#selGPS_end2 <- lapply(file.path("Individual_Results_20241022", temp_GPS2), readRDS)
#names(selGPS_end2) <- sub("ctmmSelect_", "", fs::path_ext_remove(temp_GPS2))
#rm(temp_GPS1, temp_GPS2)
#### EB68M are having issues running so manual selection of all other models
##### In addition, some of the fit functions are not working properly
##### Runs the fit functions that are still outstanding
#EB68M_end1 <- ctmmfitFun(x = selGPS_start, ID = "EB68M", sigma = 26.09349 %#% "hm^2", tau1 = 12.32579 %#% "hr", tau2 = 19.82107 %#% "min", selections = names.ds[c(1,3),1], pp = TRUE, cores.left = 18)
##### Load the fit functions that already successfully ran
#temp_EB68M <- list.files(path = "Individual_Results_20241022", pattern = "EB68M")
#EB68M_end2 <- lapply(file.path("Individual_Results_20241022", temp_EB68M), readRDS)
#names(EB68M_end2) <- names.ds$full[names.ds$short %in% sub("EB68M_", "", sub("_20241022.RDS", "", temp_EB68M))]
##### Combine the new and previously run functions
#EB68M_end <- c(EB68M_end1, EB68M_end2)
##### Do model selection to identify the best model
#summary(EB68M_end)
##### Sort the models so the "best" model is first
#EB68M_end <- EB68M_end[rownames(summary(EB68M_end))]
#saveRDS(EB68M_end, file = "ctmmSelect_EB68M.RDS")
#EB68M_end <- readRDS("ctmmSelect_EB68M.RDS")
#rm(temp_EB68M)
#EB67M_end1 <- ctmmfitFun(x = selGPS_start, ID = "EB67M", sigma = 4.76586 %#% "km^2", tau1 = 6.52887 %#% "day", tau2 = 126.41529 %#% "min", selections = names.ds[,1], pp = TRUE, cores.left = 14)
#temp_EB67M <- list.files(path = "Individual_Results_20241022", pattern = "EB67M_O")
#EB67M_end1 <- lapply(file.path("Individual_Results_20241022", temp_EB67M), readRDS)
#names(EB67M_end1) <- names.ds$full[names.ds$short %in% sub("EB67M_", "", sub("_20241107.RDS", "", temp_EB67M))]
#EB67M_end2 <- EB67M_end1[rownames(summary(EB67M_end1))]
#saveRDS(EB67M_end2, file = "ctmmSelect_EB67M_ind.RDS")
#EB83F_end1 <- ctmmfitFun(x = selGPS_start, ID = "EB83F", sigma = 2.05901 %#% "km^2", tau1 = 9.70996 %#% "day", tau2 = 7.21100 %#% "min", selections = names.ds[,1], pp = TRUE, cores.left = 14)
#temp_EB83F <- list.files(path = "Individual_Results_20241022", pattern = "EB83F_O")
#EB83F_end1 <- lapply(file.path("Individual_Results_20241022", temp_EB83F), readRDS)
#names(EB83F_end1) <- names.ds$full[names.ds$short %in% sub("EB83F_", "", sub("_20241107.RDS", "", temp_EB83F))]
#EB83F_end2 <- EB83F_end1[rownames(summary(EB83F_end1))]
#saveRDS(EB83F_end2, file = "ctmmSelect_EB83F_ind.RDS")
#### Now, bring it all back together
#selGPS_end <- c(selGPS_end1, selGPS_end2)
#selGPS <- selGPS_end[names(tGPS)]
#saveRDS(selGPS, file = paste("ctmm.select_GPS_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
selGPS <- readRDS(file.path("Data", "ctmm.select_GPS_20241108.RDS"))   # This is complete
selGPS_best <- lapply(selGPS, function(x){x[[1]]})

### Moon cats (Complete)
selMoon_start <- c(list(EB16M = list(data = tMoon_1, CTMM = ctmmMoon_1)), 
                   lapply(tMoon_2, function(x){list(data = x, CTMM = ctmmMoon_2)}))
selMoon_start <- selMoon_start[names(tMoon)]
#selMoon_end <- ctmmSelectFun(selMoon_start, pp = T, cores.left = 2)
#selMoon <- selMoon_end[names(tMoonu)]
#saveRDS(selMoon, file = paste("ctmm.select_Moon_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
selMoon <- readRDS(file.path("Data", "ctmm.select_Moon_20230616.RDS")) # This is complete
selMoon_best <- lapply(selMoon, function(x){x[[1]]})

### WTF cats (Complete)
selUnk_start <- list(EB18M = list(data = tUnk_1, CTMM = ctmmUnk_1))
#selUnk_end <- ctmmSelectFun(selUnk_start, pp = F, cores.left = 2)
#selUnk <- selUnk_end[names(tUnku)]
#saveRDS(selUnk, file = paste("ctmm.select_Unk_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
selUnk <- readRDS(file.path("Data", "ctmm.select_Unk_20230616.RDS"))   # This is complete
selUnk_best <- lapply(selUnk, function(x){x[[1]]})


# Step 7: Re-group cats based on time collar active ####
## Recombine cats into a single group for dividing up data 
selAll <- c(selVHF_best, selHix_best, selGPS_best, selMoon_best, selUnk_best)
selAll <- selAll[names(tAll2)]
names(tAll2) == names(selAll)
all(names(tAll2) == names(selAll))

## Group cats by their overlap group for future analyses (based on location and time collar was active)
### Write a function to reduce the number of outputs we are dealing with here
overlapSelect <- function(group){
  if(group == "all"){
    x1 <- meta1[meta1$Useable == "YES",]
    name = "all"
  }else if(group == "resident"){
    x1 <- meta1[meta1$Useable == "YES" & !(meta1$Note %in% c("Dispersal", "Exploratory movements")),]
    name = "res"
  }else if(group == "res_exp"){
    x1 <- meta1[meta1$Useable == "YES" & !(meta1$Note %in% c("Dispersal")),]
    name = "resexp"
  }else if(group == "large"){
    x1 <- meta1[meta1$Useable == "YES" & (meta1$Note %in% c("Dispersal", "Exploratory movements")),]
    name = "big"
  }else if(group == "explore"){
    x1 <- meta1[meta1$Useable == "YES" & (meta1$Note %in% c("Exploratory movements")),]
    name = "exp"
  }else if(group == "disperse"){
    x1 <- meta1[meta1$Useable == "YES" & (meta1$Note %in% c("Dispersal")),]
    name = "dis"
  }else{
    stop("You chose an invalid group. Try again with either c('all', 'resident', 'large', 'explore', 'disperses').")
  }
  x2 <- lapply(letters, function(a){
    a1 <- x1[grep(pattern = a, x1$OverlapGroup),]
    a2 <- split(a1$ID, f = a1$Location)
  })
  names(x2) <- letters
  x3 <- do.call(c, x2[which(sapply(x2, length) > 0)])
  
  tel <- lapply(x3, function(x){tAll2[x]})
  sel <- lapply(x3, function(x){selAll[x]})
  
  return(list(tel = tel, sel = sel, name = name))
  rm(x1, x2, x3, name, tel, sel)
  #rm(group)
}

### Include all cats, including those that dispersed or had exploratory movements
groups1 <- c(overlapSelect(group = "all"), grid = NA, pp = T, cores.left = 4)

### Remove cats that had dispersals or exploratory movements
groups2 <- c(overlapSelect(group = "resident"), grid = NA, pp = T, cores.left = 4)

### Only cats that had dispersals or exploratory movements
groups3 <- c(overlapSelect(group = "large"), grid = NA, pp = T, cores.left = 4)

### Remove dispersers (because they can cause issues with the grid calculation)
groups3a <- c(overlapSelect(group = "res_exp"), grid = NA, pp = T, cores.left = 4)

### Only dispersers (because they can cause issues with the grid calculation)
groups3b <- c(overlapSelect(group = "disperse"), grid = NA, pp = T, cores.left = 4)

# Step 8: Calculate the kdes and degrees of overlap ####
## Function for calculating all the kdes and overlaps in parallel
kdeFun <- function(tel, sel, grid=NULL, name, pp, cores.left){
  # grid should either be NULL or a list of grid objects for each item in the tel and sel objects
  
  if(is.na(grid)){
    grid <- NULL
  }
  
  if (isTRUE(pp)) {
    message("Parallel processing enabled")
    if (is.null(cores.left)) {
      cores.left <- 2
    }
    else {
      cores.left <- tryCatch(as.numeric(cores.left), error = function(e) {
        message("There was an error coercing cores.left to a number. The default of 2 cores are not utilized")
        return(2)
      }, warning = function(w) {
        message("Could not coerce cores.left to a number. The default of 2 cores are not utilized")
        return(2)
      })
    }
    cl1 <- parallel::makeCluster(parallel::detectCores() - cores.left, outfile = "out.txt")
    parallel::clusterExport(cl1, varlist = c("tel", "sel", "grid", "name"), envir = environment())
  }else {
    cl1 <- NULL
  }
  
  out1 <- pbapply::pblapply(1:length(tel), cl = cl1, function(i){
    print(paste("Starting group ", names(tel)[i], " at ", Sys.time(), sep = ""))
    x <- ctmm::akde(data = tel[[i]], CTMM = sel[[i]], grid = grid[[i]])
    print(paste("KDE calculated for ", names(tel)[i], ". Starting overlap calculation...", sep = ""))
    o <- ctmm::overlap(x, method = "Bhattacharyya", level = 0.95)
    out <- list(kde = x, overlap = o)
    saveRDS(out, file = paste("kdeoverlap_", name, "_", names(tel)[i], "_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
    print(paste("Finished group ", names(tel)[i], " at ", Sys.time(), sep = ""))
    return(out)
  })
  
  kdes <- lapply(out1, function(x){x$kde})
  overs <- lapply(out1, function(x){x$overlap})
  names(kdes) <- names(tel)
  names(overs) <- names(tel)
  
  out2 <- list(kde = kdes, overlap = overs)
  return(out2)
}

## Calculate KDEs and degrees of overlap
### Attempt all cats (doesn't work due to dispersers being weird)
#kdeoverAll <- do.call(kdeFun, groups1)
#saveRDS(kdeoverAll, file = paste("AKDE_overlap_All_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
#kdeoverAll <- readRDS(file.path("Data", "AKDE_overlap_All_20240124.RDS"))

### Only resident cats and those with exploratory movements
#kdeoverRes <- do.call(kdeFun, groups3a)
#saveRDS(kdeoverRes, file = paste("AKDE_overlap_resexp_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
#kdeoverRes <- readRDS(file.path("Individual_Results", "AKDE_overlap_resexp_20240124.RDS"))

### Only dispersing or exploring cats
#### Need to do each step separately, starting with calculation of KDE
#kdeDis1 <- ctmm::akde(data = groups3b$tel$a.HixonRanch, CTMM = groups3b$sel$a.HixonRanch, grid = kdeoverRes$kde$a.HixonRanch[[1]])
#kdeDis2 <- ctmm::akde(data = groups3b$tel$a.Hacienda, CTMM = groups3b$sel$a.Hacienda, grid = kdeoverRes$kde$a.Hacienda[[1]])
#### Now combine kdes to calculate overlap
#overs1 <- c(("HR01F" = kdeDis1), kdeoverRes$kde$a.HixonRanch)
#overs1 <- overs1[names(groups1$tel$a.HixonRanch)]
#overs2 <- c(("HB729232F" = kdeDis2), kdeoverRes$kde$a.Hacienda)
#overs2 <- overs2[names(groups1$tel$a.Hacienda)]
#### Calculate overlap for each
#overDis1 <- ctmm::overlap(overs1, method = "Bhattacharyya")
#overDis2 <- ctmm::overlap(overs2, method = "Bhattacharyya")
#### Save the individual results
#saveRDS(list(kde = overs1, overlap = overDis1), file = paste("kdeoverlap_dis_a.HixonRanch_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
#saveRDS(list(kde = overs2, overlap = overDis2), file = paste("kdeoverlap_dis_a.Hacienda_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
#### Combine into same structure as others
#temp_overDis1 <- readRDS(file.path("Individual_Results_20241022", "kdeoverlap_dis_a.HixonRanch_20240124.RDS"))
#overs1 <- temp_overDis1$kde
#overDis1 <- temp_overDis1$overlap
#kdeoverDis <- list(kde = list(a.Hacienda = overs2, a.HixonRanch = overs1), overlap = list(a.Hacienda = overDis2, a.HixonRanch = overDis1))
#rm(temp_overDis1)
#saveRDS(kdeoverDis, file = paste("AKDE_overlap_dis_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
#kdeoverDis <- readRDS(file.path("Individual_Results_20241022", "AKDE_overlap_dis_20240124.RDS"))

### Combine resident and dispersed cats
#kdeoverAll <- kdeoverRes
#kdeoverAll$kde[which(names(kdeoverAll$kde) %in% names(kdeoverDis$kde))] <- kdeoverDis$kde
#kdeoverAll$overlap[which(names(kdeoverAll$overlap) %in% names(kdeoverDis$overlap))] <- kdeoverDis$overlap
#saveRDS(kdeoverAll, file = paste("AKDE_overlap_All_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
kdeoverAll <- readRDS(file.path("Data", "AKDE_overlap_All_20241108.RDS"))

### Get only kde information
kdeAll <- kdeoverAll$kde
### Get only overlap information
overAll <- kdeoverAll$overlap

# Step 9: Estimating home range size for populations and individuals ####
## View the estimated home range size for each cat
indHR1 <- do.call(c, kdeAll)
names(indHR1) <- do.call(rbind, strsplit(names(indHR1), split = "[.]"))[,3]
indHR1 <- indHR1[meta1$ID[meta1$ID %in% names(indHR1)]]

## Save the home ranges as rasters
pbapply::pblapply(1:length(indHR1), function(i){
  name <- names(indHR1)[i]
  #ctmm::writeRaster(x = indHR1[[i]], filename = file.path("HomeRangeRasters", paste("HR_rast_", name, "_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), format = "GTiff", overwrite = TRUE)
  #ctmm::writeVector(x = indHR1[[i]], filename = file.path("HomeRangeVectors", paste("HR_95per_", name, "_", format(Sys.Date(), "%Y%m%d"), sep = "")), filetype = "ESRI Shapefile", level.UD = 0.95, overwrite = TRUE)
  #x1 <- sf::st_as_sf(ctmm::SpatialPolygonsDataFrame.UD(indHR1[[i]], convex = FALSE, level.UD = 0.95))
  #sf::st_write(obj = x1, dsn = file.path("HomeRangeVectors", paste("HR_95per_", name, "_", format(Sys.Date(), "%Y%m%d"), sep = "")), driver = "ESRI Shapefile")
})


## Convert the ctmm UDs to a data frame
indHR2 <- data.frame(do.call(rbind, lapply(1:length(indHR1), function(i){
  x1 <- data.frame(ID = names(indHR1)[i], summary(indHR1[[i]])$CI)
})))

## Plot the home ranges
### Do not attempt to plot everything at once...
#lapply(kdeAll, plot)

plot(groups1$tel$i.ElSauz$EB31M, kdeAll$i.ElSauz$EB31M)

plot(groups1$tel$a.Hacienda$HB729232F, kdeAll$a.Hacienda$HB729232F)
plot(groups1$tel$n.ElSauz$EB81F, kdeAll$n.ElSauz$EB81F)
#plot(groups1$tel$a.Haceinda, kdeAll$a.Hacienda, col = c("red", "orange", "green", "blue", "purple"))

#plot(kdeAll$a.HixonRanch)
#plot(groups1$a.HixonRanch, kdeAll$a.HixonRanch, col = rainbow(10))

# All home range analyses
inHR_all <- list(resident = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES"]], 
                 resident2 = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 females = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Sex == "F"]],
                 females2 = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 males = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Sex == "M"]],
                 males2 = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 all_1980s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1980s" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 sex_1980s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1980s" & meta1$Sex != "U" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 females_1980s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1980s" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 males_1980s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1980s" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 all_1990s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1990s" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 females_1990s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1990s" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 males_1990s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1990s" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 all_2000s = inHR_2000s <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2000s" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 females_2000s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2000s" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 males_2000s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2000s" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 all_2010s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2010s" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 females_2010s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2010s" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 males_2010s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2010s" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 all_2020s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2020s" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 females_2020s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2020s" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 males_2020s = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2020s" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 VHF = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="VHF" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 sex_VHF = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="VHF" & meta1$Sex != "U" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 females_VHF = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="VHF" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 males_VHF = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="VHF" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 all_GPS = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="GPS" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 females_GPS = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="GPS" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]],
                 males_GPS = indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="GPS" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M", "EB67M", "EB83F"))]])

inHR_out <- pbapply::pblapply(1:length(inHR_all), function(i){ctmm::meta(inHR_all[[i]], main = names(inHR_all)[i])})

homerange <- data.frame(dataset = c("residents", "residents2", "females", "females2", "males", "males2", "1980s", "1980s", "1980s", "1980s", "1990s", "1990s", "1990s", "2000s", "2000s", "2000s", "2010s", "2010s", "2010s", "2020s", "2020s", "2020s", "VHF", "VHF", "VHF", "VHF", "GPS", "GPS", "GPS"), 
                        sex = c("all", "all", "females", "females", "males", "males", "all", "sex", "females", "males", "all", "females", "males", "all", "females", "males", "all", "females", "males", "all", "females", "males", "all", "sex", "females", "males", "all", "females", "males"), 
                        sample_size = sapply(inHR_all, length), 
                        do.call(rbind, lapply(inHR_out, function(x){x[1,]})), 
                        units = "km^2", 
                        note = c(# Resident and resident2
                                 "", "WB39F, LB071M, EB67M, EB83F removed", 
                                 # Females
                                 "", "WB39F and EB83F removed", 
                                 # Males
                                 "", "LB071M and EB67M removed", 
                                 # 1980s
                                 "LB071M removed", "Cats of unknown sex removed", "", "LB071M removed", 
                                 # 1990s
                                 "WB39F removed", "WB39F removed", "", 
                                 # 2000s
                                 "", "", "", 
                                 # 2010s
                                 "", "", "", 
                                 # 2020s
                                 "EB67M and EB83F removed", "EB83F removed", "EB67M removed", 
                                 # VHF
                                 "WB39F and LB071M removed", "Cats of unknown sex removed", "WB39F removed", "LB071M removed", 
                                 # GPS
                                 "EB67M and EB83F removed", "EB83F removed", "EB67M removed"))                         

#openxlsx::write.xlsx(homerange, file = paste("HomeRanges_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = ""))

names(inHR_out) <- names(inHR_all)

# Use ratios to determine difference in home range estimates
## CI that overlaps 1 = same
HR_comps <- list("MF" = c("females2", "males2"), 
                 "all_collar" = c("sex_VHF", "all_GPS"),
                 "F_collar" = c("females_VHF", "females_GPS"), 
                 "M_collar" = c("males_VHF", "males_GPS"), 
                 "A80_90" = c("sex_1980s", "all_1990s"), 
                 "A80_00" = c("sex_1980s", "all_2000s"), 
                 "A80_10" = c("sex_1980s", "all_2010s"), 
                 "A80_20" = c("sex_1980s", "all_2020s"), 
                 "A90_00" = c("all_1990s", "all_2000s"), 
                 "A90_10" = c("all_1990s", "all_2010s"), 
                 "A90_20" = c("all_1990s", "all_2020s"), 
                 "A00_10" = c("all_2000s", "all_2010s"),
                 "A00_20" = c("all_2000s", "all_2020s"), 
                 "A10_20" = c("all_2010s", "all_2020s"), 
                 "F80_90" = c("females_1980s", "females_1990s"), 
                 "F80_00" = c("females_1980s", "females_2000s"), 
                 "F80_10" = c("females_1980s", "females_2010s"), 
                 "F80_20" = c("females_1980s", "females_2020s"), 
                 "F90_00" = c("females_1990s", "females_2000s"), 
                 "F90_10" = c("females_1990s", "females_2010s"), 
                 "F90_20" = c("females_1990s", "females_2020s"), 
                 "F00_10" = c("females_2000s", "females_2010s"),
                 "F00_20" = c("females_2000s", "females_2020s"), 
                 "F10_20" = c("females_2010s", "females_2020s"), 
                 "M80_90" = c("males_1980s", "males_1990s"), 
                 "M80_00" = c("males_1980s", "males_2000s"), 
                 "M80_10" = c("males_1980s", "males_2010s"), 
                 "M80_20" = c("males_1980s", "males_2020s"), 
                 "M90_00" = c("males_1990s", "males_2000s"), 
                 "M90_10" = c("males_1990s", "males_2010s"), 
                 "M90_20" = c("males_1990s", "males_2020s"), 
                 "M00_10" = c("males_2000s", "males_2010s"),
                 "M00_20" = c("males_2000s", "males_2020s"), 
                 "M10_20" = c("males_2010s", "males_2020s"))


inHR_comp <- lapply(HR_comps, function(x){meta(inHR_all[x], level = 0.95)})
hr_meta <- data.frame(comparison = names(inHR_comp), 
                      t(sapply(inHR_comp, function(x){
                        if(x[2,1,2] > x[1,2,2]){
                          c(low = x[2,1,1], est = x[2,1,2], high = x[2,1,3])
                          }else{
                            c(low = x[1,2,1], est = x[1,2,2], high = x[1,2,3])
                          }
})))
hr_meta$sig <- ifelse(hr_meta$low > 1 & hr_meta$high > 1, "YES", "NO")

sidak.level <- 1 - (1-0.05)^(1/length(HR_comps))

inHR_comp_sidak <- lapply(HR_comps[hr_meta$sig == "YES"], function(x){meta(inHR_all[x], level = 1 - sidak.level)})
hr_meta_sidak <- data.frame(comparison = names(inHR_comp_sidak), 
                      t(sapply(inHR_comp_sidak, function(x){
                        if(x[2,1,2] > x[1,2,2]){
                          c(low = x[2,1,1], est = x[2,1,2], high = x[2,1,3])
                        }else{
                          c(low = x[1,2,1], est = x[1,2,2], high = x[1,2,3])
                        }
                      })))
hr_meta_sidak$sig <- ifelse(hr_meta_sidak$low > 1 & hr_meta_sidak$high > 1, "YES", "NO")
View(hr_meta_sidak)

hr_meta$cor.sig[hr_meta$sig == "YES"] <- hr_meta_sidak$sig

#openxlsx::write.xlsx(hr_meta, file = paste("HomeRange_comp_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = ""))

# Step 10: Home Range Overlap ####
## View home range overlap for each group
overlaps <- lapply(overAll, function(x){
  e <- x$CI[,,"est"]
  l <- x$CI[,,"low"]
  h <- x$CI[,,"high"]
  
  if(length(x$CI[,,"est"]) > 1){
    comps <- lapply(1:ncol(e), function(i){
      do.call(rbind, lapply(i:nrow(e), function(j){data.frame(comp1 = colnames(e)[i], comp2 = row.names(e)[j])})[-1])
    })
    comps <- data.frame(do.call(rbind, comps))
    comps$low <- l[lower.tri(l)]
    comps$est <- e[lower.tri(e)]
    comps$high <- h[lower.tri(h)]
    comps$sig <- with(comps, ifelse(est > 0.137, "sig", "not"))
  }else{
    comps <- NULL
  }
  
  return(comps)
})
overlaps2 <- do.call(rbind, overlaps)
overlaps3 <- merge.data.frame(overlaps2, meta1[,c("ID", "CollarType")], by.x = "comp1", by.y = "ID")
overlaps3$Proximity <- ifelse(overlaps3$sig == "sig" & overlaps3$CollarType == "GPS", "YES", "NO")

overlaps4 <- overlaps3[overlaps3$comp1 %in% meta1$ID[meta1$Note %in% c("", "Found dead of natural causes - starvation", "died")],]
overlaps4 <- overlaps4[overlaps4$comp2 %in% meta1$ID[meta1$Note %in% c("", "Found dead of natural causes - starvation", "died")],]


overlaps4$type = paste(substr(overlaps4$comp1, start = nchar(overlaps4$comp1), stop = nchar(overlaps4$comp1)), 
                       substr(overlaps4$comp2, start = nchar(overlaps4$comp2), stop = nchar(overlaps4$comp2)), 
                       sep = "")
overlaps4$type <- factor(overlaps4$type, 
                         levels = c("FF", "MM", "FM", "MF", "MU", "UU", "UM"), 
                         labels = c("F", "M", "FM", "FM", "U", "U", "U"))

mean(overlaps4$est)
dplyr::summarise(dplyr::group_by(overlaps4, type), avg = mean(est))
dplyr::summarise(dplyr::group_by(overlaps4, type, CollarType), avg = mean(est))


## Save home range overlap for each group
openxlsx::write.xlsx(overlaps4, file = paste("Overlap_analysis_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = ""))


# Step 11: Calculate proximity of individuals to each other in time ####
### Explore the ctmm::difference, ctmm::distances, and ctmm::proximity functions

## Calculate proportion of same sex and opposite sex comparisons for VHF data
compsVHF <- do.call(rbind, lapply(1:nrow(overlaps4[overlaps4$sig == "sig" & overlaps4$CollarType == "VHF",]), function(i){
  a <- overlaps4[overlaps4$sig == "sig" & overlaps4$CollarType == "VHF",][i,1:2]
  if(substr(a[,1], start = nchar(a[,1]), stop = nchar(a[,1])) == substr(a[,2], start = nchar(a[,2]), stop = nchar(a[,2]))){
    g <- "GPS"
    if(substr(a[,1], start = nchar(a[,1]), stop = nchar(a[,1])) == "F"){
      S <- "F"
    }else{
      S <- "M"
    }
  }else{
    g <- "MF"
    S <- "MF"
  }
  
  data.frame(type = g, S, a)
}))
with(compsVHF, table(S))

## Define the groups that can be combined
comps <- do.call(c, lapply(1:nrow(overlaps4[overlaps4$Proximity == "YES",]), function(i){
  a <- overlaps4[overlaps4$Proximity == "YES",][i,1:2]
  if(substr(a[,1], start = nchar(a[,1]), stop = nchar(a[,1])) == substr(a[,2], start = nchar(a[,2]), stop = nchar(a[,2]))){
    g <- "GPS"
  }else{
    g <- "MF"
  }
  name <- paste(g, ".", a[,1], "_", substr(a[,2], 3,nchar(a[,2])), sep = "")
  a1 <- list(c(a[,1], a[,2]))
  names(a1) <- name
  return(a1)
}))

## Calculate Proximity
### Function for calculating proximity
proxFun <- function(comps, tel, sel, pp, cores.left){
  #comps <- comps_GPS[1:2]
  #tel <- tAllu
  #sel <- selAllu
  #pp <- FALSE
  #cores.left <- 4
  
  print(paste("This function started at ", Sys.time(), sep = ""))
  if(isTRUE(pp)){
    message("Parallel processing enabled")
    if (is.null(cores.left)) {
      cores.left <- 2
    }
    else {
      cores.left <- tryCatch(as.numeric(cores.left), error = function(e) {
        message("There was an error coercing cores.left to a number. The default of 2 cores are not utilized")
        return(2)
      }, warning = function(w) {
        message("Could not coerce cores.left to a number. The default of 2 cores are not utilized")
        return(2)
      })
    }
    cl1 <- parallel::makeCluster(parallel::detectCores() - cores.left, outfile = "out.txt")
    parallel::clusterExport(cl1, varlist = c("comps", "tel", "sel"), envir = environment())
  }else{
    cl1 <- NULL
  }
  
  out1 <- pbapply::pblapply(1:length(comps), cl = cl1, function(i){
    x <- comps[[i]]
    print(paste("Starting comparison between ", paste(x, collapse = " and "), sep = ""))
    dist <- ctmm::distances(tel[x], CTMM = sel[x])
    diff <- ctmm::difference(tel[x], CTMM = sel[x])
    prox <- ctmm::proximity(tel[x], CTMM = sel[x])
    out <- list(distance = dist, difference = diff, proximity = prox)
    saveRDS(out, file = paste("Proximity_", names(comps)[i], "_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
    print(paste("Finished comparison between ", paste(x, collapse = " and "), sep = ""))
    return(out)
  })
  
  if(isTRUE(pp)){
    parallel::stopCluster(cl1)
  }
  
  print(paste("This function completed at ", Sys.time(), sep = ""))
  return(out1)
  rm(cl1, out1)
  #rm(comps, tel, sel, pp, cores.left)
}
### Calculations for all overlaps
#proxAll <- proxFun(comps = comps, tel = tAll2, sel = selAll, pp = TRUE, cores.left = 4)
#proxAll_add <- proxFun(comps = comps[c(4,43)], tel = tAll2, sel = selAll, pp = TRUE, cores.left = 12)
#names(proxAll_add) <- lapply(comps[c(4,43)], function(x){paste(x, collapse="-")})
#names(proxAll) <- lapply(comps, function(x){paste(x, collapse="-")})
#proxAll <- c(proxAll, proxAll_add)
#saveRDS(proxAll, file = paste("Proximity_All_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
proxAll <- readRDS(file = file.path("Data", "Proximity_All_20241211.RDS"))

## Extracting distances, differences, and proximities
### Differences
diffsAll <- lapply(proxAll, function(x){x$difference})
### Distances
distsAll <- lapply(proxAll, function(x){x$distance})
### Proximity
proxsAll <- lapply(proxAll, function(x){x$proximity})
### Combine all together
proximityAll <- data.frame(comp = names(proxsAll), do.call(rbind, strsplit(names(proxsAll), split = "-")), type = NA, do.call(rbind, proxsAll))
colnames(proximityAll)[2:3] <- c("Ind1", "Ind2")
row.names(proximityAll) <- NULL
proximityAll$type <- with(proximityAll, ifelse(substr(Ind1, nchar(Ind1), nchar(Ind1)) == "M" & substr(Ind2, nchar(Ind2), nchar(Ind2)) == "M", "M", ifelse(substr(Ind1, nchar(Ind1), nchar(Ind1)) == "F" & substr(Ind2, nchar(Ind2), nchar(Ind2)) == "F", "F", "MF")))
proximityAll$sig <- with(proximityAll, ifelse(low < 1 & high < 1, "closer", ifelse(low > 1 & high > 1, "farther", "independent")))
### Write the output to excel
#openxlsx::write.xlsx(proximityAll, file = paste("Proximity_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = ""))

## Future steps: 
### Calculate average resident cat home ranges and compare through time

# Step 12: Create Plots ####
#library(ggplot2)

hrPlot <- function(ds, size.text, size.lab){
  #ds <- hr_time
  #ds <- hr_type
  #size.text <- 8
  #size.lab <- 15
  
  library(ggplot2)
  
  # The theme to make the graph look how we want it to
  theme.plot <- theme(text = element_text(family = "serif")) + 
    theme(plot.title = element_text(hjust = 0.5, size = size.lab*1.50, margin = margin(b = 0.5, unit = "inch")), 
          plot.background = element_rect(fill = "white", color = NULL), 
          plot.margin = unit(c(.1,.1,.1,.1), "inch")) +
    theme(axis.ticks = element_line(color = NA, linewidth = 1, linetype = "solid"), 
          axis.line = element_line(color = NA, linewidth = .1, linetype = "solid"), 
          axis.title=element_text(size=size.lab*1.15, margin = margin(t = 0.25, unit="inch")),  
          axis.title.x = element_text(vjust = -1.0, hjust = 0.50), 
          axis.title.y = element_text(angle = 90, hjust = 0.50, vjust = 2.0), 
          axis.text = element_text(size = size.lab, color = "black"), 
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.0), 
          axis.text.y = element_text(angle = 0, hjust = 0)) + 
    theme(panel.border = element_rect(fill = NA, color = "black"), 
          panel.background = element_rect(fill = NA, color = NA), 
          panel.grid.major = element_line(color = "grey75"), 
          #panel.grid.major.x = element_line(color = "grey75"), 
          #panel.grid.major.y = element_line(color = "grey75"), 
          panel.grid.minor = element_line(color = NA), 
          panel.grid.minor.x = element_line(color = NA), 
          panel.grid.minor.y = element_line(color = NA)) +  
    theme(legend.margin=margin(c(0.15,0.15,0.15,0.15), unit = "inch"), 
          legend.background = element_rect(fill = NA, color = NA), 
          legend.text=element_text(size = size.lab), 
          legend.title=element_text(size = size.lab), 
          legend.position = "right", 
          #legend.key = element_rect(color = "black", fill = NA), 
          legend.key.height = unit(0.25,"inch"), 
          legend.key.width = unit(0.25, "inch")) + 
    theme(strip.background = element_rect(fill = "gray85", color = "black"), 
          strip.text = element_text(size = size.lab), 
          strip.text.x = element_text(margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")), 
          strip.text.y = element_text(angle = -90, margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")))
  
  cl <- length(unique(ds$dataset))
  
  if(cl == 5){
    scale_col <- scale_color_manual("Decade", values = c("#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928"))
  }else if(cl == 2){
    scale_col <- scale_color_manual("Collar Type", values = c("#762A83", "#1B7837"))
  }
  
  plot1 <- ggplot(data = ds, mapping = aes(x = est, y = sex, xmin = low, xmax = high, color = dataset)) + 
    geom_errorbar(linewidth = 0.75, width = 0.25, position = position_dodge(width = 0.5)) + 
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) + 
    scale_x_continuous(expression("Estimated Home Range Size (km"^2*")"), breaks = seq(0,50,5), limits = c(0,25)) + 
    scale_y_discrete("Sex", labels = c("All cats", "Females", "Males")) + 
    scale_col + 
    theme.plot
  plot1
  
  return(plot1)
  rm(theme.plot, cl, plot1)
  #rm(ds, size.text, size.lab)
}


#### Home range estimates
# Average home range size through time for females and males
hr_time <- homerange[homerange$dataset %in% c("1980s", "1990s", "2000s", "2010s", "2020s") & homerange$sex != "sex",]

hr_time_plot <- hrPlot(ds = hr_time, size.text = 8, size.lab = 15)
hr_time_plot

# VHF vs GPS for males and females
hr_type <- homerange[homerange$dataset %in% c("VHF", "GPS") & homerange$sex != "sex",]

hr_type_plot <- hrPlot(ds = hr_type, size.text = 8, size.lab = 15)
hr_type_plot


# Males vs females
#hr_sex <- homerange[homerange$dataset %in% c("females2", "males2"),]

#hr_sex_plot <- hrPlot(ds = hr_sex, size.text = 8, size.lab = 15)
#hr_sex_plot

# Save the plots
ggsave(filename = paste("HR_time_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = ""), plot = hr_time_plot, device = "tiff", width = 6.5, height = 4, dpi = 600, compression = "lzw")
ggsave(filename = paste("HR_collar_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = ""), plot = hr_type_plot, device = "tiff", width = 6.5, height = 4, dpi = 600, compression = "lzw")


#### Overlap estimates
# 2 kdes that overlap vs 2 kdes that don't overlap
## Produced in ArcGIS

# Table of the 8 overlaps that are valid

## Supplemental figure: 
### Bar graph
  ### Collar type -> MM, FF, MF -> proportion of overlaps > 15%



#### Proximity estimates
# Dot and whisker plot for proximity with 1 being the center line
#proximityGPS_plot <- proximityGPS[c(3:4,7:8,14:16,20),]

proxPlot <- function(X, size.text, size.lab){
  #X <- proximityAll
  #size.text <- 4
  #size.lab <- 12
  
  library(ggplot2)
  
  # The theme to make the graph look how we want it to
  theme.plot <- theme(text = element_text(family = "serif")) + 
    theme(plot.title = element_text(hjust = 0.5, size = size.lab*1.50, margin = margin(b = 0.5, unit = "inch")), 
          plot.background = element_rect(fill = "white", color = NULL), 
          plot.margin = unit(c(.1,.1,.1,.1), "inch")) +
    theme(axis.ticks = element_line(color = NA, linewidth = 1, linetype = "solid"), 
          axis.line = element_line(color = NA, linewidth = .1, linetype = "solid"), 
          axis.title=element_text(size=size.lab*1.15, margin = margin(t = 0.25, unit="inch")),  
          axis.title.x = element_text(vjust = -1.0, hjust = 0.50), 
          axis.title.y = element_text(angle = 90, hjust = 0.50, vjust = 2.0), 
          axis.text = element_text(size = size.lab, color = "black"), 
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.0), 
          axis.text.y = element_text(angle = 0, hjust = 0)) + 
    theme(panel.border = element_rect(fill = NA, color = "black"), 
          panel.background = element_rect(fill = NA, color = NA), 
          panel.grid.major = element_line(color = NA), 
          #panel.grid.major.x = element_line(color = NA), 
          #panel.grid.major.y = element_line(color = "grey75"), 
          panel.grid.minor = element_line(color = NA), 
          panel.grid.minor.x = element_line(color = NA), 
          panel.grid.minor.y = element_line(color = NA)) +  
    theme(legend.margin=margin(c(0.15,0.15,0.15,0.15), unit = "inch"), 
          legend.background = element_rect(fill = NA, color = NA), 
          legend.text=element_text(size = size.lab), 
          legend.title=element_text(size = size.lab), 
          legend.position = "top", 
          #legend.key = element_rect(color = "black", fill = NA), 
          legend.key.height = unit(0.25,"inch"), 
          legend.key.width = unit(0.25, "inch")) + 
    theme(strip.background = element_rect(fill = "gray85", color = "black"), 
          strip.text = element_text(size = size.lab), 
          strip.text.x = element_text(margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")), 
          strip.text.y = element_text(angle = -90, margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")))
  
  X$type <- factor(X$type, labels = c("F" = "F-F", "M" = "M-M", "MF" = "F-M"))
  
  X2 <- split(X, f = X$type)
  X2 <- X2[c(3,2,1)]
  
  X3 <- do.call(rbind, X2)
  rownames(X3) <- NULL
  #X3$type <- as.factor(X3$type)
  X3$lab <- paste("comp", formatC(seq(1,nrow(X3)), width = 2, flag = "0"), sep = "")
  
  p <- lapply(X2, function(y){
    ggplot(y, mapping = aes(x = est, y = comp, xmin = low, xmax = high, color = sig)) + 
      geom_vline(xintercept = 1, linewidth = 1, color = "black") + 
      #geom_segment(mapping = aes(x = 1, xend = 1, y = -Inf, yend = Inf)) + 
      geom_point(position = position_dodge(width = 0.5)) + 
      geom_errorbar(width = 0.25, position = position_dodge(width = 0.5)) + 
      scale_x_continuous("Proximity Ratio", breaks = seq(0,10,0.5), limits = c(0,10), labels = c("", "Closer", "Independent", "Farther", rep("", times = 17))) + 
      scale_y_discrete("Comparison", labels = NULL) + 
      scale_color_manual("", values = c("closer" = "blue", "independent" = "black", "farther" = "red")) + 
      coord_cartesian(xlim = c(0,4)) + 
      #facet_grid(cols = vars(type)) + 
      theme.plot
  })
  p2 <- ggpubr::ggarrange(plotlist = p, nrow = 1, ncol = 3, common.legend = TRUE, labels = "AUTO")
  p2
  
  p <- ggplot(X3, mapping = aes(x = est, y = lab, xmin = low, xmax = high, color = sig)) + 
    geom_vline(xintercept = 1, linewidth = 1, color = "grey75") + 
    #geom_segment(mapping = aes(x = 1, xend = 1, y = -Inf, yend = Inf)) + 
    geom_point(position = position_dodge(width = 0.5), size = 0.75) + 
    geom_errorbar(width = 0.25, position = position_dodge(width = 0.5)) + 
    #scale_x_continuous("Proximity Ratio", breaks = seq(0,10,0.5), limits = c(0,10), labels = c("", "Closer", "Independent", "Farther", rep("", times = 17))) + 
    scale_x_continuous("Proximity Ratio", breaks = seq(0,10,1), limits = c(0,10)) + 
    scale_y_discrete("Comparison", labels = NULL) + 
    scale_color_manual("", values = c("closer" = "blue", "independent" = "black", "farther" = "red")) + 
    coord_cartesian(xlim = c(0,6)) + 
    facet_grid(cols = vars(type)) + 
    theme.plot
  p
  
  return(p)
  rm(theme.plot, X2, p, p2)
  #rm(X, size.lab, size.text)
}

prox_plot <- proxPlot(X = proximityAll, size.lab = 12, size.text = 3)
prox_plot
ggsave(filename = paste("Proximity_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = ""), plot = prox_plot, device = "tiff", width = 6.5, height = 5, dpi = 600, compression = "lzw")


#####################################################################################################################
############################################# Test data sets and old code ########################################### 
#####################################################################################################################


## Population-level home range based on year of data ####
popKDE <- do.call(rbind, lapply(1:length(kdeAll), function(i){
  x <- kdeAll[[i]]
  x1 <- ctmm::meta(x)[1,]
  data.frame(Location = names(kdeAll)[i], low = x1[1], est = x1[2], high = x1[3])
}))
popKDE

## Population-level home range based on ranch combined
ranchKDEs <- lapply(unique(meta1$Location), function(x){
  x1 <- do.call(c, kdeAll[grep(x, names(kdeAll))])
  names(x1) <- do.call(rbind, strsplit(names(x1), split = "[.]"))[,3]
  return(x1)
})
names(ranchKDEs) <- unique(meta1$Location)
ranchKDE <- lapply(ranchKDEs, ctmm::meta)
ranchKDE

## Save the home range information 
#openxlsx::write.xlsx(list(population = popKDE, individual = indHR2), file = paste("HomeRanges_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = ""))


# Home range overlaps ####
## Home range estimates (all resident individuals)
inHR_res <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES"]]
ctmm::meta(inHR_res)
### Remove WB39F and LB071M due to the reasons below
inHR_res2 <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
ctmm::meta(inHR_res2)

## Home range estimates (resident females)
inHR_F <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Sex == "F"]]
ctmm::meta(inHR_F)
### WB39F has a weird home range. She exists in two distinct clusters that she seemingly regularly moved between in 1996 and 1997, creating a larger home range than other females
plot(groups2$tel$a.Welder$WB39F, kdeAll$a.Welder$WB39F)
### Home range estimate (resident females, excluding WB39F)
inHR_F2 <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Sex == "F" & meta1$ID != "WB39F"]]
ctmm::meta(inHR_F2)

## Home range estimates (resident males)
inHR_M <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Sex == "M"]]
ctmm::meta(inHR_M)
### LB071M seems to have shifted its home range over time, accounting for its large overall home range
plot(tAll2$LB071M, indHR1$LB071M)
### Home range estimate (resident males, excluding LB071M)
inHR_M2 <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Sex == "M" & !(meta1$ID %in% c("LB071M"))]]
ctmm::meta(inHR_M2)

## Home range estimates (by decade)
### 1980s
inHR_1980s <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1980s" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_1980s_sex <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1980s" & meta1$Sex != "U" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_1980s_F <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1980s" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_1980s_M <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1980s" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
ctmm::meta(inHR_1980s)
ctmm::meta(inHR_1980s_sex)
ctmm::meta(inHR_1980s_F)
ctmm::meta(inHR_1980s_M)
### 1990s
inHR_1990s <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1990s" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_1990s_F <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1990s" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_1990s_M <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="1990s" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
ctmm::meta(inHR_1990s)
ctmm::meta(inHR_1990s_F)
ctmm::meta(inHR_1990s_M)
### 2000s
inHR_2000s <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2000s" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_2000s_F <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2000s" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_2000s_M <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2000s" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
ctmm::meta(inHR_2000s)
ctmm::meta(inHR_2000s_F)
ctmm::meta(inHR_2000s_M)
### 2010s
inHR_2010s <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2010s" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_2010s_F <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2010s" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_2010s_M <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2010s" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
ctmm::meta(inHR_2010s)
ctmm::meta(inHR_2010s_F)
ctmm::meta(inHR_2010s_M)
### 2020s
inHR_2020s <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2020s" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_2020s_F <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2020s" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_2020s_M <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$Decade=="2020s" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
ctmm::meta(inHR_2020s)
ctmm::meta(inHR_2020s_F)
ctmm::meta(inHR_2020s_M)

## Home range estimates by collar type
### VHF
inHR_VHF <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="VHF" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_VHF_sex <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="VHF" & meta1$Sex != "U" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_VHF_F <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="VHF" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_VHF_M <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="VHF" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
ctmm::meta(inHR_VHF)
ctmm::meta(inHR_VHF_sex)
ctmm::meta(inHR_VHF_F)
ctmm::meta(inHR_VHF_M)
### GPS
inHR_GPS <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="GPS" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_GPS_F <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="GPS" & meta1$Sex == "F" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
inHR_GPS_M <- indHR1[meta1$ID[!(meta1$Note %in% c("Exploratory movements", "Dispersal", "Collar Error?")) & meta1$Useable=="YES" & meta1$CollarType=="GPS" & meta1$Sex == "M" & !(meta1$ID %in% c("WB39F", "LB071M"))]]
ctmm::meta(inHR_GPS)
ctmm::meta(inHR_GPS_F)
ctmm::meta(inHR_GPS_M)




# Calculating proximity ####
### VHF Cats
#### CANNOT BE USED FOR FINE SCALE ANALYSIS. VHF DATA is not high quality enough
### GPS Cats (male/female overlap)
comps_FM <- list(es25M_26F = c("EB25M", "EB26F"), 
                 es26F_31M = c("EB26F", "EB31M"), 
                 es27M_29F = c("EB27M", "EB29F"), 
                 es27M_30F = c("EB27M", "EB30F"), 
                 es28F_31M = c("EB28F", "EB31M"), 
                 es29F_31M = c("EB29F", "EB31M"), 
                 es30F_31M = c("EB30F", "EB31M"), 
                 es32F_34M = c("EB32F", "EB34M"), 
                 es35M_41F = c("EB35M", "EB41F"), 
                 es36M_41F = c("EB36M", "EB41F"), 
                 es37M_40F = c("EB37M", "EB40F"), 
                 es37M_41F = c("EB37M", "EB41F"), 
                 es39M_41F = c("EB39M", "EB41F"), 
                 hc727054F_729230M = c("HB727054F", "HB729230M"), 
                 hc727054F_729233M = c("HB727054F", "HB729233M"), 
                 hc729230M_729232F = c("HB729230M", "HB729232F"), 
                 hc729230M_729234F = c("HB729230M", "HB729234F"), 
                 hc729232F_729233M = c("HB729232F", "HB729233M"), 
                 hc729233M_729234F = c("HB729233M", "HB729234F"), 
                 hr01F_03M = c("HR01F", "HR03M"), 
                 hr01F_06M = c("HR01F", "HR06M"), 
                 hr01F_08M = c("HR01F", "HR06M"), 
                 hr01F_10M = c("HR01F", "HR10M"), 
                 hr02F_03M = c("HR02F", "HR03M"), 
                 hr03M_04F = c("HR03M", "HR04F"), 
                 hr05F_06M = c("HR05F", "HR06M"), 
                 hr05F_08M = c("HR05F", "HR08M"), 
                 hr05F_10M = c("HR05F", "HR10M"), 
                 hr06M_07F = c("HR06M", "HR07F"), 
                 hr06M_09F = c("HR06M", "HR09F"), 
                 hr07F_08M = c("HR07F", "HR08M"),
                 hr07F_10M = c("HR07F", "HR10M"), 
                 hr09F_10M = c("HR09F", "HR10M"), 
                 kr02M_03F = c("KRB2M", "KRB3F"), 
                 kr05F_09M = c("KRB5F", "KRB9M"))
proxFM <- proxFun(comps = comps_FM, tel = tAll2, sel = selAll, pp = T, cores.left = 4)
names(proxFM) <- names(comps_FM)
saveRDS(proxFM, file = paste("Proximity_OppositeSex_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
proxFM <- readRDS(file.path("Data", "Proximity_OppositeSex_20240125.RDS"))

### GPS Cats (same sex overlap)
comps_GPS <- list(es25M_31M = c("EB25M", "EB31M"), 
                  es27M_31M = c("EB27M", "EB31M"), 
                  es29F_30F = c("EB29F", "EB30F"), 
                  es35M_36M = c("EB35M", "EB36M"), 
                  es37M_39M = c("EB37M", "EB39M"), 
                  es40F_41F = c("EB40F", "EB41F"), 
                  hc727054F_729232F = c("HB727054F", "HB729232F"), 
                  hc729230M_729233M = c("HB729230M", "HB729233M"), 
                  hr01F_02F = c("HR01F", "HR02F"), 
                  hr01F_04F = c("HR01F", "HR04F"), 
                  hr01F_05F = c("HR01F", "HR05F"), 
                  hr01F_07F = c("HR01F", "HR07F"), 
                  hr01F_09F = c("HR01F", "HR09F"), 
                  hr05F_07F = c("HR05F", "HR07F"), 
                  hr05F_09F = c("HR05F", "HR09F"),
                  hr07F_09F = c("HR07F", "HR09F"), 
                  hr06M_08M = c("HR06M", "HR08M"), 
                  hr06M_10M = c("HR06M", "HR10M"), 
                  hr08M_10M = c("HR08M", "HR10M"), 
                  la726981M_727047M = c("LB726981M", "LB727047M"))
proxGPS <- proxFun(comps = comps_GPS, tel = tAll2, sel = selAll, pp = T, cores.left = 4)
names(proxGPS) <- names(comps_GPS)
saveRDS(proxGPS, file = paste("Proximity_SameSex_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = ""))
proxGPS <- readRDS(file.path("Data", "Proximity_SameSex_20240125.RDS"))

### Differences
diffsFM <- lapply(proxFM, function(x){x$difference})
diffsGPS <- lapply(proxGPS, function(x){x$difference})

### Distances
distsFM <- lapply(proxFM, function(x){x$distance})
distsGPS <- lapply(proxGPS, function(x){x$distance})

### Proximity
proxsFM <- lapply(proxFM, function(x){x$proximity})
proxsGPS <- lapply(proxGPS, function(x){x$proximity})

proximityGPS <- do.call(rbind, proxsGPS)
proximityGPS <- data.frame(comp = row.names(proximityGPS), proximityGPS)
row.names(proximityGPS) <- NULL
proximityGPS$sig <- with(proximityGPS, ifelse(low < 1 & high < 1, "closer", ifelse(low > 1 & high > 1, "farther", "independent")))

proximityFM <- do.call(rbind, proxsFM)
proximityFM <- data.frame(comp = row.names(proximityFM), proximityFM)
row.names(proximityFM) <- NULL
proximityFM$sig <- with(proximityFM, ifelse(low < 1 & high < 1, "closer", ifelse(low > 1 & high > 1, "farther", "independent")))

#openxlsx::write.xlsx(rbind(proximityFM, proximityGPS), file = paste("Proximity_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = ""))
