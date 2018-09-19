## data cleaning file for FFSFuels Project 
## author - Jacob Levine

## Load packages
library(here)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(nlme)

## Load in data
path <- here("/Fuels_byplot_Plots_ffs_filtered_postburns.csv")
FFS_data <- read.csv(path, header = TRUE)

## remove "intermediate measurements", not going to be useful for analysis
FFS_data <- subset(FFS_data, FFS_data$Timestep != "Intermediate")

## Remove "Post_burn1", exact replica of Post_1 in all cases
FFS_data <- subset(FFS_data, Timestep != "Post_burn1")

## create numeric vector for timestep:
create_timenum <- function() {
  TimestepNum <- numeric()
  for (i in 1:nrow(FFS_data)) {
    if (FFS_data$Timestep[i] == "Pretreatment") {
      TimestepNum <- c(TimestepNum, -1)
    }
    else if (FFS_data$Timestep[i] == "Post_1") {
      TimestepNum <- c(TimestepNum, 1)
    }
    else if (FFS_data$Timestep[i] == "Post_burn2") {
      TimestepNum <- c(TimestepNum, 2)
    }
    else if (FFS_data$Timestep[i] == "Post_burn3") {
      TimestepNum <- c(TimestepNum, 3)
    }
    else if (FFS_data$Timestep[i] == "Post_7") { 
      TimestepNum <- c(TimestepNum, 7)
    }
    else if (FFS_data$Timestep[i] == "Post_8") {
      TimestepNum <- c(TimestepNum, 8)
    }
    else if (FFS_data$Timestep[i] == "Post_14") {
      TimestepNum <- c(TimestepNum, 14)
    }
  }
  return(TimestepNum)
}
TimestepNum <- create_timenum()
FFS_data$TimestepNum <- TimestepNum
FFS_data[FFS_data$Unit == 40, "Timestep"]

## put Timestep in correct order:
FFS_data$Timestep <- factor(FFS_data$Timestep, levels = c("Pretreatment", "Post_1", "Post_burn1", "Post_burn2", "Post_burn3", "Post_2", "Post_3", "Post_7", "Post_8", "Post_14"))

FFS_data[FFS_data$Timestep == "Post_burn1", "Timestep"] <- "Post_1"
FFS_data[FFS_data$Timestep == "Post_burn2", "Timestep"] <- "Post_2"
FFS_data[FFS_data$Timestep == "Post_burn3", "Timestep"] <- "Post_3"

## remove Post_2, 3, and 8, only have data from Fire treatments
FFS_data <- subset(FFS_data, FFS_data$Timestep != "Post_2" &
                     FFS_data$Timestep != "Post_3" &
                     FFS_data$Timestep != "Post_8")

## calculate deltas:
calc_deltas <- function() {
  ## initialize vectors
  fuelchange_litter <- numeric()
  fuelchange_duff <- numeric()
  fuelchange_litterduff <- numeric()
  fuelchange_1h <- numeric()
  fuelchange_10h <- numeric()
  fuelchange_100h <- numeric()
  fuelchange_1000s <- numeric()
  fuelchange_1000r <- numeric()
  fuelchange_1.100h <- numeric()
  fuelchange_1.1000h <- numeric()
  fuelchange_total <- numeric()
  fuelchange_surface <- numeric()
  
  ## fill vectors with calculated values
  for (i in 1:length(unique(FFS_data$PlotID))) {
    for (j in 1:length(unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"]))) {
      if (unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j] == "Pretreatment" |
          unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j] == "Post_burn1") {
        fuelchange_litter <- c(fuelchange_litter, NA)
        fuelchange_duff <- c(fuelchange_duff, NA)
        fuelchange_litterduff <- c(fuelchange_litterduff, NA)
        fuelchange_1h <- c(fuelchange_1h, NA)
        fuelchange_10h <- c(fuelchange_10h, NA)
        fuelchange_100h <- c(fuelchange_100h, NA)
        fuelchange_1000s <- c(fuelchange_1000s, NA)
        fuelchange_1000r <- c(fuelchange_1000r, NA)
        fuelchange_1.100h <- c(fuelchange_1.100h, NA)
        fuelchange_1.1000h <- c(fuelchange_1.1000h, NA)
        fuelchange_total <- c(fuelchange_total, NA)
        fuelchange_surface <- c(fuelchange_surface, NA)
      }
      else {
        timestepcur <- unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j]
        timestepprev <- unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j-1]
        
        if (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] & FFS_data$Timestep == timestepcur, "TimestepNum"] < FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] & FFS_data$Timestep == timestepprev, "TimestepNum"]) {
          print ("datasheet out of order, please sort in correct order of timesteps to continue")
        }
        else {
          ## litter
          Newfuelchange_litter <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                              FFS_data$Timestep == timestepcur, "fuelload_litter_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                      FFS_data$Timestep == timestepprev, "fuelload_litter_tonha"])  
          fuelchange_litter <- c(fuelchange_litter, Newfuelchange_litter)
          
          ## duff
          Newfuelchange_duff <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                            FFS_data$Timestep == timestepcur, "fuelload_duff_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                  FFS_data$Timestep == timestepprev, "fuelload_duff_tonha"]) 
          fuelchange_duff <- c(fuelchange_duff, Newfuelchange_duff)
          
          ## litterduff
          Newfuelchange_litterduff <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                  FFS_data$Timestep == timestepcur, "fuelload_litterduff_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                              FFS_data$Timestep == timestepprev, "fuelload_litterduff_tonha"])   
          fuelchange_litterduff <- c(fuelchange_litterduff, Newfuelchange_litterduff)
          
          ## 1h
          Newfuelchange_1h <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                          FFS_data$Timestep == timestepcur, "fuelload_1h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                              FFS_data$Timestep == timestepprev, "fuelload_1h_tonha"])  
          fuelchange_1h <- c(fuelchange_1h, Newfuelchange_1h)
          
          ## 10h
          Newfuelchange_10h <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                           FFS_data$Timestep == timestepcur, "fuelload_10h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                FFS_data$Timestep == timestepprev, "fuelload_10h_tonha"]) 
          fuelchange_10h <- c(fuelchange_10h, Newfuelchange_10h)
          
          ## 100h
          Newfuelchange_100h <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                            FFS_data$Timestep == timestepcur, "fuelload_100h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                  FFS_data$Timestep == timestepprev, "fuelload_100h_tonha"])   
          fuelchange_100h <- c(fuelchange_100h, Newfuelchange_100h)
          
          ## 1000s
          Newfuelchange_1000s <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                             FFS_data$Timestep == timestepcur, "fuelload_1000s_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                    FFS_data$Timestep == timestepprev, "fuelload_1000s_tonha"]) 
          fuelchange_1000s <- c(fuelchange_1000s, Newfuelchange_1000s)
          
          ## 1000r
          Newfuelchange_1000r <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                             FFS_data$Timestep == timestepcur, "fuelload_1000r_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                    FFS_data$Timestep == timestepprev, "fuelload_1000r_tonha"]) 
          fuelchange_1000r <- c(fuelchange_1000r, Newfuelchange_1000r)
          
          ## 1.100h
          Newfuelchange_1.100h <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                              FFS_data$Timestep == timestepcur, "fuelload_1.100h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                      FFS_data$Timestep == timestepprev, "fuelload_1.100h_tonha"])  
          fuelchange_1.100h <- c(fuelchange_1.100h, Newfuelchange_1.100h)
          
          ## 1.1000h
          Newfuelchange_1.1000h <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                               FFS_data$Timestep == timestepcur, "fuelload_1.1000h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                        FFS_data$Timestep == timestepprev, "fuelload_1.1000h_tonha"]) 
          
          ## total
          Newfuelchange_total <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                             FFS_data$Timestep == timestepcur, "fuelload_total_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                    FFS_data$Timestep == timestepprev, "fuelload_total_tonha"])  
          fuelchange_total <- c(fuelchange_total, Newfuelchange_total)
          
          ## surface
          Newfuelchange_surface <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                               FFS_data$Timestep == timestepcur, "fuelload_surface_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                        FFS_data$Timestep == timestepprev, "fuelload_surface_tonha"])  
          fuelchange_surface <- c(fuelchange_surface, Newfuelchange_surface)
        }
      }
    }
  }
  ## append new vectors to data:
  FFS_data$fuelchange_litter <- fuelchange_litter
  FFS_data$fuelchange_duff <- fuelchange_duff
  FFS_data$fuelchange_litterduff <- fuelchange_litterduff
  FFS_data$fuelchange_1h <- fuelchange_1h
  FFS_data$fuelchange_10h <- fuelchange_10h
  FFS_data$fuelchange_100h <- fuelchange_100h
  FFS_data$fuelchange_1000s <- fuelchange_1000s
  FFS_data$fuelchange_1000r <- fuelchange_1000r
  FFS_data$fuelchange_1.100h <- fuelchange_1.100h
  FFS_data$fuelchange_1.1000h <- fuelchange_1.1000h
  FFS_data$fuelchange_total <- fuelchange_total
  FFS_data$fuelchange_surface <- fuelchange_surface
  return(FFS_data)
  
}
FFS_data <- calc_deltas()

## calculate rates of change:
calc_rates <- function() {
  ## initialize vectors
  fuelrate_litter <- numeric()
  fuelrate_duff <- numeric()
  fuelrate_litterduff <- numeric()
  fuelrate_1h <- numeric()
  fuelrate_10h <- numeric()
  fuelrate_100h <- numeric()
  fuelrate_1000s <- numeric()
  fuelrate_1000r <- numeric()
  fuelrate_1.100h <- numeric()
  fuelrate_1.1000h <- numeric()
  fuelrate_total <- numeric()
  fuelrate_surface <- numeric()
  
  ## fill vectors with calculated values
  for (i in 1:length(unique(FFS_data$PlotID))) {
    for (j in 1:length(unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"]))) {
      if (unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j] == "Pretreatment" |
          unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j] == "Post_burn1") {
        fuelrate_litter <- c(fuelrate_litter, NA)
        fuelrate_duff <- c(fuelrate_duff, NA)
        fuelrate_litterduff <- c(fuelrate_litterduff, NA)
        fuelrate_1h <- c(fuelrate_1h, NA)
        fuelrate_10h <- c(fuelrate_10h, NA)
        fuelrate_100h <- c(fuelrate_100h, NA)
        fuelrate_1000s <- c(fuelrate_1000s, NA)
        fuelrate_1000r <- c(fuelrate_1000r, NA)
        fuelrate_1.100h <- c(fuelrate_1.100h, NA)
        fuelrate_1.1000h <- c(fuelrate_1.1000h, NA)
        fuelrate_total <- c(fuelrate_total, NA)
        fuelrate_surface <- c(fuelrate_surface, NA)
      }
      else {
        timestepcur <- unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j]
        timestepprev <- unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j-1]
        
        if (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] & FFS_data$Timestep == timestepcur, "TimestepNum"] < FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] & FFS_data$Timestep == timestepprev, "TimestepNum"]) {
          print ("datasheet out of order, please sort in correct order of timesteps to continue")
        }
        else {
          ## litter
          Newfuelrate_litter <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                            FFS_data$Timestep == timestepcur, "fuelload_litter_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                    FFS_data$Timestep == timestepprev, "fuelload_litter_tonha"]) / FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                                                                                              FFS_data$Timestep == timestepprev, "fuelload_litter_tonha"]  
          fuelrate_litter <- c(fuelrate_litter, Newfuelrate_litter)
          
          ## duff
          Newfuelrate_duff <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                          FFS_data$Timestep == timestepcur, "fuelload_duff_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                FFS_data$Timestep == timestepprev, "fuelload_duff_tonha"]) / FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                                                                                        FFS_data$Timestep == timestepprev, "fuelload_duff_tonha"]
          fuelrate_duff <- c(fuelrate_duff, Newfuelrate_duff)
          
          ## litterduff
          Newfuelrate_litterduff <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                FFS_data$Timestep == timestepcur, "fuelload_litterduff_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                            FFS_data$Timestep == timestepprev, "fuelload_litterduff_tonha"]) / FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                                                                                                          FFS_data$Timestep == timestepprev, "fuelload_litterduff_tonha"]  
          fuelrate_litterduff <- c(fuelrate_litterduff, Newfuelrate_litterduff)
          
          ## 1h
          Newfuelrate_1h <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                        FFS_data$Timestep == timestepcur, "fuelload_1h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                            FFS_data$Timestep == timestepprev, "fuelload_1h_tonha"]) / FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                                                                                  FFS_data$Timestep == timestepprev, "fuelload_1h_tonha"] 
          fuelrate_1h <- c(fuelrate_1h, Newfuelrate_1h)
          
          ## 10h
          Newfuelrate_10h <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                         FFS_data$Timestep == timestepcur, "fuelload_10h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                              FFS_data$Timestep == timestepprev, "fuelload_10h_tonha"]) / FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                                                                                     FFS_data$Timestep == timestepprev, "fuelload_10h_tonha"] 
          fuelrate_10h <- c(fuelrate_10h, Newfuelrate_10h)
          
          ## 100h
          Newfuelrate_100h <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                          FFS_data$Timestep == timestepcur, "fuelload_100h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                FFS_data$Timestep == timestepprev, "fuelload_100h_tonha"]) / FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                                                                                        FFS_data$Timestep == timestepprev, "fuelload_100h_tonha"]  
          fuelrate_100h <- c(fuelrate_100h, Newfuelrate_100h)
          
          ## 1000s
          Newfuelrate_1000s <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                           FFS_data$Timestep == timestepcur, "fuelload_1000s_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                  FFS_data$Timestep == timestepprev, "fuelload_1000s_tonha"]) / FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                                                                                           FFS_data$Timestep == timestepprev, "fuelload_1000s_tonha"]
          fuelrate_1000s <- c(fuelrate_1000s, Newfuelrate_1000s)
          
          ## 1000r
          Newfuelrate_1000r <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                           FFS_data$Timestep == timestepcur, "fuelload_1000r_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                  FFS_data$Timestep == timestepprev, "fuelload_1000r_tonha"]) / FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                                                                                           FFS_data$Timestep == timestepprev, "fuelload_1000r_tonha"]
          fuelrate_1000r <- c(fuelrate_1000r, Newfuelrate_1000r)
          
          ## 1.100h
          Newfuelrate_1.100h <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                            FFS_data$Timestep == timestepcur, "fuelload_1.100h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                    FFS_data$Timestep == timestepprev, "fuelload_1.100h_tonha"]) / FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                                                                                              FFS_data$Timestep == timestepprev, "fuelload_1.100h_tonha"] 
          fuelrate_1.100h <- c(fuelrate_1.100h, Newfuelrate_1.100h)
          
          ## 1.1000h
          Newfuelrate_1.1000h <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                             FFS_data$Timestep == timestepcur, "fuelload_1.1000h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                      FFS_data$Timestep == timestepprev, "fuelload_1.1000h_tonha"]) / FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                                                                                                 FFS_data$Timestep == timestepprev, "fuelload_1.1000h_tonha"]  
          fuelrate_1.1000h <- c(fuelrate_1.1000h, Newfuelrate_1.1000h)
          
          ## total
          Newfuelrate_total <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                           FFS_data$Timestep == timestepcur, "fuelload_total_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                  FFS_data$Timestep == timestepprev, "fuelload_total_tonha"]) / FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                                                                                           FFS_data$Timestep == timestepprev, "fuelload_total_tonha"]  
          fuelrate_total <- c(fuelrate_total, Newfuelrate_total)
          
          ## surface
          Newfuelrate_surface <- (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                             FFS_data$Timestep == timestepcur, "fuelload_surface_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                      FFS_data$Timestep == timestepprev, "fuelload_surface_tonha"]) / FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                                                                                                 FFS_data$Timestep == timestepprev, "fuelload_surface_tonha"] 
          fuelrate_surface <- c(fuelrate_surface, Newfuelrate_surface)
        }
      }
    }
  }
  ## append new vectors to data:
  FFS_data$fuelrate_litter <- fuelrate_litter
  FFS_data$fuelrate_duff <- fuelrate_duff
  FFS_data$fuelrate_litterduff <- fuelrate_litterduff
  FFS_data$fuelrate_1h <- fuelrate_1h
  FFS_data$fuelrate_10h <- fuelrate_10h
  FFS_data$fuelrate_100h <- fuelrate_100h
  FFS_data$fuelrate_1000s <- fuelrate_1000s
  FFS_data$fuelrate_1000r <- fuelrate_1000r
  FFS_data$fuelrate_1.100h <- fuelrate_1.100h
  FFS_data$fuelrate_1.1000h <- fuelrate_1.1000h
  FFS_data$fuelrate_total <- fuelrate_total
  FFS_data$fuelrate_surface <- fuelrate_surface
  return(FFS_data)
  
}
FFS_data <- calc_deltas()

## log+1 transform fuel stock data:
FFS_data$logfuelload_total <- log(FFS_data$fuelload_total_tonha + 1)
FFS_data$logfuelload_litter <- log(FFS_data$fuelload_litter_tonha + 1)
FFS_data$logfuelload_duff <- log(FFS_data$fuelload_duff_tonha + 1)
FFS_data$logfuelload_litterduff <- log(FFS_data$fuelload_litterduff_tonha + 1)
FFS_data$logfuelload_1h <- log(FFS_data$fuelload_1h_tonha + 1)
FFS_data$logfuelload_10h <- log(FFS_data$fuelload_10h_tonha + 1)
FFS_data$logfuelload_100h <- log(FFS_data$fuelload_100h_tonha + 1)
FFS_data$logfuelload_1000s <- log(FFS_data$fuelload_1000s_tonha + 1)
FFS_data$logfuelload_1000r <- log(FFS_data$fuelload_1000r_tonha + 1)
FFS_data$logfuelload_1.100h <- log(FFS_data$fuelload_1.100h_tonha + 1)
FFS_data$logfuelload_1.1000h <- log(FFS_data$fuelload_1.1000h_tonha + 1)
FFS_data$logfuelload_surface <- log(FFS_data$fuelload_surface_tonha + 1)

