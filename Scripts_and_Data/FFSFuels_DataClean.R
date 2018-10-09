## data cleaning file for FFSFuels Project 
## author - Jacob Levine

## Load packages
library(here)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(nlme)
library(lme4)
library(pscl)

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
calc_deltasBACI <- function() {
  ## initialize vectors
  fuelchangeBACI_litter <- numeric()
  fuelchangeBACI_duff <- numeric()
  fuelchangeBACI_litterduff <- numeric()
  fuelchangeBACI_1h <- numeric()
  fuelchangeBACI_10h <- numeric()
  fuelchangeBACI_100h <- numeric()
  fuelchangeBACI_1000h <- numeric()
  fuelchangeBACI_1.100h <- numeric()
  fuelchangeBACI_1.1000h <- numeric()
  fuelchangeBACI_total <- numeric()
  fuelchangeBACI_surface <- numeric()
  
  ## fill vectors with calculated values
  for (i in 1:length(unique(FFS_data$PlotID))) {
    for (j in 1:length(unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"]))) {
      if (unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j] == "Pretreatment" |
          unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j] == "Post_burn1") {
        fuelchangeBACI_litter <- c(fuelchangeBACI_litter, NA)
        fuelchangeBACI_duff <- c(fuelchangeBACI_duff, NA)
        fuelchangeBACI_litterduff <- c(fuelchangeBACI_litterduff, NA)
        fuelchangeBACI_1h <- c(fuelchangeBACI_1h, NA)
        fuelchangeBACI_10h <- c(fuelchangeBACI_10h, NA)
        fuelchangeBACI_100h <- c(fuelchangeBACI_100h, NA)
        fuelchangeBACI_1000h <- c(fuelchangeBACI_1000h, NA)
        fuelchangeBACI_1.100h <- c(fuelchangeBACI_1.100h, NA)
        fuelchangeBACI_1.1000h <- c(fuelchangeBACI_1.1000h, NA)
        fuelchangeBACI_total <- c(fuelchangeBACI_total, NA)
        fuelchangeBACI_surface <- c(fuelchangeBACI_surface, NA)
      }
      else {
          ## litter
          NewfuelchangeBACI_litter <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                              FFS_data$Timestep == "Post_7", "fuelload_litter_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                      FFS_data$Timestep == "Pretreatment", "fuelload_litter_tonha"]))
          fuelchangeBACI_litter <- c(fuelchangeBACI_litter, NewfuelchangeBACI_litter)
          
          ## duff
          NewfuelchangeBACI_duff <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                            FFS_data$Timestep == "Post_7", "fuelload_duff_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                  FFS_data$Timestep == "Pretreatment", "fuelload_duff_tonha"]))
          fuelchangeBACI_duff <- c(fuelchangeBACI_duff, NewfuelchangeBACI_duff)
          
          ## litterduff
          NewfuelchangeBACI_litterduff <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                  FFS_data$Timestep == "Post_7", "fuelload_litterduff_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                              FFS_data$Timestep == "Pretreatment", "fuelload_litterduff_tonha"]))   
          fuelchangeBACI_litterduff <- c(fuelchangeBACI_litterduff, NewfuelchangeBACI_litterduff)
          
          ## 1h
          NewfuelchangeBACI_1h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                          FFS_data$Timestep == "Post_7", "fuelload_1h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                              FFS_data$Timestep == "Pretreatment", "fuelload_1h_tonha"]))  
          fuelchangeBACI_1h <- c(fuelchangeBACI_1h, NewfuelchangeBACI_1h)
          
          ## 10h
          NewfuelchangeBACI_10h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                           FFS_data$Timestep == "Post_7", "fuelload_10h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                FFS_data$Timestep == "Pretreatment", "fuelload_10h_tonha"])) 
          fuelchangeBACI_10h <- c(fuelchangeBACI_10h, NewfuelchangeBACI_10h)
          
          ## 100h
          NewfuelchangeBACI_100h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                            FFS_data$Timestep == "Post_7", "fuelload_100h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                  FFS_data$Timestep == "Pretreatment", "fuelload_100h_tonha"]))   
          fuelchangeBACI_100h <- c(fuelchangeBACI_100h, NewfuelchangeBACI_100h)
          
          ## 1000r
          NewfuelchangeBACI_1000h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                             FFS_data$Timestep == "Post_7", "fuelload_1000h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                    FFS_data$Timestep == "Pretreatment", "fuelload_1000r_tonha"]))
          fuelchangeBACI_1000h <- c(fuelchangeBACI_1000h, NewfuelchangeBACI_1000h)
          
          ## 1.100h
          NewfuelchangeBACI_1.100h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                              FFS_data$Timestep == "Post_7", "fuelload_1.100h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                      FFS_data$Timestep == "Pretreatment", "fuelload_1.100h_tonha"]))  
          fuelchangeBACI_1.100h <- c(fuelchangeBACI_1.100h, NewfuelchangeBACI_1.100h)
          
          ## 1.1000h
          NewfuelchangeBACI_1.1000h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                               FFS_data$Timestep == "Post_7", "fuelload_1.1000h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                        FFS_data$Timestep == "Pretreatment", "fuelload_1.1000h_tonha"])) 
          fuelchangeBACI_1.1000h <- c(fuelchangeBACI_1.1000h, NewfuelchangeBACI_1.1000h)
          ## total
          NewfuelchangeBACI_total <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                             FFS_data$Timestep == "Post_7", "fuelload_total_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                    FFS_data$Timestep == "Pretreatment", "fuelload_total_tonha"]))  
          fuelchangeBACI_total <- c(fuelchangeBACI_total, NewfuelchangeBACI_total)
          
          ## surface
          NewfuelchangeBACI_surface <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                               FFS_data$Timestep == "Post_7", "fuelload_surface_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                        FFS_data$Timestep == "Pretreatment", "fuelload_surface_tonha"]))  
          fuelchangeBACI_surface <- c(fuelchangeBACI_surface, NewfuelchangeBACI_surface)
      }
    }
  }
  ## append new vectors to data:
  FFS_data$fuelchangeBACI_litter <- fuelchangeBACI_litter
  FFS_data$fuelchangeBACI_duff <- fuelchangeBACI_duff
  FFS_data$fuelchangeBACI_litterduff <- fuelchangeBACI_litterduff
  FFS_data$fuelchangeBACI_1h <- fuelchangeBACI_1h
  FFS_data$fuelchangeBACI_10h <- fuelchangeBACI_10h
  FFS_data$fuelchangeBACI_100h <- fuelchangeBACI_100h
  FFS_data$fuelchangeBACI_1000h <- fuelchangeBACI_1000h
  FFS_data$fuelchangeBACI_1.100h <- fuelchangeBACI_1.100h
  FFS_data$fuelchangeBACI_1.1000h <- fuelchangeBACI_1.1000h
  FFS_data$fuelchangeBACI_total <- fuelchangeBACI_total
  FFS_data$fuelchangeBACI_surface <- fuelchangeBACI_surface
  return(FFS_data)
}
FFS_data <- calc_deltasBACI()

calc_deltasBACI14 <- function() {
  ## initialize vectors
  fuelchangeBACI14_litter <- numeric()
  fuelchangeBACI14_duff <- numeric()
  fuelchangeBACI14_litterduff <- numeric()
  fuelchangeBACI14_1h <- numeric()
  fuelchangeBACI14_10h <- numeric()
  fuelchangeBACI14_100h <- numeric()
  fuelchangeBACI14_1000h <- numeric()
  fuelchangeBACI14_1.100h <- numeric()
  fuelchangeBACI14_1.1000h <- numeric()
  fuelchangeBACI14_total <- numeric()
  fuelchangeBACI14_surface <- numeric()
  
  ## fill vectors with calculated values
  for (i in 1:length(unique(FFS_data$PlotID))) {
    for (j in 1:length(unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"]))) {
      if (unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j] == "Pretreatment" |
          unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j] == "Post_burn1") {
        fuelchangeBACI14_litter <- c(fuelchangeBACI14_litter, NA)
        fuelchangeBACI14_duff <- c(fuelchangeBACI14_duff, NA)
        fuelchangeBACI14_litterduff <- c(fuelchangeBACI14_litterduff, NA)
        fuelchangeBACI14_1h <- c(fuelchangeBACI14_1h, NA)
        fuelchangeBACI14_10h <- c(fuelchangeBACI14_10h, NA)
        fuelchangeBACI14_100h <- c(fuelchangeBACI14_100h, NA)
        fuelchangeBACI14_1000h <- c(fuelchangeBACI14_1000h, NA)
        fuelchangeBACI14_1.100h <- c(fuelchangeBACI14_1.100h, NA)
        fuelchangeBACI14_1.1000h <- c(fuelchangeBACI14_1.1000h, NA)
        fuelchangeBACI14_total <- c(fuelchangeBACI14_total, NA)
        fuelchangeBACI14_surface <- c(fuelchangeBACI14_surface, NA)
      }
      else {
        ## litter
        NewfuelchangeBACI14_litter <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                 FFS_data$Timestep == "Post_14", "fuelload_litter_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                      FFS_data$Timestep == "Pretreatment", "fuelload_litter_tonha"]))
        fuelchangeBACI14_litter <- c(fuelchangeBACI14_litter, NewfuelchangeBACI14_litter)
        
        ## duff
        NewfuelchangeBACI14_duff <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                               FFS_data$Timestep == "Post_14", "fuelload_duff_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                  FFS_data$Timestep == "Pretreatment", "fuelload_duff_tonha"]))
        fuelchangeBACI14_duff <- c(fuelchangeBACI14_duff, NewfuelchangeBACI14_duff)
        
        ## litterduff
        NewfuelchangeBACI14_litterduff <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                     FFS_data$Timestep == "Post_14", "fuelload_litterduff_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                              FFS_data$Timestep == "Pretreatment", "fuelload_litterduff_tonha"]))   
        fuelchangeBACI14_litterduff <- c(fuelchangeBACI14_litterduff, NewfuelchangeBACI14_litterduff)
        
        ## 1h
        NewfuelchangeBACI14_1h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                             FFS_data$Timestep == "Post_14", "fuelload_1h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                              FFS_data$Timestep == "Pretreatment", "fuelload_1h_tonha"]))  
        fuelchangeBACI14_1h <- c(fuelchangeBACI14_1h, NewfuelchangeBACI14_1h)
        
        ## 10h
        NewfuelchangeBACI14_10h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                              FFS_data$Timestep == "Post_14", "fuelload_10h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                FFS_data$Timestep == "Pretreatment", "fuelload_10h_tonha"])) 
        fuelchangeBACI14_10h <- c(fuelchangeBACI14_10h, NewfuelchangeBACI14_10h)
        
        ## 100h
        NewfuelchangeBACI14_100h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                               FFS_data$Timestep == "Post_14", "fuelload_100h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                  FFS_data$Timestep == "Pretreatment", "fuelload_100h_tonha"]))   
        fuelchangeBACI14_100h <- c(fuelchangeBACI14_100h, NewfuelchangeBACI14_100h)
        
        ## 1000r
        NewfuelchangeBACI14_1000h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                FFS_data$Timestep == "Post_14", "fuelload_1000h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                    FFS_data$Timestep == "Pretreatment", "fuelload_1000r_tonha"]))
        fuelchangeBACI14_1000h <- c(fuelchangeBACI14_1000h, NewfuelchangeBACI14_1000h)
        
        ## 1.100h
        NewfuelchangeBACI14_1.100h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                 FFS_data$Timestep == "Post_14", "fuelload_1.100h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                      FFS_data$Timestep == "Pretreatment", "fuelload_1.100h_tonha"]))  
        fuelchangeBACI14_1.100h <- c(fuelchangeBACI14_1.100h, NewfuelchangeBACI14_1.100h)
        
        ## 1.1000h
        NewfuelchangeBACI14_1.1000h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                  FFS_data$Timestep == "Post_14", "fuelload_1.1000h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                        FFS_data$Timestep == "Pretreatment", "fuelload_1.1000h_tonha"])) 
        fuelchangeBACI14_1.1000h <- c(fuelchangeBACI14_1.1000h, NewfuelchangeBACI14_1.1000h)
        ## total
        NewfuelchangeBACI14_total <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                FFS_data$Timestep == "Post_14", "fuelload_total_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                    FFS_data$Timestep == "Pretreatment", "fuelload_total_tonha"]))  
        fuelchangeBACI14_total <- c(fuelchangeBACI14_total, NewfuelchangeBACI14_total)
        
        ## surface
        NewfuelchangeBACI14_surface <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                  FFS_data$Timestep == "Post_14", "fuelload_surface_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                        FFS_data$Timestep == "Pretreatment", "fuelload_surface_tonha"]))  
        fuelchangeBACI14_surface <- c(fuelchangeBACI14_surface, NewfuelchangeBACI14_surface)
      }
    }
  }
  ## append new vectors to data:
  FFS_data$fuelchangeBACI14_litter <- fuelchangeBACI14_litter
  FFS_data$fuelchangeBACI14_duff <- fuelchangeBACI14_duff
  FFS_data$fuelchangeBACI14_litterduff <- fuelchangeBACI14_litterduff
  FFS_data$fuelchangeBACI14_1h <- fuelchangeBACI14_1h
  FFS_data$fuelchangeBACI14_10h <- fuelchangeBACI14_10h
  FFS_data$fuelchangeBACI14_100h <- fuelchangeBACI14_100h
  FFS_data$fuelchangeBACI14_1000h <- fuelchangeBACI14_1000h
  FFS_data$fuelchangeBACI14_1.100h <- fuelchangeBACI14_1.100h
  FFS_data$fuelchangeBACI14_1.1000h <- fuelchangeBACI14_1.1000h
  FFS_data$fuelchangeBACI14_total <- fuelchangeBACI14_total
  FFS_data$fuelchangeBACI14_surface <- fuelchangeBACI14_surface
  return(FFS_data)
}
FFS_data <- calc_deltasBACI14()

## calculate rates of change:
calc_rates <- function() {
  ## initialize vectors
  fuelrate_litter <- numeric()
  fuelrate_duff <- numeric()
  fuelrate_litterduff <- numeric()
  fuelrate_1h <- numeric()
  fuelrate_10h <- numeric()
  fuelrate_100h <- numeric()
  fuelrate_1000h <- numeric()
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
        fuelrate_1000h <- c(fuelrate_1000h, NA)
        fuelrate_1.100h <- c(fuelrate_1.100h, NA)
        fuelrate_1.1000h <- c(fuelrate_1.1000h, NA)
        fuelrate_total <- c(fuelrate_total, NA)
        fuelrate_surface <- c(fuelrate_surface, NA)
      }
      else {
        timestepcur <- unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j]
        timestepcurnum <- unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "TimestepNum"])[j]
        timestepprev <- unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j-1]
        timestepprevnum <- unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "TimestepNum"])[j-1]
        
        if (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] & FFS_data$Timestep == timestepcur, "TimestepNum"] < FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] & FFS_data$Timestep == timestepprev, "TimestepNum"]) {
          print ("datasheet out of order, please sort in correct order of timesteps to continue")
        }
        else {
          ## litter
          Newfuelrate_litter <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                            FFS_data$Timestep == timestepcur, "fuelload_litter_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                    FFS_data$Timestep == timestepprev, "fuelload_litter_tonha"])/(timestepcurnum-timestepprevnum))
          fuelrate_litter <- c(fuelrate_litter, Newfuelrate_litter)
          
          ## duff
          Newfuelrate_duff <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                          FFS_data$Timestep == timestepcur, "fuelload_duff_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                FFS_data$Timestep == timestepprev, "fuelload_duff_tonha"])/(timestepcurnum-timestepprevnum))
          fuelrate_duff <- c(fuelrate_duff, Newfuelrate_duff)
          
          ## litterduff
          Newfuelrate_litterduff <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                FFS_data$Timestep == timestepcur, "fuelload_litterduff_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                            FFS_data$Timestep == timestepprev, "fuelload_litterduff_tonha"])/(timestepcurnum-timestepprevnum))  
          fuelrate_litterduff <- c(fuelrate_litterduff, Newfuelrate_litterduff)
          
          ## 1h
          Newfuelrate_1h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                        FFS_data$Timestep == timestepcur, "fuelload_1h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                            FFS_data$Timestep == timestepprev, "fuelload_1h_tonha"])/(timestepcurnum-timestepprevnum))   
          fuelrate_1h <- c(fuelrate_1h, Newfuelrate_1h)
          
          ## 10h
          Newfuelrate_10h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                         FFS_data$Timestep == timestepcur, "fuelload_10h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                              FFS_data$Timestep == timestepprev, "fuelload_10h_tonha"])/(timestepcurnum-timestepprevnum))   
          fuelrate_10h <- c(fuelrate_10h, Newfuelrate_10h)
          
          ## 100h
          Newfuelrate_100h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                          FFS_data$Timestep == timestepcur, "fuelload_100h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                FFS_data$Timestep == timestepprev, "fuelload_100h_tonha"])/(timestepcurnum-timestepprevnum))   
          fuelrate_100h <- c(fuelrate_100h, Newfuelrate_100h)
          
          ## 1000s
          Newfuelrate_1000h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                           FFS_data$Timestep == timestepcur, "fuelload_1000h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                  FFS_data$Timestep == timestepprev, "fuelload_1000h_tonha"])/(timestepcurnum-timestepprevnum))  
          fuelrate_1000h <- c(fuelrate_1000h, Newfuelrate_1000h)
          
          ## 1.100h
          Newfuelrate_1.100h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                            FFS_data$Timestep == timestepcur, "fuelload_1.100h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                    FFS_data$Timestep == timestepprev, "fuelload_1.100h_tonha"])/(timestepcurnum-timestepprevnum))   
          fuelrate_1.100h <- c(fuelrate_1.100h, Newfuelrate_1.100h)
          
          ## 1.1000h
          Newfuelrate_1.1000h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                             FFS_data$Timestep == timestepcur, "fuelload_1.1000h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                      FFS_data$Timestep == timestepprev, "fuelload_1.1000h_tonha"])/(timestepcurnum-timestepprevnum))    
          fuelrate_1.1000h <- c(fuelrate_1.1000h, Newfuelrate_1.1000h)
          
          ## total
          Newfuelrate_total <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                           FFS_data$Timestep == timestepcur, "fuelload_total_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                  FFS_data$Timestep == timestepprev, "fuelload_total_tonha"])/(timestepcurnum-timestepprevnum))    
          fuelrate_total <- c(fuelrate_total, Newfuelrate_total)
          
          ## surface
          Newfuelrate_surface <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                             FFS_data$Timestep == timestepcur, "fuelload_surface_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                      FFS_data$Timestep == timestepprev, "fuelload_surface_tonha"])/(timestepcurnum-timestepprevnum))   
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
  FFS_data$fuelrate_1000h <- fuelrate_1000h
  FFS_data$fuelrate_1.100h <- fuelrate_1.100h
  FFS_data$fuelrate_1.1000h <- fuelrate_1.1000h
  FFS_data$fuelrate_total <- fuelrate_total
  FFS_data$fuelrate_surface <- fuelrate_surface
  return(FFS_data)
  
}
FFS_data <- calc_rates()

## calculate deltas:
calc_deltas <- function() {
  ## initialize vectors
  fuelchange_litter <- numeric()
  fuelchange_duff <- numeric()
  fuelchange_litterduff <- numeric()
  fuelchange_1h <- numeric()
  fuelchange_10h <- numeric()
  fuelchange_100h <- numeric()
  fuelchange_1000h <- numeric()
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
        fuelchange_1000h <- c(fuelchange_1000h, NA)
        fuelchange_1.100h <- c(fuelchange_1.100h, NA)
        fuelchange_1.1000h <- c(fuelchange_1.1000h, NA)
        fuelchange_total <- c(fuelchange_total, NA)
        fuelchange_surface <- c(fuelchange_surface, NA)
      }
      else {
        timestepcur <- unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j]
        timestepcurnum <- unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "TimestepNum"])[j]
        timestepprev <- unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "Timestep"])[j-1]
        timestepprevnum <- unique(FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i], "TimestepNum"])[j-1]
        
        if (FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] & FFS_data$Timestep == timestepcur, "TimestepNum"] < FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] & FFS_data$Timestep == timestepprev, "TimestepNum"]) {
          print ("datasheet out of order, please sort in correct order of timesteps to continue")
        }
        else {
          ## litter
          Newfuelchange_litter <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                               FFS_data$Timestep == timestepcur, "fuelload_litter_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                       FFS_data$Timestep == timestepprev, "fuelload_litter_tonha"]))
          fuelchange_litter <- c(fuelchange_litter, Newfuelchange_litter)
          
          ## duff
          Newfuelchange_duff <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                             FFS_data$Timestep == timestepcur, "fuelload_duff_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                   FFS_data$Timestep == timestepprev, "fuelload_duff_tonha"]))
          fuelchange_duff <- c(fuelchange_duff, Newfuelchange_duff)
          
          ## litterduff
          Newfuelchange_litterduff <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                   FFS_data$Timestep == timestepcur, "fuelload_litterduff_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                               FFS_data$Timestep == timestepprev, "fuelload_litterduff_tonha"]))   
          fuelchange_litterduff <- c(fuelchange_litterduff, Newfuelchange_litterduff)
          
          ## 1h
          Newfuelchange_1h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                           FFS_data$Timestep == timestepcur, "fuelload_1h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                               FFS_data$Timestep == timestepprev, "fuelload_1h_tonha"]))  
          fuelchange_1h <- c(fuelchange_1h, Newfuelchange_1h)
          
          ## 10h
          Newfuelchange_10h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                            FFS_data$Timestep == timestepcur, "fuelload_10h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                 FFS_data$Timestep == timestepprev, "fuelload_10h_tonha"])) 
          fuelchange_10h <- c(fuelchange_10h, Newfuelchange_10h)
          
          ## 100h
          Newfuelchange_100h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                             FFS_data$Timestep == timestepcur, "fuelload_100h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                   FFS_data$Timestep == timestepprev, "fuelload_100h_tonha"]))   
          fuelchange_100h <- c(fuelchange_100h, Newfuelchange_100h)
          
          ## 1000r
          Newfuelchange_1000h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                              FFS_data$Timestep == timestepcur, "fuelload_1000h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                     FFS_data$Timestep == timestepprev, "fuelload_1000r_tonha"]))
          fuelchange_1000h <- c(fuelchange_1000h, Newfuelchange_1000h)
          
          ## 1.100h
          Newfuelchange_1.100h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                               FFS_data$Timestep == timestepcur, "fuelload_1.100h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                       FFS_data$Timestep == timestepprev, "fuelload_1.100h_tonha"]))  
          fuelchange_1.100h <- c(fuelchange_1.100h, Newfuelchange_1.100h)
          
          ## 1.1000h
          Newfuelchange_1.1000h <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                FFS_data$Timestep == timestepcur, "fuelload_1.1000h_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                         FFS_data$Timestep == timestepprev, "fuelload_1.1000h_tonha"])) 
          fuelchange_1.1000h <- c(fuelchange_1.1000h, Newfuelchange_1.1000h)
          ## total
          Newfuelchange_total <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                              FFS_data$Timestep == timestepcur, "fuelload_total_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                     FFS_data$Timestep == timestepprev, "fuelload_total_tonha"]))  
          fuelchange_total <- c(fuelchange_total, Newfuelchange_total)
          
          ## surface
          Newfuelchange_surface <- ((FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                FFS_data$Timestep == timestepcur, "fuelload_surface_tonha"] - FFS_data[FFS_data$PlotID == unique(FFS_data$PlotID)[i] &
                                                                                                                         FFS_data$Timestep == timestepprev, "fuelload_surface_tonha"]))  
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
  FFS_data$fuelchange_1000h <- fuelchange_1000h
  FFS_data$fuelchange_1.100h <- fuelchange_1.100h
  FFS_data$fuelchange_1.1000h <- fuelchange_1.1000h
  FFS_data$fuelchange_total <- fuelchange_total
  FFS_data$fuelchange_surface <- fuelchange_surface
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
FFS_data$logfuelload_1000h <- log(FFS_data$fuelload_1000h_tonha + 1)
FFS_data$logfuelload_1.100h <- log(FFS_data$fuelload_1.100h_tonha + 1)
FFS_data$logfuelload_1.1000h <- log(FFS_data$fuelload_1.1000h_tonha + 1)
FFS_data$logfuelload_surface <- log(FFS_data$fuelload_surface_tonha + 1)

## put treatments in right order:
FFS_data$Treatment <- factor(FFS_data$Treatment, levels = c("Control", "Burn", "Mech", "MechBurn"))
