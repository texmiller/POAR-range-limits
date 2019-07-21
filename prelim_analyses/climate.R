# script to read and calculate weather data
# currently problems with heterogeneity in the data (mesowest vs. GHCN data)
# and data sets completely missing data:
# (e.g. elreno, katy, llano, llela, ninnescah, ozona(precip)) 
setwd("C:/Users/ac79/Downloads/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography/data/climate")
options(stringsAsFactors = F)
library(tidyr)
library(dplyr)
library(magrittr)

# identify folder names
not_folders <- grep("\\.",list.files())
folder_n    <- list.files()[-not_folders]

# read climate data within each folder 
read_climate <- function(x){ 
  
  fil <- list.files(x)
  id  <- grep(".csv", fil)
  d   <- read.csv( paste0(x, "/",fil[id]), skip = 15)
  
  return(d)
  
}

# growing season percipitation and 
format_gs <- function(x){
  
  # split date information
  tmp     <- x %>%
              select(Day) %>%
              separate(Day, c("year","month","day"), sep = "-")
  # Remove "Day"
  x       <- cbind(tmp, x) %>% select(-Day) 
  # transform to numeric
  x[,]    <- lapply(x[,],function(x) x %<>% as.numeric(x)) 
  # ignore months after may, before september
  x       <- subset(x, month < 6 | month > 9)

  #growing season "split": before and after May! 
  gs      <- subset(x, month > 8) %>%
              mutate(year = year + 1) %>%
              bind_rows( subset(x, month < 5) )
  precip  <-  gs %>%
                group_by(year) %>%
                summarise(precip = sum(Precipitation, na.rm=T) )
  gd_days <- gs %>%
                group_by(year) %>%
                summarise(gd_days = (sum(Min.Temperature,Max.Temperature,na.rm=T) / 2) - 4 )
  
  out     <- left_join(precip, gd_days) %>%
              subset(year > 2009) %>%
              as.data.frame()
  
  return(as.data.frame(out))
  
}

# loop through folders
big_l       <- lapply(folder_n, read_climate)

# create growing season data
gs_climate  <- lapply(big_l,format_gs) %>%
                setNames(folder_n)
