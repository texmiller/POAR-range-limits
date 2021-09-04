library(raster)
library(tidyverse)
in_dir  <- "C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/GIS/"

# read and format coordinate data
sites     <- read.csv("POAR-range-limits/data/demography_allsites.csv") %>% 
              dplyr::select( site, Latitude, Longitude ) %>% 
              unique %>% 
              mutate( coord_type = 'Field site' )
collect   <- read.csv( paste0(in_dir, 'GPSPoarData.csv') ) %>% 
                mutate( coord_type = 'Collection' ) %>% 
                rename( site = Object.ID )
coord_df  <- bind_rows( sites, collect ) %>% 
               rename(  code = Description,
                        lat  = Latitude,
                        lon  = Longitude )

# format site coordinates for raster extraction
site_coord <-  matrix(c(coord_df$lon,
                        coord_df$lat),
                      dimnames = list(rep('value',nrow(coord_df)),
                                      c('Long','Lat') ),
                      byrow = FALSE, nrow = nrow(coord_df) )

# extract year and month from CHELSA file names
get_yr_month <- function(x){
  
  regmatches(x, 
             gregexpr("[0-9]{4}_[0-9]{2}", 
                      x) ) %>% 
    unlist %>% 
    strsplit( '_' ) %>% 
    unlist
}

# names of CHELSA files
in_dir <- 'P:/iDiv_GLOBAL_GEODATA/CHELSA_Datasets/Timeseries_1979-2013/CHELSA_prec'
file_l <- list.files( in_dir )

# year/month data frame
yr_mon_df <- lapply(file_l, get_yr_month) %>% 
               do.call( rbind, .) %>% 
               as.data.frame() %>% 
               setNames( c('year', 'month') )


# Download --------------------------------------------------------

# list of data frames with precipitation values
clim_l <- list()

# loop through the 420 files!
for( ii in 1:length(file_l) ){
  
  repP        <- raster( paste0(in_dir, '/', file_l[ii]) )
  values_clim <- raster::extract( repP, site_coord, 
                                  method = 'bilinear')
  clim_l[[ii]]<- coord_df %>% 
                  mutate( year     = yr_mon_df$year[ii],
                          month    = yr_mon_df$month[ii],
                          value    = values_clim)
  print( ii )
  
}

# store it
clim_df <- clim_l %>% bind_rows()

write.csv(clim_df, 'C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/Monitoring/monthly_climate.csv', row.names=F)
