library(raster)
library(tidyverse)
library(ggthemes)
in_dir  <- "C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/GIS/"

# read data
clim_df <- read.csv( 'C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/Monitoring/monthly_climate.csv' ) 
#              # Update "double codes"
#              mutate( code = replace( code,
#                                      full_name == 'buffalolake',
#                                      'TAMU_N') ) %>% 
#              mutate( code = replace( code,
#                                      full_name == 'Ozona',
#                                      'TAMU_S') ) %>% 
#              mutate( coord_type = replace( coord_type,
#                                            coord_type == 'origin',
#                                            'Collection' ) )  %>% 
#              mutate( coord_type = replace( coord_type,
#                                            coord_type == 'site',
#                                            'Field site' ) )

# Climatic yearly averages for all sites
clim_yr_avg  <- clim_df %>% 
  group_by( site, lon, lat, coord_type, year ) %>% 
  summarise( ppt_year = sum(value) ) %>% 
  ungroup %>% 
  group_by( site, lon, lat, coord_type ) %>% 
  summarise( ppt_avg = mean(ppt_year) ) %>% 
  ungroup  

# Climatic yearly averages per COLLECTION SITE
clim_yr_avg_collection  <- clim_df %>% 
  subset( coord_type == 'Collection' ) %>% 
  group_by( site, code, year ) %>% 
  summarise( ppt_year = sum(value) ) %>% 
  ungroup %>% 
  group_by( site, code ) %>% 
  summarise( ppt_avg_coll = mean(ppt_year) ) %>% 
  ungroup %>% 
  # Capitalize to merge with "demography" data frame later
  rename( Code = code )

# Climatic yearly averages
clim_mon_avg  <- clim_df %>% 
  group_by( site, lon, lat, coord_type, month ) %>% 
  summarise( ppt_month = mean(value) ) %>% 
  ungroup %>% 
  mutate( month = factor(month, levels = paste0(1:12)) )
  

# write down means
write.csv( clim_yr_avg,
           'C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/Monitoring/avg_yr_climate_sites.csv',
           row.names = F )
write.csv( clim_yr_avg_collection,
           'C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/Monitoring/avg_yr_climate_codes.csv',
           row.names = F )


# Plots ------------------------------------------------------------------------

# Average monthly precipitations
clim_mon_avg %>% 
  # subset( lon < -100 ) %>% 
  ggplot() +
  geom_line( aes(month,ppt_month,
                 group = site, 
                 color = lon) ) +
  scale_color_viridis_c() +
  theme_minimal() +
  theme( axis.title = element_text( size = 20),
         axis.text  = element_text( size = 15) ) +
  labs( x = 'Month',
        y = 'Monthly precipitation (mm)' ) +
  ggsave( 'results/monthly_precip.tiff',
          width = 6.3, height = 6.3, 
          compression = 'lzw')

# Precipitation by longitude
ggplot(clim_yr_avg) + 
  geom_point( aes( lon,ppt_avg, 
                   color = coord_type),
              size = 5, alpha = 0.5 ) +
  labs( x = 'Longitude',
        y = "Average precipitation (mm)") +
  ylim( 0, 1400 ) +
  theme_minimal() +
  theme( axis.title   = element_text( size = 20),
         axis.text    = element_text( size = 15),
         legend.text  = element_text( size = 15),
         legend.title = element_text( size = 15) ) +
  scale_color_colorblind() +
  guides( color = guide_legend( title="Coordinate type") ) +
  ggsave( 'results/precip_vs_longitude.tiff',
          width = 6.3, height = 6.3, 
          compression = 'lzw')
