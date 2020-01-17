# Compare size distribution of observational and experimental POAR individuals
# Read observational (GPS-Polygon) data 
# Read in experimental data
library(dplyr)
library(ggridges)
library(ggplot2)
library(ggthemes)

in_dir <- "C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Data/"

# Polygon (gps) data---------------------------------------

files_12  <- list.files( paste0(in_dir, '2012/csv') )
files_13  <- list.files( paste0(in_dir, '2013/csv') )

# read polygon files
read_poly <- function(x, year){
  
  read.csv( paste0(in_dir,year,"/csv/",x),
            stringsAsFactors = F ) %>% 
    select(Sex, GNSS_Area) %>% 
    mutate( Population = strsplit(x, "_Pol.csv")[[1]][1] )
  
}

# read polygon data
poly_12 <- lapply(files_12, read_poly, 2012) %>% bind_rows
poly_13 <- lapply(files_13, read_poly, 2013) %>% bind_rows

# put polygon data together, and format it!
poly_df <- bind_rows(poly_12, poly_13) %>% 
            rename( area = GNSS_Area ) %>% 
            mutate( id       = as.character( 1:nrow(.) ),
                    log_area = log(area),
                    type     = 'Observational' ) %>% 
            select( -Population )


# Experimental data --------------------------------------

exp_df <- read.csv( 'data/demography.csv',
                    stringsAsFactors = F ) %>% 
            # only keep records that include the area of clones
            subset( (!is.na(clone_area_t0) | !is.na(clone_area_t1)) ) %>% 
            subset( !(clone_area_t1 %in% 1) ) %>% 
            # convert from cm2 to m2  
            mutate( clone_area_t1 = clone_area_t1 / 10000,
                    clone_area_t0 = clone_area_t0 / 10000 ) %>% 
            mutate( id = paste(site, ID, Aluminum.Tag, sep='_') ) %>% 
            select( id, Sex, clone_area_t1 ) %>% 
            rename( area     = clone_area_t1 ) %>% 
            mutate( log_area = log(area),
                    type     = 'Experimental' ) %>% 
            na.omit 
   
# final file for plots
area_df <- bind_rows( poly_df, exp_df )

# plot the density
ggplot(area_df, 
       aes(x = log_area, y = type ) ) + 
  # geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  # scale_fill_viridis_c() +
  geom_density_ridges( scale = 0.8 ) +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 0.85))) +
  theme_few() +
  # theme_ridges() + 
  labs( x = expression('Area: log(m'^2*')'),
        y = 'Type of data' ) +
  theme( axis.text  = element_text( size = 7.5 ),
         axis.title = element_text( size = 7.5 ) ) +
  ggsave( 'Manuscript/Figures/log_size_obs_experiment.pdf')

# with color
ggplot(area_df, 
       aes(x = log_area, y = type,
           fill = 0.5 - abs(0.5 - stat(ecdf))) ) +
  stat_density_ridges( geom = "density_ridges_gradient", 
                       calc_ecdf = TRUE,
                       scale = 0.95 ) +
  # geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_c( name = "Tail probability", direction = -1 ) +
  # geom_density_ridges( scale = 0.8 ) +
  scale_y_discrete(expand = expand_scale(mult = c(0.1, 1))) +
  theme_few() +
  labs( x = expression('Area: log(m'^2*')'),
        y = 'Type of data' ) +
  theme( axis.text   = element_text( size = 7.5 ),
         axis.title  = element_text( size = 7.5 ),
         legend.title = element_text( size = 7.5 ) ) +
  ggsave( 'Manuscript/Figures/log_size_obs_experiment_color.pdf')



area_df %>% 
  group_by( type ) %>% 
  summarise( mean_val_log   = mean(log_area),
             median_val_log = median(log_area) ) %>% 
  ungroup %>% 
  mutate( mean_val   = exp(mean_val_log),
          median_val = exp(median_val_log) )

# identify the problematic ids -------------------------------------------

probl_id <- exp_df %>% 
              subset( log_area < -8 ) %>% 
              .$id

# read 
prbl_df <- read.csv( 'data/demography.csv',
                    stringsAsFactors = F ) %>% 
            # only keep records that include the area of clones
            subset( (!is.na(clone_area_t0) | !is.na(clone_area_t1)) ) %>% 
            # convert from cm2 to m2  
            mutate( clone_area_t1 = clone_area_t1 / 10000,
                    clone_area_t0 = clone_area_t0 / 10000 ) %>% 
            mutate( id = paste(site, ID, Aluminum.Tag, sep='_') ) %>% 
            rename( area     = clone_area_t1 ) %>% 
            mutate( log_area = log(area) ) %>% 
            subset( id %in% probl_id )


prbl_df %>% head
