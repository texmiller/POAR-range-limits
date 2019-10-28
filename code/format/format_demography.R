# Format demographic data
rm(list=ls())
setwd("D:/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography")
setwd("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography")
setwd("C:/Users/tm634/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography")
setwd("C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/Experiment/Demography")
options(stringsAsFactors=F)
library(dplyr)
library(testthat)

# read data ----------------------------------------------------------------
# 2014
f14     <- read.csv("fall2014/f2014DemoData.csv")
# 2015
s15     <- read.csv("spr15/s2015DemoData.csv")
# 2016
s16     <- read.csv("spr16/s2016DemoData.csv")
# 2017
s17     <- read.csv("spr17/s2017DemoData.csv")
# all collections in greenhouse (for debugging mistakes)
all_coll<- read.csv("D:/Dropbox/POAR--Aldo&Tom/Range limits/Genetics/allCollections.csv")
all_coll<- read.csv("C:/Users/tm9/Dropbox/POAR--Aldo&Tom/Range limits/Genetics/allCollections.csv")
all_coll<- read.csv("C:/Users/tm634/Dropbox/POAR--Aldo&Tom/Range limits/Genetics/allCollections.csv")
all_coll<- read.csv("C:/Users/ac22qawo/Dropbox/POAR--Aldo&Tom/Range limits/Genetics/allCollections.csv")

# format data for merge-----------------------------------------------------------------

# characters fields to "NA", and convert to "numeric" 
to_numeric <- function(x){
    
  vars      <- c("tillerN", "flowerN", paste0("flowerL",c(1:3)) )
  sub_x     <- select(x, vars)
  new_x     <- lapply(sub_x, 
                         function(x) {
                            id    <- grep("[a-z]",x, ignore.case = TRUE)
                            x[id] <- NA
                            x     <- as.numeric(x)
                         }
                      )
  x[,vars]  <- as.data.frame(new_x) 
  return(x)
  
}

# format data as referring to time t + 1 
t1_format <- function(x,col_id){
  
  nn          <- names(x)
  nn[col_id]  <- paste0(nn[col_id], "_t1") 
  nn          <- gsub("flowerL","flowerLength",nn) 
  nn          <- gsub("Tillering.","Tillering",nn) 
  out         <- setNames(x, nn)
 
  return(out)
  
}

# change "_t1" columns to "_t0"
t0_format <- function(x){
  nn    <- names(x)
  nn    <- gsub("_t1","_t0",nn) 
  out   <- setNames(x, nn)
  return(out)
}


# 2014 to 2015 transition --------------------------------------------------------------

# introduce NAs and convert to numeric 
s15 <- to_numeric(s15)
s16 <- to_numeric(s16)
s17 <- to_numeric(s17)

# Format site names to make data frames "mergeable"
site_n <- data.frame( site_raw = c("Amarillo","ElReno","PSU","Wichita","Woodward",
                                    "LLELA","WF","VERNON","LUBBOCK","OZONA","LLANO",
                                    "LOST PINES","KATY","SHSU"),
                      site_good = c("buffalolake","elreno","psu","ninnescah","woodward",
                                    "llela","wf","vernon","lubbock","ozona","llano",
                                    "lostP","Katy","shsu") )

# change site names
f14 <- f14 %>% mutate( site = site_n[match(site, site_n$site_raw),
                              "site_good"] )

#Column Names of f14 identical to those of s15
f14 <- setNames(f14, c("site","ID","Aluminum.Tag","Code","Sex",
                       "leafN_t0","tillerN_t0","bottomLeafN_t0") )

# 2015
s15     <- t1_format(s15,col_id=c(6:16)) 
              # trim white space in CODE names
              # mutate( Code = trimws(Code) )

# Tom's merging code --------------------------------------------------

# merge data and format 
# NOTE: MISMATCH HERE: FIND WHERE MISMATCH IS!
# Sites seem to match
setdiff(unique(f14$site),unique(s15$site))
setdiff(unique(s15$site),unique(f14$site))

#TOM: There are problems in s15 codes!
setdiff(unique(f14$Code),unique(s15$Code))
setdiff(unique(s15$Code),unique(f14$Code)) ##<-- these are the problem codes

## There is one ID in 2015 that is missing from 2014
setdiff(unique(f14$ID),unique(s15$ID))
setdiff(unique(s15$ID),unique(f14$ID)) ##<-- these are the problem IDs

#TOM: Fixing problems here in the script, not touching raw data files
#1. I think 'QL' should be 'QLP
code.problems15<-setdiff(unique(s15$Code),unique(f14$Code))
s15[s15$Code==code.problems15[1],]$Code<-"QLP"
#2. I think 'CO' should be 'COB'
s15[s15$Code==code.problems15[2],]$Code<-"COB"
#3. I think CWN should be CWM
s15[s15$Code==code.problems15[3],]$Code<-"CWM"
#4. I think SSL should be SSC (ID numbers consistent with SSC)
s15[s15$Code==code.problems15[4],]$Code<-"SSC"
#5. I think HAC should be HHC
s15[s15$Code==code.problems15[5],]$Code<-"HHC"
#6. I think CAR should be LAR, based on ID number
s15[s15$Code==code.problems15[6],]$Code<-"LAR"
#7. Here, HHC just has an extra space
s15[s15$Code==code.problems15[7],]$Code<-"HHC"
#8. Here, LAR just has an extra space
s15[s15$Code==code.problems15[8],]$Code<-"LAR"
#9. I think JJC should be HHC
s15[s15$Code==code.problems15[9],]$Code<-"HHC"
#10. COB has an extra space
s15[s15$Code==code.problems15[10],]$Code<-"COB"
#11. CWM has an extra space
s15[s15$Code==code.problems15[11],]$Code<-"CWM"
#12. I think QLF should be QLP
s15[s15$Code==code.problems15[12],]$Code<-"QLP"
#13. SLR has an extra space
s15[s15$Code==code.problems15[13],]$Code<-"SLR"

## Now fix problem ID (just one)
## This is a female LAR 396. I think it is actually LAR 296 (also female)
## Making this change
s15[s15$ID==setdiff(unique(s15$ID),unique(f14$ID)),]$ID<-296

s15_tom <- s15


# Aldo's reproduction of merge issues ------------------------------------------------------

# re-read data
s15          <- read.csv("spr15/s2015DemoData.csv") %>% 
                  to_numeric( ) %>% 
                  t1_format( col_id=c(6:16) ) %>% 
                  # remove mistakes originating from white spaces
                  mutate( Code = trimws(Code) )

# get code issues 
code_issues  <- setdiff(unique(s15$Code),unique(f14$Code))

# issues coming from codes
s15_issue_df <- s15 %>% 
                  subset( Code %in% code_issues ) %>% 
                  select( Code, ID ) %>% 
                  unique

# Test that there are not more than 1 "issue" per ID
s15_issue_df %>% 
  count( ID ) %>% 
  .$n %>% 
  unique %>% 
  expect_equal( 1 )
                  
# IDs associated with issues
issue_ids     <- s15_issue_df$ID

# corrected Codes (from f14)
correct_id_df <- f14 %>% 
                  subset( ID %in% issue_ids ) %>% 
                  select( Code, ID ) %>% 
                  rename( code_correct = Code) %>% 
                  unique

# update codes and wrong ID
s15_aldo <- s15 %>%
              # update codes
              left_join( correct_id_df ) %>% 
              # swap bad codes in 2015 with those from 2014
              mutate( Code = replace(Code,
                                     !is.na(code_correct),
                                     code_correct[!is.na(code_correct)]) ) %>% 
              select( -code_correct ) %>% 
              # replace one wrong ID
              mutate( ID = replace(ID, ID == 396, 296) )
  
# test: no difference in codes (from 2014 to 15)
unique(s15_aldo$Code) %>% 
  setdiff( unique(f14$Code) ) %>% 
  length %>% 
  expect_equal( 0 )

# test: no difference in IDs (from 2014 to 15)
unique(s15_aldo$ID) %>% 
  setdiff( unique(f14$ID) ) %>% 
  length %>% 
  expect_equal( 0 )

# test that Tom's and Aldo's files are the same
all.equal( s15_aldo %>% arrange( Block, Aluminum.Tag, ID, Code, Sex),
           s15_tom %>% arrange( Block, Aluminum.Tag, ID, Code, Sex) ) %>% 
  expect_true

# official s15 file!
s15 <- s15_tom

# Now merge, all problems fixed
d_14_15 <- merge(f14, s15)
d_14_15 <- mutate(d_14_15, year = 2015)

# SHOULD FAIL: mistakes create mismatches b/w "all=F" and "all=T" 
expect_false( identical( merge(f14, s15),
                         merge(f14, s15, all=T) 
                        )
              )


# 2015 to 2016 transition --------------------------------------------------------------

# reformat 2015
s15      <- t0_format(s15)
# reformat 2016
s16      <- t1_format(s16, c(6:14))

## Check for mismatch errors in site, ID, and code
setdiff(unique(s16$site),unique(s15$site))
setdiff(unique(s15$site),unique(s16$site))
setdiff(unique(s16$Code),unique(s15$Code))
setdiff(unique(s15$Code),unique(s16$Code))
setdiff(unique(s16$ID),unique(s15$ID)) ##<-- this 396 problems shows up again
setdiff(unique(s15$ID),unique(s16$ID))

## Fix the ID problem (396->296)
s16[s16$ID==setdiff(unique(s16$ID),unique(s15$ID)),]$ID<-296

# merge data and format 
d_15_16 <- merge(s15, s16)
d_15_16 <- mutate(d_15_16, year = 2016)

# make sure merge with "all=F" and "all=T" yield same n. of rows
expect_equal( merge(s15, s16) %>% nrow,
              merge(s15, s16, all=T) %>% nrow )


# 2016 to 2017 transition --------------------------------------------------------------

# reformat 2016
s16      <- t0_format(s16)
# reformat 2017
s17      <- t1_format(s17, c(6:14))

## Check for mismatch errors in site, ID, and code
setdiff(unique(s16$site),unique(s17$site)) ## Lubbock and WF not in 2017
                                            # That's OK because we drop these sites anyway
setdiff(unique(s17$site),unique(s16$site))
setdiff(unique(s16$Code),unique(s17$Code))
setdiff(unique(s17$Code),unique(s16$Code)) ##<-Lots of problems here
setdiff(unique(s16$ID),unique(s17$ID)) ##problems
setdiff(unique(s17$ID),unique(s16$ID)) ##problems

## Fix code problems. Same as above
code.problems17<-setdiff(unique(s17$Code),unique(s16$Code))
s17[s17$Code==code.problems17[1],]$Code<-"QLP"
s17[s17$Code==code.problems17[2],]$Code<-"COB"
s17[s17$Code==code.problems17[3],]$Code<-"CWM"
s17[s17$Code==code.problems17[4],]$Code<-"LAR"
s17[s17$Code==code.problems17[5],]$Code<-"SSC"
s17[s17$Code==code.problems17[6],]$Code<-"HHC"
s17[s17$Code==code.problems17[7],]$Code<-"HHC"
s17[s17$Code==code.problems17[8],]$Code<-"LAR"
s17[s17$Code==code.problems17[9],]$Code<-"HHC"
s17[s17$Code==code.problems17[10],]$Code<-"COB"
s17[s17$Code==code.problems17[11],]$Code<-"CWM"
s17[s17$Code==code.problems17[12],]$Code<-"QLP"
s17[s17$Code==code.problems17[13],]$Code<-"SLR"

## Fix ID problems
#1. This one is a plant from wf that has no data in 2017. 
# Not worrying about this one since wf gets dropped.
s16[s16$ID==setdiff(unique(s16$ID),unique(s17$ID)),]
#2. Fix the other ID problem
s17[s17$ID==setdiff(unique(s17$ID),unique(s16$ID))[1],]$ID<-296

# merge data and format (NOTE: )
d_16_17 <- merge(s16, s17)
d_16_17 <- mutate(d_16_17, year = 2017)

# "all=F" and "all=T" DO NOT yield same rows: 2017 misses lubbock and Witchita Falls data
expect_equal( anti_join(s16, s17) %>% 
                .$site %>% 
                unique, c('lubbock','wf') )


# Format final data frame -------------------------------------------------------

# combine data frames and calculate clone areas
d_all <- bind_rows( list(d_14_15, d_15_16, d_16_17) )

# claculate clone areas 
# When size >0, idenfity when width/length is == 0 
d_all <- d_all %>% mutate(sum_widths = MaxWidth_t0 + MaxLength_t0)
r_w   <- which(d_all$sum_widths > 0 & d_all$MaxWidth_t0 == 0 )
r_l   <- which(d_all$sum_widths > 0 & d_all$MaxLength_t0 == 0 )

# test: no overlap b/w r_W and r_l
expect_equal(length( intersect(r_w, r_l) ), 0) # expect no overlap b/w r_w and r_l 

# Introduce "1" whenever a clone area > 1, but width/length is == 0
# Do this because otherwise, clone area would be 0, even if it should not be.
d_all <- d_all %>%
          mutate( MaxWidth_t0   = replace(MaxWidth_t0,  r_w, 1)) %>%
          mutate( MaxLength_t0  = replace(MaxLength_t0, r_l, 1)) %>% 
          # Calculate the area of clones 
          mutate( clone_area_t0 = pi * MaxWidth_t0 * MaxLength_t0) %>%
          mutate( clone_area_t1 = pi * MaxWidth_t1 * MaxLength_t1)


# mistakes in the character fields

# trim white space in character columns
d_all$site <- trimws(d_all$site)
d_all$Code <- trimws(d_all$Code)
d_all$Sex  <- trimws(d_all$Sex)


# This bit of code does not make sense any more: The codes are now correct.
# "correct" contains 16 elements: the 8 "correct" ones, and the 8 mistakes shown in line 155 (code_issues)

# # Mistakes in the "Code" variable
# # Format site names to make data frames "mergeable"
# coll_codes <- data.frame( raw = unique(d_all$Code),
#                           correct = c( 'QLP','CWM','HHC','SSC','LAR','COB','SLR','LLELA',
#                                        'LAR','CWM','HHC','QLP','QLP','COB','SSC','SSC')
#                          )
# 
# # change site names
# d_all <- d_all %>% 
#           mutate(Code = coll_codes[match(Code, coll_codes$raw),
#                                    "correct"] )

# mistakes found in "Woodward/2017/collections_2017.R"

# check mistakes one by one: 5, 68, 309, 319, 377
# no mistake here
subset(d_all, ID == 68) %>% 
  select( year, ID, site, Code, Sex ) %>% 
  unique 

# Clear mistake: all ids from COB are three-digits
subset(d_all, ID == 5) %>% 
  select( year, ID, site, Code, Sex ) %>% 
  unique

# clear mistake: 377 is a Female ONLY at woodward. Likely wrong tag, or wrong coloured tag.
subset(d_all, ID == 377) %>% 
  select( year, ID, site, Code, Sex ) %>% 
  unique

# clear mistake: 309 is a male only in LLANO
# Also, this is LLELA only in elreno, in 2016 and 2017 (should be LAR)
subset(d_all, ID == 309) %>% 
  select( year, ID, site, Code, Sex ) %>% 
  unique

# clear mistake: 319 is QLP male only in LLANO (should be SLR)
subset(d_all, ID == 319) %>% 
  select( year, ID, site, Code, Sex ) %>% 
  unique


# Clear mistakes: IDs 5, 309, 319, 377
repl_orig_id <- function(coll_d, field_d, num){
  
  tmp <- subset(coll_d, ID == num)[,c("Origin","Sex")]
  tmp <- setNames(tmp, c("Code","Sex"))
  field_d[field_d$ID == num,c("Code","Sex")] <- tmp
  return(field_d) 

}
d_all <- repl_orig_id(all_coll, d_all, 5)
d_all <- repl_orig_id(all_coll, d_all, 68)
d_all <- repl_orig_id(all_coll, d_all, 309)
d_all <- repl_orig_id(all_coll, d_all, 319)
d_all <- repl_orig_id(all_coll, d_all, 377)


# final formatting for data --------------------------------------

# get lat/lon information
latlong <- read.csv("data/SiteLatLong.csv") %>% 
              # scale/center longitude
              mutate( long.scaled = scale(Longitude)[,1],
                      long.center = scale(Longitude,scale=F)[,1] )
poar    <- merge(d_all, latlong, by="site")

# Symbols for the sexes 
sex_symbol <- function(x){
  out <- x %>% 
    mutate(symb = factor(Sex, labels=c("21","24")) ) %>%
    mutate(symb = as.character(symb) ) %>%
    mutate(symb = as.integer(symb) ) %>%
    mutate(col = factor(Sex, labels=c("blue","red")) ) %>%
    mutate(col = as.character(col) )
  return(out)
}

# epigenetic changes
poar <- d_all %>% 
          # add lat/lon information
          merge( latlong, by="site" ) %>% 
          # add a unique block ID to fit as random effect
          mutate(unique.block = interaction(site,Block) ) %>%
          mutate(flow_t1 = as.numeric(flowerN_t1>0) ) %>% # add flowering (bernoulli)
          mutate(pc_fn = flowerN_t1 / tillerN_t1) %>%
          mutate(pc_fn = replace(pc_fn, pc_fn == Inf | is.nan(pc_fn), NA) ) %>%
          mutate(surv_t1 = as.numeric(tillerN_t1>0) ) %>% # add survival (bernoulli)
          mutate(delta.tiller = tillerN_t1 - tillerN_t0) %>% # change in size
          mutate(log_ratio = log(tillerN_t1 / tillerN_t0) ) %>% # log ratio
          mutate(log_ratio = replace(log_ratio, 
                                     log_ratio == -Inf | log_ratio == Inf | is.nan(log_ratio), 
                                     NA) ) %>%
          ## if clone area is zero, make it arbitrarily small
          mutate(clone_area_t0 = replace(clone_area_t0, clone_area_t0==0, 1) ) %>%
          mutate(clone_area_t1 = replace(clone_area_t1, clone_area_t1==0, 1) ) %>%
          # flag ploidy of provenance
          mutate(ploidy = "low") %>%
          mutate(ploidy = replace(ploidy, Code == "HHC" | Code == "LLELA", "high") ) %>%
          sex_symbol() %>% 
          # Remove resuscitated individuals
          subset( !(surv_t1 %in% 1 & tillerN_t0 %in% 0) )


# Write out data ----------------------------------------------------------------

# demographic data
write.csv(d_all, "C:/CODE/POAR-range-limits/data/f14_s17_data.csv", row.names = F)

# analysis data frame dropping three "bad" sites
write.csv(poar %>% 
            # drop the 3 "bad" sites (Lubbock, Wichita Falls, LLELA)
            # lubbock: only 7 individuals survived, only in 2015
            # Wichita Falls: only 7 individuals survive to 2015, 3 to 2016
            # LLELA: 22 indiv. surviva to 2015, 3 to 2016, 2 to 2017
            subset( !(site %in% c("llela", "lubbock", "wf")) ), 
          "POAR-range-limits/data/demography.csv", row.names = F)

# analysis data including three "bad" sites
write.csv(poar,"POAR-range-limits/data/demography_allsites.csv", row.names = F)
