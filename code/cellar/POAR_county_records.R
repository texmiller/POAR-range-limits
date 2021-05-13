## I searched all POSR records on TORCH on 4/12/2020
## DOwnloaded these search results: http://portal.torcherbaria.org/portal/collections/listtabledisplay.php?db=all&taxa=Poa+arachnifera&usethes=1&taxontype=2
library(tidyverse)

torch <- read.csv("TORCH_POAR_records.csv")
levels(torch$stateProvince)

torch_counties <- torch %>% 
  select(stateProvince,county) %>% 
  subset(stateProvince %in% c("Texas","TEXAS","Oklahoma","Kansas")) %>% 
  unique() %>% 
  arrange(stateProvince,county) %>% 
  mutate(POAR = 1)

##are there any duplicated county names between states?
TX_OK <- which(torch_counties$county[torch_counties$stateProvince=="Texas"|torch_counties$stateProvince=="TEXAS"] %in% torch_counties$county[torch_counties$stateProvince=="Oklahoma"])
torch_counties$county[torch_counties$stateProvince=="Texas"|torch_counties$stateProvince=="TEXAS"][TX_OK]
## there is an Ellis county TX and Ellis county OK

TX_KS <- which(torch_counties$county[torch_counties$stateProvince=="Texas"|torch_counties$stateProvince=="TEXAS"] %in% torch_counties$county[torch_counties$stateProvince=="Kansas"])
## no TX-KS overlap

OK_KS <-which(torch_counties$county[torch_counties$stateProvince=="Oklahoma"] %in% torch_counties$county[torch_counties$stateProvince=="Kansas"])
torch_counties$county[torch_counties$stateProvince=="Oklahoma"][OK_KS]
## Comanche and Kiowa counties in OK and KS


write_csv(torch_counties,"POAR_county_records.csv")
## now go through the geodatabase feature by hand to make sure these counties have historical==yes in the attribute table

torch %>% filter(stateProvince=="Kansas") %>% select(county)
