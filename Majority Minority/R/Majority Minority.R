#Setup ####

options(scipen = 999)

library(tidycensus)
library(tidyverse)
library(skmeans)
library(stringr)
library(httr)
library(readr)
library(sf)
library(maptools)
library(GGally)
library(cluster)
library(factoextra)
library(ipumsr)

census_api_key("Your API Key goes here",install=TRUE,overwrite=TRUE)

#Historical Census data files obtained from nhgis.org
Directory <- "Data\\nhgis0001_csv\\"
Data_Files<-list.files(Directory)

Data_Files_DF<-data.frame(FileName = list.files(Directory)) %>% 
  mutate(Path = str_c(Directory, Data_Files),
         Geo = if_else(str_detect(Data_Files, "block"), "Block",
                       if_else(str_detect(Data_Files, "grp"), "Block Group", "Tract")
         )
  ) %>% 
  filter(FileName != "Codebooks")

Process_Files <- function (File_Path) {
  if (str_detect(File_Path, "1990")) {
  data.frame(read_csv(File_Path)) %>% 
    mutate(across(ET2001:ET2010, as.numeric)) %>% 
    mutate(Population = select(., ET2001:ET2010) %>% rowSums(na.rm = TRUE),
           Proportion = ET2001/Population) %>% 
    filter(Population>0) %>% 
    select(GISJOIN, YEAR, Non_Hispanic_White = ET2001, Population, Proportion)
  } else if (str_detect(File_Path, "2010")) {
    data.frame(read_csv(File_Path))  %>% 
      select(GISJOIN, YEAR, H7Z001, H7Z003) %>% 
      mutate(across(H7Z001:H7Z003, as.numeric)) %>% 
      mutate(Proportion = H7Z003/H7Z001) %>% 
      filter(H7Z001>0) %>% 
      select(GISJOIN, YEAR, Non_Hispanic_White = H7Z003, Population = H7Z001, Proportion)
  } else if (str_detect(File_Path, "2000") & (str_detect(File_Path, "tract")==FALSE)) {
    data.frame(read_csv(File_Path)) %>% 
      mutate(across(FYF001:FYF014, as.numeric)) %>% 
      mutate(Population = select(., FYF001:FYF014) %>% rowSums(na.rm = TRUE),
             Proportion = FYF001/Population) %>% 
      filter(Population>0) %>% 
      select(GISJOIN, YEAR, Non_Hispanic_White = FYF001, Population, Proportion)
  } else if (str_detect(File_Path, "2000") & (str_detect(File_Path, "tract"))) {
    data.frame(read_csv(File_Path)) %>% 
      mutate(across(FMS001:FMS014, as.numeric)) %>% 
      mutate(Population = select(., FMS001:FMS014) %>% rowSums(na.rm = TRUE),
             Proportion = FMS001/Population) %>% 
      filter(Population>0) %>% 
      select(GISJOIN, YEAR, Non_Hispanic_White = FMS001, Population, Proportion)
  }
}

Files <- map2(Data_Files_DF$Path, Data_Files_DF$Geo, 
              ~Process_Files(.x) %>% 
                mutate(Geo = .y))

ACS_Grab <- function(Vars, Geogr, Years = 2019) {
  
  National <- c("state", "place", "county", "metropolitan statistical area/micropolitan statistical area", "zcta")
  By_State <- c("tract", "block group", "block")
  
  if(Geogr %in% National) {
    lapply(Years, function(Yrs) {
      get_acs(year = Yrs, geography = Geogr, variables = Vars, geometry = FALSE) %>% 
        mutate(Year = Yrs,
               High_Error = ((moe/1.645)/estimate)>.2)
    } 
    )}
  else if (Geogr %in% By_State) {
    lapply(Years, function(Yrs) {
      map_df(states, ~get_acs(year = Yrs, geography = Geogr, state = .x, variables = Vars, geometry = FALSE)) %>% 
        mutate(Year = Yrs,
               High_Error = ((moe/1.645)/estimate)>.2)}
    )
  }
}

Race_Vars <- data.frame(VariableName = c("Population", "Non-Hispanic White"),
                        variable = c("B01001_001", "B01001H_001"))

Tract_2019<- map_df(Race_Vars$variable, 
                        ~ACS_Grab(.x, Geogr = "tract",))

Tract_2019_Errors <- Tract_2019 %>% 
  filter(variable == "B01001_001",
         High_Error == TRUE | is.na(estimate) == TRUE) %>% 
  select(GEOID)

Tract_2019_Fix <- Tract_2019 %>% 
  anti_join(Tract_2019_Errors) %>% 
  select(-NAME) %>% 
  inner_join(Race_Vars) %>% 
  select(-High_Error, -variable, -moe) %>% 
  pivot_wider(names_from = VariableName, values_from = estimate) %>% 
  mutate(Proportion = `Non-Hispanic White`/Population,
         YEAR = as.character(Year),
         Geo = "Tract") %>% 
  na.omit() %>% 
  select(GISJOIN = GEOID, YEAR, Non_Hispanic_White = `Non-Hispanic White`, Population, Proportion, Geo)

All_Data<-do.call(bind_rows, Files) %>% 
  bind_rows(Tract_2019_Fix) %>% 
  mutate(Proportion_5 = if_else(Proportion == 1, .95, floor(Proportion*20)/20))

Summary_Data<-All_Data %>% 
  group_by(YEAR, Proportion_5, Geo) %>% 
  summarize(Units = n(),
            Population = sum(Population),
            Non_Hispanic_White = sum(Non_Hispanic_White))

write_csv(Summary_Data, "Data\\Summary Data.csv")

ggplot(All_Data %>% 
         filter(Geo == "Tract", YEAR == 2019), 
       aes(x = Proportion_5)) +
  geom_histogram(position = 'identity', bins = 20) +
  labs(x = "Proportion White", y = "Number of Census Tracts") +
  ggtitle("Tracts by % of White Residents") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggsave(filename = "Tracts_2019.png", device = 'png', dpi = 600)

ggplot(All_Data %>% 
         filter(Geo == "Tract"), 
       aes(x = Proportion_5)) +
  geom_histogram(position = 'identity', bins = 20) +
  labs(x = "Proportion White", y = "Number of Census Tracts") +
  facet_grid(~YEAR) +
  ggtitle("Tracts by % of White Residents") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(0.5, 1)) +
  ggsave(filename = "Tracts_All.png", device = 'png', dpi = 600)

ggplot(All_Data %>% 
         filter(Geo == "Block Group"), 
       aes(x = Proportion_5)) +
  geom_histogram(position = 'identity', bins = 20) +
  labs(x = "Proportion White", y = "Number of Census Block Groups") +
  facet_grid(~YEAR) +
  ggtitle("Block Groups by % of White Residents") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(0.5, 1)) +
  ggsave(filename = "Block_Groups_All.png", device = 'png', dpi = 600)

ggplot(All_Data %>% 
         filter(Geo == "Block"), 
       aes(x = Proportion_5)) +
  geom_histogram(position = 'identity', bins = 20) +
  labs(x = "Proportion White", y = "Number of Census Blocks") +
  facet_grid(~YEAR) +
  ggtitle("Blocks by % of White Residents") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(0.5, 1)) +
  ggsave(filename = "Blocks_All.png", device = 'png', dpi = 600)

#How many white PEOPLE are living in more diverse tracts now?

All_Whites <- All_Data %>% 
  group_by(YEAR, Geo) %>% 
  summarize(Total_NH_Whites = sum(Non_Hispanic_White))

White_Diverse <- All_Data %>% 
  group_by(Proportion_5, YEAR, Geo) %>% 
  summarize(Count = n(),
            White = sum(Non_Hispanic_White)) %>% 
  inner_join(All_Whites) %>% 
  mutate(Percent_of_Whites = White/Total_NH_Whites)

ggplot(White_Diverse %>% filter(Geo == "Tract", YEAR == "2019"), aes(x = Proportion_5, y = White)) +
  geom_col() +
  labs(x = "Proportion White", y = "Number of White Residents")  +
  ggtitle("White Residents by Tract % White") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggsave(filename = "Whites_All_Tract_2019.png", device = 'png', dpi = 600)

ggplot(White_Diverse %>% filter(Geo == "Tract"), aes(x = Proportion_5, y = White, xend = 1.05)) +
  geom_col() +
  labs(x = "Proportion White", y = "Number of White Residents") +
  facet_grid(cols = vars(YEAR)) +
  ggtitle("White Residents by Tract % White") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(0.5, 1)) +
  ggsave(filename = "Whites_All_Tract.png", device = 'png', dpi = 600)

ggplot(White_Diverse %>% filter(Geo == "Block Group"), aes(x = Proportion_5, y = White, xend = 1.03)) +
  geom_bar(stat = 'identity', position = 'identity') +
  labs(x = "Proportion White", y = "Number of White Residents") +
  facet_grid(cols = vars(YEAR)) +
  ggtitle("White Residents by Block Group % White") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(0.5, 1)) +
  ggsave(filename = "Whites_All_Block_Group.png", device = 'png', dpi = 600)

ggplot(White_Diverse %>% filter(Geo == "Block"), aes(x = Proportion_5, y = White, xend = 1.03)) +
  geom_bar(stat = 'identity', position = 'identity') +
  facet_grid(cols = vars(YEAR)) +
  labs(x = "Proportion White", y = "Number of White Residents") +
  ggtitle("White Residents by Block % White") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(0.5, 1)) +
  ggsave(filename = "Whites_All_Block.png", device = 'png', dpi = 600)

ggplot(White_Diverse %>% filter(YEAR != 2019), aes(x = Proportion_5, y = White, xend = 1.05)) +
  geom_bar(stat = 'identity', position = 'identity') +
  facet_grid(rows = vars(Geo), cols = vars(YEAR)) +
  labs(x = "Proportion White", y = "Number of White Residents") +
  ggtitle("White Residents by Geography % White") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(0.5, 1)) +
  ggsave(filename = "Whites_All.png", device = 'png', dpi = 600)