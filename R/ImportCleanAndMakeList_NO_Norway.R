##################
### NO_NORWAY  ###
##################

source("R/community_NO_Norway/loadCover.r")

#### Import Community ####
ImportCommunity_NO_Norway <- function(con){
  ## ---- load_community
  
  #load cover data and metadata
  #cover
  cover_thin_NO_Norway <- load_cover_NO_Norway(con = con)
  
  return(cover_thin_NO_Norway)
}


#get taxonomy table
ImportTaxa_NO_Norway <- function(con){
  
  # need to move all code to dplyr for consistancy
  #get taxonomy table
  taxa_NO_Norway <- tbl(con, "taxon") %>%
    collect()
  
  return(taxa_NO_Norway)
}


#### Cleaning Code ####
# Cleaning China meta community data
CleanCommunity_NO_Norway <- function(community_NO_Norway_raw, taxa_NO_Norway){
  dat2 <- community_NO_Norway_raw %>% 
    left_join(taxa_NO_Norway, by = "species") %>% 
    select(-temperature_level, -summerTemperature, -annualPrecipitation, -precipitation_level, -notbad, -authority, -family, -comment) %>% 
    rename(originSiteID = siteID, originBlockID = blockID, Treatment = TTtreat, Cover = cover, SpeciesShort = species, Year = year, SpeciesName = speciesName, Collector = recorder) %>% 
    # only select control, local control, warm/down transplant
    filter(Treatment %in% c("TTC", "TT1", "TT2")) %>% 
    mutate(Treatment = recode(Treatment, "TTC" = "Control", "TT1" = "LocalControl", "TTC" = "Warm"))
  
  return(dat2)
}


# Clean trait data
CleanTrait_NO_Norway <- function(trait_NO_Norway_raw){
  dat2 <- trait_NO_Norway_raw %>% 
    rename(SpeciesName = Taxon, Collector = Data_collected_by, Leaf_Thickness_Ave_mm = Leaf_Thicness_Ave_mm, PlantID = ID) %>%
    mutate(SpeciesName = trimws(SpeciesName),
           Year = year(Date),
           Country = "Norway",
           Site = recode(Site, "Lav" = "Lavisdalen", "Hog" = "Hogsete", "Ulv" =  "Ulvhaugen", "Vik" = "Vikesland", "Gud" = "Gudmedalen", "Ram" = "Rambera", "Arh" = "Arhelleren", "Skj" = "Skjellingahaugen", "Ves" = "Veskre", "Alr" = "Alrust", "Ovs" = "Ovstedal", "Fau" = "Fauske")) %>% 
    select(-X1, -Date, -Longitude, -Latitude, -Elevation, -Project, -Gradient) %>%
    gather(key = Trait, value = Value, -Country, -Year, -Site, -Individual_number, -SpeciesName, -PlantID, -Collector) %>% 
    filter(!is.na(Value))
  
  return(dat2)
}


CleanMeta_NO_Norway <- function(meta_NO_Norway_raw){
  dat2 <- meta_NO_Norway_raw %>% 
    mutate(Elevation = as.numeric(as.character(Elevation)),
           Gradient = "NO_Norway",
           Country = as.character("Norway"),
           YearEstablished = 2009,
           PlotSize_m2 = 0.0625,
           Site = recode(Site, "Lav" = "Lavisdalen", "Hog" = "Hogsete", "Ulv" =  "Ulvhaugen", "Vik" = "Vikesland", "Gud" = "Gudmedalen", "Ram" = "Rambera", "Arh" = "Arhelleren", "Skj" = "Skjellingahaugen", "Ves" = "Veskre", "Alr" = "Alrust", "Ovs" = "Ovstedal", "Fau" = "Fauske"))
  
  return(dat2)
}


# Cleaning meta community data
CleanMetaCommunity_NO_Norway <- function(metaCommunity_NO_Norway_raw, g){
  dat2 <- metaCommunity_NO_Norway_raw %>% 
    mutate(Gradient = paste0("NO_Norway", g),
           Country = as.character("Norway"))
  return(dat2)
}



#### IMPORT, CLEAN AND MAKE LIST #### 
ImportClean_NO_Norway <- function(g){
  sites <- case_when(g == 1 ~ c("Ulvhaugen", "Alrust", "Fauske"),
                     g == 2 ~ c("Lavisdalen", "Hogsete", "Vikesland"),
                     g == 3 ~ c("Gudmedalen", "Rambera", "Arhelleren"),
                     g == 4 ~ c("Skjellingahaugen", "Veskre", "Ovstedal")
  )
  
  ### IMPORT DATA
  meta_NO_Norway_raw = get(load(file = file_in("data/NO_Norway/meta_NO_Norway.Rdata")))
  #metaCommunity_CN_Gongga_raw = get(load(file = file_in("data/NO_Norway/metaCommunity_NO_Norway.Rdata")))
  
  #make database connection
  con <- src_sqlite(path = "data/NO_Norway/seedclim.sqlite", create = FALSE)
  community_NO_Norway_raw = ImportCommunity_NO_Norway(con)
  taxa_NO_Norway = ImportTaxa_NO_Norway(con)
  trait_NO_Norway_raw = read_csv(file = file_in("data/NO_Norway/traitdata_NO.csv"))
  
  ### CLEAN DATA SETS
  meta_NO_Norway = CleanMeta_NO_Norway(meta_NO_Norway_raw) %>% 
    filter(Site %in% sites)
  #metaCommunity_NO_Norway = CleanMetaCommunity_NO_Norway(metaCommunity_NO_Norway_raw) %>% 
  #filter(Site %in% sites)
  community_NO_Norway = CleanCommunity_NO_Norway(community_NO_Norway_raw, taxa_NO_Norway) %>% 
    filter(destSiteID %in% sites)
  trait_NO_Norway = CleanTrait_NO_Norway(trait_NO_Norway_raw) %>% 
    filter(Site %in% sites)
  
  # Make list
  NO_Norway = list(meta = meta_NO_Norway,
                   #metaCommunity = metaCommunity_NO_Norway,
                   community = community_NO_Norway,
                   taxa = taxa_NO_Norway,
                   trait = trait_NO_Norway)
  
  return(NO_Norway)
}

