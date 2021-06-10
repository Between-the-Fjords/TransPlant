##################
### NO_NORWAY  ###
##################

source("R/ImportData/community_NO_Norway/loadCover.r")

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
  taxa <- tbl(con, "taxon") %>%
    collect() 
  taxa <- unique(taxa$speciesName)
  
  return(taxa)
}


#### Cleaning Code ####
# Cleaning NO community data
CleanCommunity_NO_Norway <- function(community_NO_Norway_raw){
  dat <- community_NO_Norway_raw %>% 
    #left_join(taxa_NO_Norway, by = "species") %>% 
    select(-temperature_level, -summerTemperature, -annualPrecipitation, -precipitation_level) %>% #, -notbad, -authority, -family, -comment) #%>% 
    rename(originSiteID = siteID, originBlockID = blockID, Treatment = TTtreat, Cover = cover, SpeciesName = species, Year = year, Collector = recorder) %>% 
    mutate(destPlotID = paste(destBlockID, originBlockID, turfID, sep='_')) %>%
    # only select control, local control, warm/down transplant
    filter(Treatment %in% c("TT1", "TT2")) %>% 
    mutate(Treatment = recode(Treatment, "TT1" = "LocalControl", "TT2" = "Warm")) %>% 
    filter(Year != 2016) %>% 
    mutate(UniqueID = paste(Year, originSiteID, destSiteID, destPlotID, sep='_')) %>% 
    mutate(destPlotID = as.character(destPlotID),
           destBlockID = if (exists('destBlockID', where = .)){ as.character(destBlockID)} else {NA})

  
  dat2 <- dat %>%  
    filter(!is.na(Cover)) %>%
    group_by_at(vars(-SpeciesName, -Cover)) %>%
    summarise(SpeciesName = "Other",Cover = pmax((100 - sum(Cover)), 0)) %>% 
    bind_rows(dat) %>% 
    mutate(Total_Cover = sum(Cover), Rel_Cover = Cover / Total_Cover)
  
  comm <- dat2 %>% filter(!SpeciesName %in% c('Other')) %>% 
    filter(Cover > 0) 
  cover <- dat2 %>% filter(SpeciesName %in% c('Other')) %>% 
    select(UniqueID, destSiteID, SpeciesName, Cover, Rel_Cover) %>% group_by(UniqueID, destSiteID, SpeciesName) %>% summarize(OtherCover=sum(Cover), Rel_OtherCover=sum(Rel_Cover)) %>%
    rename(CoverClass=SpeciesName)
  return(list(comm=comm, cover=cover))
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
    filter(!is.na(Value)) %>% rename(destSiteID=Site)
  
  return(dat2)
}


CleanMeta_NO_Norway <- function(meta_NO_Norway_raw){
  meta_NO_Norway <- meta_NO_Norway_raw %>% 
    mutate(Elevation = as.numeric(as.character(Elevation)),
           Country = "Norway",
           YearEstablished = 2009,
           YearMin = 2009,
           YearMax = 2017,
           YearRange = YearMax-YearEstablished,
           PlotSize_m2 = 0.0625,
           destSiteID = recode(Site, "Lav" = "Lavisdalen", "Hog" = "Hogsete", "Ulv" =  "Ulvhaugen", "Vik" = "Vikesland", "Gud" = "Gudmedalen", "Ram" = "Rambera", "Arh" = "Arhelleren", "Skj" = "Skjellingahaugen", "Ves" = "Veskre", "Alr" = "Alrust", "Ovs" = "Ovstedal", "Fau" = "Fauske"), 
          Gradient = case_when(destSiteID %in% c("Ulvhaugen", "Alrust", "Fauske")~ "NO_Ulvhaugen" ,
                                destSiteID %in% c("Lavisdalen", "Hogsete", "Vikesland")~ "NO_Lavisdalen" ,
                                destSiteID %in% c("Gudmedalen", "Rambera", "Arhelleren")~ "NO_Gudmedalen" ,  
                                destSiteID %in% c("Skjellingahaugen", "Veskre", "Ovstedal")~ "NO_Skjellingahaugen")) %>%
          #Longitude = as.numeric(recode(destSiteID, 'Lavisdalen'=7.2759600, 'Hogsete'=7.17666, 'Ulvhaugen'=8.1234300,'Vikesland'=7.16981999, 'Gudmedalen'=7.17560999, 'Rambera'=6.63028,'Arhelleren'=6.337379999, 'Skjellingahaugen'=6.4150400, 'Veskre'=6.5146800,'Alrust'=8.7046600, 'Ovstedal'=5.9648700, 'Fauske'=9.0787600)),
          #Latitude = as.numeric(recode(destSiteID, 'Lavisdalen'=60.8230999, 'Hogsete'=60.8759999, 'Ulvhaugen'=61.0242999,'Vikesland'=60.880299, 'Gudmedalen'=60.8327999, 'Rambera'=61.086599,'Arhelleren'=60.6651999, 'Skjellingahaugen'=60.9335000, 'Veskre'=60.5444999,'Alrust'=60.8203000, 'Ovstedal'=60.690100, 'Fauske'=61.0354999))) %>%
  select(Gradient, destSiteID, Longitude, Latitude, Elevation, YearEstablished, YearMin, YearMax, YearRange, PlotSize_m2, Country)

  return(meta_NO_Norway)
}


# # Cleaning meta community data
# CleanMetaCommunity_NO_Norway <- function(metaCommunity_NO_Norway_raw, g){
#   dat2 <- metaCommunity_NO_Norway_raw %>%
#     rename(destSiteID = SiteID) %>%
#     mutate(Gradient = paste("NO_Norway", g, sep='_'),
#            Country = as.character("Norway"))
#   return(dat2)
# }



#### IMPORT, CLEAN AND MAKE LIST #### 
ImportClean_NO_Norway <- function(g){
  sites <- case_when(g == 1 ~ c("Ulvhaugen", "Alrust", "Fauske"),
                     g == 2 ~ c("Lavisdalen", "Hogsete", "Vikesland"),
                     g == 3 ~ c("Gudmedalen", "Rambera", "Arhelleren"),
                     g == 4 ~ c("Skjellingahaugen", "Veskre", "Ovstedal")
  )
  
  ### IMPORT DATA
  meta_NO_Norway_raw = get(load(file = file_in("data/NO_Norway/meta_NO_Norway.Rdata")))
  metaCommunity_NO_Norway_raw = get(load(file = file_in("data/NO_Norway/metaCommunity_NO_Norway.Rdata")))
  
  #make database connection
  con <- src_sqlite(path = file_in("data/NO_Norway/seedclim.sqlite"), create = FALSE)
  community_NO_Norway_raw = ImportCommunity_NO_Norway(con)
  taxa_NO_Norway = ImportTaxa_NO_Norway(con)
  trait_NO_Norway_raw = read_csv(file = file_in("data/NO_Norway/traitdata_NO.csv"))
  
  ### CLEAN DATA SETS
  meta_NO_Norway = CleanMeta_NO_Norway(meta_NO_Norway_raw) %>% 
    filter(destSiteID %in% sites)
  #metaCommunity_NO_Norway = CleanMetaCommunity_NO_Norway(metaCommunity_NO_Norway_raw, g) %>% 
  #filter(destSiteID %in% sites)
  trait_NO_Norway = CleanTrait_NO_Norway(trait_NO_Norway_raw) %>% 
    filter(destSiteID %in% sites)
  
  cleaned_NO_Norway = CleanCommunity_NO_Norway(community_NO_Norway_raw) 
  community_NO_Norway = cleaned_NO_Norway$comm %>% filter(destSiteID %in% sites)
  cover_NO_Norway = cleaned_NO_Norway$cover %>% filter(destSiteID %in% sites)
  
  # Make list
  NO_Norway = list(meta = meta_NO_Norway,
                   cover = cover_NO_Norway, # this is e.g. total cover of PFGs.
                   community = community_NO_Norway,
                   taxa = taxa_NO_Norway,
                   trait = trait_NO_Norway)
  
  return(NO_Norway)
}

