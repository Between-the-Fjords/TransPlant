####################
### FR_Lautaret  ###
####################

#### Import Community ####
ImportCommunity_FR_Lautaret <- function(){
  community_FR_Lautaret_raw <- read.csv("data/FR_Lautaret/FR_Lautaret_commdata/pinpoints_cleaned_reduced.csv", sep=",")
  return(community_FR_Lautaret_raw)
} 

ImportCommunity_FR_Lautaret2 <- function(){
  community_FR_Lautaret_raw2 <- read.csv("data/FR_Lautaret/FR_Lautaret_commdata/transalp.csv", sep=";")
  return(community_FR_Lautaret_raw2)
} 

#### Cleaning Code ####
# Cleaning Lautaret community data
CleanCommunity_FR_Lautaret <- function(community_FR_Lautaret_raw, community_FR_Lautaret_raw2){
  dat <- community_FR_Lautaret_raw %>% 
    mutate(Replicate = ifelse(nchar(Replicate)==1, paste0("0", Replicate), Replicate))%>%
    filter(Subplot %in% c('0','B')) %>% #Choosing 0 subplot (controls) and B subplot (treatment) as it had no added individuals into it
    mutate(plotID = paste0(Site_Treatment, "_", Replicate), 
    Cover = number_of_obs,
    SpeciesName = ifelse(is.na(word(Species, 2)), paste0(Species, " sp."), Species))%>%
    separate(Site_Treatment, c("destSiteID", "Treatment"), "_") %>% #called CP or TP (control or transplant)
    mutate(originSiteID = case_when(destSiteID == "L" & Treatment == 'TP' ~ "G", 
                                    destSiteID == "L" & Treatment == 'CP' ~ "L",
                                    destSiteID == "G" & Treatment == 'CP' ~ "G",
                                    destSiteID == "G" & Treatment == 'TP' ~ "L")) %>% 
    mutate(Treatment = case_when(destSiteID == "L" & Treatment == 'TP' ~ "Warm", 
                                    destSiteID == "L" & Treatment == 'CP' ~ "LocalControl",
                                    destSiteID == "G" & Treatment == 'CP' ~ "LocalControl",
                                    destSiteID == "G" & Treatment == 'TP' ~ "Cold"))%>%
    select(Year, destSiteID, originSiteID, plotID, Treatment, SpeciesName, Cover) %>%
    mutate(UniqueID = paste(Year, originSiteID, destSiteID, plotID, sep='_')) %>% 
    mutate(destPlotID = paste(originSiteID, destSiteID, plotID, sep='_')) %>% 
    select(-plotID) %>% 
    mutate(destPlotID = as.character(destPlotID), destBlockID = if (exists('destBlockID', where = .)) as.character(destBlockID) else NA) %>%
    distinct() %>% #one duplicated row in original dataframe
    group_by(Year, originSiteID, destSiteID, destBlockID, destPlotID, UniqueID, Treatment, SpeciesName) %>%
    summarize(Cover = sum(Cover, na.rm=T)) %>% #had one species which occured twice in a plot, summing across
    mutate(SpeciesName = ifelse(SpeciesName == "Undetermined sp.", "Undetermined", SpeciesName))%>%
    ungroup()
  
  dat2 <- community_FR_Lautaret_raw2 %>% 
    mutate(Year = 2022,
           subplot = sapply(code_plot, function(x) sub(".*_(.*)$", "\\1", x)),
           destSiteID = sapply(code_plot, function(x) sub(".*_(High|Low)_.*", "\\1", x)),
           rep  = sapply(code_plot, function(x) sub(".*_(\\d{2})_.*", "\\1", x)))%>%
    rename(SpeciesName = lb_nom, Cover = prct_recouvrement)%>%
    filter(subplot %in% c("CTRL","B"))%>%
    mutate(destSiteID = recode(destSiteID, 
                               High = "G",
                               Low = "L"))%>%
    mutate(originSiteID = case_when(destSiteID == "L" & subplot == 'CTRL' ~ "L", 
                                    destSiteID == "G" & subplot == 'CTRL' ~ "G",
                                    destSiteID == "L" & subplot == 'B' ~ "G",
                                    destSiteID == "G" & subplot == 'B' ~ "L"),
           plotID = case_when(destSiteID == "L" & subplot == 'CTRL' ~ "L_L_L_CP", 
                              destSiteID == "G" & subplot == 'CTRL' ~ "G_G_G_CP",
                              destSiteID == "L" & subplot == 'B' ~ "G_L_L_TP",
                              destSiteID == "G" & subplot == 'B' ~ "L_G_G_TP"),
           Treatment = case_when(destSiteID == "L" & subplot == 'CTRL' ~ "LocalControl", 
                                 destSiteID == "G" & subplot == 'CTRL' ~ "LocalControl",
                                 destSiteID == "L" & subplot == 'B' ~ "Warm",
                                 destSiteID == "G" & subplot == 'B' ~ "Cold"),
           SpeciesName = recode(SpeciesName, 
                                'Carex sempervirens subsp. sempervirens' = "Carex sempervirens",
                                'Pilosella officinarum' = "Pilosella",
                                'Patzkea paniculata subsp. paniculata' = "Patzkea paniculata"))%>%
    mutate(SpeciesName = ifelse(is.na(word(SpeciesName, 2)), paste0(SpeciesName, " sp."), SpeciesName))%>%
    filter(!is.na(Treatment))%>%
    mutate(destPlotID = paste0(plotID, "_",rep),
           UniqueID = paste0(Year, "_", destPlotID),
           destBlockID = NA)%>%
    select(Year, originSiteID, destSiteID, destBlockID, destPlotID, UniqueID, Treatment, SpeciesName, Cover) %>%
    distinct() %>% #one duplicated row in original dataframe
    group_by(Year, originSiteID, destSiteID, destBlockID, destPlotID, UniqueID, Treatment, SpeciesName) %>%
    summarize(Cover = sum(Cover, na.rm=T)) %>% #had one species which occured twice in a plot, summing across
    ungroup()
  
  dat <- bind_rows(dat, dat2)
  dat2 <- dat %>%  
    filter(!is.na(Cover)) %>% #no Nas, just a precaution
    group_by_at(vars(-SpeciesName, -Cover)) %>%
    summarise(SpeciesName = "Other",Cover = pmax((100 - sum(Cover)), 0)) %>% #All total cover >100, pmax rebases this to zero
    bind_rows(dat) %>% 
    mutate(Total_Cover = sum(Cover), Rel_Cover = Cover / Total_Cover)
  
  comm <- dat2 %>% filter(!SpeciesName %in% c('Other')) %>%
    filter(Cover > 0)  
  cover <- dat2 %>% filter(SpeciesName %in% c('Other')) %>% 
    select(UniqueID, SpeciesName, Cover, Rel_Cover) %>% group_by(UniqueID, SpeciesName) %>% summarize(OtherCover=sum(Cover), Rel_OtherCover=sum(Rel_Cover)) %>%
    rename(CoverClass=SpeciesName)
  return(list(comm=comm, cover=cover)) 

}

# Clean taxa list (add these to end of above)
CleanTaxa_FR_Lautaret <- function(community_FR_Lautaret){
taxa <- unique(community_FR_Lautaret$SpeciesName)
  return(taxa)
}

# Clean metadata
CleanMeta_FR_Lautaret <- function(community_FR_Lautaret){
  dat <- community_FR_Lautaret %>%
    select(destSiteID, Year) %>%
    group_by(destSiteID) %>%
    summarize(YearMin = min(Year), YearMax = max(Year)) %>%
    mutate(Elevation = as.numeric(recode(destSiteID, 'G' = 2450, 'L' = 1950)),
           Gradient = 'FR_Lautaret',
           Country = 'France',
           Longitude = as.numeric(recode(destSiteID, 'G' = 6.40048, 'L' = 6.4190699)),
           Latitude = as.numeric(recode(destSiteID, 'G' = 45.0543600, 'L' = 45.04006)),
           YearEstablished = 2017,
           PlotSize_m2 = 1) %>% 
    mutate(YearRange = (YearMax-YearEstablished)) %>% 
    select(Gradient, destSiteID, Longitude, Latitude, Elevation, YearEstablished, YearMin, YearMax, YearRange, PlotSize_m2, Country) 
  
  return(dat)
}

#Clean trait data
CleanTrait_FR_Lautaret <- function(trait_FR_Lautaret_raw){
  dat2 <- trait_FR_Lautaret_raw %>%
    rename(SpeciesName = species, Individual_number = rep, destPlotID = plot, Trait = trait, destSiteID = climate, Treatment = treatment, Value=value) %>% 
    mutate(Country = "France",
           Trait = recode(Trait, 'F_Mass' = 'W_Mass_g',  'D_Mass' = 'Dry_Mass_g', 'L_Area' = 'Leaf_Area_cm2', 'H_Repr' = 'Plant_Rep_Height_cm',  'H_Veg' = 'Plant_Veg_Height_cm'),
           Treatment = recode(Treatment, "CP" = "LocalControl", "TP" = "Warm")) %>%
    dplyr::select(Country, destSiteID, destPlotID, Treatment, SpeciesName, Individual_number, Trait, Value) %>%
    mutate(Individual_number = as.character(Individual_number), Value = as.numeric(Value)) %>%
    filter(!is.na(Value), !is.na(SpeciesName))
  return(dat2)
}

#### IMPORT, CLEAN AND MAKE LIST #### 
ImportClean_FR_Lautaret <- function(){
  
  ### IMPORT DATA
  community_FR_Lautaret_raw = ImportCommunity_FR_Lautaret()
  community_FR_Lautaret_raw2 = ImportCommunity_FR_Lautaret2()
  trait_FR_Lautaret_raw = read.table("./data/FR_Lautaret/FR_lautaret_traitdata/TransPlant_Lautaret_traits_2018_30062020.txt")
  
  ### CLEAN DATA SETS
  cleaned_FR_Lautaret = CleanCommunity_FR_Lautaret(community_FR_Lautaret_raw, community_FR_Lautaret_raw2)
  community_FR_Lautaret = cleaned_FR_Lautaret$comm
  cover_FR_Lautaret = cleaned_FR_Lautaret$cover
  meta_FR_Lautaret = CleanMeta_FR_Lautaret(community_FR_Lautaret) 
  taxa_FR_Lautaret = CleanTaxa_FR_Lautaret(community_FR_Lautaret)
  trait_FR_Lautaret = CleanTrait_FR_Lautaret(trait_FR_Lautaret_raw)
  
  
  # Make list
  FR_Lautaret = list(meta = meta_FR_Lautaret,
                   community = community_FR_Lautaret,
                   cover = cover_FR_Lautaret,
                   taxa = taxa_FR_Lautaret,
                   trait = trait_FR_Lautaret)
  
  return(FR_Lautaret)
}

