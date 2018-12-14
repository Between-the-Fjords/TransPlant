##################
### CN_GONGGA  ###
##################

#### Import Community ####
ImportCommunity_CN_Gongga <- function(){
  fl <- list.files("R/community_CN_Gongga/", full.names = TRUE)
  sapply(fl, source)
  
  ## ---- load_community
  con <- src_sqlite(path = "data/CN_Gongga/transplant.sqlite", create = FALSE)
  # need to move all code to dplyr for consistancy
  
  #load cover data and metadata
  cover_thin_CN_Gongga <- load_comm(con = con)
  
  return(cover_thin_CN_Gongga)
}


#get taxonomy table
ImportTaxa_CN_Gongga <- function(){
  fl <- list.files("R/community_CN_Gongga/", full.names = TRUE)
  sapply(fl, source)
  
  ## ---- load_community
  con <- src_sqlite(path = "data/CN_Gongga/transplant.sqlite", create = FALSE)
  
  # need to move all code to dplyr for consistancy
  taxa_CN_Gongga <- tbl(con, "taxon") %>%
  collect()
return(taxa_CN_Gongga)
}


#### Cleaning Code ####
# Clean trait data
CleanTrait_CN_Gongga <- function(dat){
  dat2 <- dat %>% 
    filter(Project %in% c("LOCAL", "0", "C")) %>% 
    mutate(Treatment = plyr::mapvalues(Project, c("C", "0", "LOCAL"), c("C", "O", "Gradient"))) %>% 
    mutate(Taxon = trimws(Taxon)) %>% 
    mutate(Year = year(Date),
           Country = "CH",
           Gradient = as.character(1),
           Project = "T") %>% 
    rename(BlockID = Location) %>%
    mutate(PlotID = paste(BlockID, Treatment, sep = "-"),
           ID = paste(Site, Treatment, Taxon, Individual_number, Leaf_number, sep = "_")) %>% 
    select(Country, Year, Site, Gradient, BlockID, PlotID, Taxon, Wet_Mass_g, Dry_Mass_g, Leaf_Thickness_Ave_mm, Leaf_Area_cm2, SLA_cm2_g, LDMC, C_percent, N_percent , CN_ratio, dN15_percent, dC13_percent, P_AVG, P_Std_Dev, P_Co_Var) %>% 
    gather(key = Trait, value = Value, -Country, -Year, -Site, -Gradient, -BlockID, -PlotID, -Taxon) %>% 
    filter(!is.na(Value))
  
  return(dat2)
}

# Cleaning community data
CleanCommunity_CN_Gongga <- function(community_CN_Gongga_raw){
  dat2 <- community_CN_Gongga_raw %>% 
    filter(TTtreat != c("OTC")) %>% 
    rename(Year = year, Treatment = TTtreat, Cover = cover, SpeciesName = speciesName) %>% 
    mutate(Gradient = "CN_Gongga",
                      Country = as.character("China"),
           Treatment = recode(Treatment, "control" = "Control", "local" = "LocalControl", "warm1" = "Warm", "cool1" = "Cold", "warm3" = "WarmLong", "cool3" = "ColdLong")) %>% 
    mutate(SpeciesName = recode(SpeciesName, "Potentilla stenophylla var. emergens" = "Potentilla stenophylla")) %>% 
    filter(!is.na(Cover), !Cover == 0)
  
  return(dat2)
}

  
# Cleaning China meta community data
CleanMetaCommunity_CN_Gongga <- function(metaCommunity_CN_Gongga_raw){
  dat2 <- metaCommunity_CN_Gongga_raw %>% 
    select(PlotID, Year, Moss, Lichen2, Litter, BareGround, Rock, Vascular, Bryophyte, Lichen, MedianHeight_cm, MedianMossHeight_cm) %>% 
    mutate(Gradient = "CN_Gongga",
           Country = as.character("China"))
  return(dat2)
}

CleanMeta_CN_Gongga <- function(dat){
  dat2 <- dat %>% 
    mutate(Elevation = as.numeric(as.character(Elevation)),
           Gradient = "CN_Gongga",
           Country = as.character("China"),
           Site = as.character(Site),
           YearEstablished = 2012,
           PlotSize_m2 = 0.0625)
  
  return(dat2)
}



#### IMPORT, CLEAN AND MAKE LIST #### 
ImportClean_CN_Gongga <- function(){
  
  ### IMPORT DATA
  meta_CN_Gongga_raw = get(load(file = file_in("data/CN_Gongga/metaCN_Gongga.Rdata")))
  metaCommunity_CN_Gongga_raw = get(load(file = file_in("data/CN_Gongga/metaCommunityCN_Gongga_2012_2016.Rdata")))
  community_CN_Gongga_raw = ImportCommunity_CN_Gongga()
  taxa_CN_Gongga = ImportTaxa_CN_Gongga()
  trait_CN_Gongga_raw = get(load(file = file_in("data/CN_Gongga/traits_2015_2016_China.Rdata")))

  ### CLEAN DATA SETS
  ## CN_Gongga
  meta_CN_Gongga = CleanMeta_CN_Gongga(meta_CN_Gongga_raw)
  metaCommunity_CN_Gongga = CleanMetaCommunity_CN_Gongga(metaCommunity_CN_Gongga_raw)
  community_CN_Gongga = CleanCommunity_CN_Gongga(community_CN_Gongga_raw)
  taxa_CN_Gongga = ImportTaxa_CN_Gongga()
  trait_CN_Gongga = CleanTrait_CN_Gongga(trait_CN_Gongga_raw)
  
  # Make list
  CN_Gongga = list(meta = meta_CN_Gongga,
                   metaCommunity = metaCommunity_CN_Gongga,
                   community = community_CN_Gongga,
                   taxa = taxa_CN_Gongga,
                   trait = trait_CN_Gongga)
  
  return(CN_Gongga)
}
