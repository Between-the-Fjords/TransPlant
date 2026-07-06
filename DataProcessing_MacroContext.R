#########################################
### Put together data for MC analyses ###
#########################################

# J. Lynn 11.2024

# if you get a "cache" error, it is due to a seed set somewhere in the making of the drake plan. run - clean(destroy = TRUE) 

#libraries
library(tidyverse); library(car); library(raster);library(ncdf4)

# start by loading in community data and taxa
#source("runsource_drakeplan.R")
#source("R/WrangleTaxaTraits/taxa_drakeplan.R")
load("~/Desktop/TransPlant/data/TransPlantNetwork_101024.RData")

# Start by setting up the data for analysis
## get start and end data 
startdat <- dat %>% dplyr::select(-UniqueID) %>% group_by(Region) %>% filter(Year == min(Year)) %>%  
  ungroup() %>% rename(year_start = Year, cover_start = Cover, rel_cov_start = Rel_Cover)
enddat <- dat %>% dplyr::select(-UniqueID) %>% group_by(Region) %>% filter(Year == max(Year)) %>% 
  ungroup() %>% rename(year_end = Year, cover_end = Cover, rel_cov_end = Rel_Cover)

comdat <- full_join(startdat, enddat, by = c("Region", "originSiteID", "originBlockID", "destSiteID",
                                             "destBlockID", "destPlotID", "Treatment", "turfID",
                                             "SpeciesName", "Elevation"))

# load in species climate distribution summaries
#sppclim <- read_csv("combinedSppClim.csv")
ori_clim <- read_csv("sppclimsum_TransPlant.csv")
ori_clim_minmax <- read_csv("sppclimsum_TransPlant_minmax.csv")
ext_clim <- read_csv("sppclimsum_TransPlant_extSpp_minmax.csv")
HP_clim <- read_csv("sppclimsum_TransPlant_HypPer_minmax.csv")

# add min max to ori_clim 
ori_clim <- ori_clim %>% dplyr::select(species, ahull, aoo_2)
ori_clim_minmax <- ori_clim_minmax %>% left_join(ori_clim, by = "species")
HP_clim <- filter(HP_clim, species == "Trommsdorffia maculata")

# combine the two datasets
ext_clim <- ext_clim %>% dplyr::select(-aoo) %>% mutate(species = case_when(species == "Plantago maritima"~ "Plantago maritima subsp. maritima",
                                                                            species == "Prunella vulgaris" ~ "Prunella vulgaris subsp. hispida",
                                                                            species == "Ranunculus acris" ~ "Ranunculus acris subsp. friesianus",
                                                                            species == "Ranunculus brotherusii" ~ "Ranunculus brotherusii var. tanguticus",
                                                                            species == "Solidago virgaurea" ~ "Solidago virgaurea subsp. minuta",
                                                                            species == "Thymus praecox" ~ "Thymus praecox subsp. polytrichus",
                                                                            TRUE ~ species))
climdat <- bind_rows(ori_clim_minmax, ext_clim, HP_clim) %>% unique()

sppclim <- climdat %>% filter(!is.na(mean_mat_full)) %>% filter(!is.na(species)) %>% unique() %>% 
  mutate(dup = duplicated(species))

# grab the stuff from cleaned species I need
cleanedspecies2 <- cleanedspecies %>% rename(species = matched_name) %>% 
  dplyr:: select(species, SpeciesName) %>% 
  mutate(SpeciesName = case_when(SpeciesName == "Capsella bursa pastoris" ~ "Capsella bursa-pastoris",
                                 TRUE ~ SpeciesName))

# combine then  slim down climate and species data
sppclim2 <- left_join(cleanedspecies2,sppclim, by = "species") %>% unique() %>% filter(!is.na(mean_mat_full)) %>%
  filter(!species == "Equisetum palustre") # one species was given multiple names

# combine the main dat set and species climate data
comdat <- left_join(comdat, sppclim2, by = "SpeciesName")

comdat_gspp <- comdat %>% filter(!is.na(mean_mat_sub))
comdat_bspp <- comdat %>% filter(is.na(mean_mat_sub)) %>% # use this to determine species that require extra gbif pulls
  filter(SpeciesName %in% c("Anemone alpina", "Anthoxantum alpinum", "Anthyllis vulneraria subsp. alpestris", 
                            "Anthyllis vuln. subsp. alpestris","Anthyllis vulneraria subsp. Alpestris",
                            "Capsella bursa pastoris", "Cerastium font. subsp. vulgare","Cerastium fontanum subsp. Vulgare", 
                            "Crocus albiflorus", "Deschampsia flexuosa","Elymus_trachycaulus", "Elyna myosuroides",
                            "Epilobium_angustifolium","Euphrasia rostkoviana", "Euphrasia rostkoviana aggr.",
                            "Euphrasia rostkoviana montana", "Festuca arundinacea","Festuca melanopsis",
                            "Festuca pratense", "Festuca pratensis", "Gentiana campestris", "Gentiana tenella",
                            "Gnaphalium supinum", "Helictotrichon alpinum", "Helictotrichon pubescens",
                            "Helictotrichon versicolor", "Hieracium glanduliferum", "Hieracium lactucella",
                            "Hieracium pilosella", "Hieratium pilosella", "Homalotrichon pubescens", "Humulus japonicus",
                            "Juncus trifidus", "juncus trifidus", "Leontodon helveticus", "Kobresia pygmae", 
                            "Leontodon helveticus", "Leontodon hisp. subsp. hispidus", "Leontodon his. subsp. Hispidus",
                            "Leontodon hispidus subsp. hispidus", "Leontodon hispidus subsp. hispidus",
                            "Ligusticum mutelina", "Ligusticum mutellina", "ligusticum mutellina", "Ligusticum mutellinoides",
                            "Lychnis flos-cuculi", "Osmorhiza_chilensis", "Oxalis acetocella", "Persicaria vivipara",
                            "Phleum alpinum subsp. rhaeticum", "Phleum rhaeticum", "Poa pra subsp. Alp", "polygonum viviparum", 
                            "Polygonum viviparum", "Potentilla bifruca", "Ran. acris subsp. Friesianus", "Senecio_corymbosus",
                            "Vaccinium_globulare", "Acinos alpinus", "Agrostis schraderiana", "Anemone demissa", "Anemone obtusiloba",
                            "Anemone rivularis", "Anthoxanthum alpinum", "Aquilegia caerulea", "Arenaria kansuensis","Arenaria roborowskii",
                            "Artemisia flaccida", "Aster bellidiastrum", "Astragalus membranaceus", "Athyrium subtriangulare",
                            "Castilleja sulpherea", "Coeloglossum viride", "Coluria longifolia", "Conyza canadensis", "Deyeuxia flavens",
                            "Elymus trachycaulus", "Gentiana alpestre", "Gentiana crenulatotruncata", "Gentiana farreri", "Glaux maratima",
                            "Helianthus anus", "Hie.pil", "Hieracium aurantiacum", "Hieracium_gracile", "Hypo.mac", "Hypochaeris maculata", 
                            "Jun.tri", "Kobresia capillifolia", "Kobresia cercostachys", "Kobresia humilis", "Kobresia pygmaea", 
                            "Kobresia robusta", "Laserpitium siler", "Ligusticum scapiforme", "Lloydia serotina", "Loiseleuria procumbens", 
                            "Mertensia fusiformis", "Microsteris_gracilis", "Minuartia sedoides", "Minuartia verna", "Morina kokonorica",
                            "Notopterygium forbesii", "Pedicularis mussotii var. lophocentra", "Perideridia_montana", "Polygonum sibiricum",
                            "Plantago maritima subsp. serpentina", "Poa crymophila", "Polygala chamaebuxus","Polygonum cyanandrum",
                            "Polygonum macrophyllum", "Potentilla anserina", "Potentilla bifurca", "Potentilla fruticosa", "Potentilla leuconata", 
                            "Ranunculus tuberosus","Prunella hispida", "Pyrethrum tatsienense", "Potentilla stenophylla", "Rhinanthus angustifolius",
                            "Rumex alpestris", "Salix brachista", "Saxifraga atrata", "Scirpus distigmaticum", "Sedum integrifolium",
                            "Senecio alpinus", "Sibbaldia pentaphylla", "Stellaria arenaria", "Stellaria uliginosa", "Stellaria umbellata", 
                            "Thalictrum aquilegifolium", "Thymus praecox polytrichus", "Trisetum spicatum", "Veronica collina", "Vicia angustifolia",
                            "Viola biflora var. rockiana", "Viola verecunda")) %>%
  mutate(species = case_when(SpeciesName == "Acinos alpinus"~ "Clinopodium alpinum",
                             SpeciesName == "Polygonum cyanandrum" ~ "Koenigia cyanandra", 
                             SpeciesName == "Polygonum macrophyllum" ~ "Bistorta macrophylla", 
                             SpeciesName == "Polygonum sibiricum" ~ "Knorringia sibirica", 
                             SpeciesName == "Potentilla anserina" ~ "Argentina anserina", 
                             SpeciesName == "Potentilla bifurca" ~ "Sibbaldianthe bifurca", 
                             SpeciesName == "Potentilla fruticosa" ~ "Dasiphora fruticosa",
                             SpeciesName == "Potentilla leuconata" ~ "Argentina leuconota", 
                             SpeciesName == "Potentilla stenophylla" ~ "Argentina stenophylla",
                             SpeciesName == "Prunella hispida" ~ "Prunella vulgaris subsp. hispida", 
                             SpeciesName == "Pyrethrum tatsienense" ~ "Tanacetum tatsienense", 
                             SpeciesName == "Ran. acris subsp. Friesianus" ~ "Ranunculus acris subsp. friesianus", 
                             SpeciesName == "Ranunculus tuberosus" ~ "Ranunculus montanus",
                             SpeciesName == "Rhinanthus angustifolius" ~ "Rhinanthus serotinus", 
                             SpeciesName == "Rumex alpestris" ~ "Rumex scutatus",  
                             SpeciesName == "Salix brachista" ~ "Salix lindleyana",  
                             SpeciesName == "Salix serpyllifolia" ~ "Salix serpillifolia", 
                             SpeciesName == "Sanguisorba minor" ~ "Sanguisorba minor",
                             SpeciesName == "Saxifraga atrata" ~  "Micranthes atrata", 
                             SpeciesName == "Scirpus distigmaticum" ~ "Trichophorum distigmaticum", 
                             SpeciesName == "Sedum integrifolium" ~ "Rhodiola integrifolia", 
                             SpeciesName == "Senecio alpinus" ~ "Jacobaea alpina",
                             SpeciesName == "Sibbaldia pentaphylla" ~ "Potentilla clandestina", 
                             #                    "Solidago virgaurea" ~ "Solidago virgaurea subsp. minuta", 
                             SpeciesName == "Stellaria arenaria" ~ "Stellaria arenarioides", 
                             SpeciesName == "Stellaria uliginosa" ~ "Stellaria alsine", 
                             SpeciesName == "Stellaria umbellata" ~ "Stellaria irrigua", 
                             SpeciesName == "Thalictrum aquilegifolium" ~ "Thalictrum aquilegiifolium", 
                             SpeciesName == "Thymus praecox polytrichus" ~ "Thymus praecox subsp. polytrichus", 
                             SpeciesName == "Trisetum spicatum" ~ "Koeleria spicata",
                             SpeciesName == "Veronica collina" ~ "Veronica triphyllos",
                             SpeciesName == "Vicia angustifolia" ~ "Vicia lathyroides", 
                             SpeciesName == "Viola biflora var. rockiana" ~ "Viola rockiana", 
                             SpeciesName == "Viola verecunda" ~ "Viola hamiltoniana",
                             SpeciesName == "Conyza canadensis" ~ "Erigeron canadensis",
                             SpeciesName == "Deyeuxia flavens" ~ "Calamagrostis tripilifera", 
                             SpeciesName == "Elymus trachycaulus" ~ "Elymus violaceus",
                             SpeciesName == "Gentiana crenulatotruncata" ~ "Gentiana prostrata var. crenulatotruncata", 
                             SpeciesName == "Gentiana farreri" ~ "Gentiana lawrencei",
                             SpeciesName == "Glaux maratima" ~ "Lysimachia maritima",
                             SpeciesName == "Helianthus anus" ~ "Helianthus annuus",
                             SpeciesName == "Hieracium aurantiacum" ~ "Pilosella aurantiaca",
                             SpeciesName == "Hieracium_gracile" ~ "Hieracium graciliusculum", 
                             SpeciesName == "Humulus japonicus" ~ "Humulus scandens", 
                             SpeciesName == "Hypochaeris maculata" ~ "Trommsdorffia maculata",
                             SpeciesName == "Hypo.mac" ~ "Trommsdorffia maculata",
                             SpeciesName == "Juncus trifidus" ~ "Oreojuncus trifidus",
                             SpeciesName == "Jun.tri" ~ "Oreojuncus trifidus",
                             SpeciesName == "Kobresia capillifolia" ~ "Carex capillifolia", 
                             SpeciesName == "Kobresia cercostachys" ~ "Carex cercostachys",
                             SpeciesName == "Kobresia humilis" ~ "Carex alatauensis",
                             SpeciesName == "Kobresia pygmae" ~ "Carex parvula",
                             SpeciesName == "Kobresia pygmaea" ~ "Carex parvula",
                             SpeciesName == "Kobresia robusta" ~ "Carex sargentiana", 
                             SpeciesName == "Laserpitium siler" ~ "Siler montanum",
                             SpeciesName == "Ligusticum scapiforme" ~ "Ligusticopsis scapiformis", 
                             SpeciesName == "Lloydia serotina"  ~ "Gagea serotina", 
                             SpeciesName == "Loiseleuria procumbens" ~ "Kalmia procumbens", 
                             SpeciesName == "Mertensia fusiformis" ~ "Mertensia brevistyla", 
                             SpeciesName == "Microsteris_gracilis" ~ "Phlox gracilis", 
                             SpeciesName == "Minuartia sedoides" ~ "Cherleria sedoides",
                             SpeciesName == "Minuartia verna" ~ "Sabulina verna", 
                             SpeciesName == "Morina kokonorica" ~ "Cryptothladia kokonorica", 
                             SpeciesName == "Notopterygium forbesii" ~ "Hansenia forbesii", 
                             SpeciesName == "Pedicularis mussotii var. lophocentra" ~ "Pedicularis mussotii", 
                             SpeciesName == "Perideridia_montana" ~ "Perideridia gairdneri", 
                             SpeciesName == "Plantago maritima subsp. serpentina" ~ "Plantago maritima", 
                             SpeciesName == "Poa crymophila" ~ "Poa alberti",
                             SpeciesName == "Polygala chamaebuxus" ~ "Polygaloides chamaebuxus", 
                             SpeciesName == "Coluria longifolia" ~ "Geum longifolium",
                             SpeciesName == "Coeloglossum viride" ~ "Dactylorhiza viridis",
                             SpeciesName == "Castilleja sulpherea" ~ "Castilleja septentrionalis", 
                             SpeciesName == "Astragalus membranaceus" ~ "Astragalus mongholicus", 
                             SpeciesName == "Athyrium subtriangulare" ~ "Athyrium spinulosum", 
                             SpeciesName == "Aster bellidiastrum" ~ "Bellidiastrum michelii", 
                             SpeciesName == "Arenaria roborowskii" ~ "Eremogone roborowskii", 
                             SpeciesName == "Aquilegia caerulea" ~ "Aquilegia coerulea",
                             SpeciesName == "Artemisia flaccida" ~ "Artemisia fulgens", 
                             SpeciesName == "Arenaria kansuensis" ~ "Eremogone kansuensis", 
                             SpeciesName == "Agrostis schraderiana" ~ "Agrostis schraderiana",
                             SpeciesName == "Anemone alpina"~ "Pulsatilla alpina",
                             SpeciesName == "Anemone demissa" ~ "Anemonastrum demissum", 
                             SpeciesName == "Anemone narcissiflora" ~ "Anemonastrum narcissiflorum",
                             SpeciesName == "Anemone obtusiloba" ~ "Anthoxanthum nipponicum",
                             SpeciesName == "Anthoxanthum alpinum" ~ "Anthoxanthum nipponicum", 
                             SpeciesName == "Anthoxantum alpinum" ~ "Anthoxanthum nipponicum", 
                             SpeciesName == "Anemone rivularis" ~ "Eriocapitella rivularis", 
                             SpeciesName == "Anemone rivularis" ~ "Eriocapitella rivularis", 
                             SpeciesName == "Anthyllis vulneraria subsp. alpestris" ~ "Anthyllis vulneraria",
                             SpeciesName == "Anthyllis vuln. subsp. alpestris" ~  "Anthyllis vulneraria",
                             SpeciesName == "Anthyllis vulneraria subsp. Alpestris" ~  "Anthyllis vulneraria",
                             SpeciesName == "Anthyllis vulneraria subsp. alpestris" ~  "Anthyllis vulneraria",
                             SpeciesName == "Bromopsis erecta" ~ "Bromus erectus",
                             SpeciesName == "Capsella bursa pastoris" ~ "Capsella bursa-pastoris",
                             SpeciesName == "Cerastium fontanum subsp. vulgare" ~ "Cerastium fontanum",
                             SpeciesName == "Cerastium fontanum subsp. Vulgare" ~ "Cerastium fontanum",
                             SpeciesName == "Cerastium font. subsp. vulgare" ~ "Cerastium fontanum",
                             SpeciesName == "Cerastium fontanumsp. triviale" ~ "Cerastium fontanum",
                             SpeciesName == "Crocus albiflorus" ~ "Crocus vernus",
                             SpeciesName == "Deschampsia flexuosa" ~ "Avenella flexuosa",
                             SpeciesName == "Elymus_trachycaulus" ~ "Elymus violaceus",
                             SpeciesName == "Elyna myosuroides" ~ "Carex myosuroides",
                             SpeciesName == "Epilobium_angustifolium" ~ "Chamaenerion angustifolium",
                             SpeciesName == "Euphrasia rostkoviana" ~ "Euphrasia officinalis",
                             SpeciesName == "Euphrasia rostkoviana aggr." ~ "Euphrasia officinalis",
                             SpeciesName == "Euphrasia rostkoviana montana" ~ "Euphrasia officinalis",
                             SpeciesName == "Festuca arundinacea" ~ "Festuca rubra",
                             SpeciesName == "Festuca melanopsis" ~ "Festuca violacea",
                             SpeciesName == "Festuca pratense" ~ "Festuca ovina",
                             SpeciesName == "Festuca pratensis" ~ "Festuca ovina",
                             SpeciesName == "Gentiana campestris" ~ "Gentianella campestris",
                             SpeciesName == "Gentiana alpestre" ~"Gentianella campestris",
                             SpeciesName == "Gentiana tenella" ~ "Comastoma tenellum",
                             SpeciesName == "Gnaphalium supinum" ~ "Omalotheca supina",
                             SpeciesName == "Helictotrichon alpinum" ~ "Helictochloa pratensis",
                             SpeciesName == "Helictotrichon pubescens" ~ "Avenula pubescens", 
                             SpeciesName == "Helictotrichon versicolor" ~ "Helictochloa versicolor",
                             SpeciesName == "Hieracium glanduliferum" ~ "Hieracium piliferum",
                             SpeciesName == "Hieracium lactucella" ~ "Pilosella lactucella",
                             SpeciesName == "Hieracium pilosella" ~ "Pilosella officinarum",
                             SpeciesName == "Hie.pil" ~ "Pilosella officinarum",
                             SpeciesName == "Hieratium pilosella" ~ "Pilosella officinarum",
                             SpeciesName == "Homalotrichon pubescens" ~ "Avenula pubescens", 
                             SpeciesName == "Juncus trifidus" ~ "Oreojuncus trifidus",
                             SpeciesName == "juncus trifidus" ~ "Oreojuncus trifidus",
                             SpeciesName == "Leontodon helveticus" ~ "Scorzoneroides helvetica",
                             SpeciesName == "Leontodon his. subsp. Hispidus" ~ "Leontodon hispidus",
                             SpeciesName == "Leontodon hisp. subsp. hispidus" ~ "Leontodon hispidus",
                             SpeciesName == "Leontodon hispidus subsp. hispidus" ~ "Leontodon hispidus",
                             SpeciesName == "Leontodon hispidus subsp. Hispidus" ~ "Leontodon hispidus",
                             SpeciesName == "Ligusticum mutelina" ~ "Mutellina adonidifolia", 
                             SpeciesName == "Ligusticum mutellina" ~ "Mutellina adonidifolia", 
                             SpeciesName == "ligusticum mutellina"~ "Mutellina adonidifolia",
                             SpeciesName == "Ligusticum mutellinoides" ~ "Pachypleurum mutellinoides",
                             SpeciesName == "Lychnis flos-cuculi" ~ "Silene flos-cuculi",
                             SpeciesName == "Osmorhiza_chilensis" ~ "Osmorhiza berteroi",
                             SpeciesName == "Oxalis acetocella" ~ "Oxalis acetosella",
                             SpeciesName == "Persicaria vivipara" ~ "Bistorta vivipara",
                             SpeciesName == "Phleum alpinum subsp. rhaeticum" ~ "Phleum alpinum",
                             SpeciesName == "Phleum rhaeticum" ~ "Phleum alpinum",
                             SpeciesName == "Poa pra subsp. Alp" ~ "Poa alpigena",
                             SpeciesName == "polygonum viviparum" ~ "Bistorta vivipara",
                             SpeciesName == "Polygonum viviparum" ~ "Bistorta vivipara",
                             SpeciesName == "Potentilla bifruca" ~ "Sibbaldianthe bifurca",
                             SpeciesName == "Ran. acris subsp. Friesianus" ~ "Ranunculus acris subsp. friesianus",
                             SpeciesName == "Senecio_corymbosus" ~ "Jacobaea incana",
                             SpeciesName == "Vaccinium_globulare" ~ "Vaccinium membranaceum")) %>%
  dplyr::select(Region:species)

# now add back in cliamte data and combine
correct_spp <- left_join(comdat_bspp, sppclim, by = "species") %>% unique()

comdat <- bind_rows(correct_spp, comdat_gspp) %>% unique()

# replace NAs with zeros
comdat$cover_start <- comdat$cover_start %>% replace_na(0)
comdat$cover_end <- comdat$cover_end %>% replace_na(0)
comdat$rel_cov_start <- comdat$rel_cov_start %>% replace_na(0)
comdat$rel_cov_end <- comdat$rel_cov_end %>% replace_na(0)

# function for cover change
diffunc <- function(init, post){log(post*100+1)-log(init*100+1)}
simpdif <- function(init, post){post*100-init*100}

comdat$dif_sub_rel <- diffunc(comdat$rel_cov_start, comdat$rel_cov_end)
comdat$simp_dif_rel <- simpdif(comdat$rel_cov_start, comdat$rel_cov_end)

ggplot(comdat, aes(x = dif_sub_rel))+ geom_histogram()+facet_wrap(~Region)

# replace NAs in year
comdat <- comdat %>% group_by(Region) %>%
  mutate(year_start = replace_na(year_start, min(year_start, na.rm = TRUE)),
         year_end = replace_na(year_end, max(year_end, na.rm = TRUE))) %>%
  ungroup()

comdat$exp_dur <- comdat$year_end-comdat$year_start

# determine colonzations and extinctions

comdat <- comdat %>% mutate(extinct = if_else(rel_cov_start != 0 & rel_cov_end == 0, 1, 0),
                             colonization = if_else(rel_cov_start == 0 & rel_cov_end != 0,1,0))

# add site lat and longitude
alldat <- tibble::lst( CH_Calanda,CH_Calanda2,CH_Lavey, CN_Damxung,CN_Gongga, CN_Heibei, DE_Grainau,
                       DE_Susalps,DE_TransAlps, FR_AlpeHuez, FR_Lautaret, IN_Kashmir,IT_MatschMazia1,
                       IT_MatschMazia2, NO_Gudmedalen, NO_Lavisdalen, NO_Skjellingahaugen, NO_Ulvhaugen,
                       SE_Abisko, US_Arizona, US_Colorado,US_Montana)

spatdat <- alldat %>% 
  map_df("meta", .id='Region') %>%
  ungroup() %>% dplyr::select(Region, destSiteID, Longitude, Latitude)

comdat <- left_join(comdat, spatdat, by = c("Region", "destSiteID"))

#### pull out climate data 
# access climate data here:
# https://chelsa-climate.org

mat <- raster("/Users/user/Desktop/CHELSA_climate_layers/CHELSA_bio1_1981-2010_V.2.1.tif") # mat
sum <- raster("/Users/user/Desktop/CHELSA_climate_layers/CHELSA_bio10_1981-2010_V.2.1.tif") # mean temp warmest quarter
ap <- raster("/Users/user/Desktop/CHELSA_climate_layers/CHELSA_bio12_1981-2010_V.2.1.tif") # AP
gsl <- raster("/Users/user/Desktop/CHELSA_climate_layers/CHELSA_gsl_1981-2010_V.2.1.tif") # growing season length
pet <- raster("/Users/user/Desktop/CHELSA_climate_layers/CHELSA_pet_penman_mean_1981-2010_V.2.1.tif") # PET

# pull the cliamte data
full_pts <- as.data.frame(cbind(comdat$Longitude,comdat$Latitude))
names(full_pts) <- c("lon", "lat")

matextfull30 <- raster::extract(mat, full_pts, method="simple", df=TRUE)
sumextfull30 <- raster::extract(sum, full_pts, method="simple", df=TRUE)
apextfull30 <- raster::extract(ap, full_pts, method="simple", df=TRUE)
gslextfull30 <- raster::extract(gsl, full_pts, method = "simple", df = TRUE)
petextfull30 <- raster::extract(pet, full_pts, method="simple", df=TRUE)

climdatfull <- left_join(matextfull30, sumextfull30,  by="ID") %>%
  left_join(., apextfull30, by="ID") %>% left_join(., petextfull30, by="ID") %>% 
  left_join(., gslextfull30, by = "ID")

# combine and rename
comdat <- comdat %>% bind_cols(., climdatfull) %>%
  dplyr::select(-c(ID)) %>%
  rename(mat_site=CHELSA_bio1_1981.2010_V.2.1) %>%
  rename(sum_temp_site=CHELSA_bio10_1981.2010_V.2.1) %>%
  rename(ap_site=CHELSA_bio12_1981.2010_V.2.1) %>% 
  rename(pet_site = CHELSA_pet_penman_mean_1981.2010_V.2.1)%>%
  rename(gsl_site = CHELSA_gsl_1981.2010_V.2.1)

# need to fix CHELSA units for all climate variables
comdat <- comdat %>% mutate(mean_mat_full = mean_mat_full*0.1 - 273.15, # makes units ºC
                            mean_ap_full = mean_ap_full *0.1, # makes units kg/m^2/year
                            mean_sum_full = mean_sum_full*0.1 - 273.15,# makes units ºC
                            mean_pet_full = mean_pet_full*0.01, # makes unites kg/m^2/month
                            median_mat_full = median_mat_full*0.1 - 273.15,
                            median_ap_full = median_ap_full*0.1,
                            median_sum_full = median_sum_full*0.1 - 273.15,
                            median_pet_full = median_pet_full*0.01,
                            min_mat_full = min_mat_full*0.1 - 273.15,
                            min5_mat_full = min5_mat_full*0.1 - 273.15,
                            max_mat_full = max_mat_full*0.1 - 273.15,
                            max95_mat_full = max95_mat_full*0.1 - 273.15,
                            min_ap_full = min_ap_full*0.1,
                            min5_ap_full = min5_ap_full*0.1,
                            max_ap_full = max_ap_full*0.1,
                            max95_ap_full = max95_ap_full*0.1,
                            min_sum_full = min_sum_full*0.1 - 273.15, 
                            min5_sum_full = min5_sum_full*0.1 - 273.15,
                            max_sum_full = max_sum_full*0.1 - 273.15,
                            max95_sum_full = max95_sum_full*0.1 - 273.15,
                            min_pet_full = min_pet_full*0.01,
                            min5_pet_full = min5_pet_full*0.01,
                            max_pet_full = max_pet_full*0.01, 
                            max95_pet_full = max95_pet_full*0.01,
                            mean_mat_sub = mean_mat_sub*0.1 - 273.15, # makes units ºC
                            mean_ap_sub = mean_ap_sub *0.1, # makes units kg/m^2/year
                            mean_sum_sub = mean_sum_sub*0.1 - 273.15,# makes units ºC
                            mean_pet_sub = mean_pet_sub*0.01, # makes unites kg/m^2/month
                            median_mat_sub = median_mat_sub*0.1 - 273.15,
                            median_ap_sub = median_ap_sub*0.1,
                            median_sum_sub = median_sum_sub*0.1 - 273.15,
                            median_pet_sub = median_pet_sub*0.01,
                            min_mat_sub = min_mat_sub*0.1 - 273.15,
                            min5_mat_sub = min5_mat_sub*0.1 - 273.15,
                            max_mat_sub = max_mat_sub*0.1 - 273.15,
                            max95_mat_sub = max95_mat_sub*0.1 - 273.15,
                            min_ap_sub = min_ap_sub*0.1,
                            min5_ap_sub = min5_ap_sub*0.1,
                            max_ap_sub = max_ap_sub*0.1,
                            max95_ap_sub = max95_ap_sub*0.1,
                            min_sum_sub = min_sum_sub*0.1 - 273.15, 
                            min5_sum_sub = min5_sum_sub*0.1 - 273.15,
                            max_sum_sub = max_sum_sub*0.1 - 273.15,
                            max95_sum_sub = max95_sum_sub*0.1 - 273.15,
                            min_pet_sub = min_pet_sub*0.01,
                            min5_pet_sub = min5_pet_sub*0.01,
                            max_pet_sub = max_pet_sub*0.01, 
                            max95_pet_sub = max95_pet_sub*0.01,
                            mat_site = mat_site*0.1 - 273.15,
                            sum_temp_site = sum_temp_site*0.1 - 273.15,
                            ap_site = ap_site*0.1,
                            pet_site = pet_site*0.01) 


# create climate difference metrics
comdat$mat_meandiff <- comdat$mean_mat_sub - comdat$mat_site
comdat$sum_meandiff <- comdat$mean_sum_sub - comdat$sum_temp_site
comdat$ap_meandiff <- (comdat$mean_ap_sub/1000) - (comdat$ap_site/1000)
comdat$gsl_meandiff <- comdat$mean_gsl_sub - comdat$gsl_site
comdat$pet_meandiff <- (comdat$mean_pet_sub) - (comdat$pet_site)

hist(comdat$pet_meandiff)

# save the output
write_csv(comdat, "MacroContextTransplantDat_2.csv")

# Checks on how climate data matches spp
sppcheck <- comdat %>% dplyr::select(Region, SpeciesName, species) %>% unique() ## fixed one species above
write_csv(sppcheck, "specieCheck.csv")

# check species without climate data
sppcheck2 <- comdat %>% dplyr::select(Region, SpeciesName, species, mean_mat_sub) %>%
  filter(is.na(mean_mat_sub)) %>% unique()
write_csv(sppcheck2, "speciesCheckNoClimdat.csv")
