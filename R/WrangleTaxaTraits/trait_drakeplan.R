#############################
### TRAIT DRAKE PLAN ########
#############################
#C. Chisholm, 13 April 2021

# Load libraries
library("drake")
library("tidyverse")
library("vegan")
library("readxl")
library("lubridate")
library("e1071")
library("DBI")
library("RSQLite")
library("visNetwork")
library("TR8")

# drake configurations
pkgconfig::set_config("drake::strings_in_dots" = "literals")

# trick
pn <- . %>% print(n = Inf)

# Source downstream trait scripts
source("R/WrangleTaxaTraits/taxa_drakeplan.R") #for TaxaPlan (which cleans species names)
source("R/WrangleTaxaTraits/try_traits.R") #for try traits
source("R/WrangleTaxaTraits/merge_traits.R") #to merge site and try traits

# Import TRY Data
ImportTRYDrakePlan <- drake_plan(

trytraits = load_wrangle_try(cleaned=cleanedspecies)

)

# Merge field-collected trait data

ImportSiteTraitDrakePlan <- drake_plan(
  sitetraits = merge_site_trait_data(sitetraitdata = tibble::lst(NO_Ulvhaugen, NO_Lavisdalen, NO_Gudmedalen, NO_Skjellingahaugen, 
                                             CH_Calanda, US_Colorado, CN_Gongga, FR_AlpeHuez, FR_Lautaret, IT_MatschMazia1, IT_MatschMazia2))
)

# Merge all trait data

 MergeTraitDrakePlan <- drake_plan(
   alltraits = merge_trait_data(traits = tibble::lst(trytraits, sitetraits))
)

 #Now there has been more rows added again! 5292 instead of 4599.
# CleanTraitDrakePlan <- drake_plan(
#   
#   #add cleaning script here, see error risk from tundra trait team as an example
# )

TraitPlan <- bind_rows(TaxaPlan, ImportTRYDrakePlan, ImportSiteTraitDrakePlan, MergeTraitDrakePlan)

conf <- drake_config(TraitPlan)
conf

make(TraitPlan)
# 
loadd(alltraits)

# Check all is good
#drake_failed()

# View dependency graph
vis_drake_graph(TraitPlan)

#Checked names! Everything w/o replacement ok
# Pulsatilla vernalis
# Cetraria islandica
# Polygonum viviparum
# Elyna myosuroides
# Vaccinium gaultherioides = Vaccinium uliginosum
# Sedum sp.
# Listera ovata = Neottia ovata
# Gentiana tenella = Gentianella tenella
# Euphrasia minima
# Nigritella nigra = Gymnadenia nigra
# Arenaria multicaulis
# Galium sp.
# Potentilla brauneana
# Anthoxantum alpinum = Anthoxanthum alpinum
# Hieracium lactucela = Pilosella lactucella
# Agrostis schraderiana = Agrostis agrostiflora
# Hieracium piloselloides
# Aster bellidiastrum
# Festuca pratense = Festuca pratensis