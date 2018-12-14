#############################
### TRANSPLANT DRAKE PLAN ###
#############################

# Load library
library("drake")
library("tidyverse")
library("readxl")
library("lubridate")
library("e1071")

# drake configurations
pkgconfig::set_config("drake::strings_in_dots" = "literals")

# trick
pn <- . %>% print(n = Inf)

# source scripts
source("R/ImportCleanAndMakeList_CN_Gongga.R")
source("R/ImportCleanAndMakeList_NO_Norway.R")


# Import Data
ImportDrakePlan <- drake_plan(
  CN_Gongga = ImportClean_CN_Gongga(),
  
  NO_Ulvhaugen = ImportClean_NO_Norway(g = 1),
  NO_Lavisdalen = ImportClean_NO_Norway(g = 2),
  NO_Gudmedalen = ImportClean_NO_Norway(g = 3),
  NO_Skjellingahaugen = ImportClean_NO_Norway(g = 4)
)

# my_plan <- bind_rows(cwm_report, cwm_datasets_plan, cwm_analyses, setting_plan)
MyPlan <- ImportDrakePlan

conf <- drake_config(MyPlan)
make(MyPlan)
loadd()
vis_drake_graph(conf, targets_only = TRUE)

system("evince trait_distributions/Rmd/trait_ordinations.pdf", wait = FALSE)
