#### Run drake plan for Transplant project ####
## C. Chisholm, chelsea.chisholm@usys.ethz.ch

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

# Source drakeplan
r_make(source = "./R/TransPlant_DrakePlan.R")

# Load data 
# site data available as CO_Site (Country_Site, ex. CH_Lavey is Switzerland, Lavey)
# target dat is combined plant community abundance data for all sites (you can also write loadd(dat) to get just this target)
loadd()

# Check all is good
drake_failed()

# View dependency graph
r_vis_drake_graph(source = "./R/TransPlant_DrakePlan.R", targets_only = TRUE)

