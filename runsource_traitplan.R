#### Run drake plan for Transplant project: TRAITS ####
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
library("TR8")

# Source drakeplan
r_make(source = "./R/WrangleTaxaTraits/trait_drakeplan.R") #NOT WORKING, SITETRAITS IS FAILING FOR SOME REASON

# Load data 
# site data available as CO_Site (Country_Site, ex. CH_Lavey is Switzerland, Lavey)
loadd(alltraits)

# Check all is good
drake_failed()

# View dependency graph
r_vis_drake_graph(source = "trait_drakeplan.R", targets_only = TRUE)
