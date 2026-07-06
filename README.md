# TransPlant Network Project

This repository contains code for cleaning and loading data from a distributed network of plant community transplants along elevational gradients. 

To obtain the combined dataframe across all sites for plant community data, run the runsource_drakeplan.R script. This will run the TransPlant_DrakePlan.R in the R subfolder, which in turn runs all the scripts in the R/ImportData subfolder and merges using the merge_community.R script. The resulting dataframe 'dat' will contain plant community abundance for each site, elevation, time (year), and plot.

To find the traits dataset, run the runsource_traitplan.R. This will source the trait drakeplan in the R/WrangleTaxaTraits subfolder. In this drake plan, sitetraits and trytraits dataframes are loaded, where site traits are field-collected data from the sites with average traits per species per elevation per site (if collected), and the try traits are an average value per species that has been gapfilled from TRY. The final combined trait object alltraits contains the field-collected sites per elevation per species and then gap-fills any values where we did not have field-collected data using the try traits.

## Authors

* **Chelsea Chisholm** - chelsea.chisholm@gmail.com
* **Dagmar Egelkraut** - Dagmar.Egelkraut@uib.no



