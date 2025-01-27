#===============================================================================
# Description: Prepare population grid to match resolution and CRS of ERA5 grid
# author: HENRY Loic, loic.henry@dauphine.psl.eu
#===============================================================================
# List directories 
dir <- list()
dir$root <- here()
dir$source <- here(dir$root, "source") # Folder for original data files
dir$data <- here(dir$root, "data") # Folder for prepared/created data files
dir$code <- here(dir$root, "code") # Folder for code analyses
dir$figures <- here(dir$root, "figures") # Folder for created figures
dir$tables <- here(dir$root, "tables") # Folder for created tables

# Create non existing directories
lapply(dir, function(i) dir.create(i, recursive = T, showWarnings = F))

# Load the required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, terra, maps, here, ncdf4, raster, climate, devtools, sf)

# Download population grid at: # https://www.earthdata.nasa.gov/data/catalog/sedac-ciesin-sedac-gpwv4-apct-wpp-2015-r11-4.11

#Pop File name
pop_file <- "gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2020_2pt5_min.tif"

# Load pop file
pop <- terra::rast(here(dir$source, pop_file))

# Load precipitation raster
raster <- terra::rast(here(dir$source, "TMPFILE.nc"))

# Step 1: Reproject the finer raster to match the CRS of the coarser raster
pop_reprojected <- terra::project(pop, crs(raster))

# Step 2: Align the extent of the finer raster to the extent of the coarser raster
aligned_pop <- terra::crop(pop_reprojected, raster, snap="near")

# Step 3: Replace NaN values with 0
aligned_pop_cleaned <- ifel(is.nan(aligned_pop), 0, aligned_pop)

# Step 4: Aggregate the finer raster to match the resolution of the coarser raster
# Use the sum function for aggregation
aggregated_pop_raster <- aggregate(aligned_pop_cleaned, 
                                    fact=6, 
                                    fun=sum)

# Save the result to a new file if needed
writeRaster(aggregated_pop_raster, here(dir$source,"aggregated_pop_raster.tif"), overwrite=TRUE)

# Visualize 
plot(raster[[1]], main="Coarse Raster")
plot(aggregated_pop_raster, main="Aggregated Fine Raster", alpha=0.75)

