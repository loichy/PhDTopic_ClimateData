#===============================================================================
# Description: Some basic weather and climate projections data in R
# author: HENRY Loic, loic.henry@dauphine.psl.eu
#===============================================================================

# Clean memory 
rm(list=ls())
gc()

# Load package
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, terra, maps, here, ncdf4, raster, climate, devtools, 
               sf, sp, rnaturalearth, Matrix)

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

#===============================================================================
# 1) Download data ------
#===============================================================================
# Use ERA5 data
# install.packages("ecmwfr")
library("ecmwfr")

# 1. Register to ECMWF to get an account and be able to request the API: https://www.ecmwf.int/

# 2. Get your API keys: https://cds.climate.copernicus.eu/how-to-api

# 3. Set a key to the keychain
wf_set_key(key = "97b8d178-a8e9-4fa1-ac36-3fae90a074e1") # Insert your key within the quotes

# you can retrieve the key using
wf_get_key()

# 4. Go to ERA5 daily stat data and select data of your choice: https://cds.climate.copernicus.eu/datasets/derived-era5-single-levels-daily-statistics?tab=download
# Accept terms and copy from the "API request" block
# 5. Paste in your script and transform in R language with ECMWF addins
request <- list(
  dataset_short_name = "derived-era5-single-levels-daily-statistics",
  product_type = "reanalysis",
  variable = "total_precipitation",
  year = "2024",
  month = "09",
  day = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"),
  daily_statistic = "daily_sum",
  time_zone = "utc+00:00",
  frequency = "1_hourly",
  area = c(52, -5, 42, 9),
  target = "TMPFILE.nc"
)
# Add ".nc" format at the target name file

# 5. Download using the ecmwfr::wf_request command
# TAKES AROUND 1 min for France and a month of daily data 
file <- ecmwfr::wf_request(
  request  = request,  # the request
  transfer = TRUE,     # download the file
  path     = here(dir$source)       # store data in current working directory
)

#===============================================================================
# 2) Import data ------
#===============================================================================
# Different potential packages to open ncdf data

# Open with ncdf package
nc_prec <- ncdf4::nc_open(file)

# Open ncdf file with terra package
terra_prec <- terra::rast(file)

# Or with raster package
raster <- raster::raster(file)
# and for for all layers:
raster_brick <- raster::brick(file)

#===============================================================================
# 3) Discover data ------
#===============================================================================
# information on the raster
nc_prec
terra_prec
raster
# dimensions: 41 * 57 with one grid being 0.25 degree (~20km)
# nlayers: 31 -> number of time layers, here 31 days
# names: these are the names of the layers (here : indices of the days)

# Plot the raster
plot(raster)
terra::plot(raster, main = "ERA-5 Reanalysis Demo (Precipitations)")

# Extract the precipitations values
prec <- ncdf4::ncvar_get(nc_prec, "tp")
dim(prec) 
# 3D array

# Can use the data frame vizualisation with 
G <- raster[] # For each cell: its corresponding temperature on that day

#===============================================================================
# 4) Aggregate data spatially ------
#===============================================================================

# 1. Load the sf object with the geographic units
# French department: https://github.com/gregoiredavid/france-geojson/blob/master/departements.geojson
dept_sf <- sf::st_read(here(dir$source,"departements.geojson"))
dept_sp <- as(dept_sf, "Spatial") # Transform as a spatial polygon object, for later aggregation
plot(dept_sp)

# Check if same CRS as the others
st_crs(dept_sp)
raster::crs(raster)
# Plot to check
png(here(dir$figures,"prec_dept.png"), width = 800, height = 600)
plot(raster[[1]], alpha=0.5)
plot(dept_sf$geometry, col = rgb(0, 0, 1, alpha = 0.25), add=T)
dev.off()

# 2. Load population raster, which will be useful for weighting grid cells
# For preparation of this rqster see RCode PreparePopulationRaster
pop_raster <- raster::raster(here(dir$source, "aggregated_pop_raster.tif"))

# 3. Crop precipitation raster to only keep cells touching French metropolitan boundaries
library(rnaturalearth)
library(sf)

# Get France boundary (country level)
france <- ne_countries(country = "France", scale = "medium", returnclass = "sf")
france_sp <- as(france, "Spatial") # Transform in a sp object
# Check CRS
st_crs(france_sp)
# Reduce the france limit extent to the one of the precipitation raster
cropped_france <- crop(france_sp, raster)
# Crop the precipitation raster
# cropped_raster <- mask(raster, cropped_france)
# cropped_raster_brick <- mask(raster_brick, cropped_france)

# 4. Compute population weights by department associated with each cell 
id <- raster::brick(raster[[1]], pop_raster)
info <- raster::extract(x = id, y = dept_sp, cellnumbers=T, weights=T) # return the values of the cells that are covered in each department, while here keeping all cells touching the borders (output also gives the share of the cell covered)
names(info) <- dept_sp@data$nom
info[[1]] # Cell: id of the cell; Value: population in the cell; Weight: share of the cell lying within the department

# Create weights, which sum to 1 within each department
# Note that we here consider weights corresponding to the share of the population in the department from that cell
tinfo <- lapply(names(info), function(i){
  print(i)
  # i = names(info)[10]
  df <- as.data.frame(info[[i]])
  df$Total.precipitation[is.na(df$Total.precipitation)] <- 0 # In case of missing weather index -> 0 (avoid operation which doesn't work)
  df$gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2020_2pt5_min[is.na(df$Total.precipitation)] <- 0 # In case of missing weather index -> 0 population so that this cell does not influence aggregation
  df$gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2020_2pt5_min[is.na(df$gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2020_2pt5_min)] <- 0 # In case of NA population -> replace with 0 to avoid operations that do not work
  df$w <- df$gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2020_2pt5_min * df$weight / sum(df$gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2020_2pt5_min * df$weight, na.rm = T) # make weight add to 1 in a given department
  df <- df[(!is.na(df$w)) | (!is.na(df$cell)),, drop=F] # drop=F, to make sure object stays as a matrix even with 1 row
  df$cell
  df$dept.order <- match(i,names(info))
  df
})
names(tinfo) <- dept_sp@data$nom

# Check weights construction
unique(sapply(tinfo, function(x) sum(x$w))) # wights add to 1? 
sapply(info, function(x) nrow(x)) #Number grid cell numbers inside department BEFORE adding weights
sapply(tinfo, function(x) nrow(x)) # Number of grid cell numbers inside department after adding weights

# Convert to dataframe that will go into our transformation matrix
dinfo <- do.call("rbind",tinfo)
dinfo <- dinfo[,c("cell","dept.order","weight","w")]
# Weight is the fraction of the cell which lies within the boundary of the departement
# W: is the population weighted fraction of the cell which lies within the boundary of the department

# 5. Create spatial aggregation matrix
# Projection "P" matrix to extract weather information by department

# We can use the dinfo object to create a matrix which extract weather information 
# for each department using matrix algebra. We need a transformation/projection
# matrix (named P) so that we can do the following:

# A <- P %*% G

# where:
# G is a matrix with all the climate data, obtained as G <- g[], where g is a raster stack
# G has a dimension n x t, where n is the number of grid cells in the climate data, and t is the number of time periods (or raster layers)
# A is a matrix with N rows (number of department) and t columns (number of layers in g)
# P is our weight/projection matrix we need to construct
# P dimensions? n rows (number of grid cells) by N columns (number of counties)

# Let's start with 3 raster layers (3 days of precipitations)
g <- stack(raster_brick[[1]],raster_brick[[2]],raster_brick[[3]]) # test with 3 layers (raster band, corresponding to 3 periods)
G <- g[]  
dim(G)  # should be number of grid cells in raster x 3 (number of layers)

# Create transformation matrix
# With population weights
P_weight <- Matrix::sparseMatrix(i = dinfo$cell,
                                 j = dinfo$dept.order,
                                 x = dinfo$w,
                                 dims = c(ncell(g), 
                                          nrow=length(unique(dinfo$dept.order))
                                          )
                                 )
# Use command sparsematrix from Matrix package as the computer will only hold the values and location of non-zero values in the memory (more efficient)
colnames(P_weight) <- dept_sp@data$nom

# Double check that all columns sum to 1
unique(colSums(P_weight))
colSums(P_weight)

# Test
A_weight <- t(P_weight) %*% G
dim(A_weight)
head(A_weight)
A_weight
dept_prec <- data.frame(as.matrix(A_weight))
dept_prec$nom <- row.names(dept_prec)

# Plot
# Transform as sf object the dataset by adding he geometry column of the department
dept_prec_sf <- dept_prec %>%
  left_join(dept_sf, by = "nom") %>% 
  st_as_sf()
# Transform the original raster as a dataframe
raster_df <- as.data.frame(g, xy = TRUE, na.rm = TRUE)

# Plot first the raster and then add the department
Plot1 <- ggplot() +
  # Add the raster with transparency (alpha controls transparency)
  geom_sf(data = dept_prec_sf, alpha=0) +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = X0), alpha = 0.75) 
  
ggsave(plot = Plot1, filename = here(dir$figures, "Fig1.png"), width = 8, height = 6, dpi = 300)


Plot2 <- ggplot() + 
  # Add the department level data of precipitations
  geom_sf(data = dept_prec_sf, aes(fill=X0), color = "black", alpha=0.75) 
ggsave(plot = Plot2, filename = here(dir$figures, "Fig2.png"), width = 8, height = 6, dpi = 300)
