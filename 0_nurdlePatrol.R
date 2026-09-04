# Script 0: Nurdle Patrol  -------------------------------------------------------------
# Purpose: reads, cleans, & pre-processes Nurdle Patrol data prior to synthesis with NOAA NCEI IBTrARC data
# Load Libraries
# Load Data
  # read full Gulf of Mexico Nurdle Patrol csv from NurdlePatrol.org selected by polygon bounding box around GOM
# Pre-process data
  # standardize date column
  # reorder columns
  # drop surveys that do not match survey protocol (i.e. surveys < 10min)
# Aggregate Proximal Surveys
  # 'proximal' = 800m (~ ave distance a person walks in 10min) NEED TO FIND A CITATION FOR THIS (9.1.25)
  # create a distance matrix on lon/lat with haversine equation from geosphere
  # Use DBSCAN to detect clusters of surveys within this proximity, retain sites outside of 800m clusters as independent sites
# Normalize cluster data
  # first set unique cluster ID to all clusters & independent sites so site = 0 MPs aren't aggregated
  # average nurdle standardized counts (weighted) for each cluster per day: total, new strandline, old strandline

# Load Libraries ---------------------------------------------------------------
install.packages('here')
pacman::p_load("lubridate", "tidyverse", "here", "janitor", "geosphere", "dbscan")
source("functions.R")

# janitor cleans colnames to all lower case, separated by '_'
# tidygeocoder reverse and forward engineers lat/lon or locations
# geosphere uses distm to get meters
# dbscan scans for proximity clusters

# Read Data -----------------------------------------------------------------
# read full Nurdle Patrol Gulf of Mexico csv from NurdlePatrol.org 
nurdle <- read.csv(file = here('NurdlePatrol','GulfofMexico','GOMNurdle_static_4.8.25.csv')) #date indicates access date


# Preprocess --------------------------------------------------------------
# 1) rename Patrol.Date to date 
nurdle <- rename(nurdle, 'date' = 'Patrol.Date')

# 2) clean column names with janitor
nurdle <- nurdle %>% clean_names()
#View(nurdle)

# 3) convert date to POSIXct
nurdle$date <- as.POSIXct(nurdle$date, format = "%Y-%m-%d")

# 4) drop nurdle collections that do not match protocol : surveys < 10min
nurdle <- nurdle %>%
  filter(collecting_time >= 10)

# 5) reorder so key columns are first
nurdle <- nurdle %>%
  relocate(id, date, longitude, latitude, standardized_amount)


# Cluster Proximal Surveys --------------------------------------------------
  # on average, person walks ~ 0.5mi in 10 min (survey time), which is ~ 800m 
  # use DBSCAN to create clusters w/in 800m of each other over all of time
  # any survey site > 800m of any single surgery site within the cluster will be detected as noise and not included in cluster

# 1) set coordinates object
coords <- nurdle %>% 
  select(longitude, latitude) 

# 2) using geosphere function distm to calculate the distance matrix between coords using spherical haversine method
dist_matrix <- distm(coords, fun = distHaversine)

# 3) set distance matrix to a dist object 
d <- as.dist(dist_matrix)

# 4) run DBSCAN w/ distance inclusion parameter (800m), cluster set to at least 2 points, 
db <- dbscan(d, eps = 800, minPts = 2)  #800m distance, cluster must contain at least two points (survey sites)

# 5) set nurdle clusters to clusters identified above, retain survey sites outside of clusters as 'noise', i.e. individual sites
clustered <- nurdle %>%
  mutate(site_id = as.integer(db$cluster)) %>% # name clusters with integers
  mutate(site_id = if_else(site_id == 0, 0L, site_id))  # noise stays as 0
#View(clustered)


# Aggregate Proximal Surveys-----------------------------------------
  # aggregate and weight (mean abundance) nurdle abundance for clusters/day while considering new and old strandlines
# 1) create unique site id so non-clusters (site_id = 0) aren't counted as same site across time & counts aren't aggregated
clustered <- clustered %>%
  mutate(
    site_id = if_else(site_id == 0,
                      paste0("lon_", round(longitude, 4), "_lat_", round(latitude, 4)),
                      paste0(site_id, "_cluster")))
#View(clustered)

# 2) anchor the center of each cluster id so it's location doesn't shift day to day depending on how many sites were surveyed in the cluster that day and where they were
centroids <- clustered %>%
  group_by(site_id) %>%
  summarise(
    cluster_lon = mean(longitude, na.rm = TRUE),
    cluster_lat = mean(latitude, na.rm = TRUE)
  )


# 3) weight abundance for new and old strand lines and a total count  (round to integer, digits = 0 10.5.25)
cluster_means <- clustered %>%
  group_by(site_id, date) %>%
  summarise(
    weighted_new_strandline = if (any(new_strand_line == "Yes")) {
      round(mean(standardized_amount[new_strand_line == "Yes"], na.rm = TRUE), 0)
    } else {
      0
    },
    weighted_old_strandline = if (any(new_strand_line == "No")) {
      round(mean(standardized_amount[new_strand_line == "No"], na.rm = TRUE), 0)
    } else {
      0
    },
    weighted_total_amount = round(mean(standardized_amount, na.rm = TRUE), 0),
    .groups = "drop"
  )

# join weighted amount and cluster centroids by site cluster
cluster_means <- cluster_means %>%
  left_join(centroids, by = "site_id")


# 4) join weighted means back to the original dataset
full_df <- clustered %>%
  left_join(cluster_means, by = c("site_id", "date"))
#View(full_df)

# 5) detect the max radius of each cluster in case needed
site_radii <- full_df %>%
  group_by(site_id, cluster_lat, cluster_lon) %>%
  summarise(
    cluster_max_radius_m = max(distHaversine(
      cbind(cluster_lon, cluster_lat),
      cbind(longitude, latitude)
    )),
    .groups = 'drop'
  )

# Joining, only bringing the radius column
full_df <- full_df %>%
  left_join(site_radii %>% select(site_id, cluster_max_radius_m), 
            by = "site_id")
#View(full_df)


# Write Files --------------------------------------------------------------
# clustered data that retains indiv. survey data
write.csv(full_df, here('output','GOM_nurdle.csv'), row.names = FALSE)





