# Script 1: Temporally & Spatially Link  -------------------------------------------------------------
  # Purpose: Synthesize pre-processed Nurdle Patrol and IBTrACS storm data sets
# Load data
# Pre-process: format & rename cols
# Temporal Link
  # Cross join data sets
  # Convert storm effect sizes to kilometers (max radial wind, quadrant winds, ROCI)
  # Wide filter: flag site storm-influenced if surveyed +/- 28d from storm start & end
  # Create sf objects from km data
  # Find date a storm was closest to a survey site & call that day == 0 of the storm in relation to that site
  # Link on 2 temporal scales
    # a) Fine phase - Flag whether a survey occurred pre|during|post a day of a storm (centered around day 0 closest storm track)  
    # b) Max phase: Flag whether a survey occurred pre|during|post the full storm start and end (complete_storm_phase)
  # Flag pre- and post-storm surveys based on temporal range (acute = 0-4d, subacute = 5-10, extended = 11-28)
  # Flag surveys as storm_influenced or control & rejoin data
# Temporal Clean
  # create linked_id (siteID_stormName)
  # ID site_storm combos that meet temoral criteria
  # TEMPORAL FILTER: drop rows that do not meet temporal criteria while retaining control data for sites that remain
  # Drop unnecessary columns
  # Rename & reorder columns for clarity
# Spatial Link
  # Detect distance survey is from temporally linked data
  # Detect proximal surveys: if survey is w/in ROCI either pre or during a storm = is proximal
  # Retain all observations (at all storm phases & control) for sites that are proximal to at least 1 linked storm
  # Retain control data for these sites
# Get Site Names
# Reduce Resolution: storms obs ~6/day, surveys 1/day, this step reduces storm obs to max values one a single day
  # clean errors (missing storm data)
# save master df

# Load Libraries ---------------------------------------------------------------
pacman::p_load('tidyverse','fuzzyjoin', 'janitor', 'geosphere','tidygeocoder', 'sf','purrr','stringr', 'here')
options(sf_use_s2 = TRUE) # uses spherical calculations
source(here("functions.R"))


# Call in Data ------------------------------------------------------------

# full clustered nurdle df w/ individual survey data
nurdle <- read.csv(here('output','GOM_nurdle.csv'))
#View(nurdle)

storms <- read.csv(here('output','GOM_storms.csv'))
#View(storms)


# Pre-processing ----------------------------------------------------------
# 1) rename date cols & set data type
nurdle <- nurdle %>%
  rename(date_nurdle = date)

storms <- storms %>%
  rename(date_storm = date)

storms$date_storm <- as.Date(storms$date_storm)
storms$full_storm_end <- as.Date(storms$full_storm_end)
storms$full_storm_start <- as.Date(storms$full_storm_start)
nurdle$date_nurdle <- as.Date(nurdle$date_nurdle)

# look for column name shares that will create a column collision during the join
names(nurdle)[names(nurdle) %in% names(storms)]

# 2) rename cols to distinguish b/w storm & survey data
storms <- storms %>%
  rename(longitude_storm = longitude, latitude_storm = latitude)

nurdle <- nurdle %>%
  rename(longitude_survey = longitude, latitude_survey = latitude)

# 3) remove unnamed storm (no wind metrics, etc)
storms <- storms %>%
  filter(name != "UNNAMED")


# Temporal Link -----------------------------------------------------------
# 1) cross join both data sets
df <- nurdle %>%
  cross_join(storms) # cross_join is part of dplyr

# detect storm year & month as column
df <- df %>%
  group_by(sid) %>%
  mutate(storm_year = year(date_storm),
         storm_month = month(date_storm)
  ) %>%
  ungroup()


# create a meters based quadrant wind conversion (needed for sf objects in meters -used for temporal and spatial links)
df <- df %>%
  mutate(
    r34_ne_m = usa_r34_ne * 1852, # r34 = winds reaching / exceeding 34knots (>34 = tropical storm, 20-34 = tropical depression)
    r34_se_m = usa_r34_se * 1852,
    r34_sw_m = usa_r34_sw * 1852,
    r34_nw_m = usa_r34_nw * 1852,
    r50_ne_m = usa_r50_ne * 1852, # r50 = winds reaching / exceeding 50 kts (tropical storm)
    r50_se_m = usa_r50_se * 1852,
    r50_sw_m = usa_r50_sw * 1852,
    r50_nw_m = usa_r50_nw * 1852,
    r64_ne_m = usa_r64_ne * 1852, # r64 = winds reaching / exceeding 64 kts (Cat. 1 hurricane and beyond)
    r64_se_m = usa_r64_se * 1852,
    r64_sw_m = usa_r64_sw * 1852,
    r64_nw_m = usa_r64_nw * 1852,
    roci_m = usa_roci * 1852,
    rmw_m = usa_rmw * 1852
  )


# 2) Filter down to storm-influenced (+/- 28d form storm) data. 
df <- df %>%
  group_by(sid, id) %>%
  filter(date_nurdle >= full_storm_start - days(28),
         date_nurdle <= full_storm_end + days(28)) %>%
         ungroup()


# 3) link on 2 temporal scales

# a) Fine phase - Flag whether a survey occurred pre|during|post a day of a storm when that storm was closest to that survey site (truest effect)    
survey_sf <- df %>%
  st_as_sf(coords = c("cluster_lon", "cluster_lat"), crs = 4326) %>%
  st_transform(3857)

track_sf <- df %>%
  st_as_sf(coords = c("longitude_storm", "latitude_storm"), crs = 4326) %>%
  st_transform(3857)

# calculate daily distance from each survey site to its paired storm track point
  # reflects distance to the strack point on each day of the storm
df <- df %>%
  mutate(daily_dist_to_track_m = as.numeric(
    st_distance(survey_sf, track_sf, by_element = TRUE)
  ))


# convert cluster and storm locations to sf objects for closest track detection
survey_pts <- df %>%
  distinct(site_id, sid, cluster_lon, cluster_lat, .keep_all = TRUE) %>%
  st_as_sf(coords = c("cluster_lon", "cluster_lat"), crs = 4326) %>%
  st_transform(3857)

storm_pts <- df %>%
  distinct(sid, date_storm, longitude_storm, latitude_storm, roci_m, .keep_all = TRUE) %>%
  st_as_sf(coords = c("longitude_storm", "latitude_storm"), crs = 4326) %>%
  st_transform(3857)


# find the nearest storm track point to each survey site (closest approach overall)
closest_track <- survey_pts %>%
  group_by(sid, site_id) %>%
  group_modify(~ {
    storm_subset <- storm_pts %>% filter(sid == .y$sid)
    .x %>%
      mutate(
        nearest_idx = st_nearest_feature(geometry, storm_subset),
        date_storm_nearest_site = storm_subset$date_storm[nearest_idx],
        min_dist_to_track_m = as.numeric(st_distance(geometry, storm_subset[nearest_idx, ], by_element = TRUE))
      )
  }) %>%
  ungroup()


# Fine link: join the date a storm was nearest a survey site (day zero) to the df and calculate the dates of surveys around day zero
df_flagged <- df %>%
  left_join(closest_track %>% select(site_id, sid, date_storm_nearest_site, min_dist_to_track_m),
            by = c("site_id", "sid"), relationship = "many-to-many") %>%
  mutate(
    fine_storm_phase = case_when(
      date_nurdle - date_storm_nearest_site < 0 ~ 'pre',
      date_nurdle - date_storm_nearest_site == 0 ~ 'during',
      date_nurdle - date_storm_nearest_site > 0 ~ 'post'))


# b) Max phase: Flag whether a survey occurred pre|during|post the full storm start and end (complete_storm_phase)
df_flagged <- df_flagged %>%
  mutate(
    complete_storm_phase = case_when(
      date_nurdle < full_storm_start ~ "pre",
      date_nurdle >= full_storm_start & date_nurdle <= full_storm_end ~ "during",
      date_nurdle > full_storm_end ~ "post"))
  


# 4) Flag post-storm surveys based on temporal range (immediate = 0-1d, acute = 2-4d, sub-acute = 5-10, extended = 11-28)
  # detect number of days fine and complete response are from storm bth pre- and post
df_flagged <- df_flagged %>%
  mutate(date_storm_nearest_site = as.Date(date_storm_nearest_site))

df_flagged <- df_flagged %>%
  mutate(
    fine_days = as.numeric(date_nurdle - date_storm_nearest_site),
    complete_days = case_when(
      date_nurdle < full_storm_start ~ as.numeric(date_nurdle - full_storm_start),
      date_nurdle > full_storm_end ~ as.numeric(date_nurdle - full_storm_end),
      TRUE ~ 0
    )
  )


# flag responses: not analyzing complete response, but retaining in case a reviewer asks
df_flagged <- df_flagged %>%
  mutate(
    fine_response = case_when( # case_when works by short-circuit evaluation so values are overwritten
      fine_days >= -4 & fine_days <= 4 ~ "acute",
      fine_days >= -10 & fine_days <= 10 ~ "subacute",
      fine_days >= -30 & fine_days <= 30 ~ "extended"),
    
    complete_response = case_when(
      complete_days >= - 4 & complete_days <= 4 ~ "acute",
      complete_days >= -10 & complete_days <= 10 ~ "subacute",
      complete_days >= -28 & complete_days <= 28 ~ "extended"))


# complete storm phase = 28d, so some sites have a fine_scale of > 28 days: drop fine_days > +/- 30days 
 # MAKE SURE THIS IS CLEAR IN METHODS bc link original on complete storm phase 3.23.26
df_flagged <- df_flagged %>%
  filter(fine_days < 31 & fine_days > -31)


# 5) Flag surveys as storm_influenced
df_flagged <- df_flagged %>%
  mutate(survey_type = "storm_influenced")


# 6) get eligible site_ids for control data post temporal filtering of storm rows
eligible_sites <- df_flagged %>%
  group_by(site_id, sid) %>%
  filter(all(c("pre", "post") %in% fine_storm_phase)) %>%
  ungroup() %>%
  distinct(site_id) %>%
  pull()
 

# 7) build control data from nurdle df - avoids data explosions from control surveys linking to storm rows
survey_dates_by_site <- df_flagged %>%
  distinct(site_id, date_nurdle)

control_df <- nurdle %>%
  filter(site_id %in% eligible_sites) %>%
  anti_join(survey_dates_by_site, by = c("site_id", "date_nurdle")) %>% # dates that are not in storm-flagged df
  mutate(survey_type = "control")


# bind - storm rows have storm cols, control rows have NAs for storm cols naturally
master_df <- bind_rows(df_flagged, control_df)



# Temporal Clean ----------------------------------------------------------

# 1) create a linked_id pairing site_id with linked storm
master_df <- master_df %>%
  group_by(site_id, sid) %>%
  mutate(linked_id = paste0(site_id, "_", name))


# 2) reorder cols
master_df <- master_df %>%
  relocate(c(linked_id, date_nurdle, name, date_storm, survey_type, fine_storm_phase, fine_response, complete_storm_phase, complete_response))


# 3) clean df col names
master_df <- clean_names(master_df)
#View(master_df)


# Spatial Link  -------------------------------------------

# 1) linking on ROCI, so must drop rows with no ROCI info 
master_df <- master_df %>%
  filter(survey_type == "control" | (is.finite(roci_m) & roci_m > 0)) # keep control data and storms w/ roci
#View(master_df)


# 2) remake survey_pts and storm_pts (sf objects) with the temporally linked data
survey_pts <- master_df %>%
  filter(survey_type == "storm_influenced") %>%  # controls have no linked_id
  distinct(linked_id, cluster_lon, cluster_lat) %>%
  st_as_sf(coords = c("cluster_lon", "cluster_lat"), crs = 4326) %>%
  st_transform(3857)


# detect distinct coordinates of each storm
storm_pts <- master_df %>%
  filter(survey_type == "storm_influenced") %>% # ignore control rows since it has no roci
  distinct(sid, name, date_storm, longitude_storm, latitude_storm, roci_m) %>%
  st_as_sf(coords = c("longitude_storm", "latitude_storm"), crs = 4326) %>%
  st_transform(3857)


# 3) using storm and cluster sf objects created above, buffer each storm pt by its own roci_m in all directions
storm_pts <- st_buffer(storm_pts, dist = storm_pts$roci_m)


# 4) spatially join site to storm-buffer (keeps only intersecting pairs)
  # first, identify sites that spatially intersect w/ storm roci
site_storm_pairs <- st_join(
  survey_pts %>% select(linked_id, geometry), # select only necessary cols
  storm_pts %>% select(sid, name, date_storm, geometry),
  join = st_intersects, # collect intersecting geometry
  left = FALSE) %>%
  st_drop_geometry() %>% # drop geometry
  filter(str_extract(linked_id, "[^_]+$") == name) %>%  # storm name in linked_id must match
  distinct(linked_id) # retain site_storm pairs that intersect (by roci)



# no linked id in control data, so need to extract site_id from linked_id in site_storm_pairs to detect same sites in control data
eligible_sites <- site_storm_pairs %>%
  mutate(site_id = str_remove(linked_id, "_[^_]+$")) %>%  # [^_]+$ means 1 or more char that isn't _, and $ anchors it to the end, removing only the storm name
  distinct(site_id) %>%
  pull()


# 5) for sites that are ever spatiotemporally linked to at least 1 storm, retain the surveys that failed spatial link of other storms as controld data for that site.
spatial_additions <- master_df %>%
  filter(survey_type == "storm-influenced",
         !linked_id %in% site_storm_pairs$linked_id, # failed spatial filter
         site_id %in% eligible_sites  # but site still eligible b
  ) %>%
  mutate(
    survey_type = "control",
    linked_id = NA_character_,
    fine_storm_phase = "control",
    fine_response = "control",
    complete_storm_phase = "control",
    complete_response = "control"
    
  )

# lastly, keep all rows that match site-storm pairs (even control data)
spatial <- master_df %>%
  filter(linked_id %in% site_storm_pairs$linked_id |
      (survey_type == "control" & site_id %in% eligible_sites))


# 6) Now, remove site_storm links that don't have observations both pre and post the storm (fine scale)
eligible_linked_ids <- spatial %>%
  filter(survey_type == 'storm_influenced') %>%
  group_by(linked_id) %>%
  filter(all(c("pre", "post") %in% fine_storm_phase)) %>%
  ungroup %>%
  distinct(linked_id) %>%
  pull()

# get site id from these linked_ids 
eligible_sites <- spatial %>%
  filter(linked_id %in% eligible_linked_ids) %>%
  distinct(site_id) %>%
  pull()

# now use this to retain control data only for sites that met the linked_id criteria
spatial <- spatial %>%
  filter(linked_id %in% eligible_linked_ids | 
           (survey_type == "control" & site_id %in% eligible_sites))



# Get site name ----------------------------------------------------------
cluster_names <- spatial %>%
  group_by(site_id) %>%
  summarise(
    # most common location name within cluster as the primary label, since clusters can contain several specific points inside (lose fine spatial resolution when clustering points)
    primary_location = names(sort(table(location[location != "" & !is.na(location)]), decreasing = TRUE))[1], # most common locaiton in cluster
    n_locations = n_distinct(location),
    state = first(state),
    country = first(country),
    cluster_lon = first(cluster_lon), # retain these cols for duplicate locations
    cluster_lat = first(cluster_lat), 
    .groups = "drop"
  ) %>%
  mutate(
    site_name = case_when(
      n_locations > 1 & state != "" ~ paste0(primary_location, " area, ", state), # generalize location (US)
      n_locations > 1 & state == "" ~ paste0(primary_location, " area, ", country), # generalize location (Mexico)
      state != "" ~ paste0(primary_location, ", ", state), # if only 1 specific location in cluster (US)
      TRUE ~ paste0(primary_location, ", ", country) # if only 1 specific location in cluster (Mexico)
    )
  )

# find cardinal direction for duplicate locations
cluster_names <- cluster_names %>%
  group_by(site_name) %>%
  mutate(
    n_duplicates = n(),
    centroid_lat = mean(cluster_lat),
    centroid_lon = mean(cluster_lon),
    cardinal = case_when(
      n_duplicates == 1 ~ "",
      cluster_lat > centroid_lat & cluster_lon > centroid_lon ~ "NE",
      cluster_lat > centroid_lat & cluster_lon < centroid_lon ~ "NW",
      cluster_lat < centroid_lat & cluster_lon > centroid_lon ~ "SE",
      cluster_lat < centroid_lat & cluster_lon < centroid_lon ~ "SW",
      cluster_lat > centroid_lat ~ "N",
      cluster_lat < centroid_lat ~ "S",
      cluster_lon > centroid_lon ~ "E",
      cluster_lon < centroid_lon ~ "W",
      TRUE  ~ ""
    ),
    site_name = if_else(cardinal != "", paste0(site_name, " ", cardinal), site_name)
  ) %>%
  ungroup() %>%
  select(site_id, site_name)


# join site names to df
spatial <- spatial %>%
  left_join(cluster_names %>% select (site_id, site_name), by = "site_id")



# Reduce Resolution -------------------------------------------------------
  # Collapse multiple storm observations per day (~6) into single observation per phase/day, this matches survey resolution (1/day)
  # Takes MAX values for storm characteristics (greatest effect that day)
  # Use this for storm characteristic analyses

spatial <- spatial %>%
  group_by(site_id, date_nurdle, sid, date_storm) %>%
  summarise(
    across(where(is.numeric), ~ max(., na.rm = FALSE)), # for numeric cols, take max value
    across(!where(is.numeric), ~ first(.)),             # for non-numeric cols, take first
    .groups = "drop"
  )


# convert storm_phase and storm_responses for control data from NA to "control"
spatial <- spatial %>%
  mutate(
  fine_storm_phase = coalesce(fine_storm_phase, "control"),
  fine_response = coalesce(fine_response, "control"),
  complete_storm_phase = coalesce(fine_storm_phase, "control"),
  complete_response = coalesce(fine_storm_phase, "control")
)

# save master df
write.csv(spatial, file = here('output','master_df.csv'), row.names = FALSE) # reworked this code to retain control data on 3.7.26 following OSM

