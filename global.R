# This script runs once when the app starts, available to ui.R and server.R

# Understanding
  # whatever doesn't depend on user input goes here: Data loading, color palettes, constants, helper functions that don't need reactivity
  # to run in console: shiny::runApp()

# to call specific line number for an error run: source("global.R")
# to parse syntax errors parse("global.R")

# to host App on free Shiny io and update
# rsconnect::deployApp(forceUpdate = TRUE)

library(shiny)
library(leaflet)
library(dplyr)
library(sf)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(scales)
library(lubridate)




# Read Data --------------------------------------------------------------------

# complete spatiotemporally linked dataset
  # fine_storm_phase: "control" = non-storm surveys; "pre"/"post = storm-influenced surveys
master_df <- read.csv("data/master_df.csv") %>%
  filter(fine_storm_phase != "during") %>%
  mutate(
    site_id = as.character(site_id),
    fine_storm_phase = factor(fine_storm_phase, levels = c("control", "pre", "post")),
    fine_response = factor(fine_response, levels = c("acute", "subacute", "extended")),
    survey_type = factor(survey_type, levels = c("control", "storm_influenced"))
  )

# A) Regional summary tables
regional_byStormPhase <- read.csv("data/regional_yearly_abund_byStormPhase.csv") %>%
  filter(fine_storm_phase %in% c("control", "pre", "post")) %>%
  mutate(fine_storm_phase = factor(fine_storm_phase, levels = c("control", "pre", "post")))

regional_byCondition <- read.csv("data/regional_abund_bySurveyCondtion.csv") %>%
  mutate(survey_type = factor(survey_type, levels = c("control", "storm_influenced")))

# B) Storm-level summary tables
storm_abund_byPhase <- read.csv("data/storm_abund_byStormPhase.csv") %>%
  filter(fine_storm_phase != "during") %>%
  mutate(fine_storm_phase = factor(fine_storm_phase, levels = c("pre", "post")),
         name = factor(name, levels = sort(unique(name))))

storm_temp_window <- read.csv("data/storm_temp_window_abund.csv") %>%
  filter(fine_storm_phase != "during") %>%
  mutate(fine_storm_phase = factor(fine_storm_phase, levels = c("pre", "post")),
         fine_response = factor(fine_response, levels = c("acute", "subacute", "extended")),
         name = factor(name, levels = sort(unique(name))))

# C) Site-level summary tables
abund_bySiteStorm <- read.csv("data/siteStorm_abund_byStormPhase.csv") %>%
  filter(fine_storm_phase != "during") %>%
  mutate(
    site_id = as.character(site_id),
    site_id = factor(site_id, levels = as.character(sort(as.numeric(unique(site_id))))),
    name = factor(name, levels = sort(unique(name))),
    fine_storm_phase = factor(fine_storm_phase, levels = c("control", "pre", "post")))


# D) Redistribution / change table — MP_delta_bySite
mp_change <- read.csv("data/MP_delta_bySite.csv") %>%
  mutate(
    site_id = as.character(site_id),
    site_id = factor(site_id, levels = as.character(sort(as.numeric(unique(site_id))))),
    cluster_lon = as.numeric(cluster_lon),
    cluster_lat = as.numeric(cluster_lat),
    name = factor(name, levels = sort(unique(name))))

# E) Survey pairs for rate/direction plots
survey_pairs_total <- read.csv("data/survey_pairs_total.csv")

# F) Storm metrics
GoM_level_storm_info <- read.csv("data/GoM-level_storm_metrics.csv")

closest_metrics <- read.csv("data/storm_metrics_closest_day_to_site.csv") %>%
  mutate(
    site_id = gsub("_cluster", "", as.character(site_id)),
    min_dist_to_track_m = round(min_dist_to_track_m, 0))

# G) Plastic manufacturers
manu_df <- read.csv("data/manufacturer_df.csv")




# Map / spatial objects ---------------------------------------------------

GoM_lat <- 27.0
GoM_lon <- -91.0
GoM_zoom <- 6


# site choices for dropdowns
site_choices <- master_df %>%
  distinct(site_id, site_name) %>%
  arrange(as.numeric(site_id)) %>%
  mutate(label = paste(site_id, "-", site_name))

site_choice_vec <- setNames(site_choices$site_id, site_choices$label)

# storm choices for dropdowns
storm_choices <- c("All Paired Storms",
                   master_df %>%
                     filter(!is.na(name)) %>%
                     distinct(name, storm_year) %>%
                     mutate(label = paste(name, storm_year)) %>%
                              arrange(label) %>%
                              pull(label))

# chronological storm order for clean plotting
storm_date_order <- master_df %>%
  filter(!is.na(name)) %>%
  distinct(name, full_storm_start) %>%
  arrange(full_storm_start)



# storm track spatial objects, remove NA (control) rows
sf_df <- master_df %>%
  filter(!is.na(longitude_storm), !is.na(latitude_storm))


# create plain df of 1 storm point / day (works better than sf objects for leaflet)
storm_points <- sf_df %>%
  group_by(name, date_storm) %>%
  summarize(
    lat = mean(latitude_storm),
    lon = mean(longitude_storm), 
    roci = max(usa_roci, na.rm = TRUE),
    usa_wind = max(usa_wind, na.rm = TRUE),
    usa_pres = min(usa_pres, na.rm = TRUE),
    usa_sshs = max(usa_sshs, na.rm = TRUE),
    storm_speed = round(mean(storm_speed, na.rm = TRUE), 1),
    usa_status = max(usa_status),
    .groups = "drop")




# create one sf LINESTRING per storm for addPolylines (from 3_EDA.R)
storm_tracks <- sf_df %>%
  arrange(name, date_storm) %>%
  group_by(name) %>%
  summarize(n_obs = n(), .groups = "drop") %>%
  filter(n_obs > 1) %>%  # only storms with multiple observations can have tracks, note Imelda is lost
  left_join(
    sf_df %>% 
      arrange(name, date_storm) %>%
      select(name, longitude_storm, latitude_storm),
    by = "name"
  ) %>%
  st_as_sf(coords = c("longitude_storm", "latitude_storm"), crs = 4326) %>%
  group_by(name) %>%
  summarize(do_union = FALSE) %>%
  st_cast("LINESTRING")


# weight tracks by storm category for plotting
storm_tracks_cat <- storm_tracks %>%
  left_join(
    sf_df %>%
      group_by(name) %>%
      summarize(max_category = max(usa_sshs, na.rm = TRUE), .groups = "drop"),
    by = "name") %>%
  mutate(track_weight = pmax(max_category + 1, 1))


# Site popup lookup table -------------------------------------------------
  # pre-compute per-site control and storm means for map popups - one join here in global.R so poup renders instantly on click
  # this is a map view pop-up only, not analytical yet, so mean is okay I think, can make final choice before hosting (8.18.26)

site_ctrl_summary <-  abund_bySiteStorm %>%
  filter(fine_storm_phase == "control") %>%
  select(site_id,
         ctrl_mean_total = mean_abundance_total, 
         ctrl_mean_new = mean_abundance_new, 
         ctrl_mean_old = mean_abundance_old)


site_storm_summary <- abund_bySiteStorm %>%
  filter(fine_storm_phase %in% c("pre", "post")) %>%
  mutate(fine_storm_phase = as.character(fine_storm_phase)) %>%
  pivot_wider(
    names_from = fine_storm_phase,
    values_from = mean_abundance_total,
    names_prefix = "storm_mean_"
  )


# join strand summaries into site_popup_lookup
site_popup_lookup <- master_df %>%
  distinct(site_id, site_name, cluster_lon, cluster_lat) %>%
  left_join(site_ctrl_summary, by = "site_id") %>%
  left_join(site_storm_summary, by = "site_id")


# Pre-computed plotting objects -----------------------------------------------
  # combined yearly df for regional trend plot
  # mirrors EDA build_trend_plot(yearly_abund_byCondition, group_var = "survey_type")
combined_yearly <- regional_byCondition %>%
  mutate(condition = as.character(survey_type),
         condition = factor(condition, levels = c("control", "storm_influenced")))



# plot labeller -----------------------------------------------------------

nurdle_labeller <- labeller(
  fine_storm_phase = c("control" = "Control", "pre" = "Pre-storm", "post" = "Post-storm"),
  survey_type = c("control" = "Control", "storm_influenced" = "Storm-influenced"),
  fine_response = c("acute" = "Acute", "subacute" = "Subacute", "extended" = "Extended"),
  strandline_type = c("new" = "New Strandline", "old" = "Old Strandline"),
  strand = c("new" = "New Strandline", "old" = "Old Strandline")
)




# Aesthetic Details -------------------------------------------------------

colors <- c(
  "accumulation"     = "#D97B7B",
  "depletion"        = "#5B7C99",
  "no_change"        = "#A8A8A8",
  "control"          = "black",
  "pre"              = "#ffcccc",
  "storm_influenced" = "#ff0000",
  "storm"            = "#ff0000",
  "post"             = "#996666",
  "acute"            = "#FFCC00",
  "subacute"         = "#FF9900",
  "extended"         = "#993300"
)
