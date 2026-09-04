# SCRIPT 7: Storm Characteristics & MP Response Correlations
# Purpose: Characterize storms, compute storm geometry relative to sites, and test whether storm metrics predict MP redistribution
  # Load Libraries
  # Read Data
  # Config
  # Section 1: Storm Duration
  # Section 2: Storm Metrics per Site-Storm
  # Section 3: Storm Intensity Analysis
  # Section 4: Storm Geometry
  # Section 5: Storm Quadrant Analysis
  # Section 6: MP response by Storm Metric (lme and plots)
  # Section 7: Correlations
    # A) Overall pre-post delta (site-storm level)
    # B) By fine_response window (acute, subacute, extended)
  # Section 8: Predictor Models
    # A) Magnitude (GMER)
    # B) Direction (CLMM)

# Load Libraries ----------------------------------------------------------
pacman::p_load('dplyr', 'tidyr','ggplot2', 'cplm', 'here', 'ordinal',  'jtools')
source(here("functions.R"))


# Read Data ---------------------------------------------------------------
master_df <- read.csv(here("output", "master_df.csv")) %>%
  mutate(date_nurdle = as.Date(date_nurdle),
         date_storm = as.Date(date_storm))

# MP change
mp_change <- read.csv(here('output','summary_tables', 'MP_delta_bySite.csv')) %>%
  mutate(site_id = as.integer(as.character(site_id)))

# response window MP abundance
response <- read.csv(here("output", "summary_tables", "site_temp_window_abund.csv"))

# per-storm paired test results -- produced by Script 6 (long format: one row per storm x metric x test)
per_storm_results <- read.csv(here("output", "storm_level", "stats", "pre_post", strand,
                                   paste0("per_storm_paired_stats_", strand, ".csv")))


# Config ------------------------------------------------------------------
strand <- "total"    # options: old, new, total
stat <- "median" # options: max, mean, median, sum -> median correct basedon skew of data
alpha = 0.05

strand_cols <- list(
  old = list(raw = "weighted_old_strandline", mean_resp = "mean_abundance_old", median_resp = "median_abundance_old",
               se_resp = "se_old", q1_resp = "q1_abundance_old", q3_resp = "q3_abundance_old"),
  new = list(raw = "weighted_new_strandline", mean_resp = "mean_abundance_new", median_resp = "median_abundance_new",
               se_resp = "se_new", q1_resp = "q1_abundance_new", q3_resp = "q3_abundance_new"),
  total = list(raw = "weighted_total_amount", mean_resp = "mean_abundance_total", median_resp = "median_abundance_total",
               se_resp = "se_total", q1_resp = "q1_abundance_total", q3_resp = "q3_abundance_total")
)

col <- strand_cols[[strand]]
resp_col <- paste0(stat, "_resp")   # resolves to "mean_resp" or "median_resp" -- this is what actually makes resp_delta stat-aware

pcs_delta_col   <- paste0("pcs_change_", stat, "_", strand)
change_type_col <- paste0("change_type_", stat, "_", strand)


# standardize column aliases for plotting
mp_change[["pcs_change"]]  <- mp_change[[pcs_delta_col]]
mp_change[["change_type"]] <- mp_change[[change_type_col]]


# pivot response data for plotting
resp_delta <- response %>%
  filter(fine_storm_phase %in% c("pre", "post")) %>%
  select(site_id, name, fine_storm_phase, fine_response,
         central_abund = all_of(col[[resp_col]]), # median
         q1_abund = all_of(col$q1_resp), #IQR
         q3_abund = all_of(col$q3_resp)) %>%
  pivot_wider(names_from = fine_storm_phase,
              values_from = c(central_abund, q1_abund, q3_abund))

# detect change type
resp_delta <- resp_delta %>%
  mutate(
    pcs_change = central_abund_post - central_abund_pre,
    change_type = case_when(
      is.na(pcs_change) ~ NA_character_,
      pcs_change > 0 ~ "accumulation",
      pcs_change < 0 ~ "depletion",
      TRUE ~ "no_change"
    ),
    fine_response = factor(fine_response, levels = c("acute", "subacute", "extended"))
  )


# defining usa_sshs integer storm levels with IBTrACS labels
sshs_labels <- c(
  "-5" = "NC",   # not classified
  "-4" = "PT",   # post-tropical
  "-3" = "EX",   # extratropical
  "-2" = "SD",   # subtropical depression
  "-1" = "TD",   # tropical depression
  "0"  = "TS",   # tropical storm
  "1"  = "Cat1",
  "2"  = "Cat2",
  "3"  = "Cat3",
  "4"  = "Cat4",
  "5"  = "Cat5"
)

# Section 1) Storm Duration ----------------------------------------------

storm_duration <- master_df %>%
  filter(survey_type == "storm_influenced") %>%
  distinct(name, full_storm_start, full_storm_end) %>%
  mutate(full_storm_start = as.Date(full_storm_start),
         full_storm_end = as.Date(full_storm_end)) %>%
  group_by(name) %>%
  mutate(duration_days = if_else(full_storm_end == full_storm_start, 1,
                                  as.numeric(full_storm_end - full_storm_start))) %>%
  arrange(desc(duration_days))

write.csv(storm_duration, here("output", "storm_level", "stats", "storm_metrics", "storm_durations.csv"), row.names = FALSE)

# Section 2) Storm Metrics per Site-Storm --------------------------

# storm metrics on day that storm was closest to a site (day 0, i.e. closest approach)
metrics_closest <- master_df %>%
  filter(survey_type == "storm_influenced") %>%
  distinct(linked_id, date_storm_nearest_site, .keep_all = TRUE) %>%
  group_by(site_id, name) %>%
  summarise(
    usa_wind = first(usa_wind),
    rmw_m = first(rmw_m),
    roci_m = first(roci_m),
    max_r_nw_m = suppressWarnings(max(c(r34_nw_m, r50_nw_m, r64_nw_m), na.rm = TRUE)), # suppress NA warnings
    max_r_ne_m = suppressWarnings(max(c(r34_ne_m, r50_ne_m, r64_ne_m), na.rm = TRUE)),
    max_r_sw_m = suppressWarnings(max(c(r34_sw_m, r50_sw_m, r64_sw_m), na.rm = TRUE)),
    max_r_se_m = suppressWarnings(max(c(r34_se_m, r50_se_m, r64_se_m), na.rm = TRUE)),
    storm_speed = first(storm_speed), # makes assumption clear, if grouping assumption breaks anywhere, first() fails instead of silent arbitrary row picks
    usa_pres = first(usa_pres),
    usa_status = first(usa_status),
    usa_sshs = first(usa_sshs),
    min_dist_to_track_m = first(min_dist_to_track_m),
    storm_dir = first(storm_dir),
    storm_speed = first(storm_speed),
    dist2land = first(dist2land),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ifelse(is.infinite(.), NA_real_, .))) # replace -Inf values with true NA


write.csv(metrics_closest, here("output", "site_level", "stats", "storm_metrics", "storm_metrics_closest_day_to_site.csv"), row.names = FALSE)



# GoM-level storm summary
overall_storm_summaries <- master_df %>%
  filter(survey_type == "storm_influenced") %>%
  group_by(name) %>%
  summarise(
    n_survey_sites = n_distinct(site_id),
    max_wind_kts = max(usa_wind, na.rm = TRUE),
    mean_wind_kts = mean(usa_wind, na.rm = TRUE),
    median_wind_kts = median(usa_wind, na.rm = TRUE),
    min_usa_pres = min(usa_pres, na.rm = TRUE),
    mean_usa_pres = mean(usa_pres, na.rm = TRUE),
    median_usa_pres = median(usa_pres, na.rm = TRUE),
    max_sshs = max(usa_sshs, na.rm = TRUE),
    max_speed_kts = max(storm_speed, na.rm = TRUE),
    mean_speed_kts = mean(storm_speed, na.rm = TRUE),
    median_speed_kts = median(storm_speed, na.rm = TRUE),
    primary_usa_sshs = names(which.max(table(usa_sshs))),
    primary_usa_status = names(which.max(table(usa_status))),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ifelse(is.infinite(.), NA_real_, .))) %>%
  left_join(storm_duration %>% 
              distinct(name, .keep_all = TRUE) %>%
              select(name, duration_days, full_storm_start, full_storm_end), by = "name") %>%
  arrange(desc(n_survey_sites))


write.csv(overall_storm_summaries, here("output", "GoM_level", "stats", "storm_metrics", "GoM_storm_metrics.csv"), row.names = FALSE)



# Section 3) Storm Intensity Analysis -------------------------------------
  # Question: How does storm category impact MP? 
  # looking at peak, primary, and closest approach intensities

per_storm_stat_results <- per_storm_results %>%
  filter(metric == stat, test == "paired_wilcox_ttest")

# ---- A) Peak intensity (storm-level, max category anywhere on track) ----
intensity_results_peak <- per_storm_stat_results %>%
  left_join(
    master_df %>%
      filter(survey_type == "storm_influenced") %>%
      group_by(name) %>%
      slice_max(usa_sshs, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(name, max_category = usa_sshs, peak_usa_status = usa_status),
    by = "name"
  ) %>%
  mutate(
    intensity_group = case_when(
      peak_usa_status == "TS" ~ "tropical_storm",
      peak_usa_status == "TD" ~ "tropical_depression",
      peak_usa_status %in% c("DB", "LO") ~ "disturbance_low",
      peak_usa_status %in% c("EX", "SD", "PT", "NC") ~ "subtropical_or_lower",
      max_category == 1 ~ "cat1",
      max_category == 2 ~ "cat2",
      max_category == 3 ~ "cat3",
      max_category == 4 ~ "cat4"
    )
  ) %>%
  filter(!is.na(intensity_group), !is.na(median_delta))

intensity_kruskal_peak <- kruskal.test(median_delta ~ intensity_group, data = intensity_results_peak)

intensity_summary_peak <- intensity_results_peak %>%
  group_by(intensity_group) %>%
  summarise(
    n_storms = n(),
    median_delta = median(median_delta, na.rm = TRUE),
    n_accumulation = sum(n_accumulation, na.rm = TRUE),
    n_depletion = sum(n_depletion, na.rm = TRUE),
    n_no_change = sum(n_no_change, na.rm = TRUE),
    pct_accumulation = n_accumulation / (n_accumulation + n_depletion + n_no_change) * 100,
    pct_depletion = n_depletion / (n_accumulation + n_depletion + n_no_change) * 100,
    pct_no_change = n_no_change / (n_accumulation + n_depletion + n_no_change) * 100,
    .groups = "drop"
  )

write.csv(intensity_summary_peak, here("output", "storm_level", "stats", "storm_metrics", paste0("GoM_intensity_summary_peak_", strand, ".csv")), row.names = FALSE)
write.csv(data.frame(strand = strand, test = "Kruskal-Wallis", statistic = intensity_kruskal_peak$statistic,
                     p_value = intensity_kruskal_peak$p.value, sig = intensity_kruskal_peak$p.value < alpha),
          here("output", "storm_level", "stats", "storm_metrics", paste0(strand, "_storm_intensity_test_peak.csv")), row.names = FALSE)


# ---- B) Primary/dominant intensity (storm-level, mode across track) ----
intensity_results_primary <- per_storm_stat_results %>%
  left_join(overall_storm_summaries %>% select(name, primary_usa_status, primary_usa_sshs), by = "name") %>%
  mutate(
    primary_usa_sshs_num = as.numeric(primary_usa_sshs),
    intensity_group = case_when(
      primary_usa_status == "TS" ~ "tropical_storm",
      primary_usa_status == "TD" ~ "tropical_depression",
      primary_usa_status %in% c("DB", "LO") ~ "disturbance_low",
      primary_usa_status %in% c("EX", "SD", "PT", "NC") ~ "subtropical_or_lower",
      primary_usa_sshs_num == 1 ~ "cat1",
      primary_usa_sshs_num == 2 ~ "cat2",
      primary_usa_sshs_num == 3 ~ "cat3",
      primary_usa_sshs_num == 4 ~ "cat4"
    )
  ) %>%
  filter(!is.na(intensity_group), !is.na(median_delta))

intensity_kruskal_primary <- kruskal.test(median_delta ~ intensity_group, data = intensity_results_primary)

intensity_summary_primary <- intensity_results_primary %>%
  group_by(intensity_group) %>%
  summarise(
    n_storms = n(),
    median_delta = median(median_delta, na.rm = TRUE),
    n_accumulation = sum(n_accumulation, na.rm = TRUE),
    n_depletion = sum(n_depletion, na.rm = TRUE),
    n_no_change = sum(n_no_change, na.rm = TRUE),
    pct_accumulation = n_accumulation / (n_accumulation + n_depletion + n_no_change) * 100,
    pct_depletion = n_depletion / (n_accumulation + n_depletion + n_no_change) * 100,
    pct_no_change = n_no_change / (n_accumulation + n_depletion + n_no_change) * 100,
    .groups = "drop"
  )

write.csv(intensity_summary_primary, here("output", "storm_level", "stats", "storm_metrics", paste0("GoM_intensity_summary_primary_", strand, ".csv")), row.names = FALSE)
write.csv(data.frame(strand = strand, test = "Kruskal-Wallis", statistic = intensity_kruskal_primary$statistic,
                     p_value = intensity_kruskal_primary$p.value, sig = intensity_kruskal_primary$p.value < alpha),
          here("output", "storm_level", "stats", "storm_metrics", paste0(strand, "_storm_intensity_test_primary.csv")), row.names = FALSE)


# ---- C) Closest-approach (site-storm interaction level) ------
intensity_results_site <- mp_change %>%
  mutate(site_id = as.integer(as.character(site_id))) %>%
  select(site_id, name, pcs_change = all_of(pcs_delta_col), change_type = all_of(change_type_col)) %>%
  left_join(metrics_closest %>% select(site_id, name, usa_status, usa_sshs), by = c("site_id", "name")) %>%
  mutate(
    intensity_group = case_when(
      usa_status == "TS" ~ "tropical_storm",
      usa_status == "TD" ~ "tropical_depression",
      usa_status %in% c("DB", "LO") ~ "disturbance_low",
      usa_status %in% c("EX", "SD", "PT", "NC") ~ "subtropical_or_lower",
      usa_sshs == 1 ~ "cat1",
      usa_sshs == 2 ~ "cat2",
      usa_sshs == 3 ~ "cat3",
      usa_sshs == 4 ~ "cat4"
    )
  ) %>%
  filter(!is.na(intensity_group), !is.na(pcs_change))

intensity_kruskal_site <- kruskal.test(pcs_change ~ intensity_group, data = intensity_results_site)

intensity_summary_site <- intensity_results_site %>%
  group_by(intensity_group) %>%
  summarise(
    n_interactions = n(),
    median_delta = median(pcs_change, na.rm = TRUE),
    n_accumulation = sum(change_type == "accumulation", na.rm = TRUE),
    n_depletion = sum(change_type == "depletion", na.rm = TRUE),
    n_no_change = sum(change_type == "no_change", na.rm = TRUE),
    pct_accumulation = n_accumulation / n() * 100,
    pct_depletion = n_depletion / n() * 100,
    pct_no_change = n_no_change / n() * 100,
    .groups = "drop"
  )

write.csv(intensity_summary_site, here("output", "storm_level", "stats", "storm_metrics", paste0("GoM_intensity_summary_", strand, "_.csv")), row.names = FALSE)
write.csv(data.frame(strand = strand, test = "Kruskal-Wallis", statistic = intensity_kruskal_site$statistic,
                     p_value = intensity_kruskal_site$p.value, sig = intensity_kruskal_site$p.value < alpha),
          here("output", "storm_level", "stats", "storm_metrics", paste0(strand, "_storm_intensity_test.csv")), row.names = FALSE)



# Section 4) Storm Geometry  --------------------------------

# at the closest approach to a site get: bearing from storm center to site, storm heading, relative to site bearing, and quadrant classification
  # quadrant defined relative to storm direction: 
    # right-front: 0-90 deg -> ahead (0) & to the right (most impactful storm quadrant in NHemi)
    # right-rear: 90-180 -> behind (180) & to the right
    # left-rear: 180-270 -> behind (180) & to the left
    # left-front: 270-360 -> ahead(360) & to left
    # Note: 0 deg = site is due North of storm center

closest <- master_df %>%
  filter(survey_type == "storm_influenced") %>%
  distinct(linked_id, date_storm, .keep_all = TRUE) %>%
  filter(date_storm == date_storm_nearest_site)
  

# convert site & storm coords to radians
storm_geometry <- closest %>%
  mutate(
    # convert to radians (* pi/180) for trig functions
    lon_storm_rad = longitude_storm * pi / 180,
    lat_storm_rad = latitude_storm * pi / 180,
    lon_site_rad = cluster_lon * pi / 180,
    lat_site_rad = cluster_lat * pi / 180,
    dlon = lon_site_rad - lon_storm_rad, # difference in longitude b/w storm center and site
    
    # use great-circle bearing formula to cal. bearing from storm center to site (0 = N, 90 = E, clockwise)
    bearing_storm_to_site = (atan2(
      sin(dlon) * cos(lat_site_rad),
      cos(lat_storm_rad) * sin(lat_site_rad) - sin(lat_storm_rad) * cos(lat_site_rad) * cos(dlon)
    ) * 180 / pi + 360) %% 360, # * 180 / pi turns result back into degrees and + 360 %% 360 ensures the result is always b/w 0-360
    
    # get relative bearing (site position relative to storm heading)
    relative_bearing = (bearing_storm_to_site - storm_dir + 360) %% 360,
    
    # quadrant classification relative to storm motion
    storm_quadrant = case_when(
      relative_bearing <  90 ~ "right_front",
      relative_bearing < 180 ~ "right_rear",
      relative_bearing < 270 ~ "left_rear",
      TRUE ~ "left_front"
    )
  ) %>%
  select(site_id, name, linked_id,
         longitude_storm, latitude_storm, cluster_lon, cluster_lat,
         storm_dir, bearing_storm_to_site, relative_bearing, storm_quadrant,
         daily_dist_to_track_m)



# Section 5) Storm Quadrant Analysis ------------------------------------------
  # does storm quadrant predict whether a site accumulates or depletes MP?

# join quadrant and MP data
quad_data <- storm_geometry %>%
  left_join(
    mp_change %>% 
      select(site_id, name, pcs_change, change_type),
            by = c("site_id", "name")) 


# build a contingency table: quadrant x change_type
quad_contingency <- table(quad_data$storm_quadrant, quad_data$change_type)


# run quad_test (functions.R) & collect results
quad_result <- data.frame(
  strand = strand,
  test = quad_test$test,
  p_value = quad_test$result$p.value,
  sig = quad_test$result$p.value < alpha
)
quad_result

write.csv(quad_result, here("output","storm_level", "stats", "storm_metrics", paste0(strand, "_storm_quadrant_test.csv")), row.names = FALSE)



# make an inferential summary table of accum|dep counts per quad
quadrant_summary <- quad_data %>%
  group_by(storm_quadrant, change_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = change_type, values_from = n, values_fill = 0) %>%
  mutate(
    strand = strand,
    total = accumulation + depletion + no_change,
    pct_accum = round(accumulation / total * 100, 1),
    pct_deplete = round(depletion / total * 100, 1),
    pct_no_change = round(no_change / total * 100, 1)
  )

write.csv(quadrant_summary, here("output", "storm_level", "stats", "storm_metrics", strand, "GoM_quadrant_summary.csv"), row.names = FALSE)


# site-level quadrant interactions
  # does a specific site consistently sit in the same quadrant relative to strom's based on geo position relative to typical Gulf storm tracks
  # or, does a site-specific quadrant exposure explain that sites accum/dep pattern

# create a descriptive summary to show which quad a site most commonly occupies,
  # cross reference w/ script 6 per-site results to look for alignments with dominant change type of a site (helps inform case study locations maybe 4.12.26)
site_quad_summary <- quad_data %>%
  group_by(site_id, storm_quadrant, change_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(strand = strand)
site_quad_summary 

write.csv(site_quad_summary, here("output", "site_level", "stats", "storm_metrics", strand, "site_quadrant_change.csv"), row.names = FALSE)


# Section 6) Plots: MP response by Storm Metric when closest to site --------------------------------

# post-storm MP change
plot_data <- metrics_closest %>%
  left_join(mp_change %>%
              mutate(site_id = as.integer(as.character(site_id))) %>%
              distinct(site_id, name, .keep_all = TRUE) %>%
              select(site_id, name, pcs_change, change_type), by = c("site_id", "name")) %>%
  filter(!is.na(usa_sshs), !is.na(change_type), !is.na(pcs_change), !is.na(usa_wind))
  

# response window MP change
plot_data2 <- metrics_closest %>%
  left_join(resp_delta %>%
              select(site_id, name, fine_response, pcs_change, change_type), by = c("site_id", "name")) %>%
filter(!is.na(usa_sshs), !is.na(change_type), !is.na(pcs_change), !is.na(usa_wind))



#### Plot 1) Directional change by Intensity at Closest Approach ####
p_intensity <- intensity_results_site %>%
  mutate(
    intensity_group = factor(intensity_group,
                             levels = c("disturbance_low", "tropical_depression", "tropical_storm", "cat1", "cat2", "cat3", "cat4"),
                             labels = c("< TD", "TD", "TS", "Cat1", "Cat2", "Cat3", "Cat4")
    )
  ) %>%
  ggplot(aes(x = intensity_group, fill = change_type)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = colors) +
  labs(title = "MP Change by Storm Intensity",
       x = "Storm Intensity Category (at closest approach)",
       y = "Number of Sites",
       fill = "Post-storm MP Change Type") +
  theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13, facet.title.size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13), plot.title = element_text(size = 16))
p_intensity

ggsave(here("output", "storm_level", "graphs", "storm_metrics", strand, paste0("closest_intensity_direction_", strand, ".pdf")),
       plot = p_intensity, width = 8, height = 6, device = cairo_pdf)


#### Plot 2) Directional change by Peak Intensity ####
p_intensity_peak <- intensity_summary_peak %>%
  select(intensity_group, accumulation = n_accumulation, depletion = n_depletion, no_change = n_no_change) %>%
  pivot_longer(cols = c(accumulation, depletion, no_change), names_to = "change_type", values_to = "n") %>%
  mutate(
    intensity_group = factor(intensity_group,
                             levels = c("disturbance_low", "subtropical_or_lower", "tropical_depression", "tropical_storm", "cat1", "cat2", "cat3", "cat4"),
                             labels = c("Disturbance", "< TD", "TD", "TS", "Cat1", "Cat2", "Cat3", "Cat4")
    )
  ) %>%
  ggplot(aes(x = intensity_group, y = n, fill = change_type)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colors) +
  labs(title = "MP Change by Peak Storm Intensity",
       x = "Storm Peak Intensity Category (max category during storm's full track)",
       y = "Number of Sites",
       fill = "Post-storm MP Change Type") +
  theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13, facet.title.size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13), plot.title = element_text(size = 16))
p_intensity_peak

ggsave(here("output", "storm_level", "graphs", "storm_metrics", strand, paste0("peak_intensity_direction_", strand, ".pdf")),
       plot = p_intensity_peak, width = 8, height = 6, device = cairo_pdf)



#### Plot 3) Direction change by Primary Intensity ####
p_intensity_primary <- intensity_summary_primary %>%
  select(intensity_group, accumulation = n_accumulation, depletion = n_depletion, no_change = n_no_change) %>%
  pivot_longer(cols = c(accumulation, depletion, no_change), names_to = "change_type", values_to = "n") %>%
  mutate(
    intensity_group = factor(intensity_group,
                             levels = c("disturbance_low", "subtropical_or_lower", "tropical_depression", "tropical_storm", "cat1", "cat2", "cat3", "cat4"),
                             labels = c("Disturbance", "< TD", "TD", "TS", "Cat1", "Cat2", "Cat3", "Cat4")
    )
  ) %>%
  ggplot(aes(x = intensity_group, y = n, fill = change_type)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colors) +
  labs(title = "MP Directional Change by Primary Storm Intensity",
       x = "Storm Primary Intensity Category (dominant status across track)",
       y = "Number of Sites",
       fill = "Post-storm MP Change Type") +
  theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13, facet.title.size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13), plot.title = element_text(size = 16))
p_intensity_primary

ggsave(here("output", "storm_level", "graphs", "storm_metrics", strand, paste0("primary_intensity_direction_", strand, ".pdf")),
       plot = p_intensity_primary, width = 8, height = 6, device = cairo_pdf)



#### Plot 4) MP change vs storm characteristic predictors scatterplots ####
predictors <- c("usa_wind", "usa_sshs", "min_dist_to_track_m",
                "roci_m", "rmw_m", "storm_speed", "usa_pres", "storm_dir",
                "bearing_storm_to_site", "relative_bearing", "dist2land")

plot_data <- plot_data %>%
  left_join(storm_geometry %>% select(site_id, name, bearing_storm_to_site, relative_bearing), by = c("site_id", "name"))

predictor_labels <- c(
  usa_wind = "Wind Speed (kt)",
  usa_sshs = "Saffir-Simpson Category",
  min_dist_to_track_m = "Min Distance to Track (m)",
  roci_m = "ROCI (m)",
  rmw_m = "Radius of Max Wind (m)",
  storm_speed = "Storm Speed (kt)",
  usa_pres = "Pressure (mb)",
  storm_dir = "Storm Direction (deg)",
  bearing_storm_to_site = "Bearing, Storm to Site (deg)",
  relative_bearing = "Relative Bearing (deg)",
  dist2land = "Distance to Land (m)"
)

for (p in predictors) {
  
  plot_df <- plot_data %>% filter(!is.na(.data[[p]]), !is.na(pcs_change))
  
  p_scatter <- ggplot(plot_df, aes(x = .data[[p]], y = pcs_change, color = change_type)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
    scale_color_manual(values = colors) +
    labs(title = paste0("MP Delta vs ", predictor_labels[[p]], ": ", strand, " strand"),
         x = predictor_labels[[p]], y = "Pcs Change (post-pre)", color = "Response") +
    theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13) +
    theme(plot.title = element_text(size = 16))
  
  ggsave(here("output", "storm_level", "graphs", "storm_metrics", strand, paste0(p, "_vs_MPdelta_", strand, ".pdf")),
         plot = p_scatter, width = 8, height = 6, device = cairo_pdf)
  
}



#### Plot 5) Directional change count by response window - general plot, not storm metric (added 4.25.26) ####

p_window_dir <- ggplot(plot_data2, aes(x = fine_response, fill = change_type)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = colors) +
  labs(title = paste0("Directional Change by Temporal Window: ", strand, " strand"),
       x = "Response Window", y = "Number of Sites", fill = "Response") +
  theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13) +
  theme(plot.title = element_text(size = 16))
p_window_dir

ggsave(here("output", "storm_level", "graphs", "storm_metrics", strand, paste0("temp_window_directional_change_", strand, ".pdf")),
       plot = p_window_dir, width = 8, height = 6, dpi = 300, device = cairo_pdf)



#### Plot 6) MP change distribution by response window - general plot, not storm metric (added 4.25.26) ####
p_window_change <- ggplot(plot_data2, aes(x = fine_response, y = pcs_change, fill = change_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = colors) +
  labs(title = paste0("MP Change Distribution by Response Window: ", strand, " strand"),
       x = "Response Window", y = "Pcs Change (post - pre)", fill = "Response") +
  theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13) +
  theme(plot.title = element_text(size = 16))
p_window_change


# add stats to plot
resp_window_stats <- read.csv(here("output", "GoM_level", "stats", "resp_window", strand,
                                   paste0("GoM_resp_window_stats_", strand, ".csv")))

# per-group SD and significance annotations for the response-window boxplot
window_annotations <- plot_data2 %>%
  group_by(fine_response) %>%
  summarise(
    sd_pcs = sd(pcs_change, na.rm = TRUE),
    y_pos = min(pcs_change, na.rm = TRUE) - 0.05 * diff(range(plot_data2$pcs_change, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(
    resp_window_stats %>% filter(metric == stat) %>%
      select(fine_response = response_level, wilcox_p),
    by = "fine_response"
  ) %>%
  mutate(
    sig_stars = case_when(
      wilcox_p < 0.001 ~ "***",
      wilcox_p < 0.01 ~ "**",
      wilcox_p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = paste0("SD = ", round(sd_pcs, 1), " ", sig_stars)
  )

p_window_change <- p_window_change +
  geom_text(data = window_annotations, aes(x = fine_response, y = y_pos, label = label),
            inherit.aes = FALSE, fontface = "italic", size = 4)
p_window_change

ggsave(here("output", "storm_level", "graphs", "storm_metrics", strand, paste0("temp_window_dir_change_", strand, ".pdf")),
       plot = p_window_change, width = 8, height = 6, dpi = 300, device = cairo_pdf)






# Section 7) Correlations ----------------------------------------
# Spearman rho -- rank-based, does not assume linear relationships (data are right-skewed w/ outliers)

# build analysis df from storm geometry & mp_change: storm metrics + geometry + MP change per site-storm
corr_base <- mp_change %>%
  select(site_id, name, pcs_change, change_type) %>%
  left_join(metrics_closest, by = c("site_id", "name")) %>%
  left_join(storm_geometry %>% select(site_id, name, bearing_storm_to_site, 
                                      relative_bearing, storm_quadrant),
            by = c("site_id", "name")) %>%
  mutate(
    abs_pcs_change = abs(pcs_change),
    change_type_ord = factor(change_type,
                             levels = c("depletion", "no_change", "accumulation"),
                             ordered = TRUE)
  )

# check for missing predictors
missing_predictors <- predictors[!predictors %in% names(corr_base)]
if (length(missing_predictors) > 0) {
  warning(paste("These predictors were not found in corr_base and will be silently excluded:",
                paste(missing_predictors, collapse = ", ")))
}

predictors <- predictors[predictors %in% names(corr_base)]


#  Multicollinearity check (Spearman) 
# build correlation matrix data
corr_matrix_data <- corr_base %>%
  select(pcs_change, all_of(predictors)) %>%
  drop_na()


# test with Spearman correlation, which is rank-based & does not assume linear relationships in the data (which we know is right skewed w/ outliers)
corr_matrix <- cor(corr_matrix_data, method = "spearman", use = "complete.obs")

# convert to long format
corr_long <- as.data.frame(corr_matrix) %>%
  tibble::rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "r")

# correlation plot aka multicolinearity check
p_corr <- ggplot(corr_long, aes(x = var1, y = var2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(r, 2)), size = 3) +
  scale_fill_gradient2(low = "#0099cc", mid = "white", high = "#ff9900",
                       midpoint = 0, limits = c(-1,1), name = "r") +
  labs(
    title = "Spearman Correlation Matrix",
    subtitle = "Post-storm MP delta vs storm characteristics at closest approach",
    x = NULL,
    y = NULL) +
  theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13), plot.title = element_text(size = 16))

p_corr


ggsave(here("output", "storm_level", "graphs", "storm_metrics", strand, paste0("GoM_corr_matrix_", strand, ".pdf")),
       plot = p_corr, width = 10, height = 8, device = cairo_pdf)





# Section 8) Predictor Models ---------------------------------------------


# A) Magnitude model: does each predictor (individually) predict |pcs_change|? 
# one predictor per model
# storm intensity metrics are highly intercorrelated, so combining them in a single model makes coefficients unstable/uninterpretable
magnitude_results <- data.frame()
for (p in predictors) {
  z_col <- paste0(p, "_z")
  corr_base[[z_col]] <- as.numeric(scale(corr_base[[p]]))
  
  f <- reformulate(c(z_col, "(1 | name)"), response = "abs_pcs_change")
  model <- tryCatch(cpglmm(f, data = corr_base), error = function(e) NULL)
  if (is.null(model)) next
  
  cf <- summary(model)@coefs
  if (!(z_col %in% rownames(cf))) next
  
  magnitude_results <- rbind(magnitude_results, data.frame(
    predictor = p,
    estimate_z = cf[z_col, "Estimate"],
    std_error_z = cf[z_col, "Std. Error"],
    statistic = cf[z_col, "t value"],
    p_value = 2 * pnorm(abs(cf[z_col, "t value"]), lower.tail = FALSE)  # normal approximation; cplm does not report df-based p-values
  ))
}
magnitude_results$p_adj <- p.adjust(magnitude_results$p_value, method = "BH")
magnitude_results$sig_adj <- magnitude_results$p_adj < alpha
write.csv(magnitude_results, here("output", "storm_level", "stats", "storm_metrics", strand,
                                  paste0("magnitude_predictor_models_", strand, ".csv")), row.names = FALSE)

# B) Direction model: does each predictor predict change_type (depletion < no_change < accumulation)? 
# ordinal mixed model (clmm) -- change_type kept as a real 3-level outcome (not collapsed to binary),
# since a meaningful cluster of "no_change" would itself be an important finding (storms not
# redistributing MP, or insufficient sampling), and the 3 levels are naturally ordered along pcs_change.
direction_results <- data.frame()

for (p in predictors) {
  z_col <- paste0(p, "_z")
  corr_base[[z_col]] <- as.numeric(scale(corr_base[[p]]))
  
  f <- reformulate(c(z_col, "(1 | site_id)", "(1 | name)"), response = "change_type_ord")
  model <- tryCatch(ordinal::clmm(f, data = corr_base), error = function(e) NULL)
  if (is.null(model)) next
  
  cf <- coef(summary(model))
  if (!(z_col %in% rownames(cf))) next
  
  direction_results <- rbind(direction_results, data.frame(
    predictor = p,
    estimate_z = cf[z_col, "Estimate"],
    std_error_z = cf[z_col, "Std. Error"],
    statistic = cf[z_col, "z value"],
    p_value = cf[z_col, "Pr(>|z|)"],
    singular = isTRUE(model$info$max.grad > 1e-2)  # clmm has no isSingular(); see note below
  ))
}

direction_results$p_adj <- p.adjust(direction_results$p_value, method = "BH")
direction_results$sig_adj <- direction_results$p_adj < alpha

write.csv(direction_results, here("output", "storm_level", "stats", "storm_metrics", strand,
                                  paste0("direction_predictor_models_", strand, ".csv")), row.names = FALSE)



# Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), here("session_info.txt"))
renv::snapshot()



