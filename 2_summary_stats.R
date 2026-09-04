# SCRIPT 2: SUMMARY STATISTICS
  # Purpose: Calculate summary statistics for clusters, storms, and MP abundances at different post-storm responses

# Load Libraries
# Read Data
# Section 1: Summary Counts
  # unique sites, unique storms, unique site-storm interactions, uniqie sites per storm
  # Number surveys per storm, per storm phase by site-storm interaction, storm phase by storm name,  per year, per year per storm condition/phase, per site per year, per month per site per year
# Section 2: Regional Abundance Summaries
  # per survey type (control, storm-influenced), per year per survey type, per site, per storm
# Section 3: Storm Phase (pre|post) Abundance Summaries
  # storm-level, site level, site-storm level
# Section 4: Temporal Window Response Abundance Summaries
  # yearly, site level, storm level
# Section 5: MP Change & Daily Abundance Deltas

# Load Libraries ----------------------------------------------------------
pacman::p_load('dplyr','tidyr','here','lubridate')
source(here("0_make_dirs.R"))

# Load functions
source(here("functions.R"))  


# Read Data ---------------------------------------------------------------
# Load master data from 0_link_storm_and_nurdle_df.R
master_df <- read.csv(here('output', 'master_df.csv')) %>%
  mutate(state = if_else(state == "", "Mexico", state)) # convert empty string to Mexico for survey site summaries on location



# Section 1: Summary survey Counts -----------------------------------------------

# 1) Unique survey sites
unique_sites <- master_df %>%
  group_by(site_id, site_name) %>%
  summarize(n_surveys = n_distinct(id), 
            latitude = first(cluster_lat),
            longitude = first(cluster_lon),
            .groups = "drop")

write.csv(unique_sites, here('output', 'summary_tables', 'distinct_sites.csv'), row.names = FALSE)


# 2) Unique storms
unique_storms <- master_df %>%
  drop_na(name) %>% # NA = control surveys
  group_by(sid, name, storm_year) %>%
  summarize(count = n(), .groups = "drop")

write.csv(unique_storms, here('output', 'summary_tables', 'distinct_storms.csv'), row.names = FALSE)


# 3) Unique site-storm interactions
unique_combos <- master_df %>%
  drop_na(name) %>%
  group_by(linked_id) %>%
  summarize(count = n(), .groups = "drop")

write.csv(unique_combos, here('output', 'summary_tables', 'distinct_combos.csv'), row.names = FALSE)


# 4) # survey sites per storm
sites_perStorm <- master_df %>%
  drop_na(name) %>%
  group_by(name) %>%
  summarise(n_clusters = n_distinct(site_id), .groups = "drop")

write.csv(sites_perStorm, here('output', 'summary_tables', 'sitesPerStorm.csv'), row.names = FALSE)


# 5) # surveys per storm phase by site-storm
surveys_perPhase_perSite <- master_df %>%
  group_by(linked_id, fine_storm_phase) %>%
  summarize(n_surveys = n_distinct(id),.groups = "drop")

write.csv(surveys_perPhase_perSite, here('output', 'summary_tables', 'surveysPerStormPhase_bySite.csv'), row.names = FALSE)


# 6) # surveys per storm phase by storm
surveys_perPhase <- master_df %>%
  group_by(name, fine_storm_phase) %>%
  summarize(n_surveys = n_distinct(id), .groups = "drop")

write.csv(surveys_perPhase, here('output','summary_tables', 'surveysPerStormPhase_byStorm.csv'), row.names = FALSE)


# 7) # surveys per year
surveys_per_year <- master_df %>%
  mutate(year = year(date_nurdle)) %>%
  group_by(year) %>%
  summarise(n_surveys = n_distinct(id), .groups = "drop") #id = unique surveys
surveys_per_year

write.csv(surveys_per_year, here('output','summary_tables', 'surveys_per_year.csv'), row.names = FALSE)


# 8) # surveys per year per storm condition/phase
surveys_per_condition_yr <- master_df %>%
  mutate(year = year(date_nurdle)) %>%
  group_by(year, fine_storm_phase) %>%
  summarise(n_surveys = n_distinct(id), .groups = "drop") #id = unique surveys
surveys_per_condition_yr

write.csv(surveys_per_condition_yr, here('output','summary_tables', 'surveys_per_condition_per_year.csv'), row.names = FALSE)


# 9) # surveys per site per year
surveys_per_site_yr <- master_df %>%
  mutate(year = year(date_nurdle)) %>%
  group_by(site_id, year) %>%
  summarise(n_surveys = n_distinct(id))

write.csv(surveys_per_site_yr, here('output','summary_tables', 'surveys_per_site_per_year.csv'), row.names = FALSE)


# 10) # surveys per month location per year (seasonality of surveying trends)
surveys_per_month <- master_df %>%
  mutate(month = month(date_nurdle)) %>%
  group_by(site_id, site_name, month) %>%
  summarise(n_surveys = n_distinct(id), .groups = "drop")

write.csv(surveys_per_month, here('output','summary_tables', 'surveys_per_month.csv'), row.names = FALSE)


# 11) Surveys per state
surveys_state <- master_df %>%
  group_by(state) %>%
  summarise(n_surveys = n_distinct(id), .groups = "drop")
surveys_state
write.csv(surveys_state, here('output','summary_tables', 'surveys_per_state.csv'), row.names = FALSE)


# 12 surveys per state per condition
surveys_state_cond <- master_df %>%
  group_by(state, fine_storm_phase) %>%
  summarise(n_surveys = n_distinct(id), .groups = "drop")
surveys_state_cond
write.csv(surveys_state_cond, here('output','summary_tables', 'surveys_per_state_condition.csv'), row.names = FALSE)


# 13 Surveys per survey type (per storm can double count surveys linked to >1 storm)
surveys_site_survey_type <- master_df %>%
  group_by(site_id, survey_type) %>%
  summarise(n_surveys = n_distinct(id))

write.csv(surveys_site_survey_type, here('output','summary_tables', 'surveys_site_survey_type.csv'), row.names = FALSE)



# Section 2) Regional Abundance Summaries ---------------------------------------------

# a) regional amount by survey condition
abund_byCondition <- master_df %>%
  distinct(id, .keep_all = TRUE) %>%
  group_by(survey_type) %>%
  summarise_abundance()

abund_byCondition
write.csv(abund_byCondition, here('output','summary_tables', 'regional_abund_bySurveyCondtion.csv'), row.names = FALSE)


# a2) regional amount by storm phase
abund_byStormPhase <- master_df %>%
  filter(fine_storm_phase != "during" & fine_storm_phase != "control") %>%
  distinct(id, .keep_all = TRUE) %>%
  group_by(fine_storm_phase) %>%
  summarise_abundance()

abund_byStormPhase
write.csv(abund_byStormPhase, here('output','summary_tables', 'regional_abund_byStormPhase.csv'), row.names = FALSE)


# b) regional amount by year - all conditions together
yearly_abund <- master_df %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(year = year(date_nurdle)) %>%
  group_by(year) %>%
  summarise_abundance()

yearly_abund
write.csv(yearly_abund, here('output','summary_tables', 'regional_yearly_abund.csv'), row.names = FALSE)


# c) regional yearly by storm phase
yearly_abund_phase <- master_df %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(year = year(date_nurdle)) %>%
  group_by(year, fine_storm_phase) %>%
  summarise_abundance()

yearly_abund_phase
write.csv(yearly_abund_phase, here('output','summary_tables', 'regional_yearly_abund_byStormPhase.csv'), row.names = FALSE)


# regional yearly by survey condition
yearly_abund_byCondition <- master_df %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(year = year(date_nurdle)) %>%
  group_by(survey_type, year) %>%
  summarise_abundance()

yearly_abund_byCondition
write.csv(yearly_abund_byCondition, here('output','summary_tables', 'regional_yearly_abund_bySurveyCondition.csv'), row.names = FALSE)


# regional by site
site_abund <- master_df %>%
  distinct(id, .keep_all = TRUE) %>%
  group_by(site_id) %>%
  summarise_abundance()

site_abund
write.csv(site_abund, here('output','summary_tables', 'regional_abund_bySite.csv'), row.names = FALSE)


# regional by storm
storm_abund <- master_df %>%
  drop_na(name) %>%
  distinct(id, .keep_all = TRUE) %>%
  group_by(name, storm_year, fine_storm_phase) %>%
  summarise_abundance()

storm_abund
write.csv(storm_abund, here('output','summary_tables', 'regional_abund_byStorm.csv'), row.names = FALSE)


# regional site per storm trends
siteStorm_abund <- master_df %>%
  drop_na(name) %>%
  distinct(id, .keep_all = TRUE) %>%
  group_by(site_id, name, storm_year) %>%
  summarise_abundance()

siteStorm_abund
write.csv(siteStorm_abund, here('output','summary_tables', 'regional_abund_bySiteStorm.csv'), row.names = FALSE)


# Section 3) Storm Phase (pre|post) Abundance Summaries ------------------------

# 1) storm-level abundance by storm phase 
storm_abund_phase <- master_df %>%
  drop_na(name) %>%  # control data has no storm name, i.e. "NA", drop these observations
  distinct(name, id, survey_type, fine_storm_phase, name, .keep_all = TRUE) %>% # handle repeated survey rows per day of storm obs.
  group_by(name, fine_storm_phase) %>% 
  summarise_abundance()

storm_abund_phase
write.csv(storm_abund_phase, here('output','summary_tables', 'storm_abund_byStormPhase.csv'), row.names = FALSE)


# 2) site-level abundance by storm phase (sites do contain control data, captured here)
  # Note: phase_sum will exceed site_abund at sites where surveys are assigned to multiple storm phases across different storm linkages (sites 13, 3, 58, 9).
  # This is expected since those surveys contribute to multiple phase categories and are not duplicates in the analytically sense.
site_abund_phase <- master_df %>%
  distinct(site_id, id, survey_type, fine_storm_phase, name, .keep_all = TRUE) %>%
  group_by(site_id, fine_storm_phase) %>% 
  summarise_abundance()

site_abund_phase
write.csv(site_abund_phase, here('output','summary_tables', 'site_abund_byStormPhase.csv'), row.names = FALSE)


# 3) site-storm-level abundance by storm phase 
siteStorm_abund_phase <- master_df %>%
  distinct(id, survey_type, fine_storm_phase, name, .keep_all = TRUE) %>% # handle repeated survey rows per day of storm obs.
  group_by(site_id, name, fine_storm_phase) %>% 
  summarise_abundance() %>%
  left_join(
    master_df %>% distinct(site_id, cluster_lat, cluster_lon), # retain site lon/lat
    by = "site_id"
  )

siteStorm_abund_phase
write.csv(siteStorm_abund_phase, here('output','summary_tables', 'siteStorm_abund_byStormPhase.csv'), row.names = FALSE)





# Section 4) Temporal Window Response Abundance Summaries  ------------------------

# 1) regional fine_response (acute|subacute|extended) by fine_storm_phase (control, pre, during, post)
yearly_temp_resp <- master_df %>%
  distinct(site_id, id, fine_storm_phase, fine_response, name, .keep_all = TRUE) %>% # handle repeated survey rows per day of storm obs.
  mutate(year = year(date_nurdle)) %>%
  group_by(year, fine_storm_phase, fine_response) %>%
  summarise_abundance()

yearly_temp_resp
write.csv(yearly_temp_resp, here('output','summary_tables', 'regional_yearly_temp_window_abund.csv'), row.names = FALSE)


# 2) Storm-level fine_response by storm phase
storm_temp_resp <- master_df %>%
  distinct(name, fine_storm_phase, fine_response, name, .keep_all = TRUE) %>% # handle repeated survey rows per day of storm obs.
  group_by(name, fine_storm_phase, fine_response) %>%
  summarise_abundance()

storm_temp_resp
write.csv(storm_temp_resp, here('output','summary_tables', 'storm_temp_window_abund.csv'), row.names = FALSE)


# 3) Site-level fine_response by storm phase (no control or during data) 
site_temp_resp <- master_df %>%
  distinct(site_id, id, fine_storm_phase, fine_response, name, .keep_all = TRUE) %>% # handle repeated survey rows per day of storm obs.
  group_by(site_id, name, fine_storm_phase, fine_response) %>%
  summarise_abundance()

site_temp_resp
write.csv(site_temp_resp, here('output','summary_tables', 'site_temp_window_abund.csv'), row.names = FALSE)





# Section 5) Daily Abundance Deltas -------------------------------------------


# 1) Calculate pre/post changes (pieces and percent) per site
  
mp_change <- siteStorm_abund_phase %>%
  filter(fine_storm_phase!= 'during') %>%
  pivot_wider(
    names_from = fine_storm_phase,
    values_from = c(mean_abundance_total, mean_abundance_new, mean_abundance_old, median_abundance_total, median_abundance_new, median_abundance_old),
    id_cols = c(site_id, name, cluster_lon, cluster_lat)
  ) %>%
  mutate(
    # Mean Piece Changes (post-pre subtraction)
    pcs_change_mean_total = round(mean_abundance_total_post - mean_abundance_total_pre, 0), 
    pcs_change_mean_new = round(mean_abundance_new_post - mean_abundance_new_pre, 0),
    pcs_change_mean_old = round(mean_abundance_old_post - mean_abundance_old_pre, 0),
    
    # Median Piece Changes (post-pre subtraction)
    pcs_change_median_total = round(median_abundance_total_post - median_abundance_total_pre, 0), 
    pcs_change_median_new = round(median_abundance_new_post - median_abundance_new_pre, 0),
    pcs_change_median_old = round(median_abundance_old_post - median_abundance_old_pre, 0)
    )


# flag change type (accumulation | depletion) based on strandline
mp_change <- mp_change %>%
  mutate(
    # mean change types (total, new, old)
    change_type_mean_total = case_when(
      pcs_change_mean_total == 0 ~ "no_change",
      pcs_change_mean_total > 0 ~ "accumulation",
      pcs_change_mean_total < 0 ~ "depletion"
    ),
    change_type_mean_new = case_when(
      pcs_change_mean_new == 0 ~ "no_change",
      pcs_change_mean_new > 0 ~ "accumulation",
      pcs_change_mean_new < 0 ~ "depletion"
      ),
    change_type_mean_old = case_when(
      pcs_change_mean_old == 0 ~ "no_change",
      pcs_change_mean_old > 0 ~ "accumulation",
      pcs_change_mean_old < 0 ~ "depletion"
    ),
    # median change types (total, new, old)
    change_type_median_total = case_when(
      pcs_change_median_total == 0 ~ "no_change",
      pcs_change_median_total > 0 ~ "accumulation",
      pcs_change_median_total < 0 ~ "depletion"
    ),
    change_type_median_new = case_when(
      pcs_change_median_new == 0 ~ "no_change",
      pcs_change_median_new > 0 ~ "accumulation",
      pcs_change_median_new < 0 ~ "depletion"
    ),
    change_type_median_old = case_when(
      pcs_change_median_old == 0 ~ "no_change",
      pcs_change_median_old > 0 ~ "accumulation",
      pcs_change_median_old < 0 ~ "depletion"
    ),
  )

write.csv(mp_change, here('output','summary_tables', 'MP_delta_bySite.csv'), row.names = FALSE)


# 2) Daily abundance Deltas
  # Combines storm-influenced & control into one pass via summarize_abundance(),
  # computes day-to-day lag deltas for both mean and median (no percent change).

# 1) Daily abundance, storm-influenced & control together
daily_abund <- master_df %>%
  distinct(site_id, date_nurdle, survey_type, name, fine_storm_phase, .keep_all = TRUE) %>%
  mutate(
    name = if_else(survey_type == "control", "control", name),
    fine_storm_phase = if_else(survey_type == "control", "control", fine_storm_phase)
  ) %>%
  group_by(site_id, date_nurdle, name, fine_storm_phase) %>%
  summarise_abundance() %>%
  left_join(
    master_df %>% distinct(site_id, cluster_lat, cluster_lon),
    by = "site_id"
  )

# 2) Day-to-day changes (lag), mean and median
daily_change <- daily_abund %>%
  arrange(site_id, name, date_nurdle) %>%
  group_by(site_id, name) %>%
  mutate(
    # mean-based lags
    mean_total_prev = lag(mean_abundance_total),
    mean_new_prev = lag(mean_abundance_new),
    mean_old_prev = lag(mean_abundance_old),
    mean_total_diff = mean_abundance_total - mean_total_prev,
    mean_new_diff = mean_abundance_new - mean_new_prev,
    mean_old_diff = mean_abundance_old - mean_old_prev,
    
    # median-based lags
    median_total_prev = lag(median_abundance_total),
    median_new_prev = lag(median_abundance_new),
    median_old_prev = lag(median_abundance_old),
    median_total_diff = median_abundance_total - median_total_prev,
    median_new_diff = median_abundance_new - median_new_prev,
    median_old_diff = median_abundance_old - median_old_prev
  ) %>%
  ungroup()

write.csv(daily_change, here('output', 'summary_tables', 'daily_delta.csv'), row.names = FALSE)
