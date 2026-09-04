# SCRIPT 4: MP Redistribution Analysis — Control vs. Storm-Influenced
# Purpose: Test whether the RATE and DIRECTION of MP change differs between storm-influenced and control periods.

# Load Libraries
# Read Data
# Config
  # set strandline to be analyzed
  # build helper functions
# Section 1: Build Survey Pairs
  # Storm Pairs: consecutive pairs within each site-storm link (linked_id)
  # Control Pairs: consecutive pairs within each site
# Section 2: Build Fine_response Matched Survey Interval Windows
# Section 3: Interval Window Distriution Checks
# Section 4: Analysis: Rate of Change & Directional Change
  # GoM-level: all storm pairs vs all control pairs (capped at observed storm survey max interval = 43 d)
  # Storm-level: all sites (control & storm pairs) linked to each storm
  # Site-level: all control & storm pairs linked to each site

# Key Question: MP movement dynamics: How fast and which direction is MP moving from one survey to the next?

# Load Libraries ----------------------------------------------------------
pacman::p_load('dplyr', 'tidyr', 'ggplot2', 'here')
source(here("functions.R"))


# Read Data ---------------------------------------------------------------
master_df <- read.csv(here("output", "master_df.csv")) %>%
  mutate(date_nurdle = as.Date(date_nurdle))


# Config ------------------------------------------------------------------

# set strand to analyze (old is our primary interest, but could analyze total and new strands)
strand <- "old" # options: old, new, total

alpha <- 0.05

# Section 1: Build Survey Pairs -------------------------------------------------

# 1) Storm Pairs
# de-dup to one row per site-storm link x survey date & detect consecutive survey info ("pairs")
storm_pairs <- master_df %>%
  filter(survey_type == "storm_influenced", fine_storm_phase != "during") %>%
  distinct(linked_id, date_nurdle, .keep_all = TRUE) %>% # linked_id ensures only linked site-storm interactions are paired
  arrange(linked_id, date_nurdle) %>%
  group_by(site_id, name) %>%
  mutate(
    next_date = lead(date_nurdle),
    next_old = lead(weighted_old_strandline),
    next_new = lead(weighted_new_strandline),
    next_total = lead(weighted_total_amount),
    next_fine_response = lead(fine_response),
    next_fine_storm_phase = lead(fine_storm_phase),
    days_between = as.numeric(next_date - date_nurdle),
    delta_old = next_old - weighted_old_strandline,
    delta_new = next_new - weighted_new_strandline,
    delta_total = next_total - weighted_total_amount,
    rate_old = delta_old / days_between,
    rate_new = delta_new / days_between,
    rate_total = delta_total / days_between
  ) 


# label directional change b/w consecutive surveys
storm_pairs <- storm_pairs %>%
  filter(!is.na(next_date), days_between > 0) %>%
  ungroup() %>%
  label_direction() %>%       # use directional helper fn
  select(linked_id, site_id, site_name, name, survey_type, date_nurdle, next_date, days_between,
         fine_storm_phase, next_fine_storm_phase, next_fine_response, delta_old, rate_old, direction_old,
         delta_new, rate_new, direction_new, delta_total, rate_total, direction_total)


# 2) Control Pairs
# no storms: pairs built within each site & detect consecutive survey info ("pairs+)
control_pairs <- master_df %>%
  filter(survey_type == "control") %>%
  distinct(site_id, date_nurdle, .keep_all = TRUE) %>%
  arrange(site_id, date_nurdle) %>%
  group_by(site_id) %>%
  mutate(
    next_date = lead(date_nurdle),
    next_old = lead(weighted_old_strandline),
    next_new = lead(weighted_new_strandline),
    next_total = lead(weighted_total_amount),
    days_between = as.numeric(next_date - date_nurdle),
    delta_old = next_old - weighted_old_strandline,
    delta_new = next_new - weighted_new_strandline,
    delta_total = next_total - weighted_total_amount,
    rate_old = delta_old / days_between,
    rate_new = delta_new / days_between,
    rate_total = delta_total / days_between
  ) 

# label directional change b/w consecutive surveys
control_pairs <- control_pairs %>%
  filter(!is.na(next_date), days_between > 0) %>%
  ungroup() %>%
  label_direction() %>%
  select(site_id, site_name, survey_type, date_nurdle, next_date, days_between,
         delta_old, rate_old, direction_old, delta_new, rate_new, direction_new,
         delta_total, rate_total, direction_total)


# 3) Bind data
survey_pairs <- bind_rows(storm_pairs, control_pairs)

write.csv(survey_pairs, here('output', 'summary_tables', paste0('survey_pairs_', strand, '.csv')), row.names = FALSE)


# Section 2: Build Fine_response Matched Survey Interval Windows ----------------------------------------
# exclusive windows: each window contains ONLY pairs within that temporal range
# control pairs interval-matched to storm window scale

pairs_acute <- bind_rows(
  storm_pairs %>% filter(next_fine_response == "acute"), # storm: 1-4d from day 0
  control_pairs %>% filter(days_between <= 8)  # matches storm ±4d scale
)

pairs_subacute <- bind_rows(
  storm_pairs %>% filter(next_fine_response %in% c("acute", "subacute")), # storm: 5-10d from day 0
  control_pairs %>% filter(days_between > 8 & days_between <= 20)  # ± 5-10d
)

pairs_extended <- bind_rows(
  storm_pairs %>% filter(next_fine_response %in% c("acute", "subacute", "extended")), # storm: 11-28d from day 0
  control_pairs %>% filter(days_between > 20 & days_between <= 56)  # ± 11 - 28d
)


temp_windows <- list(
  "acute" = pairs_acute,
  "subacute" = pairs_subacute,
  "extended" = pairs_extended
)


for (window_name in names(temp_windows)) {
  w <- temp_windows[[window_name]]
  cat(window_name, ": storm =", sum(w$survey_type == "storm_influenced", na.rm = TRUE),
      "| control =", sum(w$survey_type == "control", na.rm = TRUE), "\n")
}



# Section 3: days_between Distribution Check ------------------------------
# are survey intervals comparable between storm and control groups?

interval_summary <- survey_pairs %>%
  group_by(survey_type) %>%
  summarise(
    n = n(),
    median = median(days_between, na.rm = TRUE),
    q25 = quantile(days_between, 0.25, na.rm = TRUE),
    q75 = quantile(days_between, 0.75, na.rm = TRUE),
    mean = round(mean(days_between, na.rm = TRUE), 1),
    max = max(days_between, na.rm = TRUE),
    .groups = "drop"
  )

print(interval_summary) # no, based on mean and max values - really long survey gaps in the control surveys

# but lets test it with stats

# Mann-Whitney U (aka Wilcoxon rank-sum): are the central tendency of the 2 distrib. different? - non-parametric (non-normal) equivalent of a t-test 
  # does one group tend to have a higher value than the other (median inter-survey intervals)
interval_wilcox <- wilcox.test(days_between ~ survey_type, data = survey_pairs, exact = FALSE)

# Kolmogorov-Smirnov test: tests whether 2 distributions are sig diff w/out assuming normality by comparing cumulative distrib functions (ECDFs)
  # D stat = max diff b/w the two ECDFs (0-1). Sig p-value = 2 groups have meaningfully different sampling interval distributions (expected)
    # are the overall shapes of the distrib. different
interval_ks <- ks.test(
  survey_pairs$days_between[survey_pairs$survey_type == "storm_influenced"],
  survey_pairs$days_between[survey_pairs$survey_type == "control"]
)


# compile results
interval_tests <- data.frame(
  test = c("Wilcoxon", "KS"),
  statistic = c(interval_wilcox$statistic, interval_ks$statistic),
  p_value = c(interval_wilcox$p.value, interval_ks$p.value),
  sig = c(interval_wilcox$p.value < 0.05, interval_ks$p.value < 0.05)
)


write.csv(interval_summary, here("output", "GoM_level", "stats", "paired_surveys", "interval_distribution_check.csv"), row.names = FALSE)
write.csv(interval_tests, here("output", "GoM_level", "stats", "paired_surveys", "interval_distribution_tests.csv"), row.names = FALSE)


# medians similar but mean and max differ due to control outliers (up to 1176d), Mann
# cap control at <=56d (storm-influenced period) for overall and phase tests below




# Section 4: Analysis: Rate of Change & Directional Change (strandline level) -----------------------------------------------------

# set max survey interval days equal to storm-influenced max
max_days = 56

# prep strand column
sp <- prep_strand(storm_pairs, strand) # storm pairs
cp <- prep_strand(control_pairs, strand) %>% # cntrl pairs
  filter(days_between <= max_days)


# Analysis 1) GoM-level: all storm pairs vs all control pairs (capped at observed storm survey max interval = 43 d) rate and direction of change
gom_result <- run_rate_dir_tests(bind_rows(sp, cp), group_label = "GoM_rate_dir")

write.csv(gom_result, here("output", "GoM_level", "stats", "rate_and_dir", strand, paste0("GoM_rate_dir_", strand, ".csv")), row.names = FALSE)



# Analysis 2) Per-storm: all sites (control & storm pairs) linked to each storm
storms <- unique(sp$name)

per_storm_results <- do.call(rbind, lapply(storms, function(storm_name) {
  storm_sub <- sp[sp$name == storm_name, ]
  linked_sites <- unique(storm_sub$site_id)
  ctrl_sub <- cp[cp$site_id %in% linked_sites, ]
  run_rate_dir_tests(bind_rows(storm_sub, ctrl_sub), group_label = storm_name)
}))


# Perform Benjamini-Hochberg false discover rate (FDR) multiple test correction
# FDR estimates proportion of discoveries (sig results) that are false positives (i.e. "proportion of tests for which the null hypoth was rejected yet the null hypoth is true" - Analysis of Biological Data Whitlock & Schluter
per_storm_results <- per_storm_results %>%
  mutate(
    rate_p_adj = p.adjust(rate_p,  method = "BH"),
    rate_sig_adj = rate_p_adj < alpha,
    direction_p_adj = p.adjust(direction_p,  method = "BH"),
    direction_sig_adj = direction_p_adj < alpha
  )

write.csv(per_storm_results, here("output", "storm_level", "stats", "rate_and_dir", strand, paste0("per_storm_rate_dir_", strand, ".csv")), row.names = FALSE)


# Analysis 3) Per-site: all control & storm pairs linked to each site
sites <- unique(sp$site_id)

per_site_results <- do.call(rbind, lapply(sites, function(site_id){
  storm_sub <- sp[sp$site_id == site_id, ]
  ctrl_sub <- cp[cp$site_id == site_id, ]
  run_rate_dir_tests(bind_rows(storm_sub, ctrl_sub), group_label = site_id)
}))


# FDR adj
per_site_results <- per_site_results %>%
  mutate(
    rate_p_adj = p.adjust(rate_p,  method = "BH"),
    rate_sig_adj = rate_p_adj < alpha,
    direction_p_adj = p.adjust(direction_p,  method = "BH"),
    direction_sig_adj = direction_p_adj < alpha
  )

write.csv(per_site_results, here("output", "site_level", "stats", "rate_and_dir", strand, paste0("per_site_rate_dir_", strand, ".csv")), row.names = FALSE)

